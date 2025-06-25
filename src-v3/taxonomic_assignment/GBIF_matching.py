import pandas as pd
from pygbif import species
import os
import time
import pickle
from multiprocess import Pool, cpu_count
import logging
import re

# --- Setup logging ---
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def parse_semicolon_taxonomy(tax_string):
    """
    Helper function to parse and aggressively clean a semicolon-separated taxonomy string.
    This function handles prefixes (e.g., d__), underscores, and other non-standard formatting.
    """
    if pd.isna(tax_string) or not str(tax_string).strip():
        return []

    # Split the raw string by semicolons
    parts = str(tax_string).split(';')

    cleaned_parts = []
    for part in parts:
        # Strip leading/trailing whitespace from the part
        clean_part = part.strip()
        # Remove prefixes like d__, p__, c__, f__, g__, s__, o__
        clean_part = re.sub(r'^[dpcofgs]__', '', clean_part)
        # Replace underscores, hyphens, and slashes with spaces
        clean_part = clean_part.replace('_', ' ').replace('-', ' ').replace('/', ' ')
        # Remove any lingering "sp." or "spp." at the end of a name
        clean_part = re.sub(r'\s+spp?\.?$', '', clean_part, flags=re.IGNORECASE).strip()
        # Clean up any resulting extra whitespace
        clean_part = re.sub(r'\s+', ' ', clean_part).strip()

        # Add the part if it's not empty and not an uninformative term
        if clean_part and clean_part.lower() not in ['unassigned', 'unknown', 'no hit', '']:
            cleaned_parts.append(clean_part)
            
    return cleaned_parts

# --- Helper Functions ---

def _build_result_dict(gbif_result):
    """Helper to build the final dictionary from a GBIF result, cleaning the scientific name."""
    rank = gbif_result.get('rank', '').lower()
    
    # Logic for the scientificName field
    final_name = None
    if rank == 'species':
        # For a species match, use the clean binomial name
        final_name = gbif_result.get('species')
    else:
        # For higher ranks, use the clean rank name (e.g., 'Clausocalanus')
        # instead of the canonical name with authorship.
        final_name = gbif_result.get(rank)
    
    # Fallback if the clean name wasn't found (for unusual ranks like 'subfamily')
    if not final_name:
        final_name = gbif_result.get('scientificName')

    gbif_taxon_id = f"gbif:{gbif_result['usageKey']}"

    # Build the dictionary. We DO NOT include a 'species' key in the output dict.
    result = {
        'scientificName': final_name,
        'scientificNameID': None,
        'taxonID': gbif_taxon_id,
        'kingdom': gbif_result.get('kingdom'),
        'phylum': gbif_result.get('phylum'),
        'class': gbif_result.get('class'),
        'order': gbif_result.get('order'),
        'family': gbif_result.get('family'),
        'genus': gbif_result.get('genus'),
        'taxonRank': rank,
        'nameAccordingTo': 'GBIF',
        'match_type_debug': f"GBIF_{gbif_result.get('matchType', 'UNKNOWN')}"
    }
    return result

def _gbif_worker(args):
    """
    Worker function to match a single taxon name against the GBIF backbone.
    This function implements a tenacious "species-first", context-aware strategy.
    """
    lookup_key, skip_species = args

    if not isinstance(lookup_key, str) or not lookup_key.strip():
        return args, {'scientificName': 'No Match: Empty Input', 'match_type_debug': 'empty_input'}

    # 1. Parse and clean the taxonomy string using the comprehensive function
    # This logic is identical to the WoRMS script for consistency
    cleaned_parts = parse_semicolon_taxonomy(lookup_key)
    
    if not cleaned_parts:
        return args, {'scientificName': 'No Match: Empty Input after processing', 'match_type_debug': 'empty_input_processed'}

    # 2. Determine a kingdom hint for context
    kingdom_hint = None
    if cleaned_parts:
        first_rank_lower = cleaned_parts[0].lower()
        if first_rank_lower in ['bacteria', 'eukaryota', 'archaea', 'viruses', 'chromista', 'protozoa', 'animalia', 'plantae', 'fungi']:
            kingdom_hint = cleaned_parts[0]

    # 3. "Species First" Strategy: Tenaciously try to match the most specific part as a SPECIES
    potential_species_name = cleaned_parts[-1]
    if not skip_species and ' ' in potential_species_name:
        try:
            # Try a strict species-only match first
            gbif_result = species.name_backbone(name=potential_species_name, rank='species', kingdom=kingdom_hint, strict=True)
            if gbif_result.get('usageKey') is None:
                # If strict fails, try a fuzzy species match
                gbif_result = species.name_backbone(name=potential_species_name, rank='species', kingdom=kingdom_hint, strict=False)

            if gbif_result.get('usageKey') is not None and gbif_result.get('rank', '').lower() == 'species':
                return args, _build_result_dict(gbif_result)
        except Exception as e:
            logging.warning(f"GBIF API call (species search) failed for '{potential_species_name}': {e}")
            
    # 4. If species match failed/skipped, proceed with a more robust taxonomic walk-up.
    walk_up_parts = cleaned_parts
    if not skip_species and ' ' in potential_species_name:
        walk_up_parts = cleaned_parts[:-1]

    for taxon_to_search in reversed(walk_up_parts):
        # Skip very short, ambiguous terms.
        if not taxon_to_search or len(taxon_to_search) < 3:
            continue

        try:
            # --- PRIMARY ATTEMPT: Use name_backbone without a restrictive kingdom hint. ---
            # This is often the most effective way to find the correct higher taxon.
            # It allows GBIF's backbone to resolve the name across the entire tree of life.
            gbif_result = species.name_backbone(name=taxon_to_search, strict=False)

            # Check for a high-confidence match from the backbone service.
            # A 'NONE' matchType means GBIF found nothing conclusive.
            if gbif_result.get('matchType') != 'NONE':
                logging.info(f"Walk-up success (backbone) for '{taxon_to_search}' with type: {gbif_result.get('matchType')}")
                return args, _build_result_dict(gbif_result)

            # --- FALLBACK ATTEMPT: Use name_lookup for broader search. ---
            # This is useful if the name is ambiguous or not in the backbone as expected.
            # We look for the first 'ACCEPTED' name.
            logging.info(f"Backbone failed for '{taxon_to_search}'. Trying name_lookup fallback.")
            lookup_data = species.name_lookup(q=taxon_to_search, limit=5)
            if lookup_data and lookup_data.get('data'):
                for record in lookup_data['data']:
                    if record.get('status') == 'ACCEPTED' and record.get('rank'):
                        # Found an accepted name. Now get its full classification using its key.
                        # The species.usage() call retrieves the complete record.
                        usage_key = record.get('key')
                        if usage_key:
                            full_record = species.usage(key=usage_key)
                            if full_record and full_record.get('usageKey'):
                                logging.info(f"Walk-up success (lookup) for '{taxon_to_search}'")
                                # Manually add a match type for our debugging.
                                full_record['matchType'] = 'LOOKUP_ACCEPTED'
                                return args, _build_result_dict(full_record)

        except Exception as e:
            logging.warning(f"GBIF API call (walk-up) for '{taxon_to_search}' failed with error: {e}")
            continue

    # 5. If all attempts fail
    return args, {'scientificName': 'No GBIF Match', 'match_type_debug': 'No_GBIF_Match'}


def get_gbif_match_for_dataframe(occurrence_df, params_dict, n_proc=0):
    """
    Main function to perform taxonomic matching against the GBIF backbone for an
    entire DataFrame of occurrences.
    Args:
        occurrence_df (pd.DataFrame): The DataFrame containing occurrence data.
                                      Must have 'verbatimIdentification' and 'assay_name' columns.
        params_dict (dict): A dictionary of parameters, including:
                            - 'user_defined_assays_to_skip_species' (list of assay names)
                            - 'output_dir' (str, path to save cache file)
        n_proc (int): Number of processes to use for parallel matching. 0 means use all available.
    Returns:
        pd.DataFrame: The input DataFrame with added/updated taxonomic columns.
    """
    if occurrence_df.empty:
        logging.warning("Input DataFrame is empty. Skipping GBIF matching.")
        return occurrence_df

    # --- Caching Setup ---
    output_dir = params_dict.get('output_dir', '.')
    os.makedirs(output_dir, exist_ok=True)
    cache_file = os.path.join(output_dir, 'gbif_matches.pkl')
    
    try:
        with open(cache_file, 'rb') as f:
            cache = pickle.load(f)
        logging.info(f"Loaded {len(cache)} results from GBIF cache file: {cache_file}")
    except (FileNotFoundError, EOFError):
        cache = {}
        logging.info("No GBIF cache file found or cache is empty. Starting fresh.")

    # --- Prepare for Matching ---
    # Fix the parameter name to match what the notebook passes
    assays_to_skip_species = params_dict.get('assays_to_skip_species_match', [])
    
    # Handle case where parameter is None
    if assays_to_skip_species is None:
        assays_to_skip_species = []
    
    # Create unique combinations of verbatimIdentification and whether to skip species
    occurrence_df['skip_species_flag'] = occurrence_df['assay_name'].isin(assays_to_skip_species)
    unique_lookups = occurrence_df[['verbatimIdentification', 'skip_species_flag']].drop_duplicates()
    
    # Filter out what's already in the cache
    # The cache key is a tuple (verbatimIdentification, skip_species_flag)
    new_lookups = [
        tuple(row) for row in unique_lookups.itertuples(index=False)
        if tuple(row) not in cache
    ]

    logging.info(f"Found {len(new_lookups)} new unique taxonomic strings to match against GBIF.")

    # --- Parallel Processing ---
    if new_lookups:
        if n_proc == 0:
            n_proc = cpu_count()
        
        logging.info(f"Starting parallel GBIF matching with {n_proc} processes...")
        start_time = time.time()
        
        with Pool(processes=n_proc) as pool:
            results = pool.map(_gbif_worker, new_lookups)
        
        # Update cache with new results
        for key_tuple, result_dict in results:
            if key_tuple is not None:
                cache[key_tuple] = result_dict
            
        end_time = time.time()
        logging.info(f"Finished matching in {end_time - start_time:.2f} seconds.")

        # Save updated cache
        with open(cache_file, 'wb') as f:
            pickle.dump(cache, f)
        logging.info(f"Saved {len(cache)} total results to GBIF cache.")

    # --- Merge Results back to DataFrame ---
    logging.info("Merging GBIF results back into the main DataFrame...")

    # Create a mapping DataFrame from the cache
    cache_df_list = []
    for (verbatim_id, skip_flag), result in cache.items():
        row = {'verbatimIdentification': verbatim_id, 'skip_species_flag': skip_flag, **result}
        cache_df_list.append(row)
    
    match_results_df = pd.DataFrame(cache_df_list)

    # Merge the results back into the main df
    # First, drop any old taxonomy columns to prevent merge conflicts
    cols_to_drop = ['scientificName', 'scientificNameID', 'taxonID', 'kingdom', 'phylum', 'class', 
                    'order', 'family', 'genus', 'taxonRank', 'nameAccordingTo', 'match_type_debug']
    occurrence_df_clean = occurrence_df.drop(columns=[col for col in cols_to_drop if col in occurrence_df.columns])
    
    final_df = occurrence_df_clean.merge(
        match_results_df,
        on=['verbatimIdentification', 'skip_species_flag'],
        how='left'
    )
    
    final_df.drop(columns=['skip_species_flag'], inplace=True, errors='ignore')
    
    logging.info("GBIF matching and merging complete.")
    return final_df