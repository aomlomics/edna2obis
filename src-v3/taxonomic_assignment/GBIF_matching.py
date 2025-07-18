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
def _build_result_dict(gbif_result, match_type='UNKNOWN'):
    """Helper to build the final dictionary from a GBIF result."""
    rank = gbif_result.get('rank', '').lower()
    
    return {
        'scientificName': gbif_result.get('scientificName'),
        'taxonID': f"gbif:{gbif_result['usageKey']}",
        'kingdom': gbif_result.get('kingdom'),
        'phylum': gbif_result.get('phylum'),
        'class': gbif_result.get('class'),
        'order': gbif_result.get('order'),
        'family': gbif_result.get('family'),
        'genus': gbif_result.get('genus'),
        'taxonRank': rank,
        'confidence': gbif_result.get('confidence'),
        'nameAccordingTo': 'GBIF',
        'match_type_debug': f"GBIF_{gbif_result.get('matchType', match_type)}"
    }

def _gbif_worker(args):
    """
    Worker function to find ALL accepted matches for a taxon string.
    This function implements a tenacious "species-first", context-aware strategy.
    """
    lookup_key, skip_species = args
    
    if not isinstance(lookup_key, str) or not lookup_key.strip():
        return args, []

    cleaned_parts = parse_semicolon_taxonomy(lookup_key)
    if not cleaned_parts:
        return args, []

    if skip_species and len(cleaned_parts) >= 6:
        cleaned_parts = cleaned_parts[:-1]
        if not cleaned_parts:
            return args, []
    
    cleaned_taxonomy = ';'.join(cleaned_parts)
    
    all_accepted_matches = []
    
    # Walk up the taxonomy from most specific to least specific
    for taxon_to_search in reversed(cleaned_parts):
        if not taxon_to_search or len(taxon_to_search) < 3:
            continue
            
        try:
            # Use name_lookup which can return multiple results
            # Change the 'limit' parameter to set how many matches per unique taxon string are returned. The code will select the highest confidence match
            lookup_data = species.name_lookup(q=taxon_to_search, status='ACCEPTED', limit=5, higherTaxonKey=None)
            
            if lookup_data and lookup_data.get('results'):
                for record in lookup_data['results']:
                    if record.get('key'):
                        # Use pygbif's name_backbone to perform a proper match and get a confidence score.
                        # This requires the 'name' parameter, and is the only way to get a confidence value.
                        try:
                            full_record = species.name_backbone(name=record['scientificName'], rank=record.get('rank'), strict=True)
                            if full_record.get('usageKey'):
                                match_dict = _build_result_dict(full_record, match_type='LOOKUP_ACCEPTED')
                                match_dict['cleanedTaxonomy'] = cleaned_taxonomy
                                all_accepted_matches.append(match_dict)
                        except Exception as e_backbone:
                            logging.warning(f"GBIF API call (name_backbone) for '{record.get('scientificName')}' failed with error: {e_backbone}")
                            continue
                
                # If we found any matches, we stop the walk-up, as we've matched the most specific term possible.
                if all_accepted_matches:
                    return args, all_accepted_matches
                    
        except Exception as e:
            print(f"FAIL: General error during name_lookup for '{taxon_to_search}': {e}")
            logging.warning(f"GBIF API call (walk-up/name_lookup) for '{taxon_to_search}' failed with error: {e}")
            continue

    return args, all_accepted_matches


def get_gbif_match_for_dataframe(occurrence_df, params_dict, n_proc=0):
    """
    Main function to perform taxonomic matching against the GBIF backbone.
    This function now returns both a main DataFrame with the highest-confidence match
    and an info DataFrame with all possible matches for ambiguous cases.
    """
    if occurrence_df.empty:
        logging.warning("Input DataFrame is empty. Skipping GBIF matching.")
        return {'main_df': occurrence_df, 'info_df': pd.DataFrame()}

    output_dir = params_dict.get('output_dir', '.')
    os.makedirs(output_dir, exist_ok=True)
    cache_file = os.path.join(output_dir, 'gbif_matches_cache.pkl')
    
    try:
        with open(cache_file, 'rb') as f:
            cache = pickle.load(f)
        logging.info(f"Loaded {len(cache)} results from GBIF cache file: {cache_file}")
    except (FileNotFoundError, EOFError):
        cache = {}
        logging.info("No GBIF cache file found or cache is empty. Starting fresh.")

    assays_to_skip_species = params_dict.get('assays_to_skip_species_match', [])
    if assays_to_skip_species is None:
        assays_to_skip_species = []

    occurrence_df['skip_species_flag'] = occurrence_df['assay_name'].isin(assays_to_skip_species)
    unique_lookups = [tuple(row) for row in occurrence_df[['verbatimIdentification', 'skip_species_flag']].drop_duplicates().to_numpy()]
    
    new_lookups = [lookup for lookup in unique_lookups if lookup not in cache]

    logging.info(f"Found {len(new_lookups)} new unique taxonomic strings to match against GBIF.")

    if new_lookups:
        if n_proc == 0:
            n_proc = cpu_count()
        
        total_lookups = len(new_lookups)
        logging.info(f"Starting parallel GBIF matching for {total_lookups} lookups with {n_proc} processes...")
        start_time = time.time()
        
        processed_count = 0
        with Pool(processes=n_proc) as pool:
            # Use imap_unordered to get results as they complete, allowing for progress reporting
            for key_tuple, result_list in pool.imap_unordered(_gbif_worker, new_lookups):
                if key_tuple is not None:
                    cache[key_tuple] = result_list
                
                processed_count += 1
                # Log progress every 20 lookups or on the last one to avoid spamming the console
                if processed_count % 20 == 0 or processed_count == total_lookups:
                    print(f"  ... Progress: {processed_count}/{total_lookups} lookups completed.", flush=True)
            
        end_time = time.time()
        logging.info(f"Finished matching in {end_time - start_time:.2f} seconds.")

        with open(cache_file, 'wb') as f:
            pickle.dump(cache, f)
        logging.info(f"Saved {len(cache)} total results to GBIF cache.")

    # --- Process Results and Create DataFrames ---
    main_df_records = []
    info_df_records = []
    
    # Handle incertae sedis cases first
    incertae_sedis_record = {
        'scientificName': 'incertae sedis', 'taxonID': None, 'kingdom': None, 'phylum': None,
        'class': None, 'order': None, 'family': None, 'genus': None, 'taxonRank': None,
        'confidence': None, 'nameAccordingTo': 'GBIF'
    }

    # Map results back to the original dataframe structure
    for verbatim_id, skip_flag in unique_lookups:
        key = (verbatim_id, skip_flag)
        matches = cache.get(key, [])
        
        cleaned_taxonomy = ''
        if matches:
            cleaned_taxonomy = matches[0].get('cleanedTaxonomy', '')

        is_ambiguous = len(matches) > 1
        
        if not matches:
            # No match found, use incertae sedis
            best_match = {
                **incertae_sedis_record, 'match_type_debug': 'No_GBIF_Match',
                'cleanedTaxonomy': ';'.join(parse_semicolon_taxonomy(verbatim_id))
            }
            info_record = {'verbatimIdentification': verbatim_id, **best_match, 'ambiguous': False}
            
            main_df_records.append({'verbatimIdentification': verbatim_id, 'skip_species_flag': skip_flag, **best_match})
            info_df_records.append(info_record)
        else:
            # Matches found, determine best one and add all to info
            # Sort by confidence descending, highest first. NaNs (None) go last.
            sorted_matches = sorted(matches, key=lambda x: x['confidence'] if x['confidence'] is not None else -1, reverse=True)
            best_match = sorted_matches[0]
            
            main_df_records.append({'verbatimIdentification': verbatim_id, 'skip_species_flag': skip_flag, **best_match})
            
            for match in matches:
                info_record = {'verbatimIdentification': verbatim_id, **match, 'ambiguous': is_ambiguous}
                info_df_records.append(info_record)

    # --- Create Final DataFrames ---
    main_results_df = pd.DataFrame(main_df_records)
    info_df = pd.DataFrame(info_df_records)
    
    # Merge best match results back into the main occurrence DataFrame
    cols_to_drop = [
        'scientificName', 'taxonID', 'kingdom', 'phylum', 'class', 'order', 'family', 
        'genus', 'taxonRank', 'confidence', 'nameAccordingTo', 'match_type_debug', 'cleanedTaxonomy'
    ]
    occurrence_df_clean = occurrence_df.drop(columns=[col for col in cols_to_drop if col in occurrence_df.columns], errors='ignore')
    
    final_df = occurrence_df_clean.merge(
        main_results_df,
        on=['verbatimIdentification', 'skip_species_flag'],
        how='left'
    )
    
    final_df.drop(columns=['skip_species_flag'], inplace=True, errors='ignore')
    
    logging.info("GBIF matching and merging complete.")
    return {'main_df': final_df, 'info_df': info_df}