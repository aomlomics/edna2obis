import pandas as pd
from pygbif import species
import os
import sys
import time
import pickle
import multiprocess as mp
import logging
import re
import json
import urllib.error
import urllib.request

# --- Setup logging ---
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


def _get_mp_context():
    """
    On macOS, forking can crash with Objective-C runtime errors. Using 'spawn' avoids fork()
    and matches WoRMS matching behavior.
    """
    if sys.platform == 'darwin':
        return mp.get_context('spawn')
    return mp.get_context()


def _resolve_gbif_n_proc(n_proc):
    """Match WoRMS: 0 => min(10, cpu_count); explicit value capped at 10; always >= 1."""
    try:
        n_raw = int(n_proc) if n_proc is not None else 0
    except (TypeError, ValueError):
        n_raw = 0
    cpus = mp.cpu_count() or 1
    if n_raw <= 0:
        return max(1, min(10, cpus))
    return max(1, min(int(n_raw), 10))


# --- Constants ---
AMBIGUOUS_KINGDOMS = ['bacteria', 'plantae', 'fungi', 'animalia', 'archaea', 'protista', 'chromista']
DWC_RANKS_STD = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus']
UNACCEPTED_GBIF_STATUSES = {
    'SYNONYM',
    'HETEROTYPIC_SYNONYM',
    'HOMOTYPIC_SYNONYM',
    'INTERMEDIATE_RANK_SYNONYM',
    'PROPARTE_SYNONYM',
    'MISAPPLIED',
    'DETERMINATION_SYNONYM',
}

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
def _normalize_taxon_token(value):
    if value is None or pd.isna(value):
        return ''
    text = str(value).strip().lower()
    text = re.sub(r'\s+', ' ', text)
    return text


def _compute_lineage_consistency(parsed_names, candidate_match):
    if not parsed_names or not isinstance(candidate_match, dict):
        return 0.0, 0, 0

    ignore = {'eukaryota'}
    verbatim_tokens = []
    for token in parsed_names:
        normalized = _normalize_taxon_token(token)
        if normalized and normalized not in ignore:
            verbatim_tokens.append(normalized)

    verbatim_set = set(verbatim_tokens)
    if not verbatim_set:
        return 0.0, 0, 0

    candidate_vals = set()
    scientific_name = _normalize_taxon_token(candidate_match.get('scientificName'))
    if scientific_name:
        candidate_vals.add(scientific_name)

    for rank in DWC_RANKS_STD:
        value = _normalize_taxon_token(candidate_match.get(rank))
        if value:
            candidate_vals.add(value)

    if not candidate_vals:
        return 0.0, 0, len(verbatim_set)

    matched_count = len(verbatim_set.intersection(candidate_vals))
    total_count = len(verbatim_set)
    return matched_count / max(1, total_count), matched_count, total_count


def _extract_usage_key(record):
    if not isinstance(record, dict):
        return None
    return record.get('usageKey') or record.get('key') or record.get('nubKey')


def _extract_accepted_usage_key(record):
    if not isinstance(record, dict):
        return None
    return (
        record.get('acceptedUsageKey')
        or record.get('acceptedKey')
        or record.get('acceptedTaxonKey')
    )


def _gbif_status(record):
    if not isinstance(record, dict):
        return ''
    return str(record.get('taxonomicStatus') or record.get('status') or '').strip().upper()


def _format_environment(record):
    if not isinstance(record, dict):
        return None

    habitat = record.get('habitat')
    if not habitat:
        return None
    if isinstance(habitat, (list, tuple, set)):
        values = [str(value).strip().lower() for value in habitat if str(value).strip()]
        return ';'.join(values) if values else None
    return str(habitat).strip().lower() or None


GBIF_SPECIES_PROFILES_URL = 'https://api.gbif.org/v1/species/{key}/speciesProfiles'


def _truthy_flag(value):
    if value in (True, 1):
        return True
    if isinstance(value, str) and value.strip().lower() in {'1', 'true', 'yes', 'y'}:
        return True
    return False


def _environment_from_species_profiles(usage_key):
    """
    GBIF habitat flags (marine / freshwater / terrestrial) come from the Species Profiles
    extension, not from name_usage. Merge flags across all returned profiles (union).
    """
    if usage_key is None or (isinstance(usage_key, str) and not str(usage_key).strip()):
        return None
    try:
        key_int = int(usage_key)
    except (TypeError, ValueError):
        return None

    url = GBIF_SPECIES_PROFILES_URL.format(key=key_int)
    req = urllib.request.Request(url, headers={'User-Agent': 'edna2obis (GBIF_matching)'})
    try:
        with urllib.request.urlopen(req, timeout=45) as resp:
            payload = json.loads(resp.read().decode())
    except urllib.error.HTTPError as e_http:
        if e_http.code == 404:
            return None
        logging.warning(
            "GBIF API (speciesProfiles) for key '%s' failed: HTTP %s",
            usage_key,
            e_http.code,
        )
        return None
    except Exception as e_profiles:
        logging.warning(
            "GBIF API (speciesProfiles) for key '%s' failed with error: %s",
            usage_key,
            e_profiles,
        )
        return None

    results = payload.get('results') or []
    marine = freshwater = terrestrial = False
    for prof in results:
        if not isinstance(prof, dict):
            continue
        if _truthy_flag(prof.get('marine')):
            marine = True
        if _truthy_flag(prof.get('freshwater')):
            freshwater = True
        if _truthy_flag(prof.get('terrestrial')):
            terrestrial = True

    labels = []
    if marine:
        labels.append('marine')
    if freshwater:
        labels.append('freshwater')
    if terrestrial:
        labels.append('terrestrial')
    return ';'.join(labels) if labels else None


def _format_higher_classification(parents):
    if not parents:
        return None

    names = []
    if isinstance(parents, dict):
        parents = parents.get('results') or parents.get('parents') or []

    for parent in parents:
        if not isinstance(parent, dict):
            continue
        name = parent.get('canonicalName') or parent.get('scientificName')
        if name:
            names.append(str(name).strip())

    return '|'.join([name for name in names if name]) or None


def _get_usage_context(usage_key, usage_context_cache, include_higher_classification=False):
    if not usage_key:
        return {'environment': None, 'higherClassification': None}

    cache_key = (str(usage_key), bool(include_higher_classification))
    if cache_key in usage_context_cache:
        return usage_context_cache[cache_key]

    context = {'environment': None, 'higherClassification': None}
    usage_record = None
    try:
        usage_record = species.name_usage(key=usage_key)
    except Exception as e_usage:
        logging.warning(f"GBIF API call (name_usage) for key '{usage_key}' failed with error: {e_usage}")

    context['environment'] = _environment_from_species_profiles(usage_key)
    if context['environment'] is None and isinstance(usage_record, dict):
        context['environment'] = _format_environment(usage_record)

    if include_higher_classification:
        try:
            parents = species.name_usage(key=usage_key, data='parents')
            context['higherClassification'] = _format_higher_classification(parents)
        except Exception as e_parents:
            logging.warning(f"GBIF API call (parents) for key '{usage_key}' failed with error: {e_parents}")

    usage_context_cache[cache_key] = context
    return context


def _build_result_dict(
    gbif_result,
    match_type='UNKNOWN',
    replaced_unaccepted=False,
    unaccepted_match=False,
    usage_context=None,
):
    """Helper to build the final dictionary from a GBIF result."""
    rank = str(gbif_result.get('rank') or '').lower()
    usage_key = _extract_usage_key(gbif_result)
    usage_context = usage_context or {}
    
    return {
        'scientificName': gbif_result.get('scientificName'),
        'taxonID': f"gbif:{usage_key}" if usage_key else None,
        'kingdom': gbif_result.get('kingdom'),
        'phylum': gbif_result.get('phylum'),
        'class': gbif_result.get('class'),
        'order': gbif_result.get('order'),
        'family': gbif_result.get('family'),
        'genus': gbif_result.get('genus'),
        'taxonRank': rank,
        'confidence': gbif_result.get('confidence'),
        'nameAccordingTo': 'GBIF',
        'match_type_debug': f"GBIF_{gbif_result.get('matchType', match_type)}",
        'replaced_unaccepted': replaced_unaccepted,
        'unaccepted_match': unaccepted_match,
        'environment': usage_context.get('environment') or _format_environment(gbif_result),
        'higherClassification': usage_context.get('higherClassification'),
    }


def _score_and_tag_match(match_dict, cleaned_parts, cleaned_taxonomy):
    match_dict['cleanedTaxonomy'] = cleaned_taxonomy
    score, matched_count, total_count = _compute_lineage_consistency(cleaned_parts, match_dict)
    match_dict['assignment_score'] = score
    match_dict['ranks_matched'] = matched_count
    match_dict['ranks_provided'] = total_count
    return match_dict


def _selection_sort_key(match_dict, use_assignment_score):
    confidence = match_dict.get('confidence')
    confidence = confidence if confidence is not None else -1
    score = match_dict.get('assignment_score')
    score = score if score is not None else -1
    ranks_matched = match_dict.get('ranks_matched')
    ranks_matched = ranks_matched if ranks_matched is not None else -1

    if use_assignment_score:
        return (score, ranks_matched, confidence)
    return (confidence, score, ranks_matched)


def _coerce_bool(value):
    if isinstance(value, str):
        return value.strip().lower() in {'1', 'true', 'yes', 'y'}
    return bool(value)

def _gbif_worker(args):
    """
    Worker function to find ALL accepted matches for a taxon string.
    This function implements a tenacious "species-first", context-aware strategy.
    """
    lookup_key, skip_species, gbif_limit, include_higher_classification = args
    
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
    
    usage_context_cache = {}
    
    # Walk up the taxonomy from most specific to least specific
    for taxon_to_search in reversed(cleaned_parts):
        if not taxon_to_search or len(taxon_to_search) < 3:
            continue
            
        try:
            accepted_by_taxon_id = {}
            unaccepted_rows = []

            def add_accepted_match(record, match_type='LOOKUP_ACCEPTED', replaced_unaccepted=False):
                usage_key = _extract_usage_key(record)
                usage_context = _get_usage_context(
                    usage_key,
                    usage_context_cache,
                    include_higher_classification=include_higher_classification,
                )
                match_dict = _build_result_dict(
                    record,
                    match_type=match_type,
                    replaced_unaccepted=replaced_unaccepted,
                    unaccepted_match=False,
                    usage_context=usage_context,
                )
                match_dict = _score_and_tag_match(match_dict, cleaned_parts, cleaned_taxonomy)

                existing = accepted_by_taxon_id.get(match_dict['taxonID'])
                if existing:
                    existing['replaced_unaccepted'] = existing.get('replaced_unaccepted') or replaced_unaccepted
                    if existing.get('confidence') is None and match_dict.get('confidence') is not None:
                        existing['confidence'] = match_dict['confidence']
                    return

                accepted_by_taxon_id[match_dict['taxonID']] = match_dict

            def add_unaccepted_match(record, status):
                usage_key = _extract_usage_key(record)
                usage_context = _get_usage_context(
                    usage_key,
                    usage_context_cache,
                    include_higher_classification=include_higher_classification,
                )
                match_dict = _build_result_dict(
                    record,
                    match_type=f'LOOKUP_{status or "UNACCEPTED"}',
                    replaced_unaccepted=False,
                    unaccepted_match=True,
                    usage_context=usage_context,
                )
                match_dict = _score_and_tag_match(match_dict, cleaned_parts, cleaned_taxonomy)
                unaccepted_rows.append(match_dict)

            # For single-word kingdom names, specify the rank to avoid homonyms. (stick bug named 'Bacteria')
            api_params = {'q': taxon_to_search, 'status': 'ACCEPTED', 'limit': gbif_limit, 'higherTaxonKey': None}
            if len(cleaned_parts) == 1 and taxon_to_search.lower() in AMBIGUOUS_KINGDOMS:
                api_params['rank'] = 'kingdom'

            # Use name_lookup which can return multiple results
            lookup_data = species.name_lookup(**api_params)
            
            if lookup_data and lookup_data.get('results'):
                for record in lookup_data['results']:
                    if record.get('key'):
                        # Use pygbif's name_backbone to perform a proper match and get a confidence score.
                        # This requires the 'name' parameter, and is the only way to get a confidence value.
                        try:
                            full_record = species.name_backbone(name=record['scientificName'], rank=record.get('rank'), strict=True)
                            if full_record.get('usageKey'):
                                add_accepted_match(full_record, match_type='LOOKUP_ACCEPTED')
                        except Exception as e_backbone:
                            logging.warning(f"GBIF API call (name_backbone) for '{record.get('scientificName')}' failed with error: {e_backbone}")
                            continue

            # A second lookup without the ACCEPTED filter exposes synonym/unaccepted rows for review.
            all_status_params = {'q': taxon_to_search, 'limit': gbif_limit, 'higherTaxonKey': None}
            if len(cleaned_parts) == 1 and taxon_to_search.lower() in AMBIGUOUS_KINGDOMS:
                all_status_params['rank'] = 'kingdom'

            all_status_data = species.name_lookup(**all_status_params)
            if all_status_data and all_status_data.get('results'):
                for record in all_status_data['results']:
                    if not record.get('key'):
                        continue

                    status = _gbif_status(record)
                    if status not in UNACCEPTED_GBIF_STATUSES:
                        continue

                    add_unaccepted_match(record, status)

                    accepted_key = _extract_accepted_usage_key(record)
                    if accepted_key:
                        try:
                            accepted_record = species.name_usage(key=accepted_key)
                            if accepted_record:
                                add_accepted_match(
                                    accepted_record,
                                    match_type=f'LOOKUP_REPLACED_{status}',
                                    replaced_unaccepted=True,
                                )
                        except Exception as e_accepted:
                            logging.warning(f"GBIF API call (accepted usage) for key '{accepted_key}' failed with error: {e_accepted}")

            # If we found any matches, stop the walk-up at the most specific usable term.
            all_matches = list(accepted_by_taxon_id.values()) + unaccepted_rows
            if all_matches:
                return args, all_matches
                    
        except Exception as e:
            print(f"FAIL: General error during name_lookup for '{taxon_to_search}': {e}")
            logging.warning(f"GBIF API call (walk-up/name_lookup) for '{taxon_to_search}' failed with error: {e}")
            continue

    return args, []


def get_gbif_match_for_dataframe(occurrence_df, params_dict, n_proc=0):
    """
    Main function to perform taxonomic matching against the GBIF backbone.
    This function now returns both a main DataFrame with the highest-confidence match
    and an info DataFrame with all possible matches for ambiguous cases.
    """
    if occurrence_df.empty:
        logging.warning("Input DataFrame is empty. Skipping GBIF matching.")
        return {'main_df': occurrence_df, 'info_df': pd.DataFrame()}

    n_proc = _resolve_gbif_n_proc(n_proc)
    params_dict['gbif_n_proc_effective'] = n_proc
    logging.info("Using %s processes for GBIF name matching (0 in config = auto, max 10)", n_proc)
    if n_proc >= 8:
        logging.warning(
            "GBIF: %s parallel workers — may hit API rate limits or errors", n_proc
        )

    gbif_limit = params_dict.get('gbif_match_limit', 3)
    use_assignment_score_for_selection = _coerce_bool(
        params_dict.get('gbif_use_assignment_score_for_selection', False)
    )
    include_higher_classification = _coerce_bool(
        params_dict.get('gbif_return_higher_classification', False)
    )

    output_dir = params_dict.get('output_dir', '.')
    os.makedirs(output_dir, exist_ok=True)
    # v3: environment from GBIF speciesProfiles; bump file when match/env logic changes so stale cache is not reused.
    cache_file = os.path.join(output_dir, 'gbif_matches_cache_v3.pkl')
    
    # --- Load existing cache if it exists ---
    try:
        with open(cache_file, 'rb') as f:
            cache = pickle.load(f)
        logging.info(f"Loaded {len(cache)} results from existing GBIF cache file: {cache_file}")
    except (FileNotFoundError, EOFError):
        cache = {}
        logging.info("No GBIF cache file found or cache is empty. Starting fresh.")

    assays_to_skip_species = params_dict.get('assays_to_skip_species_match', [])
    if assays_to_skip_species is None:
        assays_to_skip_species = []

    occurrence_df['skip_species_flag'] = occurrence_df['assay_name'].isin(assays_to_skip_species)
    occurrence_df['gbif_limit'] = gbif_limit # Add the limit to the dataframe for the worker
    occurrence_df['gbif_return_higher_classification'] = include_higher_classification
    unique_lookups = [
        tuple(row)
        for row in occurrence_df[
            [
                'verbatimIdentification',
                'skip_species_flag',
                'gbif_limit',
                'gbif_return_higher_classification',
            ]
        ].drop_duplicates().to_numpy()
    ]
    
    # --- Identify which lookups are new and need to be fetched ---
    new_lookups = [lookup for lookup in unique_lookups if lookup not in cache]
    
    logging.info(f"Found {len(unique_lookups)} unique taxonomic strings in this run. {len(new_lookups)} are new and will be fetched from the GBIF API.")

    if new_lookups:
        total_lookups = len(new_lookups)
        logging.info(f"Starting parallel GBIF matching for {total_lookups} lookups with {n_proc} processes...")
        start_time = time.time()
        
        processed_count = 0
        ctx = _get_mp_context()
        with ctx.Pool(processes=n_proc) as pool:
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
        logging.info(f"Saved/updated cache with {len(new_lookups)} new results. Cache now contains {len(cache)} total results.")

    # --- Process Results and Create DataFrames ---
    main_df_records = []
    info_df_records = []
    
    # Handle incertae sedis cases first
    incertae_sedis_record = {
        'scientificName': 'incertae sedis', 'taxonID': None, 'kingdom': None, 'phylum': None,
        'class': None, 'order': None, 'family': None, 'genus': None, 'taxonRank': None,
        'confidence': None, 'nameAccordingTo': 'GBIF', 'environment': pd.NA,
        'higherClassification': pd.NA, 'replaced_unaccepted': False,
        'unaccepted_match': False, 'assignment_score': pd.NA,
        'ranks_matched': pd.NA, 'ranks_provided': pd.NA
    }

    # Map results back to the original dataframe structure
    for verbatim_id, skip_flag, limit, _include_higher_classification in unique_lookups:
        key = (verbatim_id, skip_flag, limit, _include_higher_classification)
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
            # Sort by confidence by default; optionally use lineage score first when configured.
            selection_candidates = [match for match in matches if not match.get('unaccepted_match')]
            if not selection_candidates:
                selection_candidates = matches
            sorted_matches = sorted(
                selection_candidates,
                key=lambda x: _selection_sort_key(x, use_assignment_score_for_selection),
                reverse=True,
            )
            best_match = sorted_matches[0]
            
            main_df_records.append({'verbatimIdentification': verbatim_id, 'skip_species_flag': skip_flag, **best_match})
            
            is_ambiguous = len(matches) > 1
            
            for match in matches:
                info_record = {'verbatimIdentification': verbatim_id, **match, 'ambiguous': is_ambiguous}
                info_df_records.append(info_record)

    # --- Create Final DataFrames ---
    main_results_df = pd.DataFrame(main_df_records)
    info_df = pd.DataFrame(info_df_records)
    
    # Merge best match results back into the main occurrence DataFrame
    cols_to_drop = [
        'scientificName', 'taxonID', 'kingdom', 'phylum', 'class', 'order', 'family', 
        'genus', 'taxonRank', 'confidence', 'nameAccordingTo', 'match_type_debug',
        'cleanedTaxonomy', 'environment', 'higherClassification', 'replaced_unaccepted',
        'unaccepted_match', 'assignment_score', 'ranks_matched', 'ranks_provided'
    ]
    occurrence_df_clean = occurrence_df.drop(columns=[col for col in cols_to_drop if col in occurrence_df.columns], errors='ignore')
    
    final_df = occurrence_df_clean.merge(
        main_results_df,
        on=['verbatimIdentification', 'skip_species_flag'],
        how='left'
    )
    
    final_df.drop(columns=['skip_species_flag'], inplace=True, errors='ignore')
    final_df.drop(columns=['gbif_return_higher_classification'], inplace=True, errors='ignore')
    
    logging.info("GBIF matching and merging complete.")
    return {'main_df': final_df, 'info_df': info_df}