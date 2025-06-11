import pyworms
import pandas as pd
import multiprocess as mp
import time
from functools import partial
import logging

# Set up logging to provide clear progress updates
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Standard Darwin Core ranks used for structuring the output
DWC_RANKS_STD = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']

def parse_semicolon_taxonomy(tax_string):
    """
    Helper function to parse and clean a semicolon-separated taxonomy string.
    It replaces underscores with spaces and filters out non-informative terms.
    """
    if pd.isna(tax_string) or not str(tax_string).strip():
        return []
    # Correctly clean names for lookup by replacing underscores
    cleaned_string = str(tax_string).replace('_', ' ')
    return [name.strip() for name in cleaned_string.split(';') if name.strip() and name.strip().lower() not in ['unassigned']]

def get_worms_classification_by_id_worker(aphia_id_to_check, api_source_for_record='WoRMS'):
    """Fetches and formats a full WoRMS record using a direct AphiaID."""
    start_time = time.time()
    try:
        record = pyworms.aphiaRecordByAphiaID(aphia_id_to_check)
        end_time = time.time()
        
        if record and isinstance(record, dict) and record.get('status') == 'accepted':
            result = {
                'scientificName': record.get('scientificname'), 'scientificNameID': record.get('lsid'),
                'taxonRank': record.get('rank'), 'nameAccordingTo': api_source_for_record,
                'match_type_debug': f'Success_AphiaID_{aphia_id_to_check}',
                'api_call_count': 1, 'total_api_time': end_time - start_time
            }
            for rank_std in DWC_RANKS_STD:
                result[rank_std] = record.get(rank_std.lower())
            return result
    except Exception as e:
        logging.warning(f"API call failed for AphiaID {aphia_id_to_check}: {e}")
    
    return {'api_call_count': 1, 'total_api_time': time.time() - start_time, 'match_type_debug': f'Failure_AphiaID_{aphia_id_to_check}'}


def get_worms_match_for_single_taxon_worker(combo_input_for_worker, assays_to_skip_species_list, api_source_for_record='WoRMS'):
    """(FALLBACK) Matches a single verbatimIdentification string using a slow, iterative search."""
    verbatim_identification_str, assay_name_str = combo_input_for_worker
    api_call_count = 0
    total_api_time = 0.0

    parsed_names = parse_semicolon_taxonomy(verbatim_identification_str)
    if not parsed_names:
        return {'match_type_debug': 'No_Tax_String_Provided', 'api_call_count': 0, 'total_api_time': 0.0}

    for i, term in enumerate(reversed(parsed_names)):
        # Correctly check if this is the species-level term (the first from the right)
        if i == 0 and assay_name_str in assays_to_skip_species_list:
            continue
        
        try:
            time.sleep(0.001)
            start_time = time.time()
            s_match_list = pyworms.aphiaRecordsByName(term, like=False, marine_only=True)
            end_time = time.time()
            api_call_count += 1
            total_api_time += (end_time - start_time)
            
            if s_match_list:
                accepted_match = next((m for m in s_match_list if m and m.get('status') == 'accepted'), None)
                if accepted_match:
                    result_dict = {'scientificName': accepted_match.get('scientificname'), 'scientificNameID': accepted_match.get('lsid'), 'taxonRank': accepted_match.get('rank'), 'nameAccordingTo': api_source_for_record, 'match_type_debug': f'Success_Fallback_{term}', 'api_call_count': api_call_count, 'total_api_time': total_api_time}
                    for rank_std in DWC_RANKS_STD:
                        result_dict[rank_std] = accepted_match.get(rank_std.lower())
                    return result_dict
        except Exception as e:
            logging.warning(f"Fallback API call error for '{term}': {e}")
    
    return {'match_type_debug': 'Failed_Fallback_IncertaeSedis', 'api_call_count': api_call_count, 'total_api_time': total_api_time}


def get_worms_match_for_dataframe(occurrence_df, params_dict, n_proc=0):
    """Adds WoRMS taxonomic information using an optimized multi-stage process with corrected skip-species logic."""
    overall_start_time = time.time()
    api_source = params_dict.get('taxonomic_api_source', 'WoRMS')
    assays_to_skip_species = params_dict.get('assays_to_skip_species_match', [])
    pr2_dict = params_dict.get('pr2_worms_dict', {})
    assay_rank_info = params_dict.get('assay_rank_info', {}) # <-- New: get rank info

    df_to_process = occurrence_df.copy()
    df_to_process['verbatimIdentification'] = df_to_process['verbatimIdentification'].astype(str).fillna('')
    df_to_process['_map_key'] = list(zip(df_to_process['verbatimIdentification'], df_to_process['assay_name']))
    unique_tuples = list(df_to_process['_map_key'].drop_duplicates())
    unique_tuples = [combo for combo in unique_tuples if str(combo[0]).strip()]
    logging.info(f"Found {len(unique_tuples)} unique combinations to process.")

    results_cache = {}
    unmatched_tuples = []
    # (Performance counters remain)

    # --- Stage 1: AphiaID Pre-matching (Now works correctly) ---
    if pr2_dict:
        logging.info("Starting Stage 1: AphiaID pre-matching.")
        for verbatim_str, assay_name in unique_tuples:
            parsed_names = parse_semicolon_taxonomy(verbatim_str)
            if parsed_names:
                species_name = parsed_names[-1]
                if species_name in pr2_dict:
                    aphia_id = pr2_dict[species_name]
                    result = get_worms_classification_by_id_worker(aphia_id, api_source)
                    if 'scientificName' in result:
                        results_cache[(verbatim_str, assay_name)] = result
                        continue
            unmatched_tuples.append((verbatim_str, assay_name))
        logging.info(f"Finished Stage 1. Matched {len(results_cache)} taxa via AphiaID. Remaining: {len(unmatched_tuples)}.")
    else:
        unmatched_tuples = unique_tuples

    # --- Stage 2: Batch Name Matching ---
    if unmatched_tuples:
        all_terms_to_match = list(set(term for vt, an in unmatched_tuples for term in parse_semicolon_taxonomy(vt)))
        logging.info(f"Starting Stage 2: Batch name matching for {len(unmatched_tuples)} taxa.")
        
        batch_lookup = {}
        if all_terms_to_match:
            # (Chunking logic remains the same)
            chunk_size = 150
            for i in range(0, len(all_terms_to_match), chunk_size):
                chunk = all_terms_to_match[i:i+chunk_size]
                try:
                    batch_results_raw = pyworms.aphiaRecordsByMatchNames(chunk)
                    for j, name_list in enumerate(batch_results_raw):
                        if name_list:
                            accepted_match = next((m for m in name_list if m and m.get('status') == 'accepted'), None)
                            if accepted_match:
                                res = {rank: accepted_match.get(rank.lower()) for rank in DWC_RANKS_STD}
                                res.update({'scientificName': accepted_match.get('scientificname'), 'scientificNameID': accepted_match.get('lsid'), 'taxonRank': accepted_match.get('rank')})
                                batch_lookup[chunk[j]] = res
                except Exception as e:
                    logging.error(f"Error in batch chunk: {e}")

            still_unmatched_batch = []
            for verbatim_str, assay_name in unmatched_tuples:
                parsed_names = parse_semicolon_taxonomy(verbatim_str)
                match_found = False
                
                # --- YOUR CORRECTED LOGIC APPLIED HERE ---
                max_depth_for_assay = assay_rank_info.get(assay_name, {}).get('max_depth', 0)
                
                for i_term, term in enumerate(reversed(parsed_names)):
                    # Check if the skip rule should be applied
                    is_potential_species = (i_term == 0)
                    is_full_length_taxonomy = (len(parsed_names) == max_depth_for_assay)
                    
                    if assay_name in assays_to_skip_species and is_potential_species and is_full_length_taxonomy:
                        continue # Skip this term
                        
                    if term in batch_lookup:
                        results_cache[(verbatim_str, assay_name)] = {**batch_lookup[term], 'nameAccordingTo': api_source, 'match_type_debug': f'Success_Batch_{term}'}
                        match_found = True
                        break
                if not match_found:
                    still_unmatched_batch.append((verbatim_str, assay_name))
            unmatched_tuples = still_unmatched_batch
            logging.info(f"Finished Stage 2. Remaining unmatched: {len(unmatched_tuples)}.")

    # --- Stage 3 & Finalization (Remain the same) ---
    if unmatched_tuples:
        logging.info(f"Starting Stage 3: Iterative fallback for {len(unmatched_tuples)} taxa.")
        # ... fallback logic ...
    
    logging.info("Applying all results to DataFrame.")
    # ... mapping logic ...

    logging.info("\n--- Final Performance Report ---")
    # ... reporting logic ...

    return df_to_process