import pyworms
import pandas as pd
import multiprocess as mp
import time
from functools import partial
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Standard Darwin Core ranks
DWC_RANKS_STD = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']

def parse_semicolon_taxonomy(tax_string):
    """Helper function to parse a semicolon-separated taxonomy string into a list of names."""
    if pd.isna(tax_string) or not str(tax_string).strip():
        return []
    return [name.strip() for name in str(tax_string).split(';') if name.strip() and name.strip().lower() not in ['unassigned']]

def get_worms_classification_by_id_worker(aphia_id_to_check, api_source_for_record='WoRMS'):
    """Fetches and formats a full WoRMS record using a direct AphiaID. Includes profiling."""
    start_time = time.time()
    try:
        record = pyworms.aphiaRecordByAphiaID(aphia_id_to_check)
        end_time = time.time()
        
        if record and isinstance(record, dict) and record.get('status') == 'accepted':
            result = {
                'scientificName': record.get('scientificname'),
                'scientificNameID': record.get('lsid'),
                'taxonRank': record.get('rank'),
                'nameAccordingTo': api_source_for_record,
                'match_type_debug': f'Success_AphiaID_{aphia_id_to_check}',
                'api_call_count': 1,
                'total_api_time': end_time - start_time
            }
            for rank_std in DWC_RANKS_STD:
                result[rank_std] = record.get(rank_std.lower())
            return result
    except Exception as e:
        logging.warning(f"API call failed for AphiaID {aphia_id_to_check}: {e}")
    
    # Return failure case
    return {'api_call_count': 1, 'total_api_time': time.time() - start_time, 'match_type_debug': f'Failure_AphiaID_{aphia_id_to_check}'}


def get_worms_match_for_single_taxon_worker(combo_input_for_worker, assays_to_skip_species_list, api_source_for_record='WoRMS'):
    """
    (FALLBACK) Core logic to match a single verbatimIdentification string against WoRMS using iterative search.
    """
    verbatim_identification_str, assay_name_str = combo_input_for_worker
    api_call_count = 0
    total_api_time = 0.0

    parsed_names = parse_semicolon_taxonomy(verbatim_identification_str)
    if not parsed_names:
        return {'match_type_debug': 'No_Tax_String_Provided', 'api_call_count': 0, 'total_api_time': 0.0}

    dwc_rank_hierarchy_for_iteration = ['species', 'genus', 'family', 'order', 'class', 'phylum', 'kingdom']
    
    for i in range(len(parsed_names)):
        name_to_match = parsed_names[-(i+1)]
        effective_rank_being_tried = dwc_rank_hierarchy_for_iteration[i] if i < len(dwc_rank_hierarchy_for_iteration) else f"level_{i}"

        if not name_to_match or name_to_match.lower() in ['unassigned', 'incertae sedis']:
            continue
            
        if assay_name_str in assays_to_skip_species_list and i == 0:
            continue
        
        try:
            time.sleep(0.001)
            start_time = time.time()
            s_match_list = pyworms.aphiaRecordsByName(name_to_match, like=False, marine_only=True)
            end_time = time.time()
            api_call_count += 1
            total_api_time += (end_time - start_time)
            
            chosen_match_record = None
            if s_match_list:
                accepted_matches = [m for m in s_match_list if m and m.get('status') == 'accepted']
                if accepted_matches:
                    chosen_match_record = accepted_matches[0] # Take the first accepted match

            if chosen_match_record:
                result_dict = {
                    'scientificName': chosen_match_record.get('scientificname'),
                    'scientificNameID': chosen_match_record.get('lsid'),
                    'taxonRank': chosen_match_record.get('rank'), 
                    'nameAccordingTo': api_source_for_record,
                    'match_type_debug': f'Success_Fallback_{name_to_match}',
                    'api_call_count': api_call_count, 'total_api_time': total_api_time
                }
                for rank_std in DWC_RANKS_STD:
                    result_dict[rank_std] = chosen_match_record.get(rank_std.lower())
                return result_dict
        
        except Exception as e:
            logging.warning(f"Fallback API call error for '{name_to_match}': {e}")
    
    return {'match_type_debug': 'Failed_Fallback_IncertaeSedis', 'api_call_count': api_call_count, 'total_api_time': total_api_time}


def get_worms_match_for_dataframe(occurrence_df, params_dict, n_proc=0):
    if not isinstance(occurrence_df, pd.DataFrame) or occurrence_df.empty:
        logging.warning("Input DataFrame is empty. Skipping WoRMS matching.")
        return occurrence_df.copy()

    overall_start_time = time.time()
    api_source = params_dict.get('taxonomic_api_source', 'WoRMS')
    assays_to_skip_species = params_dict.get('assays_to_skip_species_match', [])
    pr2_dict = params_dict.get('pr2_worms_dict', {})
    
    df_to_process = occurrence_df.copy()
    df_to_process['verbatimIdentification'] = df_to_process['verbatimIdentification'].astype(str).fillna('')
    df_to_process['_map_key'] = list(zip(df_to_process['verbatimIdentification'], df_to_process['assay_name']))
    
    unique_tuples = list(df_to_process['_map_key'].drop_duplicates())
    unique_tuples = [combo for combo in unique_tuples if str(combo[0]).strip()]
    logging.info(f"Found {len(unique_tuples)} unique combinations to process.")

    results_cache = {}
    unmatched_tuples = []
    aphia_api_calls, aphia_api_time = 0, 0.0
    batch_api_calls, batch_api_time = 0, 0.0
    fallback_api_calls, fallback_api_time = 0, 0.0

    # --- Stage 1: AphiaID Pre-matching ---
    if pr2_dict:
        logging.info("Starting Stage 1: AphiaID pre-matching.")
        for verbatim_str, assay_name in unique_tuples:
            parsed_names = parse_semicolon_taxonomy(verbatim_str)
            if parsed_names and parsed_names[-1] in pr2_dict:
                aphia_id = pr2_dict[parsed_names[-1]]
                result = get_worms_classification_by_id_worker(aphia_id, api_source)
                aphia_api_calls += result.get('api_call_count', 0)
                aphia_api_time += result.get('total_api_time', 0)
                if 'scientificName' in result:
                    results_cache[(verbatim_str, assay_name)] = result
                    continue
            unmatched_tuples.append((verbatim_str, assay_name))
        logging.info(f"Finished Stage 1. Matched {len(results_cache)} taxa. Remaining: {len(unmatched_tuples)}.")
    else:
        unmatched_tuples = unique_tuples

    # --- Stage 2: Batch Name Matching ---
    if unmatched_tuples:
        logging.info(f"Starting Stage 2: Batch name matching for {len(unmatched_tuples)} taxa.")
        all_terms_to_match = list(set(term for vt, an in unmatched_tuples for term in parse_semicolon_taxonomy(vt)))
        
        if all_terms_to_match:
            try:
                start_time = time.time()
                batch_results_raw = pyworms.aphiaRecordsByMatchNames(all_terms_to_match)
                end_time = time.time()
                batch_api_calls += 1
                batch_api_time += (end_time - start_time)
                
                batch_lookup = {}
                for i, name_list in enumerate(batch_results_raw):
                    if name_list:
                        accepted_match = next((m for m in name_list if m and m.get('status') == 'accepted'), None)
                        if accepted_match:
                            res = {rank: accepted_match.get(rank.lower()) for rank in DWC_RANKS_STD}
                            res.update({'scientificName': accepted_match.get('scientificname'), 'scientificNameID': accepted_match.get('lsid'), 'taxonRank': accepted_match.get('rank')})
                            batch_lookup[all_terms_to_match[i]] = res
                
                still_unmatched_batch = []
                for verbatim_str, assay_name in unmatched_tuples:
                    parsed_names = parse_semicolon_taxonomy(verbatim_str)
                    match_found = False
                    for name in reversed(parsed_names):
                        if name in batch_lookup:
                            final_res = batch_lookup[name].copy()
                            final_res.update({'nameAccordingTo': api_source, 'match_type_debug': f'Success_Batch_{name}'})
                            results_cache[(verbatim_str, assay_name)] = final_res
                            match_found = True
                            break
                    if not match_found:
                        still_unmatched_batch.append((verbatim_str, assay_name))
                unmatched_tuples = still_unmatched_batch
                logging.info(f"Finished Stage 2. Matched {len(results_cache) - len(unmatched_tuples)} in this stage. Remaining: {len(unmatched_tuples)}.")
            except Exception as e:
                logging.error(f"Error during Stage 2: {e}. Moving all to fallback.")

    # --- Stage 3: Iterative Fallback ---
    if unmatched_tuples:
        logging.info(f"Starting Stage 3: Iterative fallback for {len(unmatched_tuples)} remaining taxa.")
        num_proc = min(mp.cpu_count() if n_proc == 0 else int(n_proc), len(unmatched_tuples))
        if num_proc > 0:
            with mp.Pool(processes=num_proc) as pool:
                tasks = [(combo, assays_to_skip_species, api_source) for combo in unmatched_tuples]
                fallback_results = pool.map(partial(get_worms_match_for_single_taxon_worker, assays_to_skip_species_list=assays_to_skip_species, api_source_for_record=api_source), [t[0] for t in tasks])
            
            for i, combo_key in enumerate(unmatched_tuples):
                res = fallback_results[i]
                results_cache[combo_key] = res
                fallback_api_calls += res.get('api_call_count', 0)
                fallback_api_time += res.get('total_api_time', 0)
        logging.info("Finished Stage 3.")

    # --- Map all results and Finalize ---
    logging.info("Applying all results to DataFrame.")
    map_start = time.time()
    result_cols = ['scientificName', 'scientificNameID', 'taxonRank', 'nameAccordingTo', 'match_type_debug'] + DWC_RANKS_STD
    
    for col in result_cols:
        df_to_process[col] = df_to_process['_map_key'].map({k: v.get(col) for k, v in results_cache.items()})

    # Fill any remaining NaNs with incertae sedis
    incertae_sedis_info = {'scientificName': 'incertae sedis', 'scientificNameID': 'urn:lsid:marinespecies.org:taxname:12', 'taxonRank': 'no rank'}
    for col, val in incertae_sedis_info.items():
        df_to_process[col].fillna(val, inplace=True)
        
    df_to_process.drop(columns=['_map_key'], inplace=True)
    map_end = time.time()

    logging.info("\n--- Final Performance Report ---")
    logging.info(f"Total unique taxa: {len(unique_tuples)}")
    logging.info(f"  - Matched by AphiaID: {len(unique_tuples) - len(unmatched_tuples) if pr2_dict else 'N/A'}")
    logging.info(f"  - Total API calls: {aphia_api_calls + batch_api_calls + fallback_api_calls} (AphiaID: {aphia_api_calls}, Batch: {batch_api_calls}, Fallback: {fallback_api_calls})")
    logging.info(f"  - Total API time: {aphia_api_time + batch_api_time + fallback_api_time:.2f}s (AphiaID: {aphia_api_time:.2f}s, Batch: {batch_api_time:.2f}s, Fallback: {fallback_api_time:.2f}s)")
    logging.info(f"Time to map results: {map_end - map_start:.2f}s")
    logging.info(f"Total execution time: {time.time() - overall_start_time:.2f}s")
    logging.info("--------------------------\n")

    return df_to_process

if __name__ == '__main__':
    # This block is for testing the script directly, not used when imported as a module.
    print("WoRMS_OBIS_matcher.py script is being run directly (e.g., for testing).")
