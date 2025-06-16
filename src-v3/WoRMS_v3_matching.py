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
    try:
        record = pyworms.aphiaRecordByAphiaID(aphia_id_to_check)
        
        if record and isinstance(record, dict) and record.get('status') == 'accepted':
            result = {
                'scientificName': record.get('scientificname'), 'scientificNameID': record.get('lsid'),
                'taxonRank': record.get('rank'), 'nameAccordingTo': api_source_for_record,
                'match_type_debug': f'Success_AphiaID_{aphia_id_to_check}'
            }
            for rank_std in DWC_RANKS_STD:
                result[rank_std] = record.get(rank_std.lower())
            return aphia_id_to_check, result
    except Exception:
        # Log quietly, as some failures are expected
        pass
    
    return aphia_id_to_check, {'match_type_debug': f'Failure_AphiaID_{aphia_id_to_check}'}


def get_worms_match_for_dataframe(occurrence_df, params_dict, n_proc=0):
    """Adds WoRMS taxonomic information using an optimized multi-stage process."""
    overall_start_time = time.time()
    api_source = params_dict.get('taxonomic_api_source', 'WoRMS')
    assays_to_skip_species = params_dict.get('assays_to_skip_species_match', [])
    pr2_dict = params_dict.get('pr2_worms_dict', {})
    assay_rank_info = params_dict.get('assay_rank_info', {})

    if n_proc == 0:
        n_proc = mp.cpu_count()

    df_to_process = occurrence_df.copy()

    # --- Triage Step: Identify truly empty/unassigned verbatim strings from the start ---
    empty_verbatim_mask = (df_to_process['verbatimIdentification'].isna()) | \
                          (df_to_process['verbatimIdentification'].str.strip() == '') | \
                          (df_to_process['verbatimIdentification'].str.strip().str.lower() == 'unassigned')
    
    df_to_process['_map_key'] = ''
    df_to_process.loc[empty_verbatim_mask, '_map_key'] = 'IS_TRULY_EMPTY'
    
    non_empty_mask = ~empty_verbatim_mask
    
    # This is the fix for the ValueError. We create a Series of tuples with the
    # correct index before assigning, which is a safer operation in pandas.
    keys_to_assign = pd.Series(
        list(zip(
            df_to_process.loc[non_empty_mask, 'verbatimIdentification'],
            df_to_process.loc[non_empty_mask, 'assay_name']
        )),
        index=df_to_process[non_empty_mask].index
    )
    df_to_process.loc[non_empty_mask, '_map_key'] = keys_to_assign
    
    unique_tuples_to_process = list(df_to_process.loc[non_empty_mask, '_map_key'].drop_duplicates())
    logging.info(f"Found {len(unique_tuples_to_process)} unique, non-empty combinations to process.")

    results_cache = {}
    unmatched_tuples = []

    # --- Stage 1: Parallel AphiaID Pre-matching ---
    if pr2_dict:
        logging.info("Starting Stage 1: Parallel AphiaID pre-matching.")
        aphia_id_map = {}
        
        for verbatim_str, assay_name in unique_tuples_to_process:
            parsed_names = parse_semicolon_taxonomy(verbatim_str)
            if parsed_names:
                species_name = parsed_names[-1]
                if species_name in pr2_dict:
                    aphia_id = pr2_dict[species_name]
                    if aphia_id not in aphia_id_map:
                        aphia_id_map[aphia_id] = []
                    aphia_id_map[aphia_id].append((verbatim_str, assay_name))
                    continue
            unmatched_tuples.append((verbatim_str, assay_name))
        
        unique_aphia_ids_to_fetch = list(aphia_id_map.keys())
        if unique_aphia_ids_to_fetch:
            with mp.Pool(processes=n_proc) as pool:
                worker_func = partial(get_worms_classification_by_id_worker, api_source_for_record=api_source)
                parallel_results = pool.map(worker_func, unique_aphia_ids_to_fetch)

            for aphia_id, result in parallel_results:
                if 'scientificName' in result:
                    for combo in aphia_id_map.get(aphia_id, []):
                        results_cache[combo] = result
                else:
                    # If AphiaID lookup fails, add the original tuples to the name-based lookup list
                    unmatched_tuples.extend(aphia_id_map.get(aphia_id, []))
        
        logging.info(f"Finished Stage 1. Matched {len(results_cache)} taxa via AphiaID. Remaining: {len(unmatched_tuples)}.")
    else:
        unmatched_tuples = unique_tuples_to_process

    # --- Stage 2: Batch Name Matching ---
    if unmatched_tuples:
        all_terms_to_match = list(set(term for vt, an in unmatched_tuples for term in parse_semicolon_taxonomy(vt)))
        logging.info(f"Starting Stage 2: Batch name matching for {len(all_terms_to_match)} unique terms.")
        
        batch_lookup = {}
        if all_terms_to_match:
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
                max_depth_for_assay = assay_rank_info.get(assay_name, {}).get('max_depth', 99)
                
                for i_term, term in enumerate(reversed(parsed_names)):
                    is_potential_species = (i_term == 0)
                    is_full_length_taxonomy = (len(parsed_names) >= max_depth_for_assay)
                    
                    if assay_name in assays_to_skip_species and is_potential_species and is_full_length_taxonomy:
                        continue
                        
                    if term in batch_lookup:
                        results_cache[(verbatim_str, assay_name)] = {**batch_lookup[term], 'nameAccordingTo': api_source, 'match_type_debug': f'Success_Batch_{term}'}
                        match_found = True
                        break
                if not match_found:
                    still_unmatched_batch.append((verbatim_str, assay_name))
            unmatched_tuples = still_unmatched_batch
            logging.info(f"Finished Stage 2. Remaining unmatched: {len(unmatched_tuples)}.")

    # --- Stage 3: Handle Final Unmatched ---
    if unmatched_tuples:
        logging.info(f"Tagging {len(unmatched_tuples)} remaining lookup failures as 'No WoRMS Match'.")
        no_match_record = {
            'scientificName': 'No WoRMS Match',
            'match_type_debug': 'Failed_All_Stages_NoMatch'
        }
        # Ensure all other columns are NaN so they don't carry over old data
        for col in ['scientificNameID', 'taxonRank', 'nameAccordingTo'] + DWC_RANKS_STD:
            no_match_record[col] = None

        for combo in unmatched_tuples:
            results_cache[combo] = no_match_record

    # --- Create record for initially empty inputs ---
    results_cache['IS_TRULY_EMPTY'] = { 'scientificName': pd.NA }

    # --- APPLY ALL RESULTS TO THE DATAFRAME ---
    logging.info("Applying all results to DataFrame.")
    if results_cache:
        taxonomic_cols = ['scientificName', 'scientificNameID', 'taxonRank', 'nameAccordingTo', 'match_type_debug'] + DWC_RANKS_STD
        
        # Map the results using the unique key
        mapped_results = df_to_process['_map_key'].map(results_cache)
        
        # Create a DataFrame from the mapped results and assign it
        results_df = pd.DataFrame(mapped_results.to_list(), index=df_to_process.index)
        
        # Update the main dataframe with the new results
        df_to_process.update(results_df)

    df_to_process.drop(columns=['_map_key'], inplace=True, errors='ignore')

    overall_end_time = time.time()
    logging.info(f"\n--- WoRMS Matching Complete ---")
    logging.info(f"Total processing time: {overall_end_time - overall_start_time:.2f} seconds.")
    
    return df_to_process

if __name__ == '__main__':
    pass