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
                'match_type_debug': f'Success_AphiaID_{aphia_id_to_check}'
            }
            for rank_std in DWC_RANKS_STD:
                result[rank_std] = record.get(rank_std.lower())
            return result
    except Exception as e:
        logging.warning(f"API call failed for AphiaID {aphia_id_to_check}: {e}")
    
    return {'match_type_debug': f'Failure_AphiaID_{aphia_id_to_check}'}


def get_worms_match_for_single_taxon_worker(combo_input_for_worker, assays_to_skip_species_list, api_source_for_record='WoRMS'):
    """(FALLBACK) Matches a single verbatimIdentification string using a slow, iterative search."""
    verbatim_identification_str, assay_name_str = combo_input_for_worker

    parsed_names = parse_semicolon_taxonomy(verbatim_identification_str)
    if not parsed_names:
        return {'match_type_debug': 'No_Tax_String_Provided'}

    for i, term in enumerate(reversed(parsed_names)):
        if i == 0 and assay_name_str in assays_to_skip_species_list:
            continue
        
        try:
            time.sleep(0.001)
            s_match_list = pyworms.aphiaRecordsByName(term, like=False, marine_only=True)
            
            if s_match_list:
                accepted_match = next((m for m in s_match_list if m and m.get('status') == 'accepted'), None)
                if accepted_match:
                    result_dict = {'scientificName': accepted_match.get('scientificname'), 'scientificNameID': accepted_match.get('lsid'), 'taxonRank': accepted_match.get('rank'), 'nameAccordingTo': api_source_for_record, 'match_type_debug': f'Success_Fallback_{term}'}
                    for rank_std in DWC_RANKS_STD:
                        result_dict[rank_std] = accepted_match.get(rank_std.lower())
                    return result_dict
        except Exception as e:
            logging.warning(f"Fallback API call error for '{term}': {e}")
    
    return {'match_type_debug': 'Failed_Fallback_IncertaeSedis'}


def get_worms_match_for_dataframe(occurrence_df, params_dict, n_proc=0):
    """Adds WoRMS taxonomic information using an optimized multi-stage process."""
    overall_start_time = time.time()
    api_source = params_dict.get('taxonomic_api_source', 'WoRMS')
    assays_to_skip_species = params_dict.get('assays_to_skip_species_match', [])
    pr2_dict = params_dict.get('pr2_worms_dict', {})
    assay_rank_info = params_dict.get('assay_rank_info', {})

    df_to_process = occurrence_df.copy()
    df_to_process['verbatimIdentification'] = df_to_process['verbatimIdentification'].astype(str).fillna('')
    df_to_process['_map_key'] = list(zip(df_to_process['verbatimIdentification'], df_to_process['assay_name']))
    unique_tuples = list(df_to_process['_map_key'].drop_duplicates())
    unique_tuples = [combo for combo in unique_tuples if str(combo[0]).strip()]
    logging.info(f"Found {len(unique_tuples)} unique combinations to process.")

    results_cache = {}
    unmatched_tuples = []

    # --- Stage 1: AphiaID Pre-matching ---
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
                
                max_depth_for_assay = assay_rank_info.get(assay_name, {}).get('max_depth', 99) # Default to a high number if not found
                
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

    # --- Stage 3: Iterative Fallback for Remaining Taxa ---
    if unmatched_tuples:
        logging.info(f"Starting Stage 3: Iterative fallback for {len(unmatched_tuples)} taxa.")
        # This is a slow, single-threaded process for the remainders
        for combo in unmatched_tuples:
            result = get_worms_match_for_single_taxon_worker(combo, assays_to_skip_species, api_source)
            results_cache[combo] = result

    # --- FIX: APPLY THE RESULTS TO THE DATAFRAME ---
    logging.info("Applying all results to DataFrame.")
    if results_cache:
        # Define the columns that will be populated from the cache
        taxonomic_cols = [
            'scientificName', 'scientificNameID', 'taxonRank', 'nameAccordingTo', 
            'match_type_debug'
        ] + DWC_RANKS_STD

        # Create temporary columns by mapping the cache to the unique key
        for col in taxonomic_cols:
            df_to_process[f'_temp_{col}'] = df_to_process['_map_key'].map(
                {k: v.get(col) for k, v in results_cache.items()}
            )

        # Update the original columns from the temporary ones
        for col in taxonomic_cols:
            temp_col = f'_temp_{col}'
            # Only update rows where a new value was actually found
            update_mask = df_to_process[temp_col].notna()
            df_to_process.loc[update_mask, col] = df_to_process.loc[update_mask, temp_col]
        
        # Drop all temporary columns
        df_to_process.drop(columns=[f'_temp_{col}' for col in taxonomic_cols], inplace=True)

    # Clean up the helper key before returning
    df_to_process.drop(columns=['_map_key'], inplace=True, errors='ignore')

    overall_end_time = time.time()
    logging.info(f"\n--- WoRMS Matching Complete ---")
    logging.info(f"Total processing time: {overall_end_time - overall_start_time:.2f} seconds.")
    
    return df_to_process

if __name__ == '__main__':
    pass