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
    Helper function to parse a semicolon-separated taxonomy string.
    It cleans the string and splits it into a list of taxonomic terms.
    """
    if pd.isna(tax_string) or not str(tax_string).strip():
        return []
    # Do not replace underscores here, preserve original naming from database
    return [name.strip() for name in str(tax_string).split(';') if name.strip() and name.strip().lower() not in ['unassigned']]

def get_cleaned_species_name(species_term):
    """Cleans a single species-level term for dictionary lookup."""
    if not species_term or pd.isna(species_term):
        return None
    # Remove prefix and replace underscores for lookup
    name = str(species_term)
    if name.lower().startswith('s__'):
        name = name[3:]
    return name.replace('_', ' ').strip()

def get_worms_classification_by_id_worker(aphia_id_to_check, api_source_for_record='WoRMS'):
    """Fetches and formats a full WoRMS record using a direct AphiaID."""
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
        # Correctly check if the term is a species-level assignment
        is_species_level = str(term).lower().startswith('s__')
        
        if is_species_level and assay_name_str in assays_to_skip_species_list:
            continue
        
        # Clean the term for searching (remove prefix)
        clean_name_to_match = str(term)
        if clean_name_to_match[1:3] == '__':
            clean_name_to_match = clean_name_to_match[3:]
        
        try:
            time.sleep(0.001)
            start_time = time.time()
            s_match_list = pyworms.aphiaRecordsByName(clean_name_to_match, like=False, marine_only=True)
            end_time = time.time()
            api_call_count += 1
            total_api_time += (end_time - start_time)
            
            if s_match_list:
                accepted_match = next((m for m in s_match_list if m and m.get('status') == 'accepted'), None)
                if accepted_match:
                    result_dict = {'scientificName': accepted_match.get('scientificname'), 'scientificNameID': accepted_match.get('lsid'), 'taxonRank': accepted_match.get('rank'), 'nameAccordingTo': api_source_for_record, 'match_type_debug': f'Success_Fallback_{clean_name_to_match}', 'api_call_count': api_call_count, 'total_api_time': total_api_time}
                    for rank_std in DWC_RANKS_STD:
                        result_dict[rank_std] = accepted_match.get(rank_std.lower())
                    return result_dict
        except Exception as e:
            logging.warning(f"Fallback API call error for '{clean_name_to_match}': {e}")
    
    return {'match_type_debug': 'Failed_Fallback_IncertaeSedis', 'api_call_count': api_call_count, 'total_api_time': total_api_time}

def get_worms_match_for_dataframe(occurrence_df, params_dict, n_proc=0):
    """Adds WoRMS taxonomic information using an optimized multi-stage process."""
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
            if parsed_names:
                species_term = parsed_names[-1]
                if str(species_term).lower().startswith('s__'):
                    cleaned_species_name = get_cleaned_species_name(species_term)
                    if cleaned_species_name in pr2_dict:
                        aphia_id = pr2_dict[cleaned_species_name]
                        result = get_worms_classification_by_id_worker(aphia_id, api_source)
                        aphia_api_calls += result.get('api_call_count', 0)
                        aphia_api_time += result.get('total_api_time', 0)
                        if 'scientificName' in result:
                            results_cache[(verbatim_str, assay_name)] = result
                            continue
            unmatched_tuples.append((verbatim_str, assay_name))
        logging.info(f"Finished Stage 1. Matched {len(results_cache)} taxa via AphiaID. Remaining: {len(unmatched_tuples)}.")
    else:
        logging.info("Skipping Stage 1: No PR2 AphiaID dictionary provided.")
        unmatched_tuples = unique_tuples

    # --- Stage 2: Batch Name Matching ---
    if unmatched_tuples:
        all_terms_to_match = list(set(term for vt, an in unmatched_tuples for term in parse_semicolon_taxonomy(vt)))
        logging.info(f"Starting Stage 2: Batch name matching for {len(unmatched_tuples)} taxa, using {len(all_terms_to_match)} unique terms.")
        
        batch_lookup = {}
        if all_terms_to_match:
            chunk_size = 150
            for i in range(0, len(all_terms_to_match), chunk_size):
                chunk = all_terms_to_match[i:i+chunk_size]
                try:
                    logging.info(f"Sending batch {i//chunk_size + 1}, names {i+1}-{i+len(chunk)}/{len(all_terms_to_match)}.")
                    start_time = time.time()
                    batch_results_raw = pyworms.aphiaRecordsByMatchNames(chunk)
                    end_time = time.time()
                    batch_api_calls += 1
                    batch_api_time += (end_time - start_time)
                    
                    for j, name_list in enumerate(batch_results_raw):
                        if name_list:
                            accepted_match = next((m for m in name_list if m and m.get('status') == 'accepted'), None)
                            if accepted_match:
                                res = {rank: accepted_match.get(rank.lower()) for rank in DWC_RANKS_STD}
                                res.update({'scientificName': accepted_match.get('scientificname'), 'scientificNameID': accepted_match.get('lsid'), 'taxonRank': accepted_match.get('rank')})
                                batch_lookup[chunk[j]] = res
                except Exception as e:
                    logging.error(f"Error in batch chunk {i//chunk_size + 1}: {e}")

            still_unmatched_batch = []
            for verbatim_str, assay_name in unmatched_tuples:
                parsed_names = parse_semicolon_taxonomy(verbatim_str)
                match_found = False
                for term in reversed(parsed_names):
                    is_species_level = str(term).lower().startswith('s__')
                    if is_species_level and assay_name in assays_to_skip_species:
                        continue
                    if term in batch_lookup:
                        results_cache[(verbatim_str, assay_name)] = {**batch_lookup[term], 'nameAccordingTo': api_source, 'match_type_debug': f'Success_Batch_{term}'}
                        match_found = True
                        break
                if not match_found:
                    still_unmatched_batch.append((verbatim_str, assay_name))
            unmatched_tuples = still_unmatched_batch
            logging.info(f"Finished Stage 2. Remaining unmatched: {len(unmatched_tuples)}.")

    # --- Stage 3: Iterative Fallback ---
    if unmatched_tuples:
        logging.info(f"Starting Stage 3: Iterative fallback for {len(unmatched_tuples)} taxa.")
        num_proc = min(mp.cpu_count() if n_proc == 0 else int(n_proc), len(unmatched_tuples))
        if num_proc > 0:
            with mp.Pool(processes=num_proc) as pool:
                worker_func = partial(get_worms_match_for_single_taxon_worker, assays_to_skip_species_list=assays_to_skip_species, api_source_for_record=api_source)
                fallback_results = pool.map(worker_func, unmatched_tuples)

            for i, combo_key in enumerate(unmatched_tuples):
                res = fallback_results[i]
                results_cache[combo_key] = res
                fallback_api_calls += res.get('api_call_count', 0)
                fallback_api_time += res.get('total_api_time', 0)
        logging.info("Finished Stage 3.")

    # --- Map all results and Finalize ---
    logging.info("Applying all results to DataFrame.")
    result_cols = ['scientificName', 'scientificNameID', 'taxonRank', 'nameAccordingTo', 'match_type_debug'] + DWC_RANKS_STD
    
    for col in result_cols:
        df_to_process[col] = df_to_process['_map_key'].map({k: v.get(col) for k, v in results_cache.items()})

    incertae_sedis_info = {'scientificName': 'incertae sedis', 'scientificNameID': 'urn:lsid:marinespecies.org:taxname:12', 'taxonRank': 'no rank'}
    for col, val in incertae_sedis_info.items():
        df_to_process[col].fillna(val, inplace=True)
        
    df_to_process.drop(columns=['_map_key'], inplace=True)
    
    logging.info("\n--- Final Performance Report ---")
    logging.info(f"Total unique taxa: {len(unique_tuples)}")
    if pr2_dict:
        matched_by_aphiaid = sum(1 for v in results_cache.values() if v.get('match_type_debug', '').startswith('Success_AphiaID'))
        logging.info(f"  - Matched by AphiaID: {matched_by_aphiaid}")
    logging.info(f"  - Total API calls: {aphia_api_calls + batch_api_calls + fallback_api_calls} (AphiaID: {aphia_api_calls}, Batch: {batch_api_calls}, Fallback: {fallback_api_calls})")
    logging.info(f"  - Total API time: {aphia_api_time + batch_api_time + fallback_api_time:.2f}s")
    logging.info(f"Total execution time: {time.time() - overall_start_time:.2f}s")
    logging.info("--------------------------\n")

    return df_to_process