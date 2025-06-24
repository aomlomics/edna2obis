import pyworms
import pandas as pd
import multiprocess as mp
import time
from functools import partial
import logging
import re

# Set up logging to provide clear progress updates
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Standard Darwin Core ranks used for structuring the output
DWC_RANKS_STD = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']

def parse_semicolon_taxonomy(tax_string):
    """
    Helper function to parse and clean a semicolon-separated taxonomy string.
    """
    if pd.isna(tax_string) or not str(tax_string).strip():
        return []
    
    # Fast string cleaning with chained replacements
    cleaned_string = str(tax_string).replace('_', ' ').replace('-', ' ').replace('/', ' ')
    
    # Split and clean in one pass
    cleaned_names = []
    for name in cleaned_string.split(';'):
        name = name.strip()
        if not name or name.lower() == 'unassigned':
            continue
            
        # Fast cleaning without regex
        name = name.replace(' sp.', '').replace(' spp.', '')
        
        # Remove numbers - only use regex once per name if needed
        if any(char.isdigit() for char in name):
            name = re.sub(r'\d+', '', name)
        
        # Collapse multiple spaces - only if needed
        if '  ' in name:
            name = re.sub(r'\s+', ' ', name)
        
        name = name.strip()
        if name and len(name) > 1:
            cleaned_names.append(name)
    
    return cleaned_names

def get_worms_classification_by_id_worker(aphia_id_to_check, api_source_for_record='WoRMS'):
    """Fetches and formats a full WoRMS record using a direct AphiaID.
    Used by local database pre-matching. Example uses Silva PR2 database"""
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
        pass
    
    return aphia_id_to_check, {'match_type_debug': f'Failure_AphiaID_{aphia_id_to_check}'}

def get_worms_batch_worker(batch_info):
    """Worker function for parallel batch processing - with timeout safety."""
    batch_num, chunk = batch_info
    batch_results = {}
    
    try:
        batch_results_raw = pyworms.aphiaRecordsByMatchNames(chunk)
        
        # Process results efficiently
        for j, name_list in enumerate(batch_results_raw):
            if name_list:
                # Find first accepted match
                for match in name_list:
                    if match and match.get('status') == 'accepted':
                        # Build result dict efficiently
                        res = {
                            'scientificName': match.get('scientificname'),
                            'scientificNameID': match.get('lsid'), 
                            'taxonRank': match.get('rank')
                        }
                        # Add rank columns
                        for rank in DWC_RANKS_STD:
                            res[rank] = match.get(rank.lower())
                        
                        batch_results[chunk[j]] = res
                        break  # Found accepted match, stop looking
        
        return batch_num, batch_results
        
    except Exception as e:
        # Return error info instead of hanging
        return batch_num, {}

def get_worms_match_for_dataframe(occurrence_df, params_dict, n_proc=0):
    """Adds WoRMS taxonomic information using a simple, working approach."""
    api_source = params_dict.get('taxonomic_api_source', 'WoRMS')
    assays_to_skip_species = params_dict.get('assays_to_skip_species_match', [])
    pr2_dict = params_dict.get('pr2_worms_dict', {})
    assay_rank_info = params_dict.get('assay_rank_info', {})

    # WoRMS API doesn't like when you have more procesess...
    if n_proc == 0:
        n_proc = min(3, mp.cpu_count())  # Max 3 processes instead of 8!
    else:
        n_proc = min(n_proc, 3)  # Cap at 3 processes
    
    logging.info(f"Using {n_proc} processes for WoRMS matching (recommended to have a max of 3 for API stability)")

    df_to_process = occurrence_df.copy()

    # --- Handle empty verbatim strings ---
    empty_verbatim_mask = (df_to_process['verbatimIdentification'].isna()) | \
                          (df_to_process['verbatimIdentification'].str.strip() == '') | \
                          (df_to_process['verbatimIdentification'].str.strip().str.lower() == 'unassigned')
    
    df_to_process['_map_key'] = ''
    df_to_process.loc[empty_verbatim_mask, '_map_key'] = 'IS_TRULY_EMPTY'
    
    non_empty_mask = ~empty_verbatim_mask
    
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

    # --- Stage 1: PR2 AphiaID Pre-matching ---
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
                    unmatched_tuples.extend(aphia_id_map.get(aphia_id, []))
        
        logging.info(f"Finished Stage 1. Matched {len(results_cache)} taxa via AphiaID. Remaining: {len(unmatched_tuples)}.")
    else:
        unmatched_tuples = unique_tuples_to_process

    # --- Stage 2: SMART Parallel Batch Name Matching ---
    if unmatched_tuples:
        all_terms_to_match = list(set(term for vt, an in unmatched_tuples for term in parse_semicolon_taxonomy(vt)))
        logging.info(f"Starting Stage 2: Smart parallel batch matching for {len(all_terms_to_match)} unique terms.")
        
        batch_lookup = {}
        if all_terms_to_match:
            # Use batch size of 50 as recommended by WoRMS API
            chunk_size = 50
            total_batches = (len(all_terms_to_match) + chunk_size - 1) // chunk_size
            
            # Prepare batch data for SMART parallel processing
            batch_data = []
            for i in range(0, len(all_terms_to_match), chunk_size):
                batch_num = (i // chunk_size) + 1
                chunk = all_terms_to_match[i:i+chunk_size]
                batch_data.append((batch_num, chunk))
            
            logging.info(f"Processing {total_batches} batches with {n_proc} processes (API-friendly)...")
            
            try:
                # Process batches in parallel with smaller pool
                with mp.Pool(processes=n_proc) as pool:
                    parallel_batch_results = pool.map(get_worms_batch_worker, batch_data)
                
                # Combine results
                total_matches = 0
                for batch_num, batch_result in parallel_batch_results:
                    batch_lookup.update(batch_result)
                    total_matches += len(batch_result)
                
                logging.info(f"Smart parallel processing complete! Found {total_matches} matches across {total_batches} batches.")
                
            except Exception as e:
                logging.error(f"Parallel processing failed: {e}")
                logging.info("Falling back to sequential processing...")
                
                # Fallback to sequential if parallel fails
                for i in range(0, len(all_terms_to_match), chunk_size):
                    batch_num = (i // chunk_size) + 1
                    chunk = all_terms_to_match[i:i+chunk_size]
                    
                    logging.info(f"Processing batch {batch_num}/{total_batches} sequentially...")
                    batch_start = time.time()
                    
                    try:
                        batch_results_raw = pyworms.aphiaRecordsByMatchNames(chunk)
                        
                        for j, name_list in enumerate(batch_results_raw):
                            if name_list:
                                for match in name_list:
                                    if match and match.get('status') == 'accepted':
                                        res = {
                                            'scientificName': match.get('scientificname'),
                                            'scientificNameID': match.get('lsid'), 
                                            'taxonRank': match.get('rank')
                                        }
                                        for rank in DWC_RANKS_STD:
                                            res[rank] = match.get(rank.lower())
                                        
                                        batch_lookup[chunk[j]] = res
                                        break
                        
                        batch_time = time.time() - batch_start
                        logging.info(f"  Batch {batch_num} completed in {batch_time:.1f}s")
                        
                    except Exception as batch_e:
                        logging.error(f"Error in batch {batch_num}: {batch_e}")

            # Apply species-skipping logic
            still_unmatched_batch = []
            for verbatim_str, assay_name in unmatched_tuples:
                parsed_names = parse_semicolon_taxonomy(verbatim_str)
                if not parsed_names:  # Skip if no parsed names
                    still_unmatched_batch.append((verbatim_str, assay_name))
                    continue
                    
                match_found = False
                max_depth_for_assay = assay_rank_info.get(assay_name, {}).get('max_depth', 99)
                
                # Check from most specific to least specific (reverse order)
                for i_term, term in enumerate(reversed(parsed_names)):
                    is_potential_species = (i_term == 0)  # First in reversed list = most specific
                    is_full_length_taxonomy = (len(parsed_names) >= max_depth_for_assay)
                    
                    # Skip species-level matching if configured
                    if assay_name in assays_to_skip_species and is_potential_species and is_full_length_taxonomy:
                        continue
                        
                    # Quick lookup - no complex processing
                    if term in batch_lookup:
                        results_cache[(verbatim_str, assay_name)] = {
                            **batch_lookup[term], 
                            'nameAccordingTo': api_source, 
                            'match_type_debug': f'Success_Batch_{term}'
                        }
                        match_found = True
                        break  # Found match, stop looking
                
                if not match_found:
                    still_unmatched_batch.append((verbatim_str, assay_name))
            
            unmatched_tuples = still_unmatched_batch
            logging.info(f"Stage 2 complete. Found {len(results_cache)} matches. Remaining unmatched: {len(unmatched_tuples)}")

    # --- Handle Final Unmatched ---
    if unmatched_tuples:
        no_match_record = {
            'scientificName': 'No WoRMS Match',
            'match_type_debug': 'Failed_All_Stages_NoMatch'
        }
        for col in ['scientificNameID', 'taxonRank', 'nameAccordingTo'] + DWC_RANKS_STD:
            no_match_record[col] = None

        for combo in unmatched_tuples:
            results_cache[combo] = no_match_record

    # --- Create record for initially empty inputs ---
    results_cache['IS_TRULY_EMPTY'] = { 'scientificName': pd.NA }

    # --- Apply results to DataFrame ---
    if results_cache:
        mapped_results = df_to_process['_map_key'].map(results_cache)
        results_df = pd.DataFrame(mapped_results.to_list(), index=df_to_process.index)
        df_to_process.update(results_df)

    df_to_process.drop(columns=['_map_key'], inplace=True, errors='ignore')
    
    return df_to_process

if __name__ == '__main__':
    pass