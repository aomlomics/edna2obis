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
    """
    Worker function for parallel batch processing - with timeout safety.
    Now collects ALL accepted matches to handle ambiguous cases.
    """
    batch_num, chunk = batch_info
    batch_results = {}
    
    try:
        batch_results_raw = pyworms.aphiaRecordsByMatchNames(chunk)
        
        # Process results efficiently
        for j, name_list in enumerate(batch_results_raw):
            if name_list:
                accepted_matches = []
                # Find ALL accepted matches
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
                        
                        accepted_matches.append(res)
                
                if accepted_matches:
                    batch_results[chunk[j]] = accepted_matches  # Store list of all accepted matches
        
        return batch_num, batch_results
        
    except Exception as e:
        return batch_num, {}

def get_worms_match_for_dataframe(occurrence_df, params_dict, n_proc=0):
    """
    Adds WoRMS taxonomic information and generates a detailed info file for ambiguous matches.
    
    Returns:
        dict: A dictionary containing two DataFrames:
              - 'main_df': The occurrence DataFrame with a single best match (maintains original behavior).
              - 'info_df': A detailed DataFrame with all possible matches for ambiguous cases.
    """
    api_source = params_dict.get('taxonomic_api_source', 'WoRMS')
    assays_to_skip_species = params_dict.get('assays_to_skip_species_match', [])
    
    # Handle case where parameter is None
    if assays_to_skip_species is None:
        assays_to_skip_species = []
        
    pr2_dict = params_dict.get('pr2_worms_dict', {})
    assay_rank_info = params_dict.get('assay_rank_info', {})

    if n_proc == 0:
        n_proc = min(3, mp.cpu_count())
    else:
        n_proc = min(n_proc, 3)
    
    logging.info(f"Using {n_proc} processes for WoRMS matching (max 3 recommended)")

    df_to_process = occurrence_df.copy()

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
    info_records = []  # NEW: To build the detailed info DataFrame
    unmatched_tuples = []
    
    cases_to_handle = []
    for verbatim_str, assay_name in unique_tuples_to_process:
        cleaned_verbatim = str(verbatim_str).strip().rstrip(';').strip()
        
        is_unassigned = (not cleaned_verbatim or 
                         cleaned_verbatim.lower() in ['unassigned', 'nan', 'none', ''] or
                         pd.isna(verbatim_str))
        
        is_simple_kingdom = cleaned_verbatim.lower() in ['eukaryota']

        if is_unassigned or is_simple_kingdom:
            match_type = 'incertae_sedis_unassigned' if is_unassigned else f'incertae_sedis_simple_case_{cleaned_verbatim}'
            incertae_sedis_record = {
                'scientificName': 'incertae sedis',
                'scientificNameID': 'urn:lsid:marinespecies.org:taxname:12', 'taxonRank': None,
                'nameAccordingTo': api_source, 'match_type_debug': match_type,
                'cleanedTaxonomy': cleaned_verbatim
            }
            for col in DWC_RANKS_STD:
                incertae_sedis_record[col] = None
            
            results_cache[(verbatim_str, assay_name)] = incertae_sedis_record
            info_records.append({'verbatimIdentification': verbatim_str, **incertae_sedis_record, 'ambiguous': False})
            cases_to_handle.append((verbatim_str, assay_name))
    
    unique_tuples_to_process = [t for t in unique_tuples_to_process if t not in cases_to_handle]
    if cases_to_handle:
        logging.info(f"Assigned {len(cases_to_handle)} cases (unassigned/empty/simple) to 'incertae sedis'.")

    # --- Stage 1: PR2 AphiaID Pre-matching ---
    if pr2_dict:
        logging.info("Starting Stage 1: Parallel AphiaID pre-matching.")
        aphia_id_map = {}
        stage1_unmatched = []
        
        for verbatim_str, assay_name in unique_tuples_to_process:
            if assay_name in assays_to_skip_species:
                stage1_unmatched.append((verbatim_str, assay_name))
                continue
            
            parsed_names = parse_semicolon_taxonomy(verbatim_str)
            if parsed_names:
                species_name = parsed_names[-1]
                if species_name in pr2_dict:
                    aphia_id = pr2_dict[species_name]
                    if aphia_id not in aphia_id_map:
                        aphia_id_map[aphia_id] = []
                    aphia_id_map[aphia_id].append((verbatim_str, assay_name))
                    continue
            stage1_unmatched.append((verbatim_str, assay_name))
        
        unique_aphia_ids_to_fetch = list(aphia_id_map.keys())
        if unique_aphia_ids_to_fetch:
            with mp.Pool(processes=n_proc) as pool:
                worker_func = partial(get_worms_classification_by_id_worker, api_source_for_record=api_source)
                parallel_results = pool.map(worker_func, unique_aphia_ids_to_fetch)

            for aphia_id, result in parallel_results:
                if 'scientificName' in result:
                    for combo in aphia_id_map.get(aphia_id, []):
                        verbatim_str, _ = combo
                        parsed_names = parse_semicolon_taxonomy(verbatim_str)
                        cleaned_taxonomy = ';'.join(parsed_names) if parsed_names else verbatim_str
                        result_with_cleaned = {**result, 'cleanedTaxonomy': cleaned_taxonomy}
                        
                        results_cache[combo] = result_with_cleaned
                        info_records.append({'verbatimIdentification': verbatim_str, **result_with_cleaned, 'ambiguous': False})
                else:
                    stage1_unmatched.extend(aphia_id_map.get(aphia_id, []))
        
        unmatched_tuples = stage1_unmatched
        logging.info(f"Finished Stage 1. Matched {len(results_cache) - len(cases_to_handle)} taxa via AphiaID. Remaining: {len(unmatched_tuples)}.")
    else:
        unmatched_tuples = unique_tuples_to_process

    # --- Stage 2: SMART Parallel Batch Name Matching ---
    if unmatched_tuples:
        all_terms_to_match = list(set(term for vt, an in unmatched_tuples for term in parse_semicolon_taxonomy(vt)))
        logging.info(f"Starting Stage 2: Smart batch matching for {len(all_terms_to_match)} unique terms.")
        
        batch_lookup = {}
        if all_terms_to_match:
            chunk_size = 50
            batch_data = [( (i // chunk_size) + 1, all_terms_to_match[i:i+chunk_size] ) for i in range(0, len(all_terms_to_match), chunk_size)]
            
            logging.info(f"Processing {len(batch_data)} batches with {n_proc} processes...")
            
            with mp.Pool(processes=n_proc) as pool:
                for _, batch_result in pool.map(get_worms_batch_worker, batch_data):
                    batch_lookup.update(batch_result)
            logging.info(f"Smart parallel processing complete! Found matches for {len(batch_lookup)} terms.")

            still_unmatched_batch = []
            for verbatim_str, assay_name in unmatched_tuples:
                parsed_names = parse_semicolon_taxonomy(verbatim_str)
                if not parsed_names:
                    still_unmatched_batch.append((verbatim_str, assay_name))
                    continue
                
                max_depth = assay_rank_info.get(assay_name, {}).get('max_depth', 99)
                is_full_len = len(parsed_names) >= max_depth
                
                if assay_name in assays_to_skip_species and is_full_len:
                    parsed_names = parsed_names[:-1]
                    if not parsed_names:
                        still_unmatched_batch.append((verbatim_str, assay_name))
                        continue
                
                cleaned_taxonomy = ';'.join(parsed_names) if parsed_names else verbatim_str
                match_found = False
                
                for term in reversed(parsed_names):
                    if term in batch_lookup:
                        all_matches = batch_lookup[term]
                        is_ambiguous = len(all_matches) > 1
                        
                        # For main DF, always take the first match
                        results_cache[(verbatim_str, assay_name)] = {
                            **all_matches[0], 'nameAccordingTo': api_source,
                            'match_type_debug': f'Success_Batch_{term}', 'cleanedTaxonomy': cleaned_taxonomy
                        }
                        
                        # For info DF, add all matches
                        for match in all_matches:
                            info_records.append({
                                'verbatimIdentification': verbatim_str, **match,
                                'nameAccordingTo': api_source, 'match_type_debug': f'Success_Batch_{term}',
                                'cleanedTaxonomy': cleaned_taxonomy, 'ambiguous': is_ambiguous
                            })
                        
                        match_found = True
                        break
                
                if not match_found:
                    still_unmatched_batch.append((verbatim_str, assay_name))
            
            unmatched_tuples = still_unmatched_batch
            logging.info(f"Stage 2 complete. Total matched: {len(results_cache)}. Remaining: {len(unmatched_tuples)}")

    # --- Handle Final Unmatched ---
    if unmatched_tuples:
        for verbatim_str, assay_name in unmatched_tuples:
            parsed_names = parse_semicolon_taxonomy(verbatim_str)
            max_depth = assay_rank_info.get(assay_name, {}).get('max_depth', 99)
            is_full_len = len(parsed_names) >= max_depth
            
            if assay_name in assays_to_skip_species and is_full_len and parsed_names:
                parsed_names = parsed_names[:-1]
            
            cleaned_taxonomy = ';'.join(parsed_names) if parsed_names else verbatim_str
            
            no_match_record = {
                'scientificName': 'incertae sedis', 'scientificNameID': 'urn:lsid:marinespecies.org:taxname:12',
                'taxonRank': None, 'nameAccordingTo': api_source,
                'match_type_debug': 'Failed_All_Stages_NoMatch', 'cleanedTaxonomy': cleaned_taxonomy
            }
            for col in DWC_RANKS_STD:
                no_match_record[col] = None

            results_cache[(verbatim_str, assay_name)] = no_match_record
            info_records.append({'verbatimIdentification': verbatim_str, **no_match_record, 'ambiguous': False})

    # --- Create record for initially empty inputs (for main DF) ---
    empty_fallback_record = { 
        'scientificName': 'incertae sedis', 'scientificNameID': 'urn:lsid:marinespecies.org:taxname:12',
        'taxonRank': None, 'nameAccordingTo': api_source,
        'match_type_debug': 'incertae_sedis_truly_empty_fallback', 'cleanedTaxonomy': ''
    }
    for col in DWC_RANKS_STD:
        empty_fallback_record[col] = None
    results_cache['IS_TRULY_EMPTY'] = empty_fallback_record

    # --- Apply results to main occurrence DataFrame ---
    if results_cache:
        mapped_results = df_to_process['_map_key'].map(results_cache)
        results_df = pd.DataFrame(mapped_results.to_list(), index=df_to_process.index)
        for col in results_df.columns:
            df_to_process[col] = results_df[col]

    df_to_process.drop(columns=['_map_key'], inplace=True, errors='ignore')
    
    # --- Create the detailed info DataFrame ---
    info_df = pd.DataFrame(info_records)

    return {'main_df': df_to_process, 'info_df': info_df}

if __name__ == '__main__':
    pass