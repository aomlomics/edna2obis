import pyworms
import pandas as pd
import multiprocess as mp # Use multiprocess instead of multiprocessing for notebook compatibility
import os
import time
from functools import partial

# Standard Darwin Core ranks
DWC_RANKS_STD = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']

def parse_semicolon_taxonomy(tax_string):
    """Helper function to parse a semicolon-separated taxonomy string into a list of names."""
    if pd.isna(tax_string) or not str(tax_string).strip():
        return []
    return [name.strip() for name in str(tax_string).split(';') if name.strip()]

def get_worms_match_for_single_taxon_worker(combo_input_for_worker, assays_to_skip_species_list, api_source_for_record='WoRMS'):
    """
    Core logic to match a single taxonomy string (verbatimIdentification) against WoRMS,
    considering the assay_name for species-skipping rules.
    Returns a dictionary with WoRMS classification details and a debug log.
    """
    verbatim_identification_str, assay_name_str = combo_input_for_worker
    debug_log = [f"--- Starting Worker for: verbatim='{verbatim_identification_str}', assay='{assay_name_str}' ---"]

    # Base result if no match is found or an error occurs
    base_result_on_failure = {
        'scientificName': 'incertae sedis',
        'scientificNameID': 'urn:lsid:marinespecies.org:taxname:12', # AphiaID for Incertae Sedis
        'taxonRank': 'no rank',
        **{rank: None for rank in DWC_RANKS_STD}, # Fill standard ranks with None
        'nameAccordingTo': api_source_for_record,
        'match_type_debug': 'WorkerError_Or_InitialNoMatch'
    }

    parsed_names = parse_semicolon_taxonomy(verbatim_identification_str)
    if not parsed_names:
        debug_log.append("No parsed names found from verbatimIdentification. Returning 'incertae sedis'.")
        base_result_on_failure['match_type_debug'] = 'No_Tax_String_Provided_Or_Parsed_Empty'
        base_result_on_failure['debug_log'] = debug_log
        return base_result_on_failure
    
    debug_log.append(f"Parsed names: {parsed_names}")

    # WoRMS/DwC ranks, from most specific to most broad, to guide iteration
    # We iterate through the *parsed names* from most specific to least specific
    # and try to match them at their corresponding *effective rank*.
    dwc_rank_hierarchy_for_iteration = ['species', 'genus', 'family', 'order', 'class', 'phylum', 'kingdom']
    num_parsed_levels = len(parsed_names)

    for i in range(num_parsed_levels):
        # Iterate from the most specific part of the parsed string
        name_to_match_idx_in_parsed = num_parsed_levels - 1 - i
        name_to_match = parsed_names[name_to_match_idx_in_parsed]
        
        # Determine the 'effective rank' we are trying to match this specific name at.
        # This is based on its position in the *parsed string*.
        effective_rank_being_tried = dwc_rank_hierarchy_for_iteration[i] if i < len(dwc_rank_hierarchy_for_iteration) else f"level_{i+1}_from_specific_end"

        debug_log.append(f"Attempt {i+1}/{num_parsed_levels}: Trying term='{name_to_match}' (from parsed index {name_to_match_idx_in_parsed}). Effective rank context: '{effective_rank_being_tried}'.")

        if not name_to_match or name_to_match.lower() == 'unassigned' or name_to_match.lower() == 'incertae sedis': # Skip common placeholders
            debug_log.append(f"Skipping term '{name_to_match}' as it's a placeholder or empty.")
            continue
            
        # Check if this is the species-level term (i=0) and if the assay requires skipping species
        if assay_name_str in assays_to_skip_species_list and i == 0 : # i=0 corresponds to the most specific term (potential species)
            debug_log.append(f"Skipping API query for most specific term (potential species) '{name_to_match}' due to assay '{assay_name_str}' being in assays_to_skip_species_list.")
            continue
        
        s_match_list = None
        try:
            time.sleep(0.05) # Small delay to be kind to the API
            debug_log.append(f"Querying WoRMS API for: '{name_to_match}' (marine_only=True, like=False)")
            s_match_list = pyworms.aphiaRecordsByName(name_to_match, like=False, marine_only=True)
            debug_log.append(f"Raw WoRMS API response for '{name_to_match}': {s_match_list}")
            
            chosen_match_record = None
            if s_match_list and isinstance(s_match_list, list):
                # Prefer accepted matches whose rank aligns with the effective rank we're trying
                for m_dict in s_match_list:
                    if m_dict and isinstance(m_dict, dict) and m_dict.get('status') == 'accepted':
                        worms_rank_of_match = str(m_dict.get('rank', '')).lower()
                        if worms_rank_of_match == effective_rank_being_tried.lower():
                            chosen_match_record = m_dict
                            debug_log.append(f"Found ACCEPTED match with WoRMS rank '{worms_rank_of_match}' aligning with effective rank '{effective_rank_being_tried}' for '{name_to_match}'. Name: {m_dict.get('scientificname')}, AphiaID: {m_dict.get('AphiaID')}")
                            break 
                
                # If no rank-aligned match, take the first accepted one
                if not chosen_match_record:
                    for m_dict in s_match_list:
                        if m_dict and isinstance(m_dict, dict) and m_dict.get('status') == 'accepted':
                            chosen_match_record = m_dict
                            debug_log.append(f"No rank-aligned match. Found first available ACCEPTED match (WoRMS rank: {m_dict.get('rank')}) for '{name_to_match}'. Name: {m_dict.get('scientificname')}, AphiaID: {m_dict.get('AphiaID')}")
                            break
            
            if chosen_match_record:
                debug_log.append(f"Chosen accepted_match_record for '{name_to_match}': ID {chosen_match_record.get('AphiaID')}, Name {chosen_match_record.get('scientificname')}, Rank {chosen_match_record.get('rank')}")
                
                # Populate the result dictionary
                result_dict = {
                    'scientificName': chosen_match_record.get('scientificname'),
                    'scientificNameID': chosen_match_record.get('lsid'), # LSID
                    'taxonRank': chosen_match_record.get('rank'), 
                    'nameAccordingTo': api_source_for_record, # e.g., "WoRMS"
                    'match_type_debug': f'Success_Query_{name_to_match}_EffectiveRank_{effective_rank_being_tried}_WoRMSRank_{chosen_match_record.get("rank")}'
                }
                # Populate the standard DwC rank columns
                for rank_std in DWC_RANKS_STD: # ['kingdom', 'phylum', ..., 'species']
                    result_dict[rank_std] = chosen_match_record.get(rank_std.lower()) # WoRMS API returns ranks in lowercase
                
                debug_log.append(f"SUCCESS for verbatim='{verbatim_identification_str}'. Matched on term '{name_to_match}'. Resulting WoRMS scientificName: {result_dict.get('scientificName')}, Rank: {result_dict.get('taxonRank')}")
                result_dict['debug_log'] = debug_log
                return result_dict
            else:
                debug_log.append(f"No accepted match record derived from WoRMS response for '{name_to_match}'. Raw response was: {s_match_list}")
        
        except Exception as e:
            debug_log.append(f"WoRMS API call or processing error for '{name_to_match}': {e}")
            # Continue to try the next broader term
    
    # If loop completes without returning a match
    debug_log.append(f"No definitive match found for verbatim='{verbatim_identification_str}' after trying all parsed parts. Returning 'incertae sedis'.")
    base_result_on_failure['match_type_debug'] = 'NoMatchFound_TriedAllParts_IncertaeSedis'
    base_result_on_failure['debug_log'] = debug_log
    return base_result_on_failure

def worker_task_wrapper(task_combo_and_params):
    """ Unpacks arguments for the pool worker. """
    combo_input, assays_skip_list, api_src = task_combo_and_params
    return get_worms_match_for_single_taxon_worker(combo_input, assays_skip_list, api_src)

def get_worms_match_for_dataframe(occurrence_df, params_dict, n_proc=0):
    """
    Adds WoRMS taxonomic information to the occurrence DataFrame.
    Uses multiprocess.Pool for parallelization.
    Returns the updated DataFrame. Cache with debug logs is handled internally.
    """
    if not isinstance(occurrence_df, pd.DataFrame) or occurrence_df.empty:
        print("Input DataFrame is empty. Skipping WoRMS taxonomic matching.")
        return occurrence_df.copy()

    api_source = params_dict.get('taxonomic_api_source', 'WoRMS')
    assays_to_skip_species = params_dict.get('assays_to_skip_species_match', []) # From main params

    required_cols = ['verbatimIdentification', 'assay_name']
    for col in required_cols:
        if col not in occurrence_df.columns:
            print(f"ERROR: DataFrame must contain '{col}' column for WoRMS matching.")
            # Add empty columns for results to avoid downstream errors if processing continues
            for res_col in ['scientificName', 'scientificNameID', 'taxonRank', 'nameAccordingTo', 'match_type_debug'] + DWC_RANKS_STD:
                if res_col not in occurrence_df.columns:
                    occurrence_df[res_col] = None
            return occurrence_df
            
    df_to_process = occurrence_df.copy()
    df_to_process['verbatimIdentification'] = df_to_process['verbatimIdentification'].astype(str).fillna('')

    # Create unique combinations of (verbatimIdentification, assay_name) to query
    unique_verbatim_assay_tuples = list(df_to_process[['verbatimIdentification', 'assay_name']].drop_duplicates().itertuples(index=False, name=None))
    unique_verbatim_assay_tuples = [combo for combo in unique_verbatim_assay_tuples if str(combo[0]).strip()] # Filter out effectively empty verbatim strings
    
    print(f"Found {len(unique_verbatim_assay_tuples)} unique, non-empty (verbatimIdentification, assay_name) combinations for {api_source} matching.")

    # --- Parallel Processing with multiprocess.Pool ---
    results_from_api = {} # This will store: {(verbatim, assay): result_dict}
    
    if unique_verbatim_assay_tuples:
        num_processes_to_use = mp.cpu_count() if n_proc == 0 else int(n_proc)
        num_processes_to_use = min(num_processes_to_use, len(unique_verbatim_assay_tuples))
        if num_processes_to_use == 0 and len(unique_verbatim_assay_tuples) > 0: num_processes_to_use = 1

        print(f"Starting {api_source} queries with {num_processes_to_use} processes using multiprocess.Pool...")

        # Prepare arguments for the pool worker
        # Each item in tasks_for_pool will be: ((verbatim_str, assay_name_str), assays_to_skip_species_list, api_source)
        tasks_for_pool = [(combo, assays_to_skip_species, api_source) for combo in unique_verbatim_assay_tuples]

        with mp.Pool(processes=num_processes_to_use) as pool:
            # map will block until all results are processed
            # The worker_task_wrapper gets a single item from tasks_for_pool
            # and returns the result from get_worms_match_for_single_taxon_worker
            # The order of 'raw_results_list' will correspond to 'unique_verbatim_assay_tuples'
            raw_results_list = pool.map(worker_task_wrapper, tasks_for_pool)
        
        # Populate the results_from_api dictionary
        for i, combo_key in enumerate(unique_verbatim_assay_tuples):
            results_from_api[combo_key] = raw_results_list[i]
            # Optional: print debug logs from worker if needed, e.g., results_from_api[combo_key].get('debug_log')
            if (i + 1) % 200 == 0 or (i + 1) == len(unique_verbatim_assay_tuples):
                 print(f"  Processed {i+1}/{len(unique_verbatim_assay_tuples)} unique combinations from API/pool.")
    else:
        print("No unique combinations to process with the API.")

    # --- Map results back to the DataFrame ---
    # Create temporary columns from the tuple (verbatimIdentification, assay_name) to use for mapping
    df_to_process['_map_key'] = list(zip(df_to_process['verbatimIdentification'], df_to_process['assay_name']))
    
    # Initialize result columns
    result_cols_to_add = ['scientificName', 'scientificNameID', 'taxonRank', 'nameAccordingTo', 'match_type_debug'] + DWC_RANKS_STD
    for col in result_cols_to_add:
        df_to_process[col] = None

    # Apply cached results
    for key_tuple, result_data in results_from_api.items():
        mask = (df_to_process['_map_key'] == key_tuple)
        for res_col, res_val in result_data.items():
            if res_col in df_to_process.columns and res_col != 'debug_log': # Don't try to assign the debug_log list to a column
                df_to_process.loc[mask, res_col] = res_val
    
    df_to_process.drop(columns=['_map_key'], inplace=True)
    
    print(f"\nFinished applying {api_source} taxonomic matches to DataFrame.")
    return df_to_process

if __name__ == '__main__':
    # This block is for testing the script directly, not used when imported as a module.
    print("WoRMS_OBIS_matcher.py script is being run directly (e.g., for testing).")
    # Example usage for testing:
    # test_data = {
    #     'verbatimIdentification': [
    #         "Bacteria;Proteobacteria;Gammaproteobacteria;Alteromonadales;Alteromonadaceae;Alteromonas;Alteromonas_macleodii",
    #         "Eukaryota;Excavata;Discoba;Euglenozoa;Euglenida;Rhabdomonadales",
    #         "Bacteria;Firmicutes",
    #         "Unassigned",
    #         "Eukaryota;TSAR;Alveolata;Dinoflagellata;Dinophyceae;Gymnodiniales;Gymnodiniaceae;Gymnodinium;Gymnodinium_catenatum" # Will be species if not skipped
    #     ],
    #     'assay_name': ["16S", "18S", "16S", "16S", "18S_species_ok"], # Example assay names
    #     'occurrenceID': ["occ1", "occ2", "occ3", "occ4", "occ5"]
    # }
    # test_df = pd.DataFrame(test_data)
    # test_params = {
    #     'taxonomic_api_source': 'WoRMS',
    #     'assays_to_skip_species_match': ['16S'], # 16S will skip species, 18S will not
    #     'output_dir': './temp_test_output/' # Ensure this dir exists or is created
    # }
    # os.makedirs(test_params['output_dir'], exist_ok=True)

    # print("\n--- Running Test ---")
    # matched_df = get_worms_match_for_dataframe(test_df, test_params, n_proc=2)
    # print("\n--- Test Results (Matched DataFrame Head) ---")
    # print(matched_df.head())
    # print("\n--- Columns in Test Matched DataFrame ---")
    # print(matched_df.columns)

    # print("\n--- Specific Checks ---")
    # if not matched_df.empty:
    #     print("\nGymnodinium catenatum (should be species):")
    #     print(matched_df[matched_df['occurrenceID'] == 'occ5'][['verbatimIdentification', 'scientificName', 'taxonRank', 'match_type_debug'] + DWC_RANKS_STD])
        
    #     print("\nAlteromonas macleodii (16S, should be genus or higher due to skip):")
    #     print(matched_df[matched_df['occurrenceID'] == 'occ1'][['verbatimIdentification', 'scientificName', 'taxonRank', 'match_type_debug'] + DWC_RANKS_STD])
