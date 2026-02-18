import pyworms
import pandas as pd
import multiprocess as mp
import time
from functools import partial
import logging
import re
import sys

# Set up logging to provide clear progress updates
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Standard Darwin Core ranks used for structuring the output
DWC_RANKS_STD = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']

def _normalize_taxon_token(s):
    """
    Normalize taxonomy tokens for loose equality checks.
    Keep it conservative: we only want to match obvious identical terms.
    """
    if s is None or pd.isna(s):
        return None
    s = str(s).strip().lower()
    if not s:
        return None
    # Normalize common separators; keep letters/spaces to avoid over-matching.
    s = s.replace('_', ' ').replace('-', ' ').replace('/', ' ')
    if '  ' in s:
        s = re.sub(r'\s+', ' ', s)
    return s.strip() or None


def _compute_lineage_consistency(parsed_names, candidate_match):
    """
    Returns:
      (score, matched_count, total_count)

    Score in [0, 1]: fraction of verbatim tokens that appear somewhere in the
    candidate's returned classification (including scientificName).

    This is intentionally simple and rank-agnostic because verbatim strings
    often include non-standard ranks (e.g., tribes) and "Eukaryota" as a superkingdom.
    """
    if not parsed_names or not isinstance(candidate_match, dict):
        return 0.0, 0, 0

    # Treat Eukaryota as an optional superkingdom token (WoRMS returns Animalia/Plantae/etc as 'kingdom').
    ignore = {'eukaryota'}
    verbatim_tokens = []
    for t in parsed_names:
        nt = _normalize_taxon_token(t)
        if not nt or nt in ignore:
            continue
        verbatim_tokens.append(nt)
    verbatim_set = set(verbatim_tokens)
    if not verbatim_set:
        return 0.0, 0, 0

    candidate_vals = set()
    sci = _normalize_taxon_token(candidate_match.get('scientificName'))
    if sci:
        candidate_vals.add(sci)
    for rank in DWC_RANKS_STD:
        v = _normalize_taxon_token(candidate_match.get(rank))
        if v:
            candidate_vals.add(v)

    if not candidate_vals:
        return 0.0, 0, len(verbatim_set)

    matched = verbatim_set.intersection(candidate_vals)
    total = len(verbatim_set)
    matched_count = len(matched)
    return matched_count / max(1, total), matched_count, total


def _candidate_completeness(candidate_match):
    """
    Tie-breaker that prefers records with more populated classification fields.
    This avoids biasing toward 'species' purely because it's a more specific rank.
    """
    if not isinstance(candidate_match, dict):
        return 0
    n = 0
    if _normalize_taxon_token(candidate_match.get('scientificName')):
        n += 1
    for rank in DWC_RANKS_STD:
        if _normalize_taxon_token(candidate_match.get(rank)):
            n += 1
    return n


def _get_mp_context():
    """
    On macOS, forking can crash with Objective-C runtime errors like:
    'may have been in progress in another thread when fork() was called'.
    Using 'spawn' avoids fork() and is the safest cross-platform default.
    """
    if sys.platform == 'darwin':
        return mp.get_context('spawn')
    return mp.get_context()

def _safe_str_or_na(x):
    """Coerce to Python str or pd.NA. Avoids Mac/Linux pandas dtype 'str' error from numpy scalars/bytes."""
    if pd.isna(x) or x is None:
        return pd.NA
    return str(x).strip() if isinstance(x, str) else str(x)


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
                'match_type_debug': f'Success_AphiaID_{aphia_id_to_check}', 'name_change': False
            }
            for rank_std in DWC_RANKS_STD:
                result[rank_std] = record.get(rank_std.lower())
            return aphia_id_to_check, result
    except Exception:
        pass
    
    return aphia_id_to_check, {'match_type_debug': f'Failure_AphiaID_{aphia_id_to_check}', 'name_change': False}

def get_worms_batch_worker(batch_info):
    """
    Worker function for parallel batch processing - with timeout safety.
    Now collects ALL accepted matches to handle ambiguous cases.
    """
    batch_num, chunk = batch_info
    batch_results = {}
    
    try:
        batch_results_raw = pyworms.aphiaRecordsByMatchNames(chunk)

        for j, name_list in enumerate(batch_results_raw):
            if not name_list:
                continue

            accepted_by_lsid = {}

            for match in name_list:
                if not match:
                    continue

                status = match.get('status')
                record_to_use = None
                name_changed = False

                # Case 1: already an accepted record
                if status == 'accepted':
                    record_to_use = match

                # Case 2: unaccepted record with a valid (accepted) name behind it.
                # For example: Synagrops spinosus -> Parascombrops spinosus.
                elif status == 'unaccepted':
                    valid_aphia = (
                        match.get('valid_AphiaID')
                        or match.get('valid_aphia_id')
                        or match.get('valid_aphiaID')
                    )
                    if valid_aphia:
                        try:
                            resolved = pyworms.aphiaRecordByAphiaID(valid_aphia)
                            if resolved and resolved.get('status') == 'accepted':
                                record_to_use = resolved
                                name_changed = True
                        except Exception:
                            record_to_use = None

                # Ignore any other statuses
                if not record_to_use:
                    continue

                sci_id = record_to_use.get('lsid')
                if not sci_id:
                    continue

                existing = accepted_by_lsid.get(sci_id)
                if existing is None:
                    res = {
                        'scientificName': record_to_use.get('scientificname'),
                        'scientificNameID': record_to_use.get('lsid'),
                        'taxonRank': record_to_use.get('rank'),
                        'name_change': name_changed
                    }
                    for rank in DWC_RANKS_STD:
                        res[rank] = record_to_use.get(rank.lower())
                    accepted_by_lsid[sci_id] = res
                else:
                    if name_changed and not existing.get('name_change'):
                        existing['name_change'] = True

            accepted_matches = list(accepted_by_lsid.values())

            if accepted_matches:
                batch_results[chunk[j]] = accepted_matches
            
        return batch_num, batch_results

    except Exception:
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
    # Safe string conversion (avoids Mac/pandas dtype 'str' error)
    df_to_process['verbatimIdentification'] = df_to_process['verbatimIdentification'].apply(
        lambda x: '' if pd.isna(x) else str(x).strip()
    )
    if 'assay_name' in df_to_process.columns:
        df_to_process['assay_name'] = df_to_process['assay_name'].apply(
            lambda x: '' if pd.isna(x) else str(x).strip()
        )

    empty_verbatim_mask = (df_to_process['verbatimIdentification'].str.strip() == '') | \
                          (df_to_process['verbatimIdentification'].str.strip() == '') | \
                          (df_to_process['verbatimIdentification'].str.strip().str.lower() == 'unassigned')
    
    # IMPORTANT: This column must be able to hold tuples (verbatimIdentification, assay_name).
    # On some Mac/Linux pandas installs, assigning '' can create a strict string dtype ('str'),
    # which then rejects tuples with: "Invalid value for dtype 'str'".
    df_to_process['_map_key'] = pd.Series([None] * len(df_to_process), index=df_to_process.index, dtype='object')
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
                'cleanedTaxonomy': cleaned_verbatim, 'name_change': False,
                'assignment_score': pd.NA,
                'ranks_matched': pd.NA,
                'ranks_provided': pd.NA
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
            ctx = _get_mp_context()
            with ctx.Pool(processes=n_proc) as pool:
                worker_func = partial(get_worms_classification_by_id_worker, api_source_for_record=api_source)
                parallel_results = pool.map(worker_func, unique_aphia_ids_to_fetch)

            for aphia_id, result in parallel_results:
                if 'scientificName' in result:
                    for combo in aphia_id_map.get(aphia_id, []):
                        verbatim_str, _ = combo
                        parsed_names = parse_semicolon_taxonomy(verbatim_str)
                        cleaned_taxonomy = ';'.join(parsed_names) if parsed_names else verbatim_str
                        consistency_score, matched_count, total_count = _compute_lineage_consistency(parsed_names, result)
                        result_with_cleaned = {
                            **result,
                            'cleanedTaxonomy': cleaned_taxonomy,
                            'assignment_score': consistency_score,
                            'ranks_matched': matched_count,
                            'ranks_provided': total_count
                        }
                        
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
            total_batches = len(batch_data)
            
            logging.info(f"Processing {total_batches} batches with {n_proc} processes...")
            
            processed_batches = 0
            ctx = _get_mp_context()
            with ctx.Pool(processes=n_proc) as pool:
                # Use imap_unordered to get results as they complete, allowing for progress reporting
                for batch_num, batch_result in pool.imap_unordered(get_worms_batch_worker, batch_data):
                    batch_lookup.update(batch_result)
                    processed_batches += 1
                    logging.info(f"  ... Progress: {processed_batches}/{total_batches} batches completed (Batch ID: {batch_num}).")

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
                        
                        # Prefer exact scientificName match to the queried term when available,
                        # then choose by lineage consistency score (tie-breaker: more complete classification).
                        norm_term = _normalize_taxon_token(term)
                        exact_candidates = []
                        for m in all_matches:
                            if _normalize_taxon_token(m.get('scientificName')) == norm_term:
                                exact_candidates.append(m)
                        candidates = exact_candidates if exact_candidates else all_matches

                        scored_candidates = []
                        for m in candidates:
                            s, matched_count, total_count = _compute_lineage_consistency(parsed_names, m)
                            scored_candidates.append((s, _candidate_completeness(m), matched_count, total_count, m))
                        scored_candidates.sort(key=lambda x: (x[0], x[1]), reverse=True)
                        best_score, _, best_matched_count, best_total_count, best_match = scored_candidates[0]

                        results_cache[(verbatim_str, assay_name)] = {
                            **best_match,
                            'nameAccordingTo': api_source,
                            'match_type_debug': f'Success_Batch_{term}',
                            'cleanedTaxonomy': cleaned_taxonomy,
                            'assignment_score': best_score,
                            'ranks_matched': best_matched_count,
                            'ranks_provided': best_total_count
                        }
                        
                        # For info DF, add all matches
                        for match in all_matches:
                            s, matched_count, total_count = _compute_lineage_consistency(parsed_names, match)
                            info_records.append({
                                'verbatimIdentification': verbatim_str, **match,
                                'nameAccordingTo': api_source, 'match_type_debug': f'Success_Batch_{term}',
                                'cleanedTaxonomy': cleaned_taxonomy,
                                'assignment_score': s,
                                'ranks_matched': matched_count,
                                'ranks_provided': total_count,
                                'ambiguous': is_ambiguous
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
                'match_type_debug': 'Failed_All_Stages_NoMatch', 'cleanedTaxonomy': cleaned_taxonomy,
                'name_change': False,
                'assignment_score': pd.NA,
                'ranks_matched': pd.NA,
                'ranks_provided': pd.NA
            }
            for col in DWC_RANKS_STD:
                no_match_record[col] = None

            results_cache[(verbatim_str, assay_name)] = no_match_record
            info_records.append({'verbatimIdentification': verbatim_str, **no_match_record, 'ambiguous': False})

    # --- Create record for initially empty inputs (for main DF) ---
    empty_fallback_record = { 
        'scientificName': 'incertae sedis', 'scientificNameID': 'urn:lsid:marinespecies.org:taxname:12',
        'taxonRank': None, 'nameAccordingTo': api_source,
        'match_type_debug': 'incertae_sedis_truly_empty_fallback', 'cleanedTaxonomy': '',
        'name_change': False,
        'assignment_score': pd.NA,
        'ranks_matched': pd.NA,
        'ranks_provided': pd.NA
    }
    for col in DWC_RANKS_STD:
        empty_fallback_record[col] = None
    results_cache['IS_TRULY_EMPTY'] = empty_fallback_record

    # --- Apply results to main occurrence DataFrame ---
    if results_cache:
        mapped_results = df_to_process['_map_key'].map(results_cache)
        results_df = pd.DataFrame(mapped_results.to_list(), index=df_to_process.index)
        for col in results_df.columns:
            results_df[col] = results_df[col].apply(_safe_str_or_na)
        for col in results_df.columns:
            df_to_process[col] = results_df[col]

    df_to_process.drop(columns=['_map_key'], inplace=True, errors='ignore')
    
    # --- Create the detailed info DataFrame ---
    info_df = pd.DataFrame(info_records)
    for col in info_df.columns:
        info_df[col] = info_df[col].apply(_safe_str_or_na)

    return {'main_df': df_to_process, 'info_df': info_df}

if __name__ == '__main__':
    pass