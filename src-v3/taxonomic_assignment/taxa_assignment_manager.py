"""
Taxonomic Assignment Manager for edna2obis
Contains ALL original notebook functionality for taxonomic assignment including:
- PR2 database optimization (local reference database)
- Assay rank determination 
- WoRMS and GBIF API matching
- WoRMS-specific post-processing
- Progress capture and reporting
"""

import pandas as pd
import numpy as np
import os
import io
import sys
import logging
import importlib
import traceback
from contextlib import redirect_stdout, redirect_stderr

# Import the API-specific matching scripts
from . import WoRMS_v3_matching
from . import GBIF_matching
from create_occurrence_core.occurrence_builder import get_final_occurrence_column_order


def load_pr2_worms_dict_into_params(params, reporter=None, warn_print=None):
    """
    When use_local_reference_database is True and API is WoRMS, load PR2 (or similar) Excel
    into params['pr2_worms_dict']. Shared by assign_taxonomy and taxassign.
    warn_print: optional callable(str) when reporter is None (e.g. CLI).
    """
    def _warn(msg):
        if reporter:
            reporter.add_warning(msg)
        elif warn_print:
            warn_print(msg)

    use_local_db = params.get('use_local_reference_database', False)
    api_source = params.get('taxonomic_api_source', 'WoRMS')

    if use_local_db and api_source == 'WoRMS':
        local_db_path = params.get('local_reference_database_path')
        if local_db_path and os.path.exists(local_db_path):
            try:
                if reporter:
                    reporter.add_text(f"Loading local reference database from: {local_db_path}")
                elif warn_print:
                    warn_print(f"Loading local reference database from: {local_db_path}")

                local_df = pd.read_excel(local_db_path, index_col=None, na_values=[""])
                local_df.dropna(subset=['worms_id', 'species'], inplace=True)

                species_cleaned = local_df['species'].str.replace('_', ' ', regex=False)
                species_cleaned = species_cleaned.str.replace(' sp.', '', regex=False).str.strip()

                params['pr2_worms_dict'] = dict(zip(species_cleaned, local_df['worms_id'].astype(int)))

                if reporter:
                    reporter.add_success(
                        f"Successfully loaded local reference database with {len(params['pr2_worms_dict'])} AphiaID mappings"
                    )
                    reporter.add_text(
                        "This will significantly speed up taxonomic matching via direct AphiaID lookup!"
                    )
                elif warn_print:
                    warn_print(
                        f"Loaded local reference database: {len(params['pr2_worms_dict'])} AphiaID mappings"
                    )

            except Exception as e:
                _warn(f"Could not load local reference database: {e}")
                params['pr2_worms_dict'] = {}
        else:
            _warn(f"Local reference database enabled but file not found at: {local_db_path}")
            params['pr2_worms_dict'] = {}
    else:
        if use_local_db and api_source != 'WoRMS':
            _warn("Local reference database optimization is only available for WoRMS API")
        params['pr2_worms_dict'] = {}


def assign_taxonomy(params, data, raw_data_tables, reporter):
    """
    Perform taxonomic assignment using WoRMS or GBIF API
    EXACT implementation from original notebook with ALL functionality preserved
    """
    try:
        reporter.add_section("Taxonomic Assignment")
        
        importlib.reload(WoRMS_v3_matching)
        importlib.reload(GBIF_matching)
        
        # Load the occurrence file generated earlier
        occurrence_path = os.path.join(params.get('output_dir', '../processed-v3/'), 'occurrence.csv')
        if not os.path.exists(occurrence_path):
            reporter.add_error(f"Occurrence file not found at {occurrence_path}")
            return
            
        df_to_match = pd.read_csv(occurrence_path)
        # Coerce key columns to plain strings (avoids Mac/pandas dtype 'str' error downstream)
        for col in ['verbatimIdentification', 'assay_name']:
            if col in df_to_match.columns:
                df_to_match[col] = df_to_match[col].apply(
                    lambda x: '' if pd.isna(x) else str(x).strip()
                )
        reporter.add_text(f"Loaded occurrence data: {len(df_to_match):,} records")
        
        api_source = params.get('taxonomic_api_source', 'WoRMS')
        reporter.add_text(f"Using API source: {api_source}")
        reporter.add_text(f"Assays configured to skip species-level matching: {params.get('assays_to_skip_species_match', [])}")
        
        load_pr2_worms_dict_into_params(params, reporter)

        # Determine maximum taxonomic ranks for each assay
        # Determine maximum taxonomic ranks for each assay (silent)
        
        assay_rank_info = {}
        if raw_data_tables and data.get('analysis_data_by_assay'):
            for analysis_run_name, tables in raw_data_tables.items():
                if 'taxonomy' in tables:
                    tax_df = tables['taxonomy']
                    # Find which assay this analysis run belongs to
                    assay_name = None
                    for assay, runs in data['analysis_data_by_assay'].items():
                        if analysis_run_name in runs:
                            assay_name = assay
                            break
                    
                    if assay_name:
                        try:
                            # Find columns between the two guaranteed markers
                            start_idx = tax_df.columns.get_loc('verbatimIdentification') + 1
                            end_idx = tax_df.columns.get_loc('Confidence')
                            rank_cols = tax_df.columns[start_idx:end_idx].tolist()
                            
                            assay_rank_info[assay_name] = {'max_depth': len(rank_cols)}
                        except (KeyError, ValueError):
                            assay_rank_info[assay_name] = {'max_depth': 7}
        
        # Fallback for any missing assays
        for assay in df_to_match['assay_name'].unique():
            if assay not in assay_rank_info:
                assay_rank_info[assay] = {'max_depth': 7}
                pass
        
        params['assay_rank_info'] = assay_rank_info
        
        # Create progress capture system
        # Create a custom log handler to capture progress
        log_capture_string = io.StringIO()
        ch = logging.StreamHandler(log_capture_string)
        ch.setLevel(logging.INFO)
        formatter = logging.Formatter('%(levelname)s - %(message)s')
        ch.setFormatter(formatter)
        
        # Get the logger used by the API scripts (they likely use logging)
        logger = logging.getLogger()
        # The following line is the root cause of the logs not appearing on the console.
        # It is temporarily disabled to allow for real-time debugging in the terminal.
        # logger.addHandler(ch)
        logger.setLevel(logging.INFO)
        
        # Also capture stdout/stderr for any print statements
        captured_output = io.StringIO()
        
        reporter.add_text("Starting taxonomic matching process...")
        
        try:
            # Run the matching with both logging and stdout/stderr capture
            if api_source == 'WoRMS':
                with redirect_stdout(captured_output), redirect_stderr(captured_output):
                    # The function now returns a dict with 'main_df' and 'info_df'
                    worms_results = WoRMS_v3_matching.get_worms_match_for_dataframe(
                        occurrence_df=df_to_match,
                        params_dict=params,
                        n_proc=params.get('worms_n_proc', 0)
                    )
                matched_df = worms_results['main_df']
                # Store the detailed info df in params to pass it to the next function
                params['taxa_info_df'] = worms_results['info_df']
            elif api_source == 'GBIF':
                # The 'with' block that captures output is temporarily disabled for debugging
                # to allow real-time logs to appear on the console.
                # with redirect_stdout(captured_output), redirect_stderr(captured_output):
                gbif_results = GBIF_matching.get_gbif_match_for_dataframe(
                    occurrence_df=df_to_match,
                    params_dict=params,
                    n_proc=params.get('gbif_n_proc', 0)
                )
                matched_df = gbif_results['main_df']
                params['taxa_info_df'] = gbif_results['info_df']
            else:
                raise ValueError(f"Unknown taxonomic API source: {api_source}")
            
            # Get captured output from both sources
            stdout_output = captured_output.getvalue()
            log_output = log_capture_string.getvalue()
            
            # Combine all progress output
            all_progress = []
            if stdout_output.strip():
                all_progress.append("Console Output:")
                all_progress.append(stdout_output.strip())
            if log_output.strip():
                all_progress.append("Log Output:")
                all_progress.append(log_output.strip())
            
            if all_progress:
                progress_text = "\n".join(all_progress)
                # Add to HTML report only
                reporter.add_text("Taxonomic matching progress:")
                reporter.add_text(f"<pre>{progress_text}</pre>")
            else:
                reporter.add_text("Taxonomic matching completed")
            
        except Exception as e:
            # If there's an error, still get any captured output
            stdout_output = captured_output.getvalue()
            log_output = log_capture_string.getvalue()
            if stdout_output or log_output:
                reporter.add_text(f"Progress before error:<pre>{stdout_output}\n{log_output}</pre>")
            raise e
        finally:
            # Clean up the logger
            # logger.removeHandler(ch)
            pass
        
        if matched_df is not None and not matched_df.empty:
            # Apply post-processing for WoRMS (the manual corrections)
            if api_source == 'WoRMS':
                reporter.add_text("Applying manual taxonomic corrections...")

                match_type_debug_text = matched_df['match_type_debug'].astype('string').fillna('')
                scientific_name_text = matched_df['scientificName'].astype('string').fillna('')
                kingdom_text = matched_df['kingdom'].astype('string').fillna('')
                
                # Define all taxonomic rank columns that need to be managed
                ALL_TAX_RANKS = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'higherClassification', 'taxonRank']

                # CASE 1: Handle records where WoRMS lookup completely failed.
                # The script now assigns these to 'incertae sedis' directly, so we don't need to change them.
                # But we can log how many we have for reporting purposes.
                no_match_mask = match_type_debug_text == 'Failed_All_Stages_NoMatch'
                num_no_match = no_match_mask.sum()
                if num_no_match > 0:
                    reporter.add_text(f"Found {num_no_match:,} records assigned to 'incertae sedis' due to no WoRMS match.")
                
                # Check for pre-handled cases
                pre_handled_mask = (
                    match_type_debug_text.str.contains('incertae_sedis_simple_case', na=False) |
                    match_type_debug_text.str.contains('incertae_sedis_unassigned', na=False)
                )
                num_pre_handled = pre_handled_mask.sum()
                if num_pre_handled > 0:
                    reporter.add_text(f"Found {num_pre_handled:,} pre-handled cases (unassigned/empty/simple kingdoms) assigned to 'incertae sedis'.")

                # CASE 2: Handle specific high-level Eukaryota assignments that should be 'incertae sedis'.
                # Only reassign to incertae sedis if the kingdom is specifically 'Eukaryota' and scientificName is also 'Eukaryota'
                eukaryota_override_mask = (
                    (kingdom_text == 'Eukaryota') &
                    (scientific_name_text == 'Eukaryota') &
                    (~match_type_debug_text.str.contains('incertae_sedis_simple_case', na=False)) &
                    (~match_type_debug_text.str.contains('incertae_sedis_unassigned', na=False)))
                num_eukaryota_override = eukaryota_override_mask.sum()
                
                reporter.add_text(f"Reassigned {num_eukaryota_override:,} complex Eukaryota records to 'incertae sedis'.")
                matched_df.loc[eukaryota_override_mask, 'scientificName'] = 'incertae sedis'

                # CASE 3: Handle any remaining empty/NaN scientificName records as 'incertae sedis'.
                # This catches any edge cases that might have slipped through
                nan_mask = matched_df['scientificName'].isna()
                num_nan = nan_mask.sum()
                if num_nan > 0:
                    reporter.add_text(f"Assigned {num_nan:,} remaining empty records to 'incertae sedis'.")
                    matched_df.loc[nan_mask, 'scientificName'] = 'incertae sedis'
                    matched_df.loc[nan_mask, 'scientificNameID'] = 'urn:lsid:marinespecies.org:taxname:12'
                    # Clear all other taxonomic columns for these records
                    for rank_col in ALL_TAX_RANKS:
                         if rank_col in matched_df.columns:
                            # Use pd.NA for string-typed columns (pandas 'str' dtype is strict on Mac/Linux)
                            matched_df.loc[nan_mask, rank_col] = pd.NA
            
            elif api_source == 'GBIF':
                reporter.add_text("Applying manual taxonomic corrections for GBIF...")

                match_type_debug_text = matched_df['match_type_debug'].astype('string').fillna('')
                scientific_name_text = matched_df['scientificName'].astype('string').fillna('')
                kingdom_text = matched_df['kingdom'].astype('string').fillna('')
                
                # Define all taxonomic rank columns that need to be managed (GBIF doesn't have scientificNameID)
                ALL_TAX_RANKS_GBIF = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'taxonRank']

                # CASE 1: Handle records where GBIF lookup completely failed.
                no_match_mask = match_type_debug_text == 'No_GBIF_Match'
                num_no_match = no_match_mask.sum()
                if num_no_match > 0:
                    reporter.add_text(f"Found {num_no_match:,} records assigned to 'incertae sedis' due to no GBIF match.")
                
                # Check for pre-handled cases
                pre_handled_mask = (
                    match_type_debug_text.str.contains('incertae_sedis_simple_case', na=False) |
                    match_type_debug_text.str.contains('incertae_sedis_unassigned', na=False)
                )
                num_pre_handled = pre_handled_mask.sum()
                if num_pre_handled > 0:
                    reporter.add_text(f"Found {num_pre_handled:,} pre-handled cases (unassigned/empty/simple kingdoms) assigned to 'incertae sedis'.")

                # CASE 2: Handle specific high-level Eukaryota assignments that should be 'incertae sedis'.
                # Only reassign to incertae sedis if the kingdom is specifically 'Eukaryota' and scientificName is also 'Eukaryota'
                eukaryota_override_mask = (
                    (kingdom_text == 'Eukaryota') &
                    (scientific_name_text == 'Eukaryota') &
                    (~match_type_debug_text.str.contains('incertae_sedis_simple_case', na=False)) &
                    (~match_type_debug_text.str.contains('incertae_sedis_unassigned', na=False)))
                num_eukaryota_override = eukaryota_override_mask.sum()
                
                reporter.add_text(f"Reassigned {num_eukaryota_override:,} complex Eukaryota records to 'incertae sedis'.")
                matched_df.loc[eukaryota_override_mask, 'scientificName'] = 'incertae sedis'

                # CASE 3: Handle any remaining empty/NaN scientificName records as 'incertae sedis'.
                nan_mask = matched_df['scientificName'].isna()
                num_nan = nan_mask.sum()
                if num_nan > 0:
                    reporter.add_text(f"Assigned {num_nan:,} remaining empty records to 'incertae sedis'.")
                    matched_df.loc[nan_mask, 'scientificName'] = 'incertae sedis'
                    # No scientificNameID for GBIF
                    # Clear all other taxonomic columns for these records
                    for rank_col in ALL_TAX_RANKS_GBIF:
                         if rank_col in matched_df.columns:
                            # Use pd.NA for string-typed columns (pandas 'str' dtype is strict on Mac/Linux)
                            matched_df.loc[nan_mask, rank_col] = pd.NA
            
            # Post-matching processing and final save
            reporter.add_text("Starting post-matching processing.")
            
            # Remove temporary columns and save
            columns_to_drop = ['assay_name']
            if api_source == 'GBIF' and 'scientificNameID' in matched_df.columns:
                columns_to_drop.append('scientificNameID')
                reporter.add_text("Adjusted final columns for GBIF standard: Removed 'scientificNameID'.")
            
            matched_df = matched_df.drop(columns=columns_to_drop, errors='ignore')
            reporter.add_text(f"Removed temporary columns: {columns_to_drop}")
            
            # Now, save the final matched dataframe with a new name
            final_occurrence_path = os.path.join(params.get('output_dir', '../processed-v3/'), f'occurrence_core_{api_source.lower()}.csv')
            
            # --- Reorder columns to final desired spec, removing API-specific IDs that are not relevant ---
            final_col_order = get_final_occurrence_column_order()
            
            # Conditionally remove columns that are not applicable to the current API source
            if api_source.lower() == 'gbif':
                if 'scientificNameID' in final_col_order:
                    final_col_order.remove('scientificNameID')
            elif api_source.lower() == 'worms':
                if 'taxonID' in final_col_order:
                    final_col_order.remove('taxonID')

            # Filter the master list to only include columns that actually exist in the dataframe
            cols_to_keep = [col for col in final_col_order if col in matched_df.columns]
            
            # Select and reorder
            matched_df_final = matched_df[cols_to_keep]

            matched_df_final.to_csv(final_occurrence_path, index=False, na_rep='')
            reporter.add_success(f"Taxonomic assignment completed! Saved {len(matched_df_final):,} records to occurrence_core_{api_source.lower()}.csv")
            
            # Add summary statistics
            unique_taxa = matched_df['scientificName'].nunique()
            reporter.add_text(f"Summary: {unique_taxa:,} unique taxa identified")
            
        else:
            reporter.add_error("Taxonomic assignment returned empty results")
            
    except Exception as e:
        reporter.add_error(f"Taxonomic assignment failed: {str(e)}")
        # Include full traceback in the HTML report for cross-platform debugging
        reporter.add_text(f"<pre>{traceback.format_exc()}</pre>")


def mark_selected_match_from_main_dataframe(info_df: pd.DataFrame, main_df: pd.DataFrame, api: str) -> pd.DataFrame:
    """
    Mark selected_match True for rows matching the matcher's chosen taxon per verbatimIdentification
    (main_df), equivalent to mark_selected_* using occurrence_core but without an occurrence file.
    """
    if info_df.empty:
        return info_df
    out = info_df.copy()
    out['selected_match'] = False
    if main_df is None or main_df.empty:
        return out

    api_norm = (api or '').strip()
    if api_norm == 'WoRMS':
        key_cols = ['verbatimIdentification', 'scientificNameID']
        if not all(c in main_df.columns for c in key_cols) or not all(c in out.columns for c in key_cols):
            return out
        selected_pairs = main_df[key_cols].drop_duplicates()
    elif api_norm == 'GBIF':
        key_cols = ['verbatimIdentification', 'taxonID']
        if not all(c in main_df.columns for c in key_cols) or not all(c in out.columns for c in key_cols):
            return out
        selected_pairs = main_df[key_cols].drop_duplicates()
    else:
        return out

    for col in key_cols:
        selected_pairs[col] = selected_pairs[col].astype(str)
        out[col] = out[col].astype(str)

    out['_info_row_idx'] = out.index
    merged = pd.merge(out, selected_pairs, on=key_cols, how='inner')
    if merged.empty:
        out = out.drop(columns=['_info_row_idx'], errors='ignore')
        return out

    if api_norm == 'GBIF' and 'confidence' in merged.columns:
        merged = merged.sort_values('confidence', ascending=False)
    best_candidates = merged.drop_duplicates(subset=key_cols, keep='first')
    out.loc[best_candidates['_info_row_idx'], 'selected_match'] = True
    out = out.drop(columns=['_info_row_idx'], errors='ignore')
    return out


def limit_info_df_preserving_selected(info_df: pd.DataFrame, match_limit: int) -> pd.DataFrame:
    """Cap rows per verbatimIdentification; always keeps selected_match rows, then fills to match_limit."""
    if info_df.empty or match_limit <= 0:
        return info_df
    if 'selected_match' not in info_df.columns:
        return (
            info_df.groupby('verbatimIdentification', group_keys=False)
            .head(match_limit)
            .reset_index(drop=True)
        )

    parts = []
    for _, g in info_df.groupby('verbatimIdentification', sort=False):
        sel = g[g['selected_match'] == True]
        if sel.empty:
            parts.append(g.head(match_limit))
        else:
            rest = g[g['selected_match'] != True]
            combined = pd.concat([sel, rest], ignore_index=True)
            combined = combined.drop_duplicates(keep='first')
            parts.append(combined.head(match_limit))
    return pd.concat(parts, ignore_index=True)


def format_taxa_assignment_info_dataframe(taxa_info: pd.DataFrame, params: dict) -> pd.DataFrame:
    """
    Same column order, filled columns, and sorting as taxa_assignment_INFO.csv in the full workflow.
    Used by create_taxa_assignment_info and taxassign.py.
    """
    taxa_info = taxa_info.copy()
    api_source = params.get('taxonomic_api_source', 'WoRMS').lower()

    if 'cleanedTaxonomy' not in taxa_info.columns:
        if api_source == 'worms':
            from .WoRMS_v3_matching import parse_semicolon_taxonomy
        else:
            from .GBIF_matching import parse_semicolon_taxonomy

        cleaned_taxonomies = [
            ';'.join(parse_semicolon_taxonomy(vid)) if pd.notna(vid) else ''
            for vid in taxa_info['verbatimIdentification']
        ]
        taxa_info['cleanedTaxonomy'] = cleaned_taxonomies

    taxa_info['nameAccordingTo'] = params.get('taxonomic_api_source', 'WoRMS')

    if api_source == 'worms':
        final_column_order = [
            'verbatimIdentification', 'cleanedTaxonomy', 'ambiguous', 'name_change', 'unaccepted_match_row',
            'ranks_matched', 'ranks_provided', 'assignment_score', 'environment',
            'selected_match', 'scientificName', 'taxonRank', 'scientificNameID',
            'higherClassification', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus',
            'match_type_debug', 'nameAccordingTo'
        ]
    else:
        final_column_order = [
            'verbatimIdentification', 'cleanedTaxonomy', 'environment', 'selected_match',
            'scientificName', 'confidence', 'taxonRank', 'taxonID',
            'kingdom', 'phylum', 'class', 'order', 'family', 'genus',
            'match_type_debug', 'nameAccordingTo'
        ]

    for col in final_column_order:
        if col not in taxa_info.columns:
            taxa_info[col] = pd.NA

    taxa_info = taxa_info[final_column_order]

    if api_source == 'gbif':
        taxa_info = taxa_info.sort_values(
            ['verbatimIdentification', 'confidence'], ascending=[True, False]
        ).reset_index(drop=True)
    else:
        taxa_info = taxa_info.sort_values(
            ['verbatimIdentification', 'scientificName']
        ).reset_index(drop=True)

    return taxa_info


def create_taxa_assignment_info(params, reporter):
    """
    Create a taxa_assignment_INFO.csv file.
    For WoRMS, this now includes all ambiguous matches and an 'ambiguous' flag.
    For GBIF, it retains the original behavior of one row per unique verbatimIdentification.
    """
    try:
        reporter.add_section("Creating Taxa Assignment Info File")
        # quiet CLI; report in HTML only
        
        api_source = params.get('taxonomic_api_source', 'WoRMS').lower()
        taxa_info = pd.DataFrame()
        
        # --- Data Loading ---
        if api_source in ['worms', 'gbif'] and 'taxa_info_df' in params:
            # For WoRMS and GBIF, use the detailed info_df created during the matching process
            reporter.add_text(f"Using detailed match data from {api_source.upper()} process for info file.")
            taxa_info = params['taxa_info_df'].copy()
        else:
            # ORIGINAL PATH (Fallback if info_df is missing)
            reporter.add_text(f"Using standard method (one row per unique verbatim ID) for info file. Note: Ambiguous matches may not be shown.")
            occurrence_path = os.path.join(params.get('output_dir', '../processed-v3/'), f'occurrence_core_{api_source}.csv')
            if not os.path.exists(occurrence_path):
                reporter.add_error(f"Taxonomically matched occurrence file not found at {occurrence_path}")
                return
            
            matched_df = pd.read_csv(occurrence_path)
            reporter.add_text(f"Loaded taxonomically matched data: {len(matched_df):,} records")
            
            # This path always produces one row per verbatim ID, so ambiguous is always False.
            taxa_info = matched_df.drop_duplicates(subset=['verbatimIdentification']).copy()
            taxa_info['ambiguous'] = False
            reporter.add_text(f"Found {len(taxa_info):,} unique taxonomy strings")

        taxa_info = format_taxa_assignment_info_dataframe(taxa_info, params)

        # --- Save and Report ---
        output_dir = params.get('output_dir', '../processed-v3/')
        os.makedirs(output_dir, exist_ok=True)
        
        api_choice = params.get('taxonomic_api_source', 'WoRMS')
        output_filename = f"taxa_assignment_INFO_{api_choice}.csv"
        output_path = os.path.join(output_dir, output_filename)
        taxa_info.to_csv(output_path, index=False, na_rep='')
        
        reporter.add_success("Taxa assignment info file created successfully")
        reporter.add_text(f"Saved taxa assignment info: {len(taxa_info):,} total rows ({taxa_info['verbatimIdentification'].nunique():,} unique strings)")
        reporter.add_text(f"Output file: {output_filename}")
        
        # Verify the file was created
        if os.path.exists(output_path):
            file_size = os.path.getsize(output_path) / (1024*1024)  # Size in MB
            reporter.add_text(f"File size: {file_size:.2f} MB")
        else:
            reporter.add_error("Error: File was not created")
            
        # Add detailed explanation and table view to the report
        reporter.add_text("<h3>Detailed Taxa Assignment Information</h3>")
        reporter.add_text(
            "<p>The table below (<code>taxa_assignment_INFO.csv</code>) provides a comprehensive look at the results of the taxonomic matching process. "
            "It includes all potential matches found for each unique <code>verbatimIdentification</code> from your raw data, not just the single best match chosen for the final occurrence file. "
            "This allows for manual review and complete transparency.</p>"
            "<ul>"
            "<li><b>verbatimIdentification:</b> The original, unaltered taxonomic string from your input data.</li>"
            "<li><b>cleanedTaxonomy:</b> A standardized version of the verbatim string used for matching.</li>"
            "<li><b>scientificName:</b> The scientific name of the match returned by the taxonomic service (WoRMS or GBIF).</li>"
            "<li><b>confidence:</b> A score from 0-100 indicating GBIF's confidence in the match (GBIF only).</li>"
            "<li><b>ambiguous:</b> (WoRMS only) A flag indicating if multiple potential matches were found for the verbatim string.</li>"
            "<li><b>name_change:</b> (WoRMS only) True when WoRMS resolved an unaccepted name to its accepted valid name.</li>"
            "<li><b>unaccepted_match_row:</b> (WoRMS only) True when this row is the unaccepted (queried) taxon shown so you can compare its assignment score to the accepted name. Whether it can be selected for the occurrence core is controlled by worms_consider_unaccepted_for_selection in config.</li>"
            "<li><b>selected_match:</b> (WoRMS and GBIF) A flag indicating which of the potential matches was chosen and used in the final occurrence file.</li>"
            "<li><b>higherClassification:</b> (WoRMS only, optional) Pipe-separated higher-taxon lineage from WoRMS, used to preserve intermediate ranks without creating many sparse columns.</li>"
            "<li><b>ranks_matched:</b> (WoRMS only) Count of verbatim lineage tokens found in the matched classification.</li>"
            "<li><b>ranks_provided:</b> (WoRMS only) Total number of verbatim lineage tokens considered for scoring.</li>"
            "<li><b>assignment_score:</b> (WoRMS only) Ratio = ranks_matched / ranks_provided.</li>"
            "<li><b>environment:</b> (WoRMS only) Semicolon-separated WoRMS habitat labels for the selected match, such as marine or freshwater. Blank means WoRMS did not provide habitat data.</li>"
            "<li><b>consistency_check:</b> (GBIF only) A check to ensure the kingdom of the match is consistent with the kingdom in the verbatim string. Helps identify homonym errors.</li>"
            "</ul>"
        )
        # Show preview of the data (limited to 20 rows maximum)
        reporter.add_dataframe(taxa_info, f"Preview of {output_filename} (showing up to 20 rows)", max_rows=20)
        
        # Summary statistics
        ambiguous_count = 0
        if 'ambiguous' in taxa_info.columns:
            ambiguous_count = taxa_info[taxa_info['ambiguous'] == True]['verbatimIdentification'].nunique()
        
        if ambiguous_count > 0:
            reporter.add_text(f"Summary: Found {ambiguous_count:,} unique verbatim strings with ambiguous matches.")
        
    except Exception as e:
        reporter.add_error(f"Taxa assignment info creation failed: {str(e)}")



