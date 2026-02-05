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
        
        # OPTIONAL PR2 database optimization
        use_local_db = params.get('use_local_reference_database', False)
        if use_local_db and api_source == 'WoRMS':
            local_db_path = params.get('local_reference_database_path')
            if local_db_path and os.path.exists(local_db_path):
                try:
                    reporter.add_text(f"Loading local reference database from: {local_db_path}")
                    
                    local_df = pd.read_excel(local_db_path, index_col=None, na_values=[""])
                    
                    # Filter for rows that have a WoRMS ID and a species name
                    local_df.dropna(subset=['worms_id', 'species'], inplace=True)
                    
                    # Clean up the species names to match the format in our main data
                    # (e.g., replacing underscores, removing 'sp.', etc.)
                    species_cleaned = local_df['species'].str.replace('_', ' ', regex=False)
                    species_cleaned = species_cleaned.str.replace(' sp.', '', regex=False).str.strip()
                    
                    # Create the dictionary: {species_name: aphia_id}
                    params['pr2_worms_dict'] = dict(zip(species_cleaned, local_df['worms_id'].astype(int)))
                    
                    reporter.add_success(f"Successfully loaded local reference database with {len(params['pr2_worms_dict'])} AphiaID mappings")
                    reporter.add_text("This will significantly speed up taxonomic matching via direct AphiaID lookup!")
                    
                except Exception as e:
                    reporter.add_warning(f"Could not load local reference database: {e}")
                    params['pr2_worms_dict'] = {}
            else:
                reporter.add_warning(f"Local reference database enabled but file not found at: {local_db_path}")
                params['pr2_worms_dict'] = {}
        else:
            if use_local_db and api_source != 'WoRMS':
                reporter.add_warning("Local reference database optimization is only available for WoRMS API")
            params['pr2_worms_dict'] = {}
        
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
                
                # Define all taxonomic rank columns that need to be managed
                ALL_TAX_RANKS = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'taxonRank']

                # CASE 1: Handle records where WoRMS lookup completely failed.
                # The script now assigns these to 'incertae sedis' directly, so we don't need to change them.
                # But we can log how many we have for reporting purposes.
                no_match_mask = matched_df['match_type_debug'] == 'Failed_All_Stages_NoMatch'
                num_no_match = no_match_mask.sum()
                if num_no_match > 0:
                    reporter.add_text(f"Found {num_no_match:,} records assigned to 'incertae sedis' due to no WoRMS match.")
                
                # Check for pre-handled cases
                pre_handled_mask = (
                    matched_df['match_type_debug'].str.contains('incertae_sedis_simple_case', na=False) |
                    matched_df['match_type_debug'].str.contains('incertae_sedis_unassigned', na=False)
                )
                num_pre_handled = pre_handled_mask.sum()
                if num_pre_handled > 0:
                    reporter.add_text(f"Found {num_pre_handled:,} pre-handled cases (unassigned/empty/simple kingdoms) assigned to 'incertae sedis'.")

                # CASE 2: Handle specific high-level Eukaryota assignments that should be 'incertae sedis'.
                # Only reassign to incertae sedis if the kingdom is specifically 'Eukaryota' and scientificName is also 'Eukaryota'
                eukaryota_override_mask = (
                    (matched_df['kingdom'] == 'Eukaryota') & 
                    (matched_df['scientificName'] == 'Eukaryota') &
                    (~matched_df['match_type_debug'].str.contains('incertae_sedis_simple_case', na=False)) &
                    (~matched_df['match_type_debug'].str.contains('incertae_sedis_unassigned', na=False)))
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
                
                # Define all taxonomic rank columns that need to be managed (GBIF doesn't have scientificNameID)
                ALL_TAX_RANKS_GBIF = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'taxonRank']

                # CASE 1: Handle records where GBIF lookup completely failed.
                no_match_mask = matched_df['match_type_debug'] == 'No_GBIF_Match'
                num_no_match = no_match_mask.sum()
                if num_no_match > 0:
                    reporter.add_text(f"Found {num_no_match:,} records assigned to 'incertae sedis' due to no GBIF match.")
                
                # Check for pre-handled cases
                pre_handled_mask = (
                    matched_df['match_type_debug'].str.contains('incertae_sedis_simple_case', na=False) |
                    matched_df['match_type_debug'].str.contains('incertae_sedis_unassigned', na=False)
                )
                num_pre_handled = pre_handled_mask.sum()
                if num_pre_handled > 0:
                    reporter.add_text(f"Found {num_pre_handled:,} pre-handled cases (unassigned/empty/simple kingdoms) assigned to 'incertae sedis'.")

                # CASE 2: Handle specific high-level Eukaryota assignments that should be 'incertae sedis'.
                # Only reassign to incertae sedis if the kingdom is specifically 'Eukaryota' and scientificName is also 'Eukaryota'
                eukaryota_override_mask = (
                    (matched_df['kingdom'] == 'Eukaryota') & 
                    (matched_df['scientificName'] == 'Eukaryota') &
                    (~matched_df['match_type_debug'].str.contains('incertae_sedis_simple_case', na=False)) &
                    (~matched_df['match_type_debug'].str.contains('incertae_sedis_unassigned', na=False)))
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

        # --- Column Formatting ---
        if 'cleanedTaxonomy' not in taxa_info.columns:
            # Recreate cleanedTaxonomy if it's missing (e.g., from old GBIF path)
            if api_source == 'worms': from .WoRMS_v3_matching import parse_semicolon_taxonomy
            else: from .GBIF_matching import parse_semicolon_taxonomy
            
            cleaned_taxonomies = [';'.join(parse_semicolon_taxonomy(vid)) if pd.notna(vid) else '' for vid in taxa_info['verbatimIdentification']]
            taxa_info['cleanedTaxonomy'] = cleaned_taxonomies
            
        taxa_info['nameAccordingTo'] = params.get('taxonomic_api_source', 'WoRMS')

        # Define final column order, now including selected_match and consistency_check
        if api_source == 'worms':
            final_column_order = [
                'verbatimIdentification', 'cleanedTaxonomy', 'ambiguous', 'name_change', 'selected_match',
                'scientificName', 'taxonRank', 'scientificNameID', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus',
                'match_type_debug', 'nameAccordingTo'
            ]
        else: # GBIF, includes confidence, does not include ambiguous
            final_column_order = [
                'verbatimIdentification', 'cleanedTaxonomy', 'selected_match',
                'scientificName', 'confidence', 'taxonRank', 'taxonID', 
                'kingdom', 'phylum', 'class', 'order', 'family', 'genus',
                'match_type_debug', 'nameAccordingTo'
            ]
        
        # Ensure all required columns exist, adding any that are missing
        for col in final_column_order:
            if col not in taxa_info.columns:
                taxa_info[col] = pd.NA
        
        # Filter to only include the desired columns in the correct order
        taxa_info = taxa_info[final_column_order]
        
        # Sort for consistency and readability
        if api_source == 'gbif':
            # For GBIF, group by the original string, then show the best confidence match first
            taxa_info = taxa_info.sort_values(['verbatimIdentification', 'confidence'], ascending=[True, False]).reset_index(drop=True)
        else:
            # For WoRMS, group by original string, then by the scientific name
            taxa_info = taxa_info.sort_values(['verbatimIdentification', 'scientificName']).reset_index(drop=True)
        
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
            "<li><b>selected_match:</b> (WoRMS and GBIF) A flag indicating which of the potential matches was chosen and used in the final occurrence file.</li>"
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



