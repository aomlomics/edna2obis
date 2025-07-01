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
from contextlib import redirect_stdout, redirect_stderr

# Import the API-specific matching scripts
from . import WoRMS_v3_matching
from . import GBIF_matching


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
                    reporter.add_text(f"🧬 Loading local reference database from: {local_db_path}")
                    print(f"🧬 Loading local reference database from: {local_db_path}")
                    
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
                    print(f"✅ Successfully loaded local reference database with {len(params['pr2_worms_dict'])} AphiaID mappings")
                    
                except Exception as e:
                    reporter.add_warning(f"Could not load local reference database: {e}")
                    params['pr2_worms_dict'] = {}
                    print(f"⚠️ Could not load local reference database: {e}")
            else:
                reporter.add_warning(f"Local reference database enabled but file not found at: {local_db_path}")
                params['pr2_worms_dict'] = {}
                print(f"⚠️ Local reference database enabled but file not found at: {local_db_path}")
        else:
            if use_local_db and api_source != 'WoRMS':
                reporter.add_warning("Local reference database optimization is only available for WoRMS API")
            params['pr2_worms_dict'] = {}
        
        # Determine maximum taxonomic ranks for each assay
        reporter.add_text("Determining maximum taxonomic ranks for each assay...")
        print("Determining maximum taxonomic ranks for each assay...")
        
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
                            reporter.add_text(f"  - Assay '{assay_name}': {len(rank_cols)} rank columns")
                            print(f"  - Assay '{assay_name}': {len(rank_cols)} rank columns")
                        except (KeyError, ValueError):
                            reporter.add_text(f"  - Assay '{assay_name}': Using default 7 ranks (columns not found)")
                            assay_rank_info[assay_name] = {'max_depth': 7}
        
        # Fallback for any missing assays
        for assay in df_to_match['assay_name'].unique():
            if assay not in assay_rank_info:
                assay_rank_info[assay] = {'max_depth': 7}
                reporter.add_text(f"  - Assay '{assay}': Using default 7 ranks")
        
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
        logger.addHandler(ch)
        logger.setLevel(logging.INFO)
        
        # Also capture stdout/stderr for any print statements
        captured_output = io.StringIO()
        
        reporter.add_text("Starting taxonomic matching process...")
        print(f"🧬 Starting taxonomic assignment using {api_source} API...")
        
        try:
            # Run the matching with both logging and stdout/stderr capture
            if api_source == 'WoRMS':
                print("Running WoRMS API matching...")
                with redirect_stdout(captured_output), redirect_stderr(captured_output):
                    matched_df = WoRMS_v3_matching.get_worms_match_for_dataframe(
                        occurrence_df=df_to_match,
                        params_dict=params,
                        n_proc=params.get('worms_n_proc', 0)
                    )
            elif api_source == 'GBIF':
                print("Running GBIF API matching...")
                with redirect_stdout(captured_output), redirect_stderr(captured_output):
                    matched_df = GBIF_matching.get_gbif_match_for_dataframe(
                        occurrence_df=df_to_match,
                        params_dict=params,
                        n_proc=params.get('gbif_n_proc', 0)
                    )
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
                print("\n🐟 Taxonomic matching progress:")
                print(progress_text)
                # Add to HTML report
                reporter.add_text("Taxonomic matching progress:")
                reporter.add_text(f"<pre>{progress_text}</pre>")
            else:
                print("📊 Taxonomic matching completed (no detailed progress captured)")
                reporter.add_text("Taxonomic matching completed")
            
        except Exception as e:
            # If there's an error, still get any captured output
            stdout_output = captured_output.getvalue()
            log_output = log_capture_string.getvalue()
            if stdout_output or log_output:
                print("Progress before error:")
                if stdout_output:
                    print("Console:", stdout_output)
                if log_output:
                    print("Logs:", log_output)
                reporter.add_text(f"Progress before error:<pre>{stdout_output}\n{log_output}</pre>")
            raise e
        finally:
            # Clean up the logger
            logger.removeHandler(ch)
        
        if matched_df is not None and not matched_df.empty:
            # Apply post-processing for WoRMS (the manual corrections)
            if api_source == 'WoRMS':
                reporter.add_text("Applying manual taxonomic corrections...")
                print("🔧 Applying manual taxonomic corrections...")
                
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
                # Note: Simple "Eukaryota" cases are now handled earlier, but we may have "Eukaryota;Haptista" style cases
                eukaryota_override_mask = (
                    (matched_df['verbatimIdentification'].str.strip() == 'Eukaryota;Haptista') |
                    ((matched_df['scientificName'] == 'Eukaryota') & 
                     (~matched_df['match_type_debug'].str.contains('incertae_sedis_simple_case', na=False)) &
                     (~matched_df['match_type_debug'].str.contains('incertae_sedis_unassigned', na=False)))
                )
                num_eukaryota_override = eukaryota_override_mask.sum()
                if num_eukaryota_override > 0:
                    reporter.add_text(f"Reassigned {num_eukaryota_override:,} complex Eukaryota records to 'incertae sedis'.")
                    matched_df.loc[eukaryota_override_mask, 'scientificName'] = 'incertae sedis'
                    matched_df.loc[eukaryota_override_mask, 'scientificNameID'] = 'urn:lsid:marinespecies.org:taxname:12'
                    # Clear all other taxonomic columns for these records
                    for rank_col in ALL_TAX_RANKS:
                        if rank_col in matched_df.columns:
                            matched_df.loc[eukaryota_override_mask, rank_col] = np.nan

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
                            matched_df.loc[nan_mask, rank_col] = np.nan
            
            elif api_source == 'GBIF':
                reporter.add_text("Applying manual taxonomic corrections for GBIF...")
                print("🔧 Applying manual taxonomic corrections for GBIF...")
                
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
                eukaryota_override_mask = (
                    (matched_df['verbatimIdentification'].str.strip() == 'Eukaryota;Haptista') |
                    ((matched_df['scientificName'] == 'Eukaryota') & 
                     (~matched_df['match_type_debug'].str.contains('incertae_sedis_simple_case', na=False)) &
                     (~matched_df['match_type_debug'].str.contains('incertae_sedis_unassigned', na=False)))
                )
                num_eukaryota_override = eukaryota_override_mask.sum()
                if num_eukaryota_override > 0:
                    reporter.add_text(f"Reassigned {num_eukaryota_override:,} complex Eukaryota records to 'incertae sedis'.")
                    matched_df.loc[eukaryota_override_mask, 'scientificName'] = 'incertae sedis'
                    # No scientificNameID for GBIF
                    # Clear all other taxonomic columns for these records
                    for rank_col in ALL_TAX_RANKS_GBIF:
                        if rank_col in matched_df.columns:
                            matched_df.loc[eukaryota_override_mask, rank_col] = np.nan

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
                            matched_df.loc[nan_mask, rank_col] = np.nan
            
            # Post-matching processing and final save
            reporter.add_text("Starting post-matching processing.")
            
            # Remove temporary columns and save
            columns_to_drop = ['assay_name']
            if api_source == 'GBIF' and 'scientificNameID' in matched_df.columns:
                columns_to_drop.append('scientificNameID')
                reporter.add_text("Adjusted final columns for GBIF standard: Removed 'scientificNameID'.")
            
            matched_df = matched_df.drop(columns=columns_to_drop, errors='ignore')
            reporter.add_text(f"Removed temporary columns: {columns_to_drop}")
            
            # Define the desired final column order (same as in occurrence_builder.py but without assay_name)
            # NOTE: cleanedTaxonomy is excluded here - it should only appear in taxa_assignment_INFO.csv
            DESIRED_FINAL_COLUMNS_IN_ORDER = [
                'eventID', 'organismQuantity', 'occurrenceID', 'verbatimIdentification',
                'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 
                'scientificName', 'scientificNameID', 'match_type_debug', 
                'taxonRank', 'identificationRemarks',
                'taxonID', 'basisOfRecord', 'nameAccordingTo', 'organismQuantityType',
                'recordedBy', 'materialSampleID', 'sampleSizeValue', 'sampleSizeUnit',
                'associatedSequences', 'locationID', 'eventDate', 'minimumDepthInMeters', 'maximumDepthInMeters',
                'locality', 'decimalLatitude', 'decimalLongitude',
                'geodeticDatum', 'parentEventID', 'datasetID', 'occurrenceStatus'
            ]
            
            # For GBIF, remove scientificNameID from the desired column order
            if api_source == 'GBIF':
                DESIRED_FINAL_COLUMNS_IN_ORDER = [col for col in DESIRED_FINAL_COLUMNS_IN_ORDER if col != 'scientificNameID']
            
            # Ensure all desired columns exist in the DataFrame, add as NA if missing
            for col in DESIRED_FINAL_COLUMNS_IN_ORDER:
                if col not in matched_df.columns:
                    matched_df[col] = pd.NA
            
            # Reorder columns to match the desired order
            matched_df = matched_df.reindex(columns=DESIRED_FINAL_COLUMNS_IN_ORDER)
            reporter.add_text("Applied proper column ordering to match WoRMS output format.")
            
            # Save the taxonomically matched file
            output_filename = f"occurrence_{api_source.lower()}_matched.csv"
            output_dir = params.get('output_dir', '../processed-v3/')
            
            # Ensure output directory exists
            os.makedirs(output_dir, exist_ok=True)
            
            output_path = os.path.join(output_dir, output_filename)
            matched_df.to_csv(output_path, index=False, na_rep='')
            
            reporter.add_success(f"Taxonomic assignment completed successfully")
            reporter.add_text(f"Saved taxonomically matched data: {len(matched_df):,} records")
            reporter.add_text(f"Output file: {output_filename}")
            
            print(f"✅ Taxonomic assignment completed! Saved {len(matched_df):,} records to {output_filename}")
            
            # Add summary statistics
            unique_taxa = matched_df['scientificName'].nunique()
            reporter.add_text(f"Summary: {unique_taxa:,} unique taxa identified")
            
        else:
            reporter.add_error("Taxonomic assignment returned empty results")
            
    except Exception as e:
        reporter.add_error(f"Taxonomic assignment failed: {str(e)}")
        print(f"❌ Taxonomic assignment failed: {str(e)}")
        import traceback
        traceback.print_exc()


def create_taxa_assignment_info(params, reporter):
    """
    Create a taxa_assignment_INFO.csv file with one row per unique taxonomy string
    showing the cleaned taxonomy and API results for each unique verbatimIdentification
    """
    try:
        reporter.add_section("Creating Taxa Assignment Info File")
        print("📊 Creating taxa assignment info file...")
        
        api_source = params.get('taxonomic_api_source', 'WoRMS').lower()
        
        # Load the taxonomically matched occurrence file
        occurrence_path = os.path.join(params.get('output_dir', '../processed-v3/'), f'occurrence_{api_source}_matched.csv')
        if not os.path.exists(occurrence_path):
            reporter.add_error(f"Taxonomically matched occurrence file not found at {occurrence_path}")
            return
            
        matched_df = pd.read_csv(occurrence_path)
        reporter.add_text(f"Loaded taxonomically matched data: {len(matched_df):,} records")
        
        # Get unique taxonomy strings (one row per unique verbatimIdentification)
        unique_taxa = matched_df.drop_duplicates(subset=['verbatimIdentification']).copy()
        reporter.add_text(f"Found {len(unique_taxa):,} unique taxonomy strings")
        
        # Define columns based on API source
        if api_source == 'worms':
            taxonomic_id_col = 'scientificNameID'
            columns_to_keep = [
                'verbatimIdentification', 
                'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 
                'scientificName', 'scientificNameID', 'taxonRank', 
                'match_type_debug'
            ]
        else:  # GBIF
            taxonomic_id_col = 'taxonID'
            columns_to_keep = [
                'verbatimIdentification',
                'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 
                'scientificName', 'taxonID', 'taxonRank', 
                'match_type_debug'
            ]
        
        # Select and reorder columns (excluding cleanedTaxonomy since it's not in the final file)
        columns_to_keep_available = [col for col in columns_to_keep if col != 'cleanedTaxonomy']
        taxa_info = unique_taxa[columns_to_keep_available].copy()
        
        # Recalculate cleanedTaxonomy from verbatimIdentification
        # Import the parsing function from the appropriate matching script
        if api_source == 'worms':
            from .WoRMS_v3_matching import parse_semicolon_taxonomy
        else:  # GBIF
            from .GBIF_matching import parse_semicolon_taxonomy
        
        # Calculate cleanedTaxonomy for each verbatimIdentification
        cleaned_taxonomies = []
        for verbatim_id in taxa_info['verbatimIdentification']:
            parsed_names = parse_semicolon_taxonomy(verbatim_id)
            cleaned_taxonomy = ';'.join(parsed_names) if parsed_names else str(verbatim_id) if verbatim_id else ''
            cleaned_taxonomies.append(cleaned_taxonomy)
        
        taxa_info['cleanedTaxonomy'] = cleaned_taxonomies
        
        # Add nameAccordingTo from config (same value for all rows)
        taxa_info['nameAccordingTo'] = params.get('taxonomic_api_source', 'WoRMS')
        
        # Reorder columns to put cleanedTaxonomy right after verbatimIdentification
        if api_source == 'worms':
            final_column_order = [
                'verbatimIdentification', 'cleanedTaxonomy',
                'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 
                'scientificName', 'scientificNameID', 'taxonRank', 
                'match_type_debug', 'nameAccordingTo'
            ]
        else:  # GBIF
            final_column_order = [
                'verbatimIdentification', 'cleanedTaxonomy',
                'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 
                'scientificName', 'taxonID', 'taxonRank', 
                'match_type_debug', 'nameAccordingTo'
            ]
        
        taxa_info = taxa_info[final_column_order]
        
        # Sort by verbatimIdentification for easier review
        taxa_info = taxa_info.sort_values('verbatimIdentification')
        
        # Save the taxa assignment info file
        output_dir = params.get('output_dir', '../processed-v3/')
        os.makedirs(output_dir, exist_ok=True)
        
        output_filename = "taxa_assignment_INFO.csv"
        output_path = os.path.join(output_dir, output_filename)
        taxa_info.to_csv(output_path, index=False, na_rep='')
        
        reporter.add_success("Taxa assignment info file created successfully")
        reporter.add_text(f"Saved taxa assignment info: {len(taxa_info):,} unique taxonomy strings")
        reporter.add_text(f"Output file: {output_filename}")
        
        # Verify the file was created
        if os.path.exists(output_path):
            file_size = os.path.getsize(output_path) / (1024*1024)  # Size in MB
            reporter.add_text(f"File size: {file_size:.2f} MB")
            print(f"✅ Taxa assignment info created! Saved {len(taxa_info):,} unique taxonomy strings to {output_filename}")
        else:
            reporter.add_error("❌ Error: File was not created")
            
        # Show preview of the data
        reporter.add_text("<h4>Taxa Assignment Info Preview:</h4>")
        reporter.add_dataframe(taxa_info.head(10), "First 10 entries from taxa_assignment_INFO.csv")
        
        # Summary statistics
        empty_cleaned = taxa_info['cleanedTaxonomy'].isna().sum()
        reporter.add_text(f"Summary: {len(taxa_info) - empty_cleaned:,} taxa with cleaned taxonomy, {empty_cleaned:,} empty/unassigned")
        
    except Exception as e:
        reporter.add_error(f"Taxa assignment info creation failed: {str(e)}")
        print(f"❌ Taxa assignment info creation failed: {str(e)}")
        import traceback
        traceback.print_exc()



