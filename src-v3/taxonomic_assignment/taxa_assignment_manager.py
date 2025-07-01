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
                # The script flags these with the specific string 'No WoRMS Match'.
                # These become 'Biota'.
                no_match_mask = matched_df['scientificName'] == 'No WoRMS Match'
                num_no_match = no_match_mask.sum()
                if num_no_match > 0:
                    reporter.add_text(f"Found {num_no_match:,} records with no match in WoRMS. Assigning 'Biota'.")
                    matched_df.loc[no_match_mask, 'scientificName'] = 'Biota'
                    matched_df.loc[no_match_mask, 'scientificNameID'] = 'urn:lsid:marinespecies.org:taxname:1'
                    # Clear all other taxonomic columns for these records
                    for rank_col in ALL_TAX_RANKS:
                        if rank_col in matched_df.columns:
                            matched_df.loc[no_match_mask, rank_col] = np.nan

                # CASE 2: Handle specific high-level Eukaryota assignments that should also be 'Biota'.
                biota_override_mask = (
                    (matched_df['verbatimIdentification'].str.strip() == 'Eukaryota') |
                    (matched_df['verbatimIdentification'].str.strip() == 'Eukaryota;Haptista') |
                    (matched_df['scientificName'] == 'Eukaryota')
                )
                num_eukaryota_override = biota_override_mask.sum()
                if num_eukaryota_override > 0:
                    reporter.add_text(f"Reassigned {num_eukaryota_override:,} high-level Eukaryota records to 'Biota'.")
                    matched_df.loc[biota_override_mask, 'scientificName'] = 'Biota'
                    matched_df.loc[biota_override_mask, 'scientificNameID'] = 'urn:lsid:marinespecies.org:taxname:1'
                    # Clear all other taxonomic columns for these records
                    for rank_col in ALL_TAX_RANKS:
                        if rank_col in matched_df.columns:
                            matched_df.loc[biota_override_mask, rank_col] = np.nan

                # CASE 3: Handle any remaining empty/NaN scientificName records as 'incertae sedis'.
                # This catches inputs that were 'Unassigned' or blank from the start.
                nan_mask = matched_df['scientificName'].isna()
                num_nan = nan_mask.sum()
                if num_nan > 0:
                    reporter.add_text(f"Assigned {num_nan:,} empty records to 'incertae sedis'.")
                    matched_df.loc[nan_mask, 'scientificName'] = 'incertae sedis'
                    matched_df.loc[nan_mask, 'scientificNameID'] = 'urn:lsid:marinespecies.org:taxname:12'
                    # Clear all other taxonomic columns for these records
                    for rank_col in ALL_TAX_RANKS:
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
