"""
Occurrence Core Builder for edna2obis
Contains the complex logic for building occurrence cores from ASV data
Build a combined Occurrence Core for all Analyses
"""

import pandas as pd
import numpy as np
import os
import warnings
import traceback
import datetime


def create_occurrence_core(data, raw_data_tables, params, dwc_data, reporter):
    """
    Build a combined Occurrence Core for all Analyses
    This operation is complex, includes merging dataframes together, adding the missing fields 
    which OBIS and GBIF expect, and parsing taxonomy raw data.
    
    EXACT CODE FROM NOTEBOOK CELL 26
    """
    
    reporter.add_section("Building Combined Occurrence Core", level=2)
    
    try:
        # --- MAIN DATA PROCESSING LOOP --- Loop through each analysis run defined in params['datafiles']
        warnings.simplefilter(action='ignore', category=FutureWarning)

        reporter.add_text(f"🚀 Starting data processing for {len(params['datafiles'])} analysis run(s) to generate occurrence records.")

        # --- NEW: Load mappings and determine desired columns from the mapper ---
        import yaml
        with open('data_mapper.yaml', 'r', encoding='utf-8') as f:
            full_mapper = yaml.safe_load(f)
        
        format_prefix = "generic_" if params.get('metadata_format') == 'GENERIC' else ""
        occurrence_key = f"{format_prefix}occurrence_core"
        occurrence_map_config = full_mapper.get(occurrence_key, {})
        DESIRED_OCCURRENCE_COLUMNS_IN_ORDER = list(occurrence_map_config.keys())

        # --- FIX: Ensure internal processing columns are always present ---
        # These are needed for downstream steps but are not part of the final DwC output
        # and therefore don't belong in the mapper.
        internal_processing_cols = ['assay_name', 'match_type_debug']
        for col in internal_processing_cols:
            if col not in DESIRED_OCCURRENCE_COLUMNS_IN_ORDER:
                DESIRED_OCCURRENCE_COLUMNS_IN_ORDER.append(col)

        output_dir = params.get('output_dir', "processed-v3/")
        os.makedirs(output_dir, exist_ok=True)
        output_filename = "occurrence.csv"
        output_path = os.path.join(output_dir, output_filename)

        all_processed_occurrence_dfs = []
        successful_runs = 0
        failed_runs = 0

        project_recorded_by = "recordedBy_NotProvided"
        project_dataset_id = "DatasetID_NotProvided"

        def get_project_meta_value(project_meta_df, term_to_find, default_val=pd.NA):
            if not all(col in project_meta_df.columns for col in ['term_name', 'project_level']):
                return default_val
            term_to_find_stripped = str(term_to_find).strip()
            match = project_meta_df[project_meta_df['term_name'].astype(str).str.strip().str.lower() == term_to_find_stripped.lower()]
            if not match.empty:
                value = match['project_level'].iloc[0]
                if pd.notna(value):
                    return str(value).strip()
            return default_val

        if 'projectMetadata' in data and not data['projectMetadata'].empty:
            project_meta_df = data['projectMetadata']
            project_recorded_by = get_project_meta_value(project_meta_df, 'recordedBy', project_recorded_by)
            project_dataset_id = get_project_meta_value(project_meta_df, 'project_id', project_dataset_id)
        else:
            reporter.add_text("  ⚠️ Warning: projectMetadata is empty or not found. 'recordedBy' and 'datasetID' will use default placeholder values.")

        # --- MAIN DATA PROCESSING LOOP ---
        for analysis_run_name, file_paths_dict in params['datafiles'].items():
            reporter.add_text(f"Processing Analysis Run: {analysis_run_name}")
            try:
                # --- STEP 1: Load and Prepare Raw Taxonomy and Abundance Data ---
                if not (analysis_run_name in raw_data_tables and
                        'taxonomy' in raw_data_tables[analysis_run_name] and not raw_data_tables[analysis_run_name]['taxonomy'].empty and
                        'occurrence' in raw_data_tables[analysis_run_name] and not raw_data_tables[analysis_run_name]['occurrence'].empty):
                    reporter.add_text(f"  Skipping {analysis_run_name}: Raw taxonomy or occurrence data is missing or empty.")
                    failed_runs += 1
                    continue

                current_tax_df_raw = raw_data_tables[analysis_run_name]['taxonomy'].copy()
                current_abundance_df_raw = raw_data_tables[analysis_run_name]['occurrence'].copy()

                featureid_col_tax = current_tax_df_raw.columns[0]
                current_tax_df_raw.rename(columns={featureid_col_tax: 'featureid'}, inplace=True)
                featureid_col_abun = current_abundance_df_raw.columns[0] 
                current_abundance_df_raw.rename(columns={featureid_col_abun: 'featureid'}, inplace=True)

                # --- NEW: Dynamically identify taxonomy columns to keep from the mapper ---
                tax_cols_to_keep = {'featureid'} # Always keep featureid
                rename_map_for_tax = {}

                for dwc_term, mapping_info in occurrence_map_config.items():
                    if isinstance(mapping_info, dict) and mapping_info.get('source') == 'taxonomy':
                        faire_term = mapping_info.get('faire_term')
                        if faire_term in current_tax_df_raw.columns:
                            tax_cols_to_keep.add(faire_term)
                            if faire_term != dwc_term:
                                rename_map_for_tax[faire_term] = dwc_term
                
                # Special handling for confidence, which isn't a DwC term but is needed for remarks
                if 'Confidence' in current_tax_df_raw.columns:
                    tax_cols_to_keep.add('Confidence')
                else:
                    current_tax_df_raw['Confidence'] = pd.NA

                current_tax_df_processed = current_tax_df_raw[list(tax_cols_to_keep)].copy()
                current_tax_df_processed.rename(columns=rename_map_for_tax, inplace=True)

                # --- STEP 2: Melt Abundance and Merge with Taxonomy ---
                current_assay_occ_melted = pd.melt(
                    current_abundance_df_raw, id_vars=['featureid'],
                    var_name='samp_name', value_name='organismQuantity'
                )
                current_assay_occ_melted = current_assay_occ_melted[current_assay_occ_melted['organismQuantity'] > 0.0]
                
                current_assay_occurrence_intermediate_df = pd.merge(
                    current_assay_occ_melted, current_tax_df_processed,
                    on='featureid', how='left'
                )

                # --- STEP 3: Initialize ALL DwC Fields from the mapper ---
                for col in DESIRED_OCCURRENCE_COLUMNS_IN_ORDER:
                    if col not in current_assay_occurrence_intermediate_df.columns:
                         current_assay_occurrence_intermediate_df[col] = pd.NA
                
                # Set values for fields that are constructed or have fixed values for this process
                current_assay_occurrence_intermediate_df['taxonID'] = 'ASV:' + current_assay_occurrence_intermediate_df['featureid'].astype(str)
                current_assay_occurrence_intermediate_df['organismQuantityType'] = 'DNA sequence reads'
                current_assay_occurrence_intermediate_df['occurrenceStatus'] = 'present'
                current_assay_occurrence_intermediate_df['basisOfRecord'] = 'MaterialSample'
                current_assay_occurrence_intermediate_df['nameAccordingTo'] = 'Original Classification; WoRMS/GBIF (pending further matching)'
                
                # --- STEP 4: Assign Project-Level Metadata ---
                current_assay_occurrence_intermediate_df['recordedBy'] = project_recorded_by
                current_assay_occurrence_intermediate_df['datasetID'] = project_dataset_id

                # --- STEP 5: Merge `sampleMetadata` and map FAIRe terms to DwC terms ---
                if 'sampleMetadata' in data and not data['sampleMetadata'].empty:
                    sm_df_to_merge = data['sampleMetadata'].copy()
                    sm_df_to_merge['samp_name'] = sm_df_to_merge['samp_name'].astype(str).str.strip()
                    current_assay_occurrence_intermediate_df['samp_name'] = current_assay_occurrence_intermediate_df['samp_name'].astype(str).str.strip()

                    current_assay_occurrence_intermediate_df = pd.merge(
                        current_assay_occurrence_intermediate_df, sm_df_to_merge,
                        on='samp_name', how='left', suffixes=('', '_sm') 
                    )
                    
                    # --- NEW: Populate DwC columns from sampleMetadata based on the mapper ---
                    for dwc_col_target, mapping_info in occurrence_map_config.items():
                        if isinstance(mapping_info, dict) and mapping_info.get('source') == 'sampleMetadata':
                            faire_col_source_original = str(mapping_info.get('faire_term')).strip()
                            
                            source_col_in_df = None
                            if faire_col_source_original + '_sm' in current_assay_occurrence_intermediate_df.columns:
                                source_col_in_df = faire_col_source_original + '_sm'
                            elif faire_col_source_original in current_assay_occurrence_intermediate_df.columns:
                                source_col_in_df = faire_col_source_original
                            
                            if source_col_in_df:
                                current_assay_occurrence_intermediate_df[dwc_col_target] = current_assay_occurrence_intermediate_df[source_col_in_df]
                            elif dwc_col_target in ['locality', 'decimalLatitude', 'decimalLongitude', 'geodeticDatum', 'eventDate']:
                                 reporter.add_text(f"  DIAGNOSTIC: For DwC term '{dwc_col_target}', its mapped FAIRe term '{faire_col_source_original}' (from mapper) was NOT found as a column in the merged sample data. The column will be empty if not populated by other means.")
                else:
                    reporter.add_text(f"  Warning: 'sampleMetadata' is empty or not found. Cannot merge for DwC term population for run {analysis_run_name}.")

                # --- Construct 'locationID' with new flexible logic ---
                
                # Priority 1: Use existing 'locationID' column if it's already in the data
                if 'locationID' in current_assay_occurrence_intermediate_df.columns and not current_assay_occurrence_intermediate_df['locationID'].isna().all():
                    reporter.add_text("Found existing 'locationID' column in source data. Using it directly.")
                
                else:
                    # Priority 2: Use components from config.yaml
                    id_components = params.get('locationID_components', [])
                    
                    # Check if all specified component columns exist in the dataframe
                    missing_components = [c for c in id_components if c not in current_assay_occurrence_intermediate_df.columns]
                    
                    if id_components and not missing_components:
                        reporter.add_text(f"Constructing 'locationID' from config components: {', '.join(id_components)}")
                        
                        # Convert all component columns to string and combine
                        component_series = [current_assay_occurrence_intermediate_df[col].astype(str).fillna(f"No_{col}") for col in id_components]
                        current_assay_occurrence_intermediate_df['locationID'] = pd.Series('_'.join(map(str, t)) for t in zip(*component_series))
                    
                    else:
                        # Priority 3: Fallback to default (and warn if component columns were specified but not found)
                        if missing_components:
                            reporter.add_warning(f"Could not construct 'locationID' from config. The following columns were not found: {', '.join(missing_components)}. Falling back to default.")
                        
                        # Original default logic
                        default_components = ['line_id', 'station_id']
                        if all(c in current_assay_occurrence_intermediate_df.columns for c in default_components):
                            line_ids = current_assay_occurrence_intermediate_df['line_id'].astype(str).fillna('NoLineID')
                            station_ids = current_assay_occurrence_intermediate_df['station_id'].astype(str).fillna('NoStationID')
                            current_assay_occurrence_intermediate_df['locationID'] = line_ids + "_" + station_ids
                        else:
                            # If even the default components can't be found, set to NA
                             if 'locationID' not in current_assay_occurrence_intermediate_df.columns:
                                  current_assay_occurrence_intermediate_df['locationID'] = pd.NA


                # --- STEP 6: Merge `experimentRunMetadata` & Define `eventID`, `associatedSequences` ---
                
                # Path 1 (FAIRe-NOAA): Get assay name from the pre-processed analysisMetadata sheets.
                assay_name_for_current_run = next((an_key for an_key, runs_dict_val in data.get('analysis_data_by_assay', {}).items() if isinstance(runs_dict_val, dict) and analysis_run_name in runs_dict_val), None)
                
                # Path 2 (GENERIC FAIRe Fallback): If not found from analysis sheets, get it from the config.
                if not assay_name_for_current_run:
                    assay_name_for_current_run = params.get('datafiles', {}).get(analysis_run_name, {}).get('assay_name')
                    if assay_name_for_current_run:
                        reporter.add_text(f"  INFO: Using assay_name '{assay_name_for_current_run}' from config.yaml for run '{analysis_run_name}'.")
                
                # --- FINAL ASSIGNMENT ---
                # This ensures assay_name is assigned from the config for BOTH modes if available, providing a consistent override.
                config_assay_name = params.get('datafiles', {}).get(analysis_run_name, {}).get('assay_name')
                if config_assay_name:
                    assay_name_for_current_run = config_assay_name
                
                current_assay_occurrence_intermediate_df['assay_name'] = assay_name_for_current_run
                
                final_assay_name = assay_name_for_current_run # Use this variable consistently hereafter

                if not final_assay_name:
                    reporter.add_text(f"    ERROR: Could not determine assay_name for '{analysis_run_name}' from analysisMetadata sheets or config.yaml.")
                    current_assay_occurrence_intermediate_df['eventID'] = current_assay_occurrence_intermediate_df['eventID'].fillna(f"ERROR_eventID_for_{analysis_run_name}")
                    # Also ensure a placeholder for assay_name if it couldn't be found, though it should be an error condition
                    current_assay_occurrence_intermediate_df['assay_name'] = current_assay_occurrence_intermediate_df['assay_name'].fillna(f"UNKNOWN_ASSAY_FOR_{analysis_run_name}")
                elif 'experimentRunMetadata' in data and not data['experimentRunMetadata'].empty:
                    erm_df = data['experimentRunMetadata'].copy()
                    erm_df['samp_name'] = erm_df['samp_name'].astype(str).str.strip()
                    erm_subset = erm_df[erm_df['samp_name'].isin(current_assay_occurrence_intermediate_df['samp_name'].astype(str).str.strip())].copy()
                    if final_assay_name and final_assay_name.startswith('UNKNOWN_ASSAY_FOR_') is False and 'assay_name' in erm_subset.columns:
                        # Further filter the subset for the specific assay of the current run
                        erm_subset_assay_specific = erm_subset[erm_subset['assay_name'].astype(str).str.strip() == str(final_assay_name).strip()]
                        # If filtering results in an empty dataframe, it means no ERM records for this assay, but we can still proceed with the broader erm_subset for other potential merges.
                        if not erm_subset_assay_specific.empty:
                            erm_subset = erm_subset_assay_specific

                    # --- NEW: Get terms to merge from experimentRunMetadata from the mapper ---
                    cols_to_select_from_erm = {'samp_name'}
                    faire_terms_from_erm = {}
                    for dwc_term, mapping_info in occurrence_map_config.items():
                        if isinstance(mapping_info, dict) and mapping_info.get('source') == 'experimentRunMetadata':
                            faire_term = mapping_info.get('faire_term')
                            if faire_term in erm_subset.columns:
                                cols_to_select_from_erm.add(faire_term)
                                faire_terms_from_erm[dwc_term] = faire_term

                    if cols_to_select_from_erm and len(cols_to_select_from_erm) > 1:
                        erm_to_merge = erm_subset[list(cols_to_select_from_erm)].drop_duplicates(subset=['samp_name']).copy()
                        current_assay_occurrence_intermediate_df = pd.merge(
                            current_assay_occurrence_intermediate_df, erm_to_merge,
                            on='samp_name', how='left', suffixes=('', '_erm')
                        )

                        for dwc_term, faire_term in faire_terms_from_erm.items():
                            source_col_actual = faire_term + '_erm' if faire_term + '_erm' in current_assay_occurrence_intermediate_df.columns else faire_term
                            if source_col_actual in current_assay_occurrence_intermediate_df.columns:
                                # For eventID specifically, we fillna. For others, we might overwrite. This logic may need refinement.
                                # For now, let's assume fillna is safe for all erm-sourced columns.
                                current_assay_occurrence_intermediate_df[dwc_term] = current_assay_occurrence_intermediate_df[dwc_term].fillna(current_assay_occurrence_intermediate_df[source_col_actual])
                else:
                    # Keep placeholder if ERM not available
                    current_assay_occurrence_intermediate_df['eventID'] = current_assay_occurrence_intermediate_df['eventID'].fillna(f"ERROR_eventID_for_{analysis_run_name}")
                
                # Ensure eventID is a string for the occurrenceID construction
                if 'eventID' in current_assay_occurrence_intermediate_df.columns:
                    current_assay_occurrence_intermediate_df['eventID'] = current_assay_occurrence_intermediate_df['eventID'].fillna(f"NO_EVENT_ID_FOR_{analysis_run_name}").astype(str)
                else:
                    current_assay_occurrence_intermediate_df['eventID'] = f"NO_EVENT_ID_FOR_{analysis_run_name}"
                
                # --- STEP 7: Construct `occurrenceID` ---
                current_assay_occurrence_intermediate_df['occurrenceID'] = current_assay_occurrence_intermediate_df['eventID'] + '_occ_' + current_assay_occurrence_intermediate_df['featureid'].astype(str)

                # --- STEP 8: `identificationRemarks`, `sampleSizeValue`/`Unit`, `parentEventID` ---
                otu_seq_comp_appr_str = "Unknown sequence comparison approach"
                taxa_class_method_str = ""
                taxa_ref_db_str = "Unknown reference DB"

                if final_assay_name and analysis_run_name in data.get('analysis_data_by_assay', {}).get(final_assay_name, {}):
                    analysis_meta_df_for_run = data['analysis_data_by_assay'][final_assay_name][analysis_run_name]
                    if 'term_name' in analysis_meta_df_for_run.columns and 'values' in analysis_meta_df_for_run.columns:
                        def get_analysis_meta(term, df, default):
                            val_series = df[df['term_name'].astype(str).str.strip() == term]['values']
                            return str(val_series.iloc[0]).strip() if not val_series.empty and pd.notna(val_series.iloc[0]) else default
                        otu_seq_comp_appr_str = get_analysis_meta('otu_seq_comp_appr', analysis_meta_df_for_run, otu_seq_comp_appr_str)
                        taxa_class_method_str = get_analysis_meta('taxa_class_method', analysis_meta_df_for_run, taxa_class_method_str)
                        taxa_ref_db_str = get_analysis_meta('otu_db', analysis_meta_df_for_run, taxa_ref_db_str)
                
                confidence_value_series = pd.Series(["unknown confidence"] * len(current_assay_occurrence_intermediate_df), index=current_assay_occurrence_intermediate_df.index, dtype=object)
                if 'Confidence' in current_assay_occurrence_intermediate_df.columns:
                    confidence_value_series = current_assay_occurrence_intermediate_df['Confidence'].astype(str).fillna("unknown confidence")
                
                current_assay_occurrence_intermediate_df['identificationRemarks'] = f"{otu_seq_comp_appr_str}, confidence: " + confidence_value_series + f", against reference database: {taxa_ref_db_str}"
                
                if 'eventID' in current_assay_occurrence_intermediate_df.columns and not current_assay_occurrence_intermediate_df['eventID'].isna().all():
                    sample_size_map = current_assay_occurrence_intermediate_df.groupby('eventID')['organismQuantity'].sum().to_dict()
                    current_assay_occurrence_intermediate_df['sampleSizeValue'] = current_assay_occurrence_intermediate_df['eventID'].map(sample_size_map)
                current_assay_occurrence_intermediate_df['sampleSizeUnit'] = 'DNA sequence reads'
                
                if 'parentEventID' in dwc_data['occurrence'].index:
                    faire_parent_event_id_col_name = str(dwc_data['occurrence'].loc['parentEventID','FAIRe_term']).strip() 
                    actual_parent_event_col_name = faire_parent_event_id_col_name + "_sm" if faire_parent_event_id_col_name + "_sm" in current_assay_occurrence_intermediate_df.columns else faire_parent_event_id_col_name
                    if actual_parent_event_col_name in current_assay_occurrence_intermediate_df.columns :
                         current_assay_occurrence_intermediate_df['parentEventID'] = current_assay_occurrence_intermediate_df['parentEventID'].fillna(current_assay_occurrence_intermediate_df[actual_parent_event_col_name])

                # --- STEP 9: Final Column Selection and Order for this assay's DataFrame ---
                for col_final_desired in DESIRED_OCCURRENCE_COLUMNS_IN_ORDER:
                    if col_final_desired not in current_assay_occurrence_intermediate_df.columns:
                        current_assay_occurrence_intermediate_df[col_final_desired] = pd.NA
                
                current_assay_occurrence_final_df = current_assay_occurrence_intermediate_df[DESIRED_OCCURRENCE_COLUMNS_IN_ORDER].copy()
                
                all_processed_occurrence_dfs.append(current_assay_occurrence_final_df)
                successful_runs += 1
                reporter.add_text(f"  Successfully processed {analysis_run_name}: Generated {len(current_assay_occurrence_final_df)} records.")

            except Exception as e:
                import traceback
                reporter.add_text(f"  ❌ Error processing {analysis_run_name}: {str(e)}")
                reporter.add_text(f"  Traceback for {analysis_run_name}: {traceback.format_exc()}")
                failed_runs += 1

        # --- POST-LOOP CONCATENATION & FINALIZATION ---
        reporter.add_text(f"🏁 LOOP COMPLETED: Successful runs: {successful_runs}, Failed runs: {failed_runs}, Total DataFrames to combine: {len(all_processed_occurrence_dfs)}")

        if all_processed_occurrence_dfs:
            occ_all_final_combined = pd.concat(all_processed_occurrence_dfs, ignore_index=True, sort=False)
            
            # --- DATA QUALITY CHECKS AND DEFAULTS ---
            reporter.add_section("Data Quality Checks", level=3)
            
            # 1. Check and apply default for geodeticDatum
            if 'geodeticDatum' not in occ_all_final_combined.columns or occ_all_final_combined['geodeticDatum'].isna().all():
                occ_all_final_combined['geodeticDatum'] = 'WGS84'
                reporter.add_warning("The 'geodeticDatum' column was empty or missing. It has been filled with the default value 'WGS84'. Please verify this is appropriate for your data.")
            
            # 2. Check for empty locationID and report in the quality check section
            if 'locationID' not in occ_all_final_combined.columns or occ_all_final_combined['locationID'].isna().all() or (occ_all_final_combined['locationID'] == 'not applicable').all():
                reporter.add_warning("The 'locationID' column is empty. This is likely because the 'line_id' and/or 'station_id' columns were not found in the source sampleMetadata. The pipeline cannot generate a default for this field.")
            
            # 3. Fix eventDate formatting and warn if time is missing
            def format_event_date(d):
                # If the value is a string that can be parsed as a date
                if isinstance(d, str):
                    try:
                        # Attempt to parse it, this will recognize 'YYYY-MM-DD HH:MM:SS'
                        dt_obj = pd.to_datetime(d)
                        # If it's exactly midnight, format as date only
                        if dt_obj.time() == datetime.time(0, 0):
                            return dt_obj.strftime('%Y-%m-%d')
                    except (ValueError, TypeError):
                        # If parsing fails, it's a string like 'not applicable', so keep it as is
                        return d
                # If it's already a datetime object (less likely but possible)
                elif isinstance(d, (datetime.datetime, pd.Timestamp)):
                    if d.time() == datetime.time(0, 0):
                        return d.strftime('%Y-%m-%d')
                return d

            if 'eventDate' in occ_all_final_combined.columns:
                occ_all_final_combined['eventDate'] = occ_all_final_combined['eventDate'].apply(format_event_date)
                
                # Check if any dates are missing a time component (now that midnight is stripped)
                if any(isinstance(d, str) and len(d) == 10 for d in occ_all_final_combined['eventDate']):
                    reporter.add_warning("The 'eventDate' column contains dates without a specific time. If a more precise time is available in your source data, please include it.")


            for col_final_desired in DESIRED_OCCURRENCE_COLUMNS_IN_ORDER:
                if col_final_desired not in occ_all_final_combined.columns:
                    occ_all_final_combined[col_final_desired] = pd.NA
            
            occ_all_final_output = occ_all_final_combined.reindex(columns=DESIRED_OCCURRENCE_COLUMNS_IN_ORDER)
            
            original_rows_before_dedup = len(occ_all_final_output)
            if 'occurrenceID' in occ_all_final_output.columns and not occ_all_final_output['occurrenceID'].isna().all():
                num_duplicates = occ_all_final_output.duplicated(subset=['occurrenceID']).sum()
                if num_duplicates > 0:
                    occ_all_final_output.drop_duplicates(subset=['occurrenceID'], keep='first', inplace=True)
                    reporter.add_text(f"🔄 Dropped {num_duplicates} duplicate occurrenceID records. Final rows: {len(occ_all_final_output)}.")
                else:
                    reporter.add_text("🔄 No duplicate occurrenceID records found to drop.")
            else:
                reporter.add_text("  ⚠️ WARNING: 'occurrenceID' column not found or is all NA. Cannot effectively drop duplicates based on it.")

            try:
                occ_all_final_output.to_csv(output_path, index=False, na_rep='') 
                reporter.add_text(f"💾 Combined occurrence file '{output_filename}' saved to '{output_path}' with {len(occ_all_final_output)} records.")
                reporter.add_text(f"👀 Preview of final combined occurrence data (first 5 rows, selected columns):")
                preview_cols_subset = ['eventID', 'occurrenceID', 'assay_name', 'parentEventID', 'datasetID', 'recordedBy', 
                                       'locality', 'decimalLatitude', 'decimalLongitude', 'geodeticDatum', 
                                       'identificationRemarks', 'locationID']
                preview_cols_to_show = [col for col in preview_cols_subset if col in occ_all_final_output.columns]
                reporter.add_dataframe(occ_all_final_output[preview_cols_to_show], "Occurrence Core Preview", max_rows=5)
                
                # Add important note about assay_name column
                reporter.add_text("<h4>📝 NOTE:</h4>")
                reporter.add_text("The Occurrence Core at this step (before taxonomic assignment through WoRMS or GBIF), contains an assay_name column. This will be removed from the final Occurrence Core (after taxonomic assignment) but it is used by the taxonomic assignment code to know which assay's data you want to remove the 'species' rank from consideration. This is because some assays, like 16S for example, return non-usable assignments at species level, while, for example, 18S species assignments ARE useful.")
                
                reporter.add_success(f"Successfully created occurrence core with {len(occ_all_final_output)} records")
                
                # Return both the final output and the individual dataframes for downstream use
                return occ_all_final_output, all_processed_occurrence_dfs
                
            except Exception as e:
                error_msg = f"❌ Error saving combined occurrence file: {str(e)}"
                reporter.add_error(error_msg)
                raise Exception(error_msg)
        else:
            error_msg = f"❌ No data to combine - all analysis runs may have failed or yielded no occurrence records."
            reporter.add_error(error_msg)
            raise Exception(error_msg)
            
    except Exception as e:
        error_msg = f"Failed to create occurrence core: {str(e)}"
        reporter.add_error(error_msg)
        raise Exception(error_msg) 


def get_final_occurrence_column_order():
    """
    Returns a list defining the desired final column order for the occurrence core file.
    This provides a single, authoritative source for column ordering.
    """
    # This function is now deprecated as the column order is derived from the data_mapper.yaml
    # It is kept for reference and potential future use but is no longer called by the main process.
    return [
        'occurrenceID', 'eventID', 'verbatimIdentification', 'kingdom', 'phylum', 'class', 
        'order', 'family', 'genus', 'scientificName', 'taxonID', 'scientificNameID', 'taxonRank','parentEventID',  'datasetID', 'locationID', 'basisOfRecord', 
        'occurrenceStatus', 'organismQuantity', 'organismQuantityType', 
        'sampleSizeValue', 'sampleSizeUnit', 'recordedBy', 'materialSampleID', 
        'eventDate', 'minimumDepthInMeters', 'maximumDepthInMeters', 'locality', 
        'decimalLatitude', 'decimalLongitude', 'geodeticDatum', 
        'identificationRemarks', 'nameAccordingTo', 'associatedSequences'
        # 'match_type_debug' is intentionally excluded as it's for internal review in the INFO file.
    ] 