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

        # Define the desired final columns for occurrence.csv IN THE SPECIFIC ORDER REQUIRED
        DESIRED_OCCURRENCE_COLUMNS_IN_ORDER = [
            'eventID', 'organismQuantity', 'assay_name', 'occurrenceID', 'verbatimIdentification',
            'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 
            'scientificName', 'scientificNameID', 'match_type_debug',
            'taxonRank', 'identificationRemarks',
            'taxonID', 'basisOfRecord', 'nameAccordingTo', 'organismQuantityType',
            'recordedBy', 'materialSampleID', 'sampleSizeValue', 'sampleSizeUnit',
            'associatedSequences', 'locationID', 'eventDate', 'minimumDepthInMeters', 'maximumDepthInMeters',
            'locality', 'decimalLatitude', 'decimalLongitude',
            'geodeticDatum', 'parentEventID', 'datasetID', 'occurrenceStatus'
        ]

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

                sequence_col_dwc = 'DNA_sequence'
                sequence_col_input = 'sequence' 
                confidence_col_original_case = 'Confidence' 
                
                if 'dna_sequence' in current_tax_df_raw.columns and sequence_col_input not in current_tax_df_raw.columns:
                     current_tax_df_raw.rename(columns={'dna_sequence': sequence_col_input}, inplace=True)
                elif sequence_col_input not in current_tax_df_raw.columns:
                    current_tax_df_raw[sequence_col_input] = pd.NA
                
                # This is the name of the FAIRe column in current_tax_df_raw that holds the string to be used.
                # It will be renamed to 'verbatimIdentification' in current_tax_df_processed.
                source_column_for_verbatim_id = 'taxonomy' 

                if source_column_for_verbatim_id in current_tax_df_raw.columns:
                    verbatim_id_source_col = source_column_for_verbatim_id 
                    # print(f"  Using column '{verbatim_id_source_col}' from input as source for DwC 'verbatimIdentification'.") # Use for Debug if needed
                else:
                    # If the 'taxonomy' column is MISSING from the input file.
                    reporter.add_text(f"  CRITICAL WARNING: Column '{source_column_for_verbatim_id}' not found in input taxonomy table for '{analysis_run_name}'. DwC 'verbatimIdentification' will use a placeholder.")
                    current_tax_df_raw['verbatimIdentification_placeholder'] = f"Data from '{source_column_for_verbatim_id}' column not available in source"
                    verbatim_id_source_col = 'verbatimIdentification_placeholder'
                    # This placeholder column will be picked up by tax_cols_to_keep and then renamed.                               
                
                tax_cols_to_keep = ['featureid', sequence_col_input, verbatim_id_source_col]
                if confidence_col_original_case in current_tax_df_raw.columns:
                    tax_cols_to_keep.append(confidence_col_original_case)
                else:
                    current_tax_df_raw[confidence_col_original_case] = pd.NA 
                    tax_cols_to_keep.append(confidence_col_original_case)

                current_tax_df_processed = current_tax_df_raw[[col for col in tax_cols_to_keep if col in current_tax_df_raw.columns]].copy()
                current_tax_df_processed.rename(columns={verbatim_id_source_col: 'verbatimIdentification', 
                                                         sequence_col_input: sequence_col_dwc}, inplace=True)

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

                # --- STEP 3: Initialize ALL DwC Fields from DESIRED_OCCURRENCE_COLUMNS_IN_ORDER ---
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
                    
                    # These columns have specific assignment logic/calculation later or are core IDs
                    # They are filled by mapping only if currently NA. They all have different, specific logic to construct.
                    cols_with_specific_logic_or_origin = [
                        'datasetID', 'recordedBy', 'eventID', 'occurrenceID', 'taxonID', 
                        'organismQuantityType', 'occurrenceStatus', 'basisOfRecord', 'nameAccordingTo',
                        'parentEventID', 'associatedSequences', 
                        'sampleSizeValue', 'sampleSizeUnit', 'identificationRemarks'
                    ]

                    # Populate DwC columns using the dwc_data['occurrence'] mapping
                    for dwc_col_target, faire_row in dwc_data['occurrence'].iterrows():
                        if dwc_col_target not in DESIRED_OCCURRENCE_COLUMNS_IN_ORDER: # Ensure we only care about desired output columns
                            continue 

                        faire_col_source_original = str(faire_row['FAIRe_term']).strip()
                        source_col_in_df = None
                        
                        # Check for the column with _sm suffix first (if a clash occurred during merge with sampleMetadata)
                        if faire_col_source_original + '_sm' in current_assay_occurrence_intermediate_df.columns:
                            source_col_in_df = faire_col_source_original + '_sm'
                        # Else, check for the original FAIRe term name (if no clash)
                        elif faire_col_source_original in current_assay_occurrence_intermediate_df.columns:
                            source_col_in_df = faire_col_source_original
                        
                        if source_col_in_df:
                            # If the DwC target column has specific logic for its creation or is a core ID,
                            # only fill it from sampleMetadata if it's currently NA.
                            if dwc_col_target in cols_with_specific_logic_or_origin:
                                current_assay_occurrence_intermediate_df[dwc_col_target] = current_assay_occurrence_intermediate_df[dwc_col_target].fillna(current_assay_occurrence_intermediate_df[source_col_in_df])
                            else: 
                                # For other "standard" DwC columns (like locality, lat, lon, geodeticDatum, etc.),
                                # directly assign from the source FAIRe column.
                                current_assay_occurrence_intermediate_df[dwc_col_target] = current_assay_occurrence_intermediate_df[source_col_in_df]
                        elif dwc_col_target in ['locality', 'decimalLatitude', 'decimalLongitude', 'geodeticDatum', 'eventDate']: # Only print diagnostic for key terms if mapping is missing in checklist
                             reporter.add_text(f"  DIAGNOSTIC: For DwC term '{dwc_col_target}', its mapped FAIRe term '{faire_col_source_original}' (from checklist) was NOT found as a column in the merged sample data (checked as '{faire_col_source_original}' and '{faire_col_source_original}_sm'). The DwC column '{dwc_col_target}' will likely be empty if not populated by other means.")
                else:
                    reporter.add_text(f"  Warning: 'sampleMetadata' is empty or not found. Cannot merge for DwC term population for run {analysis_run_name}.")

                # Construct 'locationID' 
                line_id_col_sm = 'line_id_sm' if 'line_id_sm' in current_assay_occurrence_intermediate_df.columns else 'line_id'
                station_id_col_sm = 'station_id_sm' if 'station_id_sm' in current_assay_occurrence_intermediate_df.columns else 'station_id'

                if line_id_col_sm in current_assay_occurrence_intermediate_df.columns and \
                   station_id_col_sm in current_assay_occurrence_intermediate_df.columns:
                    line_ids = current_assay_occurrence_intermediate_df[line_id_col_sm].astype(str).fillna('NoLineID')
                    station_ids = current_assay_occurrence_intermediate_df[station_id_col_sm].astype(str).fillna('NoStationID')
                    current_assay_occurrence_intermediate_df['locationID'] = line_ids + "_" + station_ids
                    # print(f"    Constructed 'locationID'.") # Reduced verbosity
                else:
                    reporter.add_text(f"    Warning: Could not construct 'locationID' using '{line_id_col_sm}' or '{station_id_col_sm}'.")
                    if 'locationID' not in current_assay_occurrence_intermediate_df.columns or current_assay_occurrence_intermediate_df['locationID'].isna().all():
                         current_assay_occurrence_intermediate_df['locationID'] = "not applicable"


                # --- STEP 6: Merge `experimentRunMetadata` & Define `eventID`, `associatedSequences` ---
                assay_name_for_current_run = next((an_key for an_key, runs_dict_val in data.get('analysis_data_by_assay', {}).items() if isinstance(runs_dict_val, dict) and analysis_run_name in runs_dict_val), None)
                
                current_assay_occurrence_intermediate_df['assay_name'] = assay_name_for_current_run

                if not assay_name_for_current_run:
                    reporter.add_text(f"    ERROR: Could not determine assay_name for '{analysis_run_name}'.")
                    current_assay_occurrence_intermediate_df['eventID'] = current_assay_occurrence_intermediate_df['eventID'].fillna(f"ERROR_eventID_for_{analysis_run_name}")
                    # Also ensure a placeholder for assay_name if it couldn't be found, though it should be an error condition
                    current_assay_occurrence_intermediate_df['assay_name'] = current_assay_occurrence_intermediate_df['assay_name'].fillna(f"UNKNOWN_ASSAY_FOR_{analysis_run_name}")
                elif 'experimentRunMetadata' in data and not data['experimentRunMetadata'].empty:
                    erm_df = data['experimentRunMetadata'].copy()
                    erm_df['samp_name'] = erm_df['samp_name'].astype(str).str.strip()
                    erm_df_assay_specific = erm_df[erm_df['assay_name'].astype(str).str.strip() == str(assay_name_for_current_run).strip()]

                    if not erm_df_assay_specific.empty:
                        faire_lib_id_col = str(dwc_data['occurrence'].loc['eventID', 'FAIRe_term']).strip() if 'eventID' in dwc_data['occurrence'].index else 'lib_id'
                        faire_assoc_seq_col = str(dwc_data['occurrence'].loc['associatedSequences', 'FAIRe_term']).strip() if 'associatedSequences' in dwc_data['occurrence'].index else 'associatedSequences'

                        cols_to_select_from_erm = {'samp_name'}
                        if faire_lib_id_col in erm_df_assay_specific.columns: cols_to_select_from_erm.add(faire_lib_id_col)
                        if faire_assoc_seq_col in erm_df_assay_specific.columns: cols_to_select_from_erm.add(faire_assoc_seq_col)
                        
                        erm_to_merge = erm_df_assay_specific[list(cols_to_select_from_erm)].drop_duplicates(subset=['samp_name']).copy()
                        
                        current_assay_occurrence_intermediate_df = pd.merge(
                            current_assay_occurrence_intermediate_df, erm_to_merge,
                            on='samp_name', how='left', suffixes=('', '_erm')
                        )
                        
                        source_lib_id_col_actual = faire_lib_id_col + '_erm' if faire_lib_id_col + '_erm' in current_assay_occurrence_intermediate_df.columns else faire_lib_id_col
                        if source_lib_id_col_actual in current_assay_occurrence_intermediate_df.columns:
                            current_assay_occurrence_intermediate_df['eventID'] = current_assay_occurrence_intermediate_df['eventID'].fillna(current_assay_occurrence_intermediate_df[source_lib_id_col_actual])
                        
                        source_assoc_seq_col_actual = faire_assoc_seq_col + '_erm' if faire_assoc_seq_col + '_erm' in current_assay_occurrence_intermediate_df.columns else faire_assoc_seq_col
                        if source_assoc_seq_col_actual in current_assay_occurrence_intermediate_df.columns:
                             current_assay_occurrence_intermediate_df['associatedSequences'] = current_assay_occurrence_intermediate_df['associatedSequences'].fillna(current_assay_occurrence_intermediate_df[source_assoc_seq_col_actual])
                    else:
                        current_assay_occurrence_intermediate_df['eventID'] = current_assay_occurrence_intermediate_df['eventID'].fillna(f"NoExpMeta_eventID_for_{analysis_run_name}")
                
                current_assay_occurrence_intermediate_df['eventID'] = current_assay_occurrence_intermediate_df['eventID'].astype(str)
                
                # --- STEP 7: Construct `occurrenceID` ---
                current_assay_occurrence_intermediate_df['occurrenceID'] = current_assay_occurrence_intermediate_df['eventID'] + '_occ_' + current_assay_occurrence_intermediate_df['featureid'].astype(str)

                # --- STEP 8: `identificationRemarks`, `sampleSizeValue`/`Unit`, `parentEventID` ---
                otu_seq_comp_appr_str = "Unknown sequence comparison approach"
                taxa_class_method_str = ""
                taxa_ref_db_str = "Unknown reference DB"

                if assay_name_for_current_run and analysis_run_name in data.get('analysis_data_by_assay', {}).get(assay_name_for_current_run, {}):
                    analysis_meta_df_for_run = data['analysis_data_by_assay'][assay_name_for_current_run][analysis_run_name]
                    if 'term_name' in analysis_meta_df_for_run.columns and 'values' in analysis_meta_df_for_run.columns:
                        def get_analysis_meta(term, df, default):
                            val_series = df[df['term_name'].astype(str).str.strip() == term]['values']
                            return str(val_series.iloc[0]).strip() if not val_series.empty and pd.notna(val_series.iloc[0]) else default
                        otu_seq_comp_appr_str = get_analysis_meta('otu_seq_comp_appr', analysis_meta_df_for_run, otu_seq_comp_appr_str)
                        taxa_class_method_str = get_analysis_meta('taxa_class_method', analysis_meta_df_for_run, taxa_class_method_str)
                        taxa_ref_db_str = get_analysis_meta('otu_db', analysis_meta_df_for_run, taxa_ref_db_str)
                
                confidence_value_series = pd.Series(["unknown confidence"] * len(current_assay_occurrence_intermediate_df), index=current_assay_occurrence_intermediate_df.index, dtype=object)
                if confidence_col_original_case in current_assay_occurrence_intermediate_df.columns:
                    confidence_value_series = current_assay_occurrence_intermediate_df[confidence_col_original_case].astype(str).fillna("unknown confidence")
                
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