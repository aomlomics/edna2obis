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

from html_reporter import HTMLReporter


def create_occurrence_core(data, raw_data_tables, params, dwc_data, reporter: HTMLReporter):
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

        # --- NEW: Add pass-through columns needed by DNA Derived Extension ---
        # This prevents columns like 'concentration' from being dropped before the
        # extension builder can use them.
        dna_derived_key = f"{format_prefix}dna_derived_extension"
        dna_derived_map_config = full_mapper.get(dna_derived_key, {})
        pass_through_cols = []
        if dna_derived_map_config: # Check if the section exists
            for dwc_term, mapping_info in dna_derived_map_config.items():
                if isinstance(mapping_info, dict):
                    source = mapping_info.get('source')
                    # We need to preserve any field from sample or experiment metadata
                    if source in ['sampleMetadata', 'experimentRunMetadata']:
                        faire_term = mapping_info.get('faire_term')
                        if faire_term:
                            pass_through_cols.append(faire_term)
                            # Also keep their corresponding unit columns if they exist
                            pass_through_cols.append(f"{faire_term}_unit")

        # Add unique pass-through columns to the list of columns to keep
        for col in set(pass_through_cols):
            if col not in DESIRED_OCCURRENCE_COLUMNS_IN_ORDER:
                DESIRED_OCCURRENCE_COLUMNS_IN_ORDER.append(col)


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

        def find_and_rename_feature_id(df, df_name, reporter: HTMLReporter):
            """
            Robustly finds and renames the feature ID column in a dataframe.
            Searches for common names and falls back to the first column if none are found.
            """
            cols = df.columns
            original_cols = list(df.columns)
            
            # Create a mapping of lowercased, cleaned-up column names to original names
            cleaned_to_original_map = {str(c).lower().strip().replace('_', '-').replace(' ', ''): c for c in original_cols}
            
            feature_id_col_original = None
            
            # Common names to search for
            possible_names = ['featureid', 'feature-id', 'seq-id', 'seq_id', 'otu id', '#otuid', 'asv-id']
            
            found = False
            for name in possible_names:
                if name in cleaned_to_original_map:
                    feature_id_col_original = cleaned_to_original_map[name]
                    found = True
                    break
            
            if found:
                df.rename(columns={feature_id_col_original: 'featureid'}, inplace=True)
            else:
                # Fallback to first column if no standard name is found
                feature_id_col_original = original_cols[0]
                df.rename(columns={feature_id_col_original: 'featureid'}, inplace=True)
                reporter.add_text(f"  DIAGNOSTIC: Could not find a standard feature ID column (e.g., 'featureid', 'seq_id') in the {df_name}. Using the first column, <strong>'{feature_id_col_original}'</strong>, as the feature ID. Please ensure this is correct.")
                
            return df

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
            
            # --- NEW: Get project-level values for the new fields ---
            project_bibliographic_citation = get_project_meta_value(project_meta_df, 'bibliographicCitation', pd.NA)
            project_license = get_project_meta_value(project_meta_df, 'license', pd.NA)
            project_institution_id = get_project_meta_value(project_meta_df, 'institutionID', pd.NA)
        else:
            reporter.add_text("  ⚠️ Warning: projectMetadata is empty or not found. 'recordedBy' and 'datasetID' will use default placeholder values.")
            project_bibliographic_citation = pd.NA
            project_license = pd.NA
            project_institution_id = pd.NA

        # --- NEW: Pre-build a robust lookup map for eventID from experimentRunMetadata ---
        event_id_lookup = {}
        if 'experimentRunMetadata' in data and not data['experimentRunMetadata'].empty:
            erm_df = data['experimentRunMetadata'].copy()
            # Ensure required columns are present and cleaned
            if 'samp_name' in erm_df.columns and 'lib_id' in erm_df.columns and 'assay_name' in erm_df.columns:
                erm_df['samp_name'] = erm_df['samp_name'].astype(str).str.strip()
                erm_df['assay_name'] = erm_df['assay_name'].astype(str).str.strip()
                erm_df['lib_id'] = erm_df['lib_id'].astype(str).str.strip()
                # Create a multi-index for fast lookup
                event_id_lookup = erm_df.set_index(['samp_name', 'assay_name'])['lib_id'].to_dict()
                reporter.add_text(f"✓ Built eventID lookup map from experimentRunMetadata with {len(event_id_lookup)} entries.")
            else:
                reporter.add_warning("Could not build eventID lookup map: 'samp_name', 'lib_id', or 'assay_name' columns are missing from experimentRunMetadata.")


        # --- MAIN DATA PROCESSING LOOP ---
        # We iterate through the metadata structure, which is the authoritative source for analysis runs.
        for assay_name, analysis_runs in data.get('analysis_data_by_assay', {}).items():
            if not isinstance(analysis_runs, dict): continue

            for analysis_run_name, analysis_df in analysis_runs.items():
                reporter.add_text(f"Processing Analysis Run: {analysis_run_name}")
                try:
                    # --- Find the correct raw data for this run ---
                    # For NOAA mode, analysis_run_name comes from the Excel sheet.
                    # For GENERIC mode, analysis_run_name is the key from config.yaml.
                    # raw_data_tables is always keyed by the name from config.yaml.
                    raw_data_key = analysis_run_name
                    if params.get('metadata_format') == 'NOAA':
                        # For NOAA, we need to find which config key corresponds to this run.
                        # This is brittle, assumes a 1-to-1 mapping and that the user has set it up correctly.
                        # A better long-term solution might be to store the config key in the metadata.
                        # For now, we find the config key that has the matching assay_name.
                        # This assumes one run per assay in the config for NOAA mode.
                        raw_data_key = next((key for key, values in params['datafiles'].items() if values.get('assay_name') == assay_name), None)
                        if not raw_data_key:
                            reporter.add_warning(f"Could not find a matching entry in the 'datafiles' section of config.yaml for assay '{assay_name}'. Skipping run '{analysis_run_name}'.")
                            failed_runs += 1
                            continue

                    # --- STEP 1: Load and Prepare Raw Taxonomy and Abundance Data ---
                    if not (raw_data_key in raw_data_tables and
                            'taxonomy' in raw_data_tables[raw_data_key] and not raw_data_tables[raw_data_key]['taxonomy'].empty and
                            'occurrence' in raw_data_tables[raw_data_key] and not raw_data_tables[raw_data_key]['occurrence'].empty):
                        reporter.add_text(f"  Skipping {analysis_run_name}: Raw taxonomy or occurrence data is missing or empty (looked for key: '{raw_data_key}').")
                        failed_runs += 1
                        continue

                    current_tax_df_raw = raw_data_tables[raw_data_key]['taxonomy'].copy()
                    current_abundance_df_raw = raw_data_tables[raw_data_key]['occurrence'].copy()
                    
                    # --- NEW: Robustly find and rename feature ID ---
                    current_tax_df_raw = find_and_rename_feature_id(current_tax_df_raw, f"taxonomy file for '{analysis_run_name}'", reporter)
                    current_abundance_df_raw = find_and_rename_feature_id(current_abundance_df_raw, f"abundance file for '{analysis_run_name}'", reporter)

                    # Normalize featureid keys as strings without surrounding whitespace for safe joins
                    if 'featureid' in current_tax_df_raw.columns:
                        current_tax_df_raw['featureid'] = current_tax_df_raw['featureid'].astype(str).str.strip()
                    if 'featureid' in current_abundance_df_raw.columns:
                        current_abundance_df_raw['featureid'] = current_abundance_df_raw['featureid'].astype(str).str.strip()

                    # --- NEW, MAPPER-DRIVEN LOGIC ---
                    # 1. Identify all faire_terms sourced from 'taxonomy' and create a rename map
                    tax_source_faire_terms = {'featureid'}  # Always include featureid
                    rename_map_for_tax = {}
                    
                    DWC_RANKS_ORDER = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'scientificName']
                    all_tax_fields = DWC_RANKS_ORDER + ['verbatimIdentification', 'taxonRank']

                    for dwc_term, mapping_info in occurrence_map_config.items():
                        if isinstance(mapping_info, dict) and mapping_info.get('source') == 'taxonomy':
                            faire_term = mapping_info.get('faire_term')
                            if faire_term and faire_term in current_tax_df_raw.columns:
                                tax_source_faire_terms.add(faire_term)
                                # If the FAIRe term is different from the DwC term, map it for renaming
                                if faire_term != dwc_term:
                                    rename_map_for_tax[faire_term] = dwc_term
                            # Also, if the DwC term itself is a column (e.g., 'kingdom'), keep it.
                            elif dwc_term in current_tax_df_raw.columns:
                                 tax_source_faire_terms.add(dwc_term)


                    # 2. Filter the raw dataframe to keep only the columns we need, safely
                    final_cols_to_keep = [col for col in tax_source_faire_terms if col in current_tax_df_raw.columns]
                    current_tax_df_processed = current_tax_df_raw[final_cols_to_keep].copy()
                    
                    # 3. Rename columns from FAIRe terms to Darwin Core terms
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
                    
                    # Optional diagnostic: how many taxonomy rows matched
                    try:
                        tax_cols_present = [c for c in ['verbatimIdentification','kingdom','scientificName'] if c in current_tax_df_processed.columns]
                        if tax_cols_present:
                            matched = current_assay_occurrence_intermediate_df[tax_cols_present].notna().any(axis=1).sum()
                            total = len(current_assay_occurrence_intermediate_df)
                            reporter.add_text(f"  Join diagnostic: taxonomy fields present for {matched}/{total} rows after merge.")
                    except Exception:
                        pass

                    # --- FIX for GENERIC format: Check for a 'Taxon' column to use as verbatimIdentification ---
                    if params.get('metadata_format') == 'GENERIC':
                        # Check for 'Taxon' or 'taxonomy' (case-insensitive) - this is a common generic format variant
                        taxon_col_found = None
                        for col in current_assay_occurrence_intermediate_df.columns:
                            cleaned_col_name = str(col).lower().strip()
                            if cleaned_col_name == 'taxon' or cleaned_col_name == 'taxonomy':
                                taxon_col_found = col
                                break
                    
                        if taxon_col_found and 'verbatimIdentification' not in current_assay_occurrence_intermediate_df.columns:
                            current_assay_occurrence_intermediate_df.rename(columns={taxon_col_found: 'verbatimIdentification'}, inplace=True)
                            reporter.add_text(f"  DIAGNOSTIC: For GENERIC format, found and used column <strong>'{taxon_col_found}'</strong> to populate 'verbatimIdentification'.")


                    # --- NEW: Robustly construct verbatimIdentification ---
                    # This ensures the column exists for the taxonomic matching step.
                    if 'verbatimIdentification' not in current_assay_occurrence_intermediate_df.columns:
                        current_assay_occurrence_intermediate_df['verbatimIdentification'] = pd.NA

                    # Check if verbatimIdentification is completely empty or mostly empty
                    verbatim_empty = current_assay_occurrence_intermediate_df['verbatimIdentification'].isna().sum() / len(current_assay_occurrence_intermediate_df) > 0.9
                    
                    # If verbatimIdentification is mostly empty, construct it from individual rank columns that ARE present.
                    if verbatim_empty:
                        rank_cols_present = [col for col in DWC_RANKS_ORDER if col in current_assay_occurrence_intermediate_df.columns]
                        
                        if rank_cols_present:
                            reporter.add_text(f"  DIAGNOSTIC: verbatimIdentification is mostly empty. Constructing it from available taxonomic rank columns: {rank_cols_present}")
                            # Create a temporary DataFrame with just the rank columns, filled with empty strings for missing values
                            temp_ranks_df = current_assay_occurrence_intermediate_df[rank_cols_present].astype(str).fillna('')
                            
                            # Join the ranks with semicolons, but only include non-empty parts
                            # This creates the full taxonomic string, e.g., "Eukaryota;Dinoflagellata;..."
                            joined_series = temp_ranks_df.apply(lambda row: ';'.join(part for part in row if part and part.lower() not in ['nan', '', 'none']), axis=1)
                            
                            # Fill any remaining empty verbatimIdentification rows with the newly constructed string
                            # This uses the existing verbatim string if present, otherwise fills with the constructed one.
                            current_assay_occurrence_intermediate_df['verbatimIdentification'] = \
                                current_assay_occurrence_intermediate_df['verbatimIdentification'].fillna(joined_series)


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
                    
                    # --- NEW: Assign the project-level values we fetched earlier ---
                    current_assay_occurrence_intermediate_df['bibliographicCitation'] = project_bibliographic_citation
                    current_assay_occurrence_intermediate_df['license'] = project_license
                    current_assay_occurrence_intermediate_df['institutionID'] = project_institution_id

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


                    # --- STEP 6: Define `eventID` and `assay_name` (REWRITTEN) ---
                    # Get assay_name directly from config - it's hardcoded there!
                    final_assay_name = params.get('datafiles', {}).get(analysis_run_name, {}).get('assay_name')
                    
                    if not final_assay_name:
                        reporter.add_error(f"Could not determine assay_name for run '{analysis_run_name}' from config.yaml. Cannot assign eventID.")
                        current_assay_occurrence_intermediate_df['eventID'] = "ERROR_NO_ASSAY_NAME"
                        current_assay_occurrence_intermediate_df['assay_name'] = "ERROR_NO_ASSAY_NAME"
                    else:
                        # Assign the definitive assay_name for this entire run
                        current_assay_occurrence_intermediate_df['assay_name'] = final_assay_name
                        
                        # Use the pre-built lookup map to find the correct eventID
                        def get_event_id(row):
                            return event_id_lookup.get((row['samp_name'], final_assay_name), f"EVENT_ID_NOT_FOUND_FOR_{row['samp_name']}_{final_assay_name}")

                        if event_id_lookup:
                             current_assay_occurrence_intermediate_df['eventID'] = current_assay_occurrence_intermediate_df.apply(get_event_id, axis=1)
                        else:
                            # Fallback if the map couldn't be built
                            current_assay_occurrence_intermediate_df['eventID'] = f"NO_ERM_LOOKUP_FOR_{final_assay_name}"

                    
                    # --- STEP 7: Construct `occurrenceID` ---
                    current_assay_occurrence_intermediate_df['eventID'] = current_assay_occurrence_intermediate_df['eventID'].astype(str)
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