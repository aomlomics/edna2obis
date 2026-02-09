"""
Occurrence Core Builder for edna2obis
Contains the complex logic for building occurrence cores from ASV data
Build a combined Occurrence Core for all Analyses
"""

import pandas as pd
import numpy as np
import os
import re
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

        reporter.add_text(f"Starting data processing for {len(params['datafiles'])} analysis run(s) to generate occurrence records.")

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
            reporter.add_text("  Warning: projectMetadata is empty or not found. 'recordedBy' and 'datasetID' will use default placeholder values.")
            project_bibliographic_citation = pd.NA
            project_license = pd.NA
            project_institution_id = pd.NA

        # Note: Abundance table columns are now lib_id (not samp_name), so we no longer need
        # an event_id_lookup map. The melted column from the abundance table IS the lib_id/eventID.


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
                    # Abundance table columns may be lib_id or (legacy) samp_name
                    current_assay_occ_melted = pd.melt(
                        current_abundance_df_raw, id_vars=['featureid'],
                        var_name='lib_id', value_name='organismQuantity'
                    )
                    # Coerce to numeric so filter works even when read as strings (e.g. from Excel)
                    current_assay_occ_melted['organismQuantity'] = pd.to_numeric(current_assay_occ_melted['organismQuantity'], errors='coerce').fillna(0)
                    current_assay_occ_melted = current_assay_occ_melted[current_assay_occ_melted['organismQuantity'] > 0.0]
                    
                    # Normalize lib_id for safe joins
                    current_assay_occ_melted['lib_id'] = current_assay_occ_melted['lib_id'].astype(str).str.strip()
                    
                    # If abundance columns were samp_names (not lib_id), map to lib_id via ERM for this assay
                    final_assay_name = params.get('datafiles', {}).get(analysis_run_name, {}).get('assay_name')
                    if final_assay_name and 'experimentRunMetadata' in data and not data['experimentRunMetadata'].empty:
                        erm_df = data['experimentRunMetadata'].copy()
                        if 'lib_id' in erm_df.columns and 'samp_name' in erm_df.columns and 'assay_name' in erm_df.columns:
                            erm_for_assay = erm_df[erm_df['assay_name'].astype(str).str.strip() == str(final_assay_name).strip()]
                            valid_lib_ids = set(erm_for_assay['lib_id'].astype(str).str.strip().unique())
                            # If melted lib_id values are not in ERM lib_id, treat them as samp_name and map to lib_id
                            samp_to_lib = erm_for_assay.set_index(erm_for_assay['samp_name'].astype(str).str.strip())['lib_id'].astype(str).str.strip().to_dict()
                            def resolve_lib_id(val):
                                v = str(val).strip()
                                if v in valid_lib_ids:
                                    return v
                                return samp_to_lib.get(v, v)
                            current_assay_occ_melted['lib_id'] = current_assay_occ_melted['lib_id'].map(resolve_lib_id)
                    
                    # --- Validation: Check that all abundance lib_ids exist in experimentRunMetadata for this assay ---
                    # Get assay_name for this run
                    final_assay_name = params.get('datafiles', {}).get(analysis_run_name, {}).get('assay_name')
                    
                    # Validate lib_ids exist in ERM for this assay
                    if final_assay_name and 'experimentRunMetadata' in data and not data['experimentRunMetadata'].empty:
                        erm_df = data['experimentRunMetadata'].copy()
                        if 'lib_id' in erm_df.columns and 'assay_name' in erm_df.columns:
                            erm_for_assay = erm_df[erm_df['assay_name'].astype(str).str.strip() == str(final_assay_name).strip()]
                            if len(erm_for_assay) == 0:
                                reporter.add_warning(f"No ERM rows match assay_name '{final_assay_name}'. Check for typos or case mismatch.")
                            else:
                                valid_lib_ids = set(erm_for_assay['lib_id'].astype(str).str.strip().unique())
                                abundance_lib_ids = set(current_assay_occ_melted['lib_id'].unique())
                                missing_lib_ids = abundance_lib_ids - valid_lib_ids
                                if missing_lib_ids:
                                    reporter.add_warning(
                                        f"Found {len(missing_lib_ids)} lib_id(s) in abundance table not in ERM for assay '{final_assay_name}': "
                                        f"{', '.join(list(missing_lib_ids)[:10])}{'...' if len(missing_lib_ids) > 10 else ''}"
                                    )
                    
                    current_assay_occurrence_intermediate_df = pd.merge(
                        current_assay_occ_melted, current_tax_df_processed,
                        on='featureid', how='left'
                    )
                    

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

                    # --- STEP 5: Merge `experimentRunMetadata` first (on lib_id) to get samp_name and ERM fields ---
                    final_assay_name = params.get('datafiles', {}).get(analysis_run_name, {}).get('assay_name')
                    
                    if 'experimentRunMetadata' in data and not data['experimentRunMetadata'].empty:
                        erm_df_to_merge = data['experimentRunMetadata'].copy()
                        if 'lib_id' in erm_df_to_merge.columns:
                            # Normalize lib_id on both sides for safe merge
                            current_assay_occurrence_intermediate_df['lib_id'] = current_assay_occurrence_intermediate_df['lib_id'].astype(str).str.strip()
                            erm_df_to_merge['lib_id'] = erm_df_to_merge['lib_id'].astype(str).str.strip()
                            
                            # Drop columns from left side that exist in ERM (except lib_id) to prevent shadowing
                            erm_cols = set(erm_df_to_merge.columns)
                            left_cols_to_drop = [col for col in current_assay_occurrence_intermediate_df.columns 
                                                 if col in erm_cols and col != 'lib_id']
                            if left_cols_to_drop:
                                current_assay_occurrence_intermediate_df = current_assay_occurrence_intermediate_df.drop(columns=left_cols_to_drop)
                            
                            # Merge with ERM
                            current_assay_occurrence_intermediate_df = pd.merge(
                                current_assay_occurrence_intermediate_df,
                                erm_df_to_merge,
                                on='lib_id',
                                how='left'
                            )
                            
                            # Filter to rows for this assay only
                            if final_assay_name and 'assay_name' in current_assay_occurrence_intermediate_df.columns:
                                assay_mask = current_assay_occurrence_intermediate_df['assay_name'].astype(str).str.strip() == str(final_assay_name).strip()
                                if not assay_mask.all():
                                    current_assay_occurrence_intermediate_df = current_assay_occurrence_intermediate_df.loc[assay_mask].copy()
                            # Populate DwC columns from experimentRunMetadata based on the mapper
                            for dwc_col_target, mapping_info in occurrence_map_config.items():
                                if isinstance(mapping_info, dict) and mapping_info.get('source') == 'experimentRunMetadata':
                                    faire_col_source_original = str(mapping_info.get('faire_term')).strip()
                                    source_col_in_df = None
                                    if faire_col_source_original + '_erm' in current_assay_occurrence_intermediate_df.columns:
                                        source_col_in_df = faire_col_source_original + '_erm'
                                    elif faire_col_source_original in current_assay_occurrence_intermediate_df.columns:
                                        source_col_in_df = faire_col_source_original
                                    if source_col_in_df:
                                        current_assay_occurrence_intermediate_df[dwc_col_target] = current_assay_occurrence_intermediate_df[source_col_in_df]
                        else:
                            reporter.add_warning("experimentRunMetadata missing 'lib_id'; cannot map experiment-level fields.")
                    else:
                        reporter.add_warning(f"  Warning: 'experimentRunMetadata' is empty or not found. Cannot merge for ERM fields for run {analysis_run_name}.")

                    # --- STEP 5b: Merge `sampleMetadata` (on samp_name from ERM) and map FAIRe terms to DwC terms ---
                    if 'sampleMetadata' in data and not data['sampleMetadata'].empty:
                        sm_df_to_merge = data['sampleMetadata'].copy()
                        sm_df_to_merge['samp_name'] = sm_df_to_merge['samp_name'].astype(str).str.strip()
                        # Use samp_name from ERM merge (column is 'samp_name' from right; if conflict use 'samp_name_erm')
                        samp_key = 'samp_name_erm' if 'samp_name_erm' in current_assay_occurrence_intermediate_df.columns else 'samp_name'
                        if samp_key in current_assay_occurrence_intermediate_df.columns:
                            current_assay_occurrence_intermediate_df[samp_key] = current_assay_occurrence_intermediate_df[samp_key].astype(str).str.strip()
                            current_assay_occurrence_intermediate_df = pd.merge(
                                current_assay_occurrence_intermediate_df,
                                sm_df_to_merge,
                                left_on=samp_key,
                                right_on='samp_name',
                                how='left',
                                suffixes=('', '_sm')
                            )
                            # parentEventID = samp_name (event hierarchy; 'samp_name' is from right after merge)
                            if 'samp_name' in current_assay_occurrence_intermediate_df.columns:
                                current_assay_occurrence_intermediate_df['parentEventID'] = current_assay_occurrence_intermediate_df['samp_name']
                            elif samp_key in current_assay_occurrence_intermediate_df.columns:
                                current_assay_occurrence_intermediate_df['parentEventID'] = current_assay_occurrence_intermediate_df[samp_key]
                            # Populate DwC columns from sampleMetadata based on the mapper
                            for dwc_col_target, mapping_info in occurrence_map_config.items():
                                if isinstance(mapping_info, dict) and mapping_info.get('source') == 'sampleMetadata':
                                    faire_col_source_original = str(mapping_info.get('faire_term')).strip()
                                    source_col_in_df = None
                                    if faire_col_source_original + '_sm' in current_assay_occurrence_intermediate_df.columns:
                                        source_col_in_df = faire_col_source_original + '_sm'
                                    if source_col_in_df is None and faire_col_source_original in current_assay_occurrence_intermediate_df.columns:
                                        source_col_in_df = faire_col_source_original
                                    if source_col_in_df is not None:
                                        current_assay_occurrence_intermediate_df[dwc_col_target] = current_assay_occurrence_intermediate_df[source_col_in_df]
                                    elif dwc_col_target in ['locality', 'decimalLatitude', 'decimalLongitude', 'geodeticDatum', 'eventDate']:
                                        reporter.add_text(f"  DIAGNOSTIC: For DwC term '{dwc_col_target}', its mapped FAIRe term '{faire_col_source_original}' (from mapper) was NOT found as a column in the merged sample data. The column will be empty if not populated by other means.")
                        else:
                            reporter.add_warning(f"  Warning: 'samp_name' not found after ERM merge. Cannot merge sampleMetadata for run {analysis_run_name}.")
                    else:
                        reporter.add_text(f"  Warning: 'sampleMetadata' is empty or not found. Cannot merge for DwC term population for run {analysis_run_name}.")

                    # --- Construct 'locationID' with new flexible logic ---
                    
                    # Priority 1: Use existing 'locationID' column if it's already in the data
                    if 'locationID' in current_assay_occurrence_intermediate_df.columns and not current_assay_occurrence_intermediate_df['locationID'].isna().all():
                        reporter.add_text("Found existing 'locationID' column in source data. Using it directly.")
                    
                    else:
                        # Priority 2: Use components from config.yaml (check both col and col_sm from merge)
                        id_components = params.get('locationID_components', [])
                        def _resolve_col(name):
                            if name in current_assay_occurrence_intermediate_df.columns:
                                return name
                            if name + '_sm' in current_assay_occurrence_intermediate_df.columns:
                                return name + '_sm'
                            return None
                        resolved = [(c, _resolve_col(c)) for c in id_components]
                        missing_components = [c for c, r in resolved if r is None]
                        component_cols = [r for _, r in resolved if r is not None]
                        
                        if id_components and len(component_cols) == len(id_components):
                            reporter.add_text(f"Constructing 'locationID' from config components: {', '.join(id_components)}")
                            component_series = [current_assay_occurrence_intermediate_df[col].astype(str).fillna(f"No_{col}") for col in component_cols]
                            current_assay_occurrence_intermediate_df['locationID'] = pd.Series('_'.join(map(str, t)) for t in zip(*component_series))
                        
                        else:
                            # Priority 3: Fallback to default (and warn if component columns were specified but not found)
                            if missing_components:
                                reporter.add_warning(f"Could not construct 'locationID' from config. The following columns were not found: {', '.join(missing_components)}. Falling back to default.")
                            default_components = ['line_id', 'station_id']
                            default_resolved = [_resolve_col(c) for c in default_components]
                            if all(r is not None for r in default_resolved):
                                line_ids = current_assay_occurrence_intermediate_df[default_resolved[0]].astype(str).fillna('NoLineID')
                                station_ids = current_assay_occurrence_intermediate_df[default_resolved[1]].astype(str).fillna('NoStationID')
                                current_assay_occurrence_intermediate_df['locationID'] = line_ids + "_" + station_ids
                            else:
                                # If even the default components can't be found, set to NA
                                 if 'locationID' not in current_assay_occurrence_intermediate_df.columns:
                                      current_assay_occurrence_intermediate_df['locationID'] = pd.NA


                    # --- STEP 6: Define `eventID` and `assay_name` ---
                    # Get assay_name directly from config - it's hardcoded there!
                    if not final_assay_name:
                        reporter.add_error(f"Could not determine assay_name for run '{analysis_run_name}' from config.yaml. Cannot assign eventID.")
                        current_assay_occurrence_intermediate_df['eventID'] = "ERROR_NO_ASSAY_NAME"
                        current_assay_occurrence_intermediate_df['assay_name'] = "ERROR_NO_ASSAY_NAME"
                    else:
                        # Assign the definitive assay_name for this entire run
                        current_assay_occurrence_intermediate_df['assay_name'] = final_assay_name
                        
                        # Set eventID = lib_id (abundance table columns are lib_id)
                        current_assay_occurrence_intermediate_df['eventID'] = current_assay_occurrence_intermediate_df['lib_id']

                    # --- STEP 7: Construct `occurrenceID` ---
                    current_assay_occurrence_intermediate_df['eventID'] = current_assay_occurrence_intermediate_df['eventID'].astype(str)
                    # Namespace occurrenceID by analysis_run_name to prevent collisions when multiple analyses
                    # can share the same lib_id/eventID and featureid pairs (common in internal pipelines/DB schemas).
                    dataset_component = str(project_dataset_id).strip()
                    dataset_component = re.sub(r"\s+", "_", dataset_component)
                    dataset_component = re.sub(r"[^\w\.\-]+", "_", dataset_component)
                    analysis_component = str(analysis_run_name).strip()
                    analysis_component = re.sub(r"\s+", "_", analysis_component)
                    analysis_component = re.sub(r"[^\w\.\-]+", "_", analysis_component)
                    current_assay_occurrence_intermediate_df['occurrenceID'] = (
                        dataset_component
                        + ":"
                        + analysis_component
                        + ":"
                        + current_assay_occurrence_intermediate_df['eventID']
                        + '_occ_'
                        + current_assay_occurrence_intermediate_df['featureid'].astype(str)
                    )

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
                    reporter.add_text(f"  Processed {analysis_run_name}: {len(current_assay_occurrence_final_df)} records.")

                except Exception as e:
                    import traceback
                    reporter.add_text(f"  Error processing {analysis_run_name}: {str(e)}")
                    reporter.add_text(f"  Traceback for {analysis_run_name}: {traceback.format_exc()}")
                    failed_runs += 1

        # --- POST-LOOP CONCATENATION & FINALIZATION ---
        reporter.add_text(f"Loop completed: Successful runs: {successful_runs}, Failed runs: {failed_runs}, Total DataFrames to combine: {len(all_processed_occurrence_dfs)}")

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
                # Prefer removing fully duplicate rows (safe) over dropping by occurrenceID alone (can drop real data).
                before_rows = len(occ_all_final_output)
                occ_all_final_output = occ_all_final_output.drop_duplicates(keep='first')
                removed_full_dups = before_rows - len(occ_all_final_output)
                if removed_full_dups > 0:
                    reporter.add_text(f"Removed {removed_full_dups} fully duplicated row(s) after combining occurrence data. Rows now: {len(occ_all_final_output)}.")
                
                # occurrenceID must be unique within a published dataset (GBIF requirement; OBIS expects globally unique).
                # Do not drop rows based on occurrenceID (would cause silent data loss). Instead fail loudly with diagnostics.
                num_duplicate_ids = int(occ_all_final_output.duplicated(subset=['occurrenceID']).sum())
                if num_duplicate_ids > 0:
                    dup_ids = (
                        occ_all_final_output.loc[occ_all_final_output['occurrenceID'].duplicated(keep=False), 'occurrenceID']
                        .astype(str)
                        .value_counts()
                        .head(10)
                    )
                    reporter.add_error(
                        f"Found {num_duplicate_ids} row(s) with duplicated occurrenceID values after combining occurrence data. "
                        f"This will cause problems for GBIF/OBIS publication and can lead to empty split outputs. "
                        f"Top duplicated occurrenceID values (up to 10 shown): {dup_ids.to_dict()}"
                    )
                    reporter.add_text(
                        "Fix: ensure the components used to build occurrenceID are unique per occurrence. "
                        "Typically this means ensuring eventID/lib_id is unique across cruises/assays, and ensuring analysis_run_name is stable."
                    )
                    raise ValueError("Duplicate occurrenceID values detected in combined occurrence core.")
                else:
                    reporter.add_text("No duplicate occurrenceID values found after combining (and full-row de-duplication).")
            else:
                reporter.add_text("  WARNING: 'occurrenceID' column not found or is all NA. Cannot effectively drop duplicates based on it.")

            try:
                occ_all_final_output.to_csv(output_path, index=False, na_rep='') 
                reporter.add_text(f"Combined occurrence file '{output_filename}' saved to '{output_path}' with {len(occ_all_final_output)} records.")
                reporter.add_text("Preview of final combined occurrence data (first 5 rows, selected columns):")
                preview_cols_subset = ['eventID', 'occurrenceID', 'assay_name', 'parentEventID', 'datasetID', 'recordedBy', 
                                       'locality', 'decimalLatitude', 'decimalLongitude', 'geodeticDatum', 
                                       'identificationRemarks', 'locationID']
                preview_cols_to_show = [col for col in preview_cols_subset if col in occ_all_final_output.columns]
                reporter.add_dataframe(occ_all_final_output[preview_cols_to_show], "Occurrence Core Preview", max_rows=5)
                
                # Add important note about assay_name column
                reporter.add_text("<h4>NOTE:</h4>")
                reporter.add_text("The Occurrence Core at this step (before taxonomic assignment through WoRMS or GBIF), contains an assay_name column. This will be removed from the final Occurrence Core (after taxonomic assignment) but it is used by the taxonomic assignment code to know which assay's data you want to remove the 'species' rank from consideration. This is because some assays, like 16S for example, return non-usable assignments at species level, while, for example, 18S species assignments ARE useful.")
                
                reporter.add_success(f"Successfully created occurrence core with {len(occ_all_final_output)} records")
                
                # Return both the final output and the individual dataframes for downstream use
                return occ_all_final_output, all_processed_occurrence_dfs
                
            except Exception as e:
                error_msg = f"Error saving combined occurrence file: {str(e)}"
                reporter.add_error(error_msg)
                raise Exception(error_msg)
        else:
            error_msg = "No data to combine - all analysis runs may have failed or yielded no occurrence records."
            reporter.add_error(error_msg)
            raise Exception(error_msg)
            
    except Exception as e:
        error_msg = f"Failed to create occurrence core: {str(e)}"
        reporter.add_error(error_msg)
        raise Exception(error_msg) 


def get_final_occurrence_column_order():
    """
    Returns a list defining the desired final column order for the occurrence core file.
    This now derives from data_mapper.yaml to keep the mapper authoritative.
    """
    try:
        import yaml
        # Load mapper
        with open('data_mapper.yaml', 'r', encoding='utf-8') as f:
            mapper = yaml.safe_load(f) or {}
        # Detect metadata_format from config.yaml (default NOAA)
        metadata_format = 'NOAA'
        try:
            with open('config.yaml', 'r', encoding='utf-8') as cf:
                cfg = yaml.safe_load(cf) or {}
                val = cfg.get('metadata_format')
                if isinstance(val, str):
                    metadata_format = val.strip().upper()
        except Exception:
            pass
        format_prefix = 'generic_' if metadata_format == 'GENERIC' else ''
        occurrence_key = f"{format_prefix}occurrence_core"
        occ_map = mapper.get(occurrence_key, {}) or {}
        if isinstance(occ_map, dict) and occ_map:
            return list(occ_map.keys())
    except Exception:
        pass

    # Fallback to legacy static order
    return [
        'occurrenceID', 'eventID', 'verbatimIdentification', 'kingdom', 'phylum', 'class',
        'order', 'family', 'genus', 'scientificName', 'taxonID', 'scientificNameID', 'taxonRank', 'parentEventID', 'datasetID', 'locationID', 'basisOfRecord',
        'occurrenceStatus', 'organismQuantity', 'organismQuantityType',
        'sampleSizeValue', 'sampleSizeUnit', 'recordedBy', 'materialSampleID',
        'eventDate', 'minimumDepthInMeters', 'maximumDepthInMeters', 'locality',
        'decimalLatitude', 'decimalLongitude', 'geodeticDatum',
        'identificationRemarks', 'nameAccordingTo', 'associatedSequences'
    ]