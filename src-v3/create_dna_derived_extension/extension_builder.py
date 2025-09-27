"""
DNA Derived Extension Builder for edna2obis
Contains the complex logic for building DNA derived extensions from ASV data
Create DNA derived extension file
"""

import pandas as pd
import numpy as np
import os


def create_dna_derived_extension(params, data, raw_data_tables, dwc_data, occurrence_core, all_processed_occurrence_dfs, checklist_df, reporter):
    """Create DNA derived extension file - EXACT implementation from original notebook"""
    try:
        reporter.add_section("Creating DNA Derived Extension")
        
        # Create foundation from occurrence core
        reporter.add_text(f"Starting with the Occurrence Core's shape: {occurrence_core.shape}")
        
        # Create foundation with the columns we need, including the correct assay_name
        dna_derived_base = occurrence_core[['eventID', 'parentEventID', 'occurrenceID', 'materialSampleID', 'taxonID', 'assay_name']].copy()
        
        # Rename materialSampleID to source_mat_id for final output
        dna_derived_base = dna_derived_base.rename(columns={'materialSampleID': 'source_mat_id'})
        
        # Extract featureID from taxonID (format: 'ASV:<featureid>')
        dna_derived_base['featureID'] = dna_derived_base['taxonID'].str.replace('ASV:', '')
        
        reporter.add_text(f"DNA derived base shape: {dna_derived_base.shape}")
        
        # Merge with sampleMetadata to get sample collection information
        reporter.add_text("Merging with sample metadata...")
        
        dna_derived_with_sample = dna_derived_base.merge(
            data['sampleMetadata'], 
            left_on='parentEventID', 
            right_on='samp_name', 
            how='left'
        )
        
        if 'assay_name_x' in dna_derived_with_sample.columns:
            dna_derived_with_sample.rename(columns={'assay_name_x': 'assay_name'}, inplace=True)
            if 'assay_name_y' in dna_derived_with_sample.columns:
                dna_derived_with_sample.drop(columns=['assay_name_y'], inplace=True)
        
        reporter.add_text("After merge with sample metadata - shape: {}".format(dna_derived_with_sample.shape))
        
        # Merge with experimentRunMetadata using eventID (from occurrence) and lib_id (from experiment metadata)
        reporter.add_text("Merging with experiment run metadata...")
        dna_derived_with_experiment = dna_derived_with_sample.merge(
            data['experimentRunMetadata'],
            left_on='eventID',
            right_on='lib_id',
            how='left',
            suffixes=('', '_exp')
        )

        cols_to_drop = [col for col in dna_derived_with_experiment.columns if col.endswith('_exp')]
        if 'lib_id' in dna_derived_with_experiment.columns:
            cols_to_drop.append('lib_id')
        
        if cols_to_drop:
            dna_derived_with_experiment.drop(columns=cols_to_drop, inplace=True)

        reporter.add_text("After merge with experiment metadata - shape: {}".format(dna_derived_with_experiment.shape))
        
        # --- Dynamically add unit columns FIRST ---
        # This must run before the generic DwC mapping, to ensure unit columns are found
        # before they can be renamed by the checklist mapping.
        reporter.add_text("Dynamically generating unit columns...")
        created_unit_columns = {} # Maps a field to its new unit column name

        # Use a copy of the dataframe for this step to avoid SettingWithCopyWarning
        dna_derived_df_final = dna_derived_with_experiment.copy()

        # --- NEW: Get the list of desired columns directly from the mapper ---
        # The order of fields in the mapper will determine the final column order.
        dna_derived_mapping = dwc_data['dnaDerived']
        
        # This is a dictionary of DwC_term -> {faire_term: ..., source: ...}
        # We need to load the full mapper here, not just the simple df from main.py
        import yaml
        with open('data_mapper.yaml', 'r', encoding='utf-8') as f:
            full_mapper = yaml.safe_load(f)
        
        format_prefix = "generic_" if params.get('metadata_format') == 'GENERIC' else ""
        dna_derived_key = f"{format_prefix}dna_derived_extension"
        dna_derived_map_config = full_mapper.get(dna_derived_key, {})

        DESIRED_DNA_DERIVED_COLUMNS = list(dna_derived_map_config.keys())

        for field in DESIRED_DNA_DERIVED_COLUMNS:
            unit_col_from_faire = f"{field}_unit"
            dwc_unit_col = f"{field}Unit"

            # Strategy 1: Find a corresponding <field>_unit column from the merged data
            if unit_col_from_faire in dna_derived_df_final.columns:
                # Only create the unit column if it's not completely empty
                if not dna_derived_df_final[unit_col_from_faire].isna().all():
                    dna_derived_df_final.rename(columns={unit_col_from_faire: dwc_unit_col}, inplace=True)
                    created_unit_columns[field] = dwc_unit_col
                    reporter.add_text(f"✓ Created '{dwc_unit_col}' from existing FAIRe column '{unit_col_from_faire}'.")

            # Strategy 2: If no _unit column, look up the unit in the FAIRe checklist (legacy support)
            elif checklist_df is not None and not checklist_df.empty:
                checklist_row = checklist_df[checklist_df['term_name'] == field]
                if not checklist_row.empty and 'unit' in checklist_row.columns:
                    unit_value = checklist_row['unit'].iloc[0]
                    # Only create the column if the unit value from the checklist is meaningful
                    if pd.notna(unit_value) and str(unit_value).strip():
                        dna_derived_df_final[dwc_unit_col] = unit_value
                        created_unit_columns[field] = dwc_unit_col
                        reporter.add_text(f"✓ Created '{dwc_unit_col}' using value ('{unit_value}') from FAIRe checklist for field '{field}'.")
        
        # --- NEW UNIFIED LOOP to process fields from projectMetadata and analysisMetadata ---
        reporter.add_text("Processing bioinformatics fields based on data_mapper.yaml...")

        for dwc_term, mapping_info in dna_derived_map_config.items():
            if not isinstance(mapping_info, dict): continue

            source = mapping_info.get('source')
            faire_term = mapping_info.get('faire_term')

            # --- Logic for fields sourced from projectMetadata ---
            if source == 'projectMetadata':
                field_row = data['projectMetadata'][data['projectMetadata']['term_name'] == faire_term]
                if not field_row.empty:
                    field_row = field_row.iloc[0]
                    project_level_val = field_row.get('project_level')
                    
                    if pd.notna(project_level_val) and str(project_level_val).strip() and str(project_level_val).lower() != 'nan':
                        dna_derived_df_final[dwc_term] = project_level_val
                    else:
                        assay_value_map = {assay: field_row.get(assay) for assay in dna_derived_df_final['assay_name'].unique()}
                        dna_derived_df_final[dwc_term] = dna_derived_df_final['assay_name'].map(assay_value_map)
                else:
                    dna_derived_df_final[dwc_term] = None

            # --- Logic for fields sourced from analysisMetadata ---
            elif source == 'analysisMetadata':
                # This logic now correctly handles both NOAA and GENERIC formats because main.py
                # has already synthesized the 'analysis_data_by_assay' structure for GENERIC mode.
                analysis_run_order = list(params['datafiles'].keys())
                current_row_start = 0
                for i, df_occurrence_chunk in enumerate(all_processed_occurrence_dfs):
                    if i >= len(analysis_run_order): break
                    
                    analysis_run_name = analysis_run_order[i]
                    assay_key = next((ak for ak, runs in data['analysis_data_by_assay'].items() if analysis_run_name in runs), None)

                    if assay_key:
                        current_row_end = current_row_start + len(df_occurrence_chunk)
                        analysis_df = data['analysis_data_by_assay'][assay_key].get(analysis_run_name)
                        
                        if analysis_df is not None:
                            term_row = analysis_df[analysis_df['term_name'] == faire_term]
                            if not term_row.empty:
                                field_value = term_row.iloc[0]['values']
                                dna_derived_df_final.loc[current_row_start:current_row_end-1, dwc_term] = field_value
                        
                        current_row_start = current_row_end
        
        # --- SPECIAL LOGIC: Construct pcr_primer_reference ---
        reporter.add_text("Constructing 'pcr_primer_reference' field...")

        # Helper function to get the values for these terms, using the same projectMetadata logic
        def get_project_meta_series(faire_term, project_meta_df, assay_map_series):
            field_row = project_meta_df[project_meta_df['term_name'] == faire_term]
            if not field_row.empty:
                field_row = field_row.iloc[0]
                project_level_val = field_row.get('project_level')
                if pd.notna(project_level_val) and str(project_level_val).strip() and str(project_level_val).lower() != 'nan':
                    # Return a Series with the project-level value repeated for all rows
                    return pd.Series([project_level_val] * len(assay_map_series), index=assay_map_series.index)
                else:
                    # Map assay-specific values
                    assay_value_map = {assay: field_row.get(assay) for assay in assay_map_series.unique()}
                    return assay_map_series.map(assay_value_map)
            # Return an empty Series if the term is not found
            return pd.Series([pd.NA] * len(assay_map_series), index=assay_map_series.index)

        # Get the series of values for forward and reverse references
        fwd_series = get_project_meta_series('pcr_primer_reference_forward', data['projectMetadata'], dna_derived_df_final['assay_name'])
        rev_series = get_project_meta_series('pcr_primer_reference_reverse', data['projectMetadata'], dna_derived_df_final['assay_name'])

        # Combine them based on the user's logic
        def combine_references(fwd, rev):
            fwd_str = str(fwd).strip() if pd.notna(fwd) else ""
            rev_str = str(rev).strip() if pd.notna(rev) else ""
            
            if fwd_str and rev_str:
                if fwd_str == rev_str:
                    return fwd_str
                else:
                    return f"{fwd_str}|{rev_str}"
            elif fwd_str:
                return fwd_str
            elif rev_str:
                return rev_str
            else:
                return pd.NA

        # Apply the logic row-wise
        dna_derived_df_final['pcr_primer_reference'] = [combine_references(f, r) for f, r in zip(fwd_series, rev_series)]

        # --- (Existing logic for other sources like taxonomy, merges, etc. follows) ---
        
        # Merge DNA sequences from raw taxonomy files
        reporter.add_text("Merging DNA sequences from raw taxonomy files...")
        
        all_tax_data = []
        if raw_data_tables:
            for analysis_run_name, data_dict in raw_data_tables.items():
                if 'taxonomy' in data_dict:
                    tax_data = data_dict['taxonomy'].copy()
                    all_tax_data.append(tax_data)
        
        if all_tax_data:
            combined_tax_data = pd.concat(all_tax_data, ignore_index=True)
            
            # Find the column containing DNA sequences and feature IDs
            feature_col = None
            dna_col = None
            
            for col in combined_tax_data.columns:
                if 'featureid' in col.lower() or col.lower() == 'otu id':
                    feature_col = col
                elif 'dna_sequence' in col.lower() or 'sequence' in col.lower():
                    dna_col = col
            
            if feature_col and dna_col:
                # Remove duplicates from taxonomy data BEFORE merging
                seq_data = combined_tax_data[[feature_col, dna_col]].drop_duplicates(subset=[feature_col])
                
                # Now merge (should be one-to-one)
                dna_derived_with_sequences = dna_derived_df_final.merge(
                    seq_data, 
                    left_on='featureID', 
                    right_on=feature_col, 
                    how='left'
                )
                
                # Rename the DNA sequence column to match our desired output
                dna_derived_with_sequences = dna_derived_with_sequences.rename(columns={dna_col: 'DNA_sequence'})
                
                # Drop the extra feature column from the merge if it's different from our featureID
                if feature_col in dna_derived_with_sequences.columns and feature_col != 'featureID':
                    dna_derived_with_sequences = dna_derived_with_sequences.drop(columns=[feature_col])
            else:
                dna_derived_with_sequences = dna_derived_df_final.copy()
                dna_derived_with_sequences['DNA_sequence'] = None
        else:
            dna_derived_with_sequences = dna_derived_df_final.copy()
            dna_derived_with_sequences['DNA_sequence'] = None
        
        # Apply Darwin Core mappings and create final output
        reporter.add_text("Applying Darwin Core mappings...")
        
        # --- This section is now simplified as most columns are already named correctly ---
        # The main purpose of this step is now to rename any remaining faire_terms to their DwC equivalent
        # if they came from sampleMetadata or experimentRunMetadata merges.
        
        final_df = dna_derived_with_sequences.copy()
        
        rename_map = {}
        for dwc_term, mapping_info in dna_derived_map_config.items():
             if not isinstance(mapping_info, dict): continue
             source = mapping_info.get('source')
             faire_term = mapping_info.get('faire_term')
             if source in ['sampleMetadata', 'experimentRunMetadata'] and faire_term in final_df.columns:
                 rename_map[faire_term] = dwc_term
        
        final_df.rename(columns=rename_map, inplace=True)

        # --- Ensure we have all required columns and order them correctly ---
        # The unit columns will be placed directly after their corresponding data column.
        final_ordered_columns = []
        for col in DESIRED_DNA_DERIVED_COLUMNS:
            final_ordered_columns.append(col)
            # If we just created a unit for this column, insert it into the list
            if col in created_unit_columns:
                final_ordered_columns.append(created_unit_columns[col])

        # Ensure all columns exist in the dataframe, adding empty ones if needed.
        for col in final_ordered_columns:
            if col not in final_df.columns:
                final_df[col] = pd.NA
        
        # Select only the required columns in the correct order
        dna_derived_extension = final_df[final_ordered_columns]
        
        reporter.add_text(f"Final DNA derived extension shape: {dna_derived_extension.shape}")
        
        # Add a preview of the dataframe to the HTML report
        reporter.add_dataframe(dna_derived_extension, "DNA Derived Extension Preview", max_rows=10)

        # Save DNA derived extension to CSV file
        reporter.add_text("Saving DNA derived extension to CSV...")
        
        output_dir = params.get('output_dir', '../processed-v3/')
        os.makedirs(output_dir, exist_ok=True)
        
        output_path = os.path.join(output_dir, "dna_derived_extension.csv")
        dna_derived_extension.to_csv(output_path, index=False, encoding='utf-8-sig') # Encoding helps with special characters in units
        
        reporter.add_success("DNA derived extension created successfully")
        reporter.add_text(f"Saved DNA derived extension: {len(dna_derived_extension):,} records")
        reporter.add_text(f"Output file: dna_derived_extension.csv")
        
        # Verify the file was created
        if os.path.exists(output_path):
            file_size = os.path.getsize(output_path) / (1024*1024)  # Size in MB
            reporter.add_text(f"File size: {file_size:.2f} MB")
            print(f"✅ DNA derived extension created! Saved {len(dna_derived_extension):,} records")
        else:
            reporter.add_error("❌ Error: File was not created")
        
    except Exception as e:
        reporter.add_error(f"DNA derived extension creation failed: {str(e)}")
        print(f"❌ DNA derived extension creation failed: {str(e)}")
        import traceback
        traceback.print_exc() 