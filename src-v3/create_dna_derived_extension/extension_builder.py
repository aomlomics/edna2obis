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
        
        # Define desired columns for DNA derived extension in output order
        # Unit columns are now generated and inserted dynamically.
        DESIRED_DNA_DERIVED_COLUMNS = [
            'occurrenceID', 'eventID', 'source_mat_id', 'samp_name','env_broad_scale', 'env_local_scale', 'env_medium', 
            'samp_vol_we_dna_ext', 'samp_collect_device', 'samp_mat_process', 'size_frac', 'samp_size',
            'concentration', 'lib_layout', 'seq_meth', 'nucl_acid_ext', 'target_gene', 
            'target_subfragment', 'pcr_primer_forward', 'pcr_primer_reverse', 
            'pcr_primer_name_forward', 'pcr_primer_name_reverse', 'pcr_primer_reference', 
            'pcr_cond', 'nucl_acid_amp', 'ampliconSize', 'otu_seq_comp_appr', 'otu_db', 
            'DNA_sequence', 'otu_class_appr'
        ]
        
        reporter.add_text(f"DNA derived extension will have {len(DESIRED_DNA_DERIVED_COLUMNS)} columns (plus dynamically added unit columns).")
        
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
        
        # Perform the merge. Pandas will add suffixes _x and _y for conflicting columns like 'assay_name'.
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
            suffixes=('', '_exp') # Suffix to identify and remove duplicate columns
        )

        # Clean up any duplicated columns from the merge (e.g., samp_name_exp) and the redundant lib_id
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

            # Strategy 2: If no _unit column, look up the unit in the FAIRe checklist
            elif not checklist_df.empty:
                checklist_row = checklist_df[checklist_df['term_name'] == field]
                if not checklist_row.empty and 'unit' in checklist_row.columns:
                    unit_value = checklist_row['unit'].iloc[0]
                    # Only create the column if the unit value from the checklist is meaningful
                    if pd.notna(unit_value) and str(unit_value).strip():
                        dna_derived_df_final[dwc_unit_col] = unit_value
                        created_unit_columns[field] = dwc_unit_col
                        reporter.add_text(f"✓ Created '{dwc_unit_col}' using value ('{unit_value}') from FAIRe checklist for field '{field}'.")
        
        # Handle projectMetadata in long format using the correct assay_name column
        project_fields_needed = [
            'lib_layout', 'instrument', 'target_gene', 'target_subfragment', 
            'pcr_primer_forward', 'pcr_primer_reverse', 'pcr_primer_name_forward', 
            'pcr_primer_name_reverse', 'pcr_primer_reference_forward', 'pcr_cond', 
            'nucl_acid_amp', 'ampliconSize'
        ]
        
        # For each field, add it to our dataframe
        reporter.add_text("Processing projectMetadata fields...")
        for field in project_fields_needed:
            field_row = data['projectMetadata'][data['projectMetadata']['term_name'] == field]
            
            if not field_row.empty:
                field_row = field_row.iloc[0]
                project_level_val = field_row.get('project_level')
                
                # Check if a valid project_level value exists
                if pd.notna(project_level_val) and str(project_level_val).strip() and str(project_level_val).lower() != 'nan':
                    dna_derived_df_final[field] = project_level_val
                else:
                    # Use assay-specific values
                    assay_value_map = {assay: field_row.get(assay) for assay in dna_derived_df_final['assay_name'].unique()}
                    dna_derived_df_final[field] = dna_derived_df_final['assay_name'].map(assay_value_map)
            else:
                dna_derived_df_final[field] = None
        
        # Add analysisMetadata fields using occurrence mapping
        if all_processed_occurrence_dfs:
            analysis_run_order = list(params['datafiles'].keys())
            analysis_fields_needed = ['otu_seq_comp_appr', 'otu_db', 'otu_clust_tool']
            
            # Initialize columns
            for field in analysis_fields_needed:
                dna_derived_df_final[field] = None
            
            current_row_start = 0
            for i, df in enumerate(all_processed_occurrence_dfs):
                if i >= len(analysis_run_order):
                    break
                    
                analysis_run_name = analysis_run_order[i]
                
                # Find which assay this analysis run belongs to
                assay_key = None
                for assay_name, analysis_runs in data['analysis_data_by_assay'].items():
                    if analysis_run_name in analysis_runs:
                        assay_key = assay_name
                        break
                
                if assay_key:
                    current_row_end = current_row_start + len(df)
                    
                    # Get analysis metadata for this run
                    if analysis_run_name in data['analysis_data_by_assay'][assay_key]:
                        analysis_df = data['analysis_data_by_assay'][assay_key][analysis_run_name]
                        
                        for field in analysis_fields_needed:
                            field_row = analysis_df[analysis_df['term_name'] == field]
                            if not field_row.empty:
                                field_value = field_row.iloc[0]['values']
                                # Apply to the specific row range
                                dna_derived_df_final.loc[current_row_start:current_row_end-1, field] = field_value
                    
                    current_row_start = current_row_end
        
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
        
        # The dwc_data is a dictionary containing DataFrames
        dna_derived_mapping = dwc_data['dnaDerived']
        
        # Create mapping dictionary from FAIRe names to Darwin Core names
        field_mapping = {}
        for dwc_name, row in dna_derived_mapping.iterrows():
            faire_name = row['FAIRe_term']
            if pd.notna(faire_name) and pd.notna(dwc_name):
                field_mapping[faire_name] = faire_name # This line was wrong before
        
        # Apply the mappings
        # This was dna_derived_df_final before, but sequences need to be added first
        dna_derived_renamed = dna_derived_with_sequences.copy()
        for faire_name, dwc_name in field_mapping.items():
            if faire_name in dna_derived_renamed.columns:
                dna_derived_renamed = dna_derived_renamed.rename(columns={faire_name: dwc_name})
        
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
            if col not in dna_derived_renamed.columns:
                dna_derived_renamed[col] = pd.NA
        
        # Select only the required columns in the correct order
        dna_derived_extension = dna_derived_renamed[final_ordered_columns]
        
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