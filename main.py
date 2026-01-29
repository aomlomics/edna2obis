#!/usr/bin/env python3
"""
edna2obis Main Script
Converts eDNA sequence data from FAIRe NOAA format to Darwin Core for OBIS and GBIF submission
"""

import os
import sys
from pathlib import Path
import yaml
import traceback
import warnings
import xml.etree.ElementTree as ET

# Add the src-v3 directory to Python path for imports
sys.path.insert(0, "src-v3")

# CLI UI
from cli_output.cli_ui import (
    console, print_header, print_usage, print_separator,
    silence_output, silence_stdout, silence_stdouterr
)

# Import required libraries
import numpy as np
import pandas as pd

# Import custom modules from src-v3
from html_reporter import HTMLReporter

# Import modules from subdirectories
from create_occurrence_core.occurrence_builder import create_occurrence_core
from create_dna_derived_extension.extension_builder import create_dna_derived_extension
from taxonomic_assignment.taxa_assignment_manager import assign_taxonomy
from create_eMoF.eMoF_builder import create_emof_table
# from create_EML.EML_builder import create_eml_file


def load_config(config_path="config.yaml"):
    """Load configuration from YAML file and convert to params dict structure"""
    try:
        if not os.path.exists(config_path):
            raise FileNotFoundError(f"Config file not found: {config_path}")
        
        with open(config_path, 'r', encoding='utf-8') as f:
            config = yaml.safe_load(f)
        
        params = {}
        
        # Metadata source: one boolean chooses which paths to use.
        # true = read from single Excel file (excel_file + sheet names)
        # false = read from separate TSV files (sampleMetadata_file, experimentRunMetadata_file, projectMetadata_file)
        # You can have both sets of paths in the config; this parameter decides which is used.
        combined_metadata_file = config.get('combined_metadata_file', True)
        params['use_excel'] = bool(combined_metadata_file)
        
        if params['use_excel']:
            params['excel_file'] = config.get('excel_file')
            params['sampleMetadata'] = config.get('sampleMetadata', 'sampleMetadata')
            params['experimentRunMetadata'] = config.get('experimentRunMetadata', 'experimentRunMetadata')
            params['projectMetadata'] = config.get('projectMetadata', 'projectMetadata')
            if not params['excel_file']:
                raise ValueError(
                    "combined_metadata_file is true but 'excel_file' is missing or empty. "
                    "Provide excel_file path, or set combined_metadata_file: false and provide TSV paths."
                )
        else:
            params['sampleMetadata_file'] = config.get('sampleMetadata_file')
            params['experimentRunMetadata_file'] = config.get('experimentRunMetadata_file')
            params['projectMetadata_file'] = config.get('projectMetadata_file')
            missing = [k for k, v in [
                ('sampleMetadata_file', params['sampleMetadata_file']),
                ('experimentRunMetadata_file', params['experimentRunMetadata_file']),
                ('projectMetadata_file', params['projectMetadata_file']),
            ] if not v]
            if missing:
                raise ValueError(
                    f"combined_metadata_file is false but required paths are missing: {missing}. "
                    "Provide sampleMetadata_file, experimentRunMetadata_file, and projectMetadata_file."
                )
        
        params['FAIRe_NOAA_checklist'] = config.get('FAIRe_NOAA_checklist') # Made optional
        params['datafiles'] = config['datafiles']
        params['skip_columns'] = config.get('skip_columns', [])
        
        # Handle locationID components
        params['locationID_components'] = config.get('locationID_components', ['line_id', 'station_id'])

        # Handle control sample detection
        if 'control_sample_detection' in config and config['control_sample_detection']:
            params['control_sample_column'] = config['control_sample_detection'].get('column_name', 'samp_category')
            params['control_sample_values'] = config['control_sample_detection'].get('control_values', [])
            params['skip_sample_types'] = config['control_sample_detection'].get('control_values', [])
        else:
            params['control_sample_column'] = 'samp_category'
            params['skip_sample_types'] = config.get('skip_sample_types', [])
            params['control_sample_values'] = config.get('skip_sample_types', [])
        
        params['taxonomic_api_source'] = config['taxonomic_api_source']
        params['assays_to_skip_species_match'] = config.get('assays_to_skip_species_match', []) # Made optional
        
        # Optional parameters with defaults
        params['worms_n_proc'] = config.get('worms_n_proc', 0)
        params['gbif_n_proc'] = config.get('gbif_n_proc', 0)
        params['gbif_match_limit'] = config.get('gbif_match_limit', 3)
        params['output_dir'] = config.get('output_dir', "processed-v3/")
        
        # eMoF options
        params['emof_enabled'] = config.get('emof_enabled', True)
        params['emof_template_path'] = config.get('emof_template_path', 'raw-v3/eMoF Fields edna2obis .xlsx')
        
        # EML options (metadata loaded separately from EML_config.yaml)
        params['eml_enabled'] = config.get('eml_enabled', False)
        params['eml_config_path'] = config.get('eml_config_path', 'EML_config.yaml')
        
        # Local reference database parameters  
        params['use_local_reference_database'] = config.get('use_local_reference_database', False)
        params['local_reference_database_path'] = config.get('local_reference_database_path', '')
        
        # Output splitting by short_name (cruise/expedition)
        params['split_output_by_short_name'] = config.get('split_output_by_short_name', False)
        
        # --- NEW: Metadata format switcher ---
        params['metadata_format'] = config.get('metadata_format', 'NOAA').upper()
        
        # --- NEW: Run name for config saving and report naming ---
        params['edna2obis_run_name'] = config.get('edna2obis_run_name', 'unnamed_run')

        return params
        
    except Exception as e:
        raise Exception(f"Error loading config: {str(e)}")


def setup_pandas_display():
    """Set pandas display options"""
    pd.set_option('display.max_colwidth', 150)
    pd.set_option('display.max_columns', 50)


def save_config_for_run(params, reporter, config_path):
    """Save the config file used for this run with the run name"""
    try:
        run_name = params.get('edna2obis_run_name', 'unnamed_run')
        output_dir = params.get('output_dir', "processed-v3/")
        
        # Create the output directory if it doesn't exist
        os.makedirs(output_dir, exist_ok=True)
        
        # Define the config filename with the run name
        config_filename = f"config_{run_name}.yaml"
        config_path_to_save = os.path.join(output_dir, config_filename)
        
        # Read the original config file that was used for the run
        with open(config_path, 'r', encoding='utf-8') as f:
            config_content = f.read()
        
        # Write the config to the output directory with the run name
        with open(config_path_to_save, 'w', encoding='utf-8') as f:
            f.write(config_content)
        
        reporter.add_success(f"Configuration file saved as: {config_filename}")
        reporter.add_text(f"Saved config file: {config_path_to_save}")
        
        return config_path_to_save
        
    except Exception as e:
        error_msg = f"Failed to save config file: {str(e)}"
        reporter.add_warning(error_msg)
        return None


def split_output_files_by_short_name(params, data, reporter):
    """
    Split output files into separate subfolders based on the 'short_name' column in sampleMetadata.
    
    This allows users to submit data cruise-by-cruise to databases that require separate submissions
    connected by a parent project.
    
    Creates:
        - {output_dir}/{short_name}/ subfolder for each unique short_name
        - Filtered files with suffix: occurrence_core_{api}_{short_name}.csv, etc.
    
    Splits: occurrence_core, dna_derived_extension, eMoF
    Does NOT split: HTML report, taxa_assignment_INFO, config file, eml.xml
    """
    reporter.add_section("Splitting Output Files by short_name", level=2)
    
    try:
        output_dir = params.get('output_dir', 'processed-v3/')
        api_choice = params.get('taxonomic_api_source', 'WoRMS').lower()
        
        # Check if short_name column exists in sampleMetadata
        if 'sampleMetadata' not in data or data['sampleMetadata'].empty:
            reporter.add_error("Cannot split by short_name: sampleMetadata is missing or empty.")
            raise ValueError("sampleMetadata is missing or empty")
        
        sm_df = data['sampleMetadata']
        
        if 'short_name' not in sm_df.columns:
            reporter.add_error(
                "Cannot split by short_name: The 'short_name' column was not found in sampleMetadata. "
                "Please add a 'short_name' column to your sampleMetadata sheet, or set 'split_output_by_short_name: false' in your config."
            )
            raise ValueError("'short_name' column not found in sampleMetadata")
        
        # Check for missing short_name values
        missing_short_name = sm_df['short_name'].isna() | (sm_df['short_name'].astype(str).str.strip() == '')
        if missing_short_name.any():
            missing_samples = sm_df.loc[missing_short_name, 'samp_name'].tolist()
            reporter.add_error(
                f"Cannot split by short_name: {len(missing_samples)} sample(s) have missing or empty 'short_name' values. "
                f"All samples must have a short_name when 'split_output_by_short_name' is enabled."
            )
            reporter.add_list(missing_samples[:20], "Samples with missing short_name (first 20):")
            raise ValueError(f"{len(missing_samples)} samples have missing short_name values")
        
        # Get unique short_names
        unique_short_names = sm_df['short_name'].astype(str).str.strip().unique().tolist()
        reporter.add_text(f"Found {len(unique_short_names)} unique short_name value(s): {', '.join(unique_short_names)}")
        
        # Build lookup: samp_name -> short_name
        samp_to_short_name = dict(zip(
            sm_df['samp_name'].astype(str).str.strip(),
            sm_df['short_name'].astype(str).str.strip()
        ))
        
        # Build lookup: lib_id (eventID) -> short_name via experimentRunMetadata
        event_to_short_name = {}
        if 'experimentRunMetadata' in data and not data['experimentRunMetadata'].empty:
            erm_df = data['experimentRunMetadata']
            if 'lib_id' in erm_df.columns and 'samp_name' in erm_df.columns:
                for _, row in erm_df.iterrows():
                    lib_id = str(row['lib_id']).strip()
                    samp_name = str(row['samp_name']).strip()
                    if samp_name in samp_to_short_name:
                        event_to_short_name[lib_id] = samp_to_short_name[samp_name]
        
        reporter.add_text(f"Built eventID -> short_name lookup with {len(event_to_short_name)} entries.")
        
        # Define files to split
        files_to_split = [
            (f'occurrence_core_{api_choice}.csv', 'eventID'),
            ('dna_derived_extension.csv', 'eventID'),
        ]
        
        # Add eMoF if it exists
        emof_path = os.path.join(output_dir, 'eMoF.csv')
        if os.path.exists(emof_path):
            files_to_split.append(('eMoF.csv', 'eventID'))
        
        # Process each short_name
        for short_name in unique_short_names:
            reporter.add_text(f"<h4>Processing short_name: '{short_name}'</h4>")
            
            # Create subfolder
            short_name_dir = os.path.join(output_dir, short_name)
            os.makedirs(short_name_dir, exist_ok=True)
            reporter.add_text(f"Created subfolder: {short_name_dir}")
            
            # Get eventIDs belonging to this short_name
            event_ids_for_short_name = {
                event_id for event_id, sn in event_to_short_name.items() 
                if sn == short_name
            }
            reporter.add_text(f"Found {len(event_ids_for_short_name)} eventID(s) for this short_name.")
            
            # Filter and save each file
            for filename, filter_column in files_to_split:
                source_path = os.path.join(output_dir, filename)
                
                if not os.path.exists(source_path):
                    reporter.add_text(f"  ⚠️ Skipping {filename}: file not found")
                    continue
                
                try:
                    df = pd.read_csv(source_path, low_memory=False)
                    
                    if filter_column not in df.columns:
                        reporter.add_warning(f"  ⚠️ Skipping {filename}: column '{filter_column}' not found")
                        continue
                    
                    # Filter rows where eventID belongs to this short_name
                    df[filter_column] = df[filter_column].astype(str).str.strip()
                    filtered_df = df[df[filter_column].isin(event_ids_for_short_name)]
                    
                    # Generate output filename with short_name suffix
                    base_name, ext = os.path.splitext(filename)
                    output_filename = f"{base_name}_{short_name}{ext}"
                    output_path = os.path.join(short_name_dir, output_filename)
                    
                    # Save filtered file
                    filtered_df.to_csv(output_path, index=False, encoding='utf-8-sig')
                    
                    reporter.add_text(f"  ✓ {output_filename}: {len(filtered_df):,} rows (from {len(df):,} total)")
                    
                except Exception as e:
                    reporter.add_warning(f"  ❌ Error processing {filename}: {str(e)}")
        
        reporter.add_success(f"Successfully split output files into {len(unique_short_names)} subfolder(s).")
        
        # Summary
        summary_items = []
        for short_name in unique_short_names:
            short_name_dir = os.path.join(output_dir, short_name)
            if os.path.exists(short_name_dir):
                files_in_dir = [f for f in os.listdir(short_name_dir) if f.endswith('.csv')]
                summary_items.append(f"<strong>{short_name}/</strong>: {len(files_in_dir)} file(s)")
        reporter.add_list(summary_items, "Output folders created:")
        
    except Exception as e:
        reporter.add_error(f"Failed to split output files by short_name: {str(e)}")
        raise


def load_project_data(params, reporter):
    """Load project, sample, experimentRun, and analysis data from FAIRe Excel file or TSV files"""
    reporter.add_section("Loading Project Data", level=2)
    
    try:
        use_excel = params.get('use_excel', True)
        analysis_data_by_assay = {}
        
        if use_excel:
            # --- EXCEL FORMAT (original implementation) ---
            # Discover all sheets in the Excel file
            excel = pd.ExcelFile(params['excel_file'])
            all_sheets = excel.sheet_names
            reporter.add_text(f"Found {len(all_sheets)} sheets in Excel file: {', '.join(all_sheets)}")

            # --- Auto-detect metadata format if misconfigured ---
            try:
                has_analysis_sheets = any(str(s).startswith('analysisMetadata') for s in all_sheets)
                current_fmt = params.get('metadata_format', 'NOAA')
                if has_analysis_sheets and current_fmt != 'NOAA':
                    reporter.add_warning("Detected 'analysisMetadata' sheets in Excel. Switching metadata_format to 'NOAA' for this run.")
                    params['metadata_format'] = 'NOAA'
                elif (not has_analysis_sheets) and current_fmt == 'NOAA':
                    reporter.add_warning("No 'analysisMetadata' sheets detected. Switching metadata_format to 'GENERIC' for this run.")
                    params['metadata_format'] = 'GENERIC'
            except Exception as _fmt_e:
                reporter.add_warning(f"Could not auto-detect metadata format: {_fmt_e}")

            # --- FAIRe-NOAA FORMAT ---
            # Uses separate analysisMetadata sheets for each run.
            if params['metadata_format'] == 'NOAA':
                analysis_sheets = [sheet for sheet in all_sheets if sheet.startswith('analysisMetadata')]
                reporter.add_text(f"Found {len(analysis_sheets)} analysis metadata sheets: {', '.join(analysis_sheets)}")

                for sheet_name in analysis_sheets:
                    analysis_df = pd.read_excel(params['excel_file'], sheet_name)
                    
                    assay_name = str(analysis_df.iloc[1, 3])  # Excel cell D3
                    analysis_run_name = str(analysis_df.iloc[2, 3])  # Excel cell D4
                    
                    reporter.add_text(f"Processing sheet '{sheet_name}': assay '{assay_name}', run '{analysis_run_name}'")
                    
                    if assay_name not in analysis_data_by_assay:
                        analysis_data_by_assay[assay_name] = {}
                    analysis_data_by_assay[assay_name][analysis_run_name] = analysis_df

            # --- GENERIC FAIRe FORMAT ---
            # Synthesizes analysis metadata from projectMetadata sheet.
            elif params['metadata_format'] == 'GENERIC':
                reporter.add_text("Using GENERIC FAIRe format. Synthesizing analysis metadata from projectMetadata.")
                
                project_meta_df = pd.read_excel(params['excel_file'], params['projectMetadata'], index_col=None, na_values=[""], comment="#")
                
                # --- Build a map from assay_name to its column header (e.g., 'ssu16s...': 'assay1') ---
                assay_map = {}
                assay_name_row = project_meta_df[project_meta_df['term_name'] == 'assay_name']
                if not assay_name_row.empty:
                    assay_name_series = assay_name_row.iloc[0]
                    for col_name, assay_name_val in assay_name_series.items():
                        if col_name in ['term_name', 'project_level'] or pd.isna(assay_name_val):
                            continue
                        assay_map[assay_name_val] = col_name
                
                reporter.add_text(f"Built assay-to-column map from projectMetadata: {assay_map}")

                for run_name, run_details in params['datafiles'].items():
                    assay_name = run_details.get('assay_name')
                    if not assay_name:
                        reporter.add_warning(f"For analysis run '{run_name}', 'assay_name' is missing in config.yaml. Cannot process analysis-specific metadata.")
                        continue

                    # --- Use the map to find the correct column name ---
                    if assay_name not in assay_map:
                        reporter.add_warning(f"Assay '{assay_name}' for run '{run_name}' not found in the 'assay_name' row of projectMetadata. Please check for typos. Available assays found: {list(assay_map.keys())}")
                        continue
                    
                    # Get the actual column name (e.g., 'assay1')
                    original_assay_col_name = assay_map[assay_name]
                    
                    # Synthesize the analysis DF
                    rows = []
                    for _, row in project_meta_df.iterrows():
                        term_name = row['term_name']
                        # Assay-specific value takes precedence over project-level value
                        assay_specific_value = row[original_assay_col_name]
                        project_level_value = row['project_level']
                        
                        final_value = assay_specific_value if pd.notna(assay_specific_value) else project_level_value
                        rows.append({'term_name': term_name, 'values': final_value})
                    
                    synthesized_df = pd.DataFrame(rows)
                    
                    # Build the nested dictionary structure
                    if assay_name not in analysis_data_by_assay:
                        analysis_data_by_assay[assay_name] = {}
                    analysis_data_by_assay[assay_name][run_name] = synthesized_df
                    reporter.add_text(f"Synthesized analysis metadata for run '{run_name}' using assay '{assay_name}'.")

            # Load the main data sheets
            data = pd.read_excel(
                params['excel_file'],
                [params['projectMetadata'], params['sampleMetadata'], params['experimentRunMetadata']],
                index_col=None, na_values=[""], comment="#"
            )
            
            # Rename keys to standard terms
            data['sampleMetadata'] = data.pop(params['sampleMetadata'])
            data['experimentRunMetadata'] = data.pop(params['experimentRunMetadata'])
            data['projectMetadata'] = data.pop(params['projectMetadata'])
            
            # For backward compatibility
            if params['metadata_format'] == 'NOAA':
                analysis_sheets = [sheet for sheet in all_sheets if sheet.startswith('analysisMetadata')]
                if analysis_sheets:
                    data['analysisMetadata'] = pd.read_excel(params['excel_file'], analysis_sheets[0])
        
        else:
            # --- TSV FORMAT (new implementation) ---
            reporter.add_text("Using TSV file format. Loading separate TSV files for each metadata sheet.")
            
            # Load main metadata sheets from TSV files
            data = {}
            
            # Load projectMetadata
            project_meta_path = params['projectMetadata_file']
            reporter.add_text(f"Loading projectMetadata from: {project_meta_path}")
            data['projectMetadata'] = pd.read_csv(project_meta_path, sep='\t', na_values=[""], comment="#", low_memory=False)
            
            # Load sampleMetadata
            sample_meta_path = params['sampleMetadata_file']
            reporter.add_text(f"Loading sampleMetadata from: {sample_meta_path}")
            data['sampleMetadata'] = pd.read_csv(sample_meta_path, sep='\t', na_values=[""], comment="#", low_memory=False)
            
            # Load experimentRunMetadata
            exp_run_meta_path = params['experimentRunMetadata_file']
            reporter.add_text(f"Loading experimentRunMetadata from: {exp_run_meta_path}")
            data['experimentRunMetadata'] = pd.read_csv(exp_run_meta_path, sep='\t', na_values=[""], comment="#", low_memory=False)
            
            # Handle analysis metadata based on format
            if params['metadata_format'] == 'NOAA':
                # For NOAA format with TSV, expect analysisMetadata_file in each datafiles entry
                reporter.add_text("Using NOAA format with TSV files. Loading analysisMetadata from separate TSV files.")
                
                for run_name, run_details in params['datafiles'].items():
                    assay_name = run_details.get('assay_name')
                    analysis_meta_file = run_details.get('analysisMetadata_file')
                    
                    if not assay_name:
                        reporter.add_warning(f"For analysis run '{run_name}', 'assay_name' is missing in config.yaml. Skipping analysis metadata.")
                        continue
                    
                    if not analysis_meta_file:
                        reporter.add_warning(f"For analysis run '{run_name}', 'analysisMetadata_file' is missing in config.yaml. Skipping analysis metadata for this run.")
                        continue
                    
                    reporter.add_text(f"Loading analysisMetadata for run '{run_name}' (assay: '{assay_name}') from: {analysis_meta_file}")
                    
                    # Read the TSV file - TSV format should have term_name and values columns
                    analysis_df = pd.read_csv(analysis_meta_file, sep='\t', na_values=[""], comment="#", low_memory=False)
                    
                    # Validate expected columns exist
                    if 'term_name' not in analysis_df.columns:
                        reporter.add_warning(f"analysisMetadata file '{analysis_meta_file}' missing 'term_name' column. Found columns: {list(analysis_df.columns)}")
                        continue
                    
                    # The 'values' column might be named differently in TSV, check common names
                    value_col = None
                    for col in ['values', 'value', 'Values', 'Value']:
                        if col in analysis_df.columns:
                            value_col = col
                            break
                    
                    if not value_col:
                        reporter.add_warning(f"analysisMetadata file '{analysis_meta_file}' missing 'values' column. Found columns: {list(analysis_df.columns)}")
                        continue
                    
                    # Rename to standard 'values' column name
                    if value_col != 'values':
                        analysis_df = analysis_df.rename(columns={value_col: 'values'})
                    
                    if assay_name not in analysis_data_by_assay:
                        analysis_data_by_assay[assay_name] = {}
                    analysis_data_by_assay[assay_name][run_name] = analysis_df
                    reporter.add_text(f"Successfully loaded analysisMetadata for run '{run_name}' (assay: '{assay_name}')")
            
            elif params['metadata_format'] == 'GENERIC':
                # For GENERIC format with TSV, synthesize from projectMetadata (same as Excel)
                reporter.add_text("Using GENERIC FAIRe format with TSV files. Synthesizing analysis metadata from projectMetadata.")
                
                project_meta_df = data['projectMetadata']
                
                # --- Build a map from assay_name to its column header (e.g., 'ssu16s...': 'assay1') ---
                assay_map = {}
                assay_name_row = project_meta_df[project_meta_df['term_name'] == 'assay_name']
                if not assay_name_row.empty:
                    assay_name_series = assay_name_row.iloc[0]
                    for col_name, assay_name_val in assay_name_series.items():
                        if col_name in ['term_name', 'project_level'] or pd.isna(assay_name_val):
                            continue
                        assay_map[assay_name_val] = col_name
                
                reporter.add_text(f"Built assay-to-column map from projectMetadata: {assay_map}")

                for run_name, run_details in params['datafiles'].items():
                    assay_name = run_details.get('assay_name')
                    if not assay_name:
                        reporter.add_warning(f"For analysis run '{run_name}', 'assay_name' is missing in config.yaml. Cannot process analysis-specific metadata.")
                        continue

                    # --- Use the map to find the correct column name ---
                    if assay_name not in assay_map:
                        reporter.add_warning(f"Assay '{assay_name}' for run '{run_name}' not found in the 'assay_name' row of projectMetadata. Please check for typos. Available assays found: {list(assay_map.keys())}")
                        continue
                    
                    # Get the actual column name (e.g., 'assay1')
                    original_assay_col_name = assay_map[assay_name]
                    
                    # Synthesize the analysis DF
                    rows = []
                    for _, row in project_meta_df.iterrows():
                        term_name = row['term_name']
                        # Assay-specific value takes precedence over project-level value
                        assay_specific_value = row[original_assay_col_name]
                        project_level_value = row['project_level']
                        
                        final_value = assay_specific_value if pd.notna(assay_specific_value) else project_level_value
                        rows.append({'term_name': term_name, 'values': final_value})
                    
                    synthesized_df = pd.DataFrame(rows)
                    
                    # Build the nested dictionary structure
                    if assay_name not in analysis_data_by_assay:
                        analysis_data_by_assay[assay_name] = {}
                    analysis_data_by_assay[assay_name][run_name] = synthesized_df
                    reporter.add_text(f"Synthesized analysis metadata for run '{run_name}' using assay '{assay_name}'.")
        
        # Add analysis data to main data dictionary (common for both formats)
        data['analysis_data_by_assay'] = analysis_data_by_assay
        
        # Report summary
        reporter.add_success(f"Successfully loaded project data:")
        reporter.add_list([
            f"Sample metadata: {len(data['sampleMetadata'])} rows",
            f"Experiment run metadata: {len(data['experimentRunMetadata'])} rows", 
            f"Project metadata: {len(data['projectMetadata'])} rows",
            f"Analysis assays: {len(analysis_data_by_assay)}"
        ])
        
        # Show sample of data
        reporter.add_text("<h4>sampleMetadata preview:</h4>")
        reporter.add_text("Contextual data about the samples collected, such as when it was collected, where it was collected from, what kind of sample it is, and what were the properties of the environment or experimental condition from which the sample was taken. Each row is a distinct sample, or Event. Most of this information is recorded during sample collection. This sheet contains terms from the FAIRe NOAA data template.")
        reporter.add_dataframe(data['sampleMetadata'], "sampleMetadata (first 5 rows)", max_rows=5)
        
        reporter.add_text("<h4>experimentRunMetadata preview:</h4>")
        reporter.add_text("Contextual data about how the samples were prepared for sequencing. Includes how they were extracted, what amplicon was targeted, how they were sequenced. Each row is a separate sequencing library preparation, distinguished by a unique lib_id.")
        reporter.add_dataframe(data['experimentRunMetadata'], "experimentRunMetadata (first 2 rows)", max_rows=2)
        
        return data
        
    except Exception as e:
        error_msg = f"Failed to load project data: {str(e)}"
        reporter.add_error(error_msg)
        raise Exception(error_msg)


def load_asv_data(params, reporter):
    """Load ASV data for each analysis run"""
    reporter.add_section("Loading ASV Data", level=2)
    
    raw_data_tables = {}
    
    try:
        for analysis_run_name, file_paths in params['datafiles'].items():
            reporter.add_text(f"Processing analysis run: {analysis_run_name}")
            raw_data_tables[analysis_run_name] = {}
            
            # Load taxonomy file
            if 'taxonomy_file' in file_paths:
                tax_path = file_paths['taxonomy_file']
                try:
                    tax_df = pd.DataFrame()
                    file_extension = os.path.splitext(tax_path)[1].lower()

                    if file_extension == '.xlsx':
                        # Heuristic sheet detection: choose the sheet with expected taxonomy columns
                        xls = pd.ExcelFile(tax_path)
                        candidate_sheets = xls.sheet_names
                        # Avoid obviously non-data sheets by name
                        blacklist = {'readme', 'instructions', 'checklist', 'metadata', 'requirements'}
                        expected_cols = {
                            'seq_id', 'featureid', 'feature id', 'otu id', '#otuid',
                            'taxonomy', 'taxon', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'scientificname', 'scientificnameauthorship',
                            'dna_sequence', 'sequence', 'verbatimidentification', 'taxonrank', 'taxonid'
                        }
                        best_sheet = None
                        best_score = -1
                        best_header_row = 0
                        
                        # Find the best sheet by checking all candidate sheets
                        for sheet in candidate_sheets:
                            sheet_l = str(sheet).lower().strip()
                            if any(b in sheet_l for b in blacklist):
                                continue
                            try:
                                tmp = pd.read_excel(tax_path, sheet_name=sheet, nrows=50, header=None)
                            except Exception:
                                continue
                            # Try to find a row that looks like a header
                            for i in range(min(10, len(tmp))):
                                row_vals = [str(v).lower().strip() for v in tmp.iloc[i].tolist()]
                                row_set = set(row_vals)
                                score = len(row_set & expected_cols)
                                if ('seq_id' in row_vals or 'featureid' in row_vals or 'otu id' in row_vals or '#otuid' in row_vals) and score > best_score:
                                    best_sheet = sheet
                                    best_header_row = i
                                    best_score = score
                        
                        # If no sheet was found, use the first sheet
                        if best_sheet is None:
                            best_sheet = candidate_sheets[0]
                            best_header_row = 0
                        
                        # Read the Excel file with the determined sheet and header row
                        tax_df = pd.read_excel(tax_path, sheet_name=best_sheet, header=best_header_row)
                        # Drop any completely empty columns that may have been created by reading
                        tax_df = tax_df.loc[:, ~tax_df.columns.astype(str).str.match(r'^Unnamed', na=False)]
                        reporter.add_text(f"Note: Reading '{tax_path}' sheet '{best_sheet}' with header at row index {best_header_row}. Auto-detected based on expected taxonomy columns.")
                    
                    else:
                        # Logic for text-based files (.txt, .tsv, .csv)
                        separator = ',' if file_extension == '.csv' else '\t'
                        
                        # Determine how many commented lines to skip by reading the start of the file
                        skiprows = 0
                        with open(tax_path, 'r', encoding='utf-8', errors='ignore') as f:
                            for line in f:
                                if line.startswith('#'):
                                    skiprows += 1
                                else:
                                    break
                        
                        # Load the data, skipping any initial comment lines
                        try:
                            tax_df = pd.read_csv(tax_path, sep=separator, skiprows=skiprows, low_memory=False, encoding='utf-8')
                        except UnicodeDecodeError:
                            try:
                                tax_df = pd.read_csv(tax_path, sep=separator, skiprows=skiprows, low_memory=False, encoding='latin-1')
                            except Exception:
                                # Final fallback: sniff delimiter, use python engine, ignore commented lines and bad lines
                                tax_df = pd.read_csv(tax_path, sep=None, engine='python', comment='#', on_bad_lines='skip', encoding='latin-1')

                    # --- Perform essential input normalization ---
                    
                    # 1. Standardize the feature ID column name (prefer explicit matches)
                    # Create a mapping of lowercased column names (with spaces and underscores removed) to original names
                    cols_lower = {str(c).lower().replace('_', '').replace(' ', '').replace('-', ''): c for c in tax_df.columns if isinstance(c, str)}
                    
                    # Look for feature ID column
                    if 'featureid' not in [str(c).lower() for c in tax_df.columns]:
                        # Try to find a column that matches common feature ID names
                        if 'seqid' in cols_lower:
                            tax_df.rename(columns={cols_lower['seqid']: 'featureid'}, inplace=True)
                        elif 'featureid' in cols_lower:
                            tax_df.rename(columns={cols_lower['featureid']: 'featureid'}, inplace=True)
                        elif 'otuid' in cols_lower:
                            tax_df.rename(columns={cols_lower['otuid']: 'featureid'}, inplace=True)
                        elif '#otuid' in cols_lower:
                            tax_df.rename(columns={cols_lower['#otuid']: 'featureid'}, inplace=True)
                        # If still not found, downstream will attempt to detect
                    
                    # 2. Ensure verbatimIdentification column for taxonomy
                    # Prefer explicit 'verbatimIdentification'; fallback to 'taxonomy', then 'taxon'
                    if 'verbatimIdentification' not in tax_df.columns:
                        cols_lower_exact = {str(c).lower(): c for c in tax_df.columns if isinstance(c, str)}
                        if 'verbatimidentification' in cols_lower_exact:
                            tax_df.rename(columns={cols_lower_exact['verbatimidentification']: 'verbatimIdentification'}, inplace=True)
                        elif 'taxonomy' in cols_lower_exact:
                            tax_df.rename(columns={cols_lower_exact['taxonomy']: 'verbatimIdentification'}, inplace=True)
                        elif 'taxon' in cols_lower_exact:
                            tax_df.rename(columns={cols_lower_exact['taxon']: 'verbatimIdentification'}, inplace=True)
                    
                    # 3. Standardize DNA sequence column
                    if 'dna_sequence' not in [str(c).lower() for c in tax_df.columns]:
                        cols_lower_exact = {str(c).lower(): c for c in tax_df.columns if isinstance(c, str)}
                        if 'sequence' in cols_lower_exact and 'dna_sequence' not in tax_df.columns:
                            tax_df.rename(columns={cols_lower_exact['sequence']: 'dna_sequence'}, inplace=True)
                    
                    # 4. Standardize species/scientific name column for downstream matching
                    if 'scientificName' not in tax_df.columns:
                        cols_lower_exact = {str(c).lower(): c for c in tax_df.columns if isinstance(c, str)}
                        if 'species' in cols_lower_exact and 'scientificname' not in cols_lower_exact:
                            tax_df.rename(columns={cols_lower_exact['species']: 'scientificName'}, inplace=True)

                    # 5. Normalize featureid values for safe joins
                    if 'featureid' in tax_df.columns:
                        tax_df['featureid'] = tax_df['featureid'].astype(str).str.strip()

                    raw_data_tables[analysis_run_name]['taxonomy'] = tax_df
                    tax_shape = tax_df.shape
                    reporter.add_success(f"Loaded taxonomy file: {tax_path} (shape: {tax_shape})")
                    try:
                        reporter.add_text(f"Taxonomy columns detected: {list(tax_df.columns)[:20]}")
                    except Exception:
                        pass

                except Exception as e:
                    reporter.add_error(f"Failed to load taxonomy file {tax_path}: {e}")
                    raise
            
            # Load abundance table file  
            # --- Robustly get the path for the abundance table, checking for all possible keys ---
            abundance_key = None
            possible_keys = ['abundance_table', 'occurrence_file', 'abundance_file']
            for key in possible_keys:
                if key in file_paths:
                    abundance_key = key
                    break
            
            if abundance_key:
                abundance_path = file_paths[abundance_key]
                try:
                    df_abundance = pd.DataFrame()
                    file_extension = os.path.splitext(abundance_path)[1].lower()

                    if file_extension == '.xlsx':
                        df_abundance = pd.read_excel(abundance_path, sheet_name=0)
                        reporter.add_text(f"Note: Reading '{abundance_path}' as an Excel file. Assumes header is on the first row.")
                    else:
                        # Logic for text-based files (.txt, .tsv, .csv)
                        separator = ',' if file_extension == '.csv' else '\t'

                        # 1. Check the first line to see if we need to skip it.
                        with open(abundance_path, 'r', encoding='utf-8', errors='ignore') as f:
                            first_line = f.readline()
                        
                        # The old format has a comment line, the new one starts with the header.
                        rows_to_skip = 1 if '# Constructed from biom file' in first_line else 0

                        # 2. Load the table, skipping the comment line only if it exists.
                        try:
                            df_abundance = pd.read_csv(abundance_path,
                                                       sep=separator,
                                                       skiprows=rows_to_skip,
                                                       header=0,
                                                       low_memory=False,
                                                       encoding='utf-8')
                        except UnicodeDecodeError:
                            try:
                                df_abundance = pd.read_csv(abundance_path,
                                                           sep=separator,
                                                           skiprows=rows_to_skip,
                                                           header=0,
                                                           low_memory=False,
                                                           encoding='latin-1')
                            except Exception:
                                df_abundance = pd.read_csv(abundance_path,
                                                           sep=None,
                                                           engine='python',
                                                           comment='#',
                                                           on_bad_lines='skip',
                                                           header=0,
                                                           encoding='latin-1')
                    
                    # 3. Standardize the first column name to 'featureid'.
                    if len(df_abundance.columns) > 0:
                        first_col_name = df_abundance.columns[0]
                        # Remove leading '#' if present
                        if isinstance(first_col_name, str) and first_col_name.startswith('#'):
                            clean_col_name = first_col_name[1:]
                            df_abundance.rename(columns={first_col_name: clean_col_name}, inplace=True)
                            first_col_name = clean_col_name
                        # Unconditionally standardize the first column to 'featureid' (generic abundance uses seq_id as first column)
                        df_abundance.rename(columns={first_col_name: 'featureid'}, inplace=True)
                        # Normalize keys
                        df_abundance['featureid'] = df_abundance['featureid'].astype(str).str.strip()

                        # --- NEW CRITICAL FIX: Clean and normalize ALL sample name headers ---
                        # This prevents join failures caused by whitespace or type issues from different file formats (.xlsx, .txt)
                        header_rename_map = {col: str(col).strip() for col in df_abundance.columns if col != 'featureid'}
                        df_abundance.rename(columns=header_rename_map, inplace=True)

                    raw_data_tables[analysis_run_name]['occurrence'] = df_abundance
                    abundance_shape = df_abundance.shape
                    reporter.add_success(f"Loaded abundance table file: {abundance_path} (shape: {abundance_shape})")
                except Exception as e:
                    reporter.add_error(f"Failed to load abundance table file {abundance_path}: {e}")
                    raise
        
        # Show preview of abundance data
        # Choose the first analysis run to inspect
        if raw_data_tables:
            first_analysis = list(raw_data_tables.keys())[0]
            if 'occurrence' in raw_data_tables[first_analysis]:
                # Show first 5 rows and first 20 columns like in notebook
                preview_df = raw_data_tables[first_analysis]['occurrence'].iloc[:5, :20]
                reporter.add_dataframe(preview_df, f"Preview of Abundance Table for '{first_analysis}'")
        
        reporter.add_success(f"Successfully loaded ASV data for {len(raw_data_tables)} analysis runs")
        return raw_data_tables
        
    except Exception as e:
        error_msg = f"Failed to load ASV data: {str(e)}"
        reporter.add_error(error_msg)
        raise Exception(error_msg)


def remove_control_samples(data, raw_data_tables, params, reporter):
    """Remove control/blank samples from all datasets"""
    reporter.add_section("Removing Control Samples", level=2)
    
    try:
        # Get control detection configuration
        control_column = params.get('control_sample_column', 'samp_category')
        control_values = params.get('control_sample_values', [])
        
        if not control_values:
            reporter.add_text("No control sample values specified - skipping control sample removal")
            return data, raw_data_tables
        
        # Check if the specified column exists
        if control_column not in data['sampleMetadata'].columns:
            reporter.add_text(f"⚠️ Warning: Control detection column '{control_column}' not found in sampleMetadata. Available columns: {list(data['sampleMetadata'].columns)}")
            reporter.add_text("Skipping control sample removal")
            return data, raw_data_tables
        
        # Identify samples to remove
        samps_to_remove = data['sampleMetadata'][control_column].isin(control_values)
        samples_to_drop = data['sampleMetadata']['samp_name'][samps_to_remove].astype(str).str.strip().tolist()
        
        reporter.add_text(f"Using column '{control_column}' to detect control samples")
        reporter.add_text(f"Looking for these control values: {control_values}")
        reporter.add_text(f"Found {len(samples_to_drop)} control/blank samples to remove")
        
        # Show the complete list of samples to be dropped
        if samples_to_drop:
            reporter.add_text("You can view the list of samples to be dropped below:")
            reporter.add_list(samples_to_drop, "Complete list of samples being removed:")
        
        # Remove from sampleMetadata
        data['sampleMetadata'] = data['sampleMetadata'][~samps_to_remove]
        reporter.add_text(f"Sample metadata shape after removal: {data['sampleMetadata'].shape}")
        
        # Check remaining sample categories
        remaining_categories = data['sampleMetadata'][control_column].unique()
        reporter.add_text(f"Check the {control_column} values left in your sampleMetadata.")
        reporter.add_text(f"Remaining {control_column} values: {', '.join(map(str, remaining_categories))}")
        
        # Remove from experimentRunMetadata
        prep_samps_to_remove = data['experimentRunMetadata']['samp_name'].isin(samples_to_drop)
        data['experimentRunMetadata'] = data['experimentRunMetadata'][~prep_samps_to_remove]
        reporter.add_text(f"Experiment run metadata shape after removal: {data['experimentRunMetadata'].shape}")
        
        # Remove from ASV abundance tables
        for analysis_run_name, tables_dict in raw_data_tables.items():
            if 'occurrence' in tables_dict:
                abundance_df = tables_dict['occurrence']
                original_cols = len(abundance_df.columns)
                
                cols_to_remove = [col for col in samples_to_drop if col in abundance_df.columns]
                if cols_to_remove:
                    raw_data_tables[analysis_run_name]['occurrence'] = abundance_df.drop(columns=cols_to_remove)
                    new_cols = len(raw_data_tables[analysis_run_name]['occurrence'].columns)
                    reporter.add_text(f"Analysis '{analysis_run_name}': removed {len(cols_to_remove)} columns ({original_cols} → {new_cols})")
        
        reporter.add_success("Successfully removed control samples from all datasets")
        return data, raw_data_tables
        
    except Exception as e:
        error_msg = f"Failed to remove control samples: {str(e)}"
        reporter.add_error(error_msg)
        raise Exception(error_msg)


def drop_all_na_columns(data, reporter):
    """Drop columns with all NAs - EXACT code from notebook cell 16-18"""
    reporter.add_section("Drop columns with all NAs", level=2)
    
    try:
        dropped = pd.DataFrame()

        for sheet_name in ['sampleMetadata', 'experimentRunMetadata']:
            # Safety check: ensure the sheet exists in data and is not empty
            if sheet_name in data and not data[sheet_name].empty:
                all_na_cols = data[sheet_name].columns[data[sheet_name].isnull().all(axis=0)]
                res = pd.Series(all_na_cols, name=sheet_name)
                dropped = pd.concat([dropped, res], axis=1)
            elif sheet_name not in data:
                reporter.add_text(f"FYI: Sheet '{sheet_name}' not found in 'data' dictionary. Cannot check for all-NA columns.")
            else: # In data but empty
                reporter.add_text(f"FYI: Sheet '{sheet_name}' is empty. Cannot check for all-NA columns.")
        
        # Show which columns have only NA values
        reporter.add_text("Which columns in each sheet have only NA values?")
        reporter.add_dataframe(dropped, "Columns with all NAs by sheet")
        
        # Drops all-NA columns from 'sampleMetadata' and 'experimentRunMetadata'
        # based on the 'dropped' DataFrame.
        sheets_to_clean = ['sampleMetadata', 'experimentRunMetadata']

        for sheet_name in sheets_to_clean:
            # Check if 'dropped' has a column for this sheet AND if that column lists any actual columns to drop
            if sheet_name in dropped.columns and not dropped[sheet_name].dropna().empty:
                cols_to_drop = list(dropped[sheet_name].dropna())
                
                if sheet_name in data:
                    reporter.add_text(f"Dropping from data['{sheet_name}']: {cols_to_drop}")
                    data[sheet_name].drop(columns=cols_to_drop, inplace=True, errors='ignore')
        
        reporter.add_success("Successfully dropped all-NA columns")
        return data
        
    except Exception as e:
        error_msg = f"Failed to drop all-NA columns: {str(e)}"
        reporter.add_error(error_msg)
        raise Exception(error_msg)


def drop_empty_analysis_rows(data, reporter):
    """Drop NA rows of each analysisMetadata sheet"""
    reporter.add_section("Now lets drop NA rows of each analysisMetadata sheet", level=2)
    
    try:
        # Identify and report rows with empty 'values' in analysisMetadata sheets
        analysis_rows_to_drop_info = {} 
        expected_value_col_name = 'values' # As per the FAIRe NOAA Excel sheet (column D header)

        reporter.add_text(f"Identifying rows with empty '{expected_value_col_name}' column in analysisMetadata sheets...")

        if 'analysis_data_by_assay' in data and data['analysis_data_by_assay']:
            for assay_name, analyses_dict in data['analysis_data_by_assay'].items():
                if not isinstance(analyses_dict, dict): continue
                for run_name, analysis_df in analyses_dict.items():
                    if not isinstance(analysis_df, pd.DataFrame) or analysis_df.empty: continue

                    # Ensure the 'values' column exists
                    if expected_value_col_name not in analysis_df.columns:
                        reporter.add_text(f"  Warning: Column '{expected_value_col_name}' not found in Assay: '{assay_name}', Run: '{run_name}'. Skipping.")
                        continue
                    
                    # Identify rows where the 'values' column is NA
                    empty_values_mask = analysis_df[expected_value_col_name].isnull()
                    rows_to_drop_indices = analysis_df[empty_values_mask].index.tolist()

                    if rows_to_drop_indices:
                        # Get the list of term names for the empty rows
                        terms_to_drop = analysis_df.loc[rows_to_drop_indices, 'term_name'].tolist() if 'term_name' in analysis_df.columns else []
                        
                        analysis_rows_to_drop_info[(assay_name, run_name)] = rows_to_drop_indices

                        reporter.add_text(f"For Analysis (Assay: '{assay_name}', Run: '{run_name}'), found {len(terms_to_drop)} empty rows:")
                        reporter.add_text(f"  {terms_to_drop}")

        if not analysis_rows_to_drop_info:
            reporter.add_text(f"No rows with an empty '{expected_value_col_name}' column were found in any analysisMetadata sheet.")
        else:
            reporter.add_text("Summary: The terms listed above will be dropped from their respective sheets.")
            
            # Perform the deletion of rows with empty 'values' from analysisMetadata sheets
            reporter.add_text(f"Removing identified rows with empty 'values' from analysisMetadata sheets...")
            for (assay_name, run_name), indices_to_drop in analysis_rows_to_drop_info.items():
                if indices_to_drop: 
                    try:
                        original_count = len(data['analysis_data_by_assay'][assay_name][run_name])
                        data['analysis_data_by_assay'][assay_name][run_name].drop(index=indices_to_drop, inplace=True)
                        new_count = len(data['analysis_data_by_assay'][assay_name][run_name])
                        reporter.add_text(f"  For Assay: '{assay_name}', Run: '{run_name}': Removed {original_count - new_count} row(s).")
                    except Exception as e:
                        reporter.add_text(f"  An error occurred while dropping rows for Assay: '{assay_name}', Run: '{run_name}': {e}")
            reporter.add_text("Finished removing rows with empty 'values' from analysisMetadata sheets.")
        
        reporter.add_success("Finished removing rows with empty 'values' from analysisMetadata sheets")
        return data
        
    except Exception as e:
        error_msg = f"Failed to drop empty analysis rows: {str(e)}"
        reporter.add_error(error_msg)
        raise Exception(error_msg)


def drop_some_na_columns(data, reporter):
    """Now lets check for columns or rows with SOME missing values"""
    reporter.add_section("Now lets check for columns or rows with SOME missing values", level=2)
    
    try:
        reporter.add_text("Now let's check which columns have missing values in some of the rows. These should be filled in on the Excel sheet with the appropriate term ('not applicable', 'missing', or 'not collected'). Alternatively, you can drop the column if it is not needed for submission to OBIS.")
        
        some = pd.DataFrame()

        # Sheets to check for columns with *some* NAs
        sheets_to_examine_for_some_na = ['sampleMetadata', 'experimentRunMetadata']

        for sheet_name in sheets_to_examine_for_some_na:
            if sheet_name in data and not data[sheet_name].empty:
                cols_with_some_na = data[sheet_name].columns[data[sheet_name].isnull().any(axis=0)]
                res = pd.Series(cols_with_some_na.tolist(), name=sheet_name) # .tolist() for cleaner Series
                some = pd.concat([some, res], axis=1)

        reporter.add_text("Columns with some NAs:")
        reporter.add_dataframe(some, "Columns with some missing values")
        
        reporter.add_text("Here I'm going to drop all the columns with some missing data, as I don't need them for submission to OBIS.")
        
        # Drop columns with some missing values, but preserve any column referenced in data_mapper.yaml
        sheets_to_clean = ['sampleMetadata', 'experimentRunMetadata']

        # Build set of columns to preserve from mapper
        try:
            preserve_cols = set()
            import yaml as _yaml_loader
            with open('data_mapper.yaml','r', encoding='utf-8') as f:
                mapper = _yaml_loader.safe_load(f) or {}
            # Consider both generic and default keys
            mapping_keys = []
            # We don't know the format here; load both sections if present
            mapping_keys.extend(['occurrence_core','dna_derived_extension','generic_occurrence_core','generic_dna_derived_extension'])
            for key in mapping_keys:
                section = mapper.get(key, {}) or {}
                for _dwc, info in section.items():
                    if isinstance(info, dict):
                        ft = info.get('faire_term')
                        src = info.get('source')
                        if ft and src in ['sampleMetadata','experimentRunMetadata']:
                            preserve_cols.add(str(ft).strip())
                            preserve_cols.add(f"{str(ft).strip()}_unit")
        except Exception:
            preserve_cols = set()

        for sheet_name in sheets_to_clean:
            if sheet_name in data and not data[sheet_name].empty: # Check if DataFrame exists and is not empty
                original_columns = data[sheet_name].columns.tolist() # Get column names before dropping

                # Compute droppable columns (some NAs) excluding preserved
                cols_with_some_na = data[sheet_name].columns[data[sheet_name].isnull().any(axis=0)]
                cols_to_drop = [c for c in cols_with_some_na if c not in preserve_cols]

                if cols_to_drop:
                    data[sheet_name].drop(columns=cols_to_drop, inplace=True, errors='ignore')
                
                current_columns = data[sheet_name].columns.tolist() # Get column names after dropping
                dropped_column_names = [col for col in original_columns if col not in current_columns]
                
                if dropped_column_names:
                    reporter.add_text(f"From data['{sheet_name}']: Dropped columns: {dropped_column_names}")
                else: # Optional: if no columns were dropped
                    reporter.add_text(f"No columns dropped from data['{sheet_name}']")
            elif sheet_name not in data:
                reporter.add_text(f"FYI: Sheet '{sheet_name}' not found in 'data'. No columns dropped.")
            else: # Sheet is in data but empty
                reporter.add_text(f"FYI: Sheet '{sheet_name}' is empty. No columns dropped.")
        
        reporter.add_success("Successfully dropped columns with some missing values")
        return data
        
    except Exception as e:
        error_msg = f"Failed to drop columns with some NAs: {str(e)}"
        reporter.add_error(error_msg)
        raise Exception(error_msg)


def load_darwin_core_mappings(params, reporter):
    """Load data dictionary Excel file - EXACT code from notebook cell 23"""
    reporter.add_section("Load data dictionary Excel file", level=2)
    
    try:
        reporter.add_text("Attempting to load Darwin Core mappings from data_mapper.yaml (preferred). If not found, will fall back to the FAIRe NOAA checklist.")

        dwc_data = {}
        checklist_df = pd.DataFrame()

        # 1) Preferred: YAML-based mappings
        yaml_mapper_path = 'data_mapper.yaml'
        try:
            if os.path.exists(yaml_mapper_path):
                import yaml as _yaml_loader
                with open(yaml_mapper_path, 'r', encoding='utf-8') as f:
                    mapper = _yaml_loader.safe_load(f) or {}

                # --- NEW: Select mapping based on metadata_format ---
                format_prefix = "generic_" if params.get('metadata_format') == 'GENERIC' else ""
                occurrence_key = f"{format_prefix}occurrence_core"
                dna_derived_key = f"{format_prefix}dna_derived_extension"
                
                reporter.add_text(f"Loading mappings for '{params.get('metadata_format')}' format from keys: {occurrence_key}, {dna_derived_key}")

                # Extract sections based on the selected format
                occ_map_raw = mapper.get(occurrence_key, {}) or {}
                dna_map_raw = mapper.get(dna_derived_key, {}) or {}

                def _build_df_from_mapping(mapping_dict):
                    rows = []
                    for dwc_term, mapping_details in mapping_dict.items():
                        # Skip if the entry is not a dictionary (e.g., a comment)
                        if not isinstance(mapping_details, dict):
                            continue
                        
                        faire_term = mapping_details.get('faire_term')
                        
                        # Skip blanks/nulls to allow "omit if unmapped"
                        if faire_term is None:
                            continue
                        
                        faire_str = str(faire_term).strip()
                        if not faire_str:
                            continue
                        rows.append({'DwC_term': str(dwc_term).strip(), 'FAIRe_term': faire_str})
                    if rows:
                        return pd.DataFrame(rows).drop_duplicates().set_index('DwC_term')
                    else:
                        return pd.DataFrame(columns=['FAIRe_term']).set_index(pd.Index([], name='DwC_term'))

                dwc_data['occurrence'] = _build_df_from_mapping(occ_map_raw)
                dwc_data['dnaDerived'] = _build_df_from_mapping(dna_map_raw)

                reporter.add_text(f"Loaded Darwin Core mappings from {yaml_mapper_path}")
                reporter.add_text(f"Occurrence mappings: {len(dwc_data['occurrence'])}; DNA Derived mappings: {len(dwc_data['dnaDerived'])}")

                reporter.add_dataframe(
                    dwc_data['occurrence'].reset_index(),
                    "Darwin Core Occurrence Mappings",
                    max_rows=len(dwc_data['occurrence'])
                )
                reporter.add_dataframe(
                    dwc_data['dnaDerived'].reset_index(),
                    "DNA Derived Data Mappings",
                    max_rows=len(dwc_data['dnaDerived'])
                )

                reporter.add_success("Loaded Darwin Core mappings")
                return dwc_data, checklist_df
        except Exception as e_yaml:
            reporter.add_warning(f"Could not load YAML mappings ({yaml_mapper_path}): {e_yaml}. Will attempt NOAA checklist.")

        # 2) Fallback: FAIRe NOAA checklist
     #   reporter.add_text("Falling back to FAIRe NOAA checklist for mappings...")

        try:
            checklist_df = pd.read_excel(
                params['FAIRe_NOAA_checklist'],
                sheet_name='checklist',
                na_values=[""]
            )
        except Exception as e:
            reporter.add_text(f"Error loading 'checklist' sheet: {e}")

        # Define helper columns for edna2obis from your checklist
        col_faire_term = 'term_name'
        col_output_spec = 'edna2obis_output_file'
        col_dwc_mapping = 'dwc_term'

        occurrence_maps = []
        dna_derived_maps = []

        # Process the checklist if it loaded successfully and has the required columns
        if not checklist_df.empty and all(col in checklist_df.columns for col in [col_faire_term, col_output_spec, col_dwc_mapping]):
            for _, row in checklist_df.iterrows():
                faire_term = row[col_faire_term]
                output_file = str(row[col_output_spec]).lower() # Convert to string and lowercase
                dwc_term = row[col_dwc_mapping]

                # Add to lists if terms are valid
                if pd.notna(faire_term) and pd.notna(dwc_term) and str(faire_term).strip() and str(dwc_term).strip():
                    if 'occurrence' in output_file:
                        occurrence_maps.append({'DwC_term': dwc_term, 'FAIRe_term': faire_term})
                    if 'dnaderived' in output_file:
                        dna_derived_maps.append({'DwC_term': dwc_term, 'FAIRe_term': faire_term})
            
            # Create DataFrames, using DwC_term as index
            dwc_data['occurrence'] = pd.DataFrame(occurrence_maps).drop_duplicates().set_index('DwC_term') if occurrence_maps else \
                                     pd.DataFrame(columns=['FAIRe_term']).set_index(pd.Index([], name='DwC_term'))
            
            dwc_data['dnaDerived'] = pd.DataFrame(dna_derived_maps).drop_duplicates().set_index('DwC_term') if dna_derived_maps else \
                                     pd.DataFrame(columns=['FAIRe_term']).set_index(pd.Index([], name='DwC_term'))
        else:
            # If checklist is empty or missing columns, create empty structures for dwc_data
            if checklist_df.empty and 'FAIRe_NOAA_checklist' in params: # Avoid double error message if file load failed
                reporter.add_text(f"Checklist DataFrame is empty or required columns are missing. Creating empty DwC mappings.")
            dwc_data['occurrence'] = pd.DataFrame(columns=['FAIRe_term']).set_index(pd.Index([], name='DwC_term'))
            dwc_data['dnaDerived'] = pd.DataFrame(columns=['FAIRe_term']).set_index(pd.Index([], name='DwC_term'))

        # Print summary
        reporter.add_text(f"Darwin Core term mapping created: \n   Occurrence Core mappings: {len(dwc_data['occurrence'])} \n   DNA Derived Extension mappings: {len(dwc_data['dnaDerived'])}")
        
        # Show the mappings (display full tables, not just a preview)
        reporter.add_dataframe(
            dwc_data['occurrence'].reset_index(),
            "Darwin Core Occurrence Mappings",
            max_rows=len(dwc_data['occurrence'])
        )
        reporter.add_dataframe(
            dwc_data['dnaDerived'].reset_index(),
            "DNA Derived Data Mappings",
            max_rows=len(dwc_data['dnaDerived'])
        )
        
        reporter.add_success("Loaded Darwin Core mappings")
        return dwc_data, checklist_df
        
    except Exception as e:
        error_msg = f"Failed to load Darwin Core mappings: {str(e)}"
        reporter.add_error(error_msg)
        raise Exception(error_msg)


def main():
    """Main execution function"""
    # --- Load config first to get params for reporter initialization ---
    print_header()
    print_usage()
    print_separator("Starting Conversion")

    # Always load the default config file first to check for a custom path
    default_config_path = "config.yaml"
    final_config_path = default_config_path
    
    try:
        with open(default_config_path, 'r', encoding='utf-8') as f:
            initial_config = yaml.safe_load(f)

        # Check if a different config file is specified within the default config
        if 'run_config_path' in initial_config and initial_config['run_config_path'] and initial_config['run_config_path'] != default_config_path:
            final_config_path = initial_config['run_config_path']
            console.print(f"Using config file: {final_config_path}")

        # Load the final configuration
        params = load_config(final_config_path)

    except Exception as e:
        console.print(f"[bold red]CRITICAL ERROR[/]: Could not load configuration file. {e}")
        # We can't create a report if config fails, so exit.
        return

    # --- Initialize HTML reporter with a dynamic name ---
    global reporter
    output_dir = params.get('output_dir', "processed-v3/")
    api_choice = params.get('taxonomic_api_source', 'unknown')
    run_name = params.get('edna2obis_run_name', 'unnamed_run')
    report_filename = f"edna2obis_report_{api_choice}_{run_name}.html"
    os.makedirs(output_dir, exist_ok=True)
    report_path = os.path.join(output_dir, report_filename)
    reporter = HTMLReporter(report_path, run_name)
    
    try:
        reporter.add_section("Data Cleaning", level=1)
        reporter.add_text_with_submission_logos("Starting initial data cleaning of eDNA metadata and raw data from FAIRe NOAA format to Darwin Core for OBIS and GBIF submission")
        
        # --- Report on the loaded configuration ---
        console.print("Loading and cleaning data...")
        reporter.add_section("Loading Configuration", level=2)
        reporter.add_success("Configuration loaded successfully")
        
        # Save the config file used for this run
        save_config_for_run(params, reporter, final_config_path)
        
        config_summary = [
            f"Run Name: {params.get('edna2obis_run_name', 'unnamed_run')}",
            f"API Source: {params['taxonomic_api_source']}",
        ]
        
        # Add file format information
        if params.get('use_excel', True):
            config_summary.append(f"Excel file: {params['excel_file']}")
        else:
            config_summary.append(f"Metadata format: TSV files")
            config_summary.append(f"  - projectMetadata: {params.get('projectMetadata_file', 'N/A')}")
            config_summary.append(f"  - sampleMetadata: {params.get('sampleMetadata_file', 'N/A')}")
            config_summary.append(f"  - experimentRunMetadata: {params.get('experimentRunMetadata_file', 'N/A')}")
        
        config_summary.extend([
            f"Number of analysis runs: {len(params['datafiles'])}",
            f"Output directory: {params['output_dir']}",
            f"Local reference database: {params.get('use_local_reference_database', False)}"
        ])
        
        reporter.add_list(config_summary, "Configuration Summary:")
        
        # Set up pandas display options
        setup_pandas_display()
        
        # Load project data
        console.print("Loading project data and metadata...")
        data = load_project_data(params, reporter)
        
        # Load ASV data
        console.print("Loading ASV data...")
        raw_data_tables = load_asv_data(params, reporter)
        
        # Remove control samples
        console.print("Removing control samples...")
        data, raw_data_tables = remove_control_samples(data, raw_data_tables, params, reporter)
        
        # Drop columns with all NAs
        console.print("Dropping columns with all NAs...")
        data = drop_all_na_columns(data, reporter)
        
        # Drop NA rows of each analysisMetadata sheet
        console.print("Dropping empty analysis metadata rows...")
        data = drop_empty_analysis_rows(data, reporter)
        
        # Optional: drop columns with some missing values (quiet)
        data = drop_some_na_columns(data, reporter)
        
        # Load Darwin Core mappings
        console.print("Loading Darwin Core mappings...")
        dwc_data, checklist_df = load_darwin_core_mappings(params, reporter)
        
        # Create occurrence core
        console.print("Creating Occurrence Core...")
        occurrence_core, all_processed_occurrence_dfs = create_occurrence_core(data, raw_data_tables, params, dwc_data, reporter)
        
        # Perform taxonomic assignment
        console.print("[bold]Starting Taxonomic Assignment...[/]")
        with console.status("Running Taxonomic Assignment...", spinner="dots"):
            # Suppress logs/errors but leave stdout for spinner
            with silence_output():
                assign_taxonomy(params, data, raw_data_tables, reporter)
        console.print("[green]Finished Taxonomic Assignment.[/]")
        
        # Create taxa assignment info file
        from taxonomic_assignment.taxa_assignment_manager import create_taxa_assignment_info
        with silence_output():
            create_taxa_assignment_info(params, reporter)
        
        # After creating the GBIF info file, remove any duplicate rows
        if params.get('taxonomic_api_source') == 'GBIF':
            from taxonomic_assignment.remove_GBIF_duplicates import remove_duplicates_from_gbif_taxa_info
            with silence_output():
                remove_duplicates_from_gbif_taxa_info(params, reporter)

            from taxonomic_assignment.mark_selected_gbif_match import mark_selected_gbif_matches
            with silence_output():
                mark_selected_gbif_matches(params, reporter)

        elif params.get('taxonomic_api_source') == 'WoRMS':
            from taxonomic_assignment.mark_selected_worms_match import mark_selected_worms_matches
            # Suppress internal prints while keeping spinner visible (spinner already ended here)
            with silence_output():
                mark_selected_worms_matches(params, reporter)
            
        # Remove match_type_debug from final occurrence file (keep it only in taxa_assignment_INFO.csv)
        api_source = params.get('taxonomic_api_source', 'WoRMS').lower()
        final_occurrence_path = os.path.join(params.get('output_dir', 'processed-v3/'), f'occurrence_core_{api_source}.csv')
        try:
            if os.path.exists(final_occurrence_path):
                # Read the file, remove match_type_debug column, and save back
                final_df = pd.read_csv(final_occurrence_path)
                if 'match_type_debug' in final_df.columns:
                    final_df = final_df.drop(columns=['match_type_debug'])
                    final_df.to_csv(final_occurrence_path, index=False, na_rep='')
                    reporter.add_text(f"🔧 Removed match_type_debug from final occurrence file (kept in taxa_assignment_INFO.csv)")
        except Exception as e:
            reporter.add_text(f"⚠️ Warning: Could not remove match_type_debug from final file: {e}")
        
        # Delete intermediate occurrence.csv file
        intermediate_occurrence_path = os.path.join(params.get('output_dir', 'processed-v3/'), 'occurrence.csv')
        try:
            if os.path.exists(intermediate_occurrence_path):
                os.remove(intermediate_occurrence_path)
                reporter.add_text(f"🗑️ Deleted intermediate file: {intermediate_occurrence_path}")
        except Exception as e:
            reporter.add_text(f"⚠️ Warning: Could not delete intermediate occurrence.csv: {e}")
        
        # Create DNA derived extension
        console.print("Creating DNA derived extension...")
        with silence_output():
            create_dna_derived_extension(params, data, raw_data_tables, dwc_data, occurrence_core, all_processed_occurrence_dfs, checklist_df, reporter)
        
        # Create eMoF file (optional)
        if params.get('emof_enabled', True):
            console.print("Creating eMoF (extendedMeasurementOrFact)...")
            try:
                with silence_output():
                    emof_path = create_emof_table(params, occurrence_core, data, reporter)
                reporter.add_text(f"eMoF saved to: {emof_path}")
            except Exception as e:
                reporter.add_warning(f"eMoF creation failed: {e}")
        else:
            reporter.add_text("⏭️ Skipping eMoF creation per config (emof_enabled=false)")
        
        # Create EML file (optional)
        if params.get('eml_enabled', False):
            console.print("Creating EML (Ecological Metadata Language) file...")
            try:
                from create_EML.EML_builder import create_eml_file
                with silence_output():
                    eml_path = create_eml_file(params, data, reporter)
                reporter.add_text(f"EML saved to: {eml_path}")
            except Exception as e:
                reporter.add_warning(f"EML creation failed: {e}")
        else:
            reporter.add_text("⏭️ Skipping EML creation per config (eml_enabled=false)")
        
        # --- Final File Validation ---
        reporter.add_section("Final File Validation", level=3)
        output_dir = params.get('output_dir', 'processed-v3/')
        api_choice = params.get("taxonomic_api_source", "worms")
        
        files_to_validate = [
            f'occurrence_core_{api_choice.lower()}.csv',
            f'taxa_assignment_INFO_{api_choice}.csv',
            'dna_derived_extension.csv'
        ]
        if params.get('emof_enabled', True):
            files_to_validate.append('eMoF.csv')
        if params.get('eml_enabled', False):
            files_to_validate.append('eml.xml')

        all_empty_columns_summary = []
        removed_columns_summary = []
        for filename in files_to_validate:
            filepath = os.path.join(output_dir, filename)
            if os.path.exists(filepath):
                try:
                    if filename.lower().endswith('.xml'):
                        # Validate XML structure instead of reading as CSV
                        ET.parse(filepath)
                    else:
                        sep = '\t' if filename.lower().endswith('.tsv') else ','
                        df = pd.read_csv(filepath, sep=sep, low_memory=False)
                        
                        # Existing empty column check (for reporting)
                        empty_columns = [col for col in df.columns if df[col].isna().all()]
                        if empty_columns:
                            for col in empty_columns:
                                reporter.add_warning(f"In output file <strong>'{filename}'</strong>, the column <strong>'{col}'</strong> was found to be completely empty.")
                                all_empty_columns_summary.append(f"File: <code>{filename}</code>, Column: <code>{col}</code>")
                        
                        # --- NEW: Remove empty columns and overwrite file ---
                        if empty_columns:
                            df.drop(columns=empty_columns, inplace=True)
                            df.to_csv(filepath, index=False, sep=sep, encoding='utf-8-sig')
                            for col in empty_columns:
                                removed_columns_summary.append(f"File: <code>{filename}</code>, Removed Column: <code>{col}</code>")

                except Exception as e:
                    reporter.add_warning(f"Could not validate or clean file '{filename}': {e}")
        
        if not all_empty_columns_summary:
            reporter.add_success("Validation complete: No empty columns found in final output files.")
        else:
            # This report section now serves as a pre-cleanup warning
            reporter.add_list(all_empty_columns_summary, "<h4>Summary of All Empty Columns Found (and subsequently removed):</h4>")

        # --- Report on removed columns ---
        if removed_columns_summary:
            reporter.add_section("Empty Column Cleanup", level=3)
            reporter.add_success("Removed completely empty columns from the final output files.")
            reporter.add_list(removed_columns_summary, "<h4>Summary of Removed Columns:</h4>")

        # --- Split output files by short_name (cruise/expedition) if enabled ---
        if params.get('split_output_by_short_name', False):
            console.print("[bold]Splitting output files by short_name (cruise/expedition)...[/]")
            try:
                split_output_files_by_short_name(params, data, reporter)
                console.print("[green]Output files split successfully![/]")
            except Exception as e:
                console.print(f"[bold red]Failed to split output files: {e}[/]")
                reporter.add_error(f"Failed to split output files by short_name: {e}")
                # Re-raise to stop the pipeline since user explicitly requested splitting
                raise
        else:
            console.print("Skipping output file splitting (split_output_by_short_name=false)")
            reporter.add_text("⏭️ Skipping output file splitting (split_output_by_short_name=false)")

        # --- Final Status Check ---
        # If any warnings were logged during the run, set the final status to WARNING
        if reporter.warnings:
            reporter.set_warning()
        else:
            reporter.set_success()
            
        # Generate final report
        console.print("\n[bold green]Process completed.[/]")
        console.print("Generating HTML report...")
        
        reporter.add_section("Process Completion")
        reporter.add_success("All steps completed successfully!")
        reporter.add_text("Files generated:")
        
        output_dir = params.get('output_dir', '../processed-v3/')
        api_choice = params.get("taxonomic_api_source", "worms")
        
        # Define the list of expected final files
        files = [
            f'occurrence_core_{api_choice.lower()}.csv', 
            f'taxa_assignment_INFO_{api_choice}.csv', 
            'dna_derived_extension.csv'
        ]
        
        # Add optional files based on configuration
        if params.get('emof_enabled', True):
            files.append('eMoF.csv')
        if params.get('eml_enabled', False):
            files.append('eml.xml')
        
        # Add the saved config file
        run_name = params.get('edna2obis_run_name', 'unnamed_run')
        files.append(f'config_{run_name}.yaml')
        
        files.append(report_filename)  # Use the dynamic report filename
        
        for filename in files:
            filepath = os.path.join(output_dir, filename)
            if os.path.exists(filepath):
                size_mb = os.path.getsize(filepath) / (1024*1024)
                reporter.add_text(f"✓ {filename} ({size_mb:.2f} MB)")
                # Clean CLI success messages for final outputs
                if filename.startswith("occurrence_core_"):
                    console.print(f"[green]Created Occurrence Core:[/] {filepath}")
                elif filename == "dna_derived_extension.csv":
                    console.print(f"[green]Created DNA Derived Extension:[/] {filepath}")
                elif filename == "eMoF.csv":
                    console.print(f"[green]Created eMoF:[/] {filepath}")
                elif filename.endswith(".html"):
                    console.print(f"[green]HTML report:[/] {filepath}")
                elif filename.endswith(".xml"):
                    console.print(f"[green]EML generated:[/] {filepath}")
            else:
                reporter.add_text(f"✗ {filename} (not found)")
        
    except Exception as e:
        console.print(f"\n[bold red]Error during processing:[/] {str(e)}")
        error_msg = f"Pipeline failed: {str(e)}"
        reporter.add_error(error_msg)
        reporter.add_text("Full traceback:")
        import traceback
        reporter.add_text(traceback.format_exc())
        reporter.set_status("FAILED", error_msg)
        raise
    
    finally:
        # Save the report
        reporter.save()
        console.print("\nedna2obis conversion complete.")
        reporter.open_in_browser()


if __name__ == "__main__":
    main() 