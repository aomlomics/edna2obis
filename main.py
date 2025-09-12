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

# Add the src-v3 directory to Python path for imports
sys.path.insert(0, "src-v3")

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
from create_EML.EML_builder import create_eml_file


def load_config(config_path="config.yaml"):
    """Load configuration from YAML file and convert to params dict structure"""
    try:
        if not os.path.exists(config_path):
            raise FileNotFoundError(f"Config file not found: {config_path}")
        
        with open(config_path, 'r', encoding='utf-8') as f:
            config = yaml.safe_load(f)
        
        params = {}
        
        params['sampleMetadata'] = config['sampleMetadata']
        params['experimentRunMetadata'] = config['experimentRunMetadata'] 
        params['projectMetadata'] = config['projectMetadata']
        params['excel_file'] = config['excel_file']
        params['FAIRe_NOAA_checklist'] = config['FAIRe_NOAA_checklist']
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
        params['assays_to_skip_species_match'] = config['assays_to_skip_species_match']
        
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
        
        return params
        
    except Exception as e:
        raise Exception(f"Error loading config: {str(e)}")


def setup_pandas_display():
    """Set pandas display options"""
    pd.set_option('display.max_colwidth', 150)
    pd.set_option('display.max_columns', 50)


def load_project_data(params, reporter):
    """Load project, sample, experimentRun, and analysis data from FAIRe Excel file"""
    reporter.add_section("Loading Project Data", level=2)
    
    try:
        # Discover all sheets in the Excel file
        excel = pd.ExcelFile(params['excel_file'])
        all_sheets = excel.sheet_names
        reporter.add_text(f"Found {len(all_sheets)} sheets in Excel file: {', '.join(all_sheets)}")
        
        # Find analysis metadata sheets
        analysis_sheets = [sheet for sheet in all_sheets if sheet.startswith('analysisMetadata')]
        reporter.add_text(f"Found {len(analysis_sheets)} analysis metadata sheets: {', '.join(analysis_sheets)}")
        
        # Load the main data sheets
        data = pd.read_excel(
            params['excel_file'],
            [params['projectMetadata'], params['sampleMetadata'], params['experimentRunMetadata']],
            index_col=None, na_values=[""], comment="#"
        )
        
        # Load all analysis metadata sheets
        analysis_data_by_assay = {}
        
        for sheet_name in analysis_sheets:
            analysis_df = pd.read_excel(params['excel_file'], sheet_name)
            
            # Get assay_name and analysis_run_name from specific cells
            assay_name = str(analysis_df.iloc[1, 3])  # Excel cell D3
            analysis_run_name = str(analysis_df.iloc[2, 3])  # Excel cell D4
            
            reporter.add_text(f"Processing sheet '{sheet_name}': assay '{assay_name}', run '{analysis_run_name}'")
            
            if assay_name not in analysis_data_by_assay:
                analysis_data_by_assay[assay_name] = {}
            analysis_data_by_assay[assay_name][analysis_run_name] = analysis_df
        
        # Add analysis data to main data dictionary
        data['analysis_data_by_assay'] = analysis_data_by_assay
        
        # For backward compatibility
        if analysis_sheets:
            data['analysisMetadata'] = pd.read_excel(params['excel_file'], analysis_sheets[0])
        
        # Rename keys to standard terms
        data['sampleMetadata'] = data.pop(params['sampleMetadata'])
        data['experimentRunMetadata'] = data.pop(params['experimentRunMetadata'])
        data['projectMetadata'] = data.pop(params['projectMetadata'])
        
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
                    raw_data_tables[analysis_run_name]['taxonomy'] = pd.read_table(tax_path, sep='\t', low_memory=False)
                    tax_shape = raw_data_tables[analysis_run_name]['taxonomy'].shape
                    reporter.add_success(f"Loaded taxonomy file: {tax_path} (shape: {tax_shape})")
                except Exception as e:
                    reporter.add_error(f"Failed to load taxonomy file {tax_path}: {e}")
                    raise
            
            # Load abundance table file  
            if 'occurrence_file' in file_paths:
                abundance_path = file_paths['occurrence_file']
                try:
                    # --- Smartly handle old and new Tourmaline formats ---
                    
                    # 1. Check the first line to see if we need to skip it.
                    with open(abundance_path, 'r', encoding='utf-8', errors='ignore') as f:
                        first_line = f.readline()
                    
                    # The old format has a comment line, the new one starts with the header.
                    rows_to_skip = 1 if '# Constructed from biom file' in first_line else 0

                    # 2. Load the table, skipping the comment line only if it exists.
                    df_abundance = pd.read_table(abundance_path,
                                                 sep='\t',
                                                 skiprows=rows_to_skip,
                                                 header=0,
                                                 low_memory=False)
                    
                    # 3. Standardize the first column name to 'featureid'.
                    first_col_name = df_abundance.columns[0]
                    
                    # First, remove the leading '#' if it exists.
                    if first_col_name.startswith('#'):
                        clean_col_name = first_col_name[1:]
                        df_abundance.rename(columns={first_col_name: clean_col_name}, inplace=True)
                        first_col_name = clean_col_name # Update for the next check
                    
                    # Now, if the column is 'OTU ID', rename it to 'featureid'.
                    if first_col_name.lower() == 'otu id':
                        df_abundance.rename(columns={first_col_name: 'featureid'}, inplace=True)

                    raw_data_tables[analysis_run_name]['occurrence'] = df_abundance
                    abundance_shape = raw_data_tables[analysis_run_name]['occurrence'].shape
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
        
        # Drop columns with some missing values
        sheets_to_clean = ['sampleMetadata', 'experimentRunMetadata']

        for sheet_name in sheets_to_clean:
            if sheet_name in data and not data[sheet_name].empty: # Check if DataFrame exists and is not empty
                original_columns = data[sheet_name].columns.tolist() # Get column names before dropping
                
                data[sheet_name].dropna(axis=1, how='any', inplace=True)
                
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
        reporter.add_text("This FAIRe NOAA Checklist Excel file also contains columns for mapping FAIRe fields to the appropriate Darwin Core terms which OBIS is expecting. Currently, we are only preparing an Occurrence core file and a DNA-derived extension file, with Event information in the Occurrence file. Future versions of this workflow will prepare an extendedMeasurementOrFact file as well.")
        
        dwc_data = {}
        checklist_df = pd.DataFrame()

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
        
        # Show the mappings
        reporter.add_dataframe(dwc_data['occurrence'].reset_index(), "Darwin Core Occurrence Mappings", max_rows=15)
        reporter.add_dataframe(dwc_data['dnaDerived'].reset_index(), "DNA Derived Data Mappings", max_rows=15)
        
        reporter.add_success("Loaded Darwin Core mappings")
        return dwc_data
        
    except Exception as e:
        error_msg = f"Failed to load Darwin Core mappings: {str(e)}"
        reporter.add_error(error_msg)
        raise Exception(error_msg)


def main():
    """Main execution function"""
    # --- Load config first to get params for reporter initialization ---
    try:
        params = load_config()
    except Exception as e:
        print(f"CRITICAL ERROR: Could not load configuration file. {e}")
        # We can't create a report if config fails, so exit.
        return

    # --- Initialize HTML reporter with a dynamic name ---
    global reporter
    output_dir = params.get('output_dir', "processed-v3/")
    api_choice = params.get('taxonomic_api_source', 'unknown')
    report_filename = f"edna2obis_report_{api_choice}.html"
    os.makedirs(output_dir, exist_ok=True)
    report_path = os.path.join(output_dir, report_filename)
    reporter = HTMLReporter(report_path)
    
    print("Starting edna2obis conversion...")
    print("="*50)
    
    try:
        reporter.add_section("Data Cleaning", level=1)
        reporter.add_text_with_submission_logos("Starting initial data cleaning of eDNA metadata and raw data from FAIRe NOAA format to Darwin Core for OBIS and GBIF submission")
        
        # --- Report on the loaded configuration ---
        print("🔄 Loading and Cleaning Data...")
        reporter.add_section("Loading Configuration", level=2)
        reporter.add_success("Configuration loaded successfully")
        reporter.add_list([
            f"API Source: {params['taxonomic_api_source']}",
            f"Excel file: {params['excel_file']}", 
            f"Number of analysis runs: {len(params['datafiles'])}",
            f"Output directory: {params['output_dir']}",
            f"Local reference database: {params.get('use_local_reference_database', False)}"
        ], "Configuration Summary:")
        
        # Set up pandas display options
        print("Setting up Pandas display options...")
        setup_pandas_display()
        
        # Load project data
        print("Loading project data and metadata...")
        data = load_project_data(params, reporter)
        
        # Load ASV data
        print("Loading ASV data...")
        raw_data_tables = load_asv_data(params, reporter)
        
        # Remove control samples
        print("Removing control samples...")
        data, raw_data_tables = remove_control_samples(data, raw_data_tables, params, reporter)
        
        # Drop columns with all NAs
        print("Dropping columns with all NAs...")
        data = drop_all_na_columns(data, reporter)
        
        # Drop NA rows of each analysisMetadata sheet
        print("Dropping empty analysis metadata rows...")
        data = drop_empty_analysis_rows(data, reporter)
        
        # Drop columns with some missing values
        print("Dropping columns with some missing values...")
        data = drop_some_na_columns(data, reporter)
        
        # Load Darwin Core mappings
        print("Loading Darwin Core mappings...")
        dwc_data = load_darwin_core_mappings(params, reporter)
        
        # Create occurrence core
        print("📄 Creating Occurrence Core...")
        occurrence_core, all_processed_occurrence_dfs = create_occurrence_core(data, raw_data_tables, params, dwc_data, reporter)
        
        # Perform taxonomic assignment
        print("🐟 Performing taxonomic assignment...")
        assign_taxonomy(params, data, raw_data_tables, reporter)
        
        # Create taxa assignment info file
        from taxonomic_assignment.taxa_assignment_manager import create_taxa_assignment_info
        create_taxa_assignment_info(params, reporter)
        
        # After creating the GBIF info file, remove any duplicate rows
        if params.get('taxonomic_api_source') == 'GBIF':
            from taxonomic_assignment.remove_GBIF_duplicates import remove_duplicates_from_gbif_taxa_info
            print("Removing duplicate rows from GBIF taxa info file...")
            remove_duplicates_from_gbif_taxa_info(params)

            from taxonomic_assignment.mark_selected_gbif_match import mark_selected_gbif_matches
            print("Marking selected matches in GBIF taxa info file...")
            mark_selected_gbif_matches(params, reporter)

        elif params.get('taxonomic_api_source') == 'WoRMS':
            from taxonomic_assignment.mark_selected_worms_match import mark_selected_worms_matches
            print("Marking selected matches in WoRMS taxa info file...")
            mark_selected_worms_matches(params, reporter)
            
        # Remove match_type_debug from final occurrence file (keep it only in taxa_assignment_INFO.csv)
        api_source = params.get('taxonomic_api_source', 'WoRMS').lower()
        final_occurrence_path = os.path.join(params.get('output_dir', 'processed-v3/'), f'occurrence_{api_source}_matched.csv')
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
        print("📄 Creating DNA derived extension...")
        create_dna_derived_extension(params, data, raw_data_tables, dwc_data, occurrence_core, all_processed_occurrence_dfs, reporter)
        
        # Create eMoF file (optional)
        if params.get('emof_enabled', True):
            print("📄 Creating eMoF (extendedMeasurementOrFact)...")
            try:
                emof_path = create_emof_table(params, occurrence_core, data, reporter)
                reporter.add_text(f"eMoF saved to: {emof_path}")
            except Exception as e:
                reporter.add_warning(f"eMoF creation failed: {e}")
        else:
            reporter.add_text("⏭️ Skipping eMoF creation per config (emof_enabled=false)")
        
        # Create EML file (optional)
        if params.get('eml_enabled', False):
            print("📄 Creating EML (Ecological Metadata Language) file...")
            try:
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
            f'occurrence_{api_choice.lower()}_matched.csv',
            f'taxa_assignment_INFO_{api_choice}.csv',
            'dna_derived_extension.csv'
        ]
        if params.get('emof_enabled', True):
            files_to_validate.append('eMoF.csv')
        if params.get('eml_enabled', False):
            files_to_validate.append('eml.xml')

        all_empty_columns_summary = []
        for filename in files_to_validate:
            filepath = os.path.join(output_dir, filename)
            if os.path.exists(filepath):
                try:
                    sep = '\t' if filename.lower().endswith('.tsv') else ','
                    df = pd.read_csv(filepath, sep=sep, low_memory=False)
                    empty_columns = [col for col in df.columns if df[col].isna().all()]
                    if empty_columns:
                        for col in empty_columns:
                            reporter.add_warning(f"In output file <strong>'{filename}'</strong>, the column <strong>'{col}'</strong> was found to be completely empty.")
                            all_empty_columns_summary.append(f"File: <code>{filename}</code>, Column: <code>{col}</code>")
                except Exception as e:
                    reporter.add_warning(f"Could not validate file '{filename}': {e}")
        
        if not all_empty_columns_summary:
            reporter.add_success("Validation complete: No empty columns found in final output files.")
        else:
            reporter.add_list(all_empty_columns_summary, "<h4>Summary of All Empty Columns Found:</h4>")

        # --- Final Status Check ---
        # If any warnings were logged during the run, set the final status to WARNING
        if reporter.warnings:
            reporter.set_warning()
        else:
            reporter.set_success()
            
        # Generate final report
        print("\n🎉 Process completed!")
        print("Generating HTML report...")
        
        reporter.add_section("Process Completion")
        reporter.add_success("All steps completed successfully!")
        reporter.add_text("Files generated:")
        
        output_dir = params.get('output_dir', '../processed-v3/')
        api_choice = params.get("taxonomic_api_source", "worms")
        
        # Define the list of expected final files
        files = [
            f'occurrence_{api_choice.lower()}_matched.csv', 
            f'taxa_assignment_INFO_{api_choice}.csv', 
            'dna_derived_extension.csv'
        ]
        
        # Add optional files based on configuration
        if params.get('emof_enabled', True):
            files.append('eMoF.csv')
        if params.get('eml_enabled', False):
            files.append('eml.xml')
        
        files.append(report_filename)  # Use the dynamic report filename
        
        for filename in files:
            filepath = os.path.join(output_dir, filename)
            if os.path.exists(filepath):
                size_mb = os.path.getsize(filepath) / (1024*1024)
                reporter.add_text(f"✓ {filename} ({size_mb:.2f} MB)")
            else:
                reporter.add_text(f"✗ {filename} (not found)")
        
    except Exception as e:
        print(f"\n❌ Error during processing: {str(e)}")
        error_msg = f"Pipeline failed: {str(e)}"
        reporter.add_error(error_msg)
        reporter.add_text("Full traceback:")
        import traceback
        reporter.add_text(traceback.format_exc())
        reporter.set_status("FAILED", error_msg)
        raise
    
    finally:
        # Save the report
        print(f"HTML report saved: {reporter.filename}")
        print("\n✓ edna2obis conversion complete!")
        reporter.save()
        reporter.open_in_browser()


if __name__ == "__main__":
    main() 