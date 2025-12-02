import pandas as pd
import os

def remove_duplicates_from_gbif_taxa_info(params, reporter=None):
    """
    Reads the taxa_assignment_INFO_GBIF.csv file, removes duplicate rows,
    and overwrites the file with the deduplicated data.
    """
    try:
        api_choice = params.get('taxonomic_api_source', 'WoRMS')
        if api_choice != 'GBIF':
            return

        output_dir = params.get('output_dir', '../processed-v3/')
        filename = f"taxa_assignment_INFO_{api_choice}.csv"
        filepath = os.path.join(output_dir, filename)

        if not os.path.exists(filepath):
            return

        df = pd.read_csv(filepath)
        
        initial_rows = len(df)

        # Remove duplicate rows based on all columns
        df.drop_duplicates(inplace=True)
        
        final_rows = len(df)

        # Save the deduplicated dataframe back to the original file
        df.to_csv(filepath, index=False, na_rep='')

        if reporter:
            removed = initial_rows - final_rows
            reporter.add_text(f"Deduplicated {filename}: removed {removed} duplicate rows.")

    except Exception as e:
        if reporter:
            reporter.add_warning(f"GBIF deduplication error: {str(e)}")

if __name__ == '__main__':
    # This block allows the script to be run standalone for testing
    # You would need to provide a sample params dictionary.
    mock_params = {
        'taxonomic_api_source': 'GBIF',
        'output_dir': '../../processed-v3/'
    }
    # To test, you would need a sample 'taxa_assignment_INFO_GBIF.csv' in the 'processed-v3' directory.
    # remove_duplicates_from_gbif_taxa_info(mock_params)
    print("This script is intended to be called from the main edna2obis pipeline.")
