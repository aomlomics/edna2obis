import pandas as pd
import os

def mark_selected_worms_matches(params, reporter=None):
    """
    Reads taxa_assignment_INFO_WoRMS.csv and the final occurrence_WoRMS_matched.csv
    to mark which of the ambiguous WoRMS matches was selected.
    """
    try:
        api_choice = params.get('taxonomic_api_source', 'GBIF')
        if api_choice != 'WoRMS':
            print("Skipping WoRMS match marking as the API source is not WoRMS.")
            return

        output_dir = params.get('output_dir', '../processed-v3/')
        info_filename = f"taxa_assignment_INFO_{api_choice}.csv"
        info_filepath = os.path.join(output_dir, info_filename)

        occurrence_filename = f"occurrence_{api_choice.lower()}_matched.csv"
        occurrence_filepath = os.path.join(output_dir, occurrence_filename)

        if not os.path.exists(info_filepath):
            print(f"Error: Taxa info file not found at {info_filepath}")
            return
        if not os.path.exists(occurrence_filepath):
            print(f"Error: Occurrence file not found at {occurrence_filepath}")
            return

        print("Reading files to mark selected WoRMS matches...")
        info_df = pd.read_csv(info_filepath)
        
        if info_df.empty:
            print(f"Skipping marking for empty file: {info_filename}")
            return

        occurrence_df = pd.read_csv(occurrence_filepath, usecols=['verbatimIdentification', 'scientificNameID'])

        # Before comparison, ensure data types are compatible by converting to string.
        info_df['verbatimIdentification'] = info_df['verbatimIdentification'].astype(str)
        info_df['scientificNameID'] = info_df['scientificNameID'].astype(str)
        occurrence_df['verbatimIdentification'] = occurrence_df['verbatimIdentification'].astype(str)
        occurrence_df['scientificNameID'] = occurrence_df['scientificNameID'].astype(str)

        # Mark all as False initially
        info_df['selected_match'] = False

        # Identify the selected verbatim-taxon pairs from the occurrence file.
        selected_pairs = occurrence_df[['verbatimIdentification', 'scientificNameID']].drop_duplicates()
        print(f"Found {len(selected_pairs)} unique selected pairs in the WoRMS occurrence file.")

        # Use a merge operation to find all candidate matches in the info file.
        merged_df = pd.merge(
            info_df,
            selected_pairs,
            on=['verbatimIdentification', 'scientificNameID'],
            how='inner'
        )

        if not merged_df.empty:
            # Drop duplicates to get the unique indices of matching rows in the original info_df
            best_candidates = merged_df.drop_duplicates(subset=['verbatimIdentification', 'scientificNameID'], keep='first')
            
            # Mark the 'selected_match' column as True for the indices of the candidates.
            info_df.loc[best_candidates.index, 'selected_match'] = True

        # Move 'selected_match' column to be after 'ambiguous' for consistency
        if 'ambiguous' in info_df.columns:
            cols = info_df.columns.tolist()
            if 'selected_match' in cols:
                cols.remove('selected_match')
            
            # Find the index of 'ambiguous' and insert 'selected_match' right after it
            ambiguous_idx = cols.index('ambiguous')
            cols.insert(ambiguous_idx + 1, 'selected_match')
            info_df = info_df[cols]

        # Save the updated dataframe back to the info file
        info_df.to_csv(info_filepath, index=False, na_rep='')
        
        num_selected = info_df['selected_match'].sum()
        print(f"Successfully marked {num_selected} selected matches in {info_filename}.")

        # If a reporter object is provided, force it to update its internal view
        if reporter:
            try:
                reporter.update_dataframe_from_file(info_filename, info_filepath)
                print(f"Successfully notified reporter to update its view of '{info_filename}'.")
            except Exception as e:
                print(f"Warning: Could not update reporter view for {info_filename}. Report may be stale. Error: {e}")

    except Exception as e:
        print(f"An error occurred during WoRMS match marking: {str(e)}")
        import traceback
        traceback.print_exc()

if __name__ == '__main__':
    # This block allows the script to be run standalone for testing
    mock_params = {
        'taxonomic_api_source': 'WoRMS',
        'output_dir': '../../processed-v3/'
    }
    # To test, you would need sample files in the 'processed-v3' directory.
    # mark_selected_worms_matches(mock_params)
    print("This script is intended to be called from the main edna2obis pipeline.") 