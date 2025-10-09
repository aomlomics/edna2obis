import pandas as pd
import os

def mark_selected_gbif_matches(params, reporter=None):
    """
    Reads taxa_assignment_INFO_GBIF.csv and the final occurrence_gbif_matched.csv
    to mark which of the GBIF matches was selected.
    """
    try:
        api_choice = params.get('taxonomic_api_source', 'GBIF')
        if api_choice != 'GBIF':
            print("Skipping GBIF match marking as the API source is not GBIF.")
            return

        output_dir = params.get('output_dir', '../processed-v3/')
        info_filename = f"taxa_assignment_INFO_{api_choice}.csv"
        info_filepath = os.path.join(output_dir, info_filename)

        occurrence_filename = f"occurrence_core_{api_choice.lower()}.csv"
        occurrence_filepath = os.path.join(output_dir, occurrence_filename)

        if not os.path.exists(info_filepath):
            print(f"Error: Taxa info file not found at {info_filepath}")
            return
        if not os.path.exists(occurrence_filepath):
            print(f"Error: Occurrence file not found at {occurrence_filepath}")
            return

        print("Reading files to mark selected GBIF matches...")
        info_df = pd.read_csv(info_filepath)
        
        if info_df.empty:
            print(f"Skipping marking for empty file: {info_filename}")
            return
            
        occurrence_df = pd.read_csv(occurrence_filepath, usecols=['verbatimIdentification', 'taxonID'])

        # Before comparison, ensure data types are compatible without making assumptions
        # about the content. Convert to string to be safe.
        info_df['verbatimIdentification'] = info_df['verbatimIdentification'].astype(str)
        info_df['taxonID'] = info_df['taxonID'].astype(str)
        occurrence_df['verbatimIdentification'] = occurrence_df['verbatimIdentification'].astype(str)
        occurrence_df['taxonID'] = occurrence_df['taxonID'].astype(str)

        # Mark all as False initially
        info_df['selected_match'] = False

        # Identify the selected verbatim-taxon pairs from the occurrence file.
        # We only need the unique pairs that were chosen.
        selected_pairs = occurrence_df[['verbatimIdentification', 'taxonID']].drop_duplicates()
        print(f"Found {len(selected_pairs)} unique selected pairs in the occurrence file.")

        # Use a merge operation to find all candidate matches in the info file.
        # This is more robust than iterating or using sets.
        merged_df = pd.merge(
            info_df,
            selected_pairs,
            on=['verbatimIdentification', 'taxonID'],
            how='inner'
        )

        if not merged_df.empty:
            # From the candidates, find the single best one for each pair (highest confidence).
            # This handles any cases where the same pair might appear multiple times.
            best_candidates = merged_df.sort_values('confidence', ascending=False).drop_duplicates(
                subset=['verbatimIdentification', 'taxonID'],
                keep='first'
            )

            # Mark the 'selected_match' column as True only for the indices of the best candidates.
            info_df.loc[best_candidates.index, 'selected_match'] = True

        # Save the updated dataframe back to the info file
        info_df.to_csv(info_filepath, index=False, na_rep='')
        
        num_selected = info_df['selected_match'].sum()
        print(f"Successfully marked {num_selected} selected matches in {info_filename}.")

        # If a reporter object is provided, force it to update its internal
        # dataframe from the file we just saved. This ensures the report
        # reflects the newly added 'selected_match' column.
        if reporter:
            try:
                reporter.update_dataframe_from_file(info_filename, info_filepath)
                print(f"Successfully notified reporter to update its view of '{info_filename}'.")
            except Exception as e:
                print(f"Warning: Could not update reporter view for {info_filename}. Report may be stale. Error: {e}")

    except Exception as e:
        print(f"An error occurred during GBIF match marking: {str(e)}")
        import traceback
        traceback.print_exc()

if __name__ == '__main__':
    # This block allows the script to be run standalone for testing
    mock_params = {
        'taxonomic_api_source': 'GBIF',
        'output_dir': '../../processed-v3/'
    }
    # To test, you would need sample files in the 'processed-v3' directory.
    # mark_selected_gbif_matches(mock_params)
    print("This script is intended to be called from the main edna2obis pipeline.") 