def get_analysis_data_for_marker(data, marker, analysis_run_name=None):
    """
    Get analysis data for a specific marker
    
    Parameters:
    -----------
    data : dict
        The main data dictionary
    marker : str
        The marker name (e.g., '16S V4-V5')
    analysis_run_name : str, optional
        The specific analysis run to retrieve. If None, returns the first analysis.
        
    Returns:
    --------
    pandas.DataFrame
        The analysis metadata dataframe
    """
    assay = data['marker_to_assay'].get(marker)
    if not assay or assay not in data['analysis_data_by_assay']:
        print(f"Warning: No analysis data found for marker {marker}")
        return None
    
    analyses = data['analysis_data_by_assay'][assay]
    
    # If a specific analysis run is requested
    if analysis_run_name and analysis_run_name in analyses:
        return analyses[analysis_run_name]
    
    # Otherwise return the first analysis
    return next(iter(analyses.values()))

def get_project_info_for_assay(data, assay, field_name):
    """Get project-level information for a specific assay"""
    if not assay:
        return None
    
    # Check if there's an assay-specific column
    assay_column = f"{field_name}_{assay}"
    if assay_column in data['projectMetadata'].columns:
        return data['projectMetadata'][assay_column].iloc[0]
    
    # Otherwise check for pipe-separated values
    if field_name in data['projectMetadata'].columns:
        value = data['projectMetadata'][field_name].iloc[0]
        if '|' in str(value):
            values = str(value).split('|')
            assays = str(data['projectMetadata']['assay_name'].iloc[0]).split('|')
            assays = [a.strip() for a in assays]
            if assay in assays:
                idx = assays.index(assay)
                if idx < len(values):
                    return values[idx].strip()
    
    return None

def get_project_info_for_marker(data, marker, field_name):
    """Get project-level information for a specific marker"""
    assay = data['marker_to_assay'].get(marker)
    return get_project_info_for_assay(data, assay, field_name)