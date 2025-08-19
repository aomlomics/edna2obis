"""
Extended Measurement Or Fact (eMoF) Builder for edna2obis
Creates an eMoF Excel file based on a small template Excel and the occurrence core.

Behavior:
- Reads a template Excel (default: raw-v3/eMoF Fields edna2obis .xlsx)
  The template should have at least: measurementType, measurementValue,
  measurementTypeID, measurementValueID, measurementUnitID, measurementRemarks.
- Units are automatically detected from FAIRe metadata using the pattern: {measurementType}_unit
  For example, if measurementType = "temperature", the system looks for a "temperature_unit" column
  in sampleMetadata, experimentRunMetadata, or projectMetadata.
- Optionally, the template can include a column named 'linkTo' with values
  'occurrence' (default) or 'event'.
  - linkTo == 'occurrence': a row is expanded for every occurrence (eventID + occurrenceID)
  - linkTo == 'event': a row is expanded for every unique event (eventID only)
- Special handling: measurementType == 'assay_name'
  If the template specifies measurementValue for assay_name, rows will be created
  only for occurrences whose assay_name equals the provided measurementValue.
  If no measurementValue is provided in the template row, a per-occurrence row is
  created where measurementValue is populated from each occurrence's assay_name.

Output columns (in order):
  eventID, occurrenceID, measurementType, measurementValue, measurementUnit,
  measurementTypeID, measurementValueID, measurementUnitID, measurementRemarks

The function logs progress to the provided HTMLReporter.
"""

import os
import pandas as pd


REQUIRED_COLUMNS_IN_ORDER = [
    'eventID',
    'occurrenceID',
    'measurementType',
    'measurementValue',
    'measurementUnit',
    'measurementTypeID',
    'measurementValueID',
    'measurementUnitID',
    'measurementRemarks',
]


def _normalize_template_columns(template_df: pd.DataFrame) -> pd.DataFrame:
    """Normalize user-provided column names to the required set, keeping unknowns for possible future use."""
    col_map = {}
    for col in template_df.columns:
        key = str(col).strip()
        key_lower = key.lower()
        if key_lower in {'measurementtype', 'measurement_type'}:
            col_map[col] = 'measurementType'
        elif key_lower in {'measurementvalue', 'measurement_value'}:
            col_map[col] = 'measurementValue'
        elif key_lower in {'measurementunit', 'measurement_unit'}:
            col_map[col] = 'measurementUnit'
        elif key_lower in {'measurementtypeid', 'measurement_type_id', 'measurementtype_id'}:
            col_map[col] = 'measurementTypeID'
        elif key_lower in {'measurementvalueid', 'measurement_value_id'}:
            col_map[col] = 'measurementValueID'
        elif key_lower in {'measurementunitid', 'measurement_unit_id'}:
            col_map[col] = 'measurementUnitID'
        elif key_lower in {'measurementremarks', 'measurement_remarks'}:
            col_map[col] = 'measurementRemarks'
        elif key_lower in {'linkto', 'link_to', 'scope'}:
            col_map[col] = 'linkTo'
        else:
            col_map[col] = key  # keep as-is
    df = template_df.rename(columns=col_map).copy()

    # Ensure all required columns exist
    for col in ['measurementType', 'measurementValue', 'measurementUnit', 'measurementTypeID', 'measurementValueID', 'measurementUnitID', 'measurementRemarks']:
        if col not in df.columns:
            df[col] = pd.NA
    if 'linkTo' not in df.columns:
        df['linkTo'] = 'occurrence'
    return df


def _get_measurement_unit_from_faire_data(data: dict, measurement_type: str) -> str:
    """Look for a unit column following the pattern {measurementType}_unit in FAIRe metadata."""
    if pd.isna(measurement_type) or not str(measurement_type).strip():
        return pd.NA
    
    # Look for a column named {measurementType}_unit in FAIRe metadata
    unit_col_name = f"{measurement_type}_unit"
    
    # Check sampleMetadata first (most common place for environmental measurements)
    if 'sampleMetadata' in data and not data['sampleMetadata'].empty:
        if unit_col_name in data['sampleMetadata'].columns:
            # Get the first non-NA value
            unit_values = data['sampleMetadata'][unit_col_name].dropna()
            if not unit_values.empty:
                return str(unit_values.iloc[0])
    
    # Check experimentRunMetadata
    if 'experimentRunMetadata' in data and not data['experimentRunMetadata'].empty:
        if unit_col_name in data['experimentRunMetadata'].columns:
            unit_values = data['experimentRunMetadata'][unit_col_name].dropna()
            if not unit_values.empty:
                return str(unit_values.iloc[0])
    
    # Check projectMetadata
    if 'projectMetadata' in data and not data['projectMetadata'].empty:
        if unit_col_name in data['projectMetadata'].columns:
            unit_values = data['projectMetadata'][unit_col_name].dropna()
            if not unit_values.empty:
                return str(unit_values.iloc[0])
    
    return pd.NA


def _load_emof_template(template_path: str) -> pd.DataFrame:
    """Load the eMoF template Excel. Falls back to an empty frame if not found."""
    if not os.path.exists(template_path):
        # Return an empty template; calling code will still be able to generate assay_name rows directly
        return pd.DataFrame(columns=['measurementType', 'measurementValue', 'measurementUnit',
                                     'measurementTypeID', 'measurementValueID', 'measurementUnitID',
                                     'measurementRemarks', 'linkTo'])
    try:
        df = pd.read_excel(template_path, sheet_name=0)
        return df
    except Exception:
        # Try with engine openpyxl fallback (in case of older environments)
        df = pd.read_excel(template_path, sheet_name=0, engine='openpyxl')
        return df


def create_emof_table(params: dict, occurrence_core: pd.DataFrame, data: dict, reporter) -> str:
    """Create the eMoF Excel file and return the output path.

    params:
      - emof_template_path (optional): path to the template Excel
      - output_dir: directory where the TSV will be written
    data:
      - FAIRe metadata dictionary containing sampleMetadata, experimentRunMetadata, etc.
    """
    reporter.add_section("Creating extendedMeasurementOrFact (eMoF)")

    output_dir = params.get('output_dir', 'processed-v3/')
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, 'eMoF.xlsx')

    template_path = params.get('emof_template_path', os.path.join('raw-v3', 'eMoF Fields edna2obis .xlsx'))
    reporter.add_text(f"Using eMoF template: {template_path}")

    # Prepare bases
    if not {'eventID', 'occurrenceID'}.issubset(occurrence_core.columns):
        reporter.add_error("Occurrence core is missing 'eventID' or 'occurrenceID' required for eMoF linking.")
        raise ValueError("Occurrence core missing required identifiers")

    base_occ = occurrence_core[['eventID', 'occurrenceID']].copy()
    if 'assay_name' in occurrence_core.columns:
        base_occ['assay_name'] = occurrence_core['assay_name']
    else:
        base_occ['assay_name'] = pd.NA

    base_events = occurrence_core[['eventID']].drop_duplicates().copy()

    # Load and normalize template
    raw_template = _load_emof_template(template_path)
    template = _normalize_template_columns(raw_template)

    # If template is empty, create a synthetic row to at least export assay_name per occurrence
    if template.empty:
        template = pd.DataFrame([
            {'measurementType': 'assay_name', 'measurementValue': pd.NA, 'measurementUnit': pd.NA,
             'measurementTypeID': pd.NA, 'measurementValueID': pd.NA, 'measurementUnitID': pd.NA,
             'measurementRemarks': pd.NA, 'linkTo': 'occurrence'}
        ])
        reporter.add_warning("Template not found or empty; generating eMoF with only per-occurrence assay_name entries.")

    # Build rows
    emof_rows = []

    for _, row in template.iterrows():
        mtype = str(row.get('measurementType', '')).strip()
        if not mtype or mtype.lower() == 'nan':
            continue

        link_to = str(row.get('linkTo', 'occurrence')).strip().lower() or 'occurrence'

        # Values from template (may be NA)
        mval = row.get('measurementValue', pd.NA)
        # Look for unit in FAIRe data using {measurementType}_unit pattern
        munit = _get_measurement_unit_from_faire_data(data, mtype)
        mtype_id = row.get('measurementTypeID', pd.NA)
        mval_id = row.get('measurementValueID', pd.NA)
        munit_id = row.get('measurementUnitID', pd.NA)
        mremarks = row.get('measurementRemarks', pd.NA)

        if link_to == 'event':
            # Expand across events; occurrenceID left blank
            tmp = base_events.copy()
            tmp['occurrenceID'] = pd.NA
            tmp['measurementType'] = mtype
            tmp['measurementValue'] = mval
            tmp['measurementUnit'] = munit
            tmp['measurementTypeID'] = mtype_id
            tmp['measurementValueID'] = mval_id
            tmp['measurementUnitID'] = munit_id
            tmp['measurementRemarks'] = mremarks
            emof_rows.append(tmp[REQUIRED_COLUMNS_IN_ORDER])
            continue

        # Default: link to occurrence
        if mtype == 'assay_name':
            # If a specific measurementValue is provided, restrict to those occurrences
            occ_subset = base_occ
            if pd.notna(mval) and str(mval).strip():
                occ_subset = base_occ[base_occ['assay_name'].astype(str) == str(mval)]
                value_series = pd.Series(str(mval), index=occ_subset.index)
            else:
                value_series = occ_subset['assay_name']

            if occ_subset.empty:
                continue

            tmp = occ_subset[['eventID', 'occurrenceID']].copy()
            tmp['measurementType'] = 'assay_name'
            tmp['measurementValue'] = value_series
            tmp['measurementUnit'] = munit
            tmp['measurementTypeID'] = mtype_id
            tmp['measurementValueID'] = mval_id
            tmp['measurementUnitID'] = munit_id
            tmp['measurementRemarks'] = mremarks
            emof_rows.append(tmp[REQUIRED_COLUMNS_IN_ORDER])
        else:
            tmp = base_occ[['eventID', 'occurrenceID']].copy()
            tmp['measurementType'] = mtype
            tmp['measurementValue'] = mval
            tmp['measurementUnit'] = munit
            tmp['measurementTypeID'] = mtype_id
            tmp['measurementValueID'] = mval_id
            tmp['measurementUnitID'] = munit_id
            tmp['measurementRemarks'] = mremarks
            emof_rows.append(tmp[REQUIRED_COLUMNS_IN_ORDER])

    if not emof_rows:
        # Ensure we always output something
        result_df = pd.DataFrame(columns=REQUIRED_COLUMNS_IN_ORDER)
    else:
        result_df = pd.concat(emof_rows, ignore_index=True)

    # Ensure correct column order and replace NaNs with empty strings when saving
    for col in REQUIRED_COLUMNS_IN_ORDER:
        if col not in result_df.columns:
            result_df[col] = pd.NA
    result_df = result_df[REQUIRED_COLUMNS_IN_ORDER]

    # Save to Excel
    result_df.to_excel(output_path, index=False, engine='openpyxl')

    reporter.add_success(f"eMoF created successfully with {len(result_df):,} records")
    reporter.add_text(f"Output file: eMoF.xlsx")

    return output_path
