"""
Extended Measurement or Fact (eMoF) Builder for edna2obis

Key behavior (occurrence-based design for metabarcoding datasets):
- Occurrence-based only: one row per occurrenceID per measurement. This creates a larger file but is what GBIF requires.
- Source resolution for each measurementType: look in sampleMetadata first, else experimentRunMetadata, else error.
- Unit rules (strict, non-destructive):
  - If measurementUnit == 'provided' (case-insensitive), a per-row unit column named '<measurementType>_unit' MUST exist in the chosen source sheet and MUST be non-blank for all emitted rows; otherwise error.
  - If measurementUnit is a non-empty literal (e.g., 'm', '°C'), use it as-is for every emitted row for that measurementType.
  - If measurementUnit is blank, leave blank. Do NOT auto-fallback to any unit column. Additionally, if a '<measurementType>_unit' column exists in the source sheet while the template left unit blank, error to avoid unintended data changes.
- Value rules:
  - If template.measurementValue is non-blank: treat as categorical. Emit rows only where the source value equals the template value (trimmed, case-insensitive by default). Otherwise, do not emit for that occurrence.
  - If template.measurementValue is blank: treat as numeric/direct. Copy source value into measurementValue. Skip occurrences where the value is NA (no fabrication).
- Occurrences included: only those present in the final, taxonomically matched occurrence file (occurrence_{api}_matched.csv). This excludes controls/blank samples automatically.
"""

from __future__ import annotations

import os
from typing import Dict, List, Tuple

import pandas as pd


REQUIRED_TEMPLATE_COLUMNS_IN_ORDER: List[str] = [
    'measurementType',
    'measurementValue',
    'measurementUnit',
]

EMOF_OUTPUT_COLUMNS_IN_ORDER: List[str] = [
    'eventID',
    'occurrenceID',
    'verbatimMeasurementType',
    'measurementType',
    'measurementValue',
    'measurementUnit',
    'measurementTypeID',
    'measurementValueID',
    'measurementUnitID',
    'measurementRemarks',
]


def _normalize_str(val) -> str:
    if pd.isna(val):
        return ''
    return str(val).strip()


def _case_insensitive_equal(a, b) -> bool:
    return _normalize_str(a).lower() == _normalize_str(b).lower()


def _load_final_occurrence(params, reporter) -> pd.DataFrame:
    """Load the final, taxonomically matched occurrence file.

    Falls back to the intermediate in-memory occurrence_core DataFrame only if the
    final file is missing (should be rare), and warns the user.
    """
    output_dir = params.get('output_dir', 'processed-v3/')
    api_choice = params.get('taxonomic_api_source', 'WoRMS').lower()
    final_occurrence_path = os.path.join(output_dir, f"occurrence_core_{api_choice}.csv")

    if os.path.exists(final_occurrence_path):
        try:
            df = pd.read_csv(final_occurrence_path, low_memory=False)
            return df
        except Exception as e:
            reporter.add_warning(f"Could not read final occurrence file '{final_occurrence_path}': {e}")
    else:
        reporter.add_warning(
            f"Final occurrence file not found at '{final_occurrence_path}'. "
            f"Will attempt to use in-memory occurrence_core as fallback."
        )
    return pd.DataFrame()


def _build_event_dataframe(final_occ_df: pd.DataFrame, data: Dict[str, pd.DataFrame], reporter) -> pd.DataFrame:
    """Create a DataFrame of distinct eventIDs to include in eMoF and map to samp_name.

    Returns a DataFrame with columns:
      - eventID (lib_id)
      - samp_name (mapped from experimentRunMetadata)

    Assumes controls/blank samples have already been removed upstream from the
    provided 'data' dict (sampleMetadata, experimentRunMetadata).
    """
    required_cols = []
    if not final_occ_df.empty:
        if 'eventID' not in final_occ_df.columns:
            reporter.add_error("Final occurrence file is missing required column 'eventID'.")
            raise ValueError("Missing 'eventID' in final occurrence file")
        unique_events = (
            final_occ_df[['eventID']]
            .dropna()
            .drop_duplicates()
            .rename(columns={'eventID': 'eventID'})
        )
    else:
        # If we cannot load final file, derive eventIDs from experimentRunMetadata (non-ideal)
        reporter.add_warning(
            "Using experimentRunMetadata to derive eventIDs because final occurrence file was unavailable."
        )
        if 'experimentRunMetadata' not in data or data['experimentRunMetadata'].empty:
            reporter.add_error("experimentRunMetadata is missing or empty; cannot derive eventIDs.")
            raise ValueError("experimentRunMetadata missing for eMoF event derivation")
        erm = data['experimentRunMetadata']
        if 'lib_id' not in erm.columns:
            reporter.add_error("experimentRunMetadata missing required column 'lib_id'.")
            raise ValueError("experimentRunMetadata missing 'lib_id'")
        unique_events = erm[['lib_id']].dropna().drop_duplicates().rename(columns={'lib_id': 'eventID'})

    # Map eventID (lib_id) to samp_name using experimentRunMetadata
    if 'experimentRunMetadata' not in data or data['experimentRunMetadata'].empty:
        reporter.add_error("experimentRunMetadata is missing or empty; cannot map eventID to samp_name.")
        raise ValueError("experimentRunMetadata missing for mapping")

    erm_df = data['experimentRunMetadata'].copy()
    # Normalize lib_id and samp_name to strings
    for col in ['lib_id', 'samp_name']:
        if col in erm_df.columns:
            erm_df[col] = erm_df[col].astype(str).str.strip()

    if 'lib_id' not in erm_df.columns or 'samp_name' not in erm_df.columns:
        reporter.add_error("experimentRunMetadata must contain 'lib_id' and 'samp_name'.")
        raise ValueError("experimentRunMetadata missing 'lib_id' or 'samp_name'")

    events = unique_events.merge(erm_df[['lib_id', 'samp_name']].drop_duplicates('lib_id'),
                                 left_on='eventID', right_on='lib_id', how='left')
    events.drop(columns=['lib_id'], inplace=True)

    # Sanity check: any eventID without samp_name will limit access to sampleMetadata
    missing_samp = events['samp_name'].isna().sum()
    if missing_samp:
        reporter.add_warning(
            f"{missing_samp} eventID(s) could not be mapped to samp_name via experimentRunMetadata. "
            f"Measurements that require sampleMetadata will be missing for those events."
        )

    return events


def _resolve_source_for_measurement(meas_type: str, data: Dict[str, pd.DataFrame]) -> Tuple[str, pd.DataFrame]:
    """Determine which sheet holds the measurementType column.

    Returns a tuple (source_name, source_df) where source_name is one of
    'sampleMetadata' or 'experimentRunMetadata'. Raises on failure.
    """
    meas_type = str(meas_type)
    if 'sampleMetadata' in data and not data['sampleMetadata'].empty and meas_type in data['sampleMetadata'].columns:
        return 'sampleMetadata', data['sampleMetadata']
    if 'experimentRunMetadata' in data and not data['experimentRunMetadata'].empty and meas_type in data['experimentRunMetadata'].columns:
        return 'experimentRunMetadata', data['experimentRunMetadata']
    raise ValueError(f"measurementType '{meas_type}' not found in sampleMetadata or experimentRunMetadata")


def _prepare_join_frame_for_occurrence_measurement(
    meas_type: str,
    occurrences_df: pd.DataFrame,
    source_name: str,
    source_df: pd.DataFrame,
    data: Dict[str, pd.DataFrame],
    reporter,
) -> pd.DataFrame:
    """Build a per-measurementType join frame limited to included occurrences.

    Returns columns: occurrenceID, eventID, value (the measurementType value), optional unit column (<meas_type>_unit) if present.
    """
    # Build the base join between occurrences and the source sheet
    if source_name == 'sampleMetadata':
        if 'samp_name' not in source_df.columns:
            reporter.add_error("sampleMetadata missing required column 'samp_name'.")
            raise ValueError("sampleMetadata missing 'samp_name'")
        
        # First, we need to map occurrenceID -> eventID -> samp_name
        # Get samp_name from experimentRunMetadata using eventID
        if 'experimentRunMetadata' not in data or data['experimentRunMetadata'].empty:
            reporter.add_error("experimentRunMetadata is missing; cannot map occurrenceID to samp_name for sampleMetadata measurements.")
            raise ValueError("experimentRunMetadata missing for occurrence-to-samp_name mapping")
        
        erm_df = data['experimentRunMetadata'].copy()
        for col in ['lib_id', 'samp_name']:
            if col in erm_df.columns:
                erm_df[col] = erm_df[col].astype(str).str.strip()
        
        # Join occurrences with experimentRunMetadata to get samp_name
        occurrences_with_samp = occurrences_df.merge(
            erm_df[['lib_id', 'samp_name']].drop_duplicates('lib_id'),
            left_on='eventID',
            right_on='lib_id',
            how='left'
        )
        
        # Now join with sampleMetadata using samp_name
        join_df = occurrences_with_samp.merge(
            source_df,
            left_on='samp_name',
            right_on='samp_name',
            how='left',
        )
        
    elif source_name == 'experimentRunMetadata':
        # occurrences_df has eventID; join directly to experimentRunMetadata on lib_id
        if 'lib_id' not in source_df.columns:
            reporter.add_error("experimentRunMetadata missing required column 'lib_id'.")
            raise ValueError("experimentRunMetadata missing 'lib_id'")
        join_df = occurrences_df.merge(
            source_df,
            left_on='eventID',
            right_on='lib_id',
            how='left',
        )
    else:
        raise ValueError(f"Unknown source sheet '{source_name}'")

    # Prepare the minimal output
    out = pd.DataFrame()
    out['occurrenceID'] = join_df['occurrenceID']
    out['eventID'] = join_df['eventID']
    out['value'] = join_df[meas_type] if meas_type in join_df.columns else pd.NA

    unit_col = f"{meas_type}_unit"
    if unit_col in join_df.columns:
        out[unit_col] = join_df[unit_col]

    return out


# COMMENTED OUT: Event-based eMoF function (not used since GBIF requires occurrence-based)
# def _prepare_join_frame_for_measurement(
#     meas_type: str,
#     events_df: pd.DataFrame,
#     source_name: str,
#     source_df: pd.DataFrame,
#     reporter,
# ) -> pd.DataFrame:
#     """Build a per-measurementType join frame limited to included events.
# 
#     Returns columns: eventID, value (the measurementType value), optional unit column (<meas_type>_unit) if present.
#     """
#     # Build the base join between events and the source sheet
#     if source_name == 'sampleMetadata':
#         if 'samp_name' not in source_df.columns:
#             reporter.add_error("sampleMetadata missing required column 'samp_name'.")
#             raise ValueError("sampleMetadata missing 'samp_name'")
#         join_df = events_df.merge(
#             source_df,
#             left_on='samp_name',
#             right_on='samp_name',
#             how='left',
#         )
#     elif source_name == 'experimentRunMetadata':
#         # events_df has eventID (lib_id); join directly to experimentRunMetadata on lib_id
#         if 'lib_id' not in source_df.columns:
#             reporter.add_error("experimentRunMetadata missing required column 'lib_id'.")
#             raise ValueError("experimentRunMetadata missing 'lib_id'")
#         join_df = events_df.merge(
#             source_df,
#             left_on='eventID',
#             right_on='lib_id',
#             how='left',
#         )
#     else:
#         raise ValueError(f"Unknown source sheet '{source_name}'")
# 
#     # Prepare the minimal output
#     out = pd.DataFrame()
#     out['eventID'] = join_df['eventID']
#     out['value'] = join_df[meas_type] if meas_type in join_df.columns else pd.NA
# 
#     unit_col = f"{meas_type}_unit"
#     if unit_col in join_df.columns:
#         out[unit_col] = join_df[unit_col]
# 
#     return out


def create_emof_table(params, occurrence_core: pd.DataFrame, data: Dict[str, pd.DataFrame], reporter) -> str:
    """Create an eMoF (extendedMeasurementOrFact) CSV file.

    Creates occurrence-based eMoF: one row per occurrenceID per measurement.
    This creates a larger file but is what GBIF requires for data submission.
    """
    reporter.add_section("Creating eMoF (extendedMeasurementOrFact)")

    try:
        # ------------------------------------------------------------------
        # Load template (input_file sheet)
        # ------------------------------------------------------------------
        template_path = params.get('emof_template_path', 'raw-v3/eMoF_Fields_edna2obis.tsv')
        if not os.path.exists(template_path):
            reporter.add_error(f"eMoF template not found: {template_path}")
            raise FileNotFoundError(f"Template not found: {template_path}")

        try:
            template_df = pd.read_csv(template_path, sep='\t', dtype=object)
        except Exception as e:
            reporter.add_error(f"Failed to read eMoF template file: {e}")
            raise

        # Normalize header names (but enforce exact required set)
        template_df.columns = [str(c).strip() for c in template_df.columns]
        missing_cols = [c for c in REQUIRED_TEMPLATE_COLUMNS_IN_ORDER if c not in template_df.columns]
        if missing_cols:
            reporter.add_error(f"Template is missing required columns: {missing_cols}")
            raise ValueError("Template missing required columns")

        # Drop completely empty rows (where measurementType is NA/blank)
        template_df['measurementType'] = template_df['measurementType'].astype(object)
        template_df = template_df[template_df['measurementType'].apply(lambda x: _normalize_str(x) != '')].copy()

        reporter.add_text(f"Loaded template with {len(template_df)} configured measurement row(s)")

        # ------------------------------------------------------------------
        # Load final occurrence file to determine occurrences to include
        # ------------------------------------------------------------------
        final_occ_df = _load_final_occurrence(params, reporter)
        if final_occ_df.empty and isinstance(occurrence_core, pd.DataFrame) and not occurrence_core.empty:
            # Fallback to in-memory occurrence_core (should already have controls removed)
            reporter.add_warning("Using in-memory occurrence_core to derive occurrences (final occurrence file unavailable)")
            final_occ_df = occurrence_core.copy()

        # Use all occurrences from the final occurrence file
        if final_occ_df.empty:
            reporter.add_error("No occurrence data available for eMoF")
            raise ValueError("No occurrence data available")
        
        required_occ_cols = ['occurrenceID', 'eventID']
        missing_occ_cols = [c for c in required_occ_cols if c not in final_occ_df.columns]
        if missing_occ_cols:
            reporter.add_error(f"Final occurrence file missing required columns: {missing_occ_cols}")
            raise ValueError(f"Missing occurrence columns: {missing_occ_cols}")
        
        occurrences_df = final_occ_df[required_occ_cols].dropna().drop_duplicates()
        reporter.add_text(f"Total occurrences included (non-control): {len(occurrences_df)}")

        # ------------------------------------------------------------------
        # Build eMoF rows
        # ------------------------------------------------------------------
        emof_rows: List[Dict[str, object]] = []

        # Cache join frames per measurementType to avoid repeated merges
        prepared_frames: Dict[str, pd.DataFrame] = {}
        prepared_sources: Dict[str, Tuple[str, pd.DataFrame]] = {}

        has_verbatim_col = 'verbatimMeasurementType' in template_df.columns

        # Iterate template rows in order
        for _, trow in template_df.iterrows():
            output_meas_type = _normalize_str(trow.get('measurementType'))
            templ_value = trow.get('measurementValue')
            templ_unit = trow.get('measurementUnit')
            templ_mtid = trow.get('measurementTypeID')
            templ_mvid = trow.get('measurementValueID')
            templ_muid = trow.get('measurementUnitID')
            templ_rem = trow.get('measurementRemarks')
            verbatim_raw = trow.get('verbatimMeasurementType') if has_verbatim_col else ''

            source_field = output_meas_type
            if has_verbatim_col:
                verbatim_field = _normalize_str(trow.get('verbatimMeasurementType'))
                if verbatim_field:
                    source_field = verbatim_field

            if output_meas_type == '':
                # Skip defensive; should have been filtered out
                continue

            # Resolve data source for the measurementType
            if source_field not in prepared_sources:
                try:
                    source_name, source_df = _resolve_source_for_measurement(source_field, data)
                except Exception as e:
                    # For master-list behavior: if a verbatimMeasurementType was specified but
                    # the corresponding source field is missing in this dataset's metadata,
                    # skip emitting this measurement for this dataset without erroring out.
                    if has_verbatim_col and _normalize_str(trow.get('verbatimMeasurementType')):
                        reporter.add_text(
                            f"Skipping measurementType '{output_meas_type}': source field '{source_field}' not found in metadata for this dataset."
                        )
                        continue
                    reporter.add_error(f"ERROR: {str(e)}")
                    raise
                prepared_sources[source_field] = (source_name, source_df)
            else:
                source_name, source_df = prepared_sources[source_field]

            # Prepare the limited join frame for occurrence-based measurements
            if source_field not in prepared_frames:
                join_frame = _prepare_join_frame_for_occurrence_measurement(
                    source_field, occurrences_df, source_name, source_df, data, reporter
                )
                prepared_frames[source_field] = join_frame
            else:
                join_frame = prepared_frames[source_field]

            # Unit policy checks
            templ_unit_norm = _normalize_str(templ_unit)
            has_per_row_unit_col = f"{source_field}_unit" in join_frame.columns

            if templ_unit_norm == 'provided':
                if not has_per_row_unit_col:
                    reporter.add_error(
                        f"For measurementType '{output_meas_type}', measurementUnit='provided' but column '{source_field}_unit' "
                        f"was not found in the source sheet ('{source_name}')."
                    )
                    raise ValueError(f"Missing per-row unit column for '{source_field}'")
            else:
                if templ_unit_norm == '':
                    # Blank unit: do not fallback; but if a per-row unit column exists, error to avoid silent changes
                    if has_per_row_unit_col:
                        reporter.add_error(
                            f"For measurementType '{output_meas_type}', the template left measurementUnit blank but a per-row unit column "
                            f"'{source_field}_unit' exists in the source data. Set measurementUnit to 'provided' or a literal."
                        )
                        raise ValueError("Ambiguous unit policy detected (blank vs provided)")
                else:
                    # Literal unit: use as-is
                    pass

            # Determine categorical vs numeric/direct by template value emptiness
            is_categorical = _normalize_str(templ_value) != ''

            # Iterate each occurrence row and decide whether to emit
            for _, row in join_frame.iterrows():
                occurrence_id = _normalize_str(row.get('occurrenceID'))
                event_id = _normalize_str(row.get('eventID'))
                src_val = row.get('value')

                if pd.isna(src_val) or _normalize_str(src_val) == '':
                    # Skip occurrences without a value (non-destructive)
                    continue

                if is_categorical:
                    # Emit only when src_val matches template value (trimmed, case-insensitive)
                    if not _case_insensitive_equal(src_val, templ_value):
                        continue
                    out_value = _normalize_str(templ_value)
                    out_value_id = _normalize_str(templ_mvid)
                else:
                    # Numeric/direct: copy the source value as-is
                    out_value = src_val
                    out_value_id = ''

                # Resolve unit output per policy
                if templ_unit_norm == 'provided':
                    unit_col = f"{source_field}_unit"
                    unit_val = row.get(unit_col)
                    if pd.isna(unit_val) or _normalize_str(unit_val) == '':
                        reporter.add_error(
                            f"For measurementType '{output_meas_type}', measurementUnit='provided' but unit is blank for occurrenceID '{occurrence_id}'."
                        )
                        raise ValueError("Blank per-row unit encountered under 'provided' policy")
                    out_unit = unit_val
                elif templ_unit_norm == '':
                    out_unit = ''
                else:
                    out_unit = templ_unit_norm

                # Compose eMoF row
                emof_rows.append({
                    'eventID': event_id,
                    'occurrenceID': occurrence_id,
                    'measurementType': output_meas_type,
                    'verbatimMeasurementType': verbatim_raw,
                    'measurementValue': out_value,
                    'measurementUnit': out_unit,
                    'measurementTypeID': _normalize_str(templ_mtid),
                    'measurementValueID': out_value_id,
                    'measurementUnitID': _normalize_str(templ_muid),
                    'measurementRemarks': _normalize_str(templ_rem),
                })

        # ------------------------------------------------------------------
        # Finalize and save
        # ------------------------------------------------------------------
        if not emof_rows:
            reporter.add_warning("No eMoF rows were generated based on the template and available data.")

        emof_df = pd.DataFrame(emof_rows, columns=EMOF_OUTPUT_COLUMNS_IN_ORDER)

        # Deterministic sort: by occurrenceID, measurementType, then measurementValue
        sort_cols = [c for c in ['occurrenceID', 'measurementType', 'measurementValue'] if c in emof_df.columns]
        if sort_cols:
            # Convert sort columns to string to prevent mixed-type errors during sort
            for col in sort_cols:
                emof_df[col] = emof_df[col].astype(str)
            emof_df = emof_df.sort_values(sort_cols).reset_index(drop=True)

        output_dir = params.get('output_dir', 'processed-v3/')
        os.makedirs(output_dir, exist_ok=True)
        
        # Write only CSV file
        emof_csv_path = os.path.join(output_dir, 'eMoF.csv')
        try:
            emof_df.to_csv(emof_csv_path, index=False, encoding='utf-8-sig') # Encoding helps with special characters in units
        except Exception as e:
            reporter.add_warning(f"Could not write eMoF CSV ('{emof_csv_path}'): {e}")

        reporter.add_success(f"eMoF created successfully with {len(emof_df)} row(s) (occurrence-based)")
        reporter.add_text(f"Saved eMoF CSV to: {emof_csv_path}")
        reporter.add_text(f"File size: {len(emof_df)} rows (one per occurrence per measurement)")
        reporter.add_dataframe(emof_df.head(15), "eMoF Preview (first 15 rows)")

        return emof_csv_path

    except Exception as e:
        reporter.add_error(f"eMoF creation failed: {e}")
        raise


  