"""
Extended Measurement or Fact (eMoF) Builder for edna2obis

Key behavior (agreed design for metabarcoding datasets):
- Event-level only (for now): one row per eventID per measurement. occurrenceID is left blank.
- Source resolution for each measurementType: look in sampleMetadata first, else experimentRunMetadata, else error.
- Unit rules (strict, non-destructive):
  - If measurementUnit == 'provided' (case-insensitive), a per-row unit column named '<measurementType>_unit' MUST exist in the chosen source sheet and MUST be non-blank for all emitted rows; otherwise error.
  - If measurementUnit is a non-empty literal (e.g., 'm', '°C'), use it as-is for every emitted row for that measurementType.
  - If measurementUnit is blank, leave blank. Do NOT auto-fallback to any unit column. Additionally, if a '<measurementType>_unit' column exists in the source sheet while the template left unit blank, error to avoid unintended data changes.
- Value rules:
  - If template.measurementValue is non-blank: treat as categorical. Emit rows only where the source value equals the template value (trimmed, case-insensitive by default). Otherwise, do not emit for that event.
  - If template.measurementValue is blank: treat as numeric/direct. Copy source value into measurementValue. Skip events where the value is NA (no fabrication).
- Events included: only those present in the final, taxonomically matched occurrence file (occurrence_{api}_matched.csv). This excludes controls/blank samples automatically.
"""

from __future__ import annotations

import os
from typing import Dict, List, Tuple

import pandas as pd


REQUIRED_TEMPLATE_COLUMNS_IN_ORDER: List[str] = [
    'measurementType',
    'measurementValue',
    'measurementUnit',
    'measurementTypeID',
    'measurementValueID',
    'measurementUnitID',
    'measurementRemarks',
]

EMOF_OUTPUT_COLUMNS_IN_ORDER: List[str] = [
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
    final_occurrence_path = os.path.join(output_dir, f"occurrence_{api_choice}_matched.csv")

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


def _prepare_join_frame_for_measurement(
    meas_type: str,
    events_df: pd.DataFrame,
    source_name: str,
    source_df: pd.DataFrame,
    reporter,
) -> pd.DataFrame:
    """Build a per-measurementType join frame limited to included events.

    Returns columns: eventID, value (the measurementType value), optional unit column (<meas_type>_unit) if present.
    """
    # Build the base join between events and the source sheet
    if source_name == 'sampleMetadata':
        if 'samp_name' not in source_df.columns:
            reporter.add_error("sampleMetadata missing required column 'samp_name'.")
            raise ValueError("sampleMetadata missing 'samp_name'")
        join_df = events_df.merge(
            source_df,
            left_on='samp_name',
            right_on='samp_name',
            how='left',
        )
    elif source_name == 'experimentRunMetadata':
        # events_df has eventID (lib_id); join directly to experimentRunMetadata on lib_id
        if 'lib_id' not in source_df.columns:
            reporter.add_error("experimentRunMetadata missing required column 'lib_id'.")
            raise ValueError("experimentRunMetadata missing 'lib_id'")
        join_df = events_df.merge(
            source_df,
            left_on='eventID',
            right_on='lib_id',
            how='left',
        )
    else:
        raise ValueError(f"Unknown source sheet '{source_name}'")

    # Prepare the minimal output
    out = pd.DataFrame()
    out['eventID'] = join_df['eventID']
    out['value'] = join_df[meas_type] if meas_type in join_df.columns else pd.NA

    unit_col = f"{meas_type}_unit"
    if unit_col in join_df.columns:
        out[unit_col] = join_df[unit_col]

    return out


def create_emof_table(params, occurrence_core: pd.DataFrame, data: Dict[str, pd.DataFrame], reporter) -> str:
    """Create an eMoF (extendedMeasurementOrFact) Excel file.

    For this version, all rows are event-linked (occurrenceID left blank). The list of events is
    derived from the final, taxonomically matched occurrence file to ensure controls/blanks are excluded.
    """
    reporter.add_section("Creating eMoF (extendedMeasurementOrFact)")

    try:
        # ------------------------------------------------------------------
        # Load template (input_file sheet)
        # ------------------------------------------------------------------
        template_path = params.get('emof_template_path', 'raw-v3/eMoF Fields edna2obis .xlsx')
        if not os.path.exists(template_path):
            reporter.add_error(f"eMoF template not found: {template_path}")
            raise FileNotFoundError(f"Template not found: {template_path}")

        try:
            template_df = pd.read_excel(template_path, sheet_name='input_file', dtype=object)
        except Exception as e:
            reporter.add_error(f"Failed to read 'input_file' sheet from template: {e}")
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
        # Load final occurrence file to determine events to include
        # ------------------------------------------------------------------
        final_occ_df = _load_final_occurrence(params, reporter)
        if final_occ_df.empty and isinstance(occurrence_core, pd.DataFrame) and not occurrence_core.empty:
            # Fallback to in-memory occurrence_core (should already have controls removed)
            reporter.add_warning("Using in-memory occurrence_core to derive events (final occurrence file unavailable)")
            final_occ_df = occurrence_core.copy()

        events_df = _build_event_dataframe(final_occ_df, data, reporter)
        reporter.add_text(f"Total events included (non-control): {len(events_df)}")

        # ------------------------------------------------------------------
        # Build eMoF rows
        # ------------------------------------------------------------------
        emof_rows: List[Dict[str, object]] = []

        # Cache join frames per measurementType to avoid repeated merges
        prepared_frames: Dict[str, pd.DataFrame] = {}
        prepared_sources: Dict[str, Tuple[str, pd.DataFrame]] = {}

        # Iterate template rows in order
        for _, trow in template_df.iterrows():
            meas_type = _normalize_str(trow.get('measurementType'))
            templ_value = trow.get('measurementValue')
            templ_unit = trow.get('measurementUnit')
            templ_mtid = trow.get('measurementTypeID')
            templ_mvid = trow.get('measurementValueID')
            templ_muid = trow.get('measurementUnitID')
            templ_rem = trow.get('measurementRemarks')

            if meas_type == '':
                # Skip defensive; should have been filtered out
                continue

            # Resolve data source for the measurementType
            if meas_type not in prepared_sources:
                try:
                    source_name, source_df = _resolve_source_for_measurement(meas_type, data)
                except Exception as e:
                    reporter.add_error(str(e))
                    raise
                prepared_sources[meas_type] = (source_name, source_df)
            else:
                source_name, source_df = prepared_sources[meas_type]

            # Prepare the limited join frame (eventID + value + optional unit column)
            if meas_type not in prepared_frames:
                join_frame = _prepare_join_frame_for_measurement(meas_type, events_df, source_name, source_df, reporter)
                # Normalize value and potential unit strings
                # Do not coerce numeric here; we preserve formatting as provided (non-destructive)
                if 'value' in join_frame.columns:
                    # Keep original dtypes; but for categorical compare, we'll normalize per-comparison
                    pass
                prepared_frames[meas_type] = join_frame
            else:
                join_frame = prepared_frames[meas_type]

            # Unit policy checks
            templ_unit_norm = _normalize_str(templ_unit)
            has_per_row_unit_col = f"{meas_type}_unit" in join_frame.columns

            if templ_unit_norm == 'provided':
                if not has_per_row_unit_col:
                    reporter.add_error(
                        f"For measurementType '{meas_type}', measurementUnit='provided' but column '{meas_type}_unit' "
                        f"was not found in the source sheet ('{source_name}')."
                    )
                    raise ValueError(f"Missing per-row unit column for '{meas_type}'")
            else:
                if templ_unit_norm == '':
                    # Blank unit: do not fallback; but if a per-row unit column exists, error to avoid silent changes
                    if has_per_row_unit_col:
                        reporter.add_error(
                            f"For measurementType '{meas_type}', the template left measurementUnit blank but a per-row unit column "
                            f"'{meas_type}_unit' exists in the source data. Set measurementUnit to 'provided' or a literal."
                        )
                        raise ValueError("Ambiguous unit policy detected (blank vs provided)")
                else:
                    # Literal unit: use as-is
                    pass

            # Determine categorical vs numeric/direct by template value emptiness
            is_categorical = _normalize_str(templ_value) != ''

            # Iterate each event row and decide whether to emit
            for _, ev in join_frame.iterrows():
                event_id = _normalize_str(ev.get('eventID'))
                src_val = ev.get('value')

                if pd.isna(src_val) or _normalize_str(src_val) == '':
                    # Skip events without a value (non-destructive)
                    continue

                if is_categorical:
                    # Emit only when src_val matches template value (trimmed, case-insensitive)
                    if not _case_insensitive_equal(src_val, templ_value):
                        continue
                    out_value = _normalize_str(templ_value)
                    out_value_id = _normalize_str(templ_mvid)
                else:
                    # Numeric/direct: copy the source value as-is (stringify for Excel safety but avoid altering)
                    out_value = src_val
                    out_value_id = ''

                # Resolve unit output per policy
                if templ_unit_norm == 'provided':
                    unit_col = f"{meas_type}_unit"
                    unit_val = ev.get(unit_col)
                    if pd.isna(unit_val) or _normalize_str(unit_val) == '':
                        reporter.add_error(
                            f"For measurementType '{meas_type}', measurementUnit='provided' but unit is blank for eventID '{event_id}'."
                        )
                        raise ValueError("Blank per-row unit encountered under 'provided' policy")
                    out_unit = unit_val
                elif templ_unit_norm == '':
                    out_unit = ''
                else:
                    out_unit = templ_unit_norm

                # Compose eMoF row (occurrenceID left blank per agreed event-only design)
                emof_rows.append({
                    'eventID': event_id,
                    'occurrenceID': '',
                    'measurementType': meas_type,
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

        # Deterministic sort: by eventID, measurementType, then measurementValue
        sort_cols = [c for c in ['eventID', 'measurementType', 'measurementValue'] if c in emof_df.columns]
        if sort_cols:
            emof_df = emof_df.sort_values(sort_cols).reset_index(drop=True)

        output_dir = params.get('output_dir', 'processed-v3/')
        os.makedirs(output_dir, exist_ok=True)
        
        # Write only CSV file
        emof_csv_path = os.path.join(output_dir, 'eMoF.csv')
        try:
            emof_df.to_csv(emof_csv_path, index=False)
        except Exception as e:
            reporter.add_warning(f"Could not write eMoF CSV ('{emof_csv_path}'): {e}")

        reporter.add_success(f"eMoF created successfully with {len(emof_df)} row(s)")
        reporter.add_text(f"Saved eMoF CSV to: {emof_csv_path}")
        reporter.add_dataframe(emof_df.head(15), "eMoF Preview (first 15 rows)")

        return emof_csv_path

    except Exception as e:
        reporter.add_error(f"eMoF creation failed: {e}")
        raise


  