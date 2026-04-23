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
import csv
import heapq
import tempfile
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

# Legacy function, not used anymore. Still remains here just in case we need it later.
def _external_sort_emof_csv(path: str, chunk_size: int = 200000) -> None:
    temp_chunk_paths: List[str] = []
    temp_readers = []
    temp_handles = []

    def sort_key(row: Dict[str, object]) -> Tuple[str, str, str]:
        return (
            str(row.get('occurrenceID', '')),
            str(row.get('measurementType', '')),
            str(row.get('measurementValue', '')),
        )

    try:
        with open(path, 'r', newline='', encoding='utf-8-sig') as src_f:
            reader = csv.DictReader(src_f)
            fieldnames = list(reader.fieldnames or EMOF_OUTPUT_COLUMNS_IN_ORDER)

            chunk: List[Dict[str, object]] = []
            for row in reader:
                chunk.append(row)
                if len(chunk) >= chunk_size:
                    chunk.sort(key=sort_key)
                    with tempfile.NamedTemporaryFile('w', delete=False, newline='', encoding='utf-8', suffix='.tmp', dir=os.path.dirname(path)) as tf:
                        w = csv.DictWriter(tf, fieldnames=fieldnames)
                        w.writeheader()
                        w.writerows(chunk)
                        temp_chunk_paths.append(tf.name)
                    chunk = []

            if chunk:
                chunk.sort(key=sort_key)
                with tempfile.NamedTemporaryFile('w', delete=False, newline='', encoding='utf-8', suffix='.tmp', dir=os.path.dirname(path)) as tf:
                    w = csv.DictWriter(tf, fieldnames=fieldnames)
                    w.writeheader()
                    w.writerows(chunk)
                    temp_chunk_paths.append(tf.name)

        if not temp_chunk_paths:
            return

        out_sorted = f"{path}.sorted"
        heap = []

        for idx, p in enumerate(temp_chunk_paths):
            h = open(p, 'r', newline='', encoding='utf-8')
            temp_handles.append(h)
            r = csv.DictReader(h)
            temp_readers.append(r)
            first = next(r, None)
            if first is not None:
                heapq.heappush(heap, (sort_key(first), idx, first))

        with open(out_sorted, 'w', newline='', encoding='utf-8-sig') as out_f:
            w = csv.DictWriter(out_f, fieldnames=fieldnames)
            w.writeheader()
            while heap:
                _, idx, row = heapq.heappop(heap)
                w.writerow(row)
                nxt = next(temp_readers[idx], None)
                if nxt is not None:
                    heapq.heappush(heap, (sort_key(nxt), idx, nxt))

        os.replace(out_sorted, path)
    finally:
        for h in temp_handles:
            try:
                h.close()
            except Exception:
                pass
        for p in temp_chunk_paths:
            try:
                os.remove(p)
            except Exception:
                pass
        sorted_tmp = f"{path}.sorted"
        if os.path.exists(sorted_tmp):
            try:
                os.remove(sorted_tmp)
            except Exception:
                pass


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

        # Pre-sort occurrences by stripped occurrenceID so streaming emission is already
        # grouped by occurrenceID in the final file. This lets us skip a huge external sort.
        occurrences_df = occurrences_df.copy()
        occurrences_df['occurrenceID'] = occurrences_df['occurrenceID'].astype(str).str.strip()
        occurrences_df['eventID'] = occurrences_df['eventID'].astype(str).str.strip()
        occurrences_df = occurrences_df.sort_values('occurrenceID', kind='mergesort').reset_index(drop=True)
        N_occ = len(occurrences_df)
        occ_ids_arr = occurrences_df['occurrenceID'].values
        evt_ids_arr = occurrences_df['eventID'].values

        # Stable-sort template rows by stripped measurementType so within an occurrence,
        # measurements appear in a deterministic, readable order.
        template_df = template_df.copy()
        template_df['_mt_key'] = template_df['measurementType'].apply(_normalize_str)
        template_df = template_df.sort_values('_mt_key', kind='mergesort').drop(columns=['_mt_key'])

        has_verbatim_col = 'verbatimMeasurementType' in template_df.columns
        output_dir = params.get('output_dir', 'processed-v3/')
        os.makedirs(output_dir, exist_ok=True)
        emof_csv_path = os.path.join(output_dir, 'eMoF.csv')
        preview_rows: List[Dict[str, object]] = []
        row_count = 0

        prepared_sources: Dict[str, Tuple[str, pd.DataFrame]] = {}
        field_value_arr: Dict[str, object] = {}
        field_unit_arr: Dict[str, object] = {}
        field_has_unit: Dict[str, bool] = {}
        template_infos: List[Dict[str, object]] = []

        def _normalize_na_inplace(arr):
            """Replace NaN / pd.NA / NaT / blank-string entries in-place with None."""
            for j in range(len(arr)):
                v = arr[j]
                if v is None:
                    continue
                if isinstance(v, float):
                    if v != v:
                        arr[j] = None
                        continue
                elif isinstance(v, str):
                    if v.strip() == '':
                        arr[j] = None
                        continue
                else:
                    try:
                        if pd.isna(v):
                            arr[j] = None
                            continue
                    except Exception:
                        pass

        for _, trow in template_df.iterrows():
            output_meas_type = _normalize_str(trow.get('measurementType'))
            if output_meas_type == '':
                continue

            templ_value = trow.get('measurementValue')
            templ_unit = trow.get('measurementUnit')
            templ_mtid_norm = _normalize_str(trow.get('measurementTypeID'))
            templ_mvid_norm = _normalize_str(trow.get('measurementValueID'))
            templ_muid_norm = _normalize_str(trow.get('measurementUnitID'))
            templ_rem_norm = _normalize_str(trow.get('measurementRemarks'))
            verbatim_raw = trow.get('verbatimMeasurementType') if has_verbatim_col else ''

            source_field = output_meas_type
            if has_verbatim_col:
                verbatim_field = _normalize_str(trow.get('verbatimMeasurementType'))
                if verbatim_field:
                    source_field = verbatim_field

            # Resolve source once per source_field
            if source_field not in prepared_sources:
                try:
                    source_name, source_df = _resolve_source_for_measurement(source_field, data)
                except Exception as e:
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

            # Build aligned value/unit arrays once per source_field
            if source_field not in field_value_arr:
                jf = _prepare_join_frame_for_occurrence_measurement(
                    source_field, occurrences_df, source_name, source_df, data, reporter
                )
                if len(jf) != N_occ:
                    pos_df = pd.DataFrame({
                        'occurrenceID': occurrences_df['occurrenceID'],
                        'eventID': occurrences_df['eventID'],
                        '_pos': range(N_occ),
                    })
                    jf = jf.merge(pos_df, on=['occurrenceID', 'eventID'], how='right')
                    jf = jf.sort_values('_pos').reset_index(drop=True).drop(columns=['_pos'])
                    if len(jf) != N_occ:
                        raise ValueError(
                            f"Internal alignment error for source field '{source_field}': join frame length {len(jf)} != {N_occ}"
                        )

                val_arr = jf['value'].astype(object).values.copy()
                _normalize_na_inplace(val_arr)
                field_value_arr[source_field] = val_arr

                unit_col = f"{source_field}_unit"
                if unit_col in jf.columns:
                    unit_arr = jf[unit_col].astype(object).values.copy()
                    _normalize_na_inplace(unit_arr)
                    field_unit_arr[source_field] = unit_arr
                    field_has_unit[source_field] = True
                else:
                    field_unit_arr[source_field] = None
                    field_has_unit[source_field] = False

            has_per_row_unit_col = field_has_unit[source_field]

            # Unit policy checks (identical semantics to previous per-row logic)
            templ_unit_norm = _normalize_str(templ_unit)
            if templ_unit_norm == 'provided':
                if not has_per_row_unit_col:
                    reporter.add_error(
                        f"For measurementType '{output_meas_type}', measurementUnit='provided' but column '{source_field}_unit' "
                        f"was not found in the source sheet ('{source_name}')."
                    )
                    raise ValueError(f"Missing per-row unit column for '{source_field}'")
            elif templ_unit_norm == '':
                if has_per_row_unit_col:
                    reporter.add_error(
                        f"For measurementType '{output_meas_type}', the template left measurementUnit blank but a per-row unit column "
                        f"'{source_field}_unit' exists in the source data. Set measurementUnit to 'provided' or a literal."
                    )
                    raise ValueError("Ambiguous unit policy detected (blank vs provided)")

            templ_value_norm = _normalize_str(templ_value)
            is_categorical = templ_value_norm != ''
            templ_value_norm_lower = templ_value_norm.lower() if is_categorical else ''

            template_infos.append({
                'output_meas_type': output_meas_type,
                'verbatim_raw': verbatim_raw,
                'source_field': source_field,
                'is_categorical': is_categorical,
                'templ_value_norm': templ_value_norm,
                'templ_value_norm_lower': templ_value_norm_lower,
                'templ_unit_norm': templ_unit_norm,
                'templ_mtid_norm': templ_mtid_norm,
                'templ_mvid_norm': templ_mvid_norm,
                'templ_muid_norm': templ_muid_norm,
                'templ_rem_norm': templ_rem_norm,
            })

        # --- Stream emission: occurrence-major, already sorted. No external sort needed. ---
        with open(emof_csv_path, 'w', newline='', encoding='utf-8-sig') as emof_f:
            writer = csv.writer(emof_f)
            writer.writerow(EMOF_OUTPUT_COLUMNS_IN_ORDER)

            for i in range(N_occ):
                occ_s = occ_ids_arr[i]
                evt_s = evt_ids_arr[i]

                for ti in template_infos:
                    sf = ti['source_field']
                    src_val = field_value_arr[sf][i]
                    if src_val is None:
                        continue

                    if ti['is_categorical']:
                        if isinstance(src_val, str):
                            src_lower = src_val.strip().lower()
                        else:
                            src_lower = str(src_val).strip().lower()
                        if src_lower != ti['templ_value_norm_lower']:
                            continue
                        out_value = ti['templ_value_norm']
                        out_value_id = ti['templ_mvid_norm']
                    else:
                        out_value = src_val
                        out_value_id = ''

                    tun = ti['templ_unit_norm']
                    if tun == 'provided':
                        unit_val = field_unit_arr[sf][i]
                        if unit_val is None:
                            reporter.add_error(
                                f"For measurementType '{ti['output_meas_type']}', measurementUnit='provided' but unit is blank for occurrenceID '{occ_s}'."
                            )
                            raise ValueError("Blank per-row unit encountered under 'provided' policy")
                        out_unit = unit_val
                    elif tun == '':
                        out_unit = ''
                    else:
                        out_unit = tun

                    row_out = [
                        evt_s,
                        occ_s,
                        ti['verbatim_raw'],
                        ti['output_meas_type'],
                        out_value,
                        out_unit,
                        ti['templ_mtid_norm'],
                        out_value_id,
                        ti['templ_muid_norm'],
                        ti['templ_rem_norm'],
                    ]
                    writer.writerow(row_out)
                    row_count += 1
                    if len(preview_rows) < 15:
                        preview_rows.append(dict(zip(EMOF_OUTPUT_COLUMNS_IN_ORDER, row_out)))

        if row_count == 0:
            reporter.add_warning("No eMoF rows were generated based on the template and available data.")

        reporter.add_success(f"eMoF created successfully with {row_count} row(s) (occurrence-based)")
        reporter.add_text(f"Saved eMoF CSV to: {emof_csv_path}")
        reporter.add_text(f"File size: {row_count} rows (one per occurrence per measurement)")
        preview_df = pd.DataFrame(preview_rows, columns=EMOF_OUTPUT_COLUMNS_IN_ORDER)
        reporter.add_dataframe(preview_df, "eMoF Preview (first 15 rows)")

        return emof_csv_path

    except Exception as e:
        reporter.add_error(f"eMoF creation failed: {e}")
        raise


  