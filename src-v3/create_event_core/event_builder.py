"""
Event Core builder for edna2obis.

Builds a 2-tier Event Core from the finished Occurrence Core:
  - Sample event:  eventID = samp_name,  parentEventID = "",         eventRemarks = "sample"
  - Library event: eventID = lib_id,      parentEventID = samp_name,  eventRemarks = "library"

We derive events straight from the Occurrence Core file, reusing the exact eventID
(lib_id) and parentEventID (samp_name) values that ended up in the occurrence rows.
That way every eventID the extensions point at is guaranteed to exist here, so the
archive has no orphan records.
"""

import os
import numpy as np
import pandas as pd
import yaml


def _load_event_core_mapping(params):
    """Read the event_core column list (and per-column tier) from data_mapper.yaml."""
    with open('data_mapper.yaml', 'r', encoding='utf-8') as f:
        mapper = yaml.safe_load(f) or {}
    prefix = 'generic_' if str(params.get('metadata_format', 'NOAA')).upper() == 'GENERIC' else ''
    cfg = mapper.get(f'{prefix}event_core') or mapper.get('event_core')
    if not cfg:
        raise ValueError("No 'event_core' section found in data_mapper.yaml.")
    return cfg


def create_event_core(params, reporter=None, occurrence_filename=None):
    """Write event_core.csv next to the Occurrence Core and return its filename."""
    output_dir = params.get('output_dir', 'processed-v3/')
    api = str(params.get('taxonomic_api_source', 'worms')).lower()
    if occurrence_filename is None:
        occurrence_filename = f"occurrence_core_{api}.csv"
    occ_path = os.path.join(output_dir, occurrence_filename)
    if not os.path.exists(occ_path):
        raise FileNotFoundError(f"Occurrence Core not found for Event Core build: {occ_path}")

    event_map = _load_event_core_mapping(params)
    columns = list(event_map.keys())
    tiers = {col: str(event_map.get(col, {}).get('tier', 'sample')).strip().lower() for col in columns}

    # Read as text so Event Core values match the Occurrence Core file character-for-character.
    occ = pd.read_csv(occ_path, dtype=str, keep_default_na=False)
    for required in ('eventID', 'parentEventID'):
        if required not in occ.columns:
            raise ValueError(f"Occurrence Core is missing '{required}'; cannot build Event Core.")

    # Treat blanks as missing so groupby.first() grabs the first real value in each group.
    work = occ.replace('', np.nan)

    if reporter:
        missing_event = int(work['eventID'].isna().sum())
        missing_parent = int(work['parentEventID'].isna().sum())
        if missing_event or missing_parent:
            reporter.add_warning(
                f"Event Core: {missing_event} occurrence row(s) had an empty eventID and "
                f"{missing_parent} had an empty parentEventID; those are skipped for the affected tier."
            )

    # eventID/parentEventID/eventRemarks are set explicitly, everything else is pulled from the occurrence row.
    id_like = {'eventID', 'parentEventID', 'eventRemarks'}
    sample_value_cols = [c for c in columns if c not in id_like and tiers[c] in ('sample', 'both') and c in occ.columns]
    library_value_cols = [c for c in columns if c not in id_like and tiers[c] in ('library', 'both') and c in occ.columns]

    # Sample events: one row per samp_name (the occurrence's parentEventID).
    sample_df = work.groupby('parentEventID', sort=True)[sample_value_cols].first().reset_index()
    sample_df = sample_df.rename(columns={'parentEventID': 'eventID'})
    sample_df['parentEventID'] = ''
    sample_df['eventRemarks'] = 'sample'

    # Library events: one row per lib_id (the occurrence's eventID); parent is its samp_name.
    lib_pull = library_value_cols.copy()
    if 'parentEventID' not in lib_pull:
        lib_pull.append('parentEventID')
    library_df = work.groupby('eventID', sort=True)[lib_pull].first().reset_index()
    library_df['eventRemarks'] = 'library'

    # Line both tiers up on the same column order; anything a tier doesn't carry stays blank.
    sample_df = sample_df.reindex(columns=columns)
    library_df = library_df.reindex(columns=columns)
    event_core = pd.concat([sample_df, library_df], ignore_index=True, sort=False)
    event_core = event_core.fillna('')

    # eventID is the core id, so it has to be unique across both tiers.
    dup_ids = event_core.loc[event_core['eventID'].duplicated(keep=False), 'eventID'].unique()
    if len(dup_ids) > 0:
        raise ValueError(
            f"Event Core has {len(dup_ids)} duplicated eventID value(s) (a samp_name and a lib_id may collide). "
            f"Examples: {', '.join(map(str, dup_ids[:10]))}"
        )

    out_path = os.path.join(output_dir, 'event_core.csv')
    event_core.to_csv(out_path, index=False, na_rep='')
    if reporter:
        n_sample = int((event_core['eventRemarks'] == 'sample').sum())
        n_library = int((event_core['eventRemarks'] == 'library').sum())
        reporter.add_text(
            f"Created Event Core: {n_sample} sample event(s) + {n_library} library event(s) -> {out_path}"
        )
    return 'event_core.csv'
