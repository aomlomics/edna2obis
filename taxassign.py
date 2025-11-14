import os
import sys
import argparse
import pandas as pd

# Making sure we can import taxonomic assignment scripts 
THIS_DIR = os.path.dirname(os.path.abspath(__file__))
SRC_DIR = os.path.join(THIS_DIR, 'src-v3')
if SRC_DIR not in sys.path:
    sys.path.insert(0, SRC_DIR)


# ------------------------------
# Users: EDIT THESE PARAMETERS (optional)
# You can edit these values to set parameters if you do NOT provide CLI flags.
# Command-line flags always take precedence over these defaults.
DEFAULTS = {
    'INPUT_PATH': 'raw-v3/taxassign_example_input.tsv',
    'API': 'GBIF',                  # 'GBIF' or 'WoRMS'
    'MATCH_LIMIT': 3,               # max matches per verbatimIdentification in output
    'OUTPUT_DIR': 'processed-v3',   # used if OUTPUT_PATH is empty
    'OUTPUT_PATH': '',              # explicit output path for the result TSV; if empty, uses OUTPUT_DIR
    'N_PROC': 0,                    # parallel processes; 0 lets matcher decide
}
# ------------------------------


def _canonical_api_name(name: str) -> str:
    if not isinstance(name, str):
        return 'GBIF'
    name_clean = name.strip().lower()
    if name_clean == 'worms':
        return 'WoRMS'
    if name_clean == 'gbif':
        return 'GBIF'
    return 'GBIF'


def _read_verbatim_input(path: str) -> pd.DataFrame:
    if not os.path.exists(path):
        raise FileNotFoundError(f"Input file not found: {path}")

    # Try TSV first, then CSV as fallback
    try:
        df = pd.read_csv(path, sep='\t', dtype=str)
    except Exception:
        df = pd.read_csv(path, dtype=str)

    if df.shape[1] == 1:
        # Single column – make sure it's named verbatimIdentification
        col = df.columns[0]
        if str(col).strip().lower() != 'verbatimidentification':
            df = df.rename(columns={col: 'verbatimIdentification'})
    elif 'verbatimIdentification' not in df.columns:
        # Try case-insensitive match
        for col in df.columns:
            if str(col).strip().lower() == 'verbatimidentification':
                df = df.rename(columns={col: 'verbatimIdentification'})
                break

    if 'verbatimIdentification' not in df.columns:
        raise ValueError("Input must have a single column named 'verbatimIdentification' "
                         "or one column that will be treated as such.")

    # Clean and drop empties
    df['verbatimIdentification'] = df['verbatimIdentification'].astype(str).fillna('').str.strip()
    df = df[df['verbatimIdentification'] != ''].copy()

    # Deduplicate
    df = df.drop_duplicates(subset=['verbatimIdentification']).reset_index(drop=True)

    # Provide a dummy assay_name used by existing matching functions
    df['assay_name'] = 'taxassign'

    return df[['verbatimIdentification', 'assay_name']]


def _write_tsv(df: pd.DataFrame, out_path: str) -> None:
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    df.to_csv(out_path, sep='\t', index=False)


def run_taxassign(input_path: str,
                  api_source: str = 'GBIF',
                  match_limit: int = 3,
                  out_dir: str = 'processed-v3',
                  out_path: str | None = None,
                  n_proc: int = 0) -> str:
    api = _canonical_api_name(api_source)

    df_in = _read_verbatim_input(input_path)

    # Build params for existing matchers
    params = {
        'taxonomic_api_source': api,
        'assays_to_skip_species_match': [],
        'output_dir': out_dir,
        # GBIF uses this to cap candidate fetches; WoRMS returns all, we will post-cap
        'gbif_match_limit': int(match_limit),
        # Provide a broad depth to avoid species-trimming behavior
        'assay_rank_info': {'taxassign': {'max_depth': 99}},
    }

    if api == 'GBIF':
        from taxonomic_assignment.GBIF_matching import get_gbif_match_for_dataframe
        results = get_gbif_match_for_dataframe(df_in.copy(), params, n_proc=n_proc)
        info_df = results.get('info_df', pd.DataFrame())
    else:
        from taxonomic_assignment.WoRMS_v3_matching import get_worms_match_for_dataframe
        results = get_worms_match_for_dataframe(df_in.copy(), params, n_proc=n_proc)
        info_df = results.get('info_df', pd.DataFrame())

    if info_df is None or info_df.empty:
        info_df = pd.DataFrame({'verbatimIdentification': []})

    # Enforce top-N per verbatimIdentification for output file
    # (GBIF limits at fetch time; WoRMS may return many, so we cap here as well)
    try:
        info_df_limited = (
            info_df
            .groupby('verbatimIdentification', group_keys=False)
            .head(int(match_limit))
            .reset_index(drop=True)
        )
    except Exception:
        info_df_limited = info_df.copy()

    # Determine output path
    if out_path is None or str(out_path).strip() == '':
        out_path = os.path.join(out_dir, f"taxassign_INFO_{api}.tsv")

    _write_tsv(info_df_limited, out_path)

    print(f"Wrote {len(info_df_limited):,} rows to {out_path}")
    return out_path


def main():
    parser = argparse.ArgumentParser(
        description="Run taxonomy assignment on a list of verbatimIdentification strings and emit an INFO-style TSV."
    )
    parser.add_argument(
        '-i', '--input', default=DEFAULTS['INPUT_PATH'],
        help=f"Path to a TSV/CSV with one column named 'verbatimIdentification' (default: {DEFAULTS['INPUT_PATH']})"
    )
    parser.add_argument(
        '-a', '--api', default=DEFAULTS['API'], choices=['GBIF', 'gbif', 'WoRMS', 'worms'],
        help=f"Taxonomic API to use (default: {DEFAULTS['API']})"
    )
    parser.add_argument(
        '-n', '--limit', type=int, default=DEFAULTS['MATCH_LIMIT'],
        help=f"Maximum number of matches to include per verbatimIdentification in the output (default: {DEFAULTS['MATCH_LIMIT']})"
    )
    parser.add_argument(
        '-o', '--output', default=DEFAULTS['OUTPUT_PATH'],
        help=f"Optional explicit output path for the result TSV. If omitted, writes to {{outdir}}/taxassign_INFO_<API>.tsv (default outdir: {DEFAULTS['OUTPUT_DIR']})"
    )
    parser.add_argument(
        '--outdir', default=DEFAULTS['OUTPUT_DIR'],
        help=f"Output directory to use when --output is not specified (default: {DEFAULTS['OUTPUT_DIR']})"
    )
    parser.add_argument(
        '--n-proc', type=int, default=DEFAULTS['N_PROC'],
        help=f"Parallel processes for matching (0 lets the matcher decide a sensible default; default: {DEFAULTS['N_PROC']})"
    )

    args = parser.parse_args()

    run_taxassign(
        input_path=args.input,
        api_source=args.api,
        match_limit=args.limit,
        out_dir=args.outdir,
        out_path=args.output,
        n_proc=args.n_proc
    )


if __name__ == '__main__':
    main()


