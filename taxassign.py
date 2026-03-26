import os
import sys
import argparse
import pandas as pd
import yaml

# Ensure we can import from src-v3 before importing CLI UI
THIS_DIR = os.path.dirname(os.path.abspath(__file__))
SRC_DIR = os.path.join(THIS_DIR, 'src-v3')
if SRC_DIR not in sys.path:
    sys.path.insert(0, SRC_DIR)

# CLI UI
from cli_output.cli_ui import console, print_header, silence_output
from taxonomic_assignment.taxa_assignment_manager import (
    format_taxa_assignment_info_dataframe,
    load_pr2_worms_dict_into_params,
    limit_info_df_preserving_selected,
    mark_selected_match_from_main_dataframe,
)
from taxonomic_assignment.taxa_assignment_info_export import write_taxa_assignment_info_xlsx

# Default config next to this script (no path typing required for normal use)
DEFAULT_CONFIG_PATH = os.path.join(THIS_DIR, 'config.yaml')

# Keys merged from config.yaml when using --use-config (taxonomy / matcher only; no data paths)
TAXONOMY_CONFIG_KEYS = (
    'taxonomic_api_source',
    'assays_to_skip_species_match',
    'gbif_match_limit',
    'worms_n_proc',
    'gbif_n_proc',
    'use_local_reference_database',
    'local_reference_database_path',
    'worms_min_ranks_matched',
    'worms_max_walkup_steps',
    'worms_return_all_matches',
    'worms_return_higher_classification',
    'worms_consider_unaccepted_for_selection',
)


# ------------------------------
# Users: EDIT THESE PARAMETERS (optional)
# CLI flags override these defaults. With --use-config, CLI flags also override the
# matching keys from config.yaml (see --help epilog).
DEFAULTS = {
    'INPUT_PATH': 'raw-v3/taxassign_example_input.tsv',
    'API': 'GBIF',                  # 'GBIF' or 'WoRMS'
    'MATCH_LIMIT': 3,               # max matches per verbatimIdentification in output
    'OUTPUT_DIR': 'processed-v3',   # used if OUTPUT_PATH is empty
    'OUTPUT_PATH': '',              # explicit path for taxa_assignment_INFO .xlsx; if empty, uses OUTPUT_DIR
    'N_PROC': None,               # None = with --use-config use config worms_n_proc/gbif_n_proc; else 0
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


def _load_taxonomy_params_from_config(config_path: str) -> dict:
    if not os.path.isfile(config_path):
        raise FileNotFoundError(f"Config not found: {config_path}")
    with open(config_path, 'r', encoding='utf-8') as f:
        raw = yaml.safe_load(f) or {}
    return {k: raw[k] for k in TAXONOMY_CONFIG_KEYS if k in raw}


def _resolve_n_proc(api: str, n_proc_cli: int | None, params: dict, use_config: bool) -> int:
    if n_proc_cli is not None:
        return int(n_proc_cli)
    if use_config:
        if api == 'WoRMS':
            return int(params.get('worms_n_proc', 0) or 0)
        return int(params.get('gbif_n_proc', 0) or 0)
    return 0


def _read_verbatim_input(path: str) -> pd.DataFrame:
    if not os.path.exists(path):
        raise FileNotFoundError(f"Input file not found: {path}")

    try:
        df = pd.read_csv(path, sep='\t', dtype=str)
    except Exception:
        df = pd.read_csv(path, dtype=str)

    if df.shape[1] == 1:
        col = df.columns[0]
        if str(col).strip().lower() != 'verbatimidentification':
            df = df.rename(columns={col: 'verbatimIdentification'})
    elif 'verbatimIdentification' not in df.columns:
        for col in df.columns:
            if str(col).strip().lower() == 'verbatimidentification':
                df = df.rename(columns={col: 'verbatimIdentification'})
                break

    if 'verbatimIdentification' not in df.columns:
        raise ValueError("Input must have a single column named 'verbatimIdentification' "
                         "or one column that will be treated as such.")

    df['verbatimIdentification'] = df['verbatimIdentification'].astype(str).fillna('').str.strip()
    df = df[df['verbatimIdentification'] != ''].copy()

    df = df.drop_duplicates(subset=['verbatimIdentification']).reset_index(drop=True)

    df['assay_name'] = 'taxassign'

    return df[['verbatimIdentification', 'assay_name']]


def _normalize_taxa_output_xlsx(path: str) -> str:
    """Output is always .xlsx. If the user passed .csv, swap to .xlsx."""
    p = path.strip()
    if p.lower().endswith('.csv'):
        return p[:-4] + '.xlsx'
    if not p.lower().endswith('.xlsx'):
        return p + '.xlsx'
    return p


def _write_taxa_assignment_xlsx(df: pd.DataFrame, out_path: str) -> str:
    out_path = _normalize_taxa_output_xlsx(out_path)
    write_taxa_assignment_info_xlsx(df, out_path)
    return out_path


def run_taxassign(
    input_path: str,
    api_source: str | None = None,
    match_limit: int | None = None,
    out_dir: str = 'processed-v3',
    out_path: str | None = None,
    n_proc: int | None = None,
    use_config: bool = False,
    config_path: str | None = None,
) -> str:
    cfg_path = config_path or DEFAULT_CONFIG_PATH

    params: dict = {
        'assays_to_skip_species_match': [],
        'output_dir': out_dir,
        'assay_rank_info': {'taxassign': {'max_depth': 99}},
    }

    if use_config:
        with console.status('[bold]Loading config.yaml...[/]', spinner='dots'):
            try:
                merged = _load_taxonomy_params_from_config(cfg_path)
            except FileNotFoundError as e:
                console.print(f'[bold red]{e}[/]')
                raise SystemExit(1) from e
        params.update(merged)

    if api_source is not None:
        api = _canonical_api_name(api_source)
    elif use_config:
        api = _canonical_api_name(params.get('taxonomic_api_source', 'WoRMS'))
    else:
        api = _canonical_api_name(DEFAULTS['API'])

    if match_limit is not None:
        gbif_lim = int(match_limit)
    elif use_config:
        gbif_lim = int(params.get('gbif_match_limit', DEFAULTS['MATCH_LIMIT']))
    else:
        gbif_lim = int(DEFAULTS['MATCH_LIMIT'])

    params['taxonomic_api_source'] = api
    params['gbif_match_limit'] = gbif_lim
    params['output_dir'] = out_dir
    params['assay_rank_info'] = {'taxassign': {'max_depth': 99}}

    n_proc_use = _resolve_n_proc(api, n_proc, params, use_config)

    load_pr2_worms_dict_into_params(
        params,
        warn_print=lambda m: console.print(f'[yellow]{m}[/]'),
    )

    console.print('[bold]Taxonomy / matcher settings[/]')
    if use_config:
        console.print(f'  Config file: [cyan]{cfg_path}[/]')
        console.print(
            '  [dim]API / row limit from config unless you pass[/] [cyan]-a[/] [dim]or[/] [cyan]-n[/][dim].[/]'
        )
    else:
        console.print('  [dim]Not using config.yaml; add[/] [cyan]--use-config[/] [dim]to load matcher options from it.[/]')
    console.print(f'  API (--api):                    {api}')
    console.print(f'  Rows per name (--limit):        {gbif_lim}')
    console.print(f'  Parallel jobs (--n-proc):       {n_proc_use}')
    if api == 'WoRMS':
        console.print(f'  WoRMS min ranks matched:        {params.get("worms_min_ranks_matched", 1)}')
        console.print(f'  WoRMS max walk-up steps:        {params.get("worms_max_walkup_steps")}')
        console.print(f'  WoRMS return all matches:       {params.get("worms_return_all_matches", False)}')
        console.print(f'  WoRMS higher classification:    {params.get("worms_return_higher_classification", False)}')
        console.print(f'  WoRMS unaccepted for selection:  {params.get("worms_consider_unaccepted_for_selection", False)}')
        console.print(f'  Local reference DB (PR2):       {params.get("use_local_reference_database", False)}')
        ldb = params.get('local_reference_database_path')
        if params.get('use_local_reference_database') and ldb:
            console.print(f'  Local DB path:                  {ldb}')
        console.print(f'  PR2 name lookups loaded:        {len(params.get("pr2_worms_dict") or {})}')

    df_in = _read_verbatim_input(input_path)

    console.print('[bold]Starting Taxonomic Assignment...[/]')
    with console.status('Running Taxonomic Assignment...', spinner='dots'):
        with silence_output():
            if api == 'GBIF':
                from taxonomic_assignment.GBIF_matching import get_gbif_match_for_dataframe
                results = get_gbif_match_for_dataframe(df_in.copy(), params, n_proc=n_proc_use)
            else:
                from taxonomic_assignment.WoRMS_v3_matching import get_worms_match_for_dataframe
                results = get_worms_match_for_dataframe(df_in.copy(), params, n_proc=n_proc_use)
            info_df = results.get('info_df', pd.DataFrame())
            main_df = results.get('main_df', pd.DataFrame())

    console.print('[green]Finished Taxonomic Assignment.[/]')

    if info_df is None or info_df.empty:
        info_df = pd.DataFrame({'verbatimIdentification': []})
        main_df = pd.DataFrame()

    info_df = mark_selected_match_from_main_dataframe(info_df, main_df, api)
    try:
        info_df_limited = limit_info_df_preserving_selected(info_df, int(gbif_lim))
    except Exception:
        info_df_limited = info_df.copy()

    out_df = format_taxa_assignment_info_dataframe(info_df_limited, params)

    if out_path is None or str(out_path).strip() == '':
        out_path = os.path.join(out_dir, f'taxa_assignment_INFO_{api}.xlsx')

    out_path = _write_taxa_assignment_xlsx(out_df, out_path)

    console.print(f'[green]Wrote {len(out_df):,} rows[/] to {out_path}')
    return out_path


class _TaxassignHelpFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
    pass


def main():
    print_header()
    parser = argparse.ArgumentParser(
        prog='python taxassign.py',
        description=(
            'Match verbatimIdentification strings to WoRMS or GBIF. '
            'Set your taxonomic assignment parameters for WoRMS/GBIF in config.yaml and add --use-config to use them. Use the flags below only when you want to override the parameters set in the config.yaml file.'
        ),
        formatter_class=_TaxassignHelpFormatter,
        epilog=(
            '-------------------------------------------------------------------------------\n'
            '  CONFIG  (--use-config)  vs  COMMAND LINE\n'
            '-------------------------------------------------------------------------------\n'
            '  With --use-config, matcher tuning is read from config.yaml (next to taxassign.py).\n'
            '  Input and output paths are NEVER taken from config:\n'
            '    -i  --input     TSV/CSV of verbatimIdentification\n'
            '    -o  --output    Full path to one .xlsx file (not a folder). .csv is rewritten to .xlsx\n'
            '        --outdir    Folder only; used when you omit -o\n'
            '\n'
            '  These flags override config when you pass them (with --use-config):\n'
            '\n'
            '    -a  --api       GBIF or WoRMS (omit to use config taxonomic_api_source)\n'
            '    -n  --limit     Max candidate rows per name (omit to use config gbif_match_limit)\n'
            '        --n-proc    Parallel workers (omit to use config worms_n_proc / gbif_n_proc)\n'
            '\n'
            '\n'
            '  WoRMS-only options (walk-up, PR2 file, marine-only, etc.) have no CLI — set them in config.yaml.\n'
        ),
    )
    parser.add_argument(
        '-i', '--input', default=DEFAULTS['INPUT_PATH'],
        help=(
            f"Input TSV/CSV path (one column: verbatimIdentification). "
            f"Never read from config, always set here or via the default: {DEFAULTS['INPUT_PATH']}"
        ),
    )
    parser.add_argument(
        '-a', '--api', default=None, metavar='API',
        help=(
            'GBIF or WoRMS. With --use-config, omit this to use taxonomic_api_source from config; '
            f'without --use-config, omitting defaults to {DEFAULTS["API"]}.'
        ),
    )
    parser.add_argument(
        '-n', '--limit', type=int, default=None, metavar='N',
        help=(
            'Max candidate rows per name (best match always kept). With --use-config, omit to use gbif_match_limit from config; '
            f'without --use-config, omitting defaults to {DEFAULTS["MATCH_LIMIT"]}.'
        ),
    )
    parser.add_argument(
        '-o', '--output', default=None, metavar='FILE',
        help='Output .xlsx path (include a filename; not a folder). If you pass .csv it is saved as .xlsx. If omitted, uses --outdir and the default filename.',
    )
    parser.add_argument(
        '--outdir', default=DEFAULTS['OUTPUT_DIR'], metavar='DIR',
        help=(
            f"Folder for output when -o is omitted (writes taxa_assignment_INFO_<API>.xlsx). Default: {DEFAULTS['OUTPUT_DIR']}. Not from config."
        ),
    )
    parser.add_argument(
        '--n-proc', type=int, default=None,
        help='Number of parallel workers. With --use-config: if you omit this flag, config decides; if you pass a number, it overrides config. Without --use-config: defaults to 0.',
    )
    parser.add_argument(
        '--use-config', action='store_true',
        help=(
            f'Read WoRMS/GBIF matcher options from config.yaml next to this script '
            f'({DEFAULT_CONFIG_PATH}). '
            'Does not use config for input file path, pipeline data paths, or output folders. Config.yaml only provides taxonomic assignment options.'
        ),
    )

    args = parser.parse_args()

    if args.api is not None:
        a = str(args.api).strip().lower()
        if a not in ('gbif', 'worms'):
            parser.error('--api: expected GBIF or WoRMS')
        args.api = _canonical_api_name(args.api)

    run_taxassign(
        input_path=args.input,
        api_source=args.api,
        match_limit=args.limit,
        out_dir=args.outdir,
        out_path=args.output,
        n_proc=args.n_proc,
        use_config=args.use_config,
    )


if __name__ == '__main__':
    main()
