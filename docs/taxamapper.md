# taxamapper: taxonomic mapping only

If you only want to map taxonomy for a list of identifiers (without running the full edna2obis workflow), use `taxamapper.py`.

- Input: a TSV/CSV with a single column named `verbatimIdentification`. An example is provided at `data/example_raw/taxamapper_example_input.tsv`.
- Output: `taxa_mapping_INFO_<API>.xlsx` (same column layout as the full pipeline; sheet `taxa_mapping_INFO` only). **Important:** `selected_match` is `True` on the row that would be chosen as the official mapping (equivalent to what goes to occurrence in the full pipeline). You may get multiple candidate rows per name; the rest are for review. Candidate count is controlled by `-n` (or by `gbif_match_limit` from `config.yaml` when using `--use-config`).

Help:

```bash
python taxamapper.py --help
```

## Examples

```bash
# Example input, default API and limits from taxamapper.py, writes outputs/taxa_mapping_INFO_<API>.xlsx
python taxamapper.py

# Map taxonomies to WoRMS using an input file
python taxamapper.py -i data/example_raw/taxamapper_example_input.tsv -a WoRMS

# Up to 5 candidate rows per name, and map to GBIF backbone
python taxamapper.py -i data/example_raw/taxamapper_example_input.tsv -a GBIF -n 5

# Explicit output path
python taxamapper.py -i data/example_raw/taxamapper_example_input.tsv -a GBIF -o outputs/my_run.xlsx

# Matcher options from config.yaml
python taxamapper.py --use-config

# Parallelism
python taxamapper.py -i data/example_raw/taxamapper_example_input.tsv -a GBIF --n-proc 3
```

## Options

- `-i, --input`: TSV/CSV with column `verbatimIdentification` (default: `data/example_raw/taxamapper_example_input.tsv`)
- `-a, --api`: `GBIF` or `WoRMS` (omit with `--use-config` to use `taxonomic_api_source` from config)
- `-n, --limit`: max candidate rows per name (omit with `--use-config` to use `gbif_match_limit` from config)
- `-o, --output`: output `.xlsx` path (not a folder). If omitted, `{outdir}/taxa_mapping_INFO_<API>.xlsx`
- `--outdir`: directory when `-o` omitted (default: `outputs`)
- `--use-config`: load matcher settings from `config.yaml` beside `taxamapper.py`
- `--n-proc`: parallel workers (optional)

## In-file defaults

- You can also set defaults directly at the top of `taxamapper.py` in the `DEFAULTS` block (input path, API, match limit, outdir/output, n-proc). Any CLI flags you pass will override those defaults.

See [Taxonomic mapping](taxonomic-mapping.md) for how matches are chosen and what each output column means.
