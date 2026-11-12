# Changelog

All notable changes to edna2obis are documented here.
The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and versions follow [Semantic Versioning](https://semver.org/).
See [docs/RELEASING.md](docs/RELEASING.md) for how releases are made.

## [Unreleased]

### Changed
- **Breaking:** renamed directories to drop the version suffix: `src-v3/` → `src/`, `raw-v3/` → `data/example_raw/`, `processed-v3/` → `outputs/`. Update paths in your own config files.
- README restructured around a Quick Start; detailed reference moved to `docs/` (`input-data.md`, `outputs.md`, `taxonomic-mapping.md`, `taxamapper.md`).
- The version is now defined in `src/_version.py` and shown by `python main.py --version`, `python taxamapper.py --version`, the CLI banner and the HTML report.

### Fixed
- README: corrected example config names, eMoF template format (`.tsv`) and required columns, script path, and output file naming.
- `main.py`: eMoF template fallback now points at the bundled `data/example_raw/eMoF_Fields_edna2obis.tsv`.
- Config templates: removed `excel_file` paths that pointed at files not in the repository.

## [3.2.0] - 2026-10-26

### Added
- Event Core output (`dwc_core_type`), including short_name splitting and mixed event/occurrence eMoF linkage.
- Map preview of sample locations in the HTML report.
- HTML report states whether the archive is an Occurrence Core or Event Core.
- Video tutorial link in the README.

### Changed
- Renamed "taxonomic assignment" to "taxonomic mapping" across the repository; `taxa_mapping_INFO_<API>.xlsx` output.
- `meta.xml` is parameterized instead of hard-coding the occurrenceID linkage.
- Cleaner taxamapper CLI output and help text.

### Removed
- `identificationRemarks` from the DNA Derived Extension.

## [3.1.0] - 2026-07-17

### Added
- Bundled example datasets (small and large FAIRe NOAA, generic FAIRe).
- Configurable header rows so input files no longer need to come from FAIReSheets/Tourmaline.
- `venv` + `pip` setup for non-conda users.
- `contacts.yaml` and fuller EML generation, including per-short_name EML.
- WoRMS lineage-consistency heuristic, optional non-marine matches, and unaccepted-match review rows.
- Optional splitting of outputs by `short_name`.

### Changed
- occurrenceIDs include `project_id` and `analysis_run_name` so they are globally unique.
- Much faster eMoF creation; WoRMS worker limit raised from 3 to 10.

### Fixed
- Duplicate warnings/errors in the HTML report, EML scope attribute rejected by IPT, and `selected_match` selection bugs.

## [3.0.0-alpha] - 2025-07-08

First script-based (non-notebook) release with configuration files. Tags `v3.0.1-alpha` and `v3.0.2-alpha` contain README updates only.

## [1.0.0] - 2025-02-18

Notebook-based workflow. See the `legacy-v1` branch.

[Unreleased]: https://github.com/aomlomics/edna2obis/compare/v3.2.0...HEAD
[3.2.0]: https://github.com/aomlomics/edna2obis/compare/v3.1.0...v3.2.0
[3.1.0]: https://github.com/aomlomics/edna2obis/compare/v3.0.2-alpha...v3.1.0
[3.0.0-alpha]: https://github.com/aomlomics/edna2obis/compare/v1.0.0...v3.0.0-alpha
[1.0.0]: https://github.com/aomlomics/edna2obis/releases/tag/v1.0.0
