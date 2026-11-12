# Outputs

Each run writes its files to `<output_dir>/<edna2obis_run_name>/` (default `outputs/<edna2obis_run_name>/`).

| File | Contents |
|---|---|
| `edna2obis_report_<API>_<edna2obis_run_name>.html` | Run report: settings, warnings, errors and previews of each output |
| `config_<edna2obis_run_name>.yaml` | Copy of the config used for the run |
| `occurrence_core_worms.csv` / `occurrence_core_gbif.csv` | Occurrence Core with mapped taxonomy |
| `taxa_mapping_INFO_<API>.xlsx` | Every candidate taxonomy match; **`selected_match = TRUE`** marks the one used. See [Taxonomic mapping](taxonomic-mapping.md#taxa_mapping_info-columns) |
| `dna_derived_extension.csv` | DNA Derived Data extension |
| `eMoF.csv` | Extended Measurement or Fact extension (if `emof_enabled: true`) |
| `meta.xml` | Darwin Core Archive file mapping |
| `eml.xml` | Dataset metadata in EML (if `eml_enabled: true`) |

With `split_output_by_short_name: true`, the occurrence core, DNA derived extension, eMoF and EML are also split into one subfolder per `short_name` (cruise/expedition).

The example rows below come from the bundled small FAIRe NOAA dataset (the default `config.yaml`, WoRMS). Only a few key columns are shown; the real files contain many more.

## Occurrence Core

| occurrenceID | eventID | verbatimIdentification | scientificName | scientificNameID | taxonRank | organismQuantity | organismQuantityType | eventDate | decimalLatitude | decimalLongitude |
|---|---|---|---|---|---|---|---|---|---|---|
| noaa-sefsc-gu1901:rosel-ctd-mifish-addMFtrawlJfish-241001:GU190706-CTD11-220_MiFish_S30_occ_18109634cc2f8e156e5402bf13cf4502 | GU190706-CTD11-220_MiFish_S30 | Eukaryota;Chordata;Actinopteri;Beloniformes;Exocoetidae;Cheilopogon | Cheilopogon | urn:lsid:marinespecies.org:taxname:125691 | Genus | 4 | DNA sequence reads | 2019-07-06 | 28.662 | -85.793 |
| noaa-sefsc-gu1901:rosel-ctd-mifish-addMFtrawlJfish-241001:GU190706-CTD11-220_MiFish_S30_occ_183bc18f3e5eac45c6dd248fb86d64bf | GU190706-CTD11-220_MiFish_S30 | Eukaryota;Chordata;Actinopteri;Tetraodontiformes;Tetraodontidae;Lagocephalus;Lagocephalus laevigatus | Lagocephalus laevigatus | urn:lsid:marinespecies.org:taxname:158933 | Species | 5317 | DNA sequence reads | 2019-07-06 | 28.662 | -85.793 |

`occurrenceID` is `<project_id>:<analysis_run_name>:<eventID>_occ_<featureid>`, which keeps it unique across OBIS and GBIF.

## DNA Derived Extension

One row per occurrence, linked by `occurrenceID`.

| occurrenceID | target_gene | pcr_primer_name_forward | pcr_primer_name_reverse | seq_meth | otu_db | DNA_sequence |
|---|---|---|---|---|---|---|
| noaa-sefsc-gu1901:…_occ_18109634cc2f8e156e5402bf13cf4502 | 12S rRNA (SSU mitochondria) | MiFish-U-F | MiFish-U-R2 | Illumina MiSeq [OBI_0002003] | custom;rc_12S_efc_fish-addSI_MF_trawl_J-mifish-taxa-derepU | CACCGCGGTTATACGAGAGG… |
| noaa-sefsc-gu1901:…_occ_183bc18f3e5eac45c6dd248fb86d64bf | 12S rRNA (SSU mitochondria) | MiFish-U-F | MiFish-U-R2 | Illumina MiSeq [OBI_0002003] | custom;rc_12S_efc_fish-addSI_MF_trawl_J-mifish-taxa-derepU | CACCGCGGTTATACGATGAA… |

Which metadata fields end up in the Occurrence Core and DNA Derived Extension is controlled by [`data_mapper.yaml`](../data_mapper.yaml).

## eMoF (extended Measurement or Fact)

The eMoF file captures occurrence-level measurements linked to each occurrenceID that made it into the final occurrence file. This creates a larger file but is what GBIF requires for data submission.

- **What it contains**: One row per occurrence per configured measurement.
- **Where values come from**: `sampleMetadata` first, otherwise `experimentRunMetadata`. If a field in the template isn't in your data, it is skipped.
- **Units**:
  - If the template specifies a literal unit (e.g., `m`, `°C`), that unit is used for every row of that measurementType.
  - If the template says `provided`, a column named `<measurementType>_unit` must be present in the source sheet and non-blank for all emitted rows.
  - If the template leaves the unit blank, the output unit is blank (no auto-fallback).

### Configuring measurements

Measurements are listed in the template [`data/example_raw/eMoF_Fields_edna2obis.tsv`](../data/example_raw/eMoF_Fields_edna2obis.tsv) (tab-separated). Point `emof_template_path` in `config.yaml` at your own copy if you don't want to edit the bundled one.

- Required columns: `measurementType`, `measurementValue`, `measurementUnit`
- Optional columns: `verbatimMeasurementType`, `measurementTypeID`, `measurementValueID`, `measurementUnitID`, `measurementRemarks`

You can add environmental measurements (temperature, salinity, pH, …), sampling parameters (gear, depth, volume filtered, …), laboratory measurements (DNA concentration, …) or any other data linked to your events or occurrences. Use controlled vocabularies for the IDs; see the [OBIS eMoF formatting guide](https://manual.obis.org/format_emof.html).

### Example

| eventID | occurrenceID | verbatimMeasurementType | measurementType | measurementValue | measurementUnit | measurementTypeID | measurementUnitID | measurementRemarks |
|---|---|---|---|---|---|---|---|---|
| GU190706-CTD11-220_MiFish_S30 | noaa-sefsc-gu1901:…_occ_18109634cc2f8e156e5402bf13cf4502 | chlorophyll | Chlorophyll fluorescence | 0.038 | micrograms per litre | http://vocab.nerc.ac.uk/collection/P01/current/CPHLMOD2/ | http://vocab.nerc.ac.uk/collection/P06/current/UGPL/ | Collected by WetLabs WetStar Fluorometer |
| GU190706-CTD11-220_MiFish_S30 | noaa-sefsc-gu1901:…_occ_18109634cc2f8e156e5402bf13cf4502 | ctd_bottle_number | CTD bottle number | CTD11 | | | | |
