# edna2obis: eDNA Data Publishing Workflow

## Introduction

DNA derived data are increasingly being used to document taxon 
occurrences. To ensure these data are useful to the broadest possible 
community, [GBIF](https://www.gbif.org/) published a guide entitled "[Publishing DNA-derived 
data through biodiversity data platforms](https://docs.gbif.org/publishing-dna-derived-data/en/)." 
This guide is supported by the [DNA derived data extension](https://tools.gbif.org/dwca-validator/extension.do?id=http://rs.gbif.org/terms/1.0/DNADerivedData) 
for [Darwin Core](https://dwc.tdwg.org/), which incorporates [MIxS](https://gensc.org/mixs/) 
terms into the Darwin Core standard. 

This use case draws on both the guide and the extension to develop a workflow 
for incorporating a DNA derived data extension file into a Darwin Core
archive. 

The latest version of edna2obis (version 3) builds upon the original edna2obis, introducing new features:
- Moved from a Jupyter Notebook to script architecture (runs in one command)
- Specify parameters in the `config.yaml`, rather than in the code
- Takes the new FAIRe NOAA eDNA data format as input, which is compatible for upload to the [Ocean DNA Explorer](https://www.oceandnaexplorer.org/)
- Users can choose to perform their taxonomic assignment via [WoRMS]() or [GBIF]() APIs 
- Improved taxonomic assignment accuracy and performance, with new caching methods
- Users can specify which assays to NOT include species rank for taxonomic assignment (for example, Bacterial taxonomies often have the HOST organism as the species)
- A new output file is created, `taxa_assignment_INFO.csv`, which gives information on how the taxonomies were assigned
- Generates an HTML output report to document your edna2obis run


### Project abstract:

Seawater was collected on board the NOAA ship Ronald H. Brown as part of the fourth Gulf of Mexico Ecosystems and Carbon Cycle (GOMECC-4) cruise from September 13 to October 21, 2021. Sampling for GOMECC-4 occurred along 16 coastal-offshore transects across the entire Gulf of Mexico and an additional line at 27N latitude in the Atlantic Ocean. We also collected eDNA samples near Padre Island National Seashore (U.S. National Parks Service), a barrier island located off the coast of south Texas. Vertical CTD sampling was employed at each site to measure discrete chemical, physical, and biological properties. Water sampling for DNA filtration was conducted at 54 sites and three depths per site (surface, deep chlorophyll maximum, and near bottom) to capture horizontal and vertical gradients of bacterial, protistan, and metazoan diversity across the Gulf. The resulting ASVs, their assigned taxonomy, and the metadata associated with theircollection are the input data for the OBIS conversion scripts presented here.

### Published data
- [GBIF](https://www.gbif.org/dataset/9012def0-bd87-48a0-ac9e-e0e78dd37689)
- [OBIS](https://obis.org/dataset/210efc7c-4762-47ee-b4b5-22a0f436ef44)

## NOAA Omics FAIR eDNA-based metadata template
This code was developed to convert a custom FAIRe NOAA Google Sheet metadata template developed by NOAA Omics at AOML, and based off the [FAIRe eDNA data standard](https://fair-edna.github.io/index.html). To use the sheet for your own data, run [FAIRe2NODE](https://github.com/aomlomics/FAIReSheets/tree/FAIRe2NODE), a suite of Python scripts which will generate the FAIRe NOAA templates in Google Sheets. Here is a filled-in example:  

[FAIRe_NOAA_noaa-aoml-gomecc4_SHARING](https://docs.google.com/spreadsheets/d/1mkjfUQW3gTn3ezhMQmFDQn4EBoQ2Xv4SZeSd9sqagoU/edit?gid=0#gid=0)

## Current Repo Structure (v3.0)
```
edna2obis-3.0/
├── README.md                           
├── LICENSE                             
├── environment.yml                     # Conda environment definition
├── config.yaml                         # Configuration file
├── main.py                            # Main script
├── src-v3/                            
│   ├── html_reporter.py               # HTML report generator
│   ├── create_asv_seq_taxa_obis.sh    # ASV sequence creation script
│   ├── edna2obis_conversion_code_v3.ipynb # Original Jupyter notebook
│   ├── edna2obis_conversion_code.md   
│   ├── create_occurrence_core/        
│   │   └── occurrence_builder.py      # Main occurrence builder
│   ├── create_dna_derived_extension/   
│   │   └── extension_builder.py       # Main extension builder
│   └── taxonomic_assignment/           
│       ├── taxa_assignment_manager.py 
│       ├── WoRMS_v3_matching.py       # WoRMS API matching
│       └── GBIF_matching.py           # GBIF API matching
├── raw-v3/                           
│   ├── FAIRe_NOAA_*.xlsx              # FAIRe NOAA metadata templates
│   ├── asvTaxaFeatures_*.tsv          # ASV taxonomy features
│   ├── table_*.tsv                    # ASV abundance tables
│   └── pr2_version_*.xlsx             # PR2 taxonomy database
├── processed-v3/                      
│   ├── occurrence.csv                 # Intermediate file
│   ├── occurrence_*_matched.csv       # Occurrence Core with assigned taxonomy
│   ├── dna_derived_extension.csv      # DNA Derived Extension
│   └── taxa_assignment_INFO.csv       # Taxonomy assignment details
```

## 🚀 Setup and Installation

### Prerequisites

- [Conda](https://docs.conda.io/en/latest/miniconda.html) or [Anaconda](https://www.anaconda.com/products/distribution) installed
- [Git](https://git-scm.com/downloads) installed
- At least 8GB RAM recommended
- Internet connection required (for API calls to WoRMS/GBIF)

### Quick Start

#### 1. Clone the Repository

```bash
git clone https://github.com/baydenwillms/edna2obis-3.0.git
cd edna2obis-3.0
```

#### 2. Create Conda Environment

```bash
# Create the environment from the environment.yml file
conda env create -f environment.yml

# Activate the environment
conda activate edna2obis
```

#### 3. Configure Your Data

Edit the `config.yaml` file with your data filepaths and other parameters


Key settings to update:
- `excel_file`: Path to your FAIRe NOAA Excel file (data template)
- `datafiles`: Paths to your ASV taxonomy and occurrence files  
- `taxonomic_api_source`: Choose "WoRMS" or "GBIF"
- `output_dir`: Where to save results (default: "processed-v3/")

#### 4. Run the Pipeline

```bash
python main.py
```

The pipeline will:
- Load and clean your metadata (according to OBIS/GBIF)
- Align data to Darwin Core data standard
- Generate an Occurrence Core
- Perform taxonomic assignment via WoRMS or GBIF APIs
- Generate a DNA Derived Extension
- Create an HTML report with results from your run

### Output Files

The pipeline generates several files in your `output_dir`:

- `occurrence.csv` - Initial occurrence data (intermediate file)
- `occurrence_worms_matched.csv` / `occurrence_gbif_matched.csv` - Final Occurrence Core with taxonomy
- `taxa_assignment_INFO.csv` - Summary of HOW taxonomies were assigned
- `dna_derived_extension.csv` - DNA-Derived data extension
- `edna2obis_report.html` - HTML output report

### Troubleshooting

#### Common Issues

1. **Environment creation fails**
   ```bash
   # Try updating conda first
   conda update conda
   conda env create -f environment.yml
   ```

2. **API timeout errors**
   - Check internet connection
   - Reduce worms_n_proc or gbif_n_proc in `config.yaml`
   - The pipeline has built-in retry logic

3. **Missing data files**
   - Verify all file paths in config.yaml are correct
   - Use absolute paths if relative paths don't work

#### Getting Help

- Check the HTML report for detailed error messages
- Review the terminal output for specific error details
- **Ensure your input data follows the FAIRe NOAA format**

#### Recommended System Requirements
- Processing: 8GB+ RAM, 4+ CPU cores
- Storage: ~1GB free space for large datasets
- Network: Stable internet for API calls

## Disclaimer  
This repository is a scientific product and is not official communication of the National Oceanic and Atmospheric Administration, or the United States Department of Commerce. All NOAA GitHub project code is provided on an 'as is' basis and the user assumes responsibility for its use. Any claims against the Department of Commerce or Department of Commerce bureaus stemming from the use of this GitHub project will be governed by all applicable Federal law. Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation or favoring by the Department of Commerce. The Department of Commerce seal and logo, or the seal and logo of a DOC bureau, shall not be used in any manner to imply endorsement of any commercial product or activity by DOC or the United States Government.
