# edna2obis workflow

## 🚀 Setup and Installation

### Prerequisites

- [Conda](https://docs.conda.io/en/latest/miniconda.html) or [Anaconda](https://www.anaconda.com/products/distribution) installed
- [Git](https://git-scm.com/downloads) installed
- At least 8GB RAM recommended
- Internet connection (for API calls to WoRMS/GBIF)

### Quick Start

#### 1. Clone the Repository

```bash
git clone <repository-url>
cd edna2obis-3/edna2obis-3.0
```

#### 2. Create Conda Environment

```bash
# Create the environment from the environment.yml file
conda env create -f environment.yml

# Activate the environment
conda activate edna2obis
```

#### 3. Verify Installation

```bash
# Test that Python and key packages are available
python -c "import pandas, numpy, pyworms, pygbif, yaml; print('✅ All dependencies installed successfully!')"
```

#### 4. Configure Your Data

Edit the `config.yaml` file with your data paths and settings:

```bash
# Open config.yaml in your preferred editor
nano config.yaml  # or code config.yaml, vim config.yaml, etc.
```

Key settings to update:
- `excel_file`: Path to your FAIRe NOAA Excel file
- `datafiles`: Paths to your ASV taxonomy and occurrence files  
- `taxonomic_api_source`: Choose "WoRMS" or "GBIF"
- `output_dir`: Where to save results (default: "processed-v3/")

#### 5. Run the Pipeline

```bash
python main.py
```

The pipeline will:
- Load and clean your metadata
- Create occurrence core data
- Perform taxonomic assignment via WoRMS or GBIF APIs
- Generate Darwin Core compliant output files
- Create an HTML report with results

### Output Files

The pipeline generates several files in your `output_dir`:

- `occurrence.csv` - Initial occurrence data (intermediate file)
- `occurrence_worms_matched.csv` / `occurrence_gbif_matched.csv` - Final occurrence core with taxonomy
- `taxa_assignment_INFO.csv` - Summary of unique taxonomy assignments
- `dna_derived_extension.csv` - DNA-derived data extension
- `edna2obis_report.html` - Comprehensive HTML report

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
   - Reduce `worms_n_proc` or `gbif_n_proc` in config.yaml
   - The pipeline has built-in retry logic

3. **Memory issues**
   - Reduce batch sizes in config
   - Close other applications
   - Consider running on a machine with more RAM

4. **Missing data files**
   - Verify all file paths in config.yaml are correct
   - Use absolute paths if relative paths don't work

#### Getting Help

- Check the HTML report for detailed error messages
- Review the terminal output for specific error details
- Ensure your input data follows the FAIRe NOAA format

### API Rate Limits

- **WoRMS**: Recommended max 3 processes (`worms_n_proc: 3`)
- **GBIF**: Can handle more processes (`gbif_n_proc: 0` uses all CPU cores)
- Large datasets may take 30+ minutes depending on size and API response times

### System Requirements

- **Minimum**: 4GB RAM, 2 CPU cores
- **Recommended**: 8GB+ RAM, 4+ CPU cores
- **Storage**: ~1GB free space for large datasets
- **Network**: Stable internet for API calls

### Quick Command Reference

```bash
# Setup (one time only)
git clone <repository-url>
cd edna2obis-3/edna2obis-3.0
conda env create -f environment.yml

# Every time you want to run
conda activate edna2obis
python main.py

# To update the environment (if environment.yml changes)
conda env update -f environment.yml
```

## Introduction
**Rationale:**

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

**Project abstract:**

Seawater was collected on board the NOAA ship Ronald H. Brown as part of the fourth Gulf of Mexico Ecosystems and Carbon Cycle (GOMECC-4) cruise from September 13 to October 21, 2021. Sampling for GOMECC-4 occurred along 16 coastal-offshore transects across the entire Gulf of Mexico and an additional line at 27N latitude in the Atlantic Ocean. We also collected eDNA samples near Padre Island National Seashore (U.S. National Parks Service), a barrier island located off the coast of south Texas. Vertical CTD sampling was employed at each site to measure discrete chemical, physical, and biological properties. Water sampling for DNA filtration was conducted at 54 sites and three depths per site (surface, deep chlorophyll maximum, and near bottom) to capture horizontal and vertical gradients of bacterial, protistan, and metazoan diversity across the Gulf. The resulting ASVs, their assigned taxonomy, and the metadata associated with theircollection are the input data for the OBIS conversion scripts presented here.

## Published data
- [GBIF](https://www.gbif.org/dataset/9012def0-bd87-48a0-ac9e-e0e78dd37689)
- [OBIS](https://obis.org/dataset/210efc7c-4762-47ee-b4b5-22a0f436ef44)

## NOAA Omics MIMARKS-based metadata template
This code was developed to convert a custom Google Sheet metadata template developed by NOAA Omics at AOML. To use the sheet for your own data, copy the Google Sheet to your Google Drive. Note that we have not tested the data validation functionality when downloading the Google Sheet to use as an Excel file.  

[AOML_MIMARKS.survey.water.6.0 v1.0.2](https://docs.google.com/spreadsheets/d/1YBXFU9PuMqm7IT1tp0LTxQ1v2j0tlCWFnhSpy-EBwPw/edit?usp=sharing)

## Requirements  

All dependencies are managed via conda. See setup instructions above for complete installation guidance.

**Key dependencies:**
- Python 3.8
- pandas, numpy - Data manipulation
- pyworms, pygbif - Taxonomic API clients  
- multiprocess - Parallel processing
- openpyxl - Excel file support
- pyyaml - Configuration files

**Full dependency list:** See `environment.yml`

## Current Repo Structure (v3.0)
```
edna2obis-3.0/
├── README.md                           # This file
├── environment.yml                     # Conda environment definition
├── config.yaml                         # Configuration file
├── main.py                            # Main pipeline script
├── src-v3/                            # Source code modules
│   ├── html_reporter.py               # HTML report generator
│   ├── create_occurrence_core/        # Occurrence core creation
│   ├── create_dna_derived_extension/   # DNA extension creation  
│   └── taxonomic_assignment/           # WoRMS/GBIF matching
├── raw-v3/                            # Input data directory
├── processed-v3/                      # Output files directory
└── images/                            # Report images/logos
```

## Usage

1. **Setup environment** (one time): Follow setup instructions above
2. **Configure your data**: Edit `config.yaml` with your file paths
3. **Run pipeline**: `python main.py`
4. **View results**: Open `edna2obis_report.html` in your browser

## Output Files

- `occurrence_worms_matched.csv` / `occurrence_gbif_matched.csv` - Darwin Core occurrence data
- `taxa_assignment_INFO.csv` - Summary of unique taxonomy assignments  
- `dna_derived_extension.csv` - DNA-derived data extension
- `edna2obis_report.html` - Comprehensive processing report

## Disclaimer  
This repository is a scientific product and is not official communication of the National Oceanic and Atmospheric Administration, or the United States Department of Commerce. All NOAA GitHub project code is provided on an 'as is' basis and the user assumes responsibility for its use. Any claims against the Department of Commerce or Department of Commerce bureaus stemming from the use of this GitHub project will be governed by all applicable Federal law. Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation or favoring by the Department of Commerce. The Department of Commerce seal and logo, or the seal and logo of a DOC bureau, shall not be used in any manner to imply endorsement of any commercial product or activity by DOC or the United States Government.
