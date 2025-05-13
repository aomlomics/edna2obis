# Darwin Core Conversion of eDNA Sequence Data From the AOML_MIMARKS metadata template 

**Version:** 1.0.0

**Author:** Katherine Silliman

**Last Updated:** 2-Oct-2023

This notebook is for converting a [MIMARKS](https://fairsharing.org/FAIRsharing.zvrep1)-based data sheet to DarwinCore for submission to OBIS. It has been testing on a Mac M1 laptop running in Rosetta mode, with Python 3.11. 

[Metadata template Google Sheet](https://docs.google.com/spreadsheets/d/1jof9MBEll7Xluu8-_znLRBIP9JpyAd_5YvdioZ-REoY/edit?usp=sharing)

**Requirements:**
- Python 3
- Python 3 packages:
    - os
- External packages:
    - Bio.Entrez from biopython
    - numpy
    - pandas
    - openpyxl
    - pyworms
    - multiprocess
- Custom modules:
    - WoRMS_matching

**Resources:**
- Abarenkov K, Andersson AF, Bissett A, Finstad AG, Fossøy F, Grosjean M, Hope M, Jeppesen TS, Kõljalg U, Lundin D, Nilsson RN, Prager M, Provoost P, Schigel D, Suominen S, Svenningsen C & Frøslev TG (2023) Publishing DNA-derived data through biodiversity data platforms, v1.3. Copenhagen: GBIF Secretariat. https://doi.org/10.35035/doc-vf1a-nr22.https://doi.org/10.35035/doc-vf1a-nr22.
- [OBIS manual](https://manual.obis.org/dna_data.html)
- [TDWG Darwin Core Occurrence Core](https://dwc.tdwg.org/terms/#occurrence)
- [GBIF DNA Derived Data Extension](https://tools.gbif.org/dwca-validator/extension.do?id=http://rs.gbif.org/terms/1.0/DNADerivedData)
- https://github.com/iobis/dataset-edna

**Citation**  
Silliman K, Anderson S, Storo R, Thompson L (2023) A Case Study in Sharing Marine eDNA Metabarcoding Data to OBIS. Biodiversity Information Science and Standards 7: e111048. https://doi.org/10.3897/biss.7.111048


## Installation  

```bash
conda create -n edna2obis
conda activate edna2obis
conda install -c conda-forge notebook
conda install -c conda-forge nb_conda_kernels

conda install -c conda-forge numpy pandas
conda install -c conda-forge openpyxl

#worms conversion
conda install -c conda-forge pyworms
conda install -c conda-forge multiprocess
conda install -c conda-forge biopython
```


```python
## Imports
import os

import numpy as np
import pandas as pd

import WoRMS_matching # custom functions for querying WoRMS API
```


```python
# jupyter notebook parameters
pd.set_option('display.max_colwidth', 150)
pd.set_option('display.max_columns', 50)
```

Note that in a Jupyter Notebook, the current directory is always where the .ipynb file is being run.

## Prepare input data 

**Project data and metadata**  
This workflow assumes that you have your project metadata in an Excel sheet formatted like the template located [here](https://docs.google.com/spreadsheets/d/1jof9MBEll7Xluu8-_znLRBIP9JpyAd_5YvdioZ-REoY/edit?usp=sharing). Instructions for filling out the metadata template are located in the 'Readme' sheet.

**eDNA and taxonomy data**  
The eDNA data and assigned taxonomy should be in a specific tab-delimited format. ![asv_table format](../images/asv_table.png)

This file is generated automatically by [Tourmaline v2023.5+](https://github.com/aomlomics/tourmaline), in X location. If your data was generated with Qiime2 or a previous version of Tourmaline, you can convert the `table.qza`, `taxonomy.qza`, and `repseqs.qza` outputs to the correct format using the `create_asv_seq_taxa_obis.sh` shell script.

Example:  

``` bash
#Run this with a qiime2 environment. 
bash create_asv_seq_taxa_obis.sh -f \
../gomecc_v2_raw/table-16S-merge.qza -t ../gomecc_v2_raw/taxonomy-16S-merge.qza -r ../gomecc_v2_raw/repseqs-16S-merge.qza \
-o ../gomecc_v2_raw/gomecc-16S-asv.tsv
```


## Set configs  

Below you can set definitions for parameters used in the code. 

| Parameter           | Description                                                                                                       | Example                                                                                              |
|---------------------|-------------------------------------------------------------------------------------------------------------------|------------------------------------------------------------------------------------------------------|
| `sample_data`       | Name of sheet in project data Excel file with sample data.                                                        | "water_sample_data"                                                                                  |
| `prep_data`         | Name of sheet in project data Excel file with data about molecular preparation methods.                           | "amplicon_prep_data"                                                                                 |
| `analysis_data`     | Name of sheet in project data Excel file with data about analysis methods.                                        | "analysis_data"                                                                                      |
| `study_data`        | Name of sheet in project data Excel file with metadata about the study.                                           | "study_data"                                                                                         |
| `msmt_metadata`     | Name of sheet in project data Excel file with metadata about additional measurements. Not used in current code.   | "measurement_metadata"                                                                               |
| `excel_file`        | Path of project data Excel file.                                                                                  | "../raw/gomecc4_AOML_MIMARKS.survey.water.6.0.xlsx"                                                  |
| `md_excel`          | Path of data dictionary Excel file.                                                                               | "../raw/gomecc_AOML2DwC standards.xlsx"                                                              |
| `datafiles`         | Python dictionary, where keys are the amplicon names and the values are the paths to the cooresponding ASV table. | {'16S V4-V5': '../raw/gomecc-16S-asv.tsv', '18S V9': '../raw/gomecc-18S-asv.tsv'}                    |
| `skip_sample_types` | Python list of sample_type values to skip from OBIS submission, such as controls or blanks.                       | ['mock community','distilled water blank','extraction blank','PCR no-template control','RTSF blank'] |
| `skip_columns`      | Python list of columns to ignore when submitting to OBIS.                                                         | ['notes_sampling']                                                                                   |


```python
params = {}
params['sample_data'] = "water_sample_data"
params['prep_data']= "amplicon_prep_data"
params['analysis_data'] = "analysis_data"
params['study_data'] = "study_data"
params['msmt_metadata'] = "measurement_metadata"
params['excel_file'] = "../raw/gomecc4_AOML_MIMARKS.survey.water.6.0.xlsx"

params['datafiles'] = {'16S V4-V5': '../raw/gomecc-16S-asv.tsv',
                       '18S V9': '../raw/gomecc-18S-asv.tsv'}

params['skip_sample_types'] = ['mock community','distilled water blank','extraction blank','PCR no-template control','RTSF blank']
params['skip_columns']= ['notes_sampling']
params['md_excel'] = "../raw/gomecc_AOML2DwC standards.xlsx"

```

## Load data

Note that in a Jupyter Notebook, the current directory is always where the .ipynb file is being run.

### Load project data Excel file


```python

data = pd.read_excel(
    params['excel_file'], 
    [params['study_data'],params['sample_data'],params['prep_data'],params['analysis_data'],params['msmt_metadata']],
    index_col=None, na_values=[""], comment="#"
)
```


```python
#rename keys in data dictionary to a general term
data['sample_data'] = data.pop(params['sample_data'])
data['prep_data'] = data.pop(params['prep_data'])
data['analysis_data'] = data.pop(params['analysis_data'])
data['study_data'] = data.pop(params['study_data'])
```


```python
#remove * from headers (was required for NCBI submission, but no longer needed)
data['sample_data'].columns = data['sample_data'].columns.str.replace("*","")
```

#### sample_data  
Contextual data about the samples collected, such as when it was collected, where it was collected from, what kind of sample it is, and what were the properties of the environment or experimental condition from which the sample was taken. Each row is a distinct sample, or Event. Most of this information is recorded during sample collection. This sheet contains terms from the MIMARKS survey water 6.0 package. 


```python
data['sample_data'].head()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>sample_name</th>
      <th>serial_number</th>
      <th>cruise_id</th>
      <th>line_id</th>
      <th>station</th>
      <th>ctd_bottle_no</th>
      <th>sample_replicate</th>
      <th>source_mat_id</th>
      <th>biological_replicates</th>
      <th>extract_number</th>
      <th>sample_title</th>
      <th>bioproject_accession</th>
      <th>biosample_accession</th>
      <th>amplicon_sequenced</th>
      <th>metagenome_sequenced</th>
      <th>organism</th>
      <th>collection_date_local</th>
      <th>collection_date</th>
      <th>depth</th>
      <th>env_broad_scale</th>
      <th>env_local_scale</th>
      <th>env_medium</th>
      <th>geo_loc_name</th>
      <th>lat_lon</th>
      <th>decimalLatitude</th>
      <th>...</th>
      <th>carbonate</th>
      <th>diss_inorg_carb</th>
      <th>diss_oxygen</th>
      <th>fluor</th>
      <th>hydrogen_ion</th>
      <th>nitrate</th>
      <th>nitrite</th>
      <th>nitrate_plus_nitrite</th>
      <th>omega_arag</th>
      <th>pco2</th>
      <th>ph</th>
      <th>phosphate</th>
      <th>pressure</th>
      <th>salinity</th>
      <th>samp_store_loc</th>
      <th>samp_store_temp</th>
      <th>silicate</th>
      <th>size_frac_low</th>
      <th>size_frac_up</th>
      <th>temp</th>
      <th>tot_alkalinity</th>
      <th>tot_depth_water_col</th>
      <th>transmittance</th>
      <th>date_sheet_modified</th>
      <th>modified_by</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>GOMECC4_27N_Sta1_Deep_A</td>
      <td>GOMECC4_001</td>
      <td>GOMECC-4 (2021)</td>
      <td>27N</td>
      <td>Sta1</td>
      <td>3</td>
      <td>A</td>
      <td>GOMECC4_27N_Sta1_Deep</td>
      <td>GOMECC4_27N_Sta1_Deep_B, GOMECC4_27N_Sta1_Deep_C</td>
      <td>Plate4_52</td>
      <td>Atlantic Ocean seawater sample GOMECC4_27N_Sta1_Deep_A</td>
      <td>PRJNA887898</td>
      <td>SAMN37516091</td>
      <td>16S V4-V5 | 18S V9</td>
      <td>planned for FY24</td>
      <td>seawater metagenome</td>
      <td>2021-09-14T11:00-04:00</td>
      <td>2021-09-14T07:00</td>
      <td>618 m</td>
      <td>marine biome [ENVO:00000447]</td>
      <td>marine mesopelagic zone [ENVO:00000213]</td>
      <td>sea water [ENVO:00002149]</td>
      <td>USA: Atlantic Ocean, east of Florida (27 N)</td>
      <td>26.997 N 79.618 W</td>
      <td>26.997</td>
      <td>...</td>
      <td>88.434 µmol/kg</td>
      <td>2215.45 µmol/kg</td>
      <td>129.44 µmol/kg</td>
      <td>0.0308</td>
      <td>0.0000000142 M</td>
      <td>29.3256 µmol/kg</td>
      <td>0.00391 µmol/kg</td>
      <td>29.3295 µmol/kg</td>
      <td>1.168</td>
      <td>624 µatm</td>
      <td>7.849</td>
      <td>1.94489 µmol/kg</td>
      <td>623 dbar</td>
      <td>34.946 psu</td>
      <td>NOAA/AOML Room 248</td>
      <td>-20 °C</td>
      <td>20.3569 µmol/kg</td>
      <td>no pre-filter</td>
      <td>0.22 µm</td>
      <td>7.479 °C</td>
      <td>2318.9 µmol/kg</td>
      <td>623 m</td>
      <td>4.7221</td>
      <td>2023-10-03 13:28:31.916</td>
      <td>luke.thompson@noaa.gov</td>
    </tr>
    <tr>
      <th>1</th>
      <td>GOMECC4_27N_Sta1_Deep_B</td>
      <td>GOMECC4_002</td>
      <td>GOMECC-4 (2021)</td>
      <td>27N</td>
      <td>Sta1</td>
      <td>3</td>
      <td>B</td>
      <td>GOMECC4_27N_Sta1_Deep</td>
      <td>GOMECC4_27N_Sta1_Deep_A, GOMECC4_27N_Sta1_Deep_C</td>
      <td>Plate4_60</td>
      <td>Atlantic Ocean seawater sample GOMECC4_27N_Sta1_Deep_B</td>
      <td>PRJNA887898</td>
      <td>SAMN37516092</td>
      <td>16S V4-V5 | 18S V9</td>
      <td>planned for FY24</td>
      <td>seawater metagenome</td>
      <td>2021-09-14T11:00-04:00</td>
      <td>2021-09-14T07:00</td>
      <td>618 m</td>
      <td>marine biome [ENVO:00000447]</td>
      <td>marine mesopelagic zone [ENVO:00000213]</td>
      <td>sea water [ENVO:00002149]</td>
      <td>USA: Atlantic Ocean, east of Florida (27 N)</td>
      <td>26.997 N 79.618 W</td>
      <td>26.997</td>
      <td>...</td>
      <td>88.434 µmol/kg</td>
      <td>2215.45 µmol/kg</td>
      <td>129.44 µmol/kg</td>
      <td>0.0308</td>
      <td>0.0000000142 M</td>
      <td>29.3256 µmol/kg</td>
      <td>0.00391 µmol/kg</td>
      <td>29.3295 µmol/kg</td>
      <td>1.168</td>
      <td>624 µatm</td>
      <td>7.849</td>
      <td>1.94489 µmol/kg</td>
      <td>623 dbar</td>
      <td>34.946 psu</td>
      <td>NOAA/AOML Room 248</td>
      <td>-20 °C</td>
      <td>20.3569 µmol/kg</td>
      <td>no pre-filter</td>
      <td>0.22 µm</td>
      <td>7.479 °C</td>
      <td>2318.9 µmol/kg</td>
      <td>623 m</td>
      <td>4.7221</td>
      <td>NaT</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>2</th>
      <td>GOMECC4_27N_Sta1_Deep_C</td>
      <td>GOMECC4_003</td>
      <td>GOMECC-4 (2021)</td>
      <td>27N</td>
      <td>Sta1</td>
      <td>3</td>
      <td>C</td>
      <td>GOMECC4_27N_Sta1_Deep</td>
      <td>GOMECC4_27N_Sta1_Deep_A, GOMECC4_27N_Sta1_Deep_B</td>
      <td>Plate4_62</td>
      <td>Atlantic Ocean seawater sample GOMECC4_27N_Sta1_Deep_C</td>
      <td>PRJNA887898</td>
      <td>SAMN37516093</td>
      <td>16S V4-V5 | 18S V9</td>
      <td>planned for FY24</td>
      <td>seawater metagenome</td>
      <td>2021-09-14T11:00-04:00</td>
      <td>2021-09-14T07:00</td>
      <td>618 m</td>
      <td>marine biome [ENVO:00000447]</td>
      <td>marine mesopelagic zone [ENVO:00000213]</td>
      <td>sea water [ENVO:00002149]</td>
      <td>USA: Atlantic Ocean, east of Florida (27 N)</td>
      <td>26.997 N 79.618 W</td>
      <td>26.997</td>
      <td>...</td>
      <td>88.434 µmol/kg</td>
      <td>2215.45 µmol/kg</td>
      <td>129.44 µmol/kg</td>
      <td>0.0308</td>
      <td>0.0000000142 M</td>
      <td>29.3256 µmol/kg</td>
      <td>0.00391 µmol/kg</td>
      <td>29.3295 µmol/kg</td>
      <td>1.168</td>
      <td>624 µatm</td>
      <td>7.849</td>
      <td>1.94489 µmol/kg</td>
      <td>623 dbar</td>
      <td>34.946 psu</td>
      <td>NOAA/AOML Room 248</td>
      <td>-20 °C</td>
      <td>20.3569 µmol/kg</td>
      <td>no pre-filter</td>
      <td>0.22 µm</td>
      <td>7.479 °C</td>
      <td>2318.9 µmol/kg</td>
      <td>623 m</td>
      <td>4.7221</td>
      <td>NaT</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>3</th>
      <td>GOMECC4_27N_Sta1_DCM_A</td>
      <td>GOMECC4_004</td>
      <td>GOMECC-4 (2021)</td>
      <td>27N</td>
      <td>Sta1</td>
      <td>14</td>
      <td>A</td>
      <td>GOMECC4_27N_Sta1_DCM</td>
      <td>GOMECC4_27N_Sta1_DCM_B, GOMECC4_27N_Sta1_DCM_C</td>
      <td>Plate4_53</td>
      <td>Atlantic Ocean seawater sample GOMECC4_27N_Sta1_DCM_A</td>
      <td>PRJNA887898</td>
      <td>SAMN37516094</td>
      <td>16S V4-V5 | 18S V9</td>
      <td>planned for FY24</td>
      <td>seawater metagenome</td>
      <td>2021-09-14T11:00-04:00</td>
      <td>2021-09-14T07:00</td>
      <td>49 m</td>
      <td>marine biome [ENVO:00000447]</td>
      <td>marine photic zone [ENVO:00000209]</td>
      <td>sea water [ENVO:00002149]</td>
      <td>USA: Atlantic Ocean, east of Florida (27 N)</td>
      <td>26.997 N 79.618 W</td>
      <td>26.997</td>
      <td>...</td>
      <td>229.99 µmol/kg</td>
      <td>2033.19 µmol/kg</td>
      <td>193.443 µmol/kg</td>
      <td>0.036</td>
      <td>0.0000000094 M</td>
      <td>0 µmol/kg</td>
      <td>0 µmol/kg</td>
      <td>0 µmol/kg</td>
      <td>3.805</td>
      <td>423 µatm</td>
      <td>8.027</td>
      <td>0.0517 µmol/kg</td>
      <td>49 dbar</td>
      <td>36.325 psu</td>
      <td>NOAA/AOML Room 248</td>
      <td>-20 °C</td>
      <td>1.05635 µmol/kg</td>
      <td>no pre-filter</td>
      <td>0.22 µm</td>
      <td>28.592 °C</td>
      <td>2371 µmol/kg</td>
      <td>623 m</td>
      <td>4.665</td>
      <td>NaT</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>4</th>
      <td>GOMECC4_27N_Sta1_DCM_B</td>
      <td>GOMECC4_005</td>
      <td>GOMECC-4 (2021)</td>
      <td>27N</td>
      <td>Sta1</td>
      <td>14</td>
      <td>B</td>
      <td>GOMECC4_27N_Sta1_DCM</td>
      <td>GOMECC4_27N_Sta1_DCM_A, GOMECC4_27N_Sta1_DCM_C</td>
      <td>Plate4_46</td>
      <td>Atlantic Ocean seawater sample GOMECC4_27N_Sta1_DCM_B</td>
      <td>PRJNA887898</td>
      <td>SAMN37516095</td>
      <td>16S V4-V5 | 18S V9</td>
      <td>planned for FY24</td>
      <td>seawater metagenome</td>
      <td>2021-09-14T11:00-04:00</td>
      <td>2021-09-14T07:00</td>
      <td>49 m</td>
      <td>marine biome [ENVO:00000447]</td>
      <td>marine photic zone [ENVO:00000209]</td>
      <td>sea water [ENVO:00002149]</td>
      <td>USA: Atlantic Ocean, east of Florida (27 N)</td>
      <td>26.997 N 79.618 W</td>
      <td>26.997</td>
      <td>...</td>
      <td>229.99 µmol/kg</td>
      <td>2033.19 µmol/kg</td>
      <td>193.443 µmol/kg</td>
      <td>0.036</td>
      <td>0.0000000094 M</td>
      <td>0 µmol/kg</td>
      <td>0 µmol/kg</td>
      <td>0 µmol/kg</td>
      <td>3.805</td>
      <td>423 µatm</td>
      <td>8.027</td>
      <td>0.0517 µmol/kg</td>
      <td>49 dbar</td>
      <td>36.325 psu</td>
      <td>NOAA/AOML Room 248</td>
      <td>-20 °C</td>
      <td>1.05635 µmol/kg</td>
      <td>no pre-filter</td>
      <td>0.22 µm</td>
      <td>28.592 °C</td>
      <td>2371 µmol/kg</td>
      <td>623 m</td>
      <td>4.665</td>
      <td>NaT</td>
      <td>NaN</td>
    </tr>
  </tbody>
</table>
<p>5 rows × 79 columns</p>
</div>



#### prep_data  
Contextual data about how the samples were prepared for sequencing. Includes how they were extracted, what amplicon was targeted, how they were sequenced. Each row is a separate sequencing library preparation, distinguished by a unique library_id.


```python
data['prep_data'].head(2)
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>sample_name</th>
      <th>library_id</th>
      <th>title</th>
      <th>library_strategy</th>
      <th>library_source</th>
      <th>library_selection</th>
      <th>lib_layout</th>
      <th>platform</th>
      <th>instrument_model</th>
      <th>design_description</th>
      <th>filetype</th>
      <th>filename</th>
      <th>filename2</th>
      <th>drive_location</th>
      <th>biosample_accession</th>
      <th>sra_accession</th>
      <th>seq_method</th>
      <th>nucl_acid_ext</th>
      <th>amplicon_sequenced</th>
      <th>target_gene</th>
      <th>target_subfragment</th>
      <th>pcr_primer_forward</th>
      <th>pcr_primer_reverse</th>
      <th>pcr_primer_name_forward</th>
      <th>pcr_primer_name_reverse</th>
      <th>pcr_primer_reference</th>
      <th>pcr_cond</th>
      <th>nucl_acid_amp</th>
      <th>adapters</th>
      <th>mid_barcode</th>
      <th>date_sheet_modified</th>
      <th>modified_by</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>GOMECC4_NegativeControl_1</td>
      <td>GOMECC16S_Neg1</td>
      <td>16S amplicon metabarcoding of marine metagenome: Gulf of Mexico (USA)</td>
      <td>AMPLICON</td>
      <td>METAGENOMIC</td>
      <td>PCR</td>
      <td>paired</td>
      <td>ILLUMINA</td>
      <td>Illumina MiSeq</td>
      <td>Samples were collected and filtered onto Sterivex 0.22 um cartridge filters. DNA was extracted from Sterivex by adding lysis buffer and magnetic b...</td>
      <td>fastq</td>
      <td>GOMECC16S_Neg1_S499_L001_R1_001.fastq.gz</td>
      <td>GOMECC16S_Neg1_S499_L001_R2_001.fastq.gz</td>
      <td>NaN</td>
      <td>SAMN37516589</td>
      <td>SRR26148505</td>
      <td>Illumina MiSeq 2x250</td>
      <td>https://github.com/aomlomics/protocols/blob/main/protocol_DNA_extraction_Sterivex.md</td>
      <td>16S V4-V5</td>
      <td>16S rRNA</td>
      <td>V4-V5</td>
      <td>GTGYCAGCMGCCGCGGTAA</td>
      <td>CCGYCAATTYMTTTRAGTTT</td>
      <td>515F-Y</td>
      <td>926R</td>
      <td>10.1111/1462-2920.13023</td>
      <td>initial denaturation:95_2;denaturation:95_0.75;annealing:50_0.75;elongation:68_1.5;final elongation:68_5;25</td>
      <td>10.1111/1462-2920.13023</td>
      <td>ACACTGACGACATGGTTCTACA;TACGGTAGCAGAGACTTGGTCT</td>
      <td>missing: not provided</td>
      <td>NaT</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>1</th>
      <td>GOMECC4_27N_Sta1_DCM_A</td>
      <td>GOMECC18S_Plate4_53</td>
      <td>18S amplicon metabarcoding of marine metagenome: Gulf of Mexico (USA)</td>
      <td>AMPLICON</td>
      <td>METAGENOMIC</td>
      <td>PCR</td>
      <td>paired</td>
      <td>ILLUMINA</td>
      <td>Illumina MiSeq</td>
      <td>Samples were collected and filtered onto Sterivex 0.22 um cartridge filters. DNA was extracted from Sterivex by adding lysis buffer and magnetic b...</td>
      <td>fastq</td>
      <td>GOMECC18S_Plate4_53_S340_L001_R1_001.fastq.gz</td>
      <td>GOMECC18S_Plate4_53_S340_L001_R2_001.fastq.gz</td>
      <td>NaN</td>
      <td>SAMN37516094</td>
      <td>SRR26161153</td>
      <td>Illumina MiSeq 2x250</td>
      <td>https://github.com/aomlomics/protocols/blob/main/protocol_DNA_extraction_Sterivex.md</td>
      <td>18S V9</td>
      <td>18S rRNA</td>
      <td>V9</td>
      <td>GTACACACCGCCCGTC</td>
      <td>TGATCCTTCTGCAGGTTCACCTAC</td>
      <td>1391f</td>
      <td>EukBr</td>
      <td>10.1371/journal.pone.0006372</td>
      <td>initial denaturation:94_3;denaturation:94_0.75;annealing:65_0.25;57_0.5;elongation:72_1.5;final elongation:72_10;35</td>
      <td>10.1371/journal.pone.0006372</td>
      <td>ACACTGACGACATGGTTCTACA;TACGGTAGCAGAGACTTGGTCT</td>
      <td>missing: not provided</td>
      <td>2023-10-03 12:49:23.878</td>
      <td>katherine.silliman@noaa.gov</td>
    </tr>
  </tbody>
</table>
</div>



### Load ASV data  
There is one ASV file for each marker that was sequenced. The ASV data files have one row for each unique amplicon sequence variants (ASVs). They contain the ASV DNA sequence, a unique hash identifier the taxonomic assignment for each ASV, the confidence given that assignment by the naive-bayes classifier, and then the number of reads observed in each sample. 

This file is created automatically with [Tourmaline v.2023.5+](https://github.com/aomlomics/tourmaline), and is found in `01-taxonomy/asv_taxa_sample_table.tsv`. 

| column name    | definition                                                                                                                                                                                                                                                                                                                                                                                              |
|----------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| featureid      | A hash of the ASV sequence, used as a unique identifier for the ASV.                                                                                                                                                                                                                                                                                                                                    |
| sequence       | The DNA sequence of the ASV                                                                                                                                                                                                                                                                                                                                                                             |
| taxonomy       | The full taxonomy assigned to an ASV sequence. This string could be formatted in very different ways depending on the reference database used during classification, however it should always be in reverse rank order separated by ;. We provide examples for how to process results from a Silva classifier and the PR2 18S classifier. For other taxonomy formats, the code will need to be adapted. |
| Confidence     | This is the confidence score assigned the taxonomic classification with a naive-bayes classifier.                                                                                                                                                                                                                                                                                                       |
| sample columns | The next columns each represent a sample (or eventID), and the number of reads for that ASV observed in the sample.                                                                                                                                                                                                                                                                                     |


```python
# read in ASV tables, looping through amplicons
asv_tables = {}

for gene in params['datafiles'].keys():
    asv_tables[gene] = pd.read_table(params['datafiles'][gene])

```


```python
asv_tables.keys()
```




    dict_keys(['16S V4-V5', '18S V9'])




```python
asv_tables['16S V4-V5'].iloc[:,0:20].head()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>featureid</th>
      <th>sequence</th>
      <th>taxonomy</th>
      <th>Confidence</th>
      <th>GOMECC4_27N_Sta1_DCM_A</th>
      <th>GOMECC4_27N_Sta1_DCM_B</th>
      <th>GOMECC4_27N_Sta1_DCM_C</th>
      <th>GOMECC4_27N_Sta1_Deep_A</th>
      <th>GOMECC4_27N_Sta1_Deep_B</th>
      <th>GOMECC4_27N_Sta1_Deep_C</th>
      <th>GOMECC4_27N_Sta1_Surface_A</th>
      <th>GOMECC4_27N_Sta1_Surface_B</th>
      <th>GOMECC4_27N_Sta4_DCM_A</th>
      <th>GOMECC4_27N_Sta4_DCM_B</th>
      <th>GOMECC4_27N_Sta4_DCM_C</th>
      <th>GOMECC4_27N_Sta4_Deep_A</th>
      <th>GOMECC4_27N_Sta4_Deep_B</th>
      <th>GOMECC4_27N_Sta4_Deep_C</th>
      <th>GOMECC4_27N_Sta4_Surface_A</th>
      <th>GOMECC4_27N_Sta4_Surface_B</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>00006f0784f7dbb2f162408abb6da629</td>
      <td>TACGGAGGGTGCGAGCGTTAATCGGAATTACTGGGCGTAAAGCGCATGCAGGTGGTTTGTTAAGTCAGATGTGAAAGCCCGGGGCTCAACCTCGGAATTGCATTTGAAACTGGCAGACTAGAGTACTGTAGAGGGGGGTAGAATTT...</td>
      <td>d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Vibrionales; f__Vibrionaceae; g__Vibrio</td>
      <td>0.978926</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>25</td>
    </tr>
    <tr>
      <th>1</th>
      <td>000094731d4984ed41435a1bf65b7ef2</td>
      <td>TACAGAGAGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCGCGTAGGTGGGTATTTAAGTCGGATGTGAAATCCCCGGGCTTAACCTGGGAACTGCATCCGAAACTATTTAACTAGAGTATGGGAGAGGTAAGTAGAATTT...</td>
      <td>d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__HOC36; f__HOC36; g__HOC36; s__Candidatus_Thioglobus</td>
      <td>0.881698</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>2</th>
      <td>0001a3c11fcef1b1b8f4c72942efbbac</td>
      <td>TACGAAGGGGGCGAGCGTTGTTCGGAATTACTGGGCGTAAAGGGCGCGTAGGCGGTCTTCTAAGTTAGGCGTGAAAGCCCCGGGCTCAACCTGGGAACTGCGCTTAATACTGGAAGACTAGAAAACGGAAGAGGGTAGTGGAATTC...</td>
      <td>d__Bacteria; p__Cyanobacteria; c__Cyanobacteriia; o__Synechococcales; f__Cyanobiaceae; g__Cyanobium_PCC-6307</td>
      <td>0.762793</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>3</th>
      <td>0001ceef5162e6d689ef30418cfcc164</td>
      <td>TACAGAGGGTGCAAGCGTTGTTCGGAATCATTGGGCGTAAAGCGCGCGTAGGCGGCCAAATAAGTCTGATGTGAAGGCCCAGGGCTCAACCCTGGAAGTGCATCGGAAACTGTTTGGCTCGAGTCCCGGAGGGGGTGGTGGAATTC...</td>
      <td>d__Bacteria; p__Myxococcota; c__Myxococcia; o__Myxococcales; f__Myxococcaceae; g__P3OB-42; s__uncultured_bacterium</td>
      <td>0.997619</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>4</th>
      <td>000235534662df05bb30219a4b978dac</td>
      <td>TACGGAAGGTCCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCGCGTAGGTGGTTTTTTAAGTTGGATGTGAAAGCCCTGGGCTCAACCTAGGAACTGCATCCAAAACTAGATGACTAGAGTACGAAAGAGGGAAGTAGAATTC...</td>
      <td>d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__SAR86_clade; f__SAR86_clade; g__SAR86_clade</td>
      <td>0.999961</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
  </tbody>
</table>
</div>



### Drop samples with unwanted sample types  

Often with eDNA projects, we have control samples that are sequenced along with our survey samples. These can include filtering distilled water, using pure water instead of DNA in a PCR or DNA extraction protocol, or a mock community of known microbial taxa. Controls can help identify and mitigate contaminant DNA in our samples, but are not useful for biodiversity platforms like OBIS. You can select which sample_type values to drop with the `skip_sample_types` parameter.


```python
samps_to_remove = data['sample_data']['sample_type'].isin(params['skip_sample_types'])
#data['sample_data'][samps_to_remove]
# list of samples to drop
samples_to_drop = data['sample_data']['sample_name'][samps_to_remove]
```

You can view the list of samples to be dropped below.


```python
samples_to_drop
```




    26     GOMECC4_Blank_DIW_20210915_A
    27     GOMECC4_Blank_DIW_20210915_B
    28     GOMECC4_Blank_DIW_20210915_C
    200    GOMECC4_Blank_DIW_20210930_A
    201    GOMECC4_Blank_DIW_20210930_B
    202    GOMECC4_Blank_DIW_20210930_C
    334    GOMECC4_Blank_DIW_20211011_A
    335    GOMECC4_Blank_DIW_20211011_B
    336    GOMECC4_Blank_DIW_20211011_C
    409    GOMECC4_Blank_DIW_20211016_A
    410    GOMECC4_Blank_DIW_20211016_B
    411    GOMECC4_Blank_DIW_20211016_C
    484       GOMECC4_ExtractionBlank_1
    485      GOMECC4_ExtractionBlank_11
    486      GOMECC4_ExtractionBlank_12
    487       GOMECC4_ExtractionBlank_3
    488       GOMECC4_ExtractionBlank_5
    489       GOMECC4_ExtractionBlank_7
    490       GOMECC4_ExtractionBlank_9
    491            GOMECC4_MSUControl_1
    492            GOMECC4_MSUControl_2
    493            GOMECC4_MSUControl_3
    494            GOMECC4_MSUControl_4
    495            GOMECC4_MSUControl_5
    496            GOMECC4_MSUControl_6
    497            GOMECC4_MSUControl_7
    498       GOMECC4_NegativeControl_1
    499       GOMECC4_NegativeControl_2
    500       GOMECC4_PositiveControl_1
    501       GOMECC4_PositiveControl_2
    Name: sample_name, dtype: object




```python
# remove samples from sample_data sheet
data['sample_data'] = data['sample_data'][~samps_to_remove]
```


```python
# check the sample_type values left in your sample_data. We only want seawater.
data['sample_data']['sample_type'].unique()
```




    array(['seawater'], dtype=object)




```python
# remove samples from prep_data
prep_samps_to_remove = data['prep_data']['sample_name'].isin(samples_to_drop)
data['prep_data'] = data['prep_data'][~prep_samps_to_remove]
```

##### drop unwanted samples from ASV files



```python
for gene in params['datafiles'].keys():
    asv_tables[gene] = asv_tables[gene].drop(columns=samples_to_drop,errors='ignore')
```

### Drop columns with all NAs  

If your project data file has columns with only NAs, this code will check for those, provide their column headers for verification, then remove them.


```python
# which have all NAs?
dropped = pd.DataFrame()
for sheet in ['sample_data','prep_data','analysis_data']:
    res = pd.Series(data[sheet].columns[data[sheet].isnull().all(0)],
                name=sheet)
    dropped=pd.concat([dropped,res],axis=1)
    
```

Which columns in each sheet have only NA values?


```python
dropped
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>sample_data</th>
      <th>prep_data</th>
      <th>analysis_data</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>samp_size</td>
      <td>drive_location</td>
      <td>sop</td>
    </tr>
  </tbody>
</table>
</div>



If you are fine with leaving these columns out, proceed:


```python
for sheet in ['sample_data','prep_data','analysis_data']:
    data[sheet].dropna(axis=1, how='all',inplace=True)
```

Now let's check which columns have missing values in some of the rows. These should be filled in on the Excel sheet with the appropriate term ('not applicable', 'missing', or 'not collected'). Alternatively, you can drop the column if it is not needed for submission to OBIS.


```python
# which columns have missing data (NAs) in some rows
some = pd.DataFrame()
for sheet in ['sample_data','prep_data','analysis_data']:
    res = pd.Series(data[sheet].columns[data[sheet].isnull().any()].tolist(),
                name=sheet)
    some=pd.concat([some,res],axis=1)
```


```python
some
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>sample_data</th>
      <th>prep_data</th>
      <th>analysis_data</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>notes_bottle_metadata</td>
      <td>date_sheet_modified</td>
      <td>date_sheet_modified</td>
    </tr>
    <tr>
      <th>1</th>
      <td>date_sheet_modified</td>
      <td>modified_by</td>
      <td>modified_by</td>
    </tr>
    <tr>
      <th>2</th>
      <td>modified_by</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
  </tbody>
</table>
</div>



Here I'm going to drop all the columns with some missing data, as I don't need them for submission to OBIS.


```python
# drop columns with any missing data
for sheet in ['sample_data','prep_data','analysis_data']:
    data[sheet].dropna(axis=1, how='any',inplace=True)
```

### Load data dictionary Excel file 
This Excel file is used as a data dictionary for converting between terms used in the project data Excel file and Darwin Core terms for submission to OBIS. Currently, we are only preparing an Occurrence core file and a DNA-derived extension file, with Event information in the Occurrence file. Future versions of this workflow will prepare an extendedMeasurementOrFact file as well.


```python
# read in data dictionary excel file
dwc_data = pd.read_excel(
    params['md_excel'], 
    ['event','occurrence','dna'],
    index_col=0, na_values=[""]
)
```


```python
#example of a sheet in the data dictionary
dwc_data['event'].head()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>AOML_term</th>
      <th>AOML_file</th>
      <th>DwC_definition</th>
    </tr>
    <tr>
      <th>DwC_term</th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>eventID</th>
      <td>sample_name</td>
      <td>sample_data</td>
      <td>An identifier for the set of information associated with a dwc:Event (something that occurs at a place and time). https://dwc.tdwg.org/terms/#dwc:...</td>
    </tr>
    <tr>
      <th>eventDate</th>
      <td>collection_date_local</td>
      <td>sample_data</td>
      <td>this is the date-time when the dwc:Event was recorded. Recommended best practice is to use a date that conforms to ISO 8601-1:2019. https://dwc.td...</td>
    </tr>
    <tr>
      <th>samplingProtocol</th>
      <td>collection_method</td>
      <td>sample_data</td>
      <td>The names of, references to, or descriptions of the methods or protocols used during a dwc:Event.</td>
    </tr>
    <tr>
      <th>locationID</th>
      <td>station</td>
      <td>sample_data</td>
      <td>An identifier for the set of dcterms:Location information. May be a global unique identifier or an identifier specific to the data set.</td>
    </tr>
    <tr>
      <th>decimalLatitude</th>
      <td>decimalLatitude</td>
      <td>sample_data</td>
      <td>The geographic latitude (in decimal degrees, using the spatial reference system given in dwc:geodeticDatum) of the geographic center of a dcterms:...</td>
    </tr>
  </tbody>
</table>
</div>



## Convert to Occurrence file
In order to link the DNA-derived extension metadata to our OBIS occurrence records, we have to use the Occurrence core. For this data set, a `parentEvent` is a filtered water sample that was DNA extracted, a sequencing library from that DNA extraction is an `event`, and an `occurrence` is an ASV observed within a library. We will have an an occurence file and a DNA derived data file. Future versions will generate a measurements file.   
**Define files**


### Sampling event info 




```python
dwc_data['event']
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>AOML_term</th>
      <th>AOML_file</th>
      <th>DwC_definition</th>
    </tr>
    <tr>
      <th>DwC_term</th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>eventID</th>
      <td>sample_name</td>
      <td>sample_data</td>
      <td>An identifier for the set of information associated with a dwc:Event (something that occurs at a place and time). https://dwc.tdwg.org/terms/#dwc:...</td>
    </tr>
    <tr>
      <th>eventDate</th>
      <td>collection_date_local</td>
      <td>sample_data</td>
      <td>this is the date-time when the dwc:Event was recorded. Recommended best practice is to use a date that conforms to ISO 8601-1:2019. https://dwc.td...</td>
    </tr>
    <tr>
      <th>samplingProtocol</th>
      <td>collection_method</td>
      <td>sample_data</td>
      <td>The names of, references to, or descriptions of the methods or protocols used during a dwc:Event.</td>
    </tr>
    <tr>
      <th>locationID</th>
      <td>station</td>
      <td>sample_data</td>
      <td>An identifier for the set of dcterms:Location information. May be a global unique identifier or an identifier specific to the data set.</td>
    </tr>
    <tr>
      <th>decimalLatitude</th>
      <td>decimalLatitude</td>
      <td>sample_data</td>
      <td>The geographic latitude (in decimal degrees, using the spatial reference system given in dwc:geodeticDatum) of the geographic center of a dcterms:...</td>
    </tr>
    <tr>
      <th>decimalLongitude</th>
      <td>decimalLongitude</td>
      <td>sample_data</td>
      <td>The geographic longitude (in decimal degrees, using the spatial reference system given in dwc:geodeticDatum) of the geographic center of a dcterms...</td>
    </tr>
    <tr>
      <th>geodeticDatum</th>
      <td>none</td>
      <td>pipeline</td>
      <td>The ellipsoid, geodetic datum, or spatial reference system (SRS) upon which the geographic coordinates given in dwc:decimalLatitude and dwc:decima...</td>
    </tr>
    <tr>
      <th>countryCode</th>
      <td>none</td>
      <td>pipeline</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>minimumDepthInMeters</th>
      <td>depth</td>
      <td>sample_data</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>maximumDepthInMeters</th>
      <td>derived: depth</td>
      <td>sample_data</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>datasetID</th>
      <td>project_id_external</td>
      <td>study_data</td>
      <td>An identifier for the set of data. May be a global unique identifier or an identifier specific to a collection or institution.</td>
    </tr>
    <tr>
      <th>waterBody</th>
      <td>derived</td>
      <td>sample_data</td>
      <td>The name of the water body in which the dcterms:Location occurs.         Recommended best practice is to use a controlled vocabulary such as the G...</td>
    </tr>
    <tr>
      <th>locality</th>
      <td>geo_loc_name</td>
      <td>sample_data</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>eventRemarks</th>
      <td>derived: controls_used</td>
      <td>analysis_data</td>
      <td>Comments or notes about the dwc:Event.</td>
    </tr>
  </tbody>
</table>
</div>




```python
event_dict = dwc_data['event'].to_dict('index')
```


```python
event_dict['eventID']
```




    {'AOML_term': 'sample_name',
     'AOML_file': 'sample_data',
     'DwC_definition': 'An identifier for the set of information associated with a dwc:Event (something that occurs at a place and time). https://dwc.tdwg.org/terms/#dwc:eventID'}




```python
# check which event terms are not in sample_data sheet
for key in event_dict.keys():
    if event_dict[key]['AOML_file'] == 'sample_data':
        if event_dict[key]['AOML_term'] not in data['sample_data'].columns:
            print(key,event_dict[key])
```

    maximumDepthInMeters {'AOML_term': 'derived: depth', 'AOML_file': 'sample_data', 'DwC_definition': nan}
    waterBody {'AOML_term': 'derived', 'AOML_file': 'sample_data', 'DwC_definition': 'The name of the water body in which the dcterms:Location occurs.         Recommended best practice is to use a controlled vocabulary such as the Getty Thesaurus of Geographic Names.'}



```python
# custom add waterBody

data['sample_data'].loc[data['sample_data']['geo_loc_name'].str.contains("Atlantic Ocean"), 'waterBody']= "Atlantic Ocean"
data['sample_data'].loc[data['sample_data']['geo_loc_name'].str.contains("Gulf"), 'waterBody']= "Mexico, Gulf of"

```

    /var/folders/_f/3hfwkwps2rq9q60vkv4fnd_n9rf1vk/T/ipykernel_5945/2121370784.py:3: FutureWarning: Setting an item of incompatible dtype is deprecated and will raise in a future error of pandas. Value 'Atlantic Ocean' has dtype incompatible with float64, please explicitly cast to a compatible dtype first.
      data['sample_data'].loc[data['sample_data']['geo_loc_name'].str.contains("Atlantic Ocean"), 'waterBody']= "Atlantic Ocean"



```python
# change locationID to line_id+station
data['sample_data']['station'] = data['sample_data']['line_id']+ "_"+data['sample_data']['station'] 

```


```python
# rename sample_data columns to fit DwC standard
gen = (x for x in event_dict.keys() if event_dict[x]['AOML_file'] == 'sample_data')
rename_dict = {}
for x in gen:
    #print(x)
    rename_dict[event_dict[x]['AOML_term']] = x

event_sample = data['sample_data'].rename(columns=rename_dict)
event_sample = event_sample.drop(columns=[col for col in event_sample if col not in rename_dict.values()])

```


```python
# add maximumDepthInMeters
#remove m in depth
event_sample['minimumDepthInMeters'] = event_sample['minimumDepthInMeters'].str.strip(" m")
event_sample['maximumDepthInMeters'] = event_sample['minimumDepthInMeters']
```


```python

event_sample.head()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>eventID</th>
      <th>locationID</th>
      <th>eventDate</th>
      <th>minimumDepthInMeters</th>
      <th>locality</th>
      <th>decimalLatitude</th>
      <th>decimalLongitude</th>
      <th>samplingProtocol</th>
      <th>waterBody</th>
      <th>maximumDepthInMeters</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>GOMECC4_27N_Sta1_Deep_A</td>
      <td>27N_Sta1</td>
      <td>2021-09-14T11:00-04:00</td>
      <td>618</td>
      <td>USA: Atlantic Ocean, east of Florida (27 N)</td>
      <td>26.997</td>
      <td>-79.618</td>
      <td>CTD rosette</td>
      <td>Atlantic Ocean</td>
      <td>618</td>
    </tr>
    <tr>
      <th>1</th>
      <td>GOMECC4_27N_Sta1_Deep_B</td>
      <td>27N_Sta1</td>
      <td>2021-09-14T11:00-04:00</td>
      <td>618</td>
      <td>USA: Atlantic Ocean, east of Florida (27 N)</td>
      <td>26.997</td>
      <td>-79.618</td>
      <td>CTD rosette</td>
      <td>Atlantic Ocean</td>
      <td>618</td>
    </tr>
    <tr>
      <th>2</th>
      <td>GOMECC4_27N_Sta1_Deep_C</td>
      <td>27N_Sta1</td>
      <td>2021-09-14T11:00-04:00</td>
      <td>618</td>
      <td>USA: Atlantic Ocean, east of Florida (27 N)</td>
      <td>26.997</td>
      <td>-79.618</td>
      <td>CTD rosette</td>
      <td>Atlantic Ocean</td>
      <td>618</td>
    </tr>
    <tr>
      <th>3</th>
      <td>GOMECC4_27N_Sta1_DCM_A</td>
      <td>27N_Sta1</td>
      <td>2021-09-14T11:00-04:00</td>
      <td>49</td>
      <td>USA: Atlantic Ocean, east of Florida (27 N)</td>
      <td>26.997</td>
      <td>-79.618</td>
      <td>CTD rosette</td>
      <td>Atlantic Ocean</td>
      <td>49</td>
    </tr>
    <tr>
      <th>4</th>
      <td>GOMECC4_27N_Sta1_DCM_B</td>
      <td>27N_Sta1</td>
      <td>2021-09-14T11:00-04:00</td>
      <td>49</td>
      <td>USA: Atlantic Ocean, east of Florida (27 N)</td>
      <td>26.997</td>
      <td>-79.618</td>
      <td>CTD rosette</td>
      <td>Atlantic Ocean</td>
      <td>49</td>
    </tr>
  </tbody>
</table>
</div>




```python
# add amplicon_sequenced back 
event_sample['amplicon_sequenced'] = data['sample_data']['amplicon_sequenced']
```

Now add an event for each sequencing library, with replicate water sample as the parentEvent.  

**Future Update**: make this a for loop


```python
child_data_16S = event_sample[event_sample['amplicon_sequenced'].str.contains('16S V4-V5')].copy()
child_data_16S['parentEventID'] = child_data_16S['eventID']
child_data_16S['eventID'] = child_data_16S['eventID']+"_16S"
child_data_16S.head()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>eventID</th>
      <th>locationID</th>
      <th>eventDate</th>
      <th>minimumDepthInMeters</th>
      <th>locality</th>
      <th>decimalLatitude</th>
      <th>decimalLongitude</th>
      <th>samplingProtocol</th>
      <th>waterBody</th>
      <th>maximumDepthInMeters</th>
      <th>amplicon_sequenced</th>
      <th>parentEventID</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>GOMECC4_27N_Sta1_Deep_A_16S</td>
      <td>27N_Sta1</td>
      <td>2021-09-14T11:00-04:00</td>
      <td>618</td>
      <td>USA: Atlantic Ocean, east of Florida (27 N)</td>
      <td>26.997</td>
      <td>-79.618</td>
      <td>CTD rosette</td>
      <td>Atlantic Ocean</td>
      <td>618</td>
      <td>16S V4-V5 | 18S V9</td>
      <td>GOMECC4_27N_Sta1_Deep_A</td>
    </tr>
    <tr>
      <th>1</th>
      <td>GOMECC4_27N_Sta1_Deep_B_16S</td>
      <td>27N_Sta1</td>
      <td>2021-09-14T11:00-04:00</td>
      <td>618</td>
      <td>USA: Atlantic Ocean, east of Florida (27 N)</td>
      <td>26.997</td>
      <td>-79.618</td>
      <td>CTD rosette</td>
      <td>Atlantic Ocean</td>
      <td>618</td>
      <td>16S V4-V5 | 18S V9</td>
      <td>GOMECC4_27N_Sta1_Deep_B</td>
    </tr>
    <tr>
      <th>2</th>
      <td>GOMECC4_27N_Sta1_Deep_C_16S</td>
      <td>27N_Sta1</td>
      <td>2021-09-14T11:00-04:00</td>
      <td>618</td>
      <td>USA: Atlantic Ocean, east of Florida (27 N)</td>
      <td>26.997</td>
      <td>-79.618</td>
      <td>CTD rosette</td>
      <td>Atlantic Ocean</td>
      <td>618</td>
      <td>16S V4-V5 | 18S V9</td>
      <td>GOMECC4_27N_Sta1_Deep_C</td>
    </tr>
    <tr>
      <th>3</th>
      <td>GOMECC4_27N_Sta1_DCM_A_16S</td>
      <td>27N_Sta1</td>
      <td>2021-09-14T11:00-04:00</td>
      <td>49</td>
      <td>USA: Atlantic Ocean, east of Florida (27 N)</td>
      <td>26.997</td>
      <td>-79.618</td>
      <td>CTD rosette</td>
      <td>Atlantic Ocean</td>
      <td>49</td>
      <td>16S V4-V5 | 18S V9</td>
      <td>GOMECC4_27N_Sta1_DCM_A</td>
    </tr>
    <tr>
      <th>4</th>
      <td>GOMECC4_27N_Sta1_DCM_B_16S</td>
      <td>27N_Sta1</td>
      <td>2021-09-14T11:00-04:00</td>
      <td>49</td>
      <td>USA: Atlantic Ocean, east of Florida (27 N)</td>
      <td>26.997</td>
      <td>-79.618</td>
      <td>CTD rosette</td>
      <td>Atlantic Ocean</td>
      <td>49</td>
      <td>16S V4-V5 | 18S V9</td>
      <td>GOMECC4_27N_Sta1_DCM_B</td>
    </tr>
  </tbody>
</table>
</div>




```python
child_data_18S = event_sample[event_sample['amplicon_sequenced'].str.contains('18S V9')].copy()
child_data_18S['parentEventID'] = child_data_18S['eventID']
child_data_18S['eventID'] = child_data_18S['eventID']+"_18S"
child_data_18S.head()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>eventID</th>
      <th>locationID</th>
      <th>eventDate</th>
      <th>minimumDepthInMeters</th>
      <th>locality</th>
      <th>decimalLatitude</th>
      <th>decimalLongitude</th>
      <th>samplingProtocol</th>
      <th>waterBody</th>
      <th>maximumDepthInMeters</th>
      <th>amplicon_sequenced</th>
      <th>parentEventID</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>GOMECC4_27N_Sta1_Deep_A_18S</td>
      <td>27N_Sta1</td>
      <td>2021-09-14T11:00-04:00</td>
      <td>618</td>
      <td>USA: Atlantic Ocean, east of Florida (27 N)</td>
      <td>26.997</td>
      <td>-79.618</td>
      <td>CTD rosette</td>
      <td>Atlantic Ocean</td>
      <td>618</td>
      <td>16S V4-V5 | 18S V9</td>
      <td>GOMECC4_27N_Sta1_Deep_A</td>
    </tr>
    <tr>
      <th>1</th>
      <td>GOMECC4_27N_Sta1_Deep_B_18S</td>
      <td>27N_Sta1</td>
      <td>2021-09-14T11:00-04:00</td>
      <td>618</td>
      <td>USA: Atlantic Ocean, east of Florida (27 N)</td>
      <td>26.997</td>
      <td>-79.618</td>
      <td>CTD rosette</td>
      <td>Atlantic Ocean</td>
      <td>618</td>
      <td>16S V4-V5 | 18S V9</td>
      <td>GOMECC4_27N_Sta1_Deep_B</td>
    </tr>
    <tr>
      <th>2</th>
      <td>GOMECC4_27N_Sta1_Deep_C_18S</td>
      <td>27N_Sta1</td>
      <td>2021-09-14T11:00-04:00</td>
      <td>618</td>
      <td>USA: Atlantic Ocean, east of Florida (27 N)</td>
      <td>26.997</td>
      <td>-79.618</td>
      <td>CTD rosette</td>
      <td>Atlantic Ocean</td>
      <td>618</td>
      <td>16S V4-V5 | 18S V9</td>
      <td>GOMECC4_27N_Sta1_Deep_C</td>
    </tr>
    <tr>
      <th>3</th>
      <td>GOMECC4_27N_Sta1_DCM_A_18S</td>
      <td>27N_Sta1</td>
      <td>2021-09-14T11:00-04:00</td>
      <td>49</td>
      <td>USA: Atlantic Ocean, east of Florida (27 N)</td>
      <td>26.997</td>
      <td>-79.618</td>
      <td>CTD rosette</td>
      <td>Atlantic Ocean</td>
      <td>49</td>
      <td>16S V4-V5 | 18S V9</td>
      <td>GOMECC4_27N_Sta1_DCM_A</td>
    </tr>
    <tr>
      <th>4</th>
      <td>GOMECC4_27N_Sta1_DCM_B_18S</td>
      <td>27N_Sta1</td>
      <td>2021-09-14T11:00-04:00</td>
      <td>49</td>
      <td>USA: Atlantic Ocean, east of Florida (27 N)</td>
      <td>26.997</td>
      <td>-79.618</td>
      <td>CTD rosette</td>
      <td>Atlantic Ocean</td>
      <td>49</td>
      <td>16S V4-V5 | 18S V9</td>
      <td>GOMECC4_27N_Sta1_DCM_B</td>
    </tr>
  </tbody>
</table>
</div>




```python
# this is your full event file
all_event_data = pd.concat([child_data_16S,child_data_18S],axis=0,ignore_index=True)
```


```python
all_event_data = all_event_data.drop(columns=['amplicon_sequenced'])
```


```python
all_event_data.tail()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>eventID</th>
      <th>locationID</th>
      <th>eventDate</th>
      <th>minimumDepthInMeters</th>
      <th>locality</th>
      <th>decimalLatitude</th>
      <th>decimalLongitude</th>
      <th>samplingProtocol</th>
      <th>waterBody</th>
      <th>maximumDepthInMeters</th>
      <th>parentEventID</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>939</th>
      <td>GOMECC4_CAPECORAL_Sta141_DCM_B_18S</td>
      <td>CAPECORAL_Sta141</td>
      <td>2021-10-20T12:47-04:00</td>
      <td>59</td>
      <td>USA: Gulf of Mexico</td>
      <td>25.574</td>
      <td>-84.843</td>
      <td>CTD rosette</td>
      <td>Mexico, Gulf of</td>
      <td>59</td>
      <td>GOMECC4_CAPECORAL_Sta141_DCM_B</td>
    </tr>
    <tr>
      <th>940</th>
      <td>GOMECC4_CAPECORAL_Sta141_DCM_C_18S</td>
      <td>CAPECORAL_Sta141</td>
      <td>2021-10-20T12:47-04:00</td>
      <td>59</td>
      <td>USA: Gulf of Mexico</td>
      <td>25.574</td>
      <td>-84.843</td>
      <td>CTD rosette</td>
      <td>Mexico, Gulf of</td>
      <td>59</td>
      <td>GOMECC4_CAPECORAL_Sta141_DCM_C</td>
    </tr>
    <tr>
      <th>941</th>
      <td>GOMECC4_CAPECORAL_Sta141_Surface_A_18S</td>
      <td>CAPECORAL_Sta141</td>
      <td>2021-10-20T12:47-04:00</td>
      <td>4</td>
      <td>USA: Gulf of Mexico</td>
      <td>25.574</td>
      <td>-84.843</td>
      <td>CTD rosette</td>
      <td>Mexico, Gulf of</td>
      <td>4</td>
      <td>GOMECC4_CAPECORAL_Sta141_Surface_A</td>
    </tr>
    <tr>
      <th>942</th>
      <td>GOMECC4_CAPECORAL_Sta141_Surface_B_18S</td>
      <td>CAPECORAL_Sta141</td>
      <td>2021-10-20T12:47-04:00</td>
      <td>4</td>
      <td>USA: Gulf of Mexico</td>
      <td>25.574</td>
      <td>-84.843</td>
      <td>CTD rosette</td>
      <td>Mexico, Gulf of</td>
      <td>4</td>
      <td>GOMECC4_CAPECORAL_Sta141_Surface_B</td>
    </tr>
    <tr>
      <th>943</th>
      <td>GOMECC4_CAPECORAL_Sta141_Surface_C_18S</td>
      <td>CAPECORAL_Sta141</td>
      <td>2021-10-20T12:47-04:00</td>
      <td>4</td>
      <td>USA: Gulf of Mexico</td>
      <td>25.574</td>
      <td>-84.843</td>
      <td>CTD rosette</td>
      <td>Mexico, Gulf of</td>
      <td>4</td>
      <td>GOMECC4_CAPECORAL_Sta141_Surface_C</td>
    </tr>
  </tbody>
</table>
</div>



Which terms are still missing from the event info?


```python
for key in event_dict.keys():
    if event_dict[key]['AOML_file'] != 'sample_data':
        print(key,event_dict[key])
```

    geodeticDatum {'AOML_term': 'none', 'AOML_file': 'pipeline', 'DwC_definition': 'The ellipsoid, geodetic datum, or spatial reference system (SRS) upon which the geographic coordinates given in dwc:decimalLatitude and dwc:decimalLongitude are based.'}
    countryCode {'AOML_term': 'none', 'AOML_file': 'pipeline', 'DwC_definition': nan}
    datasetID {'AOML_term': 'project_id_external', 'AOML_file': 'study_data', 'DwC_definition': 'An identifier for the set of data. May be a global unique identifier or an identifier specific to a collection or institution.'}
    eventRemarks {'AOML_term': 'derived: controls_used', 'AOML_file': 'analysis_data', 'DwC_definition': 'Comments or notes about the dwc:Event.'}


countryCode, leave blank because it spans multiple countries. eventRemarks will be added later.


```python
#datasetID
all_event_data['datasetID'] = data['study_data']['project_id_external'].values[0]
```


```python
#geodeticDatum
all_event_data['geodeticDatum'] = "WGS84"

```


```python
all_event_data.head()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>eventID</th>
      <th>locationID</th>
      <th>eventDate</th>
      <th>minimumDepthInMeters</th>
      <th>locality</th>
      <th>decimalLatitude</th>
      <th>decimalLongitude</th>
      <th>samplingProtocol</th>
      <th>waterBody</th>
      <th>maximumDepthInMeters</th>
      <th>parentEventID</th>
      <th>datasetID</th>
      <th>geodeticDatum</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>GOMECC4_27N_Sta1_Deep_A_16S</td>
      <td>27N_Sta1</td>
      <td>2021-09-14T11:00-04:00</td>
      <td>618</td>
      <td>USA: Atlantic Ocean, east of Florida (27 N)</td>
      <td>26.997</td>
      <td>-79.618</td>
      <td>CTD rosette</td>
      <td>Atlantic Ocean</td>
      <td>618</td>
      <td>GOMECC4_27N_Sta1_Deep_A</td>
      <td>noaa-aoml-gomecc4</td>
      <td>WGS84</td>
    </tr>
    <tr>
      <th>1</th>
      <td>GOMECC4_27N_Sta1_Deep_B_16S</td>
      <td>27N_Sta1</td>
      <td>2021-09-14T11:00-04:00</td>
      <td>618</td>
      <td>USA: Atlantic Ocean, east of Florida (27 N)</td>
      <td>26.997</td>
      <td>-79.618</td>
      <td>CTD rosette</td>
      <td>Atlantic Ocean</td>
      <td>618</td>
      <td>GOMECC4_27N_Sta1_Deep_B</td>
      <td>noaa-aoml-gomecc4</td>
      <td>WGS84</td>
    </tr>
    <tr>
      <th>2</th>
      <td>GOMECC4_27N_Sta1_Deep_C_16S</td>
      <td>27N_Sta1</td>
      <td>2021-09-14T11:00-04:00</td>
      <td>618</td>
      <td>USA: Atlantic Ocean, east of Florida (27 N)</td>
      <td>26.997</td>
      <td>-79.618</td>
      <td>CTD rosette</td>
      <td>Atlantic Ocean</td>
      <td>618</td>
      <td>GOMECC4_27N_Sta1_Deep_C</td>
      <td>noaa-aoml-gomecc4</td>
      <td>WGS84</td>
    </tr>
    <tr>
      <th>3</th>
      <td>GOMECC4_27N_Sta1_DCM_A_16S</td>
      <td>27N_Sta1</td>
      <td>2021-09-14T11:00-04:00</td>
      <td>49</td>
      <td>USA: Atlantic Ocean, east of Florida (27 N)</td>
      <td>26.997</td>
      <td>-79.618</td>
      <td>CTD rosette</td>
      <td>Atlantic Ocean</td>
      <td>49</td>
      <td>GOMECC4_27N_Sta1_DCM_A</td>
      <td>noaa-aoml-gomecc4</td>
      <td>WGS84</td>
    </tr>
    <tr>
      <th>4</th>
      <td>GOMECC4_27N_Sta1_DCM_B_16S</td>
      <td>27N_Sta1</td>
      <td>2021-09-14T11:00-04:00</td>
      <td>49</td>
      <td>USA: Atlantic Ocean, east of Florida (27 N)</td>
      <td>26.997</td>
      <td>-79.618</td>
      <td>CTD rosette</td>
      <td>Atlantic Ocean</td>
      <td>49</td>
      <td>GOMECC4_27N_Sta1_DCM_B</td>
      <td>noaa-aoml-gomecc4</td>
      <td>WGS84</td>
    </tr>
  </tbody>
</table>
</div>



### Occurrence file  
Now get the occurrence info from the ASV tables, format it, then merge it with the event info.


```python
# create a dictionary to hold both markers
occ = {}
```

#### 18S


```python
asv_tables['18S V9'].iloc[0:10,0:15]
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>featureid</th>
      <th>sequence</th>
      <th>taxonomy</th>
      <th>Confidence</th>
      <th>GOMECC4_27N_Sta1_DCM_A</th>
      <th>GOMECC4_27N_Sta1_DCM_B</th>
      <th>GOMECC4_27N_Sta1_DCM_C</th>
      <th>GOMECC4_27N_Sta1_Deep_A</th>
      <th>GOMECC4_27N_Sta1_Deep_B</th>
      <th>GOMECC4_27N_Sta1_Deep_C</th>
      <th>GOMECC4_27N_Sta1_Surface_A</th>
      <th>GOMECC4_27N_Sta1_Surface_B</th>
      <th>GOMECC4_27N_Sta4_DCM_A</th>
      <th>GOMECC4_27N_Sta4_DCM_B</th>
      <th>GOMECC4_27N_Sta4_DCM_C</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>36aa75f9b28f5f831c2d631ba65c2bcb</td>
      <td>GCTACTACCGATTGAACGTTTTAGTGAGGTCCTCGGACTGTTTGCCTGGCGGATTACTCTGCCTGGCTGGCGGGAAGACGACCAAACTGTAGCGTTTAGAGGAAGTAAAAGTCGTAACAAGGTTTCC</td>
      <td>Eukaryota;Obazoa;Opisthokonta;Metazoa;Arthropoda;Crustacea;Maxillopoda;Neocalanus;Neocalanus_cristatus;</td>
      <td>0.922099</td>
      <td>1516</td>
      <td>0</td>
      <td>0</td>
      <td>6</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>4257</td>
      <td>2005</td>
      <td>0</td>
      <td>14</td>
    </tr>
    <tr>
      <th>1</th>
      <td>4e38e8ced9070952b314e1880bede1ca</td>
      <td>GCTACTACCGATTGAACGTTTTAGTGAGGTCCTCGGACTGTTTGGTAGTCGGATCACTCTGACTGCCTGGCGGGAAGACGACCAAACTGTAGCGTTTAGAGGAAGTAAAAGTCGTAACAAGGTTTCC</td>
      <td>Eukaryota;Obazoa;Opisthokonta;Metazoa;Arthropoda;Crustacea;Maxillopoda;Clausocalanus;Clausocalanus_furcatus;</td>
      <td>0.999947</td>
      <td>962</td>
      <td>316</td>
      <td>548</td>
      <td>19</td>
      <td>10</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>613</td>
      <td>561</td>
      <td>434</td>
    </tr>
    <tr>
      <th>2</th>
      <td>5d4df37251121c08397c6fbc27b06175</td>
      <td>GCTACTACCGATTGAGTGTTTTAGTGAGGTCCTCGGATTGCTTTCCTGGCGGTTAACGCTGCCTAGTTGGCGAAAAGACGACCAAACTGTAGCACTTAGAGGAAGTAAAAGTCGTAACAAGGTTTCC</td>
      <td>Eukaryota;Obazoa;Opisthokonta;Metazoa;Arthropoda;Crustacea;Maxillopoda;Sinocalanus;Sinocalanus_sinensis;</td>
      <td>0.992300</td>
      <td>0</td>
      <td>4</td>
      <td>0</td>
      <td>12</td>
      <td>5</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>9</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>3</th>
      <td>f863f671a575c6ab587e8de0190d3335</td>
      <td>GCTACTACCGATTGAACATTTTAGTGAGGTCCTCGGACTGTGAGCCAGGCGGGTCGCCCTGCCTGGTCTACGGGAAGACGACCAAACTGTAGTGTTTAGAGGAAGTAAAAGTCGTAACAAGGTTTCC</td>
      <td>Eukaryota;Obazoa;Opisthokonta;Metazoa;Arthropoda;Crustacea;Maxillopoda;Paracalanus;Paracalanus_parvus;</td>
      <td>0.998393</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>5</td>
    </tr>
    <tr>
      <th>4</th>
      <td>2a31e5c01634165da99e7381279baa75</td>
      <td>GCTACTACCGATTGGACGTTTTAGTGAGACATTTGGACTGGGTTAAGATAGTCGCAAGACTACCTTTTCTCCGGAAAGACTTTCAAACTTGAGCGTCTAGAGGAAGTAAAAGTCGTAACAAGGTTTCC</td>
      <td>Eukaryota;Obazoa;Opisthokonta;Metazoa;Arthropoda;Crustacea;Maxillopoda;Acrocalanus;Acrocalanus_sp.;</td>
      <td>0.779948</td>
      <td>1164</td>
      <td>2272</td>
      <td>2208</td>
      <td>2</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>5</th>
      <td>ecee60339b2fb88ea6d1c8d18054bed4</td>
      <td>GCTCCTACCGATTGAGTGATCCGGTGAATAATTCGGACTGCAGCAGTGTTCAGTTCCTGAACGTTGCAGCGGAAAGTTTAGTGAACCTTATCACTTAGAGGAAGGAGAAGTCGTAACAAGGTTTCC</td>
      <td>Eukaryota;TSAR;Alveolata;Dinoflagellata;Dinophyceae</td>
      <td>0.999931</td>
      <td>287</td>
      <td>414</td>
      <td>335</td>
      <td>195</td>
      <td>228</td>
      <td>298</td>
      <td>252</td>
      <td>349</td>
      <td>175</td>
      <td>102</td>
      <td>216</td>
    </tr>
    <tr>
      <th>6</th>
      <td>d70494a723d85d66aa88d2d8a975aeec</td>
      <td>GCTACTACCGATTGAATGGTTCCGTGAATTCTTGAGATCGGCGCGGGAACAACTGGCAACGGTTGATCCCGATTGCTGAGAACTTGTGTAAACGCGATCATTTAGAGGAAGTAAAAGTCGTAACAAGGTTTCC</td>
      <td>Eukaryota;Obazoa;Opisthokonta</td>
      <td>0.992451</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>4</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>5</td>
      <td>0</td>
      <td>22</td>
    </tr>
    <tr>
      <th>7</th>
      <td>fa1f1a97dd4ae7c826009186bad26384</td>
      <td>GCTCCTACCGATTGAGTGATCCGGTGAATAATTCGGACTGCAGCAATGTTTGGATCCCGAACGTTGCAGCGGAAAGTTTAGTGAACCTTATCACTTAGAGGAAGGAGAAGTCGTAACAAGGTTTCC</td>
      <td>Eukaryota;TSAR;Alveolata;Dinoflagellata;Dinophyceae;Gymnodiniales;Gymnodiniaceae</td>
      <td>0.986908</td>
      <td>250</td>
      <td>323</td>
      <td>194</td>
      <td>51</td>
      <td>59</td>
      <td>55</td>
      <td>222</td>
      <td>250</td>
      <td>230</td>
      <td>163</td>
      <td>214</td>
    </tr>
    <tr>
      <th>8</th>
      <td>bbaaf7bb4e71c80de970677779e3bf3a</td>
      <td>GCTACTACCGATTGAATGGTTTAGTGAGATCTTCGGATTGGCACAATCGCGGCCTAACGGAAGTGATGGTGCCGAAAAGTTGCTCAAACTTGATCATTTAGAGGAAGTAAAAGTCGTAACAAGGTTTCC</td>
      <td>Eukaryota;Obazoa;Opisthokonta;Metazoa;Cnidaria;Cnidaria_X;Hydrozoa;Sulculeolaria;Sulculeolaria_quadrivalvis;</td>
      <td>0.864777</td>
      <td>212</td>
      <td>50</td>
      <td>237</td>
      <td>552</td>
      <td>1278</td>
      <td>480</td>
      <td>0</td>
      <td>0</td>
      <td>26</td>
      <td>24</td>
      <td>21</td>
    </tr>
    <tr>
      <th>9</th>
      <td>7a8324bb4448b65f7adc73d70e5901da</td>
      <td>GCTACTACCGATTGAACGTTTTAGTGAGGTATTTGGACTGGGCCTTTGGAGGATTCGTTCTCCAATGTTGCTCGGGAAGACTCCCAAACTTGAGCGTTTAGAGGAAGTAAAAGTCGTAACAAGGTTTCC</td>
      <td>Eukaryota;Obazoa;Opisthokonta;Metazoa;Arthropoda;Crustacea;Maxillopoda;Delibus;Delibus_sp.;</td>
      <td>0.992088</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>15</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
  </tbody>
</table>
</div>



##### format taxonomy

How to automate this? Everyone's taxonomy might be different?


```python
#18S 
taxa_ranks_18S = ['domain','supergroup','division','subdivision','class','order','family','genus','species']

asv_tables['18S V9'][['domain','supergroup','division','subdivision','class','order','family','genus','species']] = ["","","","","","","","",""]
for index, row in asv_tables['18S V9'].iterrows():
    taxa = row['taxonomy'].split(";")
    for i in range(0,len(taxa)):
        if i < len(taxa_ranks_18S):
            asv_tables['18S V9'].loc[index,taxa_ranks_18S[i]] = taxa[i]

    
```


```python
# replace None with NA
asv_tables['18S V9'] = asv_tables['18S V9'].fillna(value=np.nan)
## Replace 'unknown', 'unassigned', etc. in species and taxonomy columns with NaN

asv_tables['18S V9'][taxa_ranks_18S] = asv_tables['18S V9'][taxa_ranks_18S].replace({'unassigned':np.nan,
                            'Unassigned':np.nan,
                              's_':np.nan,
                              'g_':np.nan,
                              'unknown':np.nan,
                              'no_hit':np.nan,
                               '':np.nan})
asv_tables['18S V9'].iloc[0:10,[0,1,2,3,4,5,6,-9,-8,-7,-6,-5,-4,-3,-2,-1]]
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>featureid</th>
      <th>sequence</th>
      <th>taxonomy</th>
      <th>Confidence</th>
      <th>GOMECC4_27N_Sta1_DCM_A</th>
      <th>GOMECC4_27N_Sta1_DCM_B</th>
      <th>GOMECC4_27N_Sta1_DCM_C</th>
      <th>domain</th>
      <th>supergroup</th>
      <th>division</th>
      <th>subdivision</th>
      <th>class</th>
      <th>order</th>
      <th>family</th>
      <th>genus</th>
      <th>species</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>36aa75f9b28f5f831c2d631ba65c2bcb</td>
      <td>GCTACTACCGATTGAACGTTTTAGTGAGGTCCTCGGACTGTTTGCCTGGCGGATTACTCTGCCTGGCTGGCGGGAAGACGACCAAACTGTAGCGTTTAGAGGAAGTAAAAGTCGTAACAAGGTTTCC</td>
      <td>Eukaryota;Obazoa;Opisthokonta;Metazoa;Arthropoda;Crustacea;Maxillopoda;Neocalanus;Neocalanus_cristatus;</td>
      <td>0.922099</td>
      <td>1516</td>
      <td>0</td>
      <td>0</td>
      <td>Eukaryota</td>
      <td>Obazoa</td>
      <td>Opisthokonta</td>
      <td>Metazoa</td>
      <td>Arthropoda</td>
      <td>Crustacea</td>
      <td>Maxillopoda</td>
      <td>Neocalanus</td>
      <td>Neocalanus cristatus</td>
    </tr>
    <tr>
      <th>1</th>
      <td>4e38e8ced9070952b314e1880bede1ca</td>
      <td>GCTACTACCGATTGAACGTTTTAGTGAGGTCCTCGGACTGTTTGGTAGTCGGATCACTCTGACTGCCTGGCGGGAAGACGACCAAACTGTAGCGTTTAGAGGAAGTAAAAGTCGTAACAAGGTTTCC</td>
      <td>Eukaryota;Obazoa;Opisthokonta;Metazoa;Arthropoda;Crustacea;Maxillopoda;Clausocalanus;Clausocalanus_furcatus;</td>
      <td>0.999947</td>
      <td>962</td>
      <td>316</td>
      <td>548</td>
      <td>Eukaryota</td>
      <td>Obazoa</td>
      <td>Opisthokonta</td>
      <td>Metazoa</td>
      <td>Arthropoda</td>
      <td>Crustacea</td>
      <td>Maxillopoda</td>
      <td>Clausocalanus</td>
      <td>Clausocalanus furcatus</td>
    </tr>
    <tr>
      <th>2</th>
      <td>5d4df37251121c08397c6fbc27b06175</td>
      <td>GCTACTACCGATTGAGTGTTTTAGTGAGGTCCTCGGATTGCTTTCCTGGCGGTTAACGCTGCCTAGTTGGCGAAAAGACGACCAAACTGTAGCACTTAGAGGAAGTAAAAGTCGTAACAAGGTTTCC</td>
      <td>Eukaryota;Obazoa;Opisthokonta;Metazoa;Arthropoda;Crustacea;Maxillopoda;Sinocalanus;Sinocalanus_sinensis;</td>
      <td>0.992300</td>
      <td>0</td>
      <td>4</td>
      <td>0</td>
      <td>Eukaryota</td>
      <td>Obazoa</td>
      <td>Opisthokonta</td>
      <td>Metazoa</td>
      <td>Arthropoda</td>
      <td>Crustacea</td>
      <td>Maxillopoda</td>
      <td>Sinocalanus</td>
      <td>Sinocalanus sinensis</td>
    </tr>
    <tr>
      <th>3</th>
      <td>f863f671a575c6ab587e8de0190d3335</td>
      <td>GCTACTACCGATTGAACATTTTAGTGAGGTCCTCGGACTGTGAGCCAGGCGGGTCGCCCTGCCTGGTCTACGGGAAGACGACCAAACTGTAGTGTTTAGAGGAAGTAAAAGTCGTAACAAGGTTTCC</td>
      <td>Eukaryota;Obazoa;Opisthokonta;Metazoa;Arthropoda;Crustacea;Maxillopoda;Paracalanus;Paracalanus_parvus;</td>
      <td>0.998393</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>Eukaryota</td>
      <td>Obazoa</td>
      <td>Opisthokonta</td>
      <td>Metazoa</td>
      <td>Arthropoda</td>
      <td>Crustacea</td>
      <td>Maxillopoda</td>
      <td>Paracalanus</td>
      <td>Paracalanus parvus</td>
    </tr>
    <tr>
      <th>4</th>
      <td>2a31e5c01634165da99e7381279baa75</td>
      <td>GCTACTACCGATTGGACGTTTTAGTGAGACATTTGGACTGGGTTAAGATAGTCGCAAGACTACCTTTTCTCCGGAAAGACTTTCAAACTTGAGCGTCTAGAGGAAGTAAAAGTCGTAACAAGGTTTCC</td>
      <td>Eukaryota;Obazoa;Opisthokonta;Metazoa;Arthropoda;Crustacea;Maxillopoda;Acrocalanus;Acrocalanus_sp.;</td>
      <td>0.779948</td>
      <td>1164</td>
      <td>2272</td>
      <td>2208</td>
      <td>Eukaryota</td>
      <td>Obazoa</td>
      <td>Opisthokonta</td>
      <td>Metazoa</td>
      <td>Arthropoda</td>
      <td>Crustacea</td>
      <td>Maxillopoda</td>
      <td>Acrocalanus</td>
      <td>Acrocalanus</td>
    </tr>
    <tr>
      <th>5</th>
      <td>ecee60339b2fb88ea6d1c8d18054bed4</td>
      <td>GCTCCTACCGATTGAGTGATCCGGTGAATAATTCGGACTGCAGCAGTGTTCAGTTCCTGAACGTTGCAGCGGAAAGTTTAGTGAACCTTATCACTTAGAGGAAGGAGAAGTCGTAACAAGGTTTCC</td>
      <td>Eukaryota;TSAR;Alveolata;Dinoflagellata;Dinophyceae</td>
      <td>0.999931</td>
      <td>287</td>
      <td>414</td>
      <td>335</td>
      <td>Eukaryota</td>
      <td>TSAR</td>
      <td>Alveolata</td>
      <td>Dinoflagellata</td>
      <td>Dinophyceae</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>6</th>
      <td>d70494a723d85d66aa88d2d8a975aeec</td>
      <td>GCTACTACCGATTGAATGGTTCCGTGAATTCTTGAGATCGGCGCGGGAACAACTGGCAACGGTTGATCCCGATTGCTGAGAACTTGTGTAAACGCGATCATTTAGAGGAAGTAAAAGTCGTAACAAGGTTTCC</td>
      <td>Eukaryota;Obazoa;Opisthokonta</td>
      <td>0.992451</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>Eukaryota</td>
      <td>Obazoa</td>
      <td>Opisthokonta</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>7</th>
      <td>fa1f1a97dd4ae7c826009186bad26384</td>
      <td>GCTCCTACCGATTGAGTGATCCGGTGAATAATTCGGACTGCAGCAATGTTTGGATCCCGAACGTTGCAGCGGAAAGTTTAGTGAACCTTATCACTTAGAGGAAGGAGAAGTCGTAACAAGGTTTCC</td>
      <td>Eukaryota;TSAR;Alveolata;Dinoflagellata;Dinophyceae;Gymnodiniales;Gymnodiniaceae</td>
      <td>0.986908</td>
      <td>250</td>
      <td>323</td>
      <td>194</td>
      <td>Eukaryota</td>
      <td>TSAR</td>
      <td>Alveolata</td>
      <td>Dinoflagellata</td>
      <td>Dinophyceae</td>
      <td>Gymnodiniales</td>
      <td>Gymnodiniaceae</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>8</th>
      <td>bbaaf7bb4e71c80de970677779e3bf3a</td>
      <td>GCTACTACCGATTGAATGGTTTAGTGAGATCTTCGGATTGGCACAATCGCGGCCTAACGGAAGTGATGGTGCCGAAAAGTTGCTCAAACTTGATCATTTAGAGGAAGTAAAAGTCGTAACAAGGTTTCC</td>
      <td>Eukaryota;Obazoa;Opisthokonta;Metazoa;Cnidaria;Cnidaria_X;Hydrozoa;Sulculeolaria;Sulculeolaria_quadrivalvis;</td>
      <td>0.864777</td>
      <td>212</td>
      <td>50</td>
      <td>237</td>
      <td>Eukaryota</td>
      <td>Obazoa</td>
      <td>Opisthokonta</td>
      <td>Metazoa</td>
      <td>Cnidaria</td>
      <td>Cnidaria X</td>
      <td>Hydrozoa</td>
      <td>Sulculeolaria</td>
      <td>Sulculeolaria quadrivalvis</td>
    </tr>
    <tr>
      <th>9</th>
      <td>7a8324bb4448b65f7adc73d70e5901da</td>
      <td>GCTACTACCGATTGAACGTTTTAGTGAGGTATTTGGACTGGGCCTTTGGAGGATTCGTTCTCCAATGTTGCTCGGGAAGACTCCCAAACTTGAGCGTTTAGAGGAAGTAAAAGTCGTAACAAGGTTTCC</td>
      <td>Eukaryota;Obazoa;Opisthokonta;Metazoa;Arthropoda;Crustacea;Maxillopoda;Delibus;Delibus_sp.;</td>
      <td>0.992088</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>Eukaryota</td>
      <td>Obazoa</td>
      <td>Opisthokonta</td>
      <td>Metazoa</td>
      <td>Arthropoda</td>
      <td>Crustacea</td>
      <td>Maxillopoda</td>
      <td>Delibus</td>
      <td>Delibus</td>
    </tr>
  </tbody>
</table>
</div>




```python
# replace _,- with space, remove sp. 

asv_tables['18S V9'][taxa_ranks_18S] = asv_tables['18S V9'][taxa_ranks_18S].replace('_',' ',regex=True)
asv_tables['18S V9'][taxa_ranks_18S] = asv_tables['18S V9'][taxa_ranks_18S].replace(' sp\.','',regex=True)
asv_tables['18S V9'][taxa_ranks_18S] = asv_tables['18S V9'][taxa_ranks_18S].replace(' spp\.','',regex=True)
asv_tables['18S V9'][taxa_ranks_18S] = asv_tables['18S V9'][taxa_ranks_18S].replace('-',' ',regex=True)
asv_tables['18S V9'][taxa_ranks_18S] = asv_tables['18S V9'][taxa_ranks_18S].replace('\/',' ',regex=True)
```


```python
asv_tables['18S V9'].shape

```




    (24067, 485)




```python
occ['18S V9'] = pd.melt(asv_tables['18S V9'],id_vars=['featureid','sequence','taxonomy','Confidence','domain','supergroup','division','subdivision','class','order','family','genus','species'],
               var_name='eventID',value_name='organismQuantity')
```


```python
occ['18S V9'].shape
```




    (11359624, 15)




```python
## Drop records where organismQuantity = 0 (absences are not meaningful for this data set)

occ['18S V9'] = occ['18S V9'][occ['18S V9']['organismQuantity'] > 0]
print(occ['18S V9'].shape)
```

    (146232, 15)



```python
occ['18S V9'].head()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>featureid</th>
      <th>sequence</th>
      <th>taxonomy</th>
      <th>Confidence</th>
      <th>domain</th>
      <th>supergroup</th>
      <th>division</th>
      <th>subdivision</th>
      <th>class</th>
      <th>order</th>
      <th>family</th>
      <th>genus</th>
      <th>species</th>
      <th>eventID</th>
      <th>organismQuantity</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>36aa75f9b28f5f831c2d631ba65c2bcb</td>
      <td>GCTACTACCGATTGAACGTTTTAGTGAGGTCCTCGGACTGTTTGCCTGGCGGATTACTCTGCCTGGCTGGCGGGAAGACGACCAAACTGTAGCGTTTAGAGGAAGTAAAAGTCGTAACAAGGTTTCC</td>
      <td>Eukaryota;Obazoa;Opisthokonta;Metazoa;Arthropoda;Crustacea;Maxillopoda;Neocalanus;Neocalanus_cristatus;</td>
      <td>0.922099</td>
      <td>Eukaryota</td>
      <td>Obazoa</td>
      <td>Opisthokonta</td>
      <td>Metazoa</td>
      <td>Arthropoda</td>
      <td>Crustacea</td>
      <td>Maxillopoda</td>
      <td>Neocalanus</td>
      <td>Neocalanus cristatus</td>
      <td>GOMECC4_27N_Sta1_DCM_A</td>
      <td>1516</td>
    </tr>
    <tr>
      <th>1</th>
      <td>4e38e8ced9070952b314e1880bede1ca</td>
      <td>GCTACTACCGATTGAACGTTTTAGTGAGGTCCTCGGACTGTTTGGTAGTCGGATCACTCTGACTGCCTGGCGGGAAGACGACCAAACTGTAGCGTTTAGAGGAAGTAAAAGTCGTAACAAGGTTTCC</td>
      <td>Eukaryota;Obazoa;Opisthokonta;Metazoa;Arthropoda;Crustacea;Maxillopoda;Clausocalanus;Clausocalanus_furcatus;</td>
      <td>0.999947</td>
      <td>Eukaryota</td>
      <td>Obazoa</td>
      <td>Opisthokonta</td>
      <td>Metazoa</td>
      <td>Arthropoda</td>
      <td>Crustacea</td>
      <td>Maxillopoda</td>
      <td>Clausocalanus</td>
      <td>Clausocalanus furcatus</td>
      <td>GOMECC4_27N_Sta1_DCM_A</td>
      <td>962</td>
    </tr>
    <tr>
      <th>4</th>
      <td>2a31e5c01634165da99e7381279baa75</td>
      <td>GCTACTACCGATTGGACGTTTTAGTGAGACATTTGGACTGGGTTAAGATAGTCGCAAGACTACCTTTTCTCCGGAAAGACTTTCAAACTTGAGCGTCTAGAGGAAGTAAAAGTCGTAACAAGGTTTCC</td>
      <td>Eukaryota;Obazoa;Opisthokonta;Metazoa;Arthropoda;Crustacea;Maxillopoda;Acrocalanus;Acrocalanus_sp.;</td>
      <td>0.779948</td>
      <td>Eukaryota</td>
      <td>Obazoa</td>
      <td>Opisthokonta</td>
      <td>Metazoa</td>
      <td>Arthropoda</td>
      <td>Crustacea</td>
      <td>Maxillopoda</td>
      <td>Acrocalanus</td>
      <td>Acrocalanus</td>
      <td>GOMECC4_27N_Sta1_DCM_A</td>
      <td>1164</td>
    </tr>
    <tr>
      <th>5</th>
      <td>ecee60339b2fb88ea6d1c8d18054bed4</td>
      <td>GCTCCTACCGATTGAGTGATCCGGTGAATAATTCGGACTGCAGCAGTGTTCAGTTCCTGAACGTTGCAGCGGAAAGTTTAGTGAACCTTATCACTTAGAGGAAGGAGAAGTCGTAACAAGGTTTCC</td>
      <td>Eukaryota;TSAR;Alveolata;Dinoflagellata;Dinophyceae</td>
      <td>0.999931</td>
      <td>Eukaryota</td>
      <td>TSAR</td>
      <td>Alveolata</td>
      <td>Dinoflagellata</td>
      <td>Dinophyceae</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>GOMECC4_27N_Sta1_DCM_A</td>
      <td>287</td>
    </tr>
    <tr>
      <th>7</th>
      <td>fa1f1a97dd4ae7c826009186bad26384</td>
      <td>GCTCCTACCGATTGAGTGATCCGGTGAATAATTCGGACTGCAGCAATGTTTGGATCCCGAACGTTGCAGCGGAAAGTTTAGTGAACCTTATCACTTAGAGGAAGGAGAAGTCGTAACAAGGTTTCC</td>
      <td>Eukaryota;TSAR;Alveolata;Dinoflagellata;Dinophyceae;Gymnodiniales;Gymnodiniaceae</td>
      <td>0.986908</td>
      <td>Eukaryota</td>
      <td>TSAR</td>
      <td>Alveolata</td>
      <td>Dinoflagellata</td>
      <td>Dinophyceae</td>
      <td>Gymnodiniales</td>
      <td>Gymnodiniaceae</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>GOMECC4_27N_Sta1_DCM_A</td>
      <td>250</td>
    </tr>
  </tbody>
</table>
</div>



Add occurenceID


```python
## Create an occurrenceID that will uniquely identify each ASV observed within a water sample

occ['18S V9']['occurrenceID'] = occ['18S V9']['featureid']
occ['18S V9']['occurrenceID'] = occ['18S V9']['eventID'] + '_occ' + occ['18S V9']['occurrenceID'].astype(str)
```


```python
occ['18S V9'].head()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>featureid</th>
      <th>sequence</th>
      <th>taxonomy</th>
      <th>Confidence</th>
      <th>domain</th>
      <th>supergroup</th>
      <th>division</th>
      <th>subdivision</th>
      <th>class</th>
      <th>order</th>
      <th>family</th>
      <th>genus</th>
      <th>species</th>
      <th>eventID</th>
      <th>organismQuantity</th>
      <th>occurrenceID</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>36aa75f9b28f5f831c2d631ba65c2bcb</td>
      <td>GCTACTACCGATTGAACGTTTTAGTGAGGTCCTCGGACTGTTTGCCTGGCGGATTACTCTGCCTGGCTGGCGGGAAGACGACCAAACTGTAGCGTTTAGAGGAAGTAAAAGTCGTAACAAGGTTTCC</td>
      <td>Eukaryota;Obazoa;Opisthokonta;Metazoa;Arthropoda;Crustacea;Maxillopoda;Neocalanus;Neocalanus_cristatus;</td>
      <td>0.922099</td>
      <td>Eukaryota</td>
      <td>Obazoa</td>
      <td>Opisthokonta</td>
      <td>Metazoa</td>
      <td>Arthropoda</td>
      <td>Crustacea</td>
      <td>Maxillopoda</td>
      <td>Neocalanus</td>
      <td>Neocalanus cristatus</td>
      <td>GOMECC4_27N_Sta1_DCM_A</td>
      <td>1516</td>
      <td>GOMECC4_27N_Sta1_DCM_A_occ36aa75f9b28f5f831c2d631ba65c2bcb</td>
    </tr>
    <tr>
      <th>1</th>
      <td>4e38e8ced9070952b314e1880bede1ca</td>
      <td>GCTACTACCGATTGAACGTTTTAGTGAGGTCCTCGGACTGTTTGGTAGTCGGATCACTCTGACTGCCTGGCGGGAAGACGACCAAACTGTAGCGTTTAGAGGAAGTAAAAGTCGTAACAAGGTTTCC</td>
      <td>Eukaryota;Obazoa;Opisthokonta;Metazoa;Arthropoda;Crustacea;Maxillopoda;Clausocalanus;Clausocalanus_furcatus;</td>
      <td>0.999947</td>
      <td>Eukaryota</td>
      <td>Obazoa</td>
      <td>Opisthokonta</td>
      <td>Metazoa</td>
      <td>Arthropoda</td>
      <td>Crustacea</td>
      <td>Maxillopoda</td>
      <td>Clausocalanus</td>
      <td>Clausocalanus furcatus</td>
      <td>GOMECC4_27N_Sta1_DCM_A</td>
      <td>962</td>
      <td>GOMECC4_27N_Sta1_DCM_A_occ4e38e8ced9070952b314e1880bede1ca</td>
    </tr>
    <tr>
      <th>4</th>
      <td>2a31e5c01634165da99e7381279baa75</td>
      <td>GCTACTACCGATTGGACGTTTTAGTGAGACATTTGGACTGGGTTAAGATAGTCGCAAGACTACCTTTTCTCCGGAAAGACTTTCAAACTTGAGCGTCTAGAGGAAGTAAAAGTCGTAACAAGGTTTCC</td>
      <td>Eukaryota;Obazoa;Opisthokonta;Metazoa;Arthropoda;Crustacea;Maxillopoda;Acrocalanus;Acrocalanus_sp.;</td>
      <td>0.779948</td>
      <td>Eukaryota</td>
      <td>Obazoa</td>
      <td>Opisthokonta</td>
      <td>Metazoa</td>
      <td>Arthropoda</td>
      <td>Crustacea</td>
      <td>Maxillopoda</td>
      <td>Acrocalanus</td>
      <td>Acrocalanus</td>
      <td>GOMECC4_27N_Sta1_DCM_A</td>
      <td>1164</td>
      <td>GOMECC4_27N_Sta1_DCM_A_occ2a31e5c01634165da99e7381279baa75</td>
    </tr>
    <tr>
      <th>5</th>
      <td>ecee60339b2fb88ea6d1c8d18054bed4</td>
      <td>GCTCCTACCGATTGAGTGATCCGGTGAATAATTCGGACTGCAGCAGTGTTCAGTTCCTGAACGTTGCAGCGGAAAGTTTAGTGAACCTTATCACTTAGAGGAAGGAGAAGTCGTAACAAGGTTTCC</td>
      <td>Eukaryota;TSAR;Alveolata;Dinoflagellata;Dinophyceae</td>
      <td>0.999931</td>
      <td>Eukaryota</td>
      <td>TSAR</td>
      <td>Alveolata</td>
      <td>Dinoflagellata</td>
      <td>Dinophyceae</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>GOMECC4_27N_Sta1_DCM_A</td>
      <td>287</td>
      <td>GOMECC4_27N_Sta1_DCM_A_occecee60339b2fb88ea6d1c8d18054bed4</td>
    </tr>
    <tr>
      <th>7</th>
      <td>fa1f1a97dd4ae7c826009186bad26384</td>
      <td>GCTCCTACCGATTGAGTGATCCGGTGAATAATTCGGACTGCAGCAATGTTTGGATCCCGAACGTTGCAGCGGAAAGTTTAGTGAACCTTATCACTTAGAGGAAGGAGAAGTCGTAACAAGGTTTCC</td>
      <td>Eukaryota;TSAR;Alveolata;Dinoflagellata;Dinophyceae;Gymnodiniales;Gymnodiniaceae</td>
      <td>0.986908</td>
      <td>Eukaryota</td>
      <td>TSAR</td>
      <td>Alveolata</td>
      <td>Dinoflagellata</td>
      <td>Dinophyceae</td>
      <td>Gymnodiniales</td>
      <td>Gymnodiniaceae</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>GOMECC4_27N_Sta1_DCM_A</td>
      <td>250</td>
      <td>GOMECC4_27N_Sta1_DCM_A_occfa1f1a97dd4ae7c826009186bad26384</td>
    </tr>
  </tbody>
</table>
</div>



#### 16S

##### 1st, format ASV file


```python
asv_tables['16S V4-V5'].iloc[0:10,0:20]
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>featureid</th>
      <th>sequence</th>
      <th>taxonomy</th>
      <th>Confidence</th>
      <th>GOMECC4_27N_Sta1_DCM_A</th>
      <th>GOMECC4_27N_Sta1_DCM_B</th>
      <th>GOMECC4_27N_Sta1_DCM_C</th>
      <th>GOMECC4_27N_Sta1_Deep_A</th>
      <th>GOMECC4_27N_Sta1_Deep_B</th>
      <th>GOMECC4_27N_Sta1_Deep_C</th>
      <th>GOMECC4_27N_Sta1_Surface_A</th>
      <th>GOMECC4_27N_Sta1_Surface_B</th>
      <th>GOMECC4_27N_Sta4_DCM_A</th>
      <th>GOMECC4_27N_Sta4_DCM_B</th>
      <th>GOMECC4_27N_Sta4_DCM_C</th>
      <th>GOMECC4_27N_Sta4_Deep_A</th>
      <th>GOMECC4_27N_Sta4_Deep_B</th>
      <th>GOMECC4_27N_Sta4_Deep_C</th>
      <th>GOMECC4_27N_Sta4_Surface_A</th>
      <th>GOMECC4_27N_Sta4_Surface_B</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>00006f0784f7dbb2f162408abb6da629</td>
      <td>TACGGAGGGTGCGAGCGTTAATCGGAATTACTGGGCGTAAAGCGCATGCAGGTGGTTTGTTAAGTCAGATGTGAAAGCCCGGGGCTCAACCTCGGAATTGCATTTGAAACTGGCAGACTAGAGTACTGTAGAGGGGGGTAGAATTT...</td>
      <td>d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Vibrionales; f__Vibrionaceae; g__Vibrio</td>
      <td>0.978926</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>25</td>
    </tr>
    <tr>
      <th>1</th>
      <td>000094731d4984ed41435a1bf65b7ef2</td>
      <td>TACAGAGAGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCGCGTAGGTGGGTATTTAAGTCGGATGTGAAATCCCCGGGCTTAACCTGGGAACTGCATCCGAAACTATTTAACTAGAGTATGGGAGAGGTAAGTAGAATTT...</td>
      <td>d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__HOC36; f__HOC36; g__HOC36; s__Candidatus_Thioglobus</td>
      <td>0.881698</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>2</th>
      <td>0001a3c11fcef1b1b8f4c72942efbbac</td>
      <td>TACGAAGGGGGCGAGCGTTGTTCGGAATTACTGGGCGTAAAGGGCGCGTAGGCGGTCTTCTAAGTTAGGCGTGAAAGCCCCGGGCTCAACCTGGGAACTGCGCTTAATACTGGAAGACTAGAAAACGGAAGAGGGTAGTGGAATTC...</td>
      <td>d__Bacteria; p__Cyanobacteria; c__Cyanobacteriia; o__Synechococcales; f__Cyanobiaceae; g__Cyanobium_PCC-6307</td>
      <td>0.762793</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>3</th>
      <td>0001ceef5162e6d689ef30418cfcc164</td>
      <td>TACAGAGGGTGCAAGCGTTGTTCGGAATCATTGGGCGTAAAGCGCGCGTAGGCGGCCAAATAAGTCTGATGTGAAGGCCCAGGGCTCAACCCTGGAAGTGCATCGGAAACTGTTTGGCTCGAGTCCCGGAGGGGGTGGTGGAATTC...</td>
      <td>d__Bacteria; p__Myxococcota; c__Myxococcia; o__Myxococcales; f__Myxococcaceae; g__P3OB-42; s__uncultured_bacterium</td>
      <td>0.997619</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>4</th>
      <td>000235534662df05bb30219a4b978dac</td>
      <td>TACGGAAGGTCCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCGCGTAGGTGGTTTTTTAAGTTGGATGTGAAAGCCCTGGGCTCAACCTAGGAACTGCATCCAAAACTAGATGACTAGAGTACGAAAGAGGGAAGTAGAATTC...</td>
      <td>d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__SAR86_clade; f__SAR86_clade; g__SAR86_clade</td>
      <td>0.999961</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>5</th>
      <td>0003aeafc4bc0522877d4804829e65b7</td>
      <td>CACCGGCATCTCGAGTGGTATCCACTTTTATTGGGCCTAAAGCATCCGTAGCCTGTTCTGTAAGTTTTCGGTTAAATCCATAAGCTCAACTTATGGGCTGCCGAAAATACTGCAGAACTAGGGAGTGGGAGAGGTAGACGGTACTC...</td>
      <td>d__Archaea; p__Crenarchaeota; c__Nitrososphaeria; o__Nitrosopumilales; f__Nitrosopumilaceae; g__Candidatus_Nitrosopelagicus; s__marine_metagenome</td>
      <td>0.779476</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>6</th>
      <td>0003b46ae196127658c07aeb11b36b1a</td>
      <td>TACGGAGGGTCCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCGCGTAGGCGGTTTAGCAAGTTGAATGTGAAAGCCCTGGGCTCAACCTAGGAACTGCATTCAAAACTACTAAGCTAGAGTACGAGAGAGGAGAGTAGAATTT...</td>
      <td>d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Thiomicrospirales; f__Thioglobaceae; g__SUP05_cluster</td>
      <td>0.810269</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>7</th>
      <td>00065c7f5701f8db77fc2c50a1204c71</td>
      <td>TACGGGAGTGGCAAGCGTTATCCGGAATTATTGGGCGTAAAGCGTTTGTAGGTGGAAAAATAAGTCTATTGTTAAATCCAGAAGCTTAACTTCTGTCAAGCGATATGAAACTATTCTTCTTGAGAATGGTAGGGGTAGAAGGAATT...</td>
      <td>d__Bacteria; p__Cyanobacteria; c__Cyanobacteriia; o__Chloroplast; f__Chloroplast; g__Chloroplast; s__Prasinoderma_coloniale</td>
      <td>1.000000</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>8</th>
      <td>0006da4e1ff162826badd8bdcfaf9dfe</td>
      <td>GACGGAGGATGCAAGTGTTATCCGGAATTATTGGGCGTAAAGCGTTTGTAGGTGGAGAAATAAGCCTATTGTTAAATCCAGGAGCTTAACTTCTGTCCAGCGATATGAAACTATTTTTCTTGAGGGTGGTAGGGGTAGAAGGAATT...</td>
      <td>d__Bacteria; p__Cyanobacteria; c__Cyanobacteriia; o__Chloroplast; f__Chloroplast; g__Chloroplast; s__Prasinoderma_coloniale</td>
      <td>1.000000</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>9</th>
      <td>0006fe6033cca4da30ea5ce9cba446f0</td>
      <td>TACGAATGCTGCAAGCGTAGTTCGGAATCACTGGGCATAAAGAGCACGTAGGCGGCCTATTAAGTCAGCTGTGAAATCCCTCGGCTTAACCGAGGAACTGCAGCTGATACTGATAGGCTTGAGTACGGGAGGGGAGAGCGGAATTC...</td>
      <td>d__Bacteria; p__Planctomycetota; c__Pla3_lineage; o__Pla3_lineage; f__Pla3_lineage; g__Pla3_lineage; s__uncultured_Planctomycetaceae</td>
      <td>0.949296</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
  </tbody>
</table>
</div>




```python
asv_tables['16S V4-V5']['taxonomy'][0]
```




    'd__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Vibrionales; f__Vibrionaceae; g__Vibrio'




```python
taxa_ranks_16S = ['domain','phylum','class','order','family','genus','species']

```


```python
asv_tables['16S V4-V5'][['domain','phylum','class','order','family','genus','species']] = asv_tables['16S V4-V5']['taxonomy'].str.split("; ",expand=True)
asv_tables['16S V4-V5'].iloc[0:10,[0,1,2,3,4,5,6,-8,-7,-6,-5,-4,-3,-2,-1]]
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>featureid</th>
      <th>sequence</th>
      <th>taxonomy</th>
      <th>Confidence</th>
      <th>GOMECC4_27N_Sta1_DCM_A</th>
      <th>GOMECC4_27N_Sta1_DCM_B</th>
      <th>GOMECC4_27N_Sta1_DCM_C</th>
      <th>GOMECC4_YUCATAN_Sta102_Surface_C</th>
      <th>domain</th>
      <th>phylum</th>
      <th>class</th>
      <th>order</th>
      <th>family</th>
      <th>genus</th>
      <th>species</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>00006f0784f7dbb2f162408abb6da629</td>
      <td>TACGGAGGGTGCGAGCGTTAATCGGAATTACTGGGCGTAAAGCGCATGCAGGTGGTTTGTTAAGTCAGATGTGAAAGCCCGGGGCTCAACCTCGGAATTGCATTTGAAACTGGCAGACTAGAGTACTGTAGAGGGGGGTAGAATTT...</td>
      <td>d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Vibrionales; f__Vibrionaceae; g__Vibrio</td>
      <td>0.978926</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>7</td>
      <td>d__Bacteria</td>
      <td>p__Proteobacteria</td>
      <td>c__Gammaproteobacteria</td>
      <td>o__Vibrionales</td>
      <td>f__Vibrionaceae</td>
      <td>g__Vibrio</td>
      <td>None</td>
    </tr>
    <tr>
      <th>1</th>
      <td>000094731d4984ed41435a1bf65b7ef2</td>
      <td>TACAGAGAGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCGCGTAGGTGGGTATTTAAGTCGGATGTGAAATCCCCGGGCTTAACCTGGGAACTGCATCCGAAACTATTTAACTAGAGTATGGGAGAGGTAAGTAGAATTT...</td>
      <td>d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__HOC36; f__HOC36; g__HOC36; s__Candidatus_Thioglobus</td>
      <td>0.881698</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>d__Bacteria</td>
      <td>p__Proteobacteria</td>
      <td>c__Gammaproteobacteria</td>
      <td>o__HOC36</td>
      <td>f__HOC36</td>
      <td>g__HOC36</td>
      <td>s__Candidatus_Thioglobus</td>
    </tr>
    <tr>
      <th>2</th>
      <td>0001a3c11fcef1b1b8f4c72942efbbac</td>
      <td>TACGAAGGGGGCGAGCGTTGTTCGGAATTACTGGGCGTAAAGGGCGCGTAGGCGGTCTTCTAAGTTAGGCGTGAAAGCCCCGGGCTCAACCTGGGAACTGCGCTTAATACTGGAAGACTAGAAAACGGAAGAGGGTAGTGGAATTC...</td>
      <td>d__Bacteria; p__Cyanobacteria; c__Cyanobacteriia; o__Synechococcales; f__Cyanobiaceae; g__Cyanobium_PCC-6307</td>
      <td>0.762793</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>d__Bacteria</td>
      <td>p__Cyanobacteria</td>
      <td>c__Cyanobacteriia</td>
      <td>o__Synechococcales</td>
      <td>f__Cyanobiaceae</td>
      <td>g__Cyanobium_PCC-6307</td>
      <td>None</td>
    </tr>
    <tr>
      <th>3</th>
      <td>0001ceef5162e6d689ef30418cfcc164</td>
      <td>TACAGAGGGTGCAAGCGTTGTTCGGAATCATTGGGCGTAAAGCGCGCGTAGGCGGCCAAATAAGTCTGATGTGAAGGCCCAGGGCTCAACCCTGGAAGTGCATCGGAAACTGTTTGGCTCGAGTCCCGGAGGGGGTGGTGGAATTC...</td>
      <td>d__Bacteria; p__Myxococcota; c__Myxococcia; o__Myxococcales; f__Myxococcaceae; g__P3OB-42; s__uncultured_bacterium</td>
      <td>0.997619</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>d__Bacteria</td>
      <td>p__Myxococcota</td>
      <td>c__Myxococcia</td>
      <td>o__Myxococcales</td>
      <td>f__Myxococcaceae</td>
      <td>g__P3OB-42</td>
      <td>s__uncultured_bacterium</td>
    </tr>
    <tr>
      <th>4</th>
      <td>000235534662df05bb30219a4b978dac</td>
      <td>TACGGAAGGTCCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCGCGTAGGTGGTTTTTTAAGTTGGATGTGAAAGCCCTGGGCTCAACCTAGGAACTGCATCCAAAACTAGATGACTAGAGTACGAAAGAGGGAAGTAGAATTC...</td>
      <td>d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__SAR86_clade; f__SAR86_clade; g__SAR86_clade</td>
      <td>0.999961</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>d__Bacteria</td>
      <td>p__Proteobacteria</td>
      <td>c__Gammaproteobacteria</td>
      <td>o__SAR86_clade</td>
      <td>f__SAR86_clade</td>
      <td>g__SAR86_clade</td>
      <td>None</td>
    </tr>
    <tr>
      <th>5</th>
      <td>0003aeafc4bc0522877d4804829e65b7</td>
      <td>CACCGGCATCTCGAGTGGTATCCACTTTTATTGGGCCTAAAGCATCCGTAGCCTGTTCTGTAAGTTTTCGGTTAAATCCATAAGCTCAACTTATGGGCTGCCGAAAATACTGCAGAACTAGGGAGTGGGAGAGGTAGACGGTACTC...</td>
      <td>d__Archaea; p__Crenarchaeota; c__Nitrososphaeria; o__Nitrosopumilales; f__Nitrosopumilaceae; g__Candidatus_Nitrosopelagicus; s__marine_metagenome</td>
      <td>0.779476</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>d__Archaea</td>
      <td>p__Crenarchaeota</td>
      <td>c__Nitrososphaeria</td>
      <td>o__Nitrosopumilales</td>
      <td>f__Nitrosopumilaceae</td>
      <td>g__Candidatus_Nitrosopelagicus</td>
      <td>s__marine_metagenome</td>
    </tr>
    <tr>
      <th>6</th>
      <td>0003b46ae196127658c07aeb11b36b1a</td>
      <td>TACGGAGGGTCCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCGCGTAGGCGGTTTAGCAAGTTGAATGTGAAAGCCCTGGGCTCAACCTAGGAACTGCATTCAAAACTACTAAGCTAGAGTACGAGAGAGGAGAGTAGAATTT...</td>
      <td>d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Thiomicrospirales; f__Thioglobaceae; g__SUP05_cluster</td>
      <td>0.810269</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>d__Bacteria</td>
      <td>p__Proteobacteria</td>
      <td>c__Gammaproteobacteria</td>
      <td>o__Thiomicrospirales</td>
      <td>f__Thioglobaceae</td>
      <td>g__SUP05_cluster</td>
      <td>None</td>
    </tr>
    <tr>
      <th>7</th>
      <td>00065c7f5701f8db77fc2c50a1204c71</td>
      <td>TACGGGAGTGGCAAGCGTTATCCGGAATTATTGGGCGTAAAGCGTTTGTAGGTGGAAAAATAAGTCTATTGTTAAATCCAGAAGCTTAACTTCTGTCAAGCGATATGAAACTATTCTTCTTGAGAATGGTAGGGGTAGAAGGAATT...</td>
      <td>d__Bacteria; p__Cyanobacteria; c__Cyanobacteriia; o__Chloroplast; f__Chloroplast; g__Chloroplast; s__Prasinoderma_coloniale</td>
      <td>1.000000</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>d__Bacteria</td>
      <td>p__Cyanobacteria</td>
      <td>c__Cyanobacteriia</td>
      <td>o__Chloroplast</td>
      <td>f__Chloroplast</td>
      <td>g__Chloroplast</td>
      <td>s__Prasinoderma_coloniale</td>
    </tr>
    <tr>
      <th>8</th>
      <td>0006da4e1ff162826badd8bdcfaf9dfe</td>
      <td>GACGGAGGATGCAAGTGTTATCCGGAATTATTGGGCGTAAAGCGTTTGTAGGTGGAGAAATAAGCCTATTGTTAAATCCAGGAGCTTAACTTCTGTCCAGCGATATGAAACTATTTTTCTTGAGGGTGGTAGGGGTAGAAGGAATT...</td>
      <td>d__Bacteria; p__Cyanobacteria; c__Cyanobacteriia; o__Chloroplast; f__Chloroplast; g__Chloroplast; s__Prasinoderma_coloniale</td>
      <td>1.000000</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>d__Bacteria</td>
      <td>p__Cyanobacteria</td>
      <td>c__Cyanobacteriia</td>
      <td>o__Chloroplast</td>
      <td>f__Chloroplast</td>
      <td>g__Chloroplast</td>
      <td>s__Prasinoderma_coloniale</td>
    </tr>
    <tr>
      <th>9</th>
      <td>0006fe6033cca4da30ea5ce9cba446f0</td>
      <td>TACGAATGCTGCAAGCGTAGTTCGGAATCACTGGGCATAAAGAGCACGTAGGCGGCCTATTAAGTCAGCTGTGAAATCCCTCGGCTTAACCGAGGAACTGCAGCTGATACTGATAGGCTTGAGTACGGGAGGGGAGAGCGGAATTC...</td>
      <td>d__Bacteria; p__Planctomycetota; c__Pla3_lineage; o__Pla3_lineage; f__Pla3_lineage; g__Pla3_lineage; s__uncultured_Planctomycetaceae</td>
      <td>0.949296</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>d__Bacteria</td>
      <td>p__Planctomycetota</td>
      <td>c__Pla3_lineage</td>
      <td>o__Pla3_lineage</td>
      <td>f__Pla3_lineage</td>
      <td>g__Pla3_lineage</td>
      <td>s__uncultured_Planctomycetaceae</td>
    </tr>
  </tbody>
</table>
</div>




```python
asv_tables['16S V4-V5']['domain'] = asv_tables['16S V4-V5']['domain'].str.replace("d__", "")
asv_tables['16S V4-V5']['phylum'] = asv_tables['16S V4-V5']['phylum'].str.replace("p__", "")
asv_tables['16S V4-V5']['class'] = asv_tables['16S V4-V5']['class'].str.replace("c__", "")
asv_tables['16S V4-V5']['order'] = asv_tables['16S V4-V5']['order'].str.replace("o__", "")
asv_tables['16S V4-V5']['family'] = asv_tables['16S V4-V5']['family'].str.replace("f__", "")
asv_tables['16S V4-V5']['genus'] = asv_tables['16S V4-V5']['genus'].str.replace("g__", "")
asv_tables['16S V4-V5']['species'] = asv_tables['16S V4-V5']['species'].str.replace("s__", "")
```


```python
# replace None with NA
asv_tables['16S V4-V5'] = asv_tables['16S V4-V5'].fillna(value=np.nan)
## Replace 'unknown', 'unassigned', etc. in species and taxonomy columns with NaN

asv_tables['16S V4-V5'][taxa_ranks_16S] = asv_tables['16S V4-V5'][taxa_ranks_16S].replace({'unassigned':np.nan,'Unassigned':np.nan,
                              's_':np.nan,
                              'g_':np.nan,
                              'unknown':np.nan,
                              'no_hit':np.nan,
                               '':np.nan})
asv_tables['16S V4-V5'].iloc[0:10,[0,1,2,3,4,5,6,-8,-7,-6,-5,-4,-3,-2,-1]]
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>featureid</th>
      <th>sequence</th>
      <th>taxonomy</th>
      <th>Confidence</th>
      <th>GOMECC4_27N_Sta1_DCM_A</th>
      <th>GOMECC4_27N_Sta1_DCM_B</th>
      <th>GOMECC4_27N_Sta1_DCM_C</th>
      <th>GOMECC4_YUCATAN_Sta102_Surface_C</th>
      <th>domain</th>
      <th>phylum</th>
      <th>class</th>
      <th>order</th>
      <th>family</th>
      <th>genus</th>
      <th>species</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>00006f0784f7dbb2f162408abb6da629</td>
      <td>TACGGAGGGTGCGAGCGTTAATCGGAATTACTGGGCGTAAAGCGCATGCAGGTGGTTTGTTAAGTCAGATGTGAAAGCCCGGGGCTCAACCTCGGAATTGCATTTGAAACTGGCAGACTAGAGTACTGTAGAGGGGGGTAGAATTT...</td>
      <td>d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Vibrionales; f__Vibrionaceae; g__Vibrio</td>
      <td>0.978926</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>7</td>
      <td>Bacteria</td>
      <td>Proteobacteria</td>
      <td>Gammaproteobacteria</td>
      <td>Vibrionales</td>
      <td>Vibrionaceae</td>
      <td>Vibrio</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>1</th>
      <td>000094731d4984ed41435a1bf65b7ef2</td>
      <td>TACAGAGAGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCGCGTAGGTGGGTATTTAAGTCGGATGTGAAATCCCCGGGCTTAACCTGGGAACTGCATCCGAAACTATTTAACTAGAGTATGGGAGAGGTAAGTAGAATTT...</td>
      <td>d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__HOC36; f__HOC36; g__HOC36; s__Candidatus_Thioglobus</td>
      <td>0.881698</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>Bacteria</td>
      <td>Proteobacteria</td>
      <td>Gammaproteobacteria</td>
      <td>HOC36</td>
      <td>HOC36</td>
      <td>HOC36</td>
      <td>Candidatus_Thioglobus</td>
    </tr>
    <tr>
      <th>2</th>
      <td>0001a3c11fcef1b1b8f4c72942efbbac</td>
      <td>TACGAAGGGGGCGAGCGTTGTTCGGAATTACTGGGCGTAAAGGGCGCGTAGGCGGTCTTCTAAGTTAGGCGTGAAAGCCCCGGGCTCAACCTGGGAACTGCGCTTAATACTGGAAGACTAGAAAACGGAAGAGGGTAGTGGAATTC...</td>
      <td>d__Bacteria; p__Cyanobacteria; c__Cyanobacteriia; o__Synechococcales; f__Cyanobiaceae; g__Cyanobium_PCC-6307</td>
      <td>0.762793</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>Bacteria</td>
      <td>Cyanobacteria</td>
      <td>Cyanobacteriia</td>
      <td>Synechococcales</td>
      <td>Cyanobiaceae</td>
      <td>Cyanobium_PCC-6307</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>3</th>
      <td>0001ceef5162e6d689ef30418cfcc164</td>
      <td>TACAGAGGGTGCAAGCGTTGTTCGGAATCATTGGGCGTAAAGCGCGCGTAGGCGGCCAAATAAGTCTGATGTGAAGGCCCAGGGCTCAACCCTGGAAGTGCATCGGAAACTGTTTGGCTCGAGTCCCGGAGGGGGTGGTGGAATTC...</td>
      <td>d__Bacteria; p__Myxococcota; c__Myxococcia; o__Myxococcales; f__Myxococcaceae; g__P3OB-42; s__uncultured_bacterium</td>
      <td>0.997619</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>Bacteria</td>
      <td>Myxococcota</td>
      <td>Myxococcia</td>
      <td>Myxococcales</td>
      <td>Myxococcaceae</td>
      <td>P3OB-42</td>
      <td>uncultured_bacterium</td>
    </tr>
    <tr>
      <th>4</th>
      <td>000235534662df05bb30219a4b978dac</td>
      <td>TACGGAAGGTCCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCGCGTAGGTGGTTTTTTAAGTTGGATGTGAAAGCCCTGGGCTCAACCTAGGAACTGCATCCAAAACTAGATGACTAGAGTACGAAAGAGGGAAGTAGAATTC...</td>
      <td>d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__SAR86_clade; f__SAR86_clade; g__SAR86_clade</td>
      <td>0.999961</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>Bacteria</td>
      <td>Proteobacteria</td>
      <td>Gammaproteobacteria</td>
      <td>SAR86_clade</td>
      <td>SAR86_clade</td>
      <td>SAR86_clade</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>5</th>
      <td>0003aeafc4bc0522877d4804829e65b7</td>
      <td>CACCGGCATCTCGAGTGGTATCCACTTTTATTGGGCCTAAAGCATCCGTAGCCTGTTCTGTAAGTTTTCGGTTAAATCCATAAGCTCAACTTATGGGCTGCCGAAAATACTGCAGAACTAGGGAGTGGGAGAGGTAGACGGTACTC...</td>
      <td>d__Archaea; p__Crenarchaeota; c__Nitrososphaeria; o__Nitrosopumilales; f__Nitrosopumilaceae; g__Candidatus_Nitrosopelagicus; s__marine_metagenome</td>
      <td>0.779476</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>Archaea</td>
      <td>Crenarchaeota</td>
      <td>Nitrososphaeria</td>
      <td>Nitrosopumilales</td>
      <td>Nitrosopumilaceae</td>
      <td>Candidatus_Nitrosopelagicus</td>
      <td>marine_metagenome</td>
    </tr>
    <tr>
      <th>6</th>
      <td>0003b46ae196127658c07aeb11b36b1a</td>
      <td>TACGGAGGGTCCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCGCGTAGGCGGTTTAGCAAGTTGAATGTGAAAGCCCTGGGCTCAACCTAGGAACTGCATTCAAAACTACTAAGCTAGAGTACGAGAGAGGAGAGTAGAATTT...</td>
      <td>d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Thiomicrospirales; f__Thioglobaceae; g__SUP05_cluster</td>
      <td>0.810269</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>Bacteria</td>
      <td>Proteobacteria</td>
      <td>Gammaproteobacteria</td>
      <td>Thiomicrospirales</td>
      <td>Thioglobaceae</td>
      <td>SUP05_cluster</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>7</th>
      <td>00065c7f5701f8db77fc2c50a1204c71</td>
      <td>TACGGGAGTGGCAAGCGTTATCCGGAATTATTGGGCGTAAAGCGTTTGTAGGTGGAAAAATAAGTCTATTGTTAAATCCAGAAGCTTAACTTCTGTCAAGCGATATGAAACTATTCTTCTTGAGAATGGTAGGGGTAGAAGGAATT...</td>
      <td>d__Bacteria; p__Cyanobacteria; c__Cyanobacteriia; o__Chloroplast; f__Chloroplast; g__Chloroplast; s__Prasinoderma_coloniale</td>
      <td>1.000000</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>Bacteria</td>
      <td>Cyanobacteria</td>
      <td>Cyanobacteriia</td>
      <td>Chloroplast</td>
      <td>Chloroplast</td>
      <td>Chloroplast</td>
      <td>Prasinoderma_coloniale</td>
    </tr>
    <tr>
      <th>8</th>
      <td>0006da4e1ff162826badd8bdcfaf9dfe</td>
      <td>GACGGAGGATGCAAGTGTTATCCGGAATTATTGGGCGTAAAGCGTTTGTAGGTGGAGAAATAAGCCTATTGTTAAATCCAGGAGCTTAACTTCTGTCCAGCGATATGAAACTATTTTTCTTGAGGGTGGTAGGGGTAGAAGGAATT...</td>
      <td>d__Bacteria; p__Cyanobacteria; c__Cyanobacteriia; o__Chloroplast; f__Chloroplast; g__Chloroplast; s__Prasinoderma_coloniale</td>
      <td>1.000000</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>Bacteria</td>
      <td>Cyanobacteria</td>
      <td>Cyanobacteriia</td>
      <td>Chloroplast</td>
      <td>Chloroplast</td>
      <td>Chloroplast</td>
      <td>Prasinoderma_coloniale</td>
    </tr>
    <tr>
      <th>9</th>
      <td>0006fe6033cca4da30ea5ce9cba446f0</td>
      <td>TACGAATGCTGCAAGCGTAGTTCGGAATCACTGGGCATAAAGAGCACGTAGGCGGCCTATTAAGTCAGCTGTGAAATCCCTCGGCTTAACCGAGGAACTGCAGCTGATACTGATAGGCTTGAGTACGGGAGGGGAGAGCGGAATTC...</td>
      <td>d__Bacteria; p__Planctomycetota; c__Pla3_lineage; o__Pla3_lineage; f__Pla3_lineage; g__Pla3_lineage; s__uncultured_Planctomycetaceae</td>
      <td>0.949296</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>Bacteria</td>
      <td>Planctomycetota</td>
      <td>Pla3_lineage</td>
      <td>Pla3_lineage</td>
      <td>Pla3_lineage</td>
      <td>Pla3_lineage</td>
      <td>uncultured_Planctomycetaceae</td>
    </tr>
  </tbody>
</table>
</div>




```python
# replace _,- with space, remove sp. 

asv_tables['16S V4-V5'][taxa_ranks_16S] = asv_tables['16S V4-V5'][taxa_ranks_16S].replace('_',' ',regex=True)
asv_tables['16S V4-V5'][taxa_ranks_16S] = asv_tables['16S V4-V5'][taxa_ranks_16S].replace(' sp\.','',regex=True)
asv_tables['16S V4-V5'][taxa_ranks_16S] = asv_tables['16S V4-V5'][taxa_ranks_16S].replace('-',' ',regex=True)
asv_tables['16S V4-V5'][taxa_ranks_16S] = asv_tables['16S V4-V5'][taxa_ranks_16S].replace(' spp\.','',regex=True)
asv_tables['16S V4-V5'][taxa_ranks_16S] = asv_tables['16S V4-V5'][taxa_ranks_16S].replace('\/',' ',regex=True)
```

##### Melt asv_tables to long format



```python
asv_tables['16S V4-V5'].shape

```




    (65048, 483)




```python
occ['16S V4-V5'] = pd.melt(asv_tables['16S V4-V5'],id_vars=['featureid','sequence','taxonomy','Confidence','domain','phylum','class','order','family','genus','species'],
               var_name='eventID',value_name='organismQuantity')
```


```python
occ['16S V4-V5'].shape
```




    (30702656, 13)




```python
## Drop records where organismQuantity = 0 (absences are not meaningful for this data set)

occ['16S V4-V5'] = occ['16S V4-V5'][occ['16S V4-V5']['organismQuantity'] > 0]
print(occ['16S V4-V5'].shape)
```

    (165158, 13)



```python
## Create an occurrenceID that will uniquely identify each ASV observed within a water sample

occ['16S V4-V5']['occurrenceID'] = occ['16S V4-V5']['featureid']
occ['16S V4-V5']['occurrenceID'] = occ['16S V4-V5']['eventID'] + '_16S_occ' + occ['16S V4-V5']['occurrenceID'].astype(str)
```


```python
occ['16S V4-V5'].head()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>featureid</th>
      <th>sequence</th>
      <th>taxonomy</th>
      <th>Confidence</th>
      <th>domain</th>
      <th>phylum</th>
      <th>class</th>
      <th>order</th>
      <th>family</th>
      <th>genus</th>
      <th>species</th>
      <th>eventID</th>
      <th>organismQuantity</th>
      <th>occurrenceID</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>182</th>
      <td>00c4c1c65d8669ed9f07abe149f9a01d</td>
      <td>TACGGAGGGGGCTAACGTTGTTCGGAATTACTGGGCGTAAAGCGCGCGTAGGCGGATTAGACAGTTGAGGGTGAAATCCCGGAGCTTAACTTCGGAACTGCCCCCAATACTACTAATCTAGAGTTCGGAAGAGGTGAGTGGAATTC...</td>
      <td>d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Parvibaculales; f__OCS116_clade; g__OCS116_clade; s__uncultured_marine</td>
      <td>0.832190</td>
      <td>Bacteria</td>
      <td>Proteobacteria</td>
      <td>Alphaproteobacteria</td>
      <td>Parvibaculales</td>
      <td>OCS116 clade</td>
      <td>OCS116 clade</td>
      <td>uncultured marine</td>
      <td>GOMECC4_27N_Sta1_DCM_A</td>
      <td>18</td>
      <td>GOMECC4_27N_Sta1_DCM_A_16S_occ00c4c1c65d8669ed9f07abe149f9a01d</td>
    </tr>
    <tr>
      <th>225</th>
      <td>00e6c13fe86364a5084987093afa1916</td>
      <td>TACGAAGGGGGCGAGCGTTGTTCGGAATTACTGGGCGTAAAGGGCGCGTAGGCGGCTCTTTAAGTTAGGCGTGAAAGCCCCGGGCTCAACCTGGGAACTGCGCTTAAGACTGGAGAGCTAGAAAACGGAAGAGGGTAGTGGAATTC...</td>
      <td>d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Puniceispirillales; f__SAR116_clade; g__SAR116_clade</td>
      <td>0.867040</td>
      <td>Bacteria</td>
      <td>Proteobacteria</td>
      <td>Alphaproteobacteria</td>
      <td>Puniceispirillales</td>
      <td>SAR116 clade</td>
      <td>SAR116 clade</td>
      <td>NaN</td>
      <td>GOMECC4_27N_Sta1_DCM_A</td>
      <td>36</td>
      <td>GOMECC4_27N_Sta1_DCM_A_16S_occ00e6c13fe86364a5084987093afa1916</td>
    </tr>
    <tr>
      <th>347</th>
      <td>015dad1fafca90944d905beb2a980bc3</td>
      <td>TACCGGCGCCTCAAGTGGTAGTCGCTTTTATTGGGCCTAAAACGTCCGTAGCCGGTCTGGTACATTCGTGGGTAAATCAACTCGCTTAACGAGTTGAATTCTGCGAGGACGGCCAGACTTGGGACCGGGAGAGGTGTGGGGTACTC...</td>
      <td>d__Archaea; p__Thermoplasmatota; c__Thermoplasmata; o__Marine_Group_II; f__Marine_Group_II; g__Marine_Group_II</td>
      <td>1.000000</td>
      <td>Archaea</td>
      <td>Thermoplasmatota</td>
      <td>Thermoplasmata</td>
      <td>Marine Group II</td>
      <td>Marine Group II</td>
      <td>Marine Group II</td>
      <td>NaN</td>
      <td>GOMECC4_27N_Sta1_DCM_A</td>
      <td>49</td>
      <td>GOMECC4_27N_Sta1_DCM_A_16S_occ015dad1fafca90944d905beb2a980bc3</td>
    </tr>
    <tr>
      <th>412</th>
      <td>019c88c6ade406f731954f38e3461564</td>
      <td>TACAGGAGGGACGAGTGTTACTCGGAATGATTAGGCGTAAAGGGTCATTTAAGCGGTCCGATAAGTTAAAAGCCAACAGTTAGAGCCTAACTCTTTCAAGCTTTTAATACTGTCAGACTAGAGTATATCAGAGAATAGTAGAATTC...</td>
      <td>d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rickettsiales; f__Mitochondria; g__Mitochondria; s__uncultured_bacterium</td>
      <td>0.952911</td>
      <td>Bacteria</td>
      <td>Proteobacteria</td>
      <td>Alphaproteobacteria</td>
      <td>Rickettsiales</td>
      <td>Mitochondria</td>
      <td>Mitochondria</td>
      <td>uncultured bacterium</td>
      <td>GOMECC4_27N_Sta1_DCM_A</td>
      <td>2</td>
      <td>GOMECC4_27N_Sta1_DCM_A_16S_occ019c88c6ade406f731954f38e3461564</td>
    </tr>
    <tr>
      <th>719</th>
      <td>02dfb0869af4bf549d290d48e66e2196</td>
      <td>TACGAGGGGTGCTAGCGTTGTCCGGAATAACTGGGCGTAAAGGGTCCGTAGGCGTTTTGCTAAGTTGATCGTTAAATCCATCGGCTTAACCGATGACATGCGATCAAAACTGGCAGAATAGAATATGTGAGGGGAATGTAGAATTC...</td>
      <td>d__Bacteria; p__Marinimicrobia_(SAR406_clade); c__Marinimicrobia_(SAR406_clade); o__Marinimicrobia_(SAR406_clade); f__Marinimicrobia_(SAR406_clade...</td>
      <td>0.818195</td>
      <td>Bacteria</td>
      <td>Marinimicrobia (SAR406 clade)</td>
      <td>Marinimicrobia (SAR406 clade)</td>
      <td>Marinimicrobia (SAR406 clade)</td>
      <td>Marinimicrobia (SAR406 clade)</td>
      <td>Marinimicrobia (SAR406 clade)</td>
      <td>uncultured bacterium</td>
      <td>GOMECC4_27N_Sta1_DCM_A</td>
      <td>3</td>
      <td>GOMECC4_27N_Sta1_DCM_A_16S_occ02dfb0869af4bf549d290d48e66e2196</td>
    </tr>
  </tbody>
</table>
</div>



##### WORMS conversion. 
Note, can't use `multiprocessing` library in a Jupyter notebook, need `multiprocess`. See [here](https://stackoverflow.com/questions/41385708/multiprocessing-example-giving-attributeerror)

OBIS currently requires taxonomy assignments that match WoRMS, however none of the commonly used metabarcoding reference databases use WoRMS as the basis of their taxonomy. This means the taxonomic ranks for any given scientific name on WoRMS may not directly compare to what is assigned. There are ongoing discussions about this problem (see [this](https://github.com/iobis/Project-team-Genetic-Data/issues/5) GitHub issue).     

Many of them, especially for microbes, include taxa that aren't on WoRMS at all. This is because the name may not have been fully and officially adopted by the scientific community (or at least not adopted by WoRMS yet). We therefore need a system for searching through the higher taxonomic ranks given, finding the lowest one that will match on WoRMS, and putting that name in the `scientificName` column. The assigned taxonomy is then recorded elsewhere.

Had some [issues with the parallelization](https://stackoverflow.com/questions/50168647/multiprocessing-causes-python-to-crash-and-gives-an-error-may-have-been-in-progr) on Mac M1. Adding 'OBJC_DISABLE_INITIALIZE_FORK_SAFETY = YES' to .bash_profile and then [This](https://github.com/python/cpython/issues/74570) fixed it.   
Try to run without the bash_profile fix LATER.


```python
os.environ["no_proxy"]="*"
```

### 16S worms

Species level IDs might be trash, [see here](https://forum.qiime2.org/t/processing-filtering-and-evaluating-the-silva-database-and-other-reference-sequence-data-with-rescript/15494), so look at genus and up.


```python
import WoRMS_matching
```


```python
import importlib
importlib.reload(WoRMS_matching)
```




    <module 'WoRMS_matching' from '/Users/katherine.silliman/Projects/NOAA/DwC/edna2obis/src/WoRMS_matching.py'>




```python
tax_16S = asv_tables['16S V4-V5'][['taxonomy','domain','phylum','class','order','family','genus','species']]
```


```python
#ignore_index is important!
tax_16S = tax_16S.drop_duplicates(ignore_index=True)
```


```python
tax_16S.shape
```




    (2729, 8)




```python
if __name__ == '__main__':
    worms_16s = WoRMS_matching.get_worms_from_scientific_name_parallel(
    tax_df = tax_16S,ordered_rank_columns=['genus','family','order','class','phylum','domain'],
    full_tax_column="taxonomy",full_tax_vI=True,n_proc=7)
```


```python
worms_16s.head()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>full_tax</th>
      <th>verbatimIdentification</th>
      <th>old_taxonRank</th>
      <th>old name</th>
      <th>scientificName</th>
      <th>scientificNameID</th>
      <th>kingdom</th>
      <th>phylum</th>
      <th>class</th>
      <th>order</th>
      <th>family</th>
      <th>genus</th>
      <th>taxonRank</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Vibrionales; f__Vibrionaceae; g__Vibrio</td>
      <td>d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Vibrionales; f__Vibrionaceae; g__Vibrio</td>
      <td>genus</td>
      <td>Vibrio</td>
      <td>Vibrio</td>
      <td>urn:lsid:marinespecies.org:taxname:480248</td>
      <td>Bacteria</td>
      <td>Proteobacteria</td>
      <td>Gammaproteobacteria</td>
      <td>Vibrionales</td>
      <td>Vibrionaceae</td>
      <td>Vibrio</td>
      <td>Genus</td>
    </tr>
    <tr>
      <th>1</th>
      <td>d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__HOC36; f__HOC36; g__HOC36; s__Candidatus_Thioglobus</td>
      <td>d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__HOC36; f__HOC36; g__HOC36; s__Candidatus_Thioglobus</td>
      <td>class</td>
      <td>Gammaproteobacteria</td>
      <td>Gammaproteobacteria</td>
      <td>urn:lsid:marinespecies.org:taxname:393018</td>
      <td>Bacteria</td>
      <td>Proteobacteria</td>
      <td>Gammaproteobacteria</td>
      <td>None</td>
      <td>None</td>
      <td>None</td>
      <td>Class</td>
    </tr>
    <tr>
      <th>2</th>
      <td>d__Bacteria; p__Cyanobacteria; c__Cyanobacteriia; o__Synechococcales; f__Cyanobiaceae; g__Cyanobium_PCC-6307</td>
      <td>d__Bacteria; p__Cyanobacteria; c__Cyanobacteriia; o__Synechococcales; f__Cyanobiaceae; g__Cyanobium_PCC-6307</td>
      <td>order</td>
      <td>Synechococcales</td>
      <td>Synechococcales</td>
      <td>urn:lsid:marinespecies.org:taxname:345514</td>
      <td>Bacteria</td>
      <td>Cyanobacteria</td>
      <td>Cyanophyceae</td>
      <td>Synechococcales</td>
      <td>None</td>
      <td>None</td>
      <td>Order</td>
    </tr>
    <tr>
      <th>3</th>
      <td>d__Bacteria; p__Myxococcota; c__Myxococcia; o__Myxococcales; f__Myxococcaceae; g__P3OB-42; s__uncultured_bacterium</td>
      <td>d__Bacteria; p__Myxococcota; c__Myxococcia; o__Myxococcales; f__Myxococcaceae; g__P3OB-42; s__uncultured_bacterium</td>
      <td>family</td>
      <td>Myxococcaceae</td>
      <td>Myxococcaceae</td>
      <td>urn:lsid:marinespecies.org:taxname:570956</td>
      <td>Bacteria</td>
      <td>Proteobacteria</td>
      <td>Deltaproteobacteria</td>
      <td>Myxococcales</td>
      <td>Myxococcaceae</td>
      <td>None</td>
      <td>Family</td>
    </tr>
    <tr>
      <th>4</th>
      <td>d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__SAR86_clade; f__SAR86_clade; g__SAR86_clade</td>
      <td>d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__SAR86_clade; f__SAR86_clade; g__SAR86_clade</td>
      <td>class</td>
      <td>Gammaproteobacteria</td>
      <td>Gammaproteobacteria</td>
      <td>urn:lsid:marinespecies.org:taxname:393018</td>
      <td>Bacteria</td>
      <td>Proteobacteria</td>
      <td>Gammaproteobacteria</td>
      <td>None</td>
      <td>None</td>
      <td>None</td>
      <td>Class</td>
    </tr>
  </tbody>
</table>
</div>




```python
worms_16s[worms_16s["scientificName"]=="No match"]
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>full_tax</th>
      <th>verbatimIdentification</th>
      <th>old_taxonRank</th>
      <th>old name</th>
      <th>scientificName</th>
      <th>scientificNameID</th>
      <th>kingdom</th>
      <th>phylum</th>
      <th>class</th>
      <th>order</th>
      <th>family</th>
      <th>genus</th>
      <th>taxonRank</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>242</th>
      <td>d__Eukaryota</td>
      <td>d__Eukaryota</td>
      <td>domain</td>
      <td>Eukaryota</td>
      <td>No match</td>
      <td>None</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
  </tbody>
</table>
</div>




```python
worms_16s.loc[worms_16s["scientificName"]=="No match",'scientificName'] = "Biota"
worms_16s.loc[worms_16s["scientificName"]=="Biota",'scientificNameID'] = "urn:lsid:marinespecies.org:taxname:1"

```


```python
worms_16s[worms_16s['scientificName'].isna() == True]
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>full_tax</th>
      <th>verbatimIdentification</th>
      <th>old_taxonRank</th>
      <th>old name</th>
      <th>scientificName</th>
      <th>scientificNameID</th>
      <th>kingdom</th>
      <th>phylum</th>
      <th>class</th>
      <th>order</th>
      <th>family</th>
      <th>genus</th>
      <th>taxonRank</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>97</th>
      <td>Unassigned</td>
      <td>Unassigned</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
  </tbody>
</table>
</div>




```python

print(worms_16s[worms_16s['scientificName'].isna() == True].shape)
worms_16s.loc[worms_16s['scientificName'].isna() == True,'scientificName'] = 'incertae sedis'
worms_16s.loc[worms_16s['scientificName'] == 'incertae sedis','scientificNameID'] =  'urn:lsid:marinespecies.org:taxname:12'
print(worms_16s[worms_16s['scientificName'].isna() == True].shape)
```

    (1, 13)
    (0, 13)



```python
worms_16s.to_csv("../processed/worms_16S_matching.tsv",sep="\t",index=False)
```


```python
worms_16s.drop(columns=['old name','old_taxonRank'],inplace=True)
worms_16s.head()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>full_tax</th>
      <th>verbatimIdentification</th>
      <th>scientificName</th>
      <th>scientificNameID</th>
      <th>kingdom</th>
      <th>phylum</th>
      <th>class</th>
      <th>order</th>
      <th>family</th>
      <th>genus</th>
      <th>taxonRank</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Vibrionales; f__Vibrionaceae; g__Vibrio</td>
      <td>d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Vibrionales; f__Vibrionaceae; g__Vibrio</td>
      <td>Vibrio</td>
      <td>urn:lsid:marinespecies.org:taxname:480248</td>
      <td>Bacteria</td>
      <td>Proteobacteria</td>
      <td>Gammaproteobacteria</td>
      <td>Vibrionales</td>
      <td>Vibrionaceae</td>
      <td>Vibrio</td>
      <td>Genus</td>
    </tr>
    <tr>
      <th>1</th>
      <td>d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__HOC36; f__HOC36; g__HOC36; s__Candidatus_Thioglobus</td>
      <td>d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__HOC36; f__HOC36; g__HOC36; s__Candidatus_Thioglobus</td>
      <td>Gammaproteobacteria</td>
      <td>urn:lsid:marinespecies.org:taxname:393018</td>
      <td>Bacteria</td>
      <td>Proteobacteria</td>
      <td>Gammaproteobacteria</td>
      <td>None</td>
      <td>None</td>
      <td>None</td>
      <td>Class</td>
    </tr>
    <tr>
      <th>2</th>
      <td>d__Bacteria; p__Cyanobacteria; c__Cyanobacteriia; o__Synechococcales; f__Cyanobiaceae; g__Cyanobium_PCC-6307</td>
      <td>d__Bacteria; p__Cyanobacteria; c__Cyanobacteriia; o__Synechococcales; f__Cyanobiaceae; g__Cyanobium_PCC-6307</td>
      <td>Synechococcales</td>
      <td>urn:lsid:marinespecies.org:taxname:345514</td>
      <td>Bacteria</td>
      <td>Cyanobacteria</td>
      <td>Cyanophyceae</td>
      <td>Synechococcales</td>
      <td>None</td>
      <td>None</td>
      <td>Order</td>
    </tr>
    <tr>
      <th>3</th>
      <td>d__Bacteria; p__Myxococcota; c__Myxococcia; o__Myxococcales; f__Myxococcaceae; g__P3OB-42; s__uncultured_bacterium</td>
      <td>d__Bacteria; p__Myxococcota; c__Myxococcia; o__Myxococcales; f__Myxococcaceae; g__P3OB-42; s__uncultured_bacterium</td>
      <td>Myxococcaceae</td>
      <td>urn:lsid:marinespecies.org:taxname:570956</td>
      <td>Bacteria</td>
      <td>Proteobacteria</td>
      <td>Deltaproteobacteria</td>
      <td>Myxococcales</td>
      <td>Myxococcaceae</td>
      <td>None</td>
      <td>Family</td>
    </tr>
    <tr>
      <th>4</th>
      <td>d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__SAR86_clade; f__SAR86_clade; g__SAR86_clade</td>
      <td>d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__SAR86_clade; f__SAR86_clade; g__SAR86_clade</td>
      <td>Gammaproteobacteria</td>
      <td>urn:lsid:marinespecies.org:taxname:393018</td>
      <td>Bacteria</td>
      <td>Proteobacteria</td>
      <td>Gammaproteobacteria</td>
      <td>None</td>
      <td>None</td>
      <td>None</td>
      <td>Class</td>
    </tr>
  </tbody>
</table>
</div>




```python
occ['16S V4-V5'].head()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>featureid</th>
      <th>sequence</th>
      <th>taxonomy</th>
      <th>Confidence</th>
      <th>domain</th>
      <th>phylum</th>
      <th>class</th>
      <th>order</th>
      <th>family</th>
      <th>genus</th>
      <th>species</th>
      <th>eventID</th>
      <th>organismQuantity</th>
      <th>occurrenceID</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>182</th>
      <td>00c4c1c65d8669ed9f07abe149f9a01d</td>
      <td>TACGGAGGGGGCTAACGTTGTTCGGAATTACTGGGCGTAAAGCGCGCGTAGGCGGATTAGACAGTTGAGGGTGAAATCCCGGAGCTTAACTTCGGAACTGCCCCCAATACTACTAATCTAGAGTTCGGAAGAGGTGAGTGGAATTC...</td>
      <td>d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Parvibaculales; f__OCS116_clade; g__OCS116_clade; s__uncultured_marine</td>
      <td>0.832190</td>
      <td>Bacteria</td>
      <td>Proteobacteria</td>
      <td>Alphaproteobacteria</td>
      <td>Parvibaculales</td>
      <td>OCS116 clade</td>
      <td>OCS116 clade</td>
      <td>uncultured marine</td>
      <td>GOMECC4_27N_Sta1_DCM_A</td>
      <td>18</td>
      <td>GOMECC4_27N_Sta1_DCM_A_16S_occ00c4c1c65d8669ed9f07abe149f9a01d</td>
    </tr>
    <tr>
      <th>225</th>
      <td>00e6c13fe86364a5084987093afa1916</td>
      <td>TACGAAGGGGGCGAGCGTTGTTCGGAATTACTGGGCGTAAAGGGCGCGTAGGCGGCTCTTTAAGTTAGGCGTGAAAGCCCCGGGCTCAACCTGGGAACTGCGCTTAAGACTGGAGAGCTAGAAAACGGAAGAGGGTAGTGGAATTC...</td>
      <td>d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Puniceispirillales; f__SAR116_clade; g__SAR116_clade</td>
      <td>0.867040</td>
      <td>Bacteria</td>
      <td>Proteobacteria</td>
      <td>Alphaproteobacteria</td>
      <td>Puniceispirillales</td>
      <td>SAR116 clade</td>
      <td>SAR116 clade</td>
      <td>NaN</td>
      <td>GOMECC4_27N_Sta1_DCM_A</td>
      <td>36</td>
      <td>GOMECC4_27N_Sta1_DCM_A_16S_occ00e6c13fe86364a5084987093afa1916</td>
    </tr>
    <tr>
      <th>347</th>
      <td>015dad1fafca90944d905beb2a980bc3</td>
      <td>TACCGGCGCCTCAAGTGGTAGTCGCTTTTATTGGGCCTAAAACGTCCGTAGCCGGTCTGGTACATTCGTGGGTAAATCAACTCGCTTAACGAGTTGAATTCTGCGAGGACGGCCAGACTTGGGACCGGGAGAGGTGTGGGGTACTC...</td>
      <td>d__Archaea; p__Thermoplasmatota; c__Thermoplasmata; o__Marine_Group_II; f__Marine_Group_II; g__Marine_Group_II</td>
      <td>1.000000</td>
      <td>Archaea</td>
      <td>Thermoplasmatota</td>
      <td>Thermoplasmata</td>
      <td>Marine Group II</td>
      <td>Marine Group II</td>
      <td>Marine Group II</td>
      <td>NaN</td>
      <td>GOMECC4_27N_Sta1_DCM_A</td>
      <td>49</td>
      <td>GOMECC4_27N_Sta1_DCM_A_16S_occ015dad1fafca90944d905beb2a980bc3</td>
    </tr>
    <tr>
      <th>412</th>
      <td>019c88c6ade406f731954f38e3461564</td>
      <td>TACAGGAGGGACGAGTGTTACTCGGAATGATTAGGCGTAAAGGGTCATTTAAGCGGTCCGATAAGTTAAAAGCCAACAGTTAGAGCCTAACTCTTTCAAGCTTTTAATACTGTCAGACTAGAGTATATCAGAGAATAGTAGAATTC...</td>
      <td>d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rickettsiales; f__Mitochondria; g__Mitochondria; s__uncultured_bacterium</td>
      <td>0.952911</td>
      <td>Bacteria</td>
      <td>Proteobacteria</td>
      <td>Alphaproteobacteria</td>
      <td>Rickettsiales</td>
      <td>Mitochondria</td>
      <td>Mitochondria</td>
      <td>uncultured bacterium</td>
      <td>GOMECC4_27N_Sta1_DCM_A</td>
      <td>2</td>
      <td>GOMECC4_27N_Sta1_DCM_A_16S_occ019c88c6ade406f731954f38e3461564</td>
    </tr>
    <tr>
      <th>719</th>
      <td>02dfb0869af4bf549d290d48e66e2196</td>
      <td>TACGAGGGGTGCTAGCGTTGTCCGGAATAACTGGGCGTAAAGGGTCCGTAGGCGTTTTGCTAAGTTGATCGTTAAATCCATCGGCTTAACCGATGACATGCGATCAAAACTGGCAGAATAGAATATGTGAGGGGAATGTAGAATTC...</td>
      <td>d__Bacteria; p__Marinimicrobia_(SAR406_clade); c__Marinimicrobia_(SAR406_clade); o__Marinimicrobia_(SAR406_clade); f__Marinimicrobia_(SAR406_clade...</td>
      <td>0.818195</td>
      <td>Bacteria</td>
      <td>Marinimicrobia (SAR406 clade)</td>
      <td>Marinimicrobia (SAR406 clade)</td>
      <td>Marinimicrobia (SAR406 clade)</td>
      <td>Marinimicrobia (SAR406 clade)</td>
      <td>Marinimicrobia (SAR406 clade)</td>
      <td>uncultured bacterium</td>
      <td>GOMECC4_27N_Sta1_DCM_A</td>
      <td>3</td>
      <td>GOMECC4_27N_Sta1_DCM_A_16S_occ02dfb0869af4bf549d290d48e66e2196</td>
    </tr>
  </tbody>
</table>
</div>



#### Merge Occurrence and worms


```python
occ['16S V4-V5'].shape
```




    (165158, 14)




```python

occ16_test = occ['16S V4-V5'].copy()
occ16_test.drop(columns=['domain','phylum','class','order','family','genus','species'],inplace=True)

occ16_test = occ16_test.merge(worms_16s, how='left', left_on ='taxonomy', right_on='full_tax')
occ16_test.drop(columns='full_tax', inplace=True)
occ16_test.head()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>featureid</th>
      <th>sequence</th>
      <th>taxonomy</th>
      <th>Confidence</th>
      <th>eventID</th>
      <th>organismQuantity</th>
      <th>occurrenceID</th>
      <th>verbatimIdentification</th>
      <th>scientificName</th>
      <th>scientificNameID</th>
      <th>kingdom</th>
      <th>phylum</th>
      <th>class</th>
      <th>order</th>
      <th>family</th>
      <th>genus</th>
      <th>taxonRank</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>00c4c1c65d8669ed9f07abe149f9a01d</td>
      <td>TACGGAGGGGGCTAACGTTGTTCGGAATTACTGGGCGTAAAGCGCGCGTAGGCGGATTAGACAGTTGAGGGTGAAATCCCGGAGCTTAACTTCGGAACTGCCCCCAATACTACTAATCTAGAGTTCGGAAGAGGTGAGTGGAATTC...</td>
      <td>d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Parvibaculales; f__OCS116_clade; g__OCS116_clade; s__uncultured_marine</td>
      <td>0.832190</td>
      <td>GOMECC4_27N_Sta1_DCM_A</td>
      <td>18</td>
      <td>GOMECC4_27N_Sta1_DCM_A_16S_occ00c4c1c65d8669ed9f07abe149f9a01d</td>
      <td>d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Parvibaculales; f__OCS116_clade; g__OCS116_clade; s__uncultured_marine</td>
      <td>Alphaproteobacteria</td>
      <td>urn:lsid:marinespecies.org:taxname:392750</td>
      <td>Bacteria</td>
      <td>Proteobacteria</td>
      <td>Alphaproteobacteria</td>
      <td>None</td>
      <td>None</td>
      <td>None</td>
      <td>Class</td>
    </tr>
    <tr>
      <th>1</th>
      <td>00e6c13fe86364a5084987093afa1916</td>
      <td>TACGAAGGGGGCGAGCGTTGTTCGGAATTACTGGGCGTAAAGGGCGCGTAGGCGGCTCTTTAAGTTAGGCGTGAAAGCCCCGGGCTCAACCTGGGAACTGCGCTTAAGACTGGAGAGCTAGAAAACGGAAGAGGGTAGTGGAATTC...</td>
      <td>d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Puniceispirillales; f__SAR116_clade; g__SAR116_clade</td>
      <td>0.867040</td>
      <td>GOMECC4_27N_Sta1_DCM_A</td>
      <td>36</td>
      <td>GOMECC4_27N_Sta1_DCM_A_16S_occ00e6c13fe86364a5084987093afa1916</td>
      <td>d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Puniceispirillales; f__SAR116_clade; g__SAR116_clade</td>
      <td>Alphaproteobacteria</td>
      <td>urn:lsid:marinespecies.org:taxname:392750</td>
      <td>Bacteria</td>
      <td>Proteobacteria</td>
      <td>Alphaproteobacteria</td>
      <td>None</td>
      <td>None</td>
      <td>None</td>
      <td>Class</td>
    </tr>
    <tr>
      <th>2</th>
      <td>015dad1fafca90944d905beb2a980bc3</td>
      <td>TACCGGCGCCTCAAGTGGTAGTCGCTTTTATTGGGCCTAAAACGTCCGTAGCCGGTCTGGTACATTCGTGGGTAAATCAACTCGCTTAACGAGTTGAATTCTGCGAGGACGGCCAGACTTGGGACCGGGAGAGGTGTGGGGTACTC...</td>
      <td>d__Archaea; p__Thermoplasmatota; c__Thermoplasmata; o__Marine_Group_II; f__Marine_Group_II; g__Marine_Group_II</td>
      <td>1.000000</td>
      <td>GOMECC4_27N_Sta1_DCM_A</td>
      <td>49</td>
      <td>GOMECC4_27N_Sta1_DCM_A_16S_occ015dad1fafca90944d905beb2a980bc3</td>
      <td>d__Archaea; p__Thermoplasmatota; c__Thermoplasmata; o__Marine_Group_II; f__Marine_Group_II; g__Marine_Group_II</td>
      <td>Thermoplasmata</td>
      <td>urn:lsid:marinespecies.org:taxname:416268</td>
      <td>Archaea</td>
      <td>Euryarchaeota</td>
      <td>Thermoplasmata</td>
      <td>None</td>
      <td>None</td>
      <td>None</td>
      <td>Class</td>
    </tr>
    <tr>
      <th>3</th>
      <td>019c88c6ade406f731954f38e3461564</td>
      <td>TACAGGAGGGACGAGTGTTACTCGGAATGATTAGGCGTAAAGGGTCATTTAAGCGGTCCGATAAGTTAAAAGCCAACAGTTAGAGCCTAACTCTTTCAAGCTTTTAATACTGTCAGACTAGAGTATATCAGAGAATAGTAGAATTC...</td>
      <td>d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rickettsiales; f__Mitochondria; g__Mitochondria; s__uncultured_bacterium</td>
      <td>0.952911</td>
      <td>GOMECC4_27N_Sta1_DCM_A</td>
      <td>2</td>
      <td>GOMECC4_27N_Sta1_DCM_A_16S_occ019c88c6ade406f731954f38e3461564</td>
      <td>d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rickettsiales; f__Mitochondria; g__Mitochondria; s__uncultured_bacterium</td>
      <td>Rickettsiales</td>
      <td>urn:lsid:marinespecies.org:taxname:570969</td>
      <td>Bacteria</td>
      <td>Proteobacteria</td>
      <td>Alphaproteobacteria</td>
      <td>Rickettsiales</td>
      <td>None</td>
      <td>None</td>
      <td>Order</td>
    </tr>
    <tr>
      <th>4</th>
      <td>02dfb0869af4bf549d290d48e66e2196</td>
      <td>TACGAGGGGTGCTAGCGTTGTCCGGAATAACTGGGCGTAAAGGGTCCGTAGGCGTTTTGCTAAGTTGATCGTTAAATCCATCGGCTTAACCGATGACATGCGATCAAAACTGGCAGAATAGAATATGTGAGGGGAATGTAGAATTC...</td>
      <td>d__Bacteria; p__Marinimicrobia_(SAR406_clade); c__Marinimicrobia_(SAR406_clade); o__Marinimicrobia_(SAR406_clade); f__Marinimicrobia_(SAR406_clade...</td>
      <td>0.818195</td>
      <td>GOMECC4_27N_Sta1_DCM_A</td>
      <td>3</td>
      <td>GOMECC4_27N_Sta1_DCM_A_16S_occ02dfb0869af4bf549d290d48e66e2196</td>
      <td>d__Bacteria; p__Marinimicrobia_(SAR406_clade); c__Marinimicrobia_(SAR406_clade); o__Marinimicrobia_(SAR406_clade); f__Marinimicrobia_(SAR406_clade...</td>
      <td>Bacteria</td>
      <td>urn:lsid:marinespecies.org:taxname:6</td>
      <td>Bacteria</td>
      <td>None</td>
      <td>None</td>
      <td>None</td>
      <td>None</td>
      <td>None</td>
      <td>Kingdom</td>
    </tr>
  </tbody>
</table>
</div>



#### identificationRemarks


```python
data['analysis_data'].head()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>amplicon_sequenced</th>
      <th>ampliconSize</th>
      <th>trim_method</th>
      <th>cluster_method</th>
      <th>pid_clustering</th>
      <th>taxa_class_method</th>
      <th>taxa_ref_db</th>
      <th>code_repo</th>
      <th>identificationReferences</th>
      <th>controls_used</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>16S V4-V5</td>
      <td>411</td>
      <td>cutadapt</td>
      <td>Tourmaline; qiime2-2021.2; dada2</td>
      <td>ASV</td>
      <td>Tourmaline; qiime2-2021.2; naive-bayes classifier</td>
      <td>Silva SSU Ref NR 99 v138.1; 515f-926r region; 10.5281/zenodo.8392695</td>
      <td>https://github.com/aomlomics/gomecc</td>
      <td>10.5281/zenodo.8392695 | https://github.com/aomlomics/tourmaline</td>
      <td>12 distilled water blanks | 2 PCR no-template controls | 7 extraction blanks | 12 2nd PCR no-template controls | 3 Zymo mock community</td>
    </tr>
    <tr>
      <th>1</th>
      <td>18S V9</td>
      <td>260</td>
      <td>cutadapt</td>
      <td>Tourmaline; qiime2-2021.2; dada2</td>
      <td>ASV</td>
      <td>Tourmaline; qiime2-2021.2; naive-bayes classifier</td>
      <td>PR2 v5.0.1; V9 1391f-1510r region; 10.5281/zenodo.8392706</td>
      <td>https://github.com/aomlomics/gomecc</td>
      <td>10.5281/zenodo.8392706 | https://pr2-database.org/ | https://github.com/aomlomics/tourmaline</td>
      <td>12 distilled water blanks | 2 PCR no-template controls | 7 extraction blanks | 7 2nd PCR no-template controls</td>
    </tr>
  </tbody>
</table>
</div>




```python
data['analysis_data'].head()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>amplicon_sequenced</th>
      <th>ampliconSize</th>
      <th>trim_method</th>
      <th>cluster_method</th>
      <th>pid_clustering</th>
      <th>taxa_class_method</th>
      <th>taxa_ref_db</th>
      <th>code_repo</th>
      <th>identificationReferences</th>
      <th>controls_used</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>16S V4-V5</td>
      <td>411</td>
      <td>cutadapt</td>
      <td>Tourmaline; qiime2-2021.2; dada2</td>
      <td>ASV</td>
      <td>Tourmaline; qiime2-2021.2; naive-bayes classifier</td>
      <td>Silva SSU Ref NR 99 v138.1; 515f-926r region; 10.5281/zenodo.8392695</td>
      <td>https://github.com/aomlomics/gomecc</td>
      <td>10.5281/zenodo.8392695 | https://github.com/aomlomics/tourmaline</td>
      <td>12 distilled water blanks | 2 PCR no-template controls | 7 extraction blanks | 12 2nd PCR no-template controls | 3 Zymo mock community</td>
    </tr>
    <tr>
      <th>1</th>
      <td>18S V9</td>
      <td>260</td>
      <td>cutadapt</td>
      <td>Tourmaline; qiime2-2021.2; dada2</td>
      <td>ASV</td>
      <td>Tourmaline; qiime2-2021.2; naive-bayes classifier</td>
      <td>PR2 v5.0.1; V9 1391f-1510r region; 10.5281/zenodo.8392706</td>
      <td>https://github.com/aomlomics/gomecc</td>
      <td>10.5281/zenodo.8392706 | https://pr2-database.org/ | https://github.com/aomlomics/tourmaline</td>
      <td>12 distilled water blanks | 2 PCR no-template controls | 7 extraction blanks | 7 2nd PCR no-template controls</td>
    </tr>
  </tbody>
</table>
</div>




```python
occ16_test['taxa_class_method'] = data['analysis_data'].loc[data['analysis_data']['amplicon_sequenced'] == '16S V4-V5','taxa_class_method'].item()
occ16_test['taxa_ref_db'] = data['analysis_data'].loc[data['analysis_data']['amplicon_sequenced'] == '16S V4-V5','taxa_ref_db'].item()

occ16_test['identificationRemarks'] = occ16_test['taxa_class_method'] +", confidence (at lowest specified taxon): "+occ16_test['Confidence'].astype(str) +", against reference database: "+occ16_test['taxa_ref_db']
```


```python
occ16_test['identificationRemarks'][0]
```




    'Tourmaline; qiime2-2021.2; naive-bayes classifier, confidence (at lowest specified taxon): 0.832189583, against reference database: Silva SSU Ref NR 99 v138.1; 515f-926r region; 10.5281/zenodo.8392695'



#### taxonID, basisOfRecord, eventID, nameAccordingTo, organismQuantityType


```python
occ16_test['taxonID'] = 'ASV:'+occ16_test['featureid']
occ16_test['basisOfRecord'] = 'MaterialSample'
occ16_test['nameAccordingTo'] = "WoRMS"
occ16_test['organismQuantityType'] = "DNA sequence reads"
occ16_test['recordedBy'] = data['study_data']['recordedBy'].values[0]
```

#### associatedSequences, materialSampleID


```python
data['prep_data'].columns
```




    Index(['sample_name', 'library_id', 'title', 'library_strategy',
           'library_source', 'library_selection', 'lib_layout', 'platform',
           'instrument_model', 'design_description', 'filetype', 'filename',
           'filename2', 'biosample_accession', 'sra_accession', 'seq_method',
           'nucl_acid_ext', 'amplicon_sequenced', 'target_gene',
           'target_subfragment', 'pcr_primer_forward', 'pcr_primer_reverse',
           'pcr_primer_name_forward', 'pcr_primer_name_reverse',
           'pcr_primer_reference', 'pcr_cond', 'nucl_acid_amp', 'adapters',
           'mid_barcode'],
          dtype='object')




```python
occ16_test = occ16_test.merge(data['prep_data'].loc[data['prep_data']['amplicon_sequenced'] == '16S V4-V5',['sample_name','sra_accession','biosample_accession']], how='left', left_on ='eventID', right_on='sample_name')
```


```python
occ16_test.head()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>featureid</th>
      <th>sequence</th>
      <th>taxonomy</th>
      <th>Confidence</th>
      <th>eventID</th>
      <th>organismQuantity</th>
      <th>occurrenceID</th>
      <th>verbatimIdentification</th>
      <th>scientificName</th>
      <th>scientificNameID</th>
      <th>kingdom</th>
      <th>phylum</th>
      <th>class</th>
      <th>order</th>
      <th>family</th>
      <th>genus</th>
      <th>taxonRank</th>
      <th>taxa_class_method</th>
      <th>taxa_ref_db</th>
      <th>identificationRemarks</th>
      <th>taxonID</th>
      <th>basisOfRecord</th>
      <th>nameAccordingTo</th>
      <th>organismQuantityType</th>
      <th>recordedBy</th>
      <th>sample_name</th>
      <th>sra_accession</th>
      <th>biosample_accession</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>00c4c1c65d8669ed9f07abe149f9a01d</td>
      <td>TACGGAGGGGGCTAACGTTGTTCGGAATTACTGGGCGTAAAGCGCGCGTAGGCGGATTAGACAGTTGAGGGTGAAATCCCGGAGCTTAACTTCGGAACTGCCCCCAATACTACTAATCTAGAGTTCGGAAGAGGTGAGTGGAATTC...</td>
      <td>d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Parvibaculales; f__OCS116_clade; g__OCS116_clade; s__uncultured_marine</td>
      <td>0.832190</td>
      <td>GOMECC4_27N_Sta1_DCM_A</td>
      <td>18</td>
      <td>GOMECC4_27N_Sta1_DCM_A_16S_occ00c4c1c65d8669ed9f07abe149f9a01d</td>
      <td>d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Parvibaculales; f__OCS116_clade; g__OCS116_clade; s__uncultured_marine</td>
      <td>Alphaproteobacteria</td>
      <td>urn:lsid:marinespecies.org:taxname:392750</td>
      <td>Bacteria</td>
      <td>Proteobacteria</td>
      <td>Alphaproteobacteria</td>
      <td>None</td>
      <td>None</td>
      <td>None</td>
      <td>Class</td>
      <td>Tourmaline; qiime2-2021.2; naive-bayes classifier</td>
      <td>Silva SSU Ref NR 99 v138.1; 515f-926r region; 10.5281/zenodo.8392695</td>
      <td>Tourmaline; qiime2-2021.2; naive-bayes classifier, confidence (at lowest specified taxon): 0.832189583, against reference database: Silva SSU Ref ...</td>
      <td>ASV:00c4c1c65d8669ed9f07abe149f9a01d</td>
      <td>MaterialSample</td>
      <td>WoRMS</td>
      <td>DNA sequence reads</td>
      <td>Luke Thompson | Katherine Silliman</td>
      <td>GOMECC4_27N_Sta1_DCM_A</td>
      <td>SRR26148187</td>
      <td>SAMN37516094</td>
    </tr>
    <tr>
      <th>1</th>
      <td>00e6c13fe86364a5084987093afa1916</td>
      <td>TACGAAGGGGGCGAGCGTTGTTCGGAATTACTGGGCGTAAAGGGCGCGTAGGCGGCTCTTTAAGTTAGGCGTGAAAGCCCCGGGCTCAACCTGGGAACTGCGCTTAAGACTGGAGAGCTAGAAAACGGAAGAGGGTAGTGGAATTC...</td>
      <td>d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Puniceispirillales; f__SAR116_clade; g__SAR116_clade</td>
      <td>0.867040</td>
      <td>GOMECC4_27N_Sta1_DCM_A</td>
      <td>36</td>
      <td>GOMECC4_27N_Sta1_DCM_A_16S_occ00e6c13fe86364a5084987093afa1916</td>
      <td>d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Puniceispirillales; f__SAR116_clade; g__SAR116_clade</td>
      <td>Alphaproteobacteria</td>
      <td>urn:lsid:marinespecies.org:taxname:392750</td>
      <td>Bacteria</td>
      <td>Proteobacteria</td>
      <td>Alphaproteobacteria</td>
      <td>None</td>
      <td>None</td>
      <td>None</td>
      <td>Class</td>
      <td>Tourmaline; qiime2-2021.2; naive-bayes classifier</td>
      <td>Silva SSU Ref NR 99 v138.1; 515f-926r region; 10.5281/zenodo.8392695</td>
      <td>Tourmaline; qiime2-2021.2; naive-bayes classifier, confidence (at lowest specified taxon): 0.867040054, against reference database: Silva SSU Ref ...</td>
      <td>ASV:00e6c13fe86364a5084987093afa1916</td>
      <td>MaterialSample</td>
      <td>WoRMS</td>
      <td>DNA sequence reads</td>
      <td>Luke Thompson | Katherine Silliman</td>
      <td>GOMECC4_27N_Sta1_DCM_A</td>
      <td>SRR26148187</td>
      <td>SAMN37516094</td>
    </tr>
    <tr>
      <th>2</th>
      <td>015dad1fafca90944d905beb2a980bc3</td>
      <td>TACCGGCGCCTCAAGTGGTAGTCGCTTTTATTGGGCCTAAAACGTCCGTAGCCGGTCTGGTACATTCGTGGGTAAATCAACTCGCTTAACGAGTTGAATTCTGCGAGGACGGCCAGACTTGGGACCGGGAGAGGTGTGGGGTACTC...</td>
      <td>d__Archaea; p__Thermoplasmatota; c__Thermoplasmata; o__Marine_Group_II; f__Marine_Group_II; g__Marine_Group_II</td>
      <td>1.000000</td>
      <td>GOMECC4_27N_Sta1_DCM_A</td>
      <td>49</td>
      <td>GOMECC4_27N_Sta1_DCM_A_16S_occ015dad1fafca90944d905beb2a980bc3</td>
      <td>d__Archaea; p__Thermoplasmatota; c__Thermoplasmata; o__Marine_Group_II; f__Marine_Group_II; g__Marine_Group_II</td>
      <td>Thermoplasmata</td>
      <td>urn:lsid:marinespecies.org:taxname:416268</td>
      <td>Archaea</td>
      <td>Euryarchaeota</td>
      <td>Thermoplasmata</td>
      <td>None</td>
      <td>None</td>
      <td>None</td>
      <td>Class</td>
      <td>Tourmaline; qiime2-2021.2; naive-bayes classifier</td>
      <td>Silva SSU Ref NR 99 v138.1; 515f-926r region; 10.5281/zenodo.8392695</td>
      <td>Tourmaline; qiime2-2021.2; naive-bayes classifier, confidence (at lowest specified taxon): 1.0, against reference database: Silva SSU Ref NR 99 v1...</td>
      <td>ASV:015dad1fafca90944d905beb2a980bc3</td>
      <td>MaterialSample</td>
      <td>WoRMS</td>
      <td>DNA sequence reads</td>
      <td>Luke Thompson | Katherine Silliman</td>
      <td>GOMECC4_27N_Sta1_DCM_A</td>
      <td>SRR26148187</td>
      <td>SAMN37516094</td>
    </tr>
    <tr>
      <th>3</th>
      <td>019c88c6ade406f731954f38e3461564</td>
      <td>TACAGGAGGGACGAGTGTTACTCGGAATGATTAGGCGTAAAGGGTCATTTAAGCGGTCCGATAAGTTAAAAGCCAACAGTTAGAGCCTAACTCTTTCAAGCTTTTAATACTGTCAGACTAGAGTATATCAGAGAATAGTAGAATTC...</td>
      <td>d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rickettsiales; f__Mitochondria; g__Mitochondria; s__uncultured_bacterium</td>
      <td>0.952911</td>
      <td>GOMECC4_27N_Sta1_DCM_A</td>
      <td>2</td>
      <td>GOMECC4_27N_Sta1_DCM_A_16S_occ019c88c6ade406f731954f38e3461564</td>
      <td>d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rickettsiales; f__Mitochondria; g__Mitochondria; s__uncultured_bacterium</td>
      <td>Rickettsiales</td>
      <td>urn:lsid:marinespecies.org:taxname:570969</td>
      <td>Bacteria</td>
      <td>Proteobacteria</td>
      <td>Alphaproteobacteria</td>
      <td>Rickettsiales</td>
      <td>None</td>
      <td>None</td>
      <td>Order</td>
      <td>Tourmaline; qiime2-2021.2; naive-bayes classifier</td>
      <td>Silva SSU Ref NR 99 v138.1; 515f-926r region; 10.5281/zenodo.8392695</td>
      <td>Tourmaline; qiime2-2021.2; naive-bayes classifier, confidence (at lowest specified taxon): 0.952910602, against reference database: Silva SSU Ref ...</td>
      <td>ASV:019c88c6ade406f731954f38e3461564</td>
      <td>MaterialSample</td>
      <td>WoRMS</td>
      <td>DNA sequence reads</td>
      <td>Luke Thompson | Katherine Silliman</td>
      <td>GOMECC4_27N_Sta1_DCM_A</td>
      <td>SRR26148187</td>
      <td>SAMN37516094</td>
    </tr>
    <tr>
      <th>4</th>
      <td>02dfb0869af4bf549d290d48e66e2196</td>
      <td>TACGAGGGGTGCTAGCGTTGTCCGGAATAACTGGGCGTAAAGGGTCCGTAGGCGTTTTGCTAAGTTGATCGTTAAATCCATCGGCTTAACCGATGACATGCGATCAAAACTGGCAGAATAGAATATGTGAGGGGAATGTAGAATTC...</td>
      <td>d__Bacteria; p__Marinimicrobia_(SAR406_clade); c__Marinimicrobia_(SAR406_clade); o__Marinimicrobia_(SAR406_clade); f__Marinimicrobia_(SAR406_clade...</td>
      <td>0.818195</td>
      <td>GOMECC4_27N_Sta1_DCM_A</td>
      <td>3</td>
      <td>GOMECC4_27N_Sta1_DCM_A_16S_occ02dfb0869af4bf549d290d48e66e2196</td>
      <td>d__Bacteria; p__Marinimicrobia_(SAR406_clade); c__Marinimicrobia_(SAR406_clade); o__Marinimicrobia_(SAR406_clade); f__Marinimicrobia_(SAR406_clade...</td>
      <td>Bacteria</td>
      <td>urn:lsid:marinespecies.org:taxname:6</td>
      <td>Bacteria</td>
      <td>None</td>
      <td>None</td>
      <td>None</td>
      <td>None</td>
      <td>None</td>
      <td>Kingdom</td>
      <td>Tourmaline; qiime2-2021.2; naive-bayes classifier</td>
      <td>Silva SSU Ref NR 99 v138.1; 515f-926r region; 10.5281/zenodo.8392695</td>
      <td>Tourmaline; qiime2-2021.2; naive-bayes classifier, confidence (at lowest specified taxon): 0.818195053, against reference database: Silva SSU Ref ...</td>
      <td>ASV:02dfb0869af4bf549d290d48e66e2196</td>
      <td>MaterialSample</td>
      <td>WoRMS</td>
      <td>DNA sequence reads</td>
      <td>Luke Thompson | Katherine Silliman</td>
      <td>GOMECC4_27N_Sta1_DCM_A</td>
      <td>SRR26148187</td>
      <td>SAMN37516094</td>
    </tr>
  </tbody>
</table>
</div>



#### eventID


```python
occ16_test['eventID'] = occ16_test['eventID']+"_16S"
```

#### sampleSize 


```python
# get sampleSize by total number of reads per sample
x = asv_tables['16S V4-V5'].sum(numeric_only=True).astype('int')
x.index = x.index+"_16S"
occ16_test['sampleSizeValue'] = occ16_test['eventID'].map(x).astype('str')
occ16_test['sampleSizeUnit'] = 'DNA sequence reads'
```


```python
occ16_test.head()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>featureid</th>
      <th>sequence</th>
      <th>taxonomy</th>
      <th>Confidence</th>
      <th>eventID</th>
      <th>organismQuantity</th>
      <th>occurrenceID</th>
      <th>verbatimIdentification</th>
      <th>scientificName</th>
      <th>scientificNameID</th>
      <th>kingdom</th>
      <th>phylum</th>
      <th>class</th>
      <th>order</th>
      <th>family</th>
      <th>genus</th>
      <th>taxonRank</th>
      <th>taxa_class_method</th>
      <th>taxa_ref_db</th>
      <th>identificationRemarks</th>
      <th>taxonID</th>
      <th>basisOfRecord</th>
      <th>nameAccordingTo</th>
      <th>organismQuantityType</th>
      <th>recordedBy</th>
      <th>sample_name</th>
      <th>sra_accession</th>
      <th>biosample_accession</th>
      <th>sampleSizeValue</th>
      <th>sampleSizeUnit</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>00c4c1c65d8669ed9f07abe149f9a01d</td>
      <td>TACGGAGGGGGCTAACGTTGTTCGGAATTACTGGGCGTAAAGCGCGCGTAGGCGGATTAGACAGTTGAGGGTGAAATCCCGGAGCTTAACTTCGGAACTGCCCCCAATACTACTAATCTAGAGTTCGGAAGAGGTGAGTGGAATTC...</td>
      <td>d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Parvibaculales; f__OCS116_clade; g__OCS116_clade; s__uncultured_marine</td>
      <td>0.832190</td>
      <td>GOMECC4_27N_Sta1_DCM_A_16S</td>
      <td>18</td>
      <td>GOMECC4_27N_Sta1_DCM_A_16S_occ00c4c1c65d8669ed9f07abe149f9a01d</td>
      <td>d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Parvibaculales; f__OCS116_clade; g__OCS116_clade; s__uncultured_marine</td>
      <td>Alphaproteobacteria</td>
      <td>urn:lsid:marinespecies.org:taxname:392750</td>
      <td>Bacteria</td>
      <td>Proteobacteria</td>
      <td>Alphaproteobacteria</td>
      <td>None</td>
      <td>None</td>
      <td>None</td>
      <td>Class</td>
      <td>Tourmaline; qiime2-2021.2; naive-bayes classifier</td>
      <td>Silva SSU Ref NR 99 v138.1; 515f-926r region; 10.5281/zenodo.8392695</td>
      <td>Tourmaline; qiime2-2021.2; naive-bayes classifier, confidence (at lowest specified taxon): 0.832189583, against reference database: Silva SSU Ref ...</td>
      <td>ASV:00c4c1c65d8669ed9f07abe149f9a01d</td>
      <td>MaterialSample</td>
      <td>WoRMS</td>
      <td>DNA sequence reads</td>
      <td>Luke Thompson | Katherine Silliman</td>
      <td>GOMECC4_27N_Sta1_DCM_A</td>
      <td>SRR26148187</td>
      <td>SAMN37516094</td>
      <td>16187</td>
      <td>DNA sequence reads</td>
    </tr>
    <tr>
      <th>1</th>
      <td>00e6c13fe86364a5084987093afa1916</td>
      <td>TACGAAGGGGGCGAGCGTTGTTCGGAATTACTGGGCGTAAAGGGCGCGTAGGCGGCTCTTTAAGTTAGGCGTGAAAGCCCCGGGCTCAACCTGGGAACTGCGCTTAAGACTGGAGAGCTAGAAAACGGAAGAGGGTAGTGGAATTC...</td>
      <td>d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Puniceispirillales; f__SAR116_clade; g__SAR116_clade</td>
      <td>0.867040</td>
      <td>GOMECC4_27N_Sta1_DCM_A_16S</td>
      <td>36</td>
      <td>GOMECC4_27N_Sta1_DCM_A_16S_occ00e6c13fe86364a5084987093afa1916</td>
      <td>d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Puniceispirillales; f__SAR116_clade; g__SAR116_clade</td>
      <td>Alphaproteobacteria</td>
      <td>urn:lsid:marinespecies.org:taxname:392750</td>
      <td>Bacteria</td>
      <td>Proteobacteria</td>
      <td>Alphaproteobacteria</td>
      <td>None</td>
      <td>None</td>
      <td>None</td>
      <td>Class</td>
      <td>Tourmaline; qiime2-2021.2; naive-bayes classifier</td>
      <td>Silva SSU Ref NR 99 v138.1; 515f-926r region; 10.5281/zenodo.8392695</td>
      <td>Tourmaline; qiime2-2021.2; naive-bayes classifier, confidence (at lowest specified taxon): 0.867040054, against reference database: Silva SSU Ref ...</td>
      <td>ASV:00e6c13fe86364a5084987093afa1916</td>
      <td>MaterialSample</td>
      <td>WoRMS</td>
      <td>DNA sequence reads</td>
      <td>Luke Thompson | Katherine Silliman</td>
      <td>GOMECC4_27N_Sta1_DCM_A</td>
      <td>SRR26148187</td>
      <td>SAMN37516094</td>
      <td>16187</td>
      <td>DNA sequence reads</td>
    </tr>
    <tr>
      <th>2</th>
      <td>015dad1fafca90944d905beb2a980bc3</td>
      <td>TACCGGCGCCTCAAGTGGTAGTCGCTTTTATTGGGCCTAAAACGTCCGTAGCCGGTCTGGTACATTCGTGGGTAAATCAACTCGCTTAACGAGTTGAATTCTGCGAGGACGGCCAGACTTGGGACCGGGAGAGGTGTGGGGTACTC...</td>
      <td>d__Archaea; p__Thermoplasmatota; c__Thermoplasmata; o__Marine_Group_II; f__Marine_Group_II; g__Marine_Group_II</td>
      <td>1.000000</td>
      <td>GOMECC4_27N_Sta1_DCM_A_16S</td>
      <td>49</td>
      <td>GOMECC4_27N_Sta1_DCM_A_16S_occ015dad1fafca90944d905beb2a980bc3</td>
      <td>d__Archaea; p__Thermoplasmatota; c__Thermoplasmata; o__Marine_Group_II; f__Marine_Group_II; g__Marine_Group_II</td>
      <td>Thermoplasmata</td>
      <td>urn:lsid:marinespecies.org:taxname:416268</td>
      <td>Archaea</td>
      <td>Euryarchaeota</td>
      <td>Thermoplasmata</td>
      <td>None</td>
      <td>None</td>
      <td>None</td>
      <td>Class</td>
      <td>Tourmaline; qiime2-2021.2; naive-bayes classifier</td>
      <td>Silva SSU Ref NR 99 v138.1; 515f-926r region; 10.5281/zenodo.8392695</td>
      <td>Tourmaline; qiime2-2021.2; naive-bayes classifier, confidence (at lowest specified taxon): 1.0, against reference database: Silva SSU Ref NR 99 v1...</td>
      <td>ASV:015dad1fafca90944d905beb2a980bc3</td>
      <td>MaterialSample</td>
      <td>WoRMS</td>
      <td>DNA sequence reads</td>
      <td>Luke Thompson | Katherine Silliman</td>
      <td>GOMECC4_27N_Sta1_DCM_A</td>
      <td>SRR26148187</td>
      <td>SAMN37516094</td>
      <td>16187</td>
      <td>DNA sequence reads</td>
    </tr>
    <tr>
      <th>3</th>
      <td>019c88c6ade406f731954f38e3461564</td>
      <td>TACAGGAGGGACGAGTGTTACTCGGAATGATTAGGCGTAAAGGGTCATTTAAGCGGTCCGATAAGTTAAAAGCCAACAGTTAGAGCCTAACTCTTTCAAGCTTTTAATACTGTCAGACTAGAGTATATCAGAGAATAGTAGAATTC...</td>
      <td>d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rickettsiales; f__Mitochondria; g__Mitochondria; s__uncultured_bacterium</td>
      <td>0.952911</td>
      <td>GOMECC4_27N_Sta1_DCM_A_16S</td>
      <td>2</td>
      <td>GOMECC4_27N_Sta1_DCM_A_16S_occ019c88c6ade406f731954f38e3461564</td>
      <td>d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rickettsiales; f__Mitochondria; g__Mitochondria; s__uncultured_bacterium</td>
      <td>Rickettsiales</td>
      <td>urn:lsid:marinespecies.org:taxname:570969</td>
      <td>Bacteria</td>
      <td>Proteobacteria</td>
      <td>Alphaproteobacteria</td>
      <td>Rickettsiales</td>
      <td>None</td>
      <td>None</td>
      <td>Order</td>
      <td>Tourmaline; qiime2-2021.2; naive-bayes classifier</td>
      <td>Silva SSU Ref NR 99 v138.1; 515f-926r region; 10.5281/zenodo.8392695</td>
      <td>Tourmaline; qiime2-2021.2; naive-bayes classifier, confidence (at lowest specified taxon): 0.952910602, against reference database: Silva SSU Ref ...</td>
      <td>ASV:019c88c6ade406f731954f38e3461564</td>
      <td>MaterialSample</td>
      <td>WoRMS</td>
      <td>DNA sequence reads</td>
      <td>Luke Thompson | Katherine Silliman</td>
      <td>GOMECC4_27N_Sta1_DCM_A</td>
      <td>SRR26148187</td>
      <td>SAMN37516094</td>
      <td>16187</td>
      <td>DNA sequence reads</td>
    </tr>
    <tr>
      <th>4</th>
      <td>02dfb0869af4bf549d290d48e66e2196</td>
      <td>TACGAGGGGTGCTAGCGTTGTCCGGAATAACTGGGCGTAAAGGGTCCGTAGGCGTTTTGCTAAGTTGATCGTTAAATCCATCGGCTTAACCGATGACATGCGATCAAAACTGGCAGAATAGAATATGTGAGGGGAATGTAGAATTC...</td>
      <td>d__Bacteria; p__Marinimicrobia_(SAR406_clade); c__Marinimicrobia_(SAR406_clade); o__Marinimicrobia_(SAR406_clade); f__Marinimicrobia_(SAR406_clade...</td>
      <td>0.818195</td>
      <td>GOMECC4_27N_Sta1_DCM_A_16S</td>
      <td>3</td>
      <td>GOMECC4_27N_Sta1_DCM_A_16S_occ02dfb0869af4bf549d290d48e66e2196</td>
      <td>d__Bacteria; p__Marinimicrobia_(SAR406_clade); c__Marinimicrobia_(SAR406_clade); o__Marinimicrobia_(SAR406_clade); f__Marinimicrobia_(SAR406_clade...</td>
      <td>Bacteria</td>
      <td>urn:lsid:marinespecies.org:taxname:6</td>
      <td>Bacteria</td>
      <td>None</td>
      <td>None</td>
      <td>None</td>
      <td>None</td>
      <td>None</td>
      <td>Kingdom</td>
      <td>Tourmaline; qiime2-2021.2; naive-bayes classifier</td>
      <td>Silva SSU Ref NR 99 v138.1; 515f-926r region; 10.5281/zenodo.8392695</td>
      <td>Tourmaline; qiime2-2021.2; naive-bayes classifier, confidence (at lowest specified taxon): 0.818195053, against reference database: Silva SSU Ref ...</td>
      <td>ASV:02dfb0869af4bf549d290d48e66e2196</td>
      <td>MaterialSample</td>
      <td>WoRMS</td>
      <td>DNA sequence reads</td>
      <td>Luke Thompson | Katherine Silliman</td>
      <td>GOMECC4_27N_Sta1_DCM_A</td>
      <td>SRR26148187</td>
      <td>SAMN37516094</td>
      <td>16187</td>
      <td>DNA sequence reads</td>
    </tr>
  </tbody>
</table>
</div>




```python
# drop unnneeded columns
occ16_test.drop(columns=['sample_name','featureid','taxonomy','Confidence','taxa_class_method','taxa_ref_db'],inplace=True)
```


```python
occ16_test['associatedSequences'] = occ16_test['sra_accession']+' | '+ occ16_test['biosample_accession']+' | '+data['study_data']['bioproject_accession'].values[0]
```


```python
occ16_test.rename(columns={'biosample_accession': 'materialSampleID',
                  'sequence': 'DNA_sequence'},inplace=True)
                   
```


```python
# drop unnneeded columns
occ16_test.drop(columns=['sra_accession'],inplace=True)
```


```python
occ16_test.columns
```




    Index(['DNA_sequence', 'eventID', 'organismQuantity', 'occurrenceID',
           'verbatimIdentification', 'scientificName', 'scientificNameID',
           'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'taxonRank',
           'identificationRemarks', 'taxonID', 'basisOfRecord', 'nameAccordingTo',
           'organismQuantityType', 'recordedBy', 'materialSampleID',
           'sampleSizeValue', 'sampleSizeUnit', 'associatedSequences'],
          dtype='object')




```python
occ16_test.head()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>DNA_sequence</th>
      <th>eventID</th>
      <th>organismQuantity</th>
      <th>occurrenceID</th>
      <th>verbatimIdentification</th>
      <th>scientificName</th>
      <th>scientificNameID</th>
      <th>kingdom</th>
      <th>phylum</th>
      <th>class</th>
      <th>order</th>
      <th>family</th>
      <th>genus</th>
      <th>taxonRank</th>
      <th>identificationRemarks</th>
      <th>taxonID</th>
      <th>basisOfRecord</th>
      <th>nameAccordingTo</th>
      <th>organismQuantityType</th>
      <th>recordedBy</th>
      <th>materialSampleID</th>
      <th>sampleSizeValue</th>
      <th>sampleSizeUnit</th>
      <th>associatedSequences</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>TACGGAGGGGGCTAACGTTGTTCGGAATTACTGGGCGTAAAGCGCGCGTAGGCGGATTAGACAGTTGAGGGTGAAATCCCGGAGCTTAACTTCGGAACTGCCCCCAATACTACTAATCTAGAGTTCGGAAGAGGTGAGTGGAATTC...</td>
      <td>GOMECC4_27N_Sta1_DCM_A_16S</td>
      <td>18</td>
      <td>GOMECC4_27N_Sta1_DCM_A_16S_occ00c4c1c65d8669ed9f07abe149f9a01d</td>
      <td>d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Parvibaculales; f__OCS116_clade; g__OCS116_clade; s__uncultured_marine</td>
      <td>Alphaproteobacteria</td>
      <td>urn:lsid:marinespecies.org:taxname:392750</td>
      <td>Bacteria</td>
      <td>Proteobacteria</td>
      <td>Alphaproteobacteria</td>
      <td>None</td>
      <td>None</td>
      <td>None</td>
      <td>Class</td>
      <td>Tourmaline; qiime2-2021.2; naive-bayes classifier, confidence (at lowest specified taxon): 0.832189583, against reference database: Silva SSU Ref ...</td>
      <td>ASV:00c4c1c65d8669ed9f07abe149f9a01d</td>
      <td>MaterialSample</td>
      <td>WoRMS</td>
      <td>DNA sequence reads</td>
      <td>Luke Thompson | Katherine Silliman</td>
      <td>SAMN37516094</td>
      <td>16187</td>
      <td>DNA sequence reads</td>
      <td>SRR26148187 | SAMN37516094 | PRJNA887898</td>
    </tr>
    <tr>
      <th>1</th>
      <td>TACGAAGGGGGCGAGCGTTGTTCGGAATTACTGGGCGTAAAGGGCGCGTAGGCGGCTCTTTAAGTTAGGCGTGAAAGCCCCGGGCTCAACCTGGGAACTGCGCTTAAGACTGGAGAGCTAGAAAACGGAAGAGGGTAGTGGAATTC...</td>
      <td>GOMECC4_27N_Sta1_DCM_A_16S</td>
      <td>36</td>
      <td>GOMECC4_27N_Sta1_DCM_A_16S_occ00e6c13fe86364a5084987093afa1916</td>
      <td>d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Puniceispirillales; f__SAR116_clade; g__SAR116_clade</td>
      <td>Alphaproteobacteria</td>
      <td>urn:lsid:marinespecies.org:taxname:392750</td>
      <td>Bacteria</td>
      <td>Proteobacteria</td>
      <td>Alphaproteobacteria</td>
      <td>None</td>
      <td>None</td>
      <td>None</td>
      <td>Class</td>
      <td>Tourmaline; qiime2-2021.2; naive-bayes classifier, confidence (at lowest specified taxon): 0.867040054, against reference database: Silva SSU Ref ...</td>
      <td>ASV:00e6c13fe86364a5084987093afa1916</td>
      <td>MaterialSample</td>
      <td>WoRMS</td>
      <td>DNA sequence reads</td>
      <td>Luke Thompson | Katherine Silliman</td>
      <td>SAMN37516094</td>
      <td>16187</td>
      <td>DNA sequence reads</td>
      <td>SRR26148187 | SAMN37516094 | PRJNA887898</td>
    </tr>
    <tr>
      <th>2</th>
      <td>TACCGGCGCCTCAAGTGGTAGTCGCTTTTATTGGGCCTAAAACGTCCGTAGCCGGTCTGGTACATTCGTGGGTAAATCAACTCGCTTAACGAGTTGAATTCTGCGAGGACGGCCAGACTTGGGACCGGGAGAGGTGTGGGGTACTC...</td>
      <td>GOMECC4_27N_Sta1_DCM_A_16S</td>
      <td>49</td>
      <td>GOMECC4_27N_Sta1_DCM_A_16S_occ015dad1fafca90944d905beb2a980bc3</td>
      <td>d__Archaea; p__Thermoplasmatota; c__Thermoplasmata; o__Marine_Group_II; f__Marine_Group_II; g__Marine_Group_II</td>
      <td>Thermoplasmata</td>
      <td>urn:lsid:marinespecies.org:taxname:416268</td>
      <td>Archaea</td>
      <td>Euryarchaeota</td>
      <td>Thermoplasmata</td>
      <td>None</td>
      <td>None</td>
      <td>None</td>
      <td>Class</td>
      <td>Tourmaline; qiime2-2021.2; naive-bayes classifier, confidence (at lowest specified taxon): 1.0, against reference database: Silva SSU Ref NR 99 v1...</td>
      <td>ASV:015dad1fafca90944d905beb2a980bc3</td>
      <td>MaterialSample</td>
      <td>WoRMS</td>
      <td>DNA sequence reads</td>
      <td>Luke Thompson | Katherine Silliman</td>
      <td>SAMN37516094</td>
      <td>16187</td>
      <td>DNA sequence reads</td>
      <td>SRR26148187 | SAMN37516094 | PRJNA887898</td>
    </tr>
    <tr>
      <th>3</th>
      <td>TACAGGAGGGACGAGTGTTACTCGGAATGATTAGGCGTAAAGGGTCATTTAAGCGGTCCGATAAGTTAAAAGCCAACAGTTAGAGCCTAACTCTTTCAAGCTTTTAATACTGTCAGACTAGAGTATATCAGAGAATAGTAGAATTC...</td>
      <td>GOMECC4_27N_Sta1_DCM_A_16S</td>
      <td>2</td>
      <td>GOMECC4_27N_Sta1_DCM_A_16S_occ019c88c6ade406f731954f38e3461564</td>
      <td>d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rickettsiales; f__Mitochondria; g__Mitochondria; s__uncultured_bacterium</td>
      <td>Rickettsiales</td>
      <td>urn:lsid:marinespecies.org:taxname:570969</td>
      <td>Bacteria</td>
      <td>Proteobacteria</td>
      <td>Alphaproteobacteria</td>
      <td>Rickettsiales</td>
      <td>None</td>
      <td>None</td>
      <td>Order</td>
      <td>Tourmaline; qiime2-2021.2; naive-bayes classifier, confidence (at lowest specified taxon): 0.952910602, against reference database: Silva SSU Ref ...</td>
      <td>ASV:019c88c6ade406f731954f38e3461564</td>
      <td>MaterialSample</td>
      <td>WoRMS</td>
      <td>DNA sequence reads</td>
      <td>Luke Thompson | Katherine Silliman</td>
      <td>SAMN37516094</td>
      <td>16187</td>
      <td>DNA sequence reads</td>
      <td>SRR26148187 | SAMN37516094 | PRJNA887898</td>
    </tr>
    <tr>
      <th>4</th>
      <td>TACGAGGGGTGCTAGCGTTGTCCGGAATAACTGGGCGTAAAGGGTCCGTAGGCGTTTTGCTAAGTTGATCGTTAAATCCATCGGCTTAACCGATGACATGCGATCAAAACTGGCAGAATAGAATATGTGAGGGGAATGTAGAATTC...</td>
      <td>GOMECC4_27N_Sta1_DCM_A_16S</td>
      <td>3</td>
      <td>GOMECC4_27N_Sta1_DCM_A_16S_occ02dfb0869af4bf549d290d48e66e2196</td>
      <td>d__Bacteria; p__Marinimicrobia_(SAR406_clade); c__Marinimicrobia_(SAR406_clade); o__Marinimicrobia_(SAR406_clade); f__Marinimicrobia_(SAR406_clade...</td>
      <td>Bacteria</td>
      <td>urn:lsid:marinespecies.org:taxname:6</td>
      <td>Bacteria</td>
      <td>None</td>
      <td>None</td>
      <td>None</td>
      <td>None</td>
      <td>None</td>
      <td>Kingdom</td>
      <td>Tourmaline; qiime2-2021.2; naive-bayes classifier, confidence (at lowest specified taxon): 0.818195053, against reference database: Silva SSU Ref ...</td>
      <td>ASV:02dfb0869af4bf549d290d48e66e2196</td>
      <td>MaterialSample</td>
      <td>WoRMS</td>
      <td>DNA sequence reads</td>
      <td>Luke Thompson | Katherine Silliman</td>
      <td>SAMN37516094</td>
      <td>16187</td>
      <td>DNA sequence reads</td>
      <td>SRR26148187 | SAMN37516094 | PRJNA887898</td>
    </tr>
  </tbody>
</table>
</div>



### merge event and occurrence


```python
all_event_data.tail()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>eventID</th>
      <th>locationID</th>
      <th>eventDate</th>
      <th>minimumDepthInMeters</th>
      <th>locality</th>
      <th>decimalLatitude</th>
      <th>decimalLongitude</th>
      <th>samplingProtocol</th>
      <th>waterBody</th>
      <th>maximumDepthInMeters</th>
      <th>parentEventID</th>
      <th>datasetID</th>
      <th>geodeticDatum</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>939</th>
      <td>GOMECC4_CAPECORAL_Sta141_DCM_B_18S</td>
      <td>CAPECORAL_Sta141</td>
      <td>2021-10-20T12:47-04:00</td>
      <td>59</td>
      <td>USA: Gulf of Mexico</td>
      <td>25.574</td>
      <td>-84.843</td>
      <td>CTD rosette</td>
      <td>Mexico, Gulf of</td>
      <td>59</td>
      <td>GOMECC4_CAPECORAL_Sta141_DCM_B</td>
      <td>noaa-aoml-gomecc4</td>
      <td>WGS84</td>
    </tr>
    <tr>
      <th>940</th>
      <td>GOMECC4_CAPECORAL_Sta141_DCM_C_18S</td>
      <td>CAPECORAL_Sta141</td>
      <td>2021-10-20T12:47-04:00</td>
      <td>59</td>
      <td>USA: Gulf of Mexico</td>
      <td>25.574</td>
      <td>-84.843</td>
      <td>CTD rosette</td>
      <td>Mexico, Gulf of</td>
      <td>59</td>
      <td>GOMECC4_CAPECORAL_Sta141_DCM_C</td>
      <td>noaa-aoml-gomecc4</td>
      <td>WGS84</td>
    </tr>
    <tr>
      <th>941</th>
      <td>GOMECC4_CAPECORAL_Sta141_Surface_A_18S</td>
      <td>CAPECORAL_Sta141</td>
      <td>2021-10-20T12:47-04:00</td>
      <td>4</td>
      <td>USA: Gulf of Mexico</td>
      <td>25.574</td>
      <td>-84.843</td>
      <td>CTD rosette</td>
      <td>Mexico, Gulf of</td>
      <td>4</td>
      <td>GOMECC4_CAPECORAL_Sta141_Surface_A</td>
      <td>noaa-aoml-gomecc4</td>
      <td>WGS84</td>
    </tr>
    <tr>
      <th>942</th>
      <td>GOMECC4_CAPECORAL_Sta141_Surface_B_18S</td>
      <td>CAPECORAL_Sta141</td>
      <td>2021-10-20T12:47-04:00</td>
      <td>4</td>
      <td>USA: Gulf of Mexico</td>
      <td>25.574</td>
      <td>-84.843</td>
      <td>CTD rosette</td>
      <td>Mexico, Gulf of</td>
      <td>4</td>
      <td>GOMECC4_CAPECORAL_Sta141_Surface_B</td>
      <td>noaa-aoml-gomecc4</td>
      <td>WGS84</td>
    </tr>
    <tr>
      <th>943</th>
      <td>GOMECC4_CAPECORAL_Sta141_Surface_C_18S</td>
      <td>CAPECORAL_Sta141</td>
      <td>2021-10-20T12:47-04:00</td>
      <td>4</td>
      <td>USA: Gulf of Mexico</td>
      <td>25.574</td>
      <td>-84.843</td>
      <td>CTD rosette</td>
      <td>Mexico, Gulf of</td>
      <td>4</td>
      <td>GOMECC4_CAPECORAL_Sta141_Surface_C</td>
      <td>noaa-aoml-gomecc4</td>
      <td>WGS84</td>
    </tr>
  </tbody>
</table>
</div>




```python
occ16_merged = occ16_test.merge(all_event_data,how='left',on='eventID')
```


```python
occ16_merged.head()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>DNA_sequence</th>
      <th>eventID</th>
      <th>organismQuantity</th>
      <th>occurrenceID</th>
      <th>verbatimIdentification</th>
      <th>scientificName</th>
      <th>scientificNameID</th>
      <th>kingdom</th>
      <th>phylum</th>
      <th>class</th>
      <th>order</th>
      <th>family</th>
      <th>genus</th>
      <th>taxonRank</th>
      <th>identificationRemarks</th>
      <th>taxonID</th>
      <th>basisOfRecord</th>
      <th>nameAccordingTo</th>
      <th>organismQuantityType</th>
      <th>recordedBy</th>
      <th>materialSampleID</th>
      <th>sampleSizeValue</th>
      <th>sampleSizeUnit</th>
      <th>associatedSequences</th>
      <th>locationID</th>
      <th>eventDate</th>
      <th>minimumDepthInMeters</th>
      <th>locality</th>
      <th>decimalLatitude</th>
      <th>decimalLongitude</th>
      <th>samplingProtocol</th>
      <th>waterBody</th>
      <th>maximumDepthInMeters</th>
      <th>parentEventID</th>
      <th>datasetID</th>
      <th>geodeticDatum</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>TACGGAGGGGGCTAACGTTGTTCGGAATTACTGGGCGTAAAGCGCGCGTAGGCGGATTAGACAGTTGAGGGTGAAATCCCGGAGCTTAACTTCGGAACTGCCCCCAATACTACTAATCTAGAGTTCGGAAGAGGTGAGTGGAATTC...</td>
      <td>GOMECC4_27N_Sta1_DCM_A_16S</td>
      <td>18</td>
      <td>GOMECC4_27N_Sta1_DCM_A_16S_occ00c4c1c65d8669ed9f07abe149f9a01d</td>
      <td>d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Parvibaculales; f__OCS116_clade; g__OCS116_clade; s__uncultured_marine</td>
      <td>Alphaproteobacteria</td>
      <td>urn:lsid:marinespecies.org:taxname:392750</td>
      <td>Bacteria</td>
      <td>Proteobacteria</td>
      <td>Alphaproteobacteria</td>
      <td>None</td>
      <td>None</td>
      <td>None</td>
      <td>Class</td>
      <td>Tourmaline; qiime2-2021.2; naive-bayes classifier, confidence (at lowest specified taxon): 0.832189583, against reference database: Silva SSU Ref ...</td>
      <td>ASV:00c4c1c65d8669ed9f07abe149f9a01d</td>
      <td>MaterialSample</td>
      <td>WoRMS</td>
      <td>DNA sequence reads</td>
      <td>Luke Thompson | Katherine Silliman</td>
      <td>SAMN37516094</td>
      <td>16187</td>
      <td>DNA sequence reads</td>
      <td>SRR26148187 | SAMN37516094 | PRJNA887898</td>
      <td>27N_Sta1</td>
      <td>2021-09-14T11:00-04:00</td>
      <td>49</td>
      <td>USA: Atlantic Ocean, east of Florida (27 N)</td>
      <td>26.997</td>
      <td>-79.618</td>
      <td>CTD rosette</td>
      <td>Atlantic Ocean</td>
      <td>49</td>
      <td>GOMECC4_27N_Sta1_DCM_A</td>
      <td>noaa-aoml-gomecc4</td>
      <td>WGS84</td>
    </tr>
    <tr>
      <th>1</th>
      <td>TACGAAGGGGGCGAGCGTTGTTCGGAATTACTGGGCGTAAAGGGCGCGTAGGCGGCTCTTTAAGTTAGGCGTGAAAGCCCCGGGCTCAACCTGGGAACTGCGCTTAAGACTGGAGAGCTAGAAAACGGAAGAGGGTAGTGGAATTC...</td>
      <td>GOMECC4_27N_Sta1_DCM_A_16S</td>
      <td>36</td>
      <td>GOMECC4_27N_Sta1_DCM_A_16S_occ00e6c13fe86364a5084987093afa1916</td>
      <td>d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Puniceispirillales; f__SAR116_clade; g__SAR116_clade</td>
      <td>Alphaproteobacteria</td>
      <td>urn:lsid:marinespecies.org:taxname:392750</td>
      <td>Bacteria</td>
      <td>Proteobacteria</td>
      <td>Alphaproteobacteria</td>
      <td>None</td>
      <td>None</td>
      <td>None</td>
      <td>Class</td>
      <td>Tourmaline; qiime2-2021.2; naive-bayes classifier, confidence (at lowest specified taxon): 0.867040054, against reference database: Silva SSU Ref ...</td>
      <td>ASV:00e6c13fe86364a5084987093afa1916</td>
      <td>MaterialSample</td>
      <td>WoRMS</td>
      <td>DNA sequence reads</td>
      <td>Luke Thompson | Katherine Silliman</td>
      <td>SAMN37516094</td>
      <td>16187</td>
      <td>DNA sequence reads</td>
      <td>SRR26148187 | SAMN37516094 | PRJNA887898</td>
      <td>27N_Sta1</td>
      <td>2021-09-14T11:00-04:00</td>
      <td>49</td>
      <td>USA: Atlantic Ocean, east of Florida (27 N)</td>
      <td>26.997</td>
      <td>-79.618</td>
      <td>CTD rosette</td>
      <td>Atlantic Ocean</td>
      <td>49</td>
      <td>GOMECC4_27N_Sta1_DCM_A</td>
      <td>noaa-aoml-gomecc4</td>
      <td>WGS84</td>
    </tr>
    <tr>
      <th>2</th>
      <td>TACCGGCGCCTCAAGTGGTAGTCGCTTTTATTGGGCCTAAAACGTCCGTAGCCGGTCTGGTACATTCGTGGGTAAATCAACTCGCTTAACGAGTTGAATTCTGCGAGGACGGCCAGACTTGGGACCGGGAGAGGTGTGGGGTACTC...</td>
      <td>GOMECC4_27N_Sta1_DCM_A_16S</td>
      <td>49</td>
      <td>GOMECC4_27N_Sta1_DCM_A_16S_occ015dad1fafca90944d905beb2a980bc3</td>
      <td>d__Archaea; p__Thermoplasmatota; c__Thermoplasmata; o__Marine_Group_II; f__Marine_Group_II; g__Marine_Group_II</td>
      <td>Thermoplasmata</td>
      <td>urn:lsid:marinespecies.org:taxname:416268</td>
      <td>Archaea</td>
      <td>Euryarchaeota</td>
      <td>Thermoplasmata</td>
      <td>None</td>
      <td>None</td>
      <td>None</td>
      <td>Class</td>
      <td>Tourmaline; qiime2-2021.2; naive-bayes classifier, confidence (at lowest specified taxon): 1.0, against reference database: Silva SSU Ref NR 99 v1...</td>
      <td>ASV:015dad1fafca90944d905beb2a980bc3</td>
      <td>MaterialSample</td>
      <td>WoRMS</td>
      <td>DNA sequence reads</td>
      <td>Luke Thompson | Katherine Silliman</td>
      <td>SAMN37516094</td>
      <td>16187</td>
      <td>DNA sequence reads</td>
      <td>SRR26148187 | SAMN37516094 | PRJNA887898</td>
      <td>27N_Sta1</td>
      <td>2021-09-14T11:00-04:00</td>
      <td>49</td>
      <td>USA: Atlantic Ocean, east of Florida (27 N)</td>
      <td>26.997</td>
      <td>-79.618</td>
      <td>CTD rosette</td>
      <td>Atlantic Ocean</td>
      <td>49</td>
      <td>GOMECC4_27N_Sta1_DCM_A</td>
      <td>noaa-aoml-gomecc4</td>
      <td>WGS84</td>
    </tr>
    <tr>
      <th>3</th>
      <td>TACAGGAGGGACGAGTGTTACTCGGAATGATTAGGCGTAAAGGGTCATTTAAGCGGTCCGATAAGTTAAAAGCCAACAGTTAGAGCCTAACTCTTTCAAGCTTTTAATACTGTCAGACTAGAGTATATCAGAGAATAGTAGAATTC...</td>
      <td>GOMECC4_27N_Sta1_DCM_A_16S</td>
      <td>2</td>
      <td>GOMECC4_27N_Sta1_DCM_A_16S_occ019c88c6ade406f731954f38e3461564</td>
      <td>d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rickettsiales; f__Mitochondria; g__Mitochondria; s__uncultured_bacterium</td>
      <td>Rickettsiales</td>
      <td>urn:lsid:marinespecies.org:taxname:570969</td>
      <td>Bacteria</td>
      <td>Proteobacteria</td>
      <td>Alphaproteobacteria</td>
      <td>Rickettsiales</td>
      <td>None</td>
      <td>None</td>
      <td>Order</td>
      <td>Tourmaline; qiime2-2021.2; naive-bayes classifier, confidence (at lowest specified taxon): 0.952910602, against reference database: Silva SSU Ref ...</td>
      <td>ASV:019c88c6ade406f731954f38e3461564</td>
      <td>MaterialSample</td>
      <td>WoRMS</td>
      <td>DNA sequence reads</td>
      <td>Luke Thompson | Katherine Silliman</td>
      <td>SAMN37516094</td>
      <td>16187</td>
      <td>DNA sequence reads</td>
      <td>SRR26148187 | SAMN37516094 | PRJNA887898</td>
      <td>27N_Sta1</td>
      <td>2021-09-14T11:00-04:00</td>
      <td>49</td>
      <td>USA: Atlantic Ocean, east of Florida (27 N)</td>
      <td>26.997</td>
      <td>-79.618</td>
      <td>CTD rosette</td>
      <td>Atlantic Ocean</td>
      <td>49</td>
      <td>GOMECC4_27N_Sta1_DCM_A</td>
      <td>noaa-aoml-gomecc4</td>
      <td>WGS84</td>
    </tr>
    <tr>
      <th>4</th>
      <td>TACGAGGGGTGCTAGCGTTGTCCGGAATAACTGGGCGTAAAGGGTCCGTAGGCGTTTTGCTAAGTTGATCGTTAAATCCATCGGCTTAACCGATGACATGCGATCAAAACTGGCAGAATAGAATATGTGAGGGGAATGTAGAATTC...</td>
      <td>GOMECC4_27N_Sta1_DCM_A_16S</td>
      <td>3</td>
      <td>GOMECC4_27N_Sta1_DCM_A_16S_occ02dfb0869af4bf549d290d48e66e2196</td>
      <td>d__Bacteria; p__Marinimicrobia_(SAR406_clade); c__Marinimicrobia_(SAR406_clade); o__Marinimicrobia_(SAR406_clade); f__Marinimicrobia_(SAR406_clade...</td>
      <td>Bacteria</td>
      <td>urn:lsid:marinespecies.org:taxname:6</td>
      <td>Bacteria</td>
      <td>None</td>
      <td>None</td>
      <td>None</td>
      <td>None</td>
      <td>None</td>
      <td>Kingdom</td>
      <td>Tourmaline; qiime2-2021.2; naive-bayes classifier, confidence (at lowest specified taxon): 0.818195053, against reference database: Silva SSU Ref ...</td>
      <td>ASV:02dfb0869af4bf549d290d48e66e2196</td>
      <td>MaterialSample</td>
      <td>WoRMS</td>
      <td>DNA sequence reads</td>
      <td>Luke Thompson | Katherine Silliman</td>
      <td>SAMN37516094</td>
      <td>16187</td>
      <td>DNA sequence reads</td>
      <td>SRR26148187 | SAMN37516094 | PRJNA887898</td>
      <td>27N_Sta1</td>
      <td>2021-09-14T11:00-04:00</td>
      <td>49</td>
      <td>USA: Atlantic Ocean, east of Florida (27 N)</td>
      <td>26.997</td>
      <td>-79.618</td>
      <td>CTD rosette</td>
      <td>Atlantic Ocean</td>
      <td>49</td>
      <td>GOMECC4_27N_Sta1_DCM_A</td>
      <td>noaa-aoml-gomecc4</td>
      <td>WGS84</td>
    </tr>
  </tbody>
</table>
</div>




```python
occ16_merged.drop(columns=['DNA_sequence']).to_csv("../processed/occurrence_16S.tsv",sep="\t",index=False)
```

### 18S worms

18S PR2 database provides WORMS IDs for species that are in worms. We will read in that file, assign known worms ids, the do a search for unannotated taxa.


```python
pr2_18S = pd.read_excel("../../../databases/18S_PR2/pr2_v5.0.0_SSU/pr2_version_5.0.0_taxonomy.xlsx",
    index_col=None, na_values=[""])
pr2_18S = pr2_18S.dropna(subset=['worms_id'])
pr2_18S['worms_id'] = pr2_18S['worms_id'].astype('int').astype('str')
pr2_18S['species'] = pr2_18S['species'].replace('_',' ',regex=True)
pr2_18S['species'] = pr2_18S['species'].replace(' sp\.','',regex=True)
pr2_18S['species'] = pr2_18S['species'].replace(' spp\.','',regex=True)
pr2_18S['species'] = pr2_18S['species'].replace('-',' ',regex=True)
pr2_18S['species'] = pr2_18S['species'].replace('\/',' ',regex=True)
```


```python
pr2_18S_dict = dict(zip(pr2_18S.species,pr2_18S.worms_id))

```


```python
(pr2_18S_dict['Aphanocapsa feldmannii'])
```




    '614894'



#### code to get record from aphia id

Had some [issues with the parallelization](https://stackoverflow.com/questions/50168647/multiprocessing-causes-python-to-crash-and-gives-an-error-may-have-been-in-progr) on Mac M1. Adding 'OBJC_DISABLE_INITIALIZE_FORK_SAFETY = YES' to .bash_profile and then [This](https://github.com/python/cpython/issues/74570) fixed it.   
Try to run without the bash_profile fix LATER.


```python
os.environ["no_proxy"]="*"
```


```python
tax_18S = asv_tables['18S V9'][['taxonomy','domain','supergroup','division','subdivision','class','order','family','genus','species']]
```


```python
tax_18S = tax_18S.drop_duplicates(ignore_index=True)
tax_18S.shape
```




    (1374, 10)




```python
if __name__ == '__main__':
    worms_18s = WoRMS_matching.get_worms_from_aphiaid_or_name_parallel(
    tax_df = tax_18S,worms_dict=pr2_18S_dict,ordered_rank_columns=['species','genus','family','order','class','subdivision','division','supergroup'],
    full_tax_column="taxonomy",full_tax_vI=True,n_proc=6)
    
```


```python
worms_18s.head()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>full_tax</th>
      <th>verbatimIdentification</th>
      <th>old_taxonRank</th>
      <th>old name</th>
      <th>scientificName</th>
      <th>scientificNameID</th>
      <th>kingdom</th>
      <th>phylum</th>
      <th>class</th>
      <th>order</th>
      <th>family</th>
      <th>genus</th>
      <th>taxonRank</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>Eukaryota;Obazoa;Opisthokonta;Fungi;Ascomycota;Pezizomycotina;Eurotiomycetes;Aspergillus;Aspergillus_penicillioides;</td>
      <td>Eukaryota;Obazoa;Opisthokonta;Fungi;Ascomycota;Pezizomycotina;Eurotiomycetes;Aspergillus;Aspergillus_penicillioides;</td>
      <td>genus</td>
      <td>Aspergillus</td>
      <td>Aspergillus</td>
      <td>urn:lsid:marinespecies.org:taxname:100211</td>
      <td>Fungi</td>
      <td>Ascomycota</td>
      <td>Eurotiomycetes</td>
      <td>Eurotiales</td>
      <td>Trichocomaceae</td>
      <td>Aspergillus</td>
      <td>Genus</td>
    </tr>
    <tr>
      <th>1</th>
      <td>Eukaryota;Cryptista;Cryptophyta;Cryptophyta_X;Cryptophyceae;Goniomonadales;Goniomonadales_X;Goniomonas;Goniomonas_sp.;</td>
      <td>Eukaryota;Cryptista;Cryptophyta;Cryptophyta_X;Cryptophyceae;Goniomonadales;Goniomonadales_X;Goniomonas;Goniomonas_sp.;</td>
      <td>species</td>
      <td>Goniomonas</td>
      <td>Goniomonas</td>
      <td>urn:lsid:marinespecies.org:taxname:106286</td>
      <td>Chromista</td>
      <td>Cryptophyta</td>
      <td>Cryptophyceae</td>
      <td>Cryptomonadales</td>
      <td>Cryptomonadaceae</td>
      <td>Goniomonas</td>
      <td>Genus</td>
    </tr>
    <tr>
      <th>2</th>
      <td>Eukaryota;TSAR;Alveolata;Ciliophora;Spirotrichea;Oligotrichida;Strombidiidae;Strombidium;Strombidium_sp.;</td>
      <td>Eukaryota;TSAR;Alveolata;Ciliophora;Spirotrichea;Oligotrichida;Strombidiidae;Strombidium;Strombidium_sp.;</td>
      <td>species</td>
      <td>Strombidium</td>
      <td>Strombidium</td>
      <td>urn:lsid:marinespecies.org:taxname:101195</td>
      <td>Chromista</td>
      <td>Ciliophora</td>
      <td>Oligotrichea</td>
      <td>Oligotrichida</td>
      <td>Strombidiidae</td>
      <td>Strombidium</td>
      <td>Genus</td>
    </tr>
    <tr>
      <th>3</th>
      <td>Eukaryota;Obazoa;Opisthokonta;Metazoa;Annelida;Annelida_X;Annelida_XX;Prionospio;Prionospio_dubia;</td>
      <td>Prionospio dubia</td>
      <td>NaN</td>
      <td>aphiaID</td>
      <td>Prionospio dubia</td>
      <td>urn:lsid:marinespecies.org:taxname:131155</td>
      <td>Animalia</td>
      <td>Annelida</td>
      <td>Polychaeta</td>
      <td>Spionida</td>
      <td>Spionidae</td>
      <td>Prionospio</td>
      <td>Species</td>
    </tr>
    <tr>
      <th>4</th>
      <td>Eukaryota;TSAR;Stramenopiles;Bigyra;Opalozoa;Opalozoa_X;MAST-12;MAST-12A;MAST-12A_sp.;</td>
      <td>Eukaryota;TSAR;Stramenopiles;Bigyra;Opalozoa;Opalozoa_X;MAST-12;MAST-12A;MAST-12A_sp.;</td>
      <td>class</td>
      <td>Opalozoa</td>
      <td>Opalozoa</td>
      <td>urn:lsid:marinespecies.org:taxname:582466</td>
      <td>Chromista</td>
      <td>Bigyra</td>
      <td>None</td>
      <td>None</td>
      <td>None</td>
      <td>None</td>
      <td>Subphylum</td>
    </tr>
  </tbody>
</table>
</div>




```python
# which taxa had absolutely no matches
worms_18s[worms_18s["scientificName"]=="No match"]['old name'].unique()
```




    array(['Haptista', 'Provora', 'Obazoa', 'TSAR', 'Cryptista:nucl',
           'Archaeplastida'], dtype=object)




```python
worms_18s[worms_18s["scientificName"]=="No match"].head(20)
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>full_tax</th>
      <th>verbatimIdentification</th>
      <th>old_taxonRank</th>
      <th>old name</th>
      <th>scientificName</th>
      <th>scientificNameID</th>
      <th>kingdom</th>
      <th>phylum</th>
      <th>class</th>
      <th>order</th>
      <th>family</th>
      <th>genus</th>
      <th>taxonRank</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>121</th>
      <td>Eukaryota;Haptista;Centroplasthelida;Centroplasthelida_X;Pterocystida</td>
      <td>Eukaryota;Haptista;Centroplasthelida;Centroplasthelida_X;Pterocystida</td>
      <td>supergroup</td>
      <td>Haptista</td>
      <td>No match</td>
      <td>None</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>4</th>
      <td>Eukaryota;Provora;Nibbleridia;Nibbleridia_X;Nibbleridea;Nibbleridida;Nibbleridae;Nibbleromonas</td>
      <td>Eukaryota;Provora;Nibbleridia;Nibbleridia_X;Nibbleridea;Nibbleridida;Nibbleridae;Nibbleromonas</td>
      <td>supergroup</td>
      <td>Provora</td>
      <td>No match</td>
      <td>None</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>48</th>
      <td>Eukaryota;Haptista;Centroplasthelida;Centroplasthelida_X;Centroplasthelida_XX;Centroplasthelida_XXX;Centroplasthelida_XXXX;Centroplasthelida_XXXXX...</td>
      <td>Eukaryota;Haptista;Centroplasthelida;Centroplasthelida_X;Centroplasthelida_XX;Centroplasthelida_XXX;Centroplasthelida_XXXX;Centroplasthelida_XXXXX...</td>
      <td>supergroup</td>
      <td>Haptista</td>
      <td>No match</td>
      <td>None</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>62</th>
      <td>Eukaryota;Provora;Nebulidia;Nebulidia_X;Nebulidea;Nebulidida;Nebulidae;Nebulomonas;Nebulomonas_marisrubri;</td>
      <td>Eukaryota;Provora;Nebulidia;Nebulidia_X;Nebulidea;Nebulidida;Nebulidae;Nebulomonas;Nebulomonas_marisrubri;</td>
      <td>supergroup</td>
      <td>Provora</td>
      <td>No match</td>
      <td>None</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>72</th>
      <td>Eukaryota;Obazoa;Breviatea;Breviatea_X;Breviatea_XX;Breviatea_XXX;Breviata-lineage;Breviata-lineage_X;Breviata-lineage_X_sp.;</td>
      <td>Eukaryota;Obazoa;Breviatea;Breviatea_X;Breviatea_XX;Breviatea_XXX;Breviata-lineage;Breviata-lineage_X;Breviata-lineage_X_sp.;</td>
      <td>supergroup</td>
      <td>Obazoa</td>
      <td>No match</td>
      <td>None</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>108</th>
      <td>Eukaryota;TSAR;Stramenopiles;Gyrista;Peronosporomycetes;Peronosporomycetes_X;Haliphthorales</td>
      <td>Eukaryota;TSAR;Stramenopiles;Gyrista;Peronosporomycetes;Peronosporomycetes_X;Haliphthorales</td>
      <td>supergroup</td>
      <td>TSAR</td>
      <td>No match</td>
      <td>None</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>120</th>
      <td>Eukaryota;TSAR;Telonemia;Telonemia_X;Telonemia_XX;Telonemia_XXX;Telonemia-Group-2</td>
      <td>Eukaryota;TSAR;Telonemia;Telonemia_X;Telonemia_XX;Telonemia_XXX;Telonemia-Group-2</td>
      <td>supergroup</td>
      <td>TSAR</td>
      <td>No match</td>
      <td>None</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>188</th>
      <td>Eukaryota;TSAR;Stramenopiles;Gyrista;Peronosporomycetes;Peronosporomycetes_X;Haliphthorales;Halocrusticida;Halocrusticida_baliensis;</td>
      <td>Eukaryota;TSAR;Stramenopiles;Gyrista;Peronosporomycetes;Peronosporomycetes_X;Haliphthorales;Halocrusticida;Halocrusticida_baliensis;</td>
      <td>supergroup</td>
      <td>TSAR</td>
      <td>No match</td>
      <td>None</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>191</th>
      <td>Eukaryota;Haptista;Centroplasthelida;Centroplasthelida_X;Pterocystida;Pterista</td>
      <td>Eukaryota;Haptista;Centroplasthelida;Centroplasthelida_X;Pterocystida;Pterista</td>
      <td>supergroup</td>
      <td>Haptista</td>
      <td>No match</td>
      <td>None</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>195</th>
      <td>Eukaryota:nucl;Cryptista:nucl;Cryptophyta:nucl;Cryptophyta_X:nucl;Cryptophyceae:nucl;Cryptomonadales:nucl;Hemiselmidaceae:nucl;Chroomonas:nucl;Chr...</td>
      <td>Eukaryota:nucl;Cryptista:nucl;Cryptophyta:nucl;Cryptophyta_X:nucl;Cryptophyceae:nucl;Cryptomonadales:nucl;Hemiselmidaceae:nucl;Chroomonas:nucl;Chr...</td>
      <td>supergroup</td>
      <td>Cryptista:nucl</td>
      <td>No match</td>
      <td>None</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>227</th>
      <td>Eukaryota;Obazoa</td>
      <td>Eukaryota;Obazoa</td>
      <td>supergroup</td>
      <td>Obazoa</td>
      <td>No match</td>
      <td>None</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>4</th>
      <td>Eukaryota;Archaeplastida</td>
      <td>Eukaryota;Archaeplastida</td>
      <td>supergroup</td>
      <td>Archaeplastida</td>
      <td>No match</td>
      <td>None</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>19</th>
      <td>Eukaryota;TSAR;Stramenopiles;Gyrista;Gyrista_X;Gyrista_XX;MAST-1;MAST-1B;MAST-1B_sp.;</td>
      <td>Eukaryota;TSAR;Stramenopiles;Gyrista;Gyrista_X;Gyrista_XX;MAST-1;MAST-1B;MAST-1B_sp.;</td>
      <td>supergroup</td>
      <td>TSAR</td>
      <td>No match</td>
      <td>None</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>46</th>
      <td>Eukaryota:nucl;Cryptista:nucl;Cryptophyta:nucl;Cryptophyta_X:nucl;Cryptophyceae:nucl;Cryptomonadales:nucl</td>
      <td>Eukaryota:nucl;Cryptista:nucl;Cryptophyta:nucl;Cryptophyta_X:nucl;Cryptophyceae:nucl;Cryptomonadales:nucl</td>
      <td>supergroup</td>
      <td>Cryptista:nucl</td>
      <td>No match</td>
      <td>None</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>60</th>
      <td>Eukaryota;TSAR</td>
      <td>Eukaryota;TSAR</td>
      <td>supergroup</td>
      <td>TSAR</td>
      <td>No match</td>
      <td>None</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>65</th>
      <td>Eukaryota;Obazoa;Opisthokonta;Metazoa;Ctenophora;Ctenophora_X;Ctenophora_XX</td>
      <td>Eukaryota;Obazoa;Opisthokonta;Metazoa;Ctenophora;Ctenophora_X;Ctenophora_XX</td>
      <td>supergroup</td>
      <td>Obazoa</td>
      <td>No match</td>
      <td>None</td>
      <td>Animalia</td>
      <td>None</td>
      <td>None</td>
      <td>None</td>
      <td>None</td>
      <td>None</td>
      <td>Kingdom</td>
    </tr>
    <tr>
      <th>81</th>
      <td>Eukaryota;TSAR;Stramenopiles;Gyrista;Gyrista_X;Gyrista_XX;MAST-1;MAST-1A;MAST-1A_sp.;</td>
      <td>Eukaryota;TSAR;Stramenopiles;Gyrista;Gyrista_X;Gyrista_XX;MAST-1;MAST-1A;MAST-1A_sp.;</td>
      <td>supergroup</td>
      <td>TSAR</td>
      <td>No match</td>
      <td>None</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>119</th>
      <td>Eukaryota;Obazoa;Opisthokonta;Opisthokonta_X;Opisthokonta_XX;Opisthokonta_XXX;Opisthokonta_XXXX;Opisthokonta_XXXXX;Opisthokonta_XXXXX_sp.;</td>
      <td>Eukaryota;Obazoa;Opisthokonta;Opisthokonta_X;Opisthokonta_XX;Opisthokonta_XXX;Opisthokonta_XXXX;Opisthokonta_XXXXX;Opisthokonta_XXXXX_sp.;</td>
      <td>supergroup</td>
      <td>Obazoa</td>
      <td>No match</td>
      <td>None</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>177</th>
      <td>Eukaryota;TSAR;Stramenopiles;Gyrista;Mediophyceae</td>
      <td>Eukaryota;TSAR;Stramenopiles;Gyrista;Mediophyceae</td>
      <td>supergroup</td>
      <td>TSAR</td>
      <td>No match</td>
      <td>None</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>203</th>
      <td>Eukaryota;TSAR;Telonemia;Telonemia_X;Telonemia_XX;Telonemia_XXX;Telonemia-Group-1</td>
      <td>Eukaryota;TSAR;Telonemia;Telonemia_X;Telonemia_XX;Telonemia_XXX;Telonemia-Group-1</td>
      <td>supergroup</td>
      <td>TSAR</td>
      <td>No match</td>
      <td>None</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
  </tbody>
</table>
</div>




```python
worms_18s.loc[worms_18s["scientificName"]=="No match",'scientificName'] = "Biota"
worms_18s.loc[worms_18s["scientificName"]=="Biota",'scientificNameID'] = "urn:lsid:marinespecies.org:taxname:1"

```


```python
worms_18s[worms_18s['scientificName'].isna() == True]
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>full_tax</th>
      <th>verbatimIdentification</th>
      <th>old_taxonRank</th>
      <th>old name</th>
      <th>scientificName</th>
      <th>scientificNameID</th>
      <th>kingdom</th>
      <th>phylum</th>
      <th>class</th>
      <th>order</th>
      <th>family</th>
      <th>genus</th>
      <th>taxonRank</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>109</th>
      <td>Unassigned</td>
      <td>Unassigned</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>120</th>
      <td>Eukaryota;Haptista</td>
      <td>Eukaryota;Haptista</td>
      <td>supergroup</td>
      <td>Haptista</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>77</th>
      <td>Eukaryota</td>
      <td>Eukaryota</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
  </tbody>
</table>
</div>




```python
worms_18s.loc[worms_18s["full_tax"]=="Eukaryota;Haptista",'scientificName'] = "Biota"
worms_18s.loc[worms_18s["full_tax"]=="Eukaryota;Haptista",'scientificNameID'] = "urn:lsid:marinespecies.org:taxname:1"
worms_18s.loc[worms_18s["full_tax"]=="Eukaryota",'scientificName'] = "Biota"
worms_18s.loc[worms_18s["full_tax"]=="Eukaryota",'scientificNameID'] = "urn:lsid:marinespecies.org:taxname:1"

```


```python
worms_18s[worms_18s['scientificName'].isna() == True]
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>full_tax</th>
      <th>verbatimIdentification</th>
      <th>old_taxonRank</th>
      <th>old name</th>
      <th>scientificName</th>
      <th>scientificNameID</th>
      <th>kingdom</th>
      <th>phylum</th>
      <th>class</th>
      <th>order</th>
      <th>family</th>
      <th>genus</th>
      <th>taxonRank</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>109</th>
      <td>Unassigned</td>
      <td>Unassigned</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
  </tbody>
</table>
</div>




```python

print(worms_18s[worms_18s['scientificName'].isna() == True].shape)
worms_18s.loc[worms_18s['scientificName'].isna() == True,'scientificName'] = 'incertae sedis'
worms_18s.loc[worms_18s['scientificName'] == 'incertae sedis','scientificNameID'] =  'urn:lsid:marinespecies.org:taxname:12'
print(worms_18s[worms_18s['scientificName'].isna() == True].shape)
```

    (1, 13)
    (0, 13)



```python
worms_18s[worms_18s["old name"]=="aphiaID"].shape
```




    (332, 13)




```python
worms_18s.to_csv("../processed/worms_18S_matching.tsv",sep="\t",index=False)
```


```python
worms_18s.drop(columns=['old name','old_taxonRank'],inplace=True)
worms_18s.head()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>full_tax</th>
      <th>verbatimIdentification</th>
      <th>scientificName</th>
      <th>scientificNameID</th>
      <th>kingdom</th>
      <th>phylum</th>
      <th>class</th>
      <th>order</th>
      <th>family</th>
      <th>genus</th>
      <th>taxonRank</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>Eukaryota;Obazoa;Opisthokonta;Fungi;Ascomycota;Pezizomycotina;Eurotiomycetes;Aspergillus;Aspergillus_penicillioides;</td>
      <td>Eukaryota;Obazoa;Opisthokonta;Fungi;Ascomycota;Pezizomycotina;Eurotiomycetes;Aspergillus;Aspergillus_penicillioides;</td>
      <td>Aspergillus</td>
      <td>urn:lsid:marinespecies.org:taxname:100211</td>
      <td>Fungi</td>
      <td>Ascomycota</td>
      <td>Eurotiomycetes</td>
      <td>Eurotiales</td>
      <td>Trichocomaceae</td>
      <td>Aspergillus</td>
      <td>Genus</td>
    </tr>
    <tr>
      <th>1</th>
      <td>Eukaryota;Cryptista;Cryptophyta;Cryptophyta_X;Cryptophyceae;Goniomonadales;Goniomonadales_X;Goniomonas;Goniomonas_sp.;</td>
      <td>Eukaryota;Cryptista;Cryptophyta;Cryptophyta_X;Cryptophyceae;Goniomonadales;Goniomonadales_X;Goniomonas;Goniomonas_sp.;</td>
      <td>Goniomonas</td>
      <td>urn:lsid:marinespecies.org:taxname:106286</td>
      <td>Chromista</td>
      <td>Cryptophyta</td>
      <td>Cryptophyceae</td>
      <td>Cryptomonadales</td>
      <td>Cryptomonadaceae</td>
      <td>Goniomonas</td>
      <td>Genus</td>
    </tr>
    <tr>
      <th>2</th>
      <td>Eukaryota;TSAR;Alveolata;Ciliophora;Spirotrichea;Oligotrichida;Strombidiidae;Strombidium;Strombidium_sp.;</td>
      <td>Eukaryota;TSAR;Alveolata;Ciliophora;Spirotrichea;Oligotrichida;Strombidiidae;Strombidium;Strombidium_sp.;</td>
      <td>Strombidium</td>
      <td>urn:lsid:marinespecies.org:taxname:101195</td>
      <td>Chromista</td>
      <td>Ciliophora</td>
      <td>Oligotrichea</td>
      <td>Oligotrichida</td>
      <td>Strombidiidae</td>
      <td>Strombidium</td>
      <td>Genus</td>
    </tr>
    <tr>
      <th>3</th>
      <td>Eukaryota;Obazoa;Opisthokonta;Metazoa;Annelida;Annelida_X;Annelida_XX;Prionospio;Prionospio_dubia;</td>
      <td>Prionospio dubia</td>
      <td>Prionospio dubia</td>
      <td>urn:lsid:marinespecies.org:taxname:131155</td>
      <td>Animalia</td>
      <td>Annelida</td>
      <td>Polychaeta</td>
      <td>Spionida</td>
      <td>Spionidae</td>
      <td>Prionospio</td>
      <td>Species</td>
    </tr>
    <tr>
      <th>4</th>
      <td>Eukaryota;TSAR;Stramenopiles;Bigyra;Opalozoa;Opalozoa_X;MAST-12;MAST-12A;MAST-12A_sp.;</td>
      <td>Eukaryota;TSAR;Stramenopiles;Bigyra;Opalozoa;Opalozoa_X;MAST-12;MAST-12A;MAST-12A_sp.;</td>
      <td>Opalozoa</td>
      <td>urn:lsid:marinespecies.org:taxname:582466</td>
      <td>Chromista</td>
      <td>Bigyra</td>
      <td>None</td>
      <td>None</td>
      <td>None</td>
      <td>None</td>
      <td>Subphylum</td>
    </tr>
  </tbody>
</table>
</div>



#### Merge Occurrence and worms


```python
occ['18S V9'].shape
```




    (146232, 16)




```python
# Get identificationRemarks
occ18_test = occ['18S V9'].copy()
occ18_test.drop(columns=['domain','supergroup','division','subdivision','class','order','family','genus','species'],inplace=True)
#occ18_test.drop(columns=['old name'],inplace=True)

occ18_test = occ18_test.merge(worms_18s, how='left', left_on ='taxonomy', right_on='full_tax')
occ18_test.drop(columns='full_tax', inplace=True)
occ18_test.head()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>featureid</th>
      <th>sequence</th>
      <th>taxonomy</th>
      <th>Confidence</th>
      <th>eventID</th>
      <th>organismQuantity</th>
      <th>occurrenceID</th>
      <th>verbatimIdentification</th>
      <th>scientificName</th>
      <th>scientificNameID</th>
      <th>kingdom</th>
      <th>phylum</th>
      <th>class</th>
      <th>order</th>
      <th>family</th>
      <th>genus</th>
      <th>taxonRank</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>36aa75f9b28f5f831c2d631ba65c2bcb</td>
      <td>GCTACTACCGATTGAACGTTTTAGTGAGGTCCTCGGACTGTTTGCCTGGCGGATTACTCTGCCTGGCTGGCGGGAAGACGACCAAACTGTAGCGTTTAGAGGAAGTAAAAGTCGTAACAAGGTTTCC</td>
      <td>Eukaryota;Obazoa;Opisthokonta;Metazoa;Arthropoda;Crustacea;Maxillopoda;Neocalanus;Neocalanus_cristatus;</td>
      <td>0.922099</td>
      <td>GOMECC4_27N_Sta1_DCM_A</td>
      <td>1516</td>
      <td>GOMECC4_27N_Sta1_DCM_A_occ36aa75f9b28f5f831c2d631ba65c2bcb</td>
      <td>Neocalanus cristatus</td>
      <td>Neocalanus cristatus</td>
      <td>urn:lsid:marinespecies.org:taxname:104470</td>
      <td>Animalia</td>
      <td>Arthropoda</td>
      <td>Copepoda</td>
      <td>Calanoida</td>
      <td>Calanidae</td>
      <td>Neocalanus</td>
      <td>Species</td>
    </tr>
    <tr>
      <th>1</th>
      <td>4e38e8ced9070952b314e1880bede1ca</td>
      <td>GCTACTACCGATTGAACGTTTTAGTGAGGTCCTCGGACTGTTTGGTAGTCGGATCACTCTGACTGCCTGGCGGGAAGACGACCAAACTGTAGCGTTTAGAGGAAGTAAAAGTCGTAACAAGGTTTCC</td>
      <td>Eukaryota;Obazoa;Opisthokonta;Metazoa;Arthropoda;Crustacea;Maxillopoda;Clausocalanus;Clausocalanus_furcatus;</td>
      <td>0.999947</td>
      <td>GOMECC4_27N_Sta1_DCM_A</td>
      <td>962</td>
      <td>GOMECC4_27N_Sta1_DCM_A_occ4e38e8ced9070952b314e1880bede1ca</td>
      <td>Clausocalanus furcatus</td>
      <td>Clausocalanus furcatus</td>
      <td>urn:lsid:marinespecies.org:taxname:104503</td>
      <td>Animalia</td>
      <td>Arthropoda</td>
      <td>Copepoda</td>
      <td>Calanoida</td>
      <td>Clausocalanidae</td>
      <td>Clausocalanus</td>
      <td>Species</td>
    </tr>
    <tr>
      <th>2</th>
      <td>2a31e5c01634165da99e7381279baa75</td>
      <td>GCTACTACCGATTGGACGTTTTAGTGAGACATTTGGACTGGGTTAAGATAGTCGCAAGACTACCTTTTCTCCGGAAAGACTTTCAAACTTGAGCGTCTAGAGGAAGTAAAAGTCGTAACAAGGTTTCC</td>
      <td>Eukaryota;Obazoa;Opisthokonta;Metazoa;Arthropoda;Crustacea;Maxillopoda;Acrocalanus;Acrocalanus_sp.;</td>
      <td>0.779948</td>
      <td>GOMECC4_27N_Sta1_DCM_A</td>
      <td>1164</td>
      <td>GOMECC4_27N_Sta1_DCM_A_occ2a31e5c01634165da99e7381279baa75</td>
      <td>Eukaryota;Obazoa;Opisthokonta;Metazoa;Arthropoda;Crustacea;Maxillopoda;Acrocalanus;Acrocalanus_sp.;</td>
      <td>Acrocalanus</td>
      <td>urn:lsid:marinespecies.org:taxname:104192</td>
      <td>Animalia</td>
      <td>Arthropoda</td>
      <td>Copepoda</td>
      <td>Calanoida</td>
      <td>Paracalanidae</td>
      <td>Acrocalanus</td>
      <td>Genus</td>
    </tr>
    <tr>
      <th>3</th>
      <td>ecee60339b2fb88ea6d1c8d18054bed4</td>
      <td>GCTCCTACCGATTGAGTGATCCGGTGAATAATTCGGACTGCAGCAGTGTTCAGTTCCTGAACGTTGCAGCGGAAAGTTTAGTGAACCTTATCACTTAGAGGAAGGAGAAGTCGTAACAAGGTTTCC</td>
      <td>Eukaryota;TSAR;Alveolata;Dinoflagellata;Dinophyceae</td>
      <td>0.999931</td>
      <td>GOMECC4_27N_Sta1_DCM_A</td>
      <td>287</td>
      <td>GOMECC4_27N_Sta1_DCM_A_occecee60339b2fb88ea6d1c8d18054bed4</td>
      <td>Eukaryota;TSAR;Alveolata;Dinoflagellata;Dinophyceae</td>
      <td>Dinophyceae</td>
      <td>urn:lsid:marinespecies.org:taxname:19542</td>
      <td>Chromista</td>
      <td>Myzozoa</td>
      <td>Dinophyceae</td>
      <td>None</td>
      <td>None</td>
      <td>None</td>
      <td>Class</td>
    </tr>
    <tr>
      <th>4</th>
      <td>fa1f1a97dd4ae7c826009186bad26384</td>
      <td>GCTCCTACCGATTGAGTGATCCGGTGAATAATTCGGACTGCAGCAATGTTTGGATCCCGAACGTTGCAGCGGAAAGTTTAGTGAACCTTATCACTTAGAGGAAGGAGAAGTCGTAACAAGGTTTCC</td>
      <td>Eukaryota;TSAR;Alveolata;Dinoflagellata;Dinophyceae;Gymnodiniales;Gymnodiniaceae</td>
      <td>0.986908</td>
      <td>GOMECC4_27N_Sta1_DCM_A</td>
      <td>250</td>
      <td>GOMECC4_27N_Sta1_DCM_A_occfa1f1a97dd4ae7c826009186bad26384</td>
      <td>Eukaryota;TSAR;Alveolata;Dinoflagellata;Dinophyceae;Gymnodiniales;Gymnodiniaceae</td>
      <td>Gymnodiniaceae</td>
      <td>urn:lsid:marinespecies.org:taxname:109410</td>
      <td>Chromista</td>
      <td>Myzozoa</td>
      <td>Dinophyceae</td>
      <td>Gymnodiniales</td>
      <td>Gymnodiniaceae</td>
      <td>None</td>
      <td>Family</td>
    </tr>
  </tbody>
</table>
</div>



#### identificationRemarks


```python
data['analysis_data'].head()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>amplicon_sequenced</th>
      <th>ampliconSize</th>
      <th>trim_method</th>
      <th>cluster_method</th>
      <th>pid_clustering</th>
      <th>taxa_class_method</th>
      <th>taxa_ref_db</th>
      <th>code_repo</th>
      <th>identificationReferences</th>
      <th>controls_used</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>16S V4-V5</td>
      <td>411</td>
      <td>cutadapt</td>
      <td>Tourmaline; qiime2-2021.2; dada2</td>
      <td>ASV</td>
      <td>Tourmaline; qiime2-2021.2; naive-bayes classifier</td>
      <td>Silva SSU Ref NR 99 v138.1; 515f-926r region; 10.5281/zenodo.8392695</td>
      <td>https://github.com/aomlomics/gomecc</td>
      <td>10.5281/zenodo.8392695 | https://github.com/aomlomics/tourmaline</td>
      <td>12 distilled water blanks | 2 PCR no-template controls | 7 extraction blanks | 12 2nd PCR no-template controls | 3 Zymo mock community</td>
    </tr>
    <tr>
      <th>1</th>
      <td>18S V9</td>
      <td>260</td>
      <td>cutadapt</td>
      <td>Tourmaline; qiime2-2021.2; dada2</td>
      <td>ASV</td>
      <td>Tourmaline; qiime2-2021.2; naive-bayes classifier</td>
      <td>PR2 v5.0.1; V9 1391f-1510r region; 10.5281/zenodo.8392706</td>
      <td>https://github.com/aomlomics/gomecc</td>
      <td>10.5281/zenodo.8392706 | https://pr2-database.org/ | https://github.com/aomlomics/tourmaline</td>
      <td>12 distilled water blanks | 2 PCR no-template controls | 7 extraction blanks | 7 2nd PCR no-template controls</td>
    </tr>
  </tbody>
</table>
</div>




```python
occ18_test['taxa_class_method'] = data['analysis_data'].loc[data['analysis_data']['amplicon_sequenced'] == '18S V9','taxa_class_method'].item()
occ18_test['taxa_ref_db'] = data['analysis_data'].loc[data['analysis_data']['amplicon_sequenced'] == '18S V9','taxa_ref_db'].item()

occ18_test['identificationRemarks'] = occ18_test['taxa_class_method'] +", confidence (at lowest specified taxon): "+occ18_test['Confidence'].astype(str) +", against reference database: "+occ18_test['taxa_ref_db']
```


```python
occ18_test['identificationRemarks'][0]
```




    'Tourmaline; qiime2-2021.2; naive-bayes classifier, confidence (at lowest specified taxon): 0.92209885, against reference database: PR2 v5.0.1; V9 1391f-1510r region; 10.5281/zenodo.8392706'



#### taxonID, basisOfRecord, eventID, nameAccordingTo, organismQuantityType


```python
occ18_test['taxonID'] = 'ASV:'+occ18_test['featureid']
occ18_test['basisOfRecord'] = 'MaterialSample'
occ18_test['nameAccordingTo'] = "WoRMS"
occ18_test['organismQuantityType'] = "DNA sequence reads"
occ18_test['recordedBy'] = data['study_data']['recordedBy'].values[0]
```

#### associatedSequences, materialSampleID


```python
data['prep_data'].columns
```




    Index(['sample_name', 'library_id', 'title', 'library_strategy',
           'library_source', 'library_selection', 'lib_layout', 'platform',
           'instrument_model', 'design_description', 'filetype', 'filename',
           'filename2', 'biosample_accession', 'sra_accession', 'seq_method',
           'nucl_acid_ext', 'amplicon_sequenced', 'target_gene',
           'target_subfragment', 'pcr_primer_forward', 'pcr_primer_reverse',
           'pcr_primer_name_forward', 'pcr_primer_name_reverse',
           'pcr_primer_reference', 'pcr_cond', 'nucl_acid_amp', 'adapters',
           'mid_barcode'],
          dtype='object')




```python
occ18_test = occ18_test.merge(data['prep_data'].loc[data['prep_data']['amplicon_sequenced'] == '18S V9',['sample_name','sra_accession','biosample_accession']], how='left', left_on ='eventID', right_on='sample_name')
```

#### eventID


```python
occ18_test['eventID'] = occ18_test['eventID']+"_18S"
```

#### sampleSize


```python
# get sampleSize by total number of reads per sample
x = asv_tables['18S V9'].sum(numeric_only=True).astype('int')
x.index = x.index+"_18S"
occ18_test['sampleSizeValue'] = occ18_test['eventID'].map(x).astype('str')
occ18_test['sampleSizeUnit'] = 'DNA sequence reads'
```


```python
# drop unnneeded columns
occ18_test.drop(columns=['sample_name','featureid','taxonomy','Confidence','taxa_class_method','taxa_ref_db'],inplace=True)
```


```python
occ18_test['associatedSequences'] = occ18_test['sra_accession']+' | '+ occ18_test['biosample_accession']+' | '+data['study_data']['bioproject_accession'].values[0]
```


```python
occ18_test.rename(columns={'biosample_accession': 'materialSampleID',
                  'sequence': 'DNA_sequence'},inplace=True)
                   
```


```python
# drop unnneeded columns
occ18_test.drop(columns=['sra_accession'],inplace=True)
```


```python
occ18_test.columns
```




    Index(['DNA_sequence', 'eventID', 'organismQuantity', 'occurrenceID',
           'verbatimIdentification', 'scientificName', 'scientificNameID',
           'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'taxonRank',
           'identificationRemarks', 'taxonID', 'basisOfRecord', 'nameAccordingTo',
           'organismQuantityType', 'recordedBy', 'materialSampleID',
           'sampleSizeValue', 'sampleSizeUnit', 'associatedSequences'],
          dtype='object')




```python
occ18_test.head()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>DNA_sequence</th>
      <th>eventID</th>
      <th>organismQuantity</th>
      <th>occurrenceID</th>
      <th>verbatimIdentification</th>
      <th>scientificName</th>
      <th>scientificNameID</th>
      <th>kingdom</th>
      <th>phylum</th>
      <th>class</th>
      <th>order</th>
      <th>family</th>
      <th>genus</th>
      <th>taxonRank</th>
      <th>identificationRemarks</th>
      <th>taxonID</th>
      <th>basisOfRecord</th>
      <th>nameAccordingTo</th>
      <th>organismQuantityType</th>
      <th>recordedBy</th>
      <th>materialSampleID</th>
      <th>sampleSizeValue</th>
      <th>sampleSizeUnit</th>
      <th>associatedSequences</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>GCTACTACCGATTGAACGTTTTAGTGAGGTCCTCGGACTGTTTGCCTGGCGGATTACTCTGCCTGGCTGGCGGGAAGACGACCAAACTGTAGCGTTTAGAGGAAGTAAAAGTCGTAACAAGGTTTCC</td>
      <td>GOMECC4_27N_Sta1_DCM_A_18S</td>
      <td>1516</td>
      <td>GOMECC4_27N_Sta1_DCM_A_occ36aa75f9b28f5f831c2d631ba65c2bcb</td>
      <td>Neocalanus cristatus</td>
      <td>Neocalanus cristatus</td>
      <td>urn:lsid:marinespecies.org:taxname:104470</td>
      <td>Animalia</td>
      <td>Arthropoda</td>
      <td>Copepoda</td>
      <td>Calanoida</td>
      <td>Calanidae</td>
      <td>Neocalanus</td>
      <td>Species</td>
      <td>Tourmaline; qiime2-2021.2; naive-bayes classifier, confidence (at lowest specified taxon): 0.92209885, against reference database: PR2 v5.0.1; V9 ...</td>
      <td>ASV:36aa75f9b28f5f831c2d631ba65c2bcb</td>
      <td>MaterialSample</td>
      <td>WoRMS</td>
      <td>DNA sequence reads</td>
      <td>Luke Thompson | Katherine Silliman</td>
      <td>SAMN37516094</td>
      <td>9838</td>
      <td>DNA sequence reads</td>
      <td>SRR26161153 | SAMN37516094 | PRJNA887898</td>
    </tr>
    <tr>
      <th>1</th>
      <td>GCTACTACCGATTGAACGTTTTAGTGAGGTCCTCGGACTGTTTGGTAGTCGGATCACTCTGACTGCCTGGCGGGAAGACGACCAAACTGTAGCGTTTAGAGGAAGTAAAAGTCGTAACAAGGTTTCC</td>
      <td>GOMECC4_27N_Sta1_DCM_A_18S</td>
      <td>962</td>
      <td>GOMECC4_27N_Sta1_DCM_A_occ4e38e8ced9070952b314e1880bede1ca</td>
      <td>Clausocalanus furcatus</td>
      <td>Clausocalanus furcatus</td>
      <td>urn:lsid:marinespecies.org:taxname:104503</td>
      <td>Animalia</td>
      <td>Arthropoda</td>
      <td>Copepoda</td>
      <td>Calanoida</td>
      <td>Clausocalanidae</td>
      <td>Clausocalanus</td>
      <td>Species</td>
      <td>Tourmaline; qiime2-2021.2; naive-bayes classifier, confidence (at lowest specified taxon): 0.999946735, against reference database: PR2 v5.0.1; V9...</td>
      <td>ASV:4e38e8ced9070952b314e1880bede1ca</td>
      <td>MaterialSample</td>
      <td>WoRMS</td>
      <td>DNA sequence reads</td>
      <td>Luke Thompson | Katherine Silliman</td>
      <td>SAMN37516094</td>
      <td>9838</td>
      <td>DNA sequence reads</td>
      <td>SRR26161153 | SAMN37516094 | PRJNA887898</td>
    </tr>
    <tr>
      <th>2</th>
      <td>GCTACTACCGATTGGACGTTTTAGTGAGACATTTGGACTGGGTTAAGATAGTCGCAAGACTACCTTTTCTCCGGAAAGACTTTCAAACTTGAGCGTCTAGAGGAAGTAAAAGTCGTAACAAGGTTTCC</td>
      <td>GOMECC4_27N_Sta1_DCM_A_18S</td>
      <td>1164</td>
      <td>GOMECC4_27N_Sta1_DCM_A_occ2a31e5c01634165da99e7381279baa75</td>
      <td>Eukaryota;Obazoa;Opisthokonta;Metazoa;Arthropoda;Crustacea;Maxillopoda;Acrocalanus;Acrocalanus_sp.;</td>
      <td>Acrocalanus</td>
      <td>urn:lsid:marinespecies.org:taxname:104192</td>
      <td>Animalia</td>
      <td>Arthropoda</td>
      <td>Copepoda</td>
      <td>Calanoida</td>
      <td>Paracalanidae</td>
      <td>Acrocalanus</td>
      <td>Genus</td>
      <td>Tourmaline; qiime2-2021.2; naive-bayes classifier, confidence (at lowest specified taxon): 0.779948049, against reference database: PR2 v5.0.1; V9...</td>
      <td>ASV:2a31e5c01634165da99e7381279baa75</td>
      <td>MaterialSample</td>
      <td>WoRMS</td>
      <td>DNA sequence reads</td>
      <td>Luke Thompson | Katherine Silliman</td>
      <td>SAMN37516094</td>
      <td>9838</td>
      <td>DNA sequence reads</td>
      <td>SRR26161153 | SAMN37516094 | PRJNA887898</td>
    </tr>
    <tr>
      <th>3</th>
      <td>GCTCCTACCGATTGAGTGATCCGGTGAATAATTCGGACTGCAGCAGTGTTCAGTTCCTGAACGTTGCAGCGGAAAGTTTAGTGAACCTTATCACTTAGAGGAAGGAGAAGTCGTAACAAGGTTTCC</td>
      <td>GOMECC4_27N_Sta1_DCM_A_18S</td>
      <td>287</td>
      <td>GOMECC4_27N_Sta1_DCM_A_occecee60339b2fb88ea6d1c8d18054bed4</td>
      <td>Eukaryota;TSAR;Alveolata;Dinoflagellata;Dinophyceae</td>
      <td>Dinophyceae</td>
      <td>urn:lsid:marinespecies.org:taxname:19542</td>
      <td>Chromista</td>
      <td>Myzozoa</td>
      <td>Dinophyceae</td>
      <td>None</td>
      <td>None</td>
      <td>None</td>
      <td>Class</td>
      <td>Tourmaline; qiime2-2021.2; naive-bayes classifier, confidence (at lowest specified taxon): 0.999930607, against reference database: PR2 v5.0.1; V9...</td>
      <td>ASV:ecee60339b2fb88ea6d1c8d18054bed4</td>
      <td>MaterialSample</td>
      <td>WoRMS</td>
      <td>DNA sequence reads</td>
      <td>Luke Thompson | Katherine Silliman</td>
      <td>SAMN37516094</td>
      <td>9838</td>
      <td>DNA sequence reads</td>
      <td>SRR26161153 | SAMN37516094 | PRJNA887898</td>
    </tr>
    <tr>
      <th>4</th>
      <td>GCTCCTACCGATTGAGTGATCCGGTGAATAATTCGGACTGCAGCAATGTTTGGATCCCGAACGTTGCAGCGGAAAGTTTAGTGAACCTTATCACTTAGAGGAAGGAGAAGTCGTAACAAGGTTTCC</td>
      <td>GOMECC4_27N_Sta1_DCM_A_18S</td>
      <td>250</td>
      <td>GOMECC4_27N_Sta1_DCM_A_occfa1f1a97dd4ae7c826009186bad26384</td>
      <td>Eukaryota;TSAR;Alveolata;Dinoflagellata;Dinophyceae;Gymnodiniales;Gymnodiniaceae</td>
      <td>Gymnodiniaceae</td>
      <td>urn:lsid:marinespecies.org:taxname:109410</td>
      <td>Chromista</td>
      <td>Myzozoa</td>
      <td>Dinophyceae</td>
      <td>Gymnodiniales</td>
      <td>Gymnodiniaceae</td>
      <td>None</td>
      <td>Family</td>
      <td>Tourmaline; qiime2-2021.2; naive-bayes classifier, confidence (at lowest specified taxon): 0.98690791, against reference database: PR2 v5.0.1; V9 ...</td>
      <td>ASV:fa1f1a97dd4ae7c826009186bad26384</td>
      <td>MaterialSample</td>
      <td>WoRMS</td>
      <td>DNA sequence reads</td>
      <td>Luke Thompson | Katherine Silliman</td>
      <td>SAMN37516094</td>
      <td>9838</td>
      <td>DNA sequence reads</td>
      <td>SRR26161153 | SAMN37516094 | PRJNA887898</td>
    </tr>
  </tbody>
</table>
</div>



### merge event and occurrence


```python
occ18_merged = occ18_test.merge(all_event_data,how='left',on='eventID')
```


```python
occ18_merged.head()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>DNA_sequence</th>
      <th>eventID</th>
      <th>organismQuantity</th>
      <th>occurrenceID</th>
      <th>verbatimIdentification</th>
      <th>scientificName</th>
      <th>scientificNameID</th>
      <th>kingdom</th>
      <th>phylum</th>
      <th>class</th>
      <th>order</th>
      <th>family</th>
      <th>genus</th>
      <th>taxonRank</th>
      <th>identificationRemarks</th>
      <th>taxonID</th>
      <th>basisOfRecord</th>
      <th>nameAccordingTo</th>
      <th>organismQuantityType</th>
      <th>recordedBy</th>
      <th>materialSampleID</th>
      <th>sampleSizeValue</th>
      <th>sampleSizeUnit</th>
      <th>associatedSequences</th>
      <th>locationID</th>
      <th>eventDate</th>
      <th>minimumDepthInMeters</th>
      <th>locality</th>
      <th>decimalLatitude</th>
      <th>decimalLongitude</th>
      <th>samplingProtocol</th>
      <th>waterBody</th>
      <th>maximumDepthInMeters</th>
      <th>parentEventID</th>
      <th>datasetID</th>
      <th>geodeticDatum</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>GCTACTACCGATTGAACGTTTTAGTGAGGTCCTCGGACTGTTTGCCTGGCGGATTACTCTGCCTGGCTGGCGGGAAGACGACCAAACTGTAGCGTTTAGAGGAAGTAAAAGTCGTAACAAGGTTTCC</td>
      <td>GOMECC4_27N_Sta1_DCM_A_18S</td>
      <td>1516</td>
      <td>GOMECC4_27N_Sta1_DCM_A_occ36aa75f9b28f5f831c2d631ba65c2bcb</td>
      <td>Neocalanus cristatus</td>
      <td>Neocalanus cristatus</td>
      <td>urn:lsid:marinespecies.org:taxname:104470</td>
      <td>Animalia</td>
      <td>Arthropoda</td>
      <td>Copepoda</td>
      <td>Calanoida</td>
      <td>Calanidae</td>
      <td>Neocalanus</td>
      <td>Species</td>
      <td>Tourmaline; qiime2-2021.2; naive-bayes classifier, confidence (at lowest specified taxon): 0.92209885, against reference database: PR2 v5.0.1; V9 ...</td>
      <td>ASV:36aa75f9b28f5f831c2d631ba65c2bcb</td>
      <td>MaterialSample</td>
      <td>WoRMS</td>
      <td>DNA sequence reads</td>
      <td>Luke Thompson | Katherine Silliman</td>
      <td>SAMN37516094</td>
      <td>9838</td>
      <td>DNA sequence reads</td>
      <td>SRR26161153 | SAMN37516094 | PRJNA887898</td>
      <td>27N_Sta1</td>
      <td>2021-09-14T11:00-04:00</td>
      <td>49</td>
      <td>USA: Atlantic Ocean, east of Florida (27 N)</td>
      <td>26.997</td>
      <td>-79.618</td>
      <td>CTD rosette</td>
      <td>Atlantic Ocean</td>
      <td>49</td>
      <td>GOMECC4_27N_Sta1_DCM_A</td>
      <td>noaa-aoml-gomecc4</td>
      <td>WGS84</td>
    </tr>
    <tr>
      <th>1</th>
      <td>GCTACTACCGATTGAACGTTTTAGTGAGGTCCTCGGACTGTTTGGTAGTCGGATCACTCTGACTGCCTGGCGGGAAGACGACCAAACTGTAGCGTTTAGAGGAAGTAAAAGTCGTAACAAGGTTTCC</td>
      <td>GOMECC4_27N_Sta1_DCM_A_18S</td>
      <td>962</td>
      <td>GOMECC4_27N_Sta1_DCM_A_occ4e38e8ced9070952b314e1880bede1ca</td>
      <td>Clausocalanus furcatus</td>
      <td>Clausocalanus furcatus</td>
      <td>urn:lsid:marinespecies.org:taxname:104503</td>
      <td>Animalia</td>
      <td>Arthropoda</td>
      <td>Copepoda</td>
      <td>Calanoida</td>
      <td>Clausocalanidae</td>
      <td>Clausocalanus</td>
      <td>Species</td>
      <td>Tourmaline; qiime2-2021.2; naive-bayes classifier, confidence (at lowest specified taxon): 0.999946735, against reference database: PR2 v5.0.1; V9...</td>
      <td>ASV:4e38e8ced9070952b314e1880bede1ca</td>
      <td>MaterialSample</td>
      <td>WoRMS</td>
      <td>DNA sequence reads</td>
      <td>Luke Thompson | Katherine Silliman</td>
      <td>SAMN37516094</td>
      <td>9838</td>
      <td>DNA sequence reads</td>
      <td>SRR26161153 | SAMN37516094 | PRJNA887898</td>
      <td>27N_Sta1</td>
      <td>2021-09-14T11:00-04:00</td>
      <td>49</td>
      <td>USA: Atlantic Ocean, east of Florida (27 N)</td>
      <td>26.997</td>
      <td>-79.618</td>
      <td>CTD rosette</td>
      <td>Atlantic Ocean</td>
      <td>49</td>
      <td>GOMECC4_27N_Sta1_DCM_A</td>
      <td>noaa-aoml-gomecc4</td>
      <td>WGS84</td>
    </tr>
    <tr>
      <th>2</th>
      <td>GCTACTACCGATTGGACGTTTTAGTGAGACATTTGGACTGGGTTAAGATAGTCGCAAGACTACCTTTTCTCCGGAAAGACTTTCAAACTTGAGCGTCTAGAGGAAGTAAAAGTCGTAACAAGGTTTCC</td>
      <td>GOMECC4_27N_Sta1_DCM_A_18S</td>
      <td>1164</td>
      <td>GOMECC4_27N_Sta1_DCM_A_occ2a31e5c01634165da99e7381279baa75</td>
      <td>Eukaryota;Obazoa;Opisthokonta;Metazoa;Arthropoda;Crustacea;Maxillopoda;Acrocalanus;Acrocalanus_sp.;</td>
      <td>Acrocalanus</td>
      <td>urn:lsid:marinespecies.org:taxname:104192</td>
      <td>Animalia</td>
      <td>Arthropoda</td>
      <td>Copepoda</td>
      <td>Calanoida</td>
      <td>Paracalanidae</td>
      <td>Acrocalanus</td>
      <td>Genus</td>
      <td>Tourmaline; qiime2-2021.2; naive-bayes classifier, confidence (at lowest specified taxon): 0.779948049, against reference database: PR2 v5.0.1; V9...</td>
      <td>ASV:2a31e5c01634165da99e7381279baa75</td>
      <td>MaterialSample</td>
      <td>WoRMS</td>
      <td>DNA sequence reads</td>
      <td>Luke Thompson | Katherine Silliman</td>
      <td>SAMN37516094</td>
      <td>9838</td>
      <td>DNA sequence reads</td>
      <td>SRR26161153 | SAMN37516094 | PRJNA887898</td>
      <td>27N_Sta1</td>
      <td>2021-09-14T11:00-04:00</td>
      <td>49</td>
      <td>USA: Atlantic Ocean, east of Florida (27 N)</td>
      <td>26.997</td>
      <td>-79.618</td>
      <td>CTD rosette</td>
      <td>Atlantic Ocean</td>
      <td>49</td>
      <td>GOMECC4_27N_Sta1_DCM_A</td>
      <td>noaa-aoml-gomecc4</td>
      <td>WGS84</td>
    </tr>
    <tr>
      <th>3</th>
      <td>GCTCCTACCGATTGAGTGATCCGGTGAATAATTCGGACTGCAGCAGTGTTCAGTTCCTGAACGTTGCAGCGGAAAGTTTAGTGAACCTTATCACTTAGAGGAAGGAGAAGTCGTAACAAGGTTTCC</td>
      <td>GOMECC4_27N_Sta1_DCM_A_18S</td>
      <td>287</td>
      <td>GOMECC4_27N_Sta1_DCM_A_occecee60339b2fb88ea6d1c8d18054bed4</td>
      <td>Eukaryota;TSAR;Alveolata;Dinoflagellata;Dinophyceae</td>
      <td>Dinophyceae</td>
      <td>urn:lsid:marinespecies.org:taxname:19542</td>
      <td>Chromista</td>
      <td>Myzozoa</td>
      <td>Dinophyceae</td>
      <td>None</td>
      <td>None</td>
      <td>None</td>
      <td>Class</td>
      <td>Tourmaline; qiime2-2021.2; naive-bayes classifier, confidence (at lowest specified taxon): 0.999930607, against reference database: PR2 v5.0.1; V9...</td>
      <td>ASV:ecee60339b2fb88ea6d1c8d18054bed4</td>
      <td>MaterialSample</td>
      <td>WoRMS</td>
      <td>DNA sequence reads</td>
      <td>Luke Thompson | Katherine Silliman</td>
      <td>SAMN37516094</td>
      <td>9838</td>
      <td>DNA sequence reads</td>
      <td>SRR26161153 | SAMN37516094 | PRJNA887898</td>
      <td>27N_Sta1</td>
      <td>2021-09-14T11:00-04:00</td>
      <td>49</td>
      <td>USA: Atlantic Ocean, east of Florida (27 N)</td>
      <td>26.997</td>
      <td>-79.618</td>
      <td>CTD rosette</td>
      <td>Atlantic Ocean</td>
      <td>49</td>
      <td>GOMECC4_27N_Sta1_DCM_A</td>
      <td>noaa-aoml-gomecc4</td>
      <td>WGS84</td>
    </tr>
    <tr>
      <th>4</th>
      <td>GCTCCTACCGATTGAGTGATCCGGTGAATAATTCGGACTGCAGCAATGTTTGGATCCCGAACGTTGCAGCGGAAAGTTTAGTGAACCTTATCACTTAGAGGAAGGAGAAGTCGTAACAAGGTTTCC</td>
      <td>GOMECC4_27N_Sta1_DCM_A_18S</td>
      <td>250</td>
      <td>GOMECC4_27N_Sta1_DCM_A_occfa1f1a97dd4ae7c826009186bad26384</td>
      <td>Eukaryota;TSAR;Alveolata;Dinoflagellata;Dinophyceae;Gymnodiniales;Gymnodiniaceae</td>
      <td>Gymnodiniaceae</td>
      <td>urn:lsid:marinespecies.org:taxname:109410</td>
      <td>Chromista</td>
      <td>Myzozoa</td>
      <td>Dinophyceae</td>
      <td>Gymnodiniales</td>
      <td>Gymnodiniaceae</td>
      <td>None</td>
      <td>Family</td>
      <td>Tourmaline; qiime2-2021.2; naive-bayes classifier, confidence (at lowest specified taxon): 0.98690791, against reference database: PR2 v5.0.1; V9 ...</td>
      <td>ASV:fa1f1a97dd4ae7c826009186bad26384</td>
      <td>MaterialSample</td>
      <td>WoRMS</td>
      <td>DNA sequence reads</td>
      <td>Luke Thompson | Katherine Silliman</td>
      <td>SAMN37516094</td>
      <td>9838</td>
      <td>DNA sequence reads</td>
      <td>SRR26161153 | SAMN37516094 | PRJNA887898</td>
      <td>27N_Sta1</td>
      <td>2021-09-14T11:00-04:00</td>
      <td>49</td>
      <td>USA: Atlantic Ocean, east of Florida (27 N)</td>
      <td>26.997</td>
      <td>-79.618</td>
      <td>CTD rosette</td>
      <td>Atlantic Ocean</td>
      <td>49</td>
      <td>GOMECC4_27N_Sta1_DCM_A</td>
      <td>noaa-aoml-gomecc4</td>
      <td>WGS84</td>
    </tr>
  </tbody>
</table>
</div>




```python
occ18_merged.drop(columns=['DNA_sequence']).to_csv("../processed/occurrence_18S.tsv",sep="\t",index=False)
```

### combine 16s and 18s occurrence


```python
occ18_merged.shape
```




    (146232, 36)




```python
occ_all = pd.concat([occ16_merged,occ18_merged],axis=0, ignore_index=True)
```


```python
occ_all['occurrenceStatus'] = 'present' 
```


```python
occ_all.shape
```




    (311390, 37)




```python
occ_all.drop(columns=['DNA_sequence']).to_csv("../processed/occurrence.csv",index=False)
```

### DNA-derived


```python
dna_dict = dwc_data['dna'].to_dict('index')
```


```python
dna_dict.keys()
```




    dict_keys(['eventID', 'occurrenceID', 'DNA_sequence', 'sop', 'nucl_acid_ext', 'samp_vol_we_dna_ext', 'samp_mat_process', 'nucl_acid_amp', 'target_gene', 'target_subfragment', 'ampliconSize', 'lib_layout', 'pcr_primer_forward', 'pcr_primer_reverse', 'pcr_primer_name_forward', 'pcr_primer_name_reverse', 'pcr_primer_reference', 'pcr_cond', 'seq_meth', 'otu_class_appr', 'otu_seq_comp_appr', 'otu_db', 'env_broad_scale', 'env_local_scale', 'env_medium', 'size_frac', 'concentration', 'concentrationUnit', 'samp_collect_device', 'source_mat_id'])



##### sample_data


```python
# check which dna file terms are in sample_data
for key in dna_dict.keys():
    if dna_dict[key]['AOML_file'] == 'sample_data':
        print(key,dna_dict[key])
```

    samp_vol_we_dna_ext {'AOML_term': 'samp_vol_we_dna_ext', 'AOML_file': 'sample_data', 'DwC_definition': 'Volume (ml) or mass (g) of total collected sample processed for DNA extraction.MIXS:0000111', 'Example': nan}
    samp_mat_process {'AOML_term': 'samp_mat_process', 'AOML_file': 'sample_data', 'DwC_definition': 'Any processing applied to the sample during or after retrieving the sample from environment. This field accepts OBI, for a browser of OBI (v 2018-02-12) terms please see http://purl.bioontology.org/ontology/OBI', 'Example': nan}
    env_broad_scale {'AOML_term': 'env_broad_scale', 'AOML_file': 'sample_data', 'DwC_definition': nan, 'Example': nan}
    env_local_scale {'AOML_term': 'env_local_scale', 'AOML_file': 'sample_data', 'DwC_definition': nan, 'Example': nan}
    env_medium {'AOML_term': 'env_medium', 'AOML_file': 'sample_data', 'DwC_definition': nan, 'Example': nan}
    size_frac {'AOML_term': 'size_frac', 'AOML_file': 'sample_data', 'DwC_definition': 'Filtering pore size used in sample preparation. Examples: 0-0.22 micrometer', 'Example': nan}
    concentration {'AOML_term': 'dna_conc', 'AOML_file': 'sample_data', 'DwC_definition': nan, 'Example': nan}
    concentrationUnit {'AOML_term': 'derived: dna_conc', 'AOML_file': 'sample_data', 'DwC_definition': nan, 'Example': nan}
    samp_collect_device {'AOML_term': 'samp_collect_device', 'AOML_file': 'sample_data', 'DwC_definition': nan, 'Example': nan}
    source_mat_id {'AOML_term': 'source_mat_id', 'AOML_file': 'sample_data', 'DwC_definition': 'used here for the nisken bottle sample', 'Example': nan}



```python
# rename sample_data columns to fit DwC standard
rename_dict = {}
gen = (x for x in dna_dict.keys() if dna_dict[x]['AOML_file'] == 'sample_data')
for x in gen:
    #print(x)
    rename_dict[dna_dict[x]['AOML_term']] = x

gen = (x for x in dna_dict.keys() if dna_dict[x]['AOML_file'] == 'prep_data')
for x in gen:
    #print(x)
    rename_dict[dna_dict[x]['AOML_term']] = x

gen = (x for x in dna_dict.keys() if dna_dict[x]['AOML_file'] == 'analysis_data')
for x in gen:
    #print(x)
    rename_dict[dna_dict[x]['AOML_term']] = x

dna_sample = data['sample_data'].rename(columns=rename_dict).copy()
dna_prep = data['prep_data'].rename(columns=rename_dict).copy()
dna_analysis = data['analysis_data'].rename(columns=rename_dict).copy()

#dna_sample = dna_sample.drop(columns=[col for col in dna_sample if col not in rename_dict.values()])
```


```python
dna_16 = dna_sample[dna_sample['amplicon_sequenced'].str.contains('16S V4-V5')].copy()
dna_16['eventID'] = dna_16['eventID']+"_16S"
dna_16.drop(columns=['amplicon_sequenced'],inplace=True)
dna_16.head()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>eventID</th>
      <th>serial_number</th>
      <th>cruise_id</th>
      <th>line_id</th>
      <th>station</th>
      <th>ctd_bottle_no</th>
      <th>sample_replicate</th>
      <th>source_mat_id</th>
      <th>biological_replicates</th>
      <th>extract_number</th>
      <th>sample_title</th>
      <th>bioproject_accession</th>
      <th>biosample_accession</th>
      <th>metagenome_sequenced</th>
      <th>organism</th>
      <th>collection_date_local</th>
      <th>collection_date</th>
      <th>depth</th>
      <th>env_broad_scale</th>
      <th>env_local_scale</th>
      <th>env_medium</th>
      <th>geo_loc_name</th>
      <th>lat_lon</th>
      <th>decimalLatitude</th>
      <th>decimalLongitude</th>
      <th>...</th>
      <th>ammonium</th>
      <th>carbonate</th>
      <th>diss_inorg_carb</th>
      <th>diss_oxygen</th>
      <th>fluor</th>
      <th>hydrogen_ion</th>
      <th>nitrate</th>
      <th>nitrite</th>
      <th>nitrate_plus_nitrite</th>
      <th>omega_arag</th>
      <th>pco2</th>
      <th>ph</th>
      <th>phosphate</th>
      <th>pressure</th>
      <th>salinity</th>
      <th>samp_store_loc</th>
      <th>samp_store_temp</th>
      <th>silicate</th>
      <th>size_frac_low</th>
      <th>size_frac_up</th>
      <th>temp</th>
      <th>tot_alkalinity</th>
      <th>tot_depth_water_col</th>
      <th>transmittance</th>
      <th>waterBody</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>GOMECC4_27N_Sta1_Deep_A_16S</td>
      <td>GOMECC4_001</td>
      <td>GOMECC-4 (2021)</td>
      <td>27N</td>
      <td>27N_Sta1</td>
      <td>3</td>
      <td>A</td>
      <td>GOMECC4_27N_Sta1_Deep</td>
      <td>GOMECC4_27N_Sta1_Deep_B, GOMECC4_27N_Sta1_Deep_C</td>
      <td>Plate4_52</td>
      <td>Atlantic Ocean seawater sample GOMECC4_27N_Sta1_Deep_A</td>
      <td>PRJNA887898</td>
      <td>SAMN37516091</td>
      <td>planned for FY24</td>
      <td>seawater metagenome</td>
      <td>2021-09-14T11:00-04:00</td>
      <td>2021-09-14T07:00</td>
      <td>618 m</td>
      <td>marine biome [ENVO:00000447]</td>
      <td>marine mesopelagic zone [ENVO:00000213]</td>
      <td>sea water [ENVO:00002149]</td>
      <td>USA: Atlantic Ocean, east of Florida (27 N)</td>
      <td>26.997 N 79.618 W</td>
      <td>26.997</td>
      <td>-79.618</td>
      <td>...</td>
      <td>0.25971 µmol/kg</td>
      <td>88.434 µmol/kg</td>
      <td>2215.45 µmol/kg</td>
      <td>129.44 µmol/kg</td>
      <td>0.0308</td>
      <td>0.0000000142 M</td>
      <td>29.3256 µmol/kg</td>
      <td>0.00391 µmol/kg</td>
      <td>29.3295 µmol/kg</td>
      <td>1.168</td>
      <td>624 µatm</td>
      <td>7.849</td>
      <td>1.94489 µmol/kg</td>
      <td>623 dbar</td>
      <td>34.946 psu</td>
      <td>NOAA/AOML Room 248</td>
      <td>-20 °C</td>
      <td>20.3569 µmol/kg</td>
      <td>no pre-filter</td>
      <td>0.22 µm</td>
      <td>7.479 °C</td>
      <td>2318.9 µmol/kg</td>
      <td>623 m</td>
      <td>4.7221</td>
      <td>Atlantic Ocean</td>
    </tr>
    <tr>
      <th>1</th>
      <td>GOMECC4_27N_Sta1_Deep_B_16S</td>
      <td>GOMECC4_002</td>
      <td>GOMECC-4 (2021)</td>
      <td>27N</td>
      <td>27N_Sta1</td>
      <td>3</td>
      <td>B</td>
      <td>GOMECC4_27N_Sta1_Deep</td>
      <td>GOMECC4_27N_Sta1_Deep_A, GOMECC4_27N_Sta1_Deep_C</td>
      <td>Plate4_60</td>
      <td>Atlantic Ocean seawater sample GOMECC4_27N_Sta1_Deep_B</td>
      <td>PRJNA887898</td>
      <td>SAMN37516092</td>
      <td>planned for FY24</td>
      <td>seawater metagenome</td>
      <td>2021-09-14T11:00-04:00</td>
      <td>2021-09-14T07:00</td>
      <td>618 m</td>
      <td>marine biome [ENVO:00000447]</td>
      <td>marine mesopelagic zone [ENVO:00000213]</td>
      <td>sea water [ENVO:00002149]</td>
      <td>USA: Atlantic Ocean, east of Florida (27 N)</td>
      <td>26.997 N 79.618 W</td>
      <td>26.997</td>
      <td>-79.618</td>
      <td>...</td>
      <td>0.25971 µmol/kg</td>
      <td>88.434 µmol/kg</td>
      <td>2215.45 µmol/kg</td>
      <td>129.44 µmol/kg</td>
      <td>0.0308</td>
      <td>0.0000000142 M</td>
      <td>29.3256 µmol/kg</td>
      <td>0.00391 µmol/kg</td>
      <td>29.3295 µmol/kg</td>
      <td>1.168</td>
      <td>624 µatm</td>
      <td>7.849</td>
      <td>1.94489 µmol/kg</td>
      <td>623 dbar</td>
      <td>34.946 psu</td>
      <td>NOAA/AOML Room 248</td>
      <td>-20 °C</td>
      <td>20.3569 µmol/kg</td>
      <td>no pre-filter</td>
      <td>0.22 µm</td>
      <td>7.479 °C</td>
      <td>2318.9 µmol/kg</td>
      <td>623 m</td>
      <td>4.7221</td>
      <td>Atlantic Ocean</td>
    </tr>
    <tr>
      <th>2</th>
      <td>GOMECC4_27N_Sta1_Deep_C_16S</td>
      <td>GOMECC4_003</td>
      <td>GOMECC-4 (2021)</td>
      <td>27N</td>
      <td>27N_Sta1</td>
      <td>3</td>
      <td>C</td>
      <td>GOMECC4_27N_Sta1_Deep</td>
      <td>GOMECC4_27N_Sta1_Deep_A, GOMECC4_27N_Sta1_Deep_B</td>
      <td>Plate4_62</td>
      <td>Atlantic Ocean seawater sample GOMECC4_27N_Sta1_Deep_C</td>
      <td>PRJNA887898</td>
      <td>SAMN37516093</td>
      <td>planned for FY24</td>
      <td>seawater metagenome</td>
      <td>2021-09-14T11:00-04:00</td>
      <td>2021-09-14T07:00</td>
      <td>618 m</td>
      <td>marine biome [ENVO:00000447]</td>
      <td>marine mesopelagic zone [ENVO:00000213]</td>
      <td>sea water [ENVO:00002149]</td>
      <td>USA: Atlantic Ocean, east of Florida (27 N)</td>
      <td>26.997 N 79.618 W</td>
      <td>26.997</td>
      <td>-79.618</td>
      <td>...</td>
      <td>0.25971 µmol/kg</td>
      <td>88.434 µmol/kg</td>
      <td>2215.45 µmol/kg</td>
      <td>129.44 µmol/kg</td>
      <td>0.0308</td>
      <td>0.0000000142 M</td>
      <td>29.3256 µmol/kg</td>
      <td>0.00391 µmol/kg</td>
      <td>29.3295 µmol/kg</td>
      <td>1.168</td>
      <td>624 µatm</td>
      <td>7.849</td>
      <td>1.94489 µmol/kg</td>
      <td>623 dbar</td>
      <td>34.946 psu</td>
      <td>NOAA/AOML Room 248</td>
      <td>-20 °C</td>
      <td>20.3569 µmol/kg</td>
      <td>no pre-filter</td>
      <td>0.22 µm</td>
      <td>7.479 °C</td>
      <td>2318.9 µmol/kg</td>
      <td>623 m</td>
      <td>4.7221</td>
      <td>Atlantic Ocean</td>
    </tr>
    <tr>
      <th>3</th>
      <td>GOMECC4_27N_Sta1_DCM_A_16S</td>
      <td>GOMECC4_004</td>
      <td>GOMECC-4 (2021)</td>
      <td>27N</td>
      <td>27N_Sta1</td>
      <td>14</td>
      <td>A</td>
      <td>GOMECC4_27N_Sta1_DCM</td>
      <td>GOMECC4_27N_Sta1_DCM_B, GOMECC4_27N_Sta1_DCM_C</td>
      <td>Plate4_53</td>
      <td>Atlantic Ocean seawater sample GOMECC4_27N_Sta1_DCM_A</td>
      <td>PRJNA887898</td>
      <td>SAMN37516094</td>
      <td>planned for FY24</td>
      <td>seawater metagenome</td>
      <td>2021-09-14T11:00-04:00</td>
      <td>2021-09-14T07:00</td>
      <td>49 m</td>
      <td>marine biome [ENVO:00000447]</td>
      <td>marine photic zone [ENVO:00000209]</td>
      <td>sea water [ENVO:00002149]</td>
      <td>USA: Atlantic Ocean, east of Florida (27 N)</td>
      <td>26.997 N 79.618 W</td>
      <td>26.997</td>
      <td>-79.618</td>
      <td>...</td>
      <td>0.32968 µmol/kg</td>
      <td>229.99 µmol/kg</td>
      <td>2033.19 µmol/kg</td>
      <td>193.443 µmol/kg</td>
      <td>0.036</td>
      <td>0.0000000094 M</td>
      <td>0 µmol/kg</td>
      <td>0 µmol/kg</td>
      <td>0 µmol/kg</td>
      <td>3.805</td>
      <td>423 µatm</td>
      <td>8.027</td>
      <td>0.0517 µmol/kg</td>
      <td>49 dbar</td>
      <td>36.325 psu</td>
      <td>NOAA/AOML Room 248</td>
      <td>-20 °C</td>
      <td>1.05635 µmol/kg</td>
      <td>no pre-filter</td>
      <td>0.22 µm</td>
      <td>28.592 °C</td>
      <td>2371 µmol/kg</td>
      <td>623 m</td>
      <td>4.665</td>
      <td>Atlantic Ocean</td>
    </tr>
    <tr>
      <th>4</th>
      <td>GOMECC4_27N_Sta1_DCM_B_16S</td>
      <td>GOMECC4_005</td>
      <td>GOMECC-4 (2021)</td>
      <td>27N</td>
      <td>27N_Sta1</td>
      <td>14</td>
      <td>B</td>
      <td>GOMECC4_27N_Sta1_DCM</td>
      <td>GOMECC4_27N_Sta1_DCM_A, GOMECC4_27N_Sta1_DCM_C</td>
      <td>Plate4_46</td>
      <td>Atlantic Ocean seawater sample GOMECC4_27N_Sta1_DCM_B</td>
      <td>PRJNA887898</td>
      <td>SAMN37516095</td>
      <td>planned for FY24</td>
      <td>seawater metagenome</td>
      <td>2021-09-14T11:00-04:00</td>
      <td>2021-09-14T07:00</td>
      <td>49 m</td>
      <td>marine biome [ENVO:00000447]</td>
      <td>marine photic zone [ENVO:00000209]</td>
      <td>sea water [ENVO:00002149]</td>
      <td>USA: Atlantic Ocean, east of Florida (27 N)</td>
      <td>26.997 N 79.618 W</td>
      <td>26.997</td>
      <td>-79.618</td>
      <td>...</td>
      <td>0.32968 µmol/kg</td>
      <td>229.99 µmol/kg</td>
      <td>2033.19 µmol/kg</td>
      <td>193.443 µmol/kg</td>
      <td>0.036</td>
      <td>0.0000000094 M</td>
      <td>0 µmol/kg</td>
      <td>0 µmol/kg</td>
      <td>0 µmol/kg</td>
      <td>3.805</td>
      <td>423 µatm</td>
      <td>8.027</td>
      <td>0.0517 µmol/kg</td>
      <td>49 dbar</td>
      <td>36.325 psu</td>
      <td>NOAA/AOML Room 248</td>
      <td>-20 °C</td>
      <td>1.05635 µmol/kg</td>
      <td>no pre-filter</td>
      <td>0.22 µm</td>
      <td>28.592 °C</td>
      <td>2371 µmol/kg</td>
      <td>623 m</td>
      <td>4.665</td>
      <td>Atlantic Ocean</td>
    </tr>
  </tbody>
</table>
<p>5 rows × 75 columns</p>
</div>




```python
dna_18 = dna_sample[dna_sample['amplicon_sequenced'].str.contains('18S V9')].copy()
dna_18['eventID'] = dna_18['eventID']+"_18S"
dna_18.drop(columns=['amplicon_sequenced'],inplace=True)
dna_18.head()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>eventID</th>
      <th>serial_number</th>
      <th>cruise_id</th>
      <th>line_id</th>
      <th>station</th>
      <th>ctd_bottle_no</th>
      <th>sample_replicate</th>
      <th>source_mat_id</th>
      <th>biological_replicates</th>
      <th>extract_number</th>
      <th>sample_title</th>
      <th>bioproject_accession</th>
      <th>biosample_accession</th>
      <th>metagenome_sequenced</th>
      <th>organism</th>
      <th>collection_date_local</th>
      <th>collection_date</th>
      <th>depth</th>
      <th>env_broad_scale</th>
      <th>env_local_scale</th>
      <th>env_medium</th>
      <th>geo_loc_name</th>
      <th>lat_lon</th>
      <th>decimalLatitude</th>
      <th>decimalLongitude</th>
      <th>...</th>
      <th>ammonium</th>
      <th>carbonate</th>
      <th>diss_inorg_carb</th>
      <th>diss_oxygen</th>
      <th>fluor</th>
      <th>hydrogen_ion</th>
      <th>nitrate</th>
      <th>nitrite</th>
      <th>nitrate_plus_nitrite</th>
      <th>omega_arag</th>
      <th>pco2</th>
      <th>ph</th>
      <th>phosphate</th>
      <th>pressure</th>
      <th>salinity</th>
      <th>samp_store_loc</th>
      <th>samp_store_temp</th>
      <th>silicate</th>
      <th>size_frac_low</th>
      <th>size_frac_up</th>
      <th>temp</th>
      <th>tot_alkalinity</th>
      <th>tot_depth_water_col</th>
      <th>transmittance</th>
      <th>waterBody</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>GOMECC4_27N_Sta1_Deep_A_18S</td>
      <td>GOMECC4_001</td>
      <td>GOMECC-4 (2021)</td>
      <td>27N</td>
      <td>27N_Sta1</td>
      <td>3</td>
      <td>A</td>
      <td>GOMECC4_27N_Sta1_Deep</td>
      <td>GOMECC4_27N_Sta1_Deep_B, GOMECC4_27N_Sta1_Deep_C</td>
      <td>Plate4_52</td>
      <td>Atlantic Ocean seawater sample GOMECC4_27N_Sta1_Deep_A</td>
      <td>PRJNA887898</td>
      <td>SAMN37516091</td>
      <td>planned for FY24</td>
      <td>seawater metagenome</td>
      <td>2021-09-14T11:00-04:00</td>
      <td>2021-09-14T07:00</td>
      <td>618 m</td>
      <td>marine biome [ENVO:00000447]</td>
      <td>marine mesopelagic zone [ENVO:00000213]</td>
      <td>sea water [ENVO:00002149]</td>
      <td>USA: Atlantic Ocean, east of Florida (27 N)</td>
      <td>26.997 N 79.618 W</td>
      <td>26.997</td>
      <td>-79.618</td>
      <td>...</td>
      <td>0.25971 µmol/kg</td>
      <td>88.434 µmol/kg</td>
      <td>2215.45 µmol/kg</td>
      <td>129.44 µmol/kg</td>
      <td>0.0308</td>
      <td>0.0000000142 M</td>
      <td>29.3256 µmol/kg</td>
      <td>0.00391 µmol/kg</td>
      <td>29.3295 µmol/kg</td>
      <td>1.168</td>
      <td>624 µatm</td>
      <td>7.849</td>
      <td>1.94489 µmol/kg</td>
      <td>623 dbar</td>
      <td>34.946 psu</td>
      <td>NOAA/AOML Room 248</td>
      <td>-20 °C</td>
      <td>20.3569 µmol/kg</td>
      <td>no pre-filter</td>
      <td>0.22 µm</td>
      <td>7.479 °C</td>
      <td>2318.9 µmol/kg</td>
      <td>623 m</td>
      <td>4.7221</td>
      <td>Atlantic Ocean</td>
    </tr>
    <tr>
      <th>1</th>
      <td>GOMECC4_27N_Sta1_Deep_B_18S</td>
      <td>GOMECC4_002</td>
      <td>GOMECC-4 (2021)</td>
      <td>27N</td>
      <td>27N_Sta1</td>
      <td>3</td>
      <td>B</td>
      <td>GOMECC4_27N_Sta1_Deep</td>
      <td>GOMECC4_27N_Sta1_Deep_A, GOMECC4_27N_Sta1_Deep_C</td>
      <td>Plate4_60</td>
      <td>Atlantic Ocean seawater sample GOMECC4_27N_Sta1_Deep_B</td>
      <td>PRJNA887898</td>
      <td>SAMN37516092</td>
      <td>planned for FY24</td>
      <td>seawater metagenome</td>
      <td>2021-09-14T11:00-04:00</td>
      <td>2021-09-14T07:00</td>
      <td>618 m</td>
      <td>marine biome [ENVO:00000447]</td>
      <td>marine mesopelagic zone [ENVO:00000213]</td>
      <td>sea water [ENVO:00002149]</td>
      <td>USA: Atlantic Ocean, east of Florida (27 N)</td>
      <td>26.997 N 79.618 W</td>
      <td>26.997</td>
      <td>-79.618</td>
      <td>...</td>
      <td>0.25971 µmol/kg</td>
      <td>88.434 µmol/kg</td>
      <td>2215.45 µmol/kg</td>
      <td>129.44 µmol/kg</td>
      <td>0.0308</td>
      <td>0.0000000142 M</td>
      <td>29.3256 µmol/kg</td>
      <td>0.00391 µmol/kg</td>
      <td>29.3295 µmol/kg</td>
      <td>1.168</td>
      <td>624 µatm</td>
      <td>7.849</td>
      <td>1.94489 µmol/kg</td>
      <td>623 dbar</td>
      <td>34.946 psu</td>
      <td>NOAA/AOML Room 248</td>
      <td>-20 °C</td>
      <td>20.3569 µmol/kg</td>
      <td>no pre-filter</td>
      <td>0.22 µm</td>
      <td>7.479 °C</td>
      <td>2318.9 µmol/kg</td>
      <td>623 m</td>
      <td>4.7221</td>
      <td>Atlantic Ocean</td>
    </tr>
    <tr>
      <th>2</th>
      <td>GOMECC4_27N_Sta1_Deep_C_18S</td>
      <td>GOMECC4_003</td>
      <td>GOMECC-4 (2021)</td>
      <td>27N</td>
      <td>27N_Sta1</td>
      <td>3</td>
      <td>C</td>
      <td>GOMECC4_27N_Sta1_Deep</td>
      <td>GOMECC4_27N_Sta1_Deep_A, GOMECC4_27N_Sta1_Deep_B</td>
      <td>Plate4_62</td>
      <td>Atlantic Ocean seawater sample GOMECC4_27N_Sta1_Deep_C</td>
      <td>PRJNA887898</td>
      <td>SAMN37516093</td>
      <td>planned for FY24</td>
      <td>seawater metagenome</td>
      <td>2021-09-14T11:00-04:00</td>
      <td>2021-09-14T07:00</td>
      <td>618 m</td>
      <td>marine biome [ENVO:00000447]</td>
      <td>marine mesopelagic zone [ENVO:00000213]</td>
      <td>sea water [ENVO:00002149]</td>
      <td>USA: Atlantic Ocean, east of Florida (27 N)</td>
      <td>26.997 N 79.618 W</td>
      <td>26.997</td>
      <td>-79.618</td>
      <td>...</td>
      <td>0.25971 µmol/kg</td>
      <td>88.434 µmol/kg</td>
      <td>2215.45 µmol/kg</td>
      <td>129.44 µmol/kg</td>
      <td>0.0308</td>
      <td>0.0000000142 M</td>
      <td>29.3256 µmol/kg</td>
      <td>0.00391 µmol/kg</td>
      <td>29.3295 µmol/kg</td>
      <td>1.168</td>
      <td>624 µatm</td>
      <td>7.849</td>
      <td>1.94489 µmol/kg</td>
      <td>623 dbar</td>
      <td>34.946 psu</td>
      <td>NOAA/AOML Room 248</td>
      <td>-20 °C</td>
      <td>20.3569 µmol/kg</td>
      <td>no pre-filter</td>
      <td>0.22 µm</td>
      <td>7.479 °C</td>
      <td>2318.9 µmol/kg</td>
      <td>623 m</td>
      <td>4.7221</td>
      <td>Atlantic Ocean</td>
    </tr>
    <tr>
      <th>3</th>
      <td>GOMECC4_27N_Sta1_DCM_A_18S</td>
      <td>GOMECC4_004</td>
      <td>GOMECC-4 (2021)</td>
      <td>27N</td>
      <td>27N_Sta1</td>
      <td>14</td>
      <td>A</td>
      <td>GOMECC4_27N_Sta1_DCM</td>
      <td>GOMECC4_27N_Sta1_DCM_B, GOMECC4_27N_Sta1_DCM_C</td>
      <td>Plate4_53</td>
      <td>Atlantic Ocean seawater sample GOMECC4_27N_Sta1_DCM_A</td>
      <td>PRJNA887898</td>
      <td>SAMN37516094</td>
      <td>planned for FY24</td>
      <td>seawater metagenome</td>
      <td>2021-09-14T11:00-04:00</td>
      <td>2021-09-14T07:00</td>
      <td>49 m</td>
      <td>marine biome [ENVO:00000447]</td>
      <td>marine photic zone [ENVO:00000209]</td>
      <td>sea water [ENVO:00002149]</td>
      <td>USA: Atlantic Ocean, east of Florida (27 N)</td>
      <td>26.997 N 79.618 W</td>
      <td>26.997</td>
      <td>-79.618</td>
      <td>...</td>
      <td>0.32968 µmol/kg</td>
      <td>229.99 µmol/kg</td>
      <td>2033.19 µmol/kg</td>
      <td>193.443 µmol/kg</td>
      <td>0.036</td>
      <td>0.0000000094 M</td>
      <td>0 µmol/kg</td>
      <td>0 µmol/kg</td>
      <td>0 µmol/kg</td>
      <td>3.805</td>
      <td>423 µatm</td>
      <td>8.027</td>
      <td>0.0517 µmol/kg</td>
      <td>49 dbar</td>
      <td>36.325 psu</td>
      <td>NOAA/AOML Room 248</td>
      <td>-20 °C</td>
      <td>1.05635 µmol/kg</td>
      <td>no pre-filter</td>
      <td>0.22 µm</td>
      <td>28.592 °C</td>
      <td>2371 µmol/kg</td>
      <td>623 m</td>
      <td>4.665</td>
      <td>Atlantic Ocean</td>
    </tr>
    <tr>
      <th>4</th>
      <td>GOMECC4_27N_Sta1_DCM_B_18S</td>
      <td>GOMECC4_005</td>
      <td>GOMECC-4 (2021)</td>
      <td>27N</td>
      <td>27N_Sta1</td>
      <td>14</td>
      <td>B</td>
      <td>GOMECC4_27N_Sta1_DCM</td>
      <td>GOMECC4_27N_Sta1_DCM_A, GOMECC4_27N_Sta1_DCM_C</td>
      <td>Plate4_46</td>
      <td>Atlantic Ocean seawater sample GOMECC4_27N_Sta1_DCM_B</td>
      <td>PRJNA887898</td>
      <td>SAMN37516095</td>
      <td>planned for FY24</td>
      <td>seawater metagenome</td>
      <td>2021-09-14T11:00-04:00</td>
      <td>2021-09-14T07:00</td>
      <td>49 m</td>
      <td>marine biome [ENVO:00000447]</td>
      <td>marine photic zone [ENVO:00000209]</td>
      <td>sea water [ENVO:00002149]</td>
      <td>USA: Atlantic Ocean, east of Florida (27 N)</td>
      <td>26.997 N 79.618 W</td>
      <td>26.997</td>
      <td>-79.618</td>
      <td>...</td>
      <td>0.32968 µmol/kg</td>
      <td>229.99 µmol/kg</td>
      <td>2033.19 µmol/kg</td>
      <td>193.443 µmol/kg</td>
      <td>0.036</td>
      <td>0.0000000094 M</td>
      <td>0 µmol/kg</td>
      <td>0 µmol/kg</td>
      <td>0 µmol/kg</td>
      <td>3.805</td>
      <td>423 µatm</td>
      <td>8.027</td>
      <td>0.0517 µmol/kg</td>
      <td>49 dbar</td>
      <td>36.325 psu</td>
      <td>NOAA/AOML Room 248</td>
      <td>-20 °C</td>
      <td>1.05635 µmol/kg</td>
      <td>no pre-filter</td>
      <td>0.22 µm</td>
      <td>28.592 °C</td>
      <td>2371 µmol/kg</td>
      <td>623 m</td>
      <td>4.665</td>
      <td>Atlantic Ocean</td>
    </tr>
  </tbody>
</table>
<p>5 rows × 75 columns</p>
</div>




```python
dna_sample = pd.concat([dna_16,dna_18],axis=0,ignore_index=True)
dna_sample.head()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>eventID</th>
      <th>serial_number</th>
      <th>cruise_id</th>
      <th>line_id</th>
      <th>station</th>
      <th>ctd_bottle_no</th>
      <th>sample_replicate</th>
      <th>source_mat_id</th>
      <th>biological_replicates</th>
      <th>extract_number</th>
      <th>sample_title</th>
      <th>bioproject_accession</th>
      <th>biosample_accession</th>
      <th>metagenome_sequenced</th>
      <th>organism</th>
      <th>collection_date_local</th>
      <th>collection_date</th>
      <th>depth</th>
      <th>env_broad_scale</th>
      <th>env_local_scale</th>
      <th>env_medium</th>
      <th>geo_loc_name</th>
      <th>lat_lon</th>
      <th>decimalLatitude</th>
      <th>decimalLongitude</th>
      <th>...</th>
      <th>ammonium</th>
      <th>carbonate</th>
      <th>diss_inorg_carb</th>
      <th>diss_oxygen</th>
      <th>fluor</th>
      <th>hydrogen_ion</th>
      <th>nitrate</th>
      <th>nitrite</th>
      <th>nitrate_plus_nitrite</th>
      <th>omega_arag</th>
      <th>pco2</th>
      <th>ph</th>
      <th>phosphate</th>
      <th>pressure</th>
      <th>salinity</th>
      <th>samp_store_loc</th>
      <th>samp_store_temp</th>
      <th>silicate</th>
      <th>size_frac_low</th>
      <th>size_frac_up</th>
      <th>temp</th>
      <th>tot_alkalinity</th>
      <th>tot_depth_water_col</th>
      <th>transmittance</th>
      <th>waterBody</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>GOMECC4_27N_Sta1_Deep_A_16S</td>
      <td>GOMECC4_001</td>
      <td>GOMECC-4 (2021)</td>
      <td>27N</td>
      <td>27N_Sta1</td>
      <td>3</td>
      <td>A</td>
      <td>GOMECC4_27N_Sta1_Deep</td>
      <td>GOMECC4_27N_Sta1_Deep_B, GOMECC4_27N_Sta1_Deep_C</td>
      <td>Plate4_52</td>
      <td>Atlantic Ocean seawater sample GOMECC4_27N_Sta1_Deep_A</td>
      <td>PRJNA887898</td>
      <td>SAMN37516091</td>
      <td>planned for FY24</td>
      <td>seawater metagenome</td>
      <td>2021-09-14T11:00-04:00</td>
      <td>2021-09-14T07:00</td>
      <td>618 m</td>
      <td>marine biome [ENVO:00000447]</td>
      <td>marine mesopelagic zone [ENVO:00000213]</td>
      <td>sea water [ENVO:00002149]</td>
      <td>USA: Atlantic Ocean, east of Florida (27 N)</td>
      <td>26.997 N 79.618 W</td>
      <td>26.997</td>
      <td>-79.618</td>
      <td>...</td>
      <td>0.25971 µmol/kg</td>
      <td>88.434 µmol/kg</td>
      <td>2215.45 µmol/kg</td>
      <td>129.44 µmol/kg</td>
      <td>0.0308</td>
      <td>0.0000000142 M</td>
      <td>29.3256 µmol/kg</td>
      <td>0.00391 µmol/kg</td>
      <td>29.3295 µmol/kg</td>
      <td>1.168</td>
      <td>624 µatm</td>
      <td>7.849</td>
      <td>1.94489 µmol/kg</td>
      <td>623 dbar</td>
      <td>34.946 psu</td>
      <td>NOAA/AOML Room 248</td>
      <td>-20 °C</td>
      <td>20.3569 µmol/kg</td>
      <td>no pre-filter</td>
      <td>0.22 µm</td>
      <td>7.479 °C</td>
      <td>2318.9 µmol/kg</td>
      <td>623 m</td>
      <td>4.7221</td>
      <td>Atlantic Ocean</td>
    </tr>
    <tr>
      <th>1</th>
      <td>GOMECC4_27N_Sta1_Deep_B_16S</td>
      <td>GOMECC4_002</td>
      <td>GOMECC-4 (2021)</td>
      <td>27N</td>
      <td>27N_Sta1</td>
      <td>3</td>
      <td>B</td>
      <td>GOMECC4_27N_Sta1_Deep</td>
      <td>GOMECC4_27N_Sta1_Deep_A, GOMECC4_27N_Sta1_Deep_C</td>
      <td>Plate4_60</td>
      <td>Atlantic Ocean seawater sample GOMECC4_27N_Sta1_Deep_B</td>
      <td>PRJNA887898</td>
      <td>SAMN37516092</td>
      <td>planned for FY24</td>
      <td>seawater metagenome</td>
      <td>2021-09-14T11:00-04:00</td>
      <td>2021-09-14T07:00</td>
      <td>618 m</td>
      <td>marine biome [ENVO:00000447]</td>
      <td>marine mesopelagic zone [ENVO:00000213]</td>
      <td>sea water [ENVO:00002149]</td>
      <td>USA: Atlantic Ocean, east of Florida (27 N)</td>
      <td>26.997 N 79.618 W</td>
      <td>26.997</td>
      <td>-79.618</td>
      <td>...</td>
      <td>0.25971 µmol/kg</td>
      <td>88.434 µmol/kg</td>
      <td>2215.45 µmol/kg</td>
      <td>129.44 µmol/kg</td>
      <td>0.0308</td>
      <td>0.0000000142 M</td>
      <td>29.3256 µmol/kg</td>
      <td>0.00391 µmol/kg</td>
      <td>29.3295 µmol/kg</td>
      <td>1.168</td>
      <td>624 µatm</td>
      <td>7.849</td>
      <td>1.94489 µmol/kg</td>
      <td>623 dbar</td>
      <td>34.946 psu</td>
      <td>NOAA/AOML Room 248</td>
      <td>-20 °C</td>
      <td>20.3569 µmol/kg</td>
      <td>no pre-filter</td>
      <td>0.22 µm</td>
      <td>7.479 °C</td>
      <td>2318.9 µmol/kg</td>
      <td>623 m</td>
      <td>4.7221</td>
      <td>Atlantic Ocean</td>
    </tr>
    <tr>
      <th>2</th>
      <td>GOMECC4_27N_Sta1_Deep_C_16S</td>
      <td>GOMECC4_003</td>
      <td>GOMECC-4 (2021)</td>
      <td>27N</td>
      <td>27N_Sta1</td>
      <td>3</td>
      <td>C</td>
      <td>GOMECC4_27N_Sta1_Deep</td>
      <td>GOMECC4_27N_Sta1_Deep_A, GOMECC4_27N_Sta1_Deep_B</td>
      <td>Plate4_62</td>
      <td>Atlantic Ocean seawater sample GOMECC4_27N_Sta1_Deep_C</td>
      <td>PRJNA887898</td>
      <td>SAMN37516093</td>
      <td>planned for FY24</td>
      <td>seawater metagenome</td>
      <td>2021-09-14T11:00-04:00</td>
      <td>2021-09-14T07:00</td>
      <td>618 m</td>
      <td>marine biome [ENVO:00000447]</td>
      <td>marine mesopelagic zone [ENVO:00000213]</td>
      <td>sea water [ENVO:00002149]</td>
      <td>USA: Atlantic Ocean, east of Florida (27 N)</td>
      <td>26.997 N 79.618 W</td>
      <td>26.997</td>
      <td>-79.618</td>
      <td>...</td>
      <td>0.25971 µmol/kg</td>
      <td>88.434 µmol/kg</td>
      <td>2215.45 µmol/kg</td>
      <td>129.44 µmol/kg</td>
      <td>0.0308</td>
      <td>0.0000000142 M</td>
      <td>29.3256 µmol/kg</td>
      <td>0.00391 µmol/kg</td>
      <td>29.3295 µmol/kg</td>
      <td>1.168</td>
      <td>624 µatm</td>
      <td>7.849</td>
      <td>1.94489 µmol/kg</td>
      <td>623 dbar</td>
      <td>34.946 psu</td>
      <td>NOAA/AOML Room 248</td>
      <td>-20 °C</td>
      <td>20.3569 µmol/kg</td>
      <td>no pre-filter</td>
      <td>0.22 µm</td>
      <td>7.479 °C</td>
      <td>2318.9 µmol/kg</td>
      <td>623 m</td>
      <td>4.7221</td>
      <td>Atlantic Ocean</td>
    </tr>
    <tr>
      <th>3</th>
      <td>GOMECC4_27N_Sta1_DCM_A_16S</td>
      <td>GOMECC4_004</td>
      <td>GOMECC-4 (2021)</td>
      <td>27N</td>
      <td>27N_Sta1</td>
      <td>14</td>
      <td>A</td>
      <td>GOMECC4_27N_Sta1_DCM</td>
      <td>GOMECC4_27N_Sta1_DCM_B, GOMECC4_27N_Sta1_DCM_C</td>
      <td>Plate4_53</td>
      <td>Atlantic Ocean seawater sample GOMECC4_27N_Sta1_DCM_A</td>
      <td>PRJNA887898</td>
      <td>SAMN37516094</td>
      <td>planned for FY24</td>
      <td>seawater metagenome</td>
      <td>2021-09-14T11:00-04:00</td>
      <td>2021-09-14T07:00</td>
      <td>49 m</td>
      <td>marine biome [ENVO:00000447]</td>
      <td>marine photic zone [ENVO:00000209]</td>
      <td>sea water [ENVO:00002149]</td>
      <td>USA: Atlantic Ocean, east of Florida (27 N)</td>
      <td>26.997 N 79.618 W</td>
      <td>26.997</td>
      <td>-79.618</td>
      <td>...</td>
      <td>0.32968 µmol/kg</td>
      <td>229.99 µmol/kg</td>
      <td>2033.19 µmol/kg</td>
      <td>193.443 µmol/kg</td>
      <td>0.036</td>
      <td>0.0000000094 M</td>
      <td>0 µmol/kg</td>
      <td>0 µmol/kg</td>
      <td>0 µmol/kg</td>
      <td>3.805</td>
      <td>423 µatm</td>
      <td>8.027</td>
      <td>0.0517 µmol/kg</td>
      <td>49 dbar</td>
      <td>36.325 psu</td>
      <td>NOAA/AOML Room 248</td>
      <td>-20 °C</td>
      <td>1.05635 µmol/kg</td>
      <td>no pre-filter</td>
      <td>0.22 µm</td>
      <td>28.592 °C</td>
      <td>2371 µmol/kg</td>
      <td>623 m</td>
      <td>4.665</td>
      <td>Atlantic Ocean</td>
    </tr>
    <tr>
      <th>4</th>
      <td>GOMECC4_27N_Sta1_DCM_B_16S</td>
      <td>GOMECC4_005</td>
      <td>GOMECC-4 (2021)</td>
      <td>27N</td>
      <td>27N_Sta1</td>
      <td>14</td>
      <td>B</td>
      <td>GOMECC4_27N_Sta1_DCM</td>
      <td>GOMECC4_27N_Sta1_DCM_A, GOMECC4_27N_Sta1_DCM_C</td>
      <td>Plate4_46</td>
      <td>Atlantic Ocean seawater sample GOMECC4_27N_Sta1_DCM_B</td>
      <td>PRJNA887898</td>
      <td>SAMN37516095</td>
      <td>planned for FY24</td>
      <td>seawater metagenome</td>
      <td>2021-09-14T11:00-04:00</td>
      <td>2021-09-14T07:00</td>
      <td>49 m</td>
      <td>marine biome [ENVO:00000447]</td>
      <td>marine photic zone [ENVO:00000209]</td>
      <td>sea water [ENVO:00002149]</td>
      <td>USA: Atlantic Ocean, east of Florida (27 N)</td>
      <td>26.997 N 79.618 W</td>
      <td>26.997</td>
      <td>-79.618</td>
      <td>...</td>
      <td>0.32968 µmol/kg</td>
      <td>229.99 µmol/kg</td>
      <td>2033.19 µmol/kg</td>
      <td>193.443 µmol/kg</td>
      <td>0.036</td>
      <td>0.0000000094 M</td>
      <td>0 µmol/kg</td>
      <td>0 µmol/kg</td>
      <td>0 µmol/kg</td>
      <td>3.805</td>
      <td>423 µatm</td>
      <td>8.027</td>
      <td>0.0517 µmol/kg</td>
      <td>49 dbar</td>
      <td>36.325 psu</td>
      <td>NOAA/AOML Room 248</td>
      <td>-20 °C</td>
      <td>1.05635 µmol/kg</td>
      <td>no pre-filter</td>
      <td>0.22 µm</td>
      <td>28.592 °C</td>
      <td>2371 µmol/kg</td>
      <td>623 m</td>
      <td>4.665</td>
      <td>Atlantic Ocean</td>
    </tr>
  </tbody>
</table>
<p>5 rows × 75 columns</p>
</div>




```python
prep_16 = dna_prep[dna_prep['amplicon_sequenced'].str.contains('16S V4-V5')].copy()
prep_16['eventID'] = prep_16['eventID']+"_16S"
prep_16.head()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>eventID</th>
      <th>library_id</th>
      <th>title</th>
      <th>library_strategy</th>
      <th>library_source</th>
      <th>library_selection</th>
      <th>lib_layout</th>
      <th>platform</th>
      <th>instrument_model</th>
      <th>design_description</th>
      <th>filetype</th>
      <th>filename</th>
      <th>filename2</th>
      <th>biosample_accession</th>
      <th>sra_accession</th>
      <th>seq_meth</th>
      <th>nucl_acid_ext</th>
      <th>amplicon_sequenced</th>
      <th>target_gene</th>
      <th>target_subfragment</th>
      <th>pcr_primer_forward</th>
      <th>pcr_primer_reverse</th>
      <th>pcr_primer_name_forward</th>
      <th>pcr_primer_name_reverse</th>
      <th>pcr_primer_reference</th>
      <th>pcr_cond</th>
      <th>nucl_acid_amp</th>
      <th>adapters</th>
      <th>mid_barcode</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>4</th>
      <td>GOMECC4_BROWNSVILLE_Sta66_DCM_B_16S</td>
      <td>GOMECC16S_Plate1_1</td>
      <td>16S amplicon metabarcoding of marine metagenome: Gulf of Mexico (USA)</td>
      <td>AMPLICON</td>
      <td>METAGENOMIC</td>
      <td>PCR</td>
      <td>paired</td>
      <td>ILLUMINA</td>
      <td>Illumina MiSeq</td>
      <td>Samples were collected and filtered onto Sterivex 0.22 um cartridge filters. DNA was extracted from Sterivex by adding lysis buffer and magnetic b...</td>
      <td>fastq</td>
      <td>GOMECC16S_Plate1_1_S1_L001_R1_001.fastq.gz</td>
      <td>GOMECC16S_Plate1_1_S1_L001_R2_001.fastq.gz</td>
      <td>SAMN37516307</td>
      <td>SRR26148474</td>
      <td>Illumina MiSeq 2x250</td>
      <td>https://github.com/aomlomics/protocols/blob/main/protocol_DNA_extraction_Sterivex.md</td>
      <td>16S V4-V5</td>
      <td>16S rRNA</td>
      <td>V4-V5</td>
      <td>GTGYCAGCMGCCGCGGTAA</td>
      <td>CCGYCAATTYMTTTRAGTTT</td>
      <td>515F-Y</td>
      <td>926R</td>
      <td>10.1111/1462-2920.13023</td>
      <td>initial denaturation:95_2;denaturation:95_0.75;annealing:50_0.75;elongation:68_1.5;final elongation:68_5;25</td>
      <td>10.1111/1462-2920.13023</td>
      <td>ACACTGACGACATGGTTCTACA;TACGGTAGCAGAGACTTGGTCT</td>
      <td>missing: not provided</td>
    </tr>
    <tr>
      <th>6</th>
      <td>GOMECC4_GALVESTON_Sta54_DCM_B_16S</td>
      <td>GOMECC16S_Plate1_10</td>
      <td>16S amplicon metabarcoding of marine metagenome: Gulf of Mexico (USA)</td>
      <td>AMPLICON</td>
      <td>METAGENOMIC</td>
      <td>PCR</td>
      <td>paired</td>
      <td>ILLUMINA</td>
      <td>Illumina MiSeq</td>
      <td>Samples were collected and filtered onto Sterivex 0.22 um cartridge filters. DNA was extracted from Sterivex by adding lysis buffer and magnetic b...</td>
      <td>fastq</td>
      <td>GOMECC16S_Plate1_10_S10_L001_R1_001.fastq.gz</td>
      <td>GOMECC16S_Plate1_10_S10_L001_R2_001.fastq.gz</td>
      <td>SAMN37516268</td>
      <td>SRR26148413</td>
      <td>Illumina MiSeq 2x250</td>
      <td>https://github.com/aomlomics/protocols/blob/main/protocol_DNA_extraction_Sterivex.md</td>
      <td>16S V4-V5</td>
      <td>16S rRNA</td>
      <td>V4-V5</td>
      <td>GTGYCAGCMGCCGCGGTAA</td>
      <td>CCGYCAATTYMTTTRAGTTT</td>
      <td>515F-Y</td>
      <td>926R</td>
      <td>10.1111/1462-2920.13023</td>
      <td>initial denaturation:95_2;denaturation:95_0.75;annealing:50_0.75;elongation:68_1.5;final elongation:68_5;25</td>
      <td>10.1111/1462-2920.13023</td>
      <td>ACACTGACGACATGGTTCTACA;TACGGTAGCAGAGACTTGGTCT</td>
      <td>missing: not provided</td>
    </tr>
    <tr>
      <th>8</th>
      <td>GOMECC4_GALVESTON_Sta54_Deep_A_16S</td>
      <td>GOMECC16S_Plate1_11</td>
      <td>16S amplicon metabarcoding of marine metagenome: Gulf of Mexico (USA)</td>
      <td>AMPLICON</td>
      <td>METAGENOMIC</td>
      <td>PCR</td>
      <td>paired</td>
      <td>ILLUMINA</td>
      <td>Illumina MiSeq</td>
      <td>Samples were collected and filtered onto Sterivex 0.22 um cartridge filters. DNA was extracted from Sterivex by adding lysis buffer and magnetic b...</td>
      <td>fastq</td>
      <td>GOMECC16S_Plate1_11_S11_L001_R1_001.fastq.gz</td>
      <td>GOMECC16S_Plate1_11_S11_L001_R2_001.fastq.gz</td>
      <td>SAMN37516264</td>
      <td>SRR26148140</td>
      <td>Illumina MiSeq 2x250</td>
      <td>https://github.com/aomlomics/protocols/blob/main/protocol_DNA_extraction_Sterivex.md</td>
      <td>16S V4-V5</td>
      <td>16S rRNA</td>
      <td>V4-V5</td>
      <td>GTGYCAGCMGCCGCGGTAA</td>
      <td>CCGYCAATTYMTTTRAGTTT</td>
      <td>515F-Y</td>
      <td>926R</td>
      <td>10.1111/1462-2920.13023</td>
      <td>initial denaturation:95_2;denaturation:95_0.75;annealing:50_0.75;elongation:68_1.5;final elongation:68_5;25</td>
      <td>10.1111/1462-2920.13023</td>
      <td>ACACTGACGACATGGTTCTACA;TACGGTAGCAGAGACTTGGTCT</td>
      <td>missing: not provided</td>
    </tr>
    <tr>
      <th>10</th>
      <td>GOMECC4_GALVESTON_Sta49_Deep_A_16S</td>
      <td>GOMECC16S_Plate1_12</td>
      <td>16S amplicon metabarcoding of marine metagenome: Gulf of Mexico (USA)</td>
      <td>AMPLICON</td>
      <td>METAGENOMIC</td>
      <td>PCR</td>
      <td>paired</td>
      <td>ILLUMINA</td>
      <td>Illumina MiSeq</td>
      <td>Samples were collected and filtered onto Sterivex 0.22 um cartridge filters. DNA was extracted from Sterivex by adding lysis buffer and magnetic b...</td>
      <td>fastq</td>
      <td>GOMECC16S_Plate1_12_S12_L001_R1_001.fastq.gz</td>
      <td>GOMECC16S_Plate1_12_S12_L001_R2_001.fastq.gz</td>
      <td>SAMN37516246</td>
      <td>SRR26148197</td>
      <td>Illumina MiSeq 2x250</td>
      <td>https://github.com/aomlomics/protocols/blob/main/protocol_DNA_extraction_Sterivex.md</td>
      <td>16S V4-V5</td>
      <td>16S rRNA</td>
      <td>V4-V5</td>
      <td>GTGYCAGCMGCCGCGGTAA</td>
      <td>CCGYCAATTYMTTTRAGTTT</td>
      <td>515F-Y</td>
      <td>926R</td>
      <td>10.1111/1462-2920.13023</td>
      <td>initial denaturation:95_2;denaturation:95_0.75;annealing:50_0.75;elongation:68_1.5;final elongation:68_5;25</td>
      <td>10.1111/1462-2920.13023</td>
      <td>ACACTGACGACATGGTTCTACA;TACGGTAGCAGAGACTTGGTCT</td>
      <td>missing: not provided</td>
    </tr>
    <tr>
      <th>12</th>
      <td>GOMECC4_BROWNSVILLE_Sta66_DCM_C_16S</td>
      <td>GOMECC16S_Plate1_13</td>
      <td>16S amplicon metabarcoding of marine metagenome: Gulf of Mexico (USA)</td>
      <td>AMPLICON</td>
      <td>METAGENOMIC</td>
      <td>PCR</td>
      <td>paired</td>
      <td>ILLUMINA</td>
      <td>Illumina MiSeq</td>
      <td>Samples were collected and filtered onto Sterivex 0.22 um cartridge filters. DNA was extracted from Sterivex by adding lysis buffer and magnetic b...</td>
      <td>fastq</td>
      <td>GOMECC16S_Plate1_13_S13_L001_R1_001.fastq.gz</td>
      <td>GOMECC16S_Plate1_13_S13_L001_R2_001.fastq.gz</td>
      <td>SAMN37516308</td>
      <td>SRR26148464</td>
      <td>Illumina MiSeq 2x250</td>
      <td>https://github.com/aomlomics/protocols/blob/main/protocol_DNA_extraction_Sterivex.md</td>
      <td>16S V4-V5</td>
      <td>16S rRNA</td>
      <td>V4-V5</td>
      <td>GTGYCAGCMGCCGCGGTAA</td>
      <td>CCGYCAATTYMTTTRAGTTT</td>
      <td>515F-Y</td>
      <td>926R</td>
      <td>10.1111/1462-2920.13023</td>
      <td>initial denaturation:95_2;denaturation:95_0.75;annealing:50_0.75;elongation:68_1.5;final elongation:68_5;25</td>
      <td>10.1111/1462-2920.13023</td>
      <td>ACACTGACGACATGGTTCTACA;TACGGTAGCAGAGACTTGGTCT</td>
      <td>missing: not provided</td>
    </tr>
  </tbody>
</table>
</div>




```python
prep_18 = dna_prep[dna_prep['amplicon_sequenced'].str.contains('18S V9')].copy()
prep_18['eventID'] = prep_18['eventID']+"_18S"
prep_18.head()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>eventID</th>
      <th>library_id</th>
      <th>title</th>
      <th>library_strategy</th>
      <th>library_source</th>
      <th>library_selection</th>
      <th>lib_layout</th>
      <th>platform</th>
      <th>instrument_model</th>
      <th>design_description</th>
      <th>filetype</th>
      <th>filename</th>
      <th>filename2</th>
      <th>biosample_accession</th>
      <th>sra_accession</th>
      <th>seq_meth</th>
      <th>nucl_acid_ext</th>
      <th>amplicon_sequenced</th>
      <th>target_gene</th>
      <th>target_subfragment</th>
      <th>pcr_primer_forward</th>
      <th>pcr_primer_reverse</th>
      <th>pcr_primer_name_forward</th>
      <th>pcr_primer_name_reverse</th>
      <th>pcr_primer_reference</th>
      <th>pcr_cond</th>
      <th>nucl_acid_amp</th>
      <th>adapters</th>
      <th>mid_barcode</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>1</th>
      <td>GOMECC4_27N_Sta1_DCM_A_18S</td>
      <td>GOMECC18S_Plate4_53</td>
      <td>18S amplicon metabarcoding of marine metagenome: Gulf of Mexico (USA)</td>
      <td>AMPLICON</td>
      <td>METAGENOMIC</td>
      <td>PCR</td>
      <td>paired</td>
      <td>ILLUMINA</td>
      <td>Illumina MiSeq</td>
      <td>Samples were collected and filtered onto Sterivex 0.22 um cartridge filters. DNA was extracted from Sterivex by adding lysis buffer and magnetic b...</td>
      <td>fastq</td>
      <td>GOMECC18S_Plate4_53_S340_L001_R1_001.fastq.gz</td>
      <td>GOMECC18S_Plate4_53_S340_L001_R2_001.fastq.gz</td>
      <td>SAMN37516094</td>
      <td>SRR26161153</td>
      <td>Illumina MiSeq 2x250</td>
      <td>https://github.com/aomlomics/protocols/blob/main/protocol_DNA_extraction_Sterivex.md</td>
      <td>18S V9</td>
      <td>18S rRNA</td>
      <td>V9</td>
      <td>GTACACACCGCCCGTC</td>
      <td>TGATCCTTCTGCAGGTTCACCTAC</td>
      <td>1391f</td>
      <td>EukBr</td>
      <td>10.1371/journal.pone.0006372</td>
      <td>initial denaturation:94_3;denaturation:94_0.75;annealing:65_0.25;57_0.5;elongation:72_1.5;final elongation:72_10;35</td>
      <td>10.1371/journal.pone.0006372</td>
      <td>ACACTGACGACATGGTTCTACA;TACGGTAGCAGAGACTTGGTCT</td>
      <td>missing: not provided</td>
    </tr>
    <tr>
      <th>3</th>
      <td>GOMECC4_27N_Sta1_DCM_B_18S</td>
      <td>GOMECC18S_Plate4_46</td>
      <td>18S amplicon metabarcoding of marine metagenome: Gulf of Mexico (USA)</td>
      <td>AMPLICON</td>
      <td>METAGENOMIC</td>
      <td>PCR</td>
      <td>paired</td>
      <td>ILLUMINA</td>
      <td>Illumina MiSeq</td>
      <td>Samples were collected and filtered onto Sterivex 0.22 um cartridge filters. DNA was extracted from Sterivex by adding lysis buffer and magnetic b...</td>
      <td>fastq</td>
      <td>GOMECC18S_Plate4_46_S333_L001_R1_001.fastq.gz</td>
      <td>GOMECC18S_Plate4_46_S333_L001_R2_001.fastq.gz</td>
      <td>SAMN37516095</td>
      <td>SRR26161138</td>
      <td>Illumina MiSeq 2x250</td>
      <td>https://github.com/aomlomics/protocols/blob/main/protocol_DNA_extraction_Sterivex.md</td>
      <td>18S V9</td>
      <td>18S rRNA</td>
      <td>V9</td>
      <td>GTACACACCGCCCGTC</td>
      <td>TGATCCTTCTGCAGGTTCACCTAC</td>
      <td>1391f</td>
      <td>EukBr</td>
      <td>10.1371/journal.pone.0006372</td>
      <td>initial denaturation:94_3;denaturation:94_0.75;annealing:65_0.25;57_0.5;elongation:72_1.5;final elongation:72_10;35</td>
      <td>10.1371/journal.pone.0006372</td>
      <td>ACACTGACGACATGGTTCTACA;TACGGTAGCAGAGACTTGGTCT</td>
      <td>missing: not provided</td>
    </tr>
    <tr>
      <th>5</th>
      <td>GOMECC4_27N_Sta1_DCM_C_18S</td>
      <td>GOMECC18S_Plate4_54</td>
      <td>18S amplicon metabarcoding of marine metagenome: Gulf of Mexico (USA)</td>
      <td>AMPLICON</td>
      <td>METAGENOMIC</td>
      <td>PCR</td>
      <td>paired</td>
      <td>ILLUMINA</td>
      <td>Illumina MiSeq</td>
      <td>Samples were collected and filtered onto Sterivex 0.22 um cartridge filters. DNA was extracted from Sterivex by adding lysis buffer and magnetic b...</td>
      <td>fastq</td>
      <td>GOMECC18S_Plate4_54_S341_L001_R1_001.fastq.gz</td>
      <td>GOMECC18S_Plate4_54_S341_L001_R2_001.fastq.gz</td>
      <td>SAMN37516096</td>
      <td>SRR26160919</td>
      <td>Illumina MiSeq 2x250</td>
      <td>https://github.com/aomlomics/protocols/blob/main/protocol_DNA_extraction_Sterivex.md</td>
      <td>18S V9</td>
      <td>18S rRNA</td>
      <td>V9</td>
      <td>GTACACACCGCCCGTC</td>
      <td>TGATCCTTCTGCAGGTTCACCTAC</td>
      <td>1391f</td>
      <td>EukBr</td>
      <td>10.1371/journal.pone.0006372</td>
      <td>initial denaturation:94_3;denaturation:94_0.75;annealing:65_0.25;57_0.5;elongation:72_1.5;final elongation:72_10;35</td>
      <td>10.1371/journal.pone.0006372</td>
      <td>ACACTGACGACATGGTTCTACA;TACGGTAGCAGAGACTTGGTCT</td>
      <td>missing: not provided</td>
    </tr>
    <tr>
      <th>7</th>
      <td>GOMECC4_27N_Sta1_Deep_A_18S</td>
      <td>GOMECC18S_Plate4_52</td>
      <td>18S amplicon metabarcoding of marine metagenome: Gulf of Mexico (USA)</td>
      <td>AMPLICON</td>
      <td>METAGENOMIC</td>
      <td>PCR</td>
      <td>paired</td>
      <td>ILLUMINA</td>
      <td>Illumina MiSeq</td>
      <td>Samples were collected and filtered onto Sterivex 0.22 um cartridge filters. DNA was extracted from Sterivex by adding lysis buffer and magnetic b...</td>
      <td>fastq</td>
      <td>GOMECC18S_Plate4_52_S339_L001_R1_001.fastq.gz</td>
      <td>GOMECC18S_Plate4_52_S339_L001_R2_001.fastq.gz</td>
      <td>SAMN37516091</td>
      <td>SRR26160709</td>
      <td>Illumina MiSeq 2x250</td>
      <td>https://github.com/aomlomics/protocols/blob/main/protocol_DNA_extraction_Sterivex.md</td>
      <td>18S V9</td>
      <td>18S rRNA</td>
      <td>V9</td>
      <td>GTACACACCGCCCGTC</td>
      <td>TGATCCTTCTGCAGGTTCACCTAC</td>
      <td>1391f</td>
      <td>EukBr</td>
      <td>10.1371/journal.pone.0006372</td>
      <td>initial denaturation:94_3;denaturation:94_0.75;annealing:65_0.25;57_0.5;elongation:72_1.5;final elongation:72_10;35</td>
      <td>10.1371/journal.pone.0006372</td>
      <td>ACACTGACGACATGGTTCTACA;TACGGTAGCAGAGACTTGGTCT</td>
      <td>missing: not provided</td>
    </tr>
    <tr>
      <th>9</th>
      <td>GOMECC4_27N_Sta1_Deep_B_18S</td>
      <td>GOMECC18S_Plate4_60</td>
      <td>18S amplicon metabarcoding of marine metagenome: Gulf of Mexico (USA)</td>
      <td>AMPLICON</td>
      <td>METAGENOMIC</td>
      <td>PCR</td>
      <td>paired</td>
      <td>ILLUMINA</td>
      <td>Illumina MiSeq</td>
      <td>Samples were collected and filtered onto Sterivex 0.22 um cartridge filters. DNA was extracted from Sterivex by adding lysis buffer and magnetic b...</td>
      <td>fastq</td>
      <td>GOMECC18S_Plate4_60_S347_L001_R1_001.fastq.gz</td>
      <td>GOMECC18S_Plate4_60_S347_L001_R2_001.fastq.gz</td>
      <td>SAMN37516092</td>
      <td>SRR26160970</td>
      <td>Illumina MiSeq 2x250</td>
      <td>https://github.com/aomlomics/protocols/blob/main/protocol_DNA_extraction_Sterivex.md</td>
      <td>18S V9</td>
      <td>18S rRNA</td>
      <td>V9</td>
      <td>GTACACACCGCCCGTC</td>
      <td>TGATCCTTCTGCAGGTTCACCTAC</td>
      <td>1391f</td>
      <td>EukBr</td>
      <td>10.1371/journal.pone.0006372</td>
      <td>initial denaturation:94_3;denaturation:94_0.75;annealing:65_0.25;57_0.5;elongation:72_1.5;final elongation:72_10;35</td>
      <td>10.1371/journal.pone.0006372</td>
      <td>ACACTGACGACATGGTTCTACA;TACGGTAGCAGAGACTTGGTCT</td>
      <td>missing: not provided</td>
    </tr>
  </tbody>
</table>
</div>




```python
dna_prep = pd.concat([prep_16,prep_18],axis=0,ignore_index=True)
dna_prep.head()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>eventID</th>
      <th>library_id</th>
      <th>title</th>
      <th>library_strategy</th>
      <th>library_source</th>
      <th>library_selection</th>
      <th>lib_layout</th>
      <th>platform</th>
      <th>instrument_model</th>
      <th>design_description</th>
      <th>filetype</th>
      <th>filename</th>
      <th>filename2</th>
      <th>biosample_accession</th>
      <th>sra_accession</th>
      <th>seq_meth</th>
      <th>nucl_acid_ext</th>
      <th>amplicon_sequenced</th>
      <th>target_gene</th>
      <th>target_subfragment</th>
      <th>pcr_primer_forward</th>
      <th>pcr_primer_reverse</th>
      <th>pcr_primer_name_forward</th>
      <th>pcr_primer_name_reverse</th>
      <th>pcr_primer_reference</th>
      <th>pcr_cond</th>
      <th>nucl_acid_amp</th>
      <th>adapters</th>
      <th>mid_barcode</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>GOMECC4_BROWNSVILLE_Sta66_DCM_B_16S</td>
      <td>GOMECC16S_Plate1_1</td>
      <td>16S amplicon metabarcoding of marine metagenome: Gulf of Mexico (USA)</td>
      <td>AMPLICON</td>
      <td>METAGENOMIC</td>
      <td>PCR</td>
      <td>paired</td>
      <td>ILLUMINA</td>
      <td>Illumina MiSeq</td>
      <td>Samples were collected and filtered onto Sterivex 0.22 um cartridge filters. DNA was extracted from Sterivex by adding lysis buffer and magnetic b...</td>
      <td>fastq</td>
      <td>GOMECC16S_Plate1_1_S1_L001_R1_001.fastq.gz</td>
      <td>GOMECC16S_Plate1_1_S1_L001_R2_001.fastq.gz</td>
      <td>SAMN37516307</td>
      <td>SRR26148474</td>
      <td>Illumina MiSeq 2x250</td>
      <td>https://github.com/aomlomics/protocols/blob/main/protocol_DNA_extraction_Sterivex.md</td>
      <td>16S V4-V5</td>
      <td>16S rRNA</td>
      <td>V4-V5</td>
      <td>GTGYCAGCMGCCGCGGTAA</td>
      <td>CCGYCAATTYMTTTRAGTTT</td>
      <td>515F-Y</td>
      <td>926R</td>
      <td>10.1111/1462-2920.13023</td>
      <td>initial denaturation:95_2;denaturation:95_0.75;annealing:50_0.75;elongation:68_1.5;final elongation:68_5;25</td>
      <td>10.1111/1462-2920.13023</td>
      <td>ACACTGACGACATGGTTCTACA;TACGGTAGCAGAGACTTGGTCT</td>
      <td>missing: not provided</td>
    </tr>
    <tr>
      <th>1</th>
      <td>GOMECC4_GALVESTON_Sta54_DCM_B_16S</td>
      <td>GOMECC16S_Plate1_10</td>
      <td>16S amplicon metabarcoding of marine metagenome: Gulf of Mexico (USA)</td>
      <td>AMPLICON</td>
      <td>METAGENOMIC</td>
      <td>PCR</td>
      <td>paired</td>
      <td>ILLUMINA</td>
      <td>Illumina MiSeq</td>
      <td>Samples were collected and filtered onto Sterivex 0.22 um cartridge filters. DNA was extracted from Sterivex by adding lysis buffer and magnetic b...</td>
      <td>fastq</td>
      <td>GOMECC16S_Plate1_10_S10_L001_R1_001.fastq.gz</td>
      <td>GOMECC16S_Plate1_10_S10_L001_R2_001.fastq.gz</td>
      <td>SAMN37516268</td>
      <td>SRR26148413</td>
      <td>Illumina MiSeq 2x250</td>
      <td>https://github.com/aomlomics/protocols/blob/main/protocol_DNA_extraction_Sterivex.md</td>
      <td>16S V4-V5</td>
      <td>16S rRNA</td>
      <td>V4-V5</td>
      <td>GTGYCAGCMGCCGCGGTAA</td>
      <td>CCGYCAATTYMTTTRAGTTT</td>
      <td>515F-Y</td>
      <td>926R</td>
      <td>10.1111/1462-2920.13023</td>
      <td>initial denaturation:95_2;denaturation:95_0.75;annealing:50_0.75;elongation:68_1.5;final elongation:68_5;25</td>
      <td>10.1111/1462-2920.13023</td>
      <td>ACACTGACGACATGGTTCTACA;TACGGTAGCAGAGACTTGGTCT</td>
      <td>missing: not provided</td>
    </tr>
    <tr>
      <th>2</th>
      <td>GOMECC4_GALVESTON_Sta54_Deep_A_16S</td>
      <td>GOMECC16S_Plate1_11</td>
      <td>16S amplicon metabarcoding of marine metagenome: Gulf of Mexico (USA)</td>
      <td>AMPLICON</td>
      <td>METAGENOMIC</td>
      <td>PCR</td>
      <td>paired</td>
      <td>ILLUMINA</td>
      <td>Illumina MiSeq</td>
      <td>Samples were collected and filtered onto Sterivex 0.22 um cartridge filters. DNA was extracted from Sterivex by adding lysis buffer and magnetic b...</td>
      <td>fastq</td>
      <td>GOMECC16S_Plate1_11_S11_L001_R1_001.fastq.gz</td>
      <td>GOMECC16S_Plate1_11_S11_L001_R2_001.fastq.gz</td>
      <td>SAMN37516264</td>
      <td>SRR26148140</td>
      <td>Illumina MiSeq 2x250</td>
      <td>https://github.com/aomlomics/protocols/blob/main/protocol_DNA_extraction_Sterivex.md</td>
      <td>16S V4-V5</td>
      <td>16S rRNA</td>
      <td>V4-V5</td>
      <td>GTGYCAGCMGCCGCGGTAA</td>
      <td>CCGYCAATTYMTTTRAGTTT</td>
      <td>515F-Y</td>
      <td>926R</td>
      <td>10.1111/1462-2920.13023</td>
      <td>initial denaturation:95_2;denaturation:95_0.75;annealing:50_0.75;elongation:68_1.5;final elongation:68_5;25</td>
      <td>10.1111/1462-2920.13023</td>
      <td>ACACTGACGACATGGTTCTACA;TACGGTAGCAGAGACTTGGTCT</td>
      <td>missing: not provided</td>
    </tr>
    <tr>
      <th>3</th>
      <td>GOMECC4_GALVESTON_Sta49_Deep_A_16S</td>
      <td>GOMECC16S_Plate1_12</td>
      <td>16S amplicon metabarcoding of marine metagenome: Gulf of Mexico (USA)</td>
      <td>AMPLICON</td>
      <td>METAGENOMIC</td>
      <td>PCR</td>
      <td>paired</td>
      <td>ILLUMINA</td>
      <td>Illumina MiSeq</td>
      <td>Samples were collected and filtered onto Sterivex 0.22 um cartridge filters. DNA was extracted from Sterivex by adding lysis buffer and magnetic b...</td>
      <td>fastq</td>
      <td>GOMECC16S_Plate1_12_S12_L001_R1_001.fastq.gz</td>
      <td>GOMECC16S_Plate1_12_S12_L001_R2_001.fastq.gz</td>
      <td>SAMN37516246</td>
      <td>SRR26148197</td>
      <td>Illumina MiSeq 2x250</td>
      <td>https://github.com/aomlomics/protocols/blob/main/protocol_DNA_extraction_Sterivex.md</td>
      <td>16S V4-V5</td>
      <td>16S rRNA</td>
      <td>V4-V5</td>
      <td>GTGYCAGCMGCCGCGGTAA</td>
      <td>CCGYCAATTYMTTTRAGTTT</td>
      <td>515F-Y</td>
      <td>926R</td>
      <td>10.1111/1462-2920.13023</td>
      <td>initial denaturation:95_2;denaturation:95_0.75;annealing:50_0.75;elongation:68_1.5;final elongation:68_5;25</td>
      <td>10.1111/1462-2920.13023</td>
      <td>ACACTGACGACATGGTTCTACA;TACGGTAGCAGAGACTTGGTCT</td>
      <td>missing: not provided</td>
    </tr>
    <tr>
      <th>4</th>
      <td>GOMECC4_BROWNSVILLE_Sta66_DCM_C_16S</td>
      <td>GOMECC16S_Plate1_13</td>
      <td>16S amplicon metabarcoding of marine metagenome: Gulf of Mexico (USA)</td>
      <td>AMPLICON</td>
      <td>METAGENOMIC</td>
      <td>PCR</td>
      <td>paired</td>
      <td>ILLUMINA</td>
      <td>Illumina MiSeq</td>
      <td>Samples were collected and filtered onto Sterivex 0.22 um cartridge filters. DNA was extracted from Sterivex by adding lysis buffer and magnetic b...</td>
      <td>fastq</td>
      <td>GOMECC16S_Plate1_13_S13_L001_R1_001.fastq.gz</td>
      <td>GOMECC16S_Plate1_13_S13_L001_R2_001.fastq.gz</td>
      <td>SAMN37516308</td>
      <td>SRR26148464</td>
      <td>Illumina MiSeq 2x250</td>
      <td>https://github.com/aomlomics/protocols/blob/main/protocol_DNA_extraction_Sterivex.md</td>
      <td>16S V4-V5</td>
      <td>16S rRNA</td>
      <td>V4-V5</td>
      <td>GTGYCAGCMGCCGCGGTAA</td>
      <td>CCGYCAATTYMTTTRAGTTT</td>
      <td>515F-Y</td>
      <td>926R</td>
      <td>10.1111/1462-2920.13023</td>
      <td>initial denaturation:95_2;denaturation:95_0.75;annealing:50_0.75;elongation:68_1.5;final elongation:68_5;25</td>
      <td>10.1111/1462-2920.13023</td>
      <td>ACACTGACGACATGGTTCTACA;TACGGTAGCAGAGACTTGGTCT</td>
      <td>missing: not provided</td>
    </tr>
  </tbody>
</table>
</div>




```python
# merge prep and sample
dna = dna_sample.merge(dna_prep, how='outer', on='eventID')
dna = dna.merge(dna_analysis,how='outer',on='amplicon_sequenced')
```


```python
rename_dict.values()
```




    dict_values(['samp_vol_we_dna_ext', 'samp_mat_process', 'env_broad_scale', 'env_local_scale', 'env_medium', 'size_frac', 'concentration', 'concentrationUnit', 'samp_collect_device', 'source_mat_id', 'eventID', 'nucl_acid_ext', 'nucl_acid_amp', 'target_gene', 'target_subfragment', 'lib_layout', 'pcr_primer_forward', 'pcr_primer_reverse', 'pcr_primer_name_forward', 'pcr_primer_name_reverse', 'pcr_primer_reference', 'pcr_cond', 'seq_meth', 'sop', 'ampliconSize', 'otu_class_appr', 'otu_seq_comp_appr', 'otu_db'])




```python
#which columns are not in the list of values for dna-derived extension?
[col for col in dna if col not in rename_dict.values()]
```




    ['serial_number',
     'cruise_id',
     'line_id',
     'station',
     'ctd_bottle_no',
     'sample_replicate',
     'biological_replicates',
     'extract_number',
     'sample_title',
     'bioproject_accession',
     'biosample_accession_x',
     'metagenome_sequenced',
     'organism',
     'collection_date_local',
     'collection_date',
     'depth',
     'geo_loc_name',
     'lat_lon',
     'decimalLatitude',
     'decimalLongitude',
     'sample_type',
     'collection_method',
     'basisOfRecord',
     'cluster_16s',
     'cluster_18s',
     'notes_sampling',
     'line_position',
     'offshore_inshore_200m_isobath',
     'depth_category',
     'ocean_acidification_status',
     'seascape_class',
     'seascape_probability',
     'seascape_window',
     'dna_sample_number',
     'dna_yield',
     'extraction_plate_name',
     'extraction_well_number',
     'extraction_well_position',
     'ship_crs_expocode',
     'woce_sect',
     'ammonium',
     'carbonate',
     'diss_inorg_carb',
     'diss_oxygen',
     'fluor',
     'hydrogen_ion',
     'nitrate',
     'nitrite',
     'nitrate_plus_nitrite',
     'omega_arag',
     'pco2',
     'ph',
     'phosphate',
     'pressure',
     'salinity',
     'samp_store_loc',
     'samp_store_temp',
     'silicate',
     'size_frac_low',
     'size_frac_up',
     'temp',
     'tot_alkalinity',
     'tot_depth_water_col',
     'transmittance',
     'waterBody',
     'library_id',
     'title',
     'library_strategy',
     'library_source',
     'library_selection',
     'platform',
     'instrument_model',
     'design_description',
     'filetype',
     'filename',
     'filename2',
     'biosample_accession_y',
     'sra_accession',
     'amplicon_sequenced',
     'adapters',
     'mid_barcode',
     'trim_method',
     'cluster_method',
     'pid_clustering',
     'code_repo',
     'identificationReferences',
     'controls_used']




```python
dna = dna.drop(columns=[col for col in dna if col not in rename_dict.values()])
```


```python
dna.tail()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>eventID</th>
      <th>source_mat_id</th>
      <th>env_broad_scale</th>
      <th>env_local_scale</th>
      <th>env_medium</th>
      <th>samp_vol_we_dna_ext</th>
      <th>samp_collect_device</th>
      <th>samp_mat_process</th>
      <th>size_frac</th>
      <th>concentration</th>
      <th>lib_layout</th>
      <th>seq_meth</th>
      <th>nucl_acid_ext</th>
      <th>target_gene</th>
      <th>target_subfragment</th>
      <th>pcr_primer_forward</th>
      <th>pcr_primer_reverse</th>
      <th>pcr_primer_name_forward</th>
      <th>pcr_primer_name_reverse</th>
      <th>pcr_primer_reference</th>
      <th>pcr_cond</th>
      <th>nucl_acid_amp</th>
      <th>ampliconSize</th>
      <th>otu_seq_comp_appr</th>
      <th>otu_db</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>939</th>
      <td>GOMECC4_CAPECORAL_Sta141_DCM_B_18S</td>
      <td>GOMECC4_CAPECORAL_Sta141_DCM</td>
      <td>marine biome [ENVO:00000447]</td>
      <td>marine photic zone [ENVO:00000209]</td>
      <td>sea water [ENVO:00002149]</td>
      <td>2040 ml</td>
      <td>Niskin bottle</td>
      <td>Pumped through Sterivex filter (0.22-µm) using peristaltic pµmp</td>
      <td>0.22 µm</td>
      <td>1.634 ng/µl</td>
      <td>paired</td>
      <td>Illumina MiSeq 2x250</td>
      <td>https://github.com/aomlomics/protocols/blob/main/protocol_DNA_extraction_Sterivex.md</td>
      <td>18S rRNA</td>
      <td>V9</td>
      <td>GTACACACCGCCCGTC</td>
      <td>TGATCCTTCTGCAGGTTCACCTAC</td>
      <td>1391f</td>
      <td>EukBr</td>
      <td>10.1371/journal.pone.0006372</td>
      <td>initial denaturation:94_3;denaturation:94_0.75;annealing:65_0.25;57_0.5;elongation:72_1.5;final elongation:72_10;35</td>
      <td>10.1371/journal.pone.0006372</td>
      <td>260</td>
      <td>Tourmaline; qiime2-2021.2; naive-bayes classifier</td>
      <td>PR2 v5.0.1; V9 1391f-1510r region; 10.5281/zenodo.8392706</td>
    </tr>
    <tr>
      <th>940</th>
      <td>GOMECC4_CAPECORAL_Sta141_DCM_C_18S</td>
      <td>GOMECC4_CAPECORAL_Sta141_DCM</td>
      <td>marine biome [ENVO:00000447]</td>
      <td>marine photic zone [ENVO:00000209]</td>
      <td>sea water [ENVO:00002149]</td>
      <td>2080 ml</td>
      <td>Niskin bottle</td>
      <td>Pumped through Sterivex filter (0.22-µm) using peristaltic pµmp</td>
      <td>0.22 µm</td>
      <td>2.307 ng/µl</td>
      <td>paired</td>
      <td>Illumina MiSeq 2x250</td>
      <td>https://github.com/aomlomics/protocols/blob/main/protocol_DNA_extraction_Sterivex.md</td>
      <td>18S rRNA</td>
      <td>V9</td>
      <td>GTACACACCGCCCGTC</td>
      <td>TGATCCTTCTGCAGGTTCACCTAC</td>
      <td>1391f</td>
      <td>EukBr</td>
      <td>10.1371/journal.pone.0006372</td>
      <td>initial denaturation:94_3;denaturation:94_0.75;annealing:65_0.25;57_0.5;elongation:72_1.5;final elongation:72_10;35</td>
      <td>10.1371/journal.pone.0006372</td>
      <td>260</td>
      <td>Tourmaline; qiime2-2021.2; naive-bayes classifier</td>
      <td>PR2 v5.0.1; V9 1391f-1510r region; 10.5281/zenodo.8392706</td>
    </tr>
    <tr>
      <th>941</th>
      <td>GOMECC4_CAPECORAL_Sta141_Surface_A_18S</td>
      <td>GOMECC4_CAPECORAL_Sta141_Surface</td>
      <td>marine biome [ENVO:00000447]</td>
      <td>marine photic zone [ENVO:00000209]</td>
      <td>sea water [ENVO:00002149]</td>
      <td>2100 ml</td>
      <td>Niskin bottle</td>
      <td>Pumped through Sterivex filter (0.22-µm) using peristaltic pµmp</td>
      <td>0.22 µm</td>
      <td>1.286 ng/µl</td>
      <td>paired</td>
      <td>Illumina MiSeq 2x250</td>
      <td>https://github.com/aomlomics/protocols/blob/main/protocol_DNA_extraction_Sterivex.md</td>
      <td>18S rRNA</td>
      <td>V9</td>
      <td>GTACACACCGCCCGTC</td>
      <td>TGATCCTTCTGCAGGTTCACCTAC</td>
      <td>1391f</td>
      <td>EukBr</td>
      <td>10.1371/journal.pone.0006372</td>
      <td>initial denaturation:94_3;denaturation:94_0.75;annealing:65_0.25;57_0.5;elongation:72_1.5;final elongation:72_10;35</td>
      <td>10.1371/journal.pone.0006372</td>
      <td>260</td>
      <td>Tourmaline; qiime2-2021.2; naive-bayes classifier</td>
      <td>PR2 v5.0.1; V9 1391f-1510r region; 10.5281/zenodo.8392706</td>
    </tr>
    <tr>
      <th>942</th>
      <td>GOMECC4_CAPECORAL_Sta141_Surface_B_18S</td>
      <td>GOMECC4_CAPECORAL_Sta141_Surface</td>
      <td>marine biome [ENVO:00000447]</td>
      <td>marine photic zone [ENVO:00000209]</td>
      <td>sea water [ENVO:00002149]</td>
      <td>2000 ml</td>
      <td>Niskin bottle</td>
      <td>Pumped through Sterivex filter (0.22-µm) using peristaltic pµmp</td>
      <td>0.22 µm</td>
      <td>1.831 ng/µl</td>
      <td>paired</td>
      <td>Illumina MiSeq 2x250</td>
      <td>https://github.com/aomlomics/protocols/blob/main/protocol_DNA_extraction_Sterivex.md</td>
      <td>18S rRNA</td>
      <td>V9</td>
      <td>GTACACACCGCCCGTC</td>
      <td>TGATCCTTCTGCAGGTTCACCTAC</td>
      <td>1391f</td>
      <td>EukBr</td>
      <td>10.1371/journal.pone.0006372</td>
      <td>initial denaturation:94_3;denaturation:94_0.75;annealing:65_0.25;57_0.5;elongation:72_1.5;final elongation:72_10;35</td>
      <td>10.1371/journal.pone.0006372</td>
      <td>260</td>
      <td>Tourmaline; qiime2-2021.2; naive-bayes classifier</td>
      <td>PR2 v5.0.1; V9 1391f-1510r region; 10.5281/zenodo.8392706</td>
    </tr>
    <tr>
      <th>943</th>
      <td>GOMECC4_CAPECORAL_Sta141_Surface_C_18S</td>
      <td>GOMECC4_CAPECORAL_Sta141_Surface</td>
      <td>marine biome [ENVO:00000447]</td>
      <td>marine photic zone [ENVO:00000209]</td>
      <td>sea water [ENVO:00002149]</td>
      <td>2000 ml</td>
      <td>Niskin bottle</td>
      <td>Pumped through Sterivex filter (0.22-µm) using peristaltic pµmp</td>
      <td>0.22 µm</td>
      <td>1.849 ng/µl</td>
      <td>paired</td>
      <td>Illumina MiSeq 2x250</td>
      <td>https://github.com/aomlomics/protocols/blob/main/protocol_DNA_extraction_Sterivex.md</td>
      <td>18S rRNA</td>
      <td>V9</td>
      <td>GTACACACCGCCCGTC</td>
      <td>TGATCCTTCTGCAGGTTCACCTAC</td>
      <td>1391f</td>
      <td>EukBr</td>
      <td>10.1371/journal.pone.0006372</td>
      <td>initial denaturation:94_3;denaturation:94_0.75;annealing:65_0.25;57_0.5;elongation:72_1.5;final elongation:72_10;35</td>
      <td>10.1371/journal.pone.0006372</td>
      <td>260</td>
      <td>Tourmaline; qiime2-2021.2; naive-bayes classifier</td>
      <td>PR2 v5.0.1; V9 1391f-1510r region; 10.5281/zenodo.8392706</td>
    </tr>
  </tbody>
</table>
</div>



#### merge with occurrenceID, DNA_sequence


```python
dna.shape
```




    (944, 25)




```python
dna_occ = dna.merge(occ_all[['eventID','occurrenceID','DNA_sequence']],how='left',on='eventID')
```


```python
dna_occ.shape
```




    (311390, 27)




```python
dna_occ.head()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>eventID</th>
      <th>source_mat_id</th>
      <th>env_broad_scale</th>
      <th>env_local_scale</th>
      <th>env_medium</th>
      <th>samp_vol_we_dna_ext</th>
      <th>samp_collect_device</th>
      <th>samp_mat_process</th>
      <th>size_frac</th>
      <th>concentration</th>
      <th>lib_layout</th>
      <th>seq_meth</th>
      <th>nucl_acid_ext</th>
      <th>target_gene</th>
      <th>target_subfragment</th>
      <th>pcr_primer_forward</th>
      <th>pcr_primer_reverse</th>
      <th>pcr_primer_name_forward</th>
      <th>pcr_primer_name_reverse</th>
      <th>pcr_primer_reference</th>
      <th>pcr_cond</th>
      <th>nucl_acid_amp</th>
      <th>ampliconSize</th>
      <th>otu_seq_comp_appr</th>
      <th>otu_db</th>
      <th>occurrenceID</th>
      <th>DNA_sequence</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>GOMECC4_27N_Sta1_Deep_A_16S</td>
      <td>GOMECC4_27N_Sta1_Deep</td>
      <td>marine biome [ENVO:00000447]</td>
      <td>marine mesopelagic zone [ENVO:00000213]</td>
      <td>sea water [ENVO:00002149]</td>
      <td>1920 ml</td>
      <td>Niskin bottle</td>
      <td>Pumped through Sterivex filter (0.22-µm) using peristaltic pµmp</td>
      <td>0.22 µm</td>
      <td>0.08038 ng/µl</td>
      <td>paired</td>
      <td>Illumina MiSeq 2x250</td>
      <td>https://github.com/aomlomics/protocols/blob/main/protocol_DNA_extraction_Sterivex.md</td>
      <td>16S rRNA</td>
      <td>V4-V5</td>
      <td>GTGYCAGCMGCCGCGGTAA</td>
      <td>CCGYCAATTYMTTTRAGTTT</td>
      <td>515F-Y</td>
      <td>926R</td>
      <td>10.1111/1462-2920.13023</td>
      <td>initial denaturation:95_2;denaturation:95_0.75;annealing:50_0.75;elongation:68_1.5;final elongation:68_5;25</td>
      <td>10.1111/1462-2920.13023</td>
      <td>411</td>
      <td>Tourmaline; qiime2-2021.2; naive-bayes classifier</td>
      <td>Silva SSU Ref NR 99 v138.1; 515f-926r region; 10.5281/zenodo.8392695</td>
      <td>GOMECC4_27N_Sta1_Deep_A_16S_occ009257b156ab4a9dd2f0b0dd33100b7e</td>
      <td>TACGAGGGGTGCTAGCGTTGTCCGGAATTACTGGGCGTAAAGGGTTCGTAGGCGTCTTGCCAAGTTGATCGTTAAAGCCACCGGCTTAACCGGTGATCTGCGATCAAAACTGGCGAGATAGAATATGTGAGGGGAATGTGGAATTC...</td>
    </tr>
    <tr>
      <th>1</th>
      <td>GOMECC4_27N_Sta1_Deep_A_16S</td>
      <td>GOMECC4_27N_Sta1_Deep</td>
      <td>marine biome [ENVO:00000447]</td>
      <td>marine mesopelagic zone [ENVO:00000213]</td>
      <td>sea water [ENVO:00002149]</td>
      <td>1920 ml</td>
      <td>Niskin bottle</td>
      <td>Pumped through Sterivex filter (0.22-µm) using peristaltic pµmp</td>
      <td>0.22 µm</td>
      <td>0.08038 ng/µl</td>
      <td>paired</td>
      <td>Illumina MiSeq 2x250</td>
      <td>https://github.com/aomlomics/protocols/blob/main/protocol_DNA_extraction_Sterivex.md</td>
      <td>16S rRNA</td>
      <td>V4-V5</td>
      <td>GTGYCAGCMGCCGCGGTAA</td>
      <td>CCGYCAATTYMTTTRAGTTT</td>
      <td>515F-Y</td>
      <td>926R</td>
      <td>10.1111/1462-2920.13023</td>
      <td>initial denaturation:95_2;denaturation:95_0.75;annealing:50_0.75;elongation:68_1.5;final elongation:68_5;25</td>
      <td>10.1111/1462-2920.13023</td>
      <td>411</td>
      <td>Tourmaline; qiime2-2021.2; naive-bayes classifier</td>
      <td>Silva SSU Ref NR 99 v138.1; 515f-926r region; 10.5281/zenodo.8392695</td>
      <td>GOMECC4_27N_Sta1_Deep_A_16S_occ01398067b1d323b7f992a6764fa69e97</td>
      <td>TACGGAGGGTGCAAGCGTTGTTCGGAATTATTGGGCGTAAAGCGGATGTAGGCGGTCTGTCAAGTCGGATGTGAAATCCCTGGGCTCAACCCAGGAACTGCATTCGAAACTGTCAGACTAGAGTCTCGGAGGGGGTGGCGGAATTC...</td>
    </tr>
    <tr>
      <th>2</th>
      <td>GOMECC4_27N_Sta1_Deep_A_16S</td>
      <td>GOMECC4_27N_Sta1_Deep</td>
      <td>marine biome [ENVO:00000447]</td>
      <td>marine mesopelagic zone [ENVO:00000213]</td>
      <td>sea water [ENVO:00002149]</td>
      <td>1920 ml</td>
      <td>Niskin bottle</td>
      <td>Pumped through Sterivex filter (0.22-µm) using peristaltic pµmp</td>
      <td>0.22 µm</td>
      <td>0.08038 ng/µl</td>
      <td>paired</td>
      <td>Illumina MiSeq 2x250</td>
      <td>https://github.com/aomlomics/protocols/blob/main/protocol_DNA_extraction_Sterivex.md</td>
      <td>16S rRNA</td>
      <td>V4-V5</td>
      <td>GTGYCAGCMGCCGCGGTAA</td>
      <td>CCGYCAATTYMTTTRAGTTT</td>
      <td>515F-Y</td>
      <td>926R</td>
      <td>10.1111/1462-2920.13023</td>
      <td>initial denaturation:95_2;denaturation:95_0.75;annealing:50_0.75;elongation:68_1.5;final elongation:68_5;25</td>
      <td>10.1111/1462-2920.13023</td>
      <td>411</td>
      <td>Tourmaline; qiime2-2021.2; naive-bayes classifier</td>
      <td>Silva SSU Ref NR 99 v138.1; 515f-926r region; 10.5281/zenodo.8392695</td>
      <td>GOMECC4_27N_Sta1_Deep_A_16S_occ01770ea2fb7f041c787e5a481888c27e</td>
      <td>TACGGAGGATCCAAGCGTTATCCGGATTTATTGGGTTTAAAGGGTCCGCAGGCGGACTATTAAGTCAGTGGTGAAAGTCTGCAGCTTAACTGTAGAATTGCCATTGAAACTGATAGTCTTGAGTGTGGTTGAAGTGGGCGGAATAT...</td>
    </tr>
    <tr>
      <th>3</th>
      <td>GOMECC4_27N_Sta1_Deep_A_16S</td>
      <td>GOMECC4_27N_Sta1_Deep</td>
      <td>marine biome [ENVO:00000447]</td>
      <td>marine mesopelagic zone [ENVO:00000213]</td>
      <td>sea water [ENVO:00002149]</td>
      <td>1920 ml</td>
      <td>Niskin bottle</td>
      <td>Pumped through Sterivex filter (0.22-µm) using peristaltic pµmp</td>
      <td>0.22 µm</td>
      <td>0.08038 ng/µl</td>
      <td>paired</td>
      <td>Illumina MiSeq 2x250</td>
      <td>https://github.com/aomlomics/protocols/blob/main/protocol_DNA_extraction_Sterivex.md</td>
      <td>16S rRNA</td>
      <td>V4-V5</td>
      <td>GTGYCAGCMGCCGCGGTAA</td>
      <td>CCGYCAATTYMTTTRAGTTT</td>
      <td>515F-Y</td>
      <td>926R</td>
      <td>10.1111/1462-2920.13023</td>
      <td>initial denaturation:95_2;denaturation:95_0.75;annealing:50_0.75;elongation:68_1.5;final elongation:68_5;25</td>
      <td>10.1111/1462-2920.13023</td>
      <td>411</td>
      <td>Tourmaline; qiime2-2021.2; naive-bayes classifier</td>
      <td>Silva SSU Ref NR 99 v138.1; 515f-926r region; 10.5281/zenodo.8392695</td>
      <td>GOMECC4_27N_Sta1_Deep_A_16S_occ017dbdc8b62705bdf3f93218ac93a030</td>
      <td>TACTAGGGGTGCAAGCGTTGTCCGGAATTACTGGGCGTAAAGGGTGCGTAGGCGTCTACGTAAGTTGTTTGTTAAATCCATCGGCTTAACCGATGATCTGCAAACAAAACTGCATAGATAGAGTTTGGAAGAGGAAAGTGGAATTC...</td>
    </tr>
    <tr>
      <th>4</th>
      <td>GOMECC4_27N_Sta1_Deep_A_16S</td>
      <td>GOMECC4_27N_Sta1_Deep</td>
      <td>marine biome [ENVO:00000447]</td>
      <td>marine mesopelagic zone [ENVO:00000213]</td>
      <td>sea water [ENVO:00002149]</td>
      <td>1920 ml</td>
      <td>Niskin bottle</td>
      <td>Pumped through Sterivex filter (0.22-µm) using peristaltic pµmp</td>
      <td>0.22 µm</td>
      <td>0.08038 ng/µl</td>
      <td>paired</td>
      <td>Illumina MiSeq 2x250</td>
      <td>https://github.com/aomlomics/protocols/blob/main/protocol_DNA_extraction_Sterivex.md</td>
      <td>16S rRNA</td>
      <td>V4-V5</td>
      <td>GTGYCAGCMGCCGCGGTAA</td>
      <td>CCGYCAATTYMTTTRAGTTT</td>
      <td>515F-Y</td>
      <td>926R</td>
      <td>10.1111/1462-2920.13023</td>
      <td>initial denaturation:95_2;denaturation:95_0.75;annealing:50_0.75;elongation:68_1.5;final elongation:68_5;25</td>
      <td>10.1111/1462-2920.13023</td>
      <td>411</td>
      <td>Tourmaline; qiime2-2021.2; naive-bayes classifier</td>
      <td>Silva SSU Ref NR 99 v138.1; 515f-926r region; 10.5281/zenodo.8392695</td>
      <td>GOMECC4_27N_Sta1_Deep_A_16S_occ069f375524db7812103fe73fdefb7d2b</td>
      <td>TACGTAGGAGGCTAGCGTTGTCCGGATTTACTGGGCGTAAAGGGAGCGCAGGTGGCTGAGTTCGTCCGTGGTGCAAGCTCCAGGCCTAACCTGGAGAGGTCTACGGATACTGCTCGGCTTGAGGGCGGTAGAGGAGCACGGAATTC...</td>
    </tr>
  </tbody>
</table>
</div>




```python
dna_occ['concentration'] = dna_occ['concentration'].str.strip(" ng/µl")
dna_occ['concentrationUnit'] = "ng/µl"
```


```python
# check if all DwC terms are in dna file
for key in dna_dict.keys():
    if key not in dna_occ.columns:
        print(key,dna_dict[key])
```

    sop {'AOML_term': 'sop', 'AOML_file': 'analysis_data', 'DwC_definition': 'Standard operating procedures used in assembly and/or annotation of genomes, metagenomes or environmental sequences. Or A reference to a well documented protocol, e.g. using protocols.io', 'Example': nan}
    otu_class_appr {'AOML_term': 'derived: cluster_method, pid_clustering', 'AOML_file': 'analysis_data', 'DwC_definition': 'Approach/algorithm when defining OTUs or ASVs, include version and parameters separated by semicolons', 'Example': '"dada2; 1.14.0; ASV"'}



```python
data['analysis_data']['cluster_method'][0]
```




    'Tourmaline; qiime2-2021.2; dada2'




```python
dna_occ['seq_meth'] = 'Illumina MiSeq 2x250'
dna_occ['otu_class_appr']= data['analysis_data']['cluster_method'][0]+"; "+data['analysis_data']['pid_clustering'][0]
```


```python
# check if all DwC terms are in dna file
for key in dna_dict.keys():
    if key not in dna_occ.columns:
        print(key,dna_dict[key])
```

    sop {'AOML_term': 'sop', 'AOML_file': 'analysis_data', 'DwC_definition': 'Standard operating procedures used in assembly and/or annotation of genomes, metagenomes or environmental sequences. Or A reference to a well documented protocol, e.g. using protocols.io', 'Example': nan}



```python
dna_occ.head()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>eventID</th>
      <th>source_mat_id</th>
      <th>env_broad_scale</th>
      <th>env_local_scale</th>
      <th>env_medium</th>
      <th>samp_vol_we_dna_ext</th>
      <th>samp_collect_device</th>
      <th>samp_mat_process</th>
      <th>size_frac</th>
      <th>concentration</th>
      <th>lib_layout</th>
      <th>seq_meth</th>
      <th>nucl_acid_ext</th>
      <th>target_gene</th>
      <th>target_subfragment</th>
      <th>pcr_primer_forward</th>
      <th>pcr_primer_reverse</th>
      <th>pcr_primer_name_forward</th>
      <th>pcr_primer_name_reverse</th>
      <th>pcr_primer_reference</th>
      <th>pcr_cond</th>
      <th>nucl_acid_amp</th>
      <th>ampliconSize</th>
      <th>otu_seq_comp_appr</th>
      <th>otu_db</th>
      <th>occurrenceID</th>
      <th>DNA_sequence</th>
      <th>concentrationUnit</th>
      <th>otu_class_appr</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>GOMECC4_27N_Sta1_Deep_A_16S</td>
      <td>GOMECC4_27N_Sta1_Deep</td>
      <td>marine biome [ENVO:00000447]</td>
      <td>marine mesopelagic zone [ENVO:00000213]</td>
      <td>sea water [ENVO:00002149]</td>
      <td>1920 ml</td>
      <td>Niskin bottle</td>
      <td>Pumped through Sterivex filter (0.22-µm) using peristaltic pµmp</td>
      <td>0.22 µm</td>
      <td>0.08038</td>
      <td>paired</td>
      <td>Illumina MiSeq 2x250</td>
      <td>https://github.com/aomlomics/protocols/blob/main/protocol_DNA_extraction_Sterivex.md</td>
      <td>16S rRNA</td>
      <td>V4-V5</td>
      <td>GTGYCAGCMGCCGCGGTAA</td>
      <td>CCGYCAATTYMTTTRAGTTT</td>
      <td>515F-Y</td>
      <td>926R</td>
      <td>10.1111/1462-2920.13023</td>
      <td>initial denaturation:95_2;denaturation:95_0.75;annealing:50_0.75;elongation:68_1.5;final elongation:68_5;25</td>
      <td>10.1111/1462-2920.13023</td>
      <td>411</td>
      <td>Tourmaline; qiime2-2021.2; naive-bayes classifier</td>
      <td>Silva SSU Ref NR 99 v138.1; 515f-926r region; 10.5281/zenodo.8392695</td>
      <td>GOMECC4_27N_Sta1_Deep_A_16S_occ009257b156ab4a9dd2f0b0dd33100b7e</td>
      <td>TACGAGGGGTGCTAGCGTTGTCCGGAATTACTGGGCGTAAAGGGTTCGTAGGCGTCTTGCCAAGTTGATCGTTAAAGCCACCGGCTTAACCGGTGATCTGCGATCAAAACTGGCGAGATAGAATATGTGAGGGGAATGTGGAATTC...</td>
      <td>ng/µl</td>
      <td>Tourmaline; qiime2-2021.2; dada2; ASV</td>
    </tr>
    <tr>
      <th>1</th>
      <td>GOMECC4_27N_Sta1_Deep_A_16S</td>
      <td>GOMECC4_27N_Sta1_Deep</td>
      <td>marine biome [ENVO:00000447]</td>
      <td>marine mesopelagic zone [ENVO:00000213]</td>
      <td>sea water [ENVO:00002149]</td>
      <td>1920 ml</td>
      <td>Niskin bottle</td>
      <td>Pumped through Sterivex filter (0.22-µm) using peristaltic pµmp</td>
      <td>0.22 µm</td>
      <td>0.08038</td>
      <td>paired</td>
      <td>Illumina MiSeq 2x250</td>
      <td>https://github.com/aomlomics/protocols/blob/main/protocol_DNA_extraction_Sterivex.md</td>
      <td>16S rRNA</td>
      <td>V4-V5</td>
      <td>GTGYCAGCMGCCGCGGTAA</td>
      <td>CCGYCAATTYMTTTRAGTTT</td>
      <td>515F-Y</td>
      <td>926R</td>
      <td>10.1111/1462-2920.13023</td>
      <td>initial denaturation:95_2;denaturation:95_0.75;annealing:50_0.75;elongation:68_1.5;final elongation:68_5;25</td>
      <td>10.1111/1462-2920.13023</td>
      <td>411</td>
      <td>Tourmaline; qiime2-2021.2; naive-bayes classifier</td>
      <td>Silva SSU Ref NR 99 v138.1; 515f-926r region; 10.5281/zenodo.8392695</td>
      <td>GOMECC4_27N_Sta1_Deep_A_16S_occ01398067b1d323b7f992a6764fa69e97</td>
      <td>TACGGAGGGTGCAAGCGTTGTTCGGAATTATTGGGCGTAAAGCGGATGTAGGCGGTCTGTCAAGTCGGATGTGAAATCCCTGGGCTCAACCCAGGAACTGCATTCGAAACTGTCAGACTAGAGTCTCGGAGGGGGTGGCGGAATTC...</td>
      <td>ng/µl</td>
      <td>Tourmaline; qiime2-2021.2; dada2; ASV</td>
    </tr>
    <tr>
      <th>2</th>
      <td>GOMECC4_27N_Sta1_Deep_A_16S</td>
      <td>GOMECC4_27N_Sta1_Deep</td>
      <td>marine biome [ENVO:00000447]</td>
      <td>marine mesopelagic zone [ENVO:00000213]</td>
      <td>sea water [ENVO:00002149]</td>
      <td>1920 ml</td>
      <td>Niskin bottle</td>
      <td>Pumped through Sterivex filter (0.22-µm) using peristaltic pµmp</td>
      <td>0.22 µm</td>
      <td>0.08038</td>
      <td>paired</td>
      <td>Illumina MiSeq 2x250</td>
      <td>https://github.com/aomlomics/protocols/blob/main/protocol_DNA_extraction_Sterivex.md</td>
      <td>16S rRNA</td>
      <td>V4-V5</td>
      <td>GTGYCAGCMGCCGCGGTAA</td>
      <td>CCGYCAATTYMTTTRAGTTT</td>
      <td>515F-Y</td>
      <td>926R</td>
      <td>10.1111/1462-2920.13023</td>
      <td>initial denaturation:95_2;denaturation:95_0.75;annealing:50_0.75;elongation:68_1.5;final elongation:68_5;25</td>
      <td>10.1111/1462-2920.13023</td>
      <td>411</td>
      <td>Tourmaline; qiime2-2021.2; naive-bayes classifier</td>
      <td>Silva SSU Ref NR 99 v138.1; 515f-926r region; 10.5281/zenodo.8392695</td>
      <td>GOMECC4_27N_Sta1_Deep_A_16S_occ01770ea2fb7f041c787e5a481888c27e</td>
      <td>TACGGAGGATCCAAGCGTTATCCGGATTTATTGGGTTTAAAGGGTCCGCAGGCGGACTATTAAGTCAGTGGTGAAAGTCTGCAGCTTAACTGTAGAATTGCCATTGAAACTGATAGTCTTGAGTGTGGTTGAAGTGGGCGGAATAT...</td>
      <td>ng/µl</td>
      <td>Tourmaline; qiime2-2021.2; dada2; ASV</td>
    </tr>
    <tr>
      <th>3</th>
      <td>GOMECC4_27N_Sta1_Deep_A_16S</td>
      <td>GOMECC4_27N_Sta1_Deep</td>
      <td>marine biome [ENVO:00000447]</td>
      <td>marine mesopelagic zone [ENVO:00000213]</td>
      <td>sea water [ENVO:00002149]</td>
      <td>1920 ml</td>
      <td>Niskin bottle</td>
      <td>Pumped through Sterivex filter (0.22-µm) using peristaltic pµmp</td>
      <td>0.22 µm</td>
      <td>0.08038</td>
      <td>paired</td>
      <td>Illumina MiSeq 2x250</td>
      <td>https://github.com/aomlomics/protocols/blob/main/protocol_DNA_extraction_Sterivex.md</td>
      <td>16S rRNA</td>
      <td>V4-V5</td>
      <td>GTGYCAGCMGCCGCGGTAA</td>
      <td>CCGYCAATTYMTTTRAGTTT</td>
      <td>515F-Y</td>
      <td>926R</td>
      <td>10.1111/1462-2920.13023</td>
      <td>initial denaturation:95_2;denaturation:95_0.75;annealing:50_0.75;elongation:68_1.5;final elongation:68_5;25</td>
      <td>10.1111/1462-2920.13023</td>
      <td>411</td>
      <td>Tourmaline; qiime2-2021.2; naive-bayes classifier</td>
      <td>Silva SSU Ref NR 99 v138.1; 515f-926r region; 10.5281/zenodo.8392695</td>
      <td>GOMECC4_27N_Sta1_Deep_A_16S_occ017dbdc8b62705bdf3f93218ac93a030</td>
      <td>TACTAGGGGTGCAAGCGTTGTCCGGAATTACTGGGCGTAAAGGGTGCGTAGGCGTCTACGTAAGTTGTTTGTTAAATCCATCGGCTTAACCGATGATCTGCAAACAAAACTGCATAGATAGAGTTTGGAAGAGGAAAGTGGAATTC...</td>
      <td>ng/µl</td>
      <td>Tourmaline; qiime2-2021.2; dada2; ASV</td>
    </tr>
    <tr>
      <th>4</th>
      <td>GOMECC4_27N_Sta1_Deep_A_16S</td>
      <td>GOMECC4_27N_Sta1_Deep</td>
      <td>marine biome [ENVO:00000447]</td>
      <td>marine mesopelagic zone [ENVO:00000213]</td>
      <td>sea water [ENVO:00002149]</td>
      <td>1920 ml</td>
      <td>Niskin bottle</td>
      <td>Pumped through Sterivex filter (0.22-µm) using peristaltic pµmp</td>
      <td>0.22 µm</td>
      <td>0.08038</td>
      <td>paired</td>
      <td>Illumina MiSeq 2x250</td>
      <td>https://github.com/aomlomics/protocols/blob/main/protocol_DNA_extraction_Sterivex.md</td>
      <td>16S rRNA</td>
      <td>V4-V5</td>
      <td>GTGYCAGCMGCCGCGGTAA</td>
      <td>CCGYCAATTYMTTTRAGTTT</td>
      <td>515F-Y</td>
      <td>926R</td>
      <td>10.1111/1462-2920.13023</td>
      <td>initial denaturation:95_2;denaturation:95_0.75;annealing:50_0.75;elongation:68_1.5;final elongation:68_5;25</td>
      <td>10.1111/1462-2920.13023</td>
      <td>411</td>
      <td>Tourmaline; qiime2-2021.2; naive-bayes classifier</td>
      <td>Silva SSU Ref NR 99 v138.1; 515f-926r region; 10.5281/zenodo.8392695</td>
      <td>GOMECC4_27N_Sta1_Deep_A_16S_occ069f375524db7812103fe73fdefb7d2b</td>
      <td>TACGTAGGAGGCTAGCGTTGTCCGGATTTACTGGGCGTAAAGGGAGCGCAGGTGGCTGAGTTCGTCCGTGGTGCAAGCTCCAGGCCTAACCTGGAGAGGTCTACGGATACTGCTCGGCTTGAGGGCGGTAGAGGAGCACGGAATTC...</td>
      <td>ng/µl</td>
      <td>Tourmaline; qiime2-2021.2; dada2; ASV</td>
    </tr>
  </tbody>
</table>
</div>




```python
dna_occ.to_csv("../processed/dna-derived.csv",index=False)
```
