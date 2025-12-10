# Marine Xenobiotics Database (In Development)

## Introduction

MXD (provisory name) is a database to store potential genes for xenobiotics biorremediation, designed for easy visualization and usage for all marine microbiologists, regardless of knowledge in informatics and database handling. 

All the data presented in MXD originated from data mining of previous marine metagenomics studies, aiming a worldwide access to taxonomical and function of marine microbial communities. All the DNA sequences retrieved from the studies were submitted to a curated pipeline containing filters of quality control, metagenomics binning, taxonomic accessing and functional annotation. The processed data is then divided into a relational database, containing all the information of the studies, oceanographic and climatic variables and genomic functions and clusters. 

## Code and directory structure

The scripts are all located in a Scripts folder, with daughter directories related to specific parts of the projects, each one of these parts are run by a sub main script, which are called by the general main script, the "run.sh". 

-MXD/
  ├── Scripts/
  │   ├── Subproject1/
  │   │   ├── script1
  │   │   ├── script2
  │   │   ├── script3
  │   │   └── runsub1.sh
  │   ├── Subproject2/
  │   │   ├── script1
  │   │   ├── script2
  │   │   ├── script3
  │   │   └── runsub2.sh
  │   └── Processing_Pipeline/
  ├── Data/
  ├── Layout/
  └── run.sh

In this way it's easier to know what each script and functions are doing in the database scheme and to debug if needed. This code are designed to be run and generate all the data retrieved from the database in a machine or a cluster, for a range of usages.

After padronization of the retrieved data, all will pass through the main pipeline of the project and the final results will be sent to de Data directory, which stores all the raw data and the processed ones, the database entities and other useful information.

-MXD/
  ├── Scripts/
  ├── Data/
  │   ├── Entities/
  │   ├── Raw/
  │   └── Reports/
  ├── Layout/
  └── run.sh

The Layout directory refers to the front and backend of the software, using the data stored in the Data folder to organize the data and the filters in a intuitive way for the user. One can just download the source code and run the database in the machine or run from a server.

-MXD/
  ├── Scripts/
  ├── Data/
  ├── Layout/
  │   ├── frontend/
  │   └── backend/
  └── run.sh

## Data processing

The data retrieved from the mining process will pass through an universal processing pipeline, which contains the quality control, binning, annotation, xenobiotics comparison and data organizing.

### Filtering

The data will be submitted to a first filtering step, consisting in an seqkit analysis about the summary information of the produced bins. The higher outliers of the number of sequences distribution will be removed, since the more fragmented is the MAGs less trustful is the protein annotation.

### Quality acessing

The quality acessing will be performed by CheckM2 (see documentation), only the Medium quality or above bins will be selected to pass.

### Taxonomic annotation

In progress

### Protein annotation

The protein annotation will be made by a combination of annotators, the result will be compared and the best one will be chosen by each bin analyzed.

### Xenobiotics comparison

In progress


