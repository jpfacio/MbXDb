# Marine Xenobiotics Database

## Introduction

MXD (provisory name) is a database to store potential genes for xenobiotics biorremediation, designed for easy visualization and usage for all marine microbiologists, regardless of knowledge in informatics and database handling. 

All the data presented in MXD originated from data mining of previous marine metagenomics studies, aiming a worldwide access to taxonomical and function of marine microbial communities. All the DNA sequences retrieved from the studies were submitted to a curated pipeline containing filters of quality control, metagenomics binning, taxonomic accessing and functional annotation. The processed data is then divided into a relational database, containing all the information of the studies, oceanographic and climatic variables and genomic functions and clusters. 

## Code and directory structure

The scripts are all located in a Scripts folder, with daughter directories related to specific parts of the projects, each one of these parts are run by a sub main script, which are called by the general main script, the "run.sh". 

-MXD
    -Scripts
        -Subproject1
            -script1
            -script2
            -script3
            -runsub1.sh
        -Subproject2
            -script1
            -script2
            -script3
            -runsub2.sh
        ...
        -Processing_Pipeline
    -Data
    -run.sh

In this way it's easier to know what each script and functions are doing in the database scheme and to debug if needed. This code are designed to be run and generate all the data retrieved from the database in a machine or a cluster, for a range of usages.

After padronization of the retrieved data, all will pass through the main pipelina of the project and the final results will be sent to de Data directory.



