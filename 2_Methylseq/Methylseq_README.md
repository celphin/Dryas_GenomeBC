# Methylseq - guide to scripts and notes

## Scripts:

### 1. make_input_file.sh
***Function:*** 
Creates list of files and file names formatted for input into methylseq config file.

***Inputs:***
* None, but adjust names/locations of following for temporary/intermediate files 
* file_dir - file that will contain all 
* sample_list=" ~/scratch/testfilecopy/samples_list.csv"

***Running:***

sh make_input_file.sh

***Outputs:***
* input_files.csv - file ready to submit to Next flow's methylseq

### 2. rename_files.sh 
***Function:*** 
Uses newnames.csv to rename files into desire new names.

***Inputs:***
* newnames.csv - file of format newname, old name

***Running:***

sh rename_files.sh

## Notes Files:

### Notes1_SeedlingsPhenoMeth_Methylseq_Sept2022.txt

Prepares files for, configures, and runs Nextflow's methylseq on Mature flowers.

### Notes1a_WildWGBS_Methylseq_Mapping_Methcalling_Jan2022


Prepares files for, configures, and runs Nextflow's methylseq on all Wild Plants. Some post processing.


## To do:
* Insert Seedling notes
* Clean up all notes
* Write up clear instructions for methylseq current version ???
* Split up both notes to just be methylseq