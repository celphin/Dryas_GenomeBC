# Interproscan - guide to scripts and notes

## Note:
This will likely change very much. Once we get interproscan running on a more recent version we can change this all into a single script (aside from pre/post processing). Maybe reformat this file to just be a map nucleotide -> fasta file, will be cleaner. 
## Scripts:

### 1. prep_interproscan.sh
***Function:*** 
Maps nucleotide fasta file to protien file, and replaces special characters. 
May not be neccessary.   

***Required modules:***
* StdEnv/2020
* gcc/9.3.0
* r-bundle-bioconductor/3.16 

***Inputs:***
* ref_protein - reference genome with protiens
* input_dir - directory with BLAST output
* output_dir - whatever directory will be inputted to interproscan

***Running:***

sh prep_interproscan.sh

***Outputs:***

* Protien fasta files ready to be inputted into interproscan

## Notes Files:

### MS_Interproscan_Notes ###
* Prepares files for submission to interproscan 
* Creates file with list of all files to run through interproscan
* Runs it on entire fasta files that are relatively small via loop 
* Runs interproscan on larger files by splitting them and running them in a loop
* Messy notes on cleaning interproscan output


### MS_RNA_Interproscan
Requires: RNAseq DER data, Interproscan output
* Runs interproscan on DER genes
* Crosses over interproscan output for DMR and DER data
* Subsets go distributions

## Todo for Folder:
1. Modify prep file - mapping to protien is unnecessary. If no special characters in Oxyria nucleotide genome - submit it straight to interproscan, this file is unnecessary.  
2. Modify loop in directory to be a single script just using array jobs 
3. For larger files replace with array script 
4. Once new version is running replace this all to be done in one single script for all the files
5. Clean up notes on output cleaning
6. Once new interproscan running - make better notes on crossing DER and DMRS 