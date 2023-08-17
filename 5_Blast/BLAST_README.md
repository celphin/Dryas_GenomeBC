# BLAST - guide to scripts and notes

## Scripts:

### 1. blast_prep_extended.sh
***Function:*** 
Extends DMR regions by 500 bp on each side. If extended position is negative, it is changed to zero. Then the positions are mapped to the fasta file. 

***Required modules:***
* StdEnv/2020
* bedtools/2.30.0

***Inputs:***
* directory - directory where DMR bed file is located
* ref_genome - name of reference genome .fa file

***Parameters***
* blastname - name of bed file to be prepared for blasting

***Running:***

sh blast_prep_extended.sh blastname

***Outputs:***
* blastname.fasta - file of DMRS mapped to genome 

## Notes Files:

### MS_Blast_Notes ###
Requires: all the files from metilene folder
* Copies over metilene output - creates several intersections (REMOVE eventually - put to a different file)
* Prepares all bedfiles for submitting to BLAST
* Runs Blast for all files to: reference genome, rosaceae, arabidopsis
* Cleans reference blasts to only leave matches within the same chromosome

### MS_Process_Blast_Output
Requires: ncbi blast output for filtering
* Note: Messy notes - need to finish
* Filters out only gene results
* Leaves top 3 matches

## Todo for Folder:
1. Decide what to do with processing BLAST ncbi, clean up code. Maybe run only on list we have identified as significant.
2. Remove DMRS - and update file copy in this case