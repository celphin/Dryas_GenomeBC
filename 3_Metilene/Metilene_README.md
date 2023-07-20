# Metilene - guide to scripts and notes

## Scripts:

### 1. metilene_prep.sh
***Function:*** 
Given a methylseq output directory, copies over files as specified, and renames them according to specified convention. The bedgraph files are then sorted and formatted for metilene input. The file names are then put into lists according to category. The correctly formatted bedgraphs are then put into a single file that is ready to be inputted into metilene. 

***Required modules:***
* StdEnv/2020
* bedtools/2.30.0

***Inputs:***
* h1 - group 1
* h2 - group 2
* methylseq_output_dir - name of location of all the methylseq output bedgraphs
* input_dir - name of directory where all the bedgraphs, name of files, and going into the metilene input file will be stored
* in_metilene - name of file to be inputted into metilene

***Running:***

sh metilene_prep.sh

***Modifications:***
* section called "copying files into directory" - modify so only desired files are copied over from specified directory. 
* section called "renaming files" - modify names to match format of either h1_filename or h2_filename (not essential) 
* detailed modifications specified in each notes file

***Outputs:***
* h1_list, h2_list - list of all the files created by group
* in_metilene file - file ready to be inputted to metilene or for use in plotting

### 2. metilene_run.sh

***Function:*** 
Runs metilene on the specified input file, with desired parameters

***Required modules:***
* StdEnv/2020
* bedtools/2.30.0

***Inputs: ***
* h1 - group 1
* h2 - group 2
* metilene_dir - name of directory where metilene is downloaded
* input_dir - name of directory where metilene input file is
* output_name - name of output metilene file
* in_metilene - name of input file 
* threads - number of threads used (no need to modify), but make sure allocation accounts for this

***Parameters:***
* maxdist - maximum distance between two cpgs
* mincpgs - minimum number cpgs within range
* mindiff - minium difference in methylation 

***Runnning:***

sh metilene_run.sh maxdist mincpgs mindiff

***Outputs***
* bedfile of the name output_name containing all the DMRS matching specifications in parameters

### 3. metilene_filter_qval.sh

***Function:*** 
Uses metilene's output processing program on the specified input file, to filter DMRs by a specified qval threshold.

***Required modules:***
* nixpkgs/16.09 
* gcc/7.3.0
*  r/3.6.0
* gdal
* udunits
* python
* R_LIBS_USER=/home/username/R/x86_64-pc-linux-gnu-library/3.6/

***Inputs:***
* h1 - group 1
* h2 - group 2
* metilene_dir - name of directory where metilene is downloaded
* input_dir - name of directory where metilene input file is
* output_name - name of output metilene file
* mincpgs=10 (don't modify this is for running qvalue filtering)

***Parameters:***
* maxdist - maximum distance between two cpgs
* mincpgs - minimum number cpgs within range
* mindiff - minimum difference in methylation
* minmeandiff 
* q val - q value

***Runnning:***

sh metilene_prep.sh maxdist mincpgs mindiff minmeandiff qval

***Outputs***
* .bed, .out, and .pdf file containing all the DMRS filtered according to specifications in parameters

## Notes Files:

### MS_Parent_Metilene_Notes:
Requests allocations, specifies modifications, loads modules, and runs the three scripts in the folder on:
* All wild plants for Warming versus control
* All true parent plants, matched too seedling ID, for warming versus control
* All wild plants for high Arctic sites versus low Arctic sites
* All wild plants for Dryas Octapetala versus Dryas Integrifolia

### MS_Wild_Metilene_Site_Specific_Notes
Requests allocations, specifies modifications, loads modules, and runs the three scripts in the folder on:
* Svalbard sites
* Sweden sites
* Nunavut sites
* Alaska sites

Currently (moved to bedgraph_intersect notes):
* Creates an intersection of all 4 sites
* Creates intersections of all groups of 3 sites

### MS_Seedling_Metilene_Notes
Requests allocations, specifies modifications, loads modules, and runs the three scripts in the folder on:
* All seedling plants for Warming versus Control
* All seedling plants for high Arctic sites versus low Arctic sites

### MS_Phenology_Metilene_Notes
Requests allocations, specifies modifications, loads modules, and runs the three scripts in the folder on:
* Samples for which data was collected for mature flowers versus senesence



## Todo for Folder:
1. Create file for metilene download and R download
2. Potentially make all scripts to sbatch runnable with slurm specifications and modules premade into the files?
3. Remove redundant cross over with bedgraph intersection file in site specific notes ??
4. Double check that all information is complete