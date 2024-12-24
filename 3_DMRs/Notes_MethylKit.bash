###########################
# methylKit
# https://github.com/al2na/methylKit
# Dec 2024
##############################

# data

cd /home/celphin/scratch/Dryas/CpG/

tmux new-session -s methylkit
tmux attach-session -t methylkit

#-----------------------------
# ideal input data - methylation call file
##         chrBase   chr    base strand coverage freqC  freqT
## 1 chr21.9764539 chr21 9764539      R       12 25.00  75.00
## 2 chr21.9764513 chr21 9764513      R       12  0.00 100.00
## 3 chr21.9820622 chr21 9820622      F       13  0.00 100.00
## 4 chr21.9837545 chr21 9837545      F       11  0.00 100.00
## 5 chr21.9849022 chr21 9849022      F      124 72.58  27.42

# https://felixkrueger.github.io/Bismark/bismark/methylation_extraction/
cd /home/celphin/scratch/Dryas/CpG/stranded_CpG_report/

# <chromosome> <position> <strand> <count methylated> <count unmethylated> <C-context> <trinucleotide context>
# Do1_00107       9       +       0       1       CG      CGG
# Do1_00107       10      -       0       0       CG      CGA
# Do1_00107       138     +       0       5       CG      CGT
# Do1_00107       139     -       0       2       CG      CGT
# Do1_00107       248     +       1       7       CG      CGC
# Do1_00107       249     -       0       8       CG      CGT
# Do1_00107       573     +       0       2       CG      CGC
# Do1_00107       574     -       5       2       CG      CGA

gunzip *.gz

# add overall coverage
awk -F'\t' '{print $0 "\t" $4 + $5}' OFS='\t' W4.G10.W2d1_DRY9W_50_185_R1_val_1_bismark_bt2_pe.deduplicated.CpG_report.txt > \
W4.G10.W2d1_DRY9W_50_185_R1_val_1_bismark_bt2_pe.deduplicated.CpG_report1.txt

awk -F'\t' '{print $0 "\t" $4 + $5}' OFS='\t' C2.H12.C2g6_ALAS0C_18_229_R1_val_1_bismark_bt2_pe.deduplicated.CpG_report.txt > \
C2.H12.C2g6_ALAS0C_18_229_R1_val_1_bismark_bt2_pe.deduplicated.CpG_report1.txt

awk -F'\t' '{print $0 "\t" $4 + $5}' OFS='\t' C2.G12.C2a6_LATD1C_4_223_R1_val_1_bismark_bt2_pe.deduplicated.CpG_report.txt > \
C2.G12.C2a6_LATD1C_4_223_R1_val_1_bismark_bt2_pe.deduplicated.CpG_report1.txt

awk -F'\t' '{print $0 "\t" $4 + $5}' OFS='\t' W2.H11.W2a4_LATD4W_9_207_R1_val_1_bismark_bt2_pe.deduplicated.CpG_report.txt > \
W2.H11.W2a4_LATD4W_9_207_R1_val_1_bismark_bt2_pe.deduplicated.CpG_report1.txt

#---------------------------------
# run on all

directory="/home/celphin/scratch/Dryas/CpG/stranded_CpG_report"  # Replace with your directory path

# Loop through all .txt files in the directory
for file in "$directory"/*R1_val_1_bismark_bt2_pe.deduplicated.CpG_report.txt; do
  # Run the awk command on each file and save the output to a new file
  awk -F'\t' '{print $0 "\t" $4 + $5}' OFS='\t' "$file" > "${file%.txt}1.txt"
done


#-----------------
# to make stranded cytosine reports for CHH
cd /home/celphin/scratch/Dryas/CHG_CHH/

#!/bin/bash
#SBATCH --account=def-cronk
#SBATCH --time=0-5:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=120000M

module load StdEnv/2023 minimap2/2.28 samtools/1.20 perl/5.36.1

/home/celphin/scratch/Dryas/CHG_CHH/Bismark/bismark_methylation_extractor --help

# add overall coverage
awk -F'\t' '{print $0 "\t" $4 + $5}' OFS='\t' input_file.txt > output_file.txt

#############################
# Get ref genome chromosome lengths

# Sort the bedGraph file by chromosome (column 1) and then find the maximum value per chromosome
awk '{if ($3 > max_end[$1]) max_end[$1] = $3} END {for (chrom in max_end) print chrom, max_end[chrom]}' Chilliwack1.F112574.bedGraph | sort > Dryas_seqinfo.txt

# add 1000 to each

###################################
module load StdEnv/2023
module load r/4.4.0

R

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("methylKit")

###################################
# to run 
# https://bioconductor.org/packages/release/bioc/vignettes/methylKit/inst/doc/methylKit.html
# https://bioconductor.org/packages/release/bioc/vignettes/methylKit/inst/doc/methylKit.html#1_Introduction

## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(dpi = 75)
knitr::opts_chunk$set(cache = FALSE)

##########################
# define genome assembly

library(GenomicRanges)

# read in data
ref_data0 <- read.table("/home/celphin/scratch/Dryas/CHG_CHH/data/bedgraph/Dryas_seqinfo.txt", header = FALSE, sep = " ")
ref_data <- head(ref_data0, -1)

# Define the chromosomes and their lengths (custom genome)
chromosomes <- ref_data$V1
lengths <- ref_data$V2  # Lengths for each chromosome

new_lengths <- lengths + 10000

# Create a Seqinfo object with chromosome information
Dryas_seqinfo <- Seqinfo(seqnames = chromosomes, seqlengths = new_lengths,  genome = "Dryas_genome" )

# View the Seqinfo object
Dryas_seqinfo

###############################
# load data

library(methylKit)

# Set the working directory to the folder containing the text files
setwd("/home/celphin/scratch/Dryas/CpG/stranded_CpG_report")

# List all text files in the folder with a specific pattern 
files <- list.files(pattern = "*CpG_report.txt$")

# Define sample ids (replace with your actual sample names)
sample_ids <- sub("_R1_val_1_bismark_bt2_pe\\.deduplicated\\.CpG_report\\.txt$", "", files)

covariate <- sub("^[^_]*_([A-Za-z]{3}).*", "\\1", sample_ids)

# Define treatment conditions (0 = control, 1 = treated)
treatment <- c(rep(0, length(files)/2), rep(1, length(files)/2))

# Extract first letter (C or W)
first_letter <- substr(sample_ids, 1, 1)

# Replace C with 0 and W with 1
treatment <- gsub("C", "0", first_letter)
treatment <- gsub("W", "1", treatment)
treatment <- gsub("L", "0", treatment)

# Coverage header: <chromosome> <start position> <end position> <methylation percentage> <count methylated> <count unmethylated>
# https://rdrr.io/bioc/methylKit/man/methRead-methods.html

d.f <- as.data.frame(cbind(sample_ids, treatment, covariate))

# Read all the files into methylRawList using methRead
myobj <- methRead(
  location = as.list(files), 
  sample.id = as.list(d.f$sample_ids), 
  assembly = "Dryas",
  treatment = as.numeric(d.f$treatment),
  dbdir = getwd(),
  context = "CpG", 
  dbtype = "tabix",
  resolution = "base", 
  mincov = 10,
  pipeline = "bismarkCytosineReport" 
)

# creating directory /lustre04/scratch/celphin/Dryas/CpG/stranded_CpG_report/methylDB_2024-12-17_9vt
# Reading file.
# compressing the file with bgzip...
# making tabix index...
# Reading file.

#------------------
# Print the object to check the result
print(myobj)

save(myobj, file = "MethylKit_data.RData")

load("MethylKit_data.RData")
## -----------------------------------------------------------------------------
getMethylationStats(myobj[[2]],plot=FALSE,both.strands=FALSE)

# methylation statistics per base
# summary:
   # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
   # 0.00    0.00   36.36   43.41   90.00  100.00
# percentiles:
       # 0%       10%       20%       30%       40%       50%       60%       70%
  # 0.00000   0.00000   0.00000   0.00000   0.00000  36.36364  71.42857  84.61538
      # 80%       90%       95%       99%     99.5%     99.9%      100%
 # 91.66667 100.00000 100.00000 100.00000 100.00000 100.00000 100.00000

## -----------------------------------------------------------------------------
getMethylationStats(myobj[[2]],plot=TRUE,both.strands=FALSE)

## -----------------------------------------------------------------------------
getCoverageStats(myobj[[2]],plot=FALSE,both.strands=FALSE)

# read coverage statistics per base
# summary:
   # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
  # 10.00   12.00   14.00   24.06   17.00 8144.00
# percentiles:
   # 0%   10%   20%   30%   40%   50%   60%   70%   80%   90%   95%   99% 99.5%
   # 10    10    11    12    13    14    15    16    18    21    27    71   258
# 99.9%  100%
 # 3901  8144

## -----------------------------------------------------------------------------
getCoverageStats(myobj[[2]],plot=TRUE,both.strands=FALSE)

## -----------------------------------------------------------------------------
filtered.myobj=filterByCoverage(myobj,lo.count=10,lo.perc=NULL,
                                      hi.count=NULL,hi.perc=99.9)

save(filtered.myobj, file = "FilteredMethylKit_data.RData")
load("FilteredMethylKit_data.RData")

## -----------------------------------------------------------------------------
# destrand=TRUE for CpG but not CHG/CHH

meth=unite(filtered.myobj, destrand=TRUE)

## -----------------------------------------------------------------------------
head(meth)
save(meth, file = "united_MethylKit_data.RData")
load("united_MethylKit_data.RData")

## ----eval=FALSE---------------------------------------------------------------
#  # creates a methylBase object,
#  # where only CpGs covered with at least 1 sample per group will be returned
#  
#  # there were two groups defined by the treatment vector,
#  # given during the creation of myobj: treatment=c(1,1,0,0)
#  meth.min=unite(myobj,min.per.group=1L)

## -----------------------------------------------------------------------------
getCorrelation(meth,plot=FALSE)

# Save the plot to a PDF file
pdf("Correlation_plot.pdf", width = 10, height = 8)
  print(getCorrelation(meth,plot=TRUE))
dev.off()
# will not plot due to margins

## -----------------------------------------------------------------------------
# Save the plot to a PDF file
pdf("Dendrogram_plot.pdf", width = 10, height = 8)
clusterSamples(meth, dist="correlation", method="ward.D2", plot=TRUE)
dev.off()

## ----message=FALSE------------------------------------------------------------
hc = clusterSamples(meth, dist="correlation", method="ward.D2", plot=FALSE)
save(hc, file = "dendrogram.RData")
load("dendrogram.RData")

## -----------------------------------------------------------------------------
pdf("PCA_Scree_plot.pdf", width = 10, height = 8)
PCASamples(meth, screeplot=TRUE)
dev.off()

# too large cannot make plot
# Error: cannot allocate vector of size 365.9 Mb


## -----------------------------------------------------------------------------
# https://rdrr.io/github/al2na/methylKit/man/PCASamples-methods.html

pdf("Wild_PCA_PC1_PC2.pdf", width = 10, height = 8)
PCASamples(meth)
dev.off()

pdf("Wild_PCA_PC3_PC4.pdf", width = 10, height = 8)
PCASamples(meth, comp=c(3,4))
dev.off()

pdf("Wild_PCA_PC5_PC6.pdf", width = 10, height = 8)
PCASamples(meth, comp=c(5,6))
dev.off()

#------------------------------
# make nicer PCA
load("united_MethylKit_data.RData")

allDataPCA <- PCASamples(meth, obj.return = TRUE) #Run a PCA analysis on percent methylation for all samples. methylKit uses prcomp to create the PCA matrix
summary(allDataPCA)


d.f$mix <- paste0(d.f$treatment, "_", d.f$covariate)

library(vegan)
pdf("Total_PCA_PC1_PC2_ordiplot.pdf", width = 10, height = 8)
fig.allDataPCA <- ordiplot(allDataPCA, choices = c(1, 2), type = "none", display = "sites", cex = 0.5, xlim = c(-400, 200), xlab = "", ylab = "", xaxt = "n", yaxt = "n") #Use ordiplot to create base biplot. Do not add any points
#points(fig.allDataPCA, "sites", col = c(rep("blue", times = 50), rep("red", times = 50)), pch = c(rep(16, times = 50), rep(17, times = 50)), cex = 3) #Add each sample. Darker samples are ambient, lighter samples are elevated pCO2

ordiellipse(allDataPCA, d.f$mix, show.groups = "1_ALA", col = "red") #Add confidence ellipse around the samples in elevated pCO2
ordiellipse(allDataPCA, d.f$mix, show.groups = "0_ALA", col = "blue") #Add confidence ellipse around the samples in ambient pCO2

ordiellipse(allDataPCA, d.f$mix, show.groups = "1_LAT", col = "red") #Add confidence ellipse around the samples in elevated pCO2
ordiellipse(allDataPCA, d.f$mix, show.groups = "0_LAT", col = "blue") #Add confidence ellipse around the samples in ambient pCO2

ordiellipse(allDataPCA, d.f$mix, show.groups = "1_MEA", col = "red") #Add confidence ellipse around the samples in elevated pCO2
ordiellipse(allDataPCA, d.f$mix, show.groups = "0_MEA", col = "blue") #Add confidence ellipse around the samples in ambient pCO2

ordiellipse(allDataPCA, d.f$mix, show.groups = "1_SVA", col = "red") #Add confidence ellipse around the samples in elevated pCO2
ordiellipse(allDataPCA, d.f$mix, show.groups = "0_SVA", col = "blue") #Add confidence ellipse around the samples in ambient pCO2


axis(side =  1, at = seq(-400, 200, 200), col = "grey80", cex.axis = 1.7) #Add x-axis
mtext(side = 1, text = "PC 1", line = 3, cex = 1.5) #Add x-axis label

axis(side =  2, labels = TRUE, col = "grey80", cex.axis = 1.7) #Add y-axis
mtext(side = 2, text = "PC 2", line = 3, cex = 1.5) #Add y-axis label

dev.off()


############################################
# relaunch
salloc -c1 --time 2:00:00 --mem 120000m --account def-cronk

library(methylKit)

# Set the working directory to the folder containing the text files
setwd("/home/celphin/scratch/Dryas/CpG/stranded_CpG_report")

# List all text files in the folder with a specific pattern 
files <- list.files(pattern = "*CpG_report.txt$")

# Define sample ids (replace with your actual sample names)
sample_ids <- sub("_R1_val_1_bismark_bt2_pe\\.deduplicated\\.CpG_report\\.txt$", "", files)

covariate <- sub("^[^_]*_([A-Za-z]{3}).*", "\\1", sample_ids)

# Define treatment conditions (0 = control, 1 = treated)
treatment <- c(rep(0, length(files)/2), rep(1, length(files)/2))

# Extract first letter (C or W)
first_letter <- substr(sample_ids, 1, 1)

# Replace C with 0 and W with 1
treatment <- gsub("C", "0", first_letter)
treatment <- gsub("W", "1", treatment)
treatment <- gsub("L", "0", treatment)

d.f <- as.data.frame(cbind(sample_ids, treatment, covariate))

load("MethylKit_data.RData")
load("FilteredMethylKit_data.RData")
load("united_MethylKit_data.RData")

#-------------------------------
# check association of random factors (e.g. site with PCA)

sampleAnnotation=data.frame(batch_id=covariate)

as=assocComp(mBase=meth,sampleAnnotation)
as

# $association
                  # PC1          PC2          PC3          PC4       PC5
# batch_id 1.650368e-11 4.705047e-16 3.497492e-16 8.506409e-16 0.2043553
                  # PC6       PC7       PC8        PC9      PC10      PC11
# batch_id 0.0007716853 0.2305742 0.7289346 0.02385734 0.1009379 0.7571989
                 # PC12         PC13         PC14         PC15      PC16
# batch_id 7.553552e-12 6.998494e-09 1.251649e-05 7.383688e-05 0.2411195
                # PC17     PC18      PC19         PC20      PC21      PC22
# batch_id 0.003351971 0.153593 0.5935674 0.0002988655 0.6161314 0.1887075
              # PC23      PC24      PC25      PC26     PC27      PC28      PC29
# batch_id 0.3119683 0.5665974 0.9517783 0.3087377 0.908184 0.2831733 0.1650447
             # PC30      PC31      PC32      PC33      PC34     PC35      PC36
# batch_id 0.336769 0.2951237 0.4239164 0.1858227 0.7689914 0.143244 0.2648281
              # PC37       PC38      PC39      PC40      PC41      PC42
# batch_id 0.5775444 0.06314453 0.8409517 0.3585235 0.5721591 0.4390139
               # PC43      PC44      PC45      PC46      PC47      PC48      PC49
# batch_id 0.08366977 0.0952027 0.7736787 0.3366397 0.6792229 0.8235588 0.2667659
              # PC50      PC51      PC52     PC53      PC54      PC55      PC56
# batch_id 0.1941657 0.5197683 0.7241807 0.470371 0.5915692 0.2434096 0.8763353
              # PC57      PC58      PC59      PC60      PC61      PC62      PC63
# batch_id 0.5967877 0.4952361 0.8211694 0.3116649 0.8741406 0.8010502 0.3033271
              # PC64      PC65      PC66      PC67      PC68    PC69      PC70
# batch_id 0.3962544 0.3620019 0.7913348 0.1318804 0.8732317 0.21293 0.3690799
             # PC71     PC72      PC73      PC74     PC75      PC76       PC77
# batch_id 0.862756 0.234529 0.9992219 0.6002552 0.201419 0.4470045 0.08001973
               # PC78      PC79      PC80       PC81      PC82      PC83
# batch_id 0.02370628 0.4678206 0.2934271 0.09474636 0.1573261 0.2425966
                # PC84      PC85        PC86       PC87      PC88      PC89
# batch_id 0.004182407 0.9060866 0.008126297 0.07838719 0.7665706 0.5992398
              # PC90     PC91      PC92      PC93      PC94      PC95      PC96
# batch_id 0.6206974 0.654191 0.9913613 0.9567689 0.3804815 0.9857119 0.8529406
              # PC97      PC98      PC99     PC100     PC101     PC102    PC103
# batch_id 0.5608475 0.9414458 0.7792603 0.8434796 0.6284271 0.9789979 0.776521


#-----------------------------
# construct a new object by removing the first pricipal component
# from percent methylation value matrix

newObj=removeComp(meth,comp=c(1,2,3,4,6), dbtype="tabix")
# ran out of memory

# PCASamples(newObj)

## -----------------------------------------------------------------------------
# # correct for batch effects

# mat=percMethylation(meth)

# # do some changes in the matrix
# # this is just a toy example
# # ideally you want to correct the matrix
# # for batch effects
# mat[mat==100]=80
 
# # reconstruct the methylBase from the corrected matrix
# newobj=reconstruct(mat,meth)

## -----------------------------------------------------------------------------
# DMC and DMRs
# https://rdrr.io/github/al2na/methylKit/man/calculateDiffMeth-methods.html

myDiff=calculateDiffMeth(meth)
save(myDiff, file = "MethDiffCs.RData")
load("MethDiffCs.RData")

# 1: glm.fit: algorithm did not converge
# 2: glm.fit: algorithm did not converge

## -----------------------------------------------------------------------------
# get hyper methylated bases
myDiff5p.hyper=getMethylDiff(myDiff,difference=5,qvalue=0.01,type="hyper")

# get hypo methylated bases
myDiff5p.hypo=getMethylDiff(myDiff,difference=5,qvalue=0.01,type="hypo")

# get all differentially methylated bases
myDiff5p=getMethylDiff(myDiff,difference=5,qvalue=0.01)

nrow(getData(myDiff5p))
nrow(getData(myDiff5p.hypo))
nrow(getData(myDiff5p.hyper))

# [1] 45176
# [1] 21298
# [1] 23878

#-------------------
myDiff10p.hyper=getMethylDiff(myDiff,difference=10,qvalue=0.01,type="hyper")

# get hypo methylated bases
myDiff10p.hypo=getMethylDiff(myDiff,difference=10,qvalue=0.01,type="hypo")

# get all differentially methylated bases
myDiff10p=getMethylDiff(myDiff,difference=10,qvalue=0.01)

nrow(getData(myDiff10p))
nrow(getData(myDiff10p.hypo))
nrow(getData(myDiff10p.hyper))

# [1] 7091
# [1] 3302
# [1] 3789

## -----------------------------------------------------------------------------
diffMethPerChr(myDiff,plot=TRUE,qvalue.cutoff=0.01, meth.cutoff=7)

pdf("MethDiff_scaffolds.pdf", width = 10, height = 8)
diffMethPerChr(myDiff,plot=TRUE,qvalue.cutoff=0.01, meth.cutoff=7)
dev.off()

#-------------------------
# Over dispersion correction

sim.methylBase1<-dataSim(replicates=6,sites=1000,
                         treatment=c(rep(1,3),rep(0,3)),
                        sample.ids=c(paste0("test",1:3),paste0("ctrl",1:3))
                        )

my.diffMeth<-calculateDiffMeth(sim.methylBase1[1:,],
                                overdispersion="MN",test="Chisq",mc.cores=1)



#-------------------------------
# add site as a covariate

covariates=data.frame(site=covariate)
sim.methylBase<-dataSim(replicates=6,sites=1000,
                        treatment=c(rep(1,3),rep(0,3)),
                        covariates=covariates,
                        sample.ids=c(paste0("test",1:3),paste0("ctrl",1:3))
                        )
my.diffMeth3<-calculateDiffMeth(sim.methylBase,
                                covariates=covariates,
                                overdispersion="MN",test="Chisq",mc.cores=1)

## ----------------------------------------------------------------
# regions vs per base

library(methylKit)

# Set the working directory to the folder containing the text files
setwd("/home/celphin/scratch/Dryas/CpG/stranded_CpG_report")

# List all text files in the folder with a specific pattern 
files <- list.files(pattern = "*CpG_report.txt$")

# Define sample ids (replace with your actual sample names)
sample_ids <- sub("_R1_val_1_bismark_bt2_pe\\.deduplicated\\.CpG_report\\.txt$", "", files)

covariate <- sub("^[^_]*_([A-Za-z]{3}).*", "\\1", sample_ids)

# Define treatment conditions (0 = control, 1 = treated)
treatment <- c(rep(0, length(files)/2), rep(1, length(files)/2))

# Extract first letter (C or W)
first_letter <- substr(sample_ids, 1, 1)

# Replace C with 0 and W with 1
treatment <- gsub("C", "0", first_letter)
treatment <- gsub("W", "1", treatment)
treatment <- gsub("L", "0", treatment)

# Coverage header: <chromosome> <start position> <end position> <methylation percentage> <count methylated> <count unmethylated>
# https://rdrr.io/bioc/methylKit/man/methRead-methods.html

d.f <- as.data.frame(cbind(sample_ids, treatment, covariate))

# Read all the files into methylRawList using methRead
myobj_lowCov <- methRead(
  location = as.list(files), 
  sample.id = as.list(d.f$sample_ids), 
  assembly = "Dryas",
  treatment = as.numeric(d.f$treatment),
  dbdir = getwd(),
  context = "CpG", 
  dbtype = "tabix",
  resolution = "base", 
  mincov = 3,
  pipeline = "bismarkCytosineReport" 
)
# creating directory /lustre04/scratch/celphin/Dryas/CpG/stranded_CpG_report/methylDB_2024-12-18_8PF

# Print the object to check the result
print(myobj_lowCov)
save(myobj_lowCov, file = "MethylKit_lowcov_data.RData")
load("MethylKit_lowcov_data.RData")

# make regions
lowtiles = tileMethylCounts(myobj_lowCov,win.size=300,step.size=300,cov.bases = 10)
head(lowtiles[[1]],3)
save(lowtiles, file = "MethylKit_tiles300_lowcov.RData")
load("MethylKit_tiles300_lowcov.RData")

#------------------------
methtiles=unite(lowtiles, destrand=FALSE)
save(methtiles, file = "united_lowtiles300.RData")
load("united_lowtiles300.RData")

#---------------------------
# calculate DMRs

myDifftiles=calculateDiffMeth(methtiles)
save(myDifftiles, file = "MethDiffRegions_300.RData")
load("MethDiffRegions_300.RData")

# get hyper methylated bases
myDifftiles10p.hyper=getMethylDiff(myDifftiles,difference=10,qvalue=0.01,type="hyper")

# get hypo methylated bases
myDifftiles10p.hypo=getMethylDiff(myDifftiles,difference=10,qvalue=0.01,type="hypo")

# get all differentially methylated bases
myDifftiles10p=getMethylDiff(myDifftiles,difference=10,qvalue=0.01)

## -----------------------------------------------------------------------------
diffMethPerChr(myDifftiles,plot=FALSE,qvalue.cutoff=0.01, meth.cutoff=10)

nrow(getData(myDifftiles10p))
nrow(getData(myDifftiles10p.hypo))
nrow(getData(myDifftiles10p.hyper))

#1000
# [1] 82
# [1] 36
# [1] 46

# 300p
# [1] 211
# [1] 95
# [1] 116


Wild_W_C_lowcovDMRs <- getData(myDifftiles10p)
write.table(Wild_W_C_lowcovDMRs, "Methylkit_Wild_W_C_lowcovDMRs.txt", sep = "\t", quote = FALSE, row.names = FALSE)


############################
# run with higher cov data for now
load("MethylKit_data.RData")
# https://rdrr.io/github/al2na/methylKit/man/tileMethylCounts-methods.html
tiles = tileMethylCounts(myobj,win.size=1000,step.size=1000,cov.bases = 10)
head(tiles[[1]],3)

#------------------------
methtiles=unite(tiles, destrand=FALSE)
save(methtiles, file = "united_tiles.RData")
load("united_tiles.RData")

#---------------------------
# calculate DMRs

myDifftiles=calculateDiffMeth(methtiles)
save(myDifftiles, file = "MethDiffRegions.RData")
load("MethDiffRegions.RData")


# get hyper methylated bases
myDifftiles10p.hyper=getMethylDiff(myDifftiles,difference=10,qvalue=0.01,type="hyper")

# get hypo methylated bases
myDifftiles10p.hypo=getMethylDiff(myDifftiles,difference=10,qvalue=0.01,type="hypo")

# get all differentially methylated bases
myDifftiles10p=getMethylDiff(myDifftiles,difference=10,qvalue=0.01)

## -----------------------------------------------------------------------------
diffMethPerChr(myDifftiles,plot=FALSE,qvalue.cutoff=0.01, meth.cutoff=10)

nrow(getData(myDifftiles10p))
nrow(getData(myDifftiles10p.hypo))
nrow(getData(myDifftiles10p.hyper))

Wild_W_C_DMRs <- getData(myDifftiles10p)
write.table(Wild_W_C_DMRs, "Methylkit_Wild_W_C_DMRs.txt", sep = "\t", quote = FALSE, row.names = FALSE)

             # chr    start      end strand        pvalue        qvalue meth.diff
# 1  Do1_01_a00001  4157001  4158000      * 2.901521e-234 5.128805e-232  10.02695
# 2  Do1_01_a00001  9180001  9181000      * 1.443590e-162 1.896554e-160 -10.86751
# 3  Do1_01_a00007  1469001  1470000      * 2.075587e-159 2.666795e-157 -10.96074
# 4  Do1_02_a00003  9872001  9873000      *  3.432551e-83  2.143752e-81  11.30865
# 5  Do1_02_a00004   851001   852000      * 4.274134e-160 5.515871e-158  10.45121
# 6  Do1_02_a00004  2304001  2305000      *  0.000000e+00  0.000000e+00  13.52473
# 7  Do1_02_a00004  4585001  4586000      *  9.979348e-63  4.299203e-61  10.02458
# 8  Do1_03_a00001  5153001  5154000      * 9.251670e-176 1.316257e-173  10.56226
# 9  Do1_03_a00002  5353001  5354000      *  0.000000e+00  0.000000e+00 -10.45337
# 10 Do1_03_a00002  5354001  5355000      * 4.940656e-324 1.022716e-321 -10.29152
# 11 Do1_03_a00002 10566001 10567000      * 1.467586e-196 2.264728e-194 -13.03470
# 12 Do1_03_a00004   878001   879000      * 6.744355e-189 1.019196e-186  10.55043
# 13 Do1_04_a00001  7587001  7588000      *  0.000000e+00  0.000000e+00  10.97409
# 14 Do1_04_a00002  1290001  1291000      * 1.834590e-147 2.148890e-145  10.68644
# 15 Do1_04_a00002  1291001  1292000      * 3.859869e-211 6.219686e-209  11.18434
# 16 Do1_04_a00002  2404001  2405000      *  0.000000e+00  0.000000e+00 -10.56695
# 17 Do1_04_a00002  2662001  2663000      * 4.367915e-157 5.491117e-155 -11.78146
# 18 Do1_04_a00005   588001   589000      * 5.261911e-176 7.522943e-174  10.16907
# 19 Do1_05_a00001 11338001 11339000      * 6.406090e-151 7.720620e-149  10.53493
# 20 Do1_05_a00001 13285001 13286000      * 2.900313e-226 5.035123e-224 -12.29982
# 21 Do1_05_a00003   840001   841000      * 1.489464e-164 1.974612e-162  10.81079
# 22 Do1_06_a00001  2923001  2924000      * 2.204060e-160 2.857034e-158 -10.17589
# 23 Do1_06_a00001  3893001  3894000      * 7.169829e-263 1.340475e-260  13.47419
# 24 Do1_06_a00001 11467001 11468000      * 4.797479e-155 5.879096e-153 -10.86645
# 25 Do1_06_a00002  2718001  2719000      *  0.000000e+00  0.000000e+00  12.08077
# 26 Do1_06_a00002  8275001  8276000      * 3.863094e-224 6.627668e-222  14.69934
# 27 Do1_07_a00002  3726001  3727000      * 1.756893e-143 1.978426e-141  10.39217
# 28 Do1_07_a00004  3326001  3327000      * 7.726279e-218 1.280360e-215  12.87097
# 29 Do1_07_a00004  4633001  4634000      *  3.857163e-96  2.929620e-94  10.08782
# 30 Do1_07_a00004  6548001  6549000      *  0.000000e+00  0.000000e+00  11.77154
# 31    Do1_a00028   195001   196000      *  0.000000e+00  0.000000e+00  13.06707

#-------------------------
# output the bedgraph files for DMC and DMR

awk '{print $1, $2, $3, $NF}' Methylkit_Wild_W_C_DMRs.txt > Methylkit_Wild_W_C_DMRs.bedGraph
sed -i 's/ /\t/g' Methylkit_Wild_W_C_DMRs.bedGraph
sed -i 's/chr\tstart\tend\tmeth.diff/track type=bedGraph/g' Methylkit_Wild_W_C_DMRs.bedGraph

# intersect with metilene output

module load StdEnv/2023 bedtools/2.31.0
bedtools intersect -u -a /home/celphin/scratch/Dryas/CpG/stranded_CpG_report/Methylkit_Wild_W_C_DMRs.bedGraph \
-b /home/celphin/scratch/Dryas/CpG/Wild_W_C_Metilene/Wild_W_C_150_5_4_0.9_qval.1e-5.bedgraph

# Do1_01_a00001   9180001 9181000 -10.8675079933135
# Do1_01_a00007   1469001 1470000 -10.9607368835372
# Do1_02_a00004   2304001 2305000 13.524726691852
# Do1_03_a00002   5353001 5354000 -10.4533739880911
# Do1_03_a00002   5354001 5355000 -10.2915198843922
# Do1_04_a00001   7587001 7588000 10.9740874201469
# Do1_04_a00002   2662001 2663000 -11.7814644387254
# Do1_04_a00005   588001  589000  10.1690690889332
# Do1_05_a00001   11338001        11339000        10.5349261806715
# Do1_05_a00001   13285001        13286000        -12.299816197244
# Do1_05_a00003   840001  841000  10.8107852722255
# Do1_06_a00001   3893001 3894000 13.4741932005561
# Do1_06_a00002   2718001 2719000 12.0807706327395
# Do1_06_a00002   8275001 8276000 14.699340662143
# Do1_07_a00004   6548001 6549000 11.7715366309823
# Do1_a00028      195001  196000  13.0670652384072

bedtools intersect -u -a /home/celphin/scratch/Dryas/CpG/stranded_CpG_report/Methylkit_Wild_W_C_DMRs.bedGraph \
-b /home/celphin/scratch/Dryas/CpG/Wild_W_C_Metilene/Wild_W_C_70_5_4_0.9_qval.0.001.bedgraph

# Do1_01_a00001   9180001 9181000 -10.8675079933135
# Do1_03_a00002   5353001 5354000 -10.4533739880911
# Do1_03_a00002   5354001 5355000 -10.2915198843922
# Do1_03_a00002   10566001        10567000        -13.034702563676
# Do1_03_a00004   878001  879000  10.5504250815761
# Do1_04_a00002   2662001 2663000 -11.7814644387254
# Do1_04_a00005   588001  589000  10.1690690889332
# Do1_05_a00001   11338001        11339000        10.5349261806715
# Do1_05_a00001   13285001        13286000        -12.299816197244
# Do1_05_a00003   840001  841000  10.8107852722255
# Do1_06_a00001   3893001 3894000 13.4741932005561
# Do1_06_a00002   2718001 2719000 12.0807706327395
# Do1_07_a00004   6548001 6549000 11.7715366309823
# Do1_a00028      195001  196000  13.0670652384072

bedtools intersect -u -a /home/celphin/scratch/Dryas/CpG/stranded_CpG_report/Methylkit_Wild_W_C_DMRs.bedGraph \
-b /home/celphin/scratch/Dryas/CpG/Wild_W_C_Metilene_random/Wild_W_C_150_5_4_0.9_qval.1e-5.bedgraph

# Do1_01_a00001   9180001 9181000 -10.8675079933135
# Do1_01_a00007   1469001 1470000 -10.9607368835372
# Do1_02_a00004   2304001 2305000 13.524726691852
# Do1_03_a00002   5353001 5354000 -10.4533739880911 Do1_03_a00002G00984
# Do1_03_a00002   5354001 5355000 -10.2915198843922
# Do1_04_a00001   7587001 7588000 10.9740874201469
# Do1_04_a00002   2662001 2663000 -11.7814644387254
# Do1_04_a00005   588001  589000  10.1690690889332
# Do1_05_a00001   11338001        11339000        10.5349261806715
# Do1_05_a00001   13285001        13286000        -12.299816197244
# Do1_05_a00003   840001  841000  10.8107852722255
# Do1_06_a00001   3893001 3894000 13.4741932005561
# Do1_06_a00002   2718001 2719000 12.0807706327395
# Do1_06_a00002   8275001 8276000 14.699340662143
# Do1_07_a00004   6548001 6549000 11.7715366309823
# Do1_a00028      195001  196000  13.0670652384072


#####################################
# try running again for each site independently - set meth diff to 25 or 50%

cd /home/celphin/scratch/Dryas/CpG/

tmux new-session -s methylkit
tmux attach-session -t methylkit

module load StdEnv/2023
module load r/4.4.0

R

## ----------------------------------------------------------------
# DMRs for Alaska

library(methylKit)

# Set the working directory to the folder containing the text files
setwd("/home/celphin/scratch/Dryas/CpG/stranded_CpG_report")

# List all text files in the folder with a specific pattern 
files <- list.files(pattern = "*CpG_report.txt$")

# Define sample ids (replace with your actual sample names)
sample_ids <- sub("_R1_val_1_bismark_bt2_pe\\.deduplicated\\.CpG_report\\.txt$", "", files)

covariate <- sub("^[^_]*_([A-Za-z]{3}).*", "\\1", sample_ids)

# Define treatment conditions (0 = control, 1 = treated)
treatment <- c(rep(0, length(files)/2), rep(1, length(files)/2))

# Extract first letter (C or W)
first_letter <- substr(sample_ids, 1, 1)

# Replace C with 0 and W with 1
treatment <- gsub("C", "0", first_letter)
treatment <- gsub("W", "1", treatment)
treatment <- gsub("L", "0", treatment)

# Coverage header: <chromosome> <start position> <end position> <methylation percentage> <count methylated> <count unmethylated>
# https://rdrr.io/bioc/methylKit/man/methRead-methods.html

d.f <- as.data.frame(cbind(sample_ids, treatment, covariate))

ALAS_d.f <- d.f[which(d.f$covariate=="ALA"),]

load("MethylKit_lowcov_data.RData")

ALAS_myobj <- reorganize(
  myobj_lowCov,
  as.list(ALAS_d.f$sample_ids),
  as.numeric(ALAS_d.f$treatment),
  chunk.size = 1e+06,
  save.db = TRUE,
  suffix = "ALAS",
  dbtype = "tabix"
)

# Print the object to check the result
print(ALAS_myobj)
save(ALAS_myobj, file = "ALAS_MethylKit_lowcov_data.RData")
load("ALAS_MethylKit_lowcov_data.RData")

meth=unite(ALAS_myobj, destrand=FALSE)
save(meth, file = "ALAS_united.RData")
load("ALAS_united.RData")
pdf("ALAS_PCA_PC1_PC2.pdf", width = 10, height = 8)
PCASamples(meth)
dev.off()
pdf("ALAS_Dendrogram_plot.pdf", width = 10, height = 8)
clusterSamples(meth, dist="correlation", method="ward.D2", plot=TRUE)
dev.off()

#------------------------------
# make nicer PCA

allDataPCA <- PCASamples(meth, obj.return = TRUE) #Run a PCA analysis on percent methylation for all samples. methylKit uses prcomp to create the PCA matrix
summary(allDataPCA)

library(vegan)
pdf("ALAS_PCA_PC1_PC2_ordiplot.pdf", width = 10, height = 8)
fig.allDataPCA <- ordiplot(allDataPCA, choices = c(1, 2), type = "none", display = "sites", cex = 0.5, xlim = c(-400, 200), xlab = "", ylab = "", xaxt = "n", yaxt = "n") #Use ordiplot to create base biplot. Do not add any points
points(fig.allDataPCA, "sites", col = c(rep("blue", times = 10), rep("red", times = 10)), pch = c(rep(16, times = 10), rep(17, times = 10)), cex = 3) #Add each sample. Darker samples are ambient, lighter samples are elevated pCO2

ordiellipse(allDataPCA, ALAS_d.f$treatment, show.groups = "1", col = "red") #Add confidence ellipse around the samples in elevated pCO2
ordiellipse(allDataPCA, ALAS_d.f$treatment, show.groups = "0", col = "blue") #Add confidence ellipse around the samples in ambient pCO2

axis(side =  1, at = seq(-400, 200, 200), col = "grey80", cex.axis = 1.7) #Add x-axis
mtext(side = 1, text = "PC 1", line = 3, cex = 1.5) #Add x-axis label

axis(side =  2, labels = TRUE, col = "grey80", cex.axis = 1.7) #Add y-axis
mtext(side = 2, text = "PC 2", line = 3, cex = 1.5) #Add y-axis label

mtext(side = 3, line = -5, adj = c(-100,0), text = "    a. All CpG Loci", cex = 1.5)

dev.off()

#-----------------------
# make regions
lowtiles = tileMethylCounts(ALAS_myobj,win.size=1000,step.size=1000,cov.bases = 10)
head(lowtiles[[1]],3)
save(lowtiles, file = "ALAS_MethylKit_tiles_lowcov.RData")
load("ALAS_MethylKit_tiles_lowcov.RData")

#------------------------
methtiles=unite(lowtiles, destrand=FALSE)
save(methtiles, file = "ALAS_united_lowtiles.RData")
load("ALAS_united_lowtiles.RData")

#---------------------------
# calculate DMRs

myDifftiles=calculateDiffMeth(methtiles)
save(myDifftiles, file = "ALAS_MethDiffRegions.RData")
load("ALAS_MethDiffRegions.RData")

# get hyper methylated bases
myDifftiles10p.hyper=getMethylDiff(myDifftiles,difference=10,qvalue=0.01,type="hyper")

# get hypo methylated bases
myDifftiles10p.hypo=getMethylDiff(myDifftiles,difference=10,qvalue=0.01,type="hypo")

# get all differentially methylated bases
myDifftiles10p=getMethylDiff(myDifftiles,difference=10,qvalue=0.01)

## -----------------------------------------------------------------------------
diffMethPerChr(myDifftiles,plot=FALSE,qvalue.cutoff=0.01, meth.cutoff=10)

nrow(getData(myDifftiles10p))
nrow(getData(myDifftiles10p.hypo))
nrow(getData(myDifftiles10p.hyper))

# [1] 5085
# [1] 2609
# [1] 2476

#--------------------------
# get hyper methylated bases
myDifftiles25p.hyper=getMethylDiff(myDifftiles,difference=25,qvalue=0.01,type="hyper")

# get hypo methylated bases
myDifftiles25p.hypo=getMethylDiff(myDifftiles,difference=25,qvalue=0.01,type="hypo")

# get all differentially methylated bases
myDifftiles25p=getMethylDiff(myDifftiles,difference=25,qvalue=0.01)

## -----------------------------------------------------------------------------
diffMethPerChr(myDifftiles,plot=FALSE,qvalue.cutoff=0.01, meth.cutoff=25)

nrow(getData(myDifftiles25p))
nrow(getData(myDifftiles25p.hypo))
nrow(getData(myDifftiles25p.hyper))

# [1] 138
# [1] 71
# [1] 67

ALAS_W_C_lowcovDMRs <- getData(myDifftiles25p)
write.table(ALAS_W_C_lowcovDMRs, "Methylkit_ALAS_W_C_25plowcovDMRs.txt", sep = "\t", quote = FALSE, row.names = FALSE)

########################################
# DMRs for Svalbard

library(methylKit)

# Set the working directory to the folder containing the text files
setwd("/home/celphin/scratch/Dryas/CpG/stranded_CpG_report")

# List all text files in the folder with a specific pattern 
files <- list.files(pattern = "*CpG_report.txt$")

# Define sample ids (replace with your actual sample names)
sample_ids <- sub("_R1_val_1_bismark_bt2_pe\\.deduplicated\\.CpG_report\\.txt$", "", files)

covariate <- sub("^[^_]*_([A-Za-z]{3}).*", "\\1", sample_ids)

# Define treatment conditions (0 = control, 1 = treated)
treatment <- c(rep(0, length(files)/2), rep(1, length(files)/2))

# Extract first letter (C or W)
first_letter <- substr(sample_ids, 1, 1)

# Replace C with 0 and W with 1
treatment <- gsub("C", "0", first_letter)
treatment <- gsub("W", "1", treatment)
treatment <- gsub("L", "0", treatment)

# Coverage header: <chromosome> <start position> <end position> <methylation percentage> <count methylated> <count unmethylated>
# https://rdrr.io/bioc/methylKit/man/methRead-methods.html

d.f <- as.data.frame(cbind(sample_ids, treatment, covariate))

SVAL_d.f <- d.f[which(d.f$covariate=="SVA"),]

load("MethylKit_lowcov_data.RData")

SVAL_myobj <- reorganize(
  myobj_lowCov,
  as.list(SVAL_d.f$sample_ids),
  as.numeric(SVAL_d.f$treatment),
  chunk.size = 1e+06,
  save.db = TRUE,
  suffix = "SVAL",
  dbtype = "tabix"
)

# Print the object to check the result
print(SVAL_myobj)
save(SVAL_myobj, file = "SVAL_MethylKit_lowcov_data.RData")
load("SVAL_MethylKit_lowcov_data.RData")

meth=unite(SVAL_myobj, destrand=FALSE)
save(meth, file = "SVAL_united.RData")

load("SVAL_united.RData")
pdf("SVAL_PCA_PC1_PC2.pdf", width = 10, height = 8)
PCASamples(meth)
dev.off()
pdf("SVAL_Dendrogram_plot.pdf", width = 10, height = 8)
clusterSamples(meth, dist="correlation", method="ward.D2", plot=TRUE)
dev.off()

#------------------------------
# make nicer PCA

allDataPCA <- PCASamples(meth, obj.return = TRUE) #Run a PCA analysis on percent methylation for all samples. methylKit uses prcomp to create the PCA matrix
summary(allDataPCA)

library(vegan)
pdf("SVAL_PCA_PC1_PC2_ordiplot.pdf", width = 10, height = 8)
fig.allDataPCA <- ordiplot(allDataPCA, choices = c(1, 2), type = "none", display = "sites", cex = 0.5, xlim = c(-400, 200), xlab = "", ylab = "", xaxt = "n", yaxt = "n") #Use ordiplot to create base biplot. Do not add any points
points(fig.allDataPCA, "sites", col = c(rep("blue", times = 6), rep("red", times = 6)), pch = c(rep(16, times = 6), rep(17, times = 6)), cex = 3) 

ordiellipse(allDataPCA, SVAL_d.f$treatment, show.groups = "1", col = "red") #Add confidence ellipse around the samples in warming
ordiellipse(allDataPCA, SVAL_d.f$treatment, show.groups = "0", col = "blue") #Add confidence ellipse around the samples in control

axis(side =  1, at = seq(-400, 200, 200), col = "grey80", cex.axis = 1.7) #Add x-axis
mtext(side = 1, text = "PC 1", line = 3, cex = 1.5) #Add x-axis label

axis(side =  2, labels = TRUE, col = "grey80", cex.axis = 1.7) #Add y-axis
mtext(side = 2, text = "PC 2", line = 3, cex = 1.5) #Add y-axis label

mtext(side = 3, line = -5, adj = c(-100,0), text = "    a. All CpG Loci", cex = 1.5)

dev.off()
#--------------------------------
# make regions
lowtiles = tileMethylCounts(SVAL_myobj,win.size=1000,step.size=1000,cov.bases = 10)
head(lowtiles[[1]],3)
save(lowtiles, file = "SVAL_MethylKit_tiles_lowcov.RData")
load("SVAL_MethylKit_tiles_lowcov.RData")

#------------------------
methtiles=unite(lowtiles, destrand=FALSE)
save(methtiles, file = "SVAL_united_lowtiles.RData")
load("SVAL_united_lowtiles.RData")

#---------------------------
# calculate DMRs

myDifftiles=calculateDiffMeth(methtiles)
save(myDifftiles, file = "SVAL_MethDiffRegions.RData")
load("SVAL_MethDiffRegions.RData")

# get hyper methylated bases
myDifftiles10p.hyper=getMethylDiff(myDifftiles,difference=10,qvalue=0.01,type="hyper")

# get hypo methylated bases
myDifftiles10p.hypo=getMethylDiff(myDifftiles,difference=10,qvalue=0.01,type="hypo")

# get all differentially methylated bases
myDifftiles10p=getMethylDiff(myDifftiles,difference=10,qvalue=0.01)

## -----------------------------------------------------------------------------
diffMethPerChr(myDifftiles,plot=FALSE,qvalue.cutoff=0.01, meth.cutoff=10)

nrow(getData(myDifftiles10p))
nrow(getData(myDifftiles10p.hypo))
nrow(getData(myDifftiles10p.hyper))

# [1] 6713
# [1] 2742
# [1] 3971

#--------------------------
# get hyper methylated bases
myDifftiles25p.hyper=getMethylDiff(myDifftiles,difference=25,qvalue=0.01,type="hyper")

# get hypo methylated bases
myDifftiles25p.hypo=getMethylDiff(myDifftiles,difference=25,qvalue=0.01,type="hypo")

# get all differentially methylated bases
myDifftiles25p=getMethylDiff(myDifftiles,difference=25,qvalue=0.01)

## -----------------------------------------------------------------------------
diffMethPerChr(myDifftiles,plot=FALSE,qvalue.cutoff=0.01, meth.cutoff=25)

nrow(getData(myDifftiles25p))
nrow(getData(myDifftiles25p.hypo))
nrow(getData(myDifftiles25p.hyper))

# [1] 254
# [1] 114
# [1] 140

SVAL_W_C_lowcovDMRs <- getData(myDifftiles25p)
write.table(SVAL_W_C_lowcovDMRs, "Methylkit_SVAL_W_C_25plowcovDMRs.txt", sep = "\t", quote = FALSE, row.names = FALSE)

############################################
# DMRs for Sweden

library(methylKit)

# Set the working directory to the folder containing the text files
setwd("/home/celphin/scratch/Dryas/CpG/stranded_CpG_report")

# List all text files in the folder with a specific pattern 
files <- list.files(pattern = "*CpG_report.txt$")

# Define sample ids (replace with your actual sample names)
sample_ids <- sub("_R1_val_1_bismark_bt2_pe\\.deduplicated\\.CpG_report\\.txt$", "", files)

covariate <- sub("^[^_]*_([A-Za-z]{3}).*", "\\1", sample_ids)

# Define treatment conditions (0 = control, 1 = treated)
treatment <- c(rep(0, length(files)/2), rep(1, length(files)/2))

# Extract first letter (C or W)
first_letter <- substr(sample_ids, 1, 1)

# Replace C with 0 and W with 1
treatment <- gsub("C", "0", first_letter)
treatment <- gsub("W", "1", treatment)
treatment <- gsub("L", "0", treatment)

# Coverage header: <chromosome> <start position> <end position> <methylation percentage> <count methylated> <count unmethylated>
# https://rdrr.io/bioc/methylKit/man/methRead-methods.html

d.f <- as.data.frame(cbind(sample_ids, treatment, covariate))

LAT_d.f <- d.f[which(d.f$covariate=="LAT"),]

load("MethylKit_lowcov_data.RData")

LAT_myobj <- reorganize(
  myobj_lowCov,
  as.list(LAT_d.f$sample_ids),
  as.numeric(LAT_d.f$treatment),
  chunk.size = 1e+06,
  save.db = TRUE,
  suffix = "LAT",
  dbtype = "tabix"
)

# Print the object to check the result
print(LAT_myobj)
save(LAT_myobj, file = "LAT_MethylKit_lowcov_data.RData")
load("LAT_MethylKit_lowcov_data.RData")

meth=unite(LAT_myobj, destrand=FALSE)
save(meth, file = "LAT_united.RData")

load("LAT_united.RData")
pdf("LAT_PCA_PC1_PC2.pdf", width = 10, height = 8)
PCASamples(meth)
dev.off()
pdf("LAT_Dendrogram_plot.pdf", width = 10, height = 8)
clusterSamples(meth, dist="correlation", method="ward.D2", plot=TRUE)
dev.off()

#------------------------------
# make nicer PCA

allDataPCA <- PCASamples(meth, obj.return = TRUE) #Run a PCA analysis on percent methylation for all samples. methylKit uses prcomp to create the PCA matrix
summary(allDataPCA)

library(vegan)
pdf("LAT_PCA_PC1_PC2_ordiplot.pdf", width = 10, height = 8)
fig.allDataPCA <- ordiplot(allDataPCA, choices = c(1, 2), type = "none", display = "sites", cex = 0.5, xlim = c(-400, 200), xlab = "", ylab = "", xaxt = "n", yaxt = "n") #Use ordiplot to create base biplot. Do not add any points
points(fig.allDataPCA, "sites", col = c(rep("blue", times = 10), rep("red", times = 10)), pch = c(rep(16, times = 10), rep(17, times = 10)), cex = 3) #Add each sample. Darker samples are ambient, lighter samples are elevated pCO2

ordiellipse(allDataPCA, LAT_d.f$treatment, show.groups = "1", col = "red") #Add confidence ellipse around the samples in elevated pCO2
ordiellipse(allDataPCA, LAT_d.f$treatment, show.groups = "0", col = "blue") #Add confidence ellipse around the samples in ambient pCO2

axis(side =  1, at = seq(-400, 200, 200), col = "grey80", cex.axis = 1.7) #Add x-axis
mtext(side = 1, text = "PC 1", line = 3, cex = 1.5) #Add x-axis label

axis(side =  2, labels = TRUE, col = "grey80", cex.axis = 1.7) #Add y-axis
mtext(side = 2, text = "PC 2", line = 3, cex = 1.5) #Add y-axis label

mtext(side = 3, line = -5, adj = c(-100,0), text = "    a. All CpG Loci", cex = 1.5)

dev.off()
#---------------------------
# make regions
lowtiles = tileMethylCounts(LAT_myobj,win.size=1000,step.size=1000,cov.bases = 10)
head(lowtiles[[1]],3)
save(lowtiles, file = "LAT_MethylKit_tiles_lowcov.RData")
load("LAT_MethylKit_tiles_lowcov.RData")

#------------------------
methtiles=unite(lowtiles, destrand=FALSE)
save(methtiles, file = "LAT_united_lowtiles.RData")
load("LAT_united_lowtiles.RData")

#---------------------------
# calculate DMRs

myDifftiles=calculateDiffMeth(methtiles)
save(myDifftiles, file = "LAT_MethDiffRegions.RData")
load("LAT_MethDiffRegions.RData")

# get hyper methylated bases
myDifftiles10p.hyper=getMethylDiff(myDifftiles,difference=10,qvalue=0.01,type="hyper")

# get hypo methylated bases
myDifftiles10p.hypo=getMethylDiff(myDifftiles,difference=10,qvalue=0.01,type="hypo")

# get all differentially methylated bases
myDifftiles10p=getMethylDiff(myDifftiles,difference=10,qvalue=0.01)

## -----------------------------------------------------------------------------
diffMethPerChr(myDifftiles,plot=FALSE,qvalue.cutoff=0.01, meth.cutoff=10)

nrow(getData(myDifftiles10p))
nrow(getData(myDifftiles10p.hypo))
nrow(getData(myDifftiles10p.hyper))

# [1] 6943
# [1] 3019
# [1] 3924


#--------------------------
# get hyper methylated bases
myDifftiles25p.hyper=getMethylDiff(myDifftiles,difference=25,qvalue=0.01,type="hyper")

# get hypo methylated bases
myDifftiles25p.hypo=getMethylDiff(myDifftiles,difference=25,qvalue=0.01,type="hypo")

# get all differentially methylated bases
myDifftiles25p=getMethylDiff(myDifftiles,difference=25,qvalue=0.01)

## -----------------------------------------------------------------------------
diffMethPerChr(myDifftiles,plot=FALSE,qvalue.cutoff=0.01, meth.cutoff=25)

nrow(getData(myDifftiles25p))
nrow(getData(myDifftiles25p.hypo))
nrow(getData(myDifftiles25p.hyper))

# [1] 327
# [1] 144
# [1] 183

LAT_W_C_lowcovDMRs <- getData(myDifftiles25p)
write.table(LAT_W_C_lowcovDMRs, "Methylkit_LAT_W_C_25plowcovDMRs.txt", sep = "\t", quote = FALSE, row.names = FALSE)

########################################
# DMRs for Nunavut

salloc -c1 --time 2:00:00 --mem 120000m --account def-cronk

R 

library(methylKit)

# Set the working directory to the folder containing the text files
setwd("/home/celphin/scratch/Dryas/CpG/stranded_CpG_report")

# List all text files in the folder with a specific pattern 
files <- list.files(pattern = "*CpG_report.txt$")

# Define sample ids (replace with your actual sample names)
sample_ids <- sub("_R1_val_1_bismark_bt2_pe\\.deduplicated\\.CpG_report\\.txt$", "", files)

covariate <- sub("^[^_]*_([A-Za-z]{3}).*", "\\1", sample_ids)

# Define treatment conditions (0 = control, 1 = treated)
treatment <- c(rep(0, length(files)/2), rep(1, length(files)/2))

# Extract first letter (C or W)
first_letter <- substr(sample_ids, 1, 1)

# Replace C with 0 and W with 1
treatment <- gsub("C", "0", first_letter)
treatment <- gsub("W", "1", treatment)
treatment <- gsub("L", "0", treatment)

# Coverage header: <chromosome> <start position> <end position> <methylation percentage> <count methylated> <count unmethylated>
# https://rdrr.io/bioc/methylKit/man/methRead-methods.html

d.f <- as.data.frame(cbind(sample_ids, treatment, covariate))

Alex_d.f <- d.f[which(d.f$covariate=="CAS" | d.f$covariate=="WIL" |d.f$covariate=="Dry" |d.f$covariate=="MEA" |d.f$covariate=="FER"),]

load("MethylKit_lowcov_data.RData")

Alex_myobj <- reorganize(
  myobj_lowCov,
  as.list(Alex_d.f$sample_ids),
  as.numeric(Alex_d.f$treatment),
  chunk.size = 1e+06,
  save.db = TRUE,
  suffix = "Alex",
  dbtype = "tabix"
)

# Print the object to check the result
print(Alex_myobj)
save(Alex_myobj, file = "Alex_MethylKit_lowcov_data.RData")
load("Alex_MethylKit_lowcov_data.RData")

meth=unite(Alex_myobj, destrand=FALSE)
save(meth, file = "Alex_united.RData")

#-------------------------------

load("Alex_united.RData")
pdf("Alex_PCA_PC1_PC2.pdf", width = 10, height = 8)
PCASamples(meth)
dev.off()
pdf("Alex_PCA_PC3_PC4.pdf", width = 10, height = 8)
PCASamples(meth, comp=c(3,4))
dev.off()

pdf("Alex_Dendrogram_plot.pdf", width = 10, height = 8)
clusterSamples(meth, dist="correlation", method="ward.D2", plot=TRUE)
dev.off()

#------------------------------
# make nicer PCA

allDataPCA <- PCASamples(meth, obj.return = TRUE) #Run a PCA analysis on percent methylation for all samples. methylKit uses prcomp to create the PCA matrix
summary(allDataPCA)

library(vegan)
pdf("Alex_PCA_PC1_PC2_ordiplot.pdf", width = 10, height = 8)
fig.allDataPCA <- ordiplot(allDataPCA, choices = c(1, 2), type = "none", display = "sites", cex = 0.5, xlim = c(-400, 200), xlab = "", ylab = "", xaxt = "n", yaxt = "n") #Use ordiplot to create base biplot. Do not add any points
points(fig.allDataPCA, "sites", col = c(rep("blue", times = 25), rep("red", times = 25)), pch = c(rep(16, times = 25), rep(17, times = 25)), cex = 3) #Add each sample. Darker samples are ambient, lighter samples are elevated pCO2

ordiellipse(allDataPCA, Alex_d.f$treatment, show.groups = "1", col = "red") #Add confidence ellipse around the samples in elevated pCO2
ordiellipse(allDataPCA, Alex_d.f$treatment, show.groups = "0", col = "blue") #Add confidence ellipse around the samples in ambient pCO2

axis(side =  1, at = seq(-400, 200, 200), col = "grey80", cex.axis = 1.7) #Add x-axis
mtext(side = 1, text = "PC 1", line = 3, cex = 1.5) #Add x-axis label

axis(side =  2, labels = TRUE, col = "grey80", cex.axis = 1.7) #Add y-axis
mtext(side = 2, text = "PC 2", line = 3, cex = 1.5) #Add y-axis label

mtext(side = 3, line = -5, adj = c(-100,0), text = "    a. All CpG Loci", cex = 1.5)

dev.off()


#------------------------------
# make regions
lowtiles = tileMethylCounts(Alex_myobj,win.size=1000,step.size=1000,cov.bases = 10)
head(lowtiles[[1]],3)
save(lowtiles, file = "Alex_MethylKit_tiles_lowcov.RData")
load("Alex_MethylKit_tiles_lowcov.RData")

#------------------------
methtiles=unite(lowtiles, destrand=FALSE)
save(methtiles, file = "Alex_united_lowtiles.RData")
load("Alex_united_lowtiles.RData")

#---------------------------
# calculate DMRs

myDifftiles=calculateDiffMeth(methtiles)
save(myDifftiles, file = "Alex_MethDiffRegions.RData")
load("Alex_MethDiffRegions.RData")

# get hyper methylated bases
myDifftiles10p.hyper=getMethylDiff(myDifftiles,difference=10,qvalue=0.01,type="hyper")

# get hypo methylated bases
myDifftiles10p.hypo=getMethylDiff(myDifftiles,difference=10,qvalue=0.01,type="hypo")

# get all differentially methylated bases
myDifftiles10p=getMethylDiff(myDifftiles,difference=10,qvalue=0.01)

## -----------------------------------------------------------------------------
diffMethPerChr(myDifftiles,plot=FALSE,qvalue.cutoff=0.01, meth.cutoff=10)

nrow(getData(myDifftiles10p))
nrow(getData(myDifftiles10p.hypo))
nrow(getData(myDifftiles10p.hyper))

# [1] 677
# [1] 301
# [1] 376


#--------------------------
# get hyper methylated bases
myDifftiles25p.hyper=getMethylDiff(myDifftiles,difference=25,qvalue=0.01,type="hyper")

# get hypo methylated bases
myDifftiles25p.hypo=getMethylDiff(myDifftiles,difference=25,qvalue=0.01,type="hypo")

# get all differentially methylated bases
myDifftiles25p=getMethylDiff(myDifftiles,difference=25,qvalue=0.01)

## -----------------------------------------------------------------------------
diffMethPerChr(myDifftiles,plot=FALSE,qvalue.cutoff=0.01, meth.cutoff=25)

nrow(getData(myDifftiles25p))
nrow(getData(myDifftiles25p.hypo))
nrow(getData(myDifftiles25p.hyper))

# 0

Alex_W_C_lowcovDMRs <- getData(myDifftiles10p)
write.table(Alex_W_C_lowcovDMRs, "Methylkit_Alex_W_C_10plowcovDMRs.txt", sep = "\t", quote = FALSE, row.names = FALSE)

############################################
# intersection of sites

cd /home/celphin/scratch/Dryas/CpG/stranded_CpG_report/DMRs
mv ../*DMRs.txt .

for file in *.txt; do
    base_name=$(basename "$file" .txt)
    awk '{print $1, $2, $3, $NF}' "$file" > "${base_name}.bedGraph"
    sed -i 's/ /\t/g' "${base_name}.bedGraph"
    sed -i 's/chr\tstart\tend\tmeth.diff/track type=bedGraph/g' "${base_name}.bedGraph"
done

# intersect with metilene output

module load StdEnv/2023 bedtools/2.31.0
bedtools intersect -u -a Methylkit_LAT_W_C_25plowcovDMRs.bedGraph \
-b Methylkit_SVAL_W_C_25plowcovDMRs.bedGraph > intersect_SVAL_LAT_W_C_25plowcovDMRs.bedGraph

bedtools intersect -u -a Methylkit_LAT_W_C_25plowcovDMRs.bedGraph \
-b Methylkit_ALAS_W_C_25plowcovDMRs.bedGraph > intersect_ALAS_LAT_W_C_25plowcovDMRs.bedGraph

bedtools intersect -u -a Methylkit_LAT_W_C_25plowcovDMRs.bedGraph \
-b Methylkit_Alex_W_C_10plowcovDMRs.bedGraph > intersect_Alex_LAT_W_C_25plowcovDMRs.bedGraph

bedtools intersect -u -a Methylkit_SVAL_W_C_25plowcovDMRs.bedGraph \
-b Methylkit_Alex_W_C_10plowcovDMRs.bedGraph > intersect_SVAL_Alex_W_C_25plowcovDMRs.bedGraph

bedtools intersect -u -a intersect_ALAS_LAT_W_C_25plowcovDMRs.bedGraph \
-b intersect_SVAL_Alex_W_C_25plowcovDMRs.bedGraph > intersect_ALL_plowcovDMRs.bedGraph

bedtools intersect -u -a intersect_Alex_LAT_W_C_25plowcovDMRs.bedGraph \
-b intersect_SVAL_Alex_W_C_25plowcovDMRs.bedGraph > intersect_Alex_LAT_SVAL_plowcovDMRs.bedGraph


wc -l *.bedGraph
    13 intersect_Alex_LAT_W_C_25plowcovDMRs.bedGraph
    13 intersect_SVAL_Alex_W_C_25plowcovDMRs.bedGraph
	 2 intersect_Alex_LAT_SVAL_plowcovDMRs.bedGraph
	 0 intersect_ALL_plowcovDMRs.bedGraph
     4 intersect_ALAS_LAT_W_C_25plowcovDMRs.bedGraph
     9 intersect_SVAL_LAT_W_C_25plowcovDMRs.bedGraph
   139 Methylkit_ALAS_W_C_25plowcovDMRs.bedGraph
  5086 Methylkit_ALAS_W_C_lowcovDMRs.bedGraph
   678 Methylkit_Alex_W_C_10plowcovDMRs.bedGraph
   328 Methylkit_LAT_W_C_25plowcovDMRs.bedGraph
   255 Methylkit_SVAL_W_C_25plowcovDMRs.bedGraph
    32 Methylkit_Wild_W_C_DMRs.bedGraph
   212 Methylkit_Wild_W_C_lowcovDMRs.bedGraph


##########################################
# get feature overlap with SNPEff and TEs

tmux new-session -s snpEff
tmux attach-session -t snpEff

mkdir /home/celphin/scratch/Dryas/snpEff/Dryas_DMRs/methylkit/
cp /home/celphin/scratch/Dryas/CpG/stranded_CpG_report/DMRs/*bedGraph /home/celphin/scratch/Dryas/snpEff/Dryas_DMRs/methylkit/
cd /home/celphin/scratch/Dryas/snpEff/Dryas_DMRs/methylkit/

FILES=$(ls /home/celphin/scratch/Dryas/snpEff/Dryas_DMRs/methylkit/*.bedGraph)
echo ${FILES}

for file in ${FILES}; do
    # Run the command on the file and create the .bed file
    grep -v track "$file" | awk '{print $1 "\t" $2 "\t" $3}' > "${file%.bedGraph}.bed"
    echo "Processed $file"
done

# run SNPEff
module load StdEnv/2023 java/21.0.1
cd /home/celphin/scratch/Dryas/snpEff/Dryas_DMRs/sites/
FILES=$(ls /home/celphin/scratch/Dryas/snpEff/Dryas_DMRs/methylkit/*.bed)
echo ${FILES}

cd /home/celphin/scratch/Dryas/snpEff

for file in ${FILES}; do
java -Xmx8g -jar snpEff.jar -i bed OldDoct "$file" > "$file".out
echo "Processed $file"
done


# extract the immediate feature types and count them
cd /home/celphin/scratch/Dryas/snpEff/Dryas_DMRs/methylkit/
FILES=$(ls *.out)
echo ${FILES}

for file in ${FILES}; do
awk -F'[;:]' '{print $1 $2}' "$file" > "$file"1
echo 
echo "$file Upstream"
grep "Upstream" "$file"1 | wc -l 
echo "$file Downstream"
grep "Downstream" "$file"1 | wc -l 
echo "$file Gene"
grep "Gene" "$file"1 | wc -l 
echo "$file Intergenic"
grep "Intergenic" "$file"1 | wc -l 
echo "$file Intron"
grep "Intron" "$file"1 | wc -l 
echo "$file Exon"
grep "Exon" "$file"1 | wc -l 
done

#--------------------------

# Look at transposons
# https://github.com/RobertsLab/project-gigas-oa-meth/blob/master/code/10-Genomic-Location-of-DML.ipynb 

module load StdEnv/2023 bedtools/2.31.0

cd /home/celphin/scratch/Dryas/snpEff/Dryas_DMRs/methylkit/
FILES=$(ls *.bed)
echo ${FILES}

for file in ${FILES}; do
intersectBed \
-u \
-a ${file} \
-b /home/celphin/scratch/Dryas/snpEff/data/OldDoct/OldDoct.genome.fa.mod.EDTA.TEanno.gff3 \
> ${file}-All-TE.bed

#head ${file}-All-TE.bed
wc -l ${file}-All-TE.bed

done

wc -l *-All-TE.bed


#-------------------------
# could subet the TEs by type and rerun for each type

awk '{print $3}' /home/celphin/scratch/Dryas/snpEff/data/OldDoct/OldDoct.genome.fa.mod.EDTA.TEanno.gff3 | sort | uniq

# CACTA_TIR_transposon
# Copia_LTR_retrotransposon
# Gypsy_LTR_retrotransposon
# hAT_TIR_transposon
# helitron
# identity
# long_terminal_repeat
# LTR_retrotransposon
# Mutator_TIR_transposon
# Nov
# PIF_Harbinger_TIR_transposon
# repeat_region
# sequence_ontology
# site
# target_site_duplication
# Tc1_Mariner_TIR_transposon

cd /home/celphin/scratch/Dryas/snpEff/Dryas_DMRs/
# Extract unique values from the third column (or any other column)
unique_values=$(awk '{print $3}' /home/celphin/scratch/Dryas/snpEff/data/OldDoct/OldDoct.genome.fa.mod.EDTA.TEanno.gff3 | sort | uniq)

# Loop through each unique value and extract lines to a new file
for value in $unique_values; do
    grep "$value" /home/celphin/scratch/Dryas/snpEff/data/OldDoct/OldDoct.genome.fa.mod.EDTA.TEanno.gff3 > "${value}.gff3"
    echo "Extracted lines for $value into ${value}.gff3"
done

mkdir gff3
mv *.gff3 gff3/


#-----------------
module load StdEnv/2023 bedtools/2.31.0
unique_values=$(awk '{print $3}' /home/celphin/scratch/Dryas/snpEff/data/OldDoct/OldDoct.genome.fa.mod.EDTA.TEanno.gff3 | sort | uniq)

cd /home/celphin/scratch/Dryas/snpEff/Dryas_DMRs/
FILES=$(ls ./methylkit/*-All-TE.bed)
echo ${FILES}

for value in $unique_values; do
for file in ${FILES}; do
intersectBed \
-u \
-a ${file} \
-b ./gff3/${value}.gff3 \
> ${file}-${value}.bed

wc -l ${file}-${value}.bed

done
done


#################################
# move to other reference with CrossMap

# Map to RagTag reference
# to run
cd /home/celphin/scratch/Dryas/CrossMap/file_for_converting

# copy data to this folder
cp /home/celphin/scratch/Dryas/snpEff/Dryas_DMRs/methylkit/*.bedGraph ./methylkit/

# Bedgraphs CpG per site_treatment
CrossMap bed ../Dryas_RagTag_Old.chain C_ALAS.bedGraph RagTag_C_ALAS.bedGraph
CrossMap bed ../Dryas_RagTag_Old.chain C_CASS.bedGraph RagTag_C_CASS.bedGraph

##################################
# GO enrichment










#######################################
# Randomize methylkit

cd /home/celphin/scratch/Dryas/CpG/

tmux new-session -s methylkit
tmux attach-session -t methylkit

module load StdEnv/2023
module load r/4.4.0

R

## ----------------------------------------------------------------
# DMRs for Alaska

library(methylKit)

# Set the working directory to the folder containing the text files
setwd("/home/celphin/scratch/Dryas/CpG/stranded_CpG_report")

# List all text files in the folder with a specific pattern 
files <- list.files(pattern = "*CpG_report.txt$")

# Define sample ids (replace with your actual sample names)
sample_ids <- sub("_R1_val_1_bismark_bt2_pe\\.deduplicated\\.CpG_report\\.txt$", "", files)

covariate <- sub("^[^_]*_([A-Za-z]{3}).*", "\\1", sample_ids)

d.f <- as.data.frame(cbind(sample_ids, covariate))

ALAS_rand_d.f <- d.f[which(d.f$covariate=="ALA"),]

# Extract first letter (C or W)
first_letter <- substr(ALAS_rand_d.f$sample_ids, 1, 1)

# Replace C with 0 and W with 1
treatment <- gsub("C", "0", first_letter)
treatment <- gsub("W", "1", treatment)

# randomly rearrange treatment
treatment0 <- sample(treatment)

ALAS_rand_d.f$treatment0 <- treatment0 

load("ALAS_united_lowtiles.RData")

methtiles_rand <- reorganize(
  methtiles,
  as.list(ALAS_rand_d.f$sample_ids),
  as.numeric(ALAS_rand_d.f$treatment0),
  chunk.size = 1e+06,
  save.db = TRUE,
  suffix = "ALAS_rand",
  dbtype = "tabix"
)


#---------------------------
# calculate DMRs

myDifftiles=calculateDiffMeth(methtiles_rand)
save(myDifftiles, file = "ALAS_rand_MethDiffRegions.RData")
load("ALAS_rand_MethDiffRegions.RData")

# get hyper methylated bases
myDifftiles10p.hyper=getMethylDiff(myDifftiles,difference=10,qvalue=0.01,type="hyper")

# get hypo methylated bases
myDifftiles10p.hypo=getMethylDiff(myDifftiles,difference=10,qvalue=0.01,type="hypo")

# get all differentially methylated bases
myDifftiles10p=getMethylDiff(myDifftiles,difference=10,qvalue=0.01)

## -----------------------------------------------------------------------------
diffMethPerChr(myDifftiles,plot=FALSE,qvalue.cutoff=0.01, meth.cutoff=10)

nrow(getData(myDifftiles10p))
nrow(getData(myDifftiles10p.hypo))
nrow(getData(myDifftiles10p.hyper))

# [1] 5085
# [1] 2609
# [1] 2476

# rand
# [1] 3405
# [1] 1488
# [1] 1917

#--------------------------
# get hyper methylated bases
myDifftiles25p.hyper=getMethylDiff(myDifftiles,difference=25,qvalue=0.01,type="hyper")

# get hypo methylated bases
myDifftiles25p.hypo=getMethylDiff(myDifftiles,difference=25,qvalue=0.01,type="hypo")

# get all differentially methylated bases
myDifftiles25p=getMethylDiff(myDifftiles,difference=25,qvalue=0.01)

## -----------------------------------------------------------------------------
diffMethPerChr(myDifftiles,plot=FALSE,qvalue.cutoff=0.01, meth.cutoff=25)

nrow(getData(myDifftiles25p))
nrow(getData(myDifftiles25p.hypo))
nrow(getData(myDifftiles25p.hyper))

# [1] 138
# [1] 71
# [1] 67

# [1] 56
# [1] 19
# [1] 37


ALAS_rand_W_C_lowcovDMRs <- getData(myDifftiles25p)
write.table(ALAS_rand_W_C_lowcovDMRs, "Methylkit_ALAS_rand_W_C_25plowcovDMRs.txt", sep = "\t", quote = FALSE, row.names = FALSE)

########################################
# DMRs for Svalbard

library(methylKit)

# Set the working directory to the folder containing the text files
setwd("/home/celphin/scratch/Dryas/CpG/stranded_CpG_report")

# List all text files in the folder with a specific pattern 
files <- list.files(pattern = "*CpG_report.txt$")

# Define sample ids (replace with your actual sample names)
sample_ids <- sub("_R1_val_1_bismark_bt2_pe\\.deduplicated\\.CpG_report\\.txt$", "", files)

covariate <- sub("^[^_]*_([A-Za-z]{3}).*", "\\1", sample_ids)

d.f <- as.data.frame(cbind(sample_ids, covariate))

SVAL_rand_d.f <- d.f[which(d.f$covariate=="SVA"),]

# Extract first letter (C or W)
first_letter <- substr(SVAL_rand_d.f$sample_ids, 1, 1)

# Replace C with 0 and W with 1
treatment <- gsub("C", "0", first_letter)
treatment <- gsub("W", "1", treatment)

# randomly rearrange treatment
treatment0 <- sample(treatment)

SVAL_rand_d.f$treatment0 <- treatment0 

load("SVAL_united_lowtiles.RData")

methtiles_rand <- reorganize(
  methtiles,
  as.list(SVAL_rand_d.f$sample_ids),
  as.numeric(SVAL_rand_d.f$treatment0),
  chunk.size = 1e+06,
  save.db = TRUE,
  suffix = "SVAL_rand",
  dbtype = "tabix"
)

#---------------------------
# calculate DMRs

myDifftiles=calculateDiffMeth(methtiles_rand)
save(myDifftiles, file = "SVAL_rand_MethDiffRegions.RData")
load("SVAL_rand_MethDiffRegions.RData")

# get hyper methylated bases
myDifftiles10p.hyper=getMethylDiff(myDifftiles,difference=10,qvalue=0.01,type="hyper")

# get hypo methylated bases
myDifftiles10p.hypo=getMethylDiff(myDifftiles,difference=10,qvalue=0.01,type="hypo")

# get all differentially methylated bases
myDifftiles10p=getMethylDiff(myDifftiles,difference=10,qvalue=0.01)

## -----------------------------------------------------------------------------
diffMethPerChr(myDifftiles,plot=FALSE,qvalue.cutoff=0.01, meth.cutoff=10)

nrow(getData(myDifftiles10p))
nrow(getData(myDifftiles10p.hypo))
nrow(getData(myDifftiles10p.hyper))

# [1] 6588
# [1] 2696
# [1] 3892

#--------------------------
# get hyper methylated bases
myDifftiles25p.hyper=getMethylDiff(myDifftiles,difference=25,qvalue=0.01,type="hyper")

# get hypo methylated bases
myDifftiles25p.hypo=getMethylDiff(myDifftiles,difference=25,qvalue=0.01,type="hypo")

# get all differentially methylated bases
myDifftiles25p=getMethylDiff(myDifftiles,difference=25,qvalue=0.01)

## -----------------------------------------------------------------------------
diffMethPerChr(myDifftiles,plot=FALSE,qvalue.cutoff=0.01, meth.cutoff=25)

nrow(getData(myDifftiles25p))
nrow(getData(myDifftiles25p.hypo))
nrow(getData(myDifftiles25p.hyper))

# [1] 239
# [1] 90
# [1] 149


SVAL_rand_W_C_lowcovDMRs <- getData(myDifftiles25p)
write.table(SVAL_rand_W_C_lowcovDMRs, "Methylkit_SVAL_rand_W_C_25plowcovDMRs.txt", sep = "\t", quote = FALSE, row.names = FALSE)

#######################################
# Sweden

library(methylKit)

# Set the working directory to the folder containing the text files
setwd("/home/celphin/scratch/Dryas/CpG/stranded_CpG_report")

# List all text files in the folder with a specific pattern 
files <- list.files(pattern = "*CpG_report.txt$")

# Define sample ids (replace with your actual sample names)
sample_ids <- sub("_R1_val_1_bismark_bt2_pe\\.deduplicated\\.CpG_report\\.txt$", "", files)

covariate <- sub("^[^_]*_([A-Za-z]{3}).*", "\\1", sample_ids)

d.f <- as.data.frame(cbind(sample_ids, covariate))

LAT_rand_d.f <- d.f[which(d.f$covariate=="LAT"),]

# Extract first letter (C or W)
first_letter <- substr(LAT_rand_d.f$sample_ids, 1, 1)

# Replace C with 0 and W with 1
treatment <- gsub("C", "0", first_letter)
treatment <- gsub("W", "1", treatment)

# randomly rearrange treatment
treatment0 <- sample(treatment)

LAT_rand_d.f$treatment0 <- treatment0 

load("LAT_united_lowtiles.RData")

methtiles_rand <- reorganize(
  methtiles,
  as.list(LAT_rand_d.f$sample_ids),
  as.numeric(LAT_rand_d.f$treatment0),
  chunk.size = 1e+06,
  save.db = TRUE,
  suffix = "LAT_rand",
  dbtype = "tabix"
)


#---------------------------
# calculate DMRs

myDifftiles=calculateDiffMeth(methtiles_rand)
save(myDifftiles, file = "LAT_rand_MethDiffRegions.RData")
load("LAT_rand_MethDiffRegions.RData")

# get hyper methylated bases
myDifftiles10p.hyper=getMethylDiff(myDifftiles,difference=10,qvalue=0.01,type="hyper")

# get hypo methylated bases
myDifftiles10p.hypo=getMethylDiff(myDifftiles,difference=10,qvalue=0.01,type="hypo")

# get all differentially methylated bases
myDifftiles10p=getMethylDiff(myDifftiles,difference=10,qvalue=0.01)

## -----------------------------------------------------------------------------
diffMethPerChr(myDifftiles,plot=FALSE,qvalue.cutoff=0.01, meth.cutoff=10)

nrow(getData(myDifftiles10p))
nrow(getData(myDifftiles10p.hypo))
nrow(getData(myDifftiles10p.hyper))

# [1] 1597
# [1] 585
# [1] 1012

#--------------------------
# get hyper methylated bases
myDifftiles25p.hyper=getMethylDiff(myDifftiles,difference=25,qvalue=0.01,type="hyper")

# get hypo methylated bases
myDifftiles25p.hypo=getMethylDiff(myDifftiles,difference=25,qvalue=0.01,type="hypo")

# get all differentially methylated bases
myDifftiles25p=getMethylDiff(myDifftiles,difference=25,qvalue=0.01)

## -----------------------------------------------------------------------------
diffMethPerChr(myDifftiles,plot=FALSE,qvalue.cutoff=0.01, meth.cutoff=25)

nrow(getData(myDifftiles25p))
nrow(getData(myDifftiles25p.hypo))
nrow(getData(myDifftiles25p.hyper))

# [1] 20
# [1] 5
# [1] 15

LAT_rand_W_C_lowcovDMRs <- getData(myDifftiles25p)
write.table(LAT_rand_W_C_lowcovDMRs, "Methylkit_LAT_rand_W_C_25plowcovDMRs.txt", sep = "\t", quote = FALSE, row.names = FALSE)

####################################
#Nunavut

library(methylKit)

# Set the working directory to the folder containing the text files
setwd("/home/celphin/scratch/Dryas/CpG/stranded_CpG_report")

# List all text files in the folder with a specific pattern 
files <- list.files(pattern = "*CpG_report.txt$")

# Define sample ids (replace with your actual sample names)
sample_ids <- sub("_R1_val_1_bismark_bt2_pe\\.deduplicated\\.CpG_report\\.txt$", "", files)

covariate <- sub("^[^_]*_([A-Za-z]{3}).*", "\\1", sample_ids)

d.f <- as.data.frame(cbind(sample_ids, covariate))

Alex_rand_d.f <- d.f[which(d.f$covariate=="CAS" | d.f$covariate=="WIL" |d.f$covariate=="Dry" |d.f$covariate=="MEA" |d.f$covariate=="FER"),]

# Extract first letter (C or W)
first_letter <- substr(Alex_rand_d.f$sample_ids, 1, 1)

# Replace C with 0 and W with 1
treatment <- gsub("C", "0", first_letter)
treatment <- gsub("W", "1", treatment)

# randomly rearrange treatment
treatment0 <- sample(treatment)

Alex_rand_d.f$treatment0 <- treatment0 

load("Alex_united_lowtiles.RData")

methtiles_rand <- reorganize(
  methtiles,
  as.list(Alex_rand_d.f$sample_ids),
  as.numeric(Alex_rand_d.f$treatment0),
  chunk.size = 1e+06,
  save.db = TRUE,
  suffix = "Alex_rand",
  dbtype = "tabix"
)


#---------------------------
# calculate DMRs

myDifftiles=calculateDiffMeth(methtiles_rand)
save(myDifftiles, file = "Alex_rand_MethDiffRegions.RData")
load("Alex_rand_MethDiffRegions.RData")

# get hyper methylated bases
myDifftiles10p.hyper=getMethylDiff(myDifftiles,difference=10,qvalue=0.01,type="hyper")

# get hypo methylated bases
myDifftiles10p.hypo=getMethylDiff(myDifftiles,difference=10,qvalue=0.01,type="hypo")

# get all differentially methylated bases
myDifftiles10p=getMethylDiff(myDifftiles,difference=10,qvalue=0.01)

## -----------------------------------------------------------------------------
diffMethPerChr(myDifftiles,plot=FALSE,qvalue.cutoff=0.01, meth.cutoff=10)

nrow(getData(myDifftiles10p))
nrow(getData(myDifftiles10p.hypo))
nrow(getData(myDifftiles10p.hyper))

# [1] 707
# [1] 448
# [1] 259

#--------------------------
# get hyper methylated bases
myDifftiles25p.hyper=getMethylDiff(myDifftiles,difference=25,qvalue=0.01,type="hyper")

# get hypo methylated bases
myDifftiles25p.hypo=getMethylDiff(myDifftiles,difference=25,qvalue=0.01,type="hypo")

# get all differentially methylated bases
myDifftiles25p=getMethylDiff(myDifftiles,difference=25,qvalue=0.01)

## -----------------------------------------------------------------------------
diffMethPerChr(myDifftiles,plot=FALSE,qvalue.cutoff=0.01, meth.cutoff=25)

nrow(getData(myDifftiles25p))
nrow(getData(myDifftiles25p.hypo))
nrow(getData(myDifftiles25p.hyper))

# 3
# 3
# 0

Alex_rand_W_C_lowcovDMRs <- getData(myDifftiles25p)
write.table(Alex_rand_W_C_lowcovDMRs, "Methylkit_Alex_rand_W_C_25plowcovDMRs.txt", sep = "\t", quote = FALSE, row.names = FALSE)

########################################