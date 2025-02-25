###########################
# methylKit
# https://github.com/al2na/methylKit
# https://bioconductor.org/packages/release/bioc/vignettes/methylKit/inst/doc/methylKit.html#4_Annotating_differentially_methylated_bases_or_regions
# Dec 2024
##############################

# data

cd /lustre04/scratch/celphin/Dryas_large_folders/CpG

tmux new-session -s methwarm4
tmux attach-session -t methwarm4

#-----------------------------
# ideal input data - methylation call file
##         chrBase   chr    base strand coverage freqC  freqT
## 1 chr21.9764539 chr21 9764539      R       12 25.00  75.00
## 2 chr21.9764513 chr21 9764513      R       12  0.00 100.00
## 3 chr21.9820622 chr21 9820622      F       13  0.00 100.00
## 4 chr21.9837545 chr21 9837545      F       11  0.00 100.00
## 5 chr21.9849022 chr21 9849022      F      124 72.58  27.42

# https://felixkrueger.github.io/Bismark/bismark/methylation_extraction/
cd /lustre04/scratch/celphin/Dryas_large_folders/CpG/stranded_CpG_report/

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

salloc -c1 --time 3:00:00 --mem 120000m --account def-henryg

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
# https://bioconductor.org/packages/release/bioc/vignettes/methylKit/inst/doc/methylKit.html#36_Finding_differentially_methylated_bases_or_regions


###############################

library(methylKit)

# Set the working directory to the folder containing the text files
setwd("/lustre04/scratch/celphin/Dryas_large_folders/CpG/Wild_stranded_CpG_report")

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
#print(myobj_lowCov)
save(myobj_lowCov, file = "MethylKit_lowcov_data.RData")
load("MethylKit_lowcov_data.RData")

#------------------------
#make regions 300bp
lowtiles = tileMethylCounts(myobj_lowCov,win.size=300,step.size=300,cov.bases = 10)
head(lowtiles[[1]],3)
save(lowtiles, file = "MethylKit_tiles300_lowcov.RData")
load("MethylKit_tiles300_lowcov.RData")


methtiles=unite(lowtiles, destrand=FALSE)
save(methtiles, file = "united_lowtiles300.RData")
load("united_lowtiles300.RData")

#------------------------
# make regions 1000bp
lowtiles = tileMethylCounts(myobj_lowCov,win.size=1000,step.size=1000,cov.bases = 10)
head(lowtiles[[1]],3)
save(lowtiles, file = "MethylKit_tiles_lowcov.RData")
load("MethylKit_tiles_lowcov.RData")

methtiles=unite(lowtiles, destrand=FALSE)
save(methtiles, file = "united_lowtiles.RData")
load("united_lowtiles.RData")

#---------------------------
# calculate DMRs

myDifftiles=calculateDiffMeth(methtiles, covariates=d.f$covariates, overdispersion="MN", test="Chisq")

#---------------------------
# get hyper methylated bases
myDifftiles10p.hyper=getMethylDiff(myDifftiles,difference=10,qvalue=0.01,type="hyper")

# get hypo methylated bases
myDifftiles10p.hypo=getMethylDiff(myDifftiles,difference=10,qvalue=0.01,type="hypo")

# get all differentially methylated bases
myDifftiles10p=getMethylDiff(myDifftiles,difference=10,qvalue=0.01)


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

# with covariates, overdisp, chisq  - none

Wild_W_C_DMRs <- getData(myDifftiles10p)
write.table(Wild_W_C_DMRs, "Methylkit_Wild_W_C_10_DMRs.txt", sep = "\t", quote = FALSE, row.names = FALSE)

#--------------------------------
# get hyper methylated bases
myDifftiles25p.hyper=getMethylDiff(myDifftiles,difference=25,qvalue=0.01,type="hyper")

# get hypo methylated bases
myDifftiles25p.hypo=getMethylDiff(myDifftiles,difference=25,qvalue=0.01,type="hypo")

# get all differentially methylated bases
myDifftiles25p=getMethylDiff(myDifftiles,difference=25,qvalue=0.01)


diffMethPerChr(myDifftiles,plot=FALSE,qvalue.cutoff=0.01, meth.cutoff=25)

nrow(getData(myDifftiles25p))
nrow(getData(myDifftiles25p.hypo))
nrow(getData(myDifftiles25p.hyper))

Wild_W_C_DMRs <- getData(myDifftiles25p)
write.table(Wild_W_C_DMRs, "Methylkit_Wild_W_C_25_DMRs.txt", sep = "\t", quote = FALSE, row.names = FALSE)


#####################################
# try running again for each site independently - set meth diff to 10 or 25%

cd /lustre04/scratch/celphin/Dryas_large_folders/CpG/Wild_stranded_CpG_report

tmux new-session -s methylkit
tmux attach-session -t methylkit

module load StdEnv/2023
module load r/4.4.0

R

## ----------------------------------------------------------------
# DMRs for Alaska

library(methylKit)

# Set the working directory to the folder containing the text files
setwd("/lustre04/scratch/celphin/Dryas_large_folders/CpG/Wild_stranded_CpG_report")

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

#-----------------------
#plotting
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

myDifftiles=calculateDiffMeth(methtiles, overdispersion="MN", test="Chisq")

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

# overdip and chisq
# [1] 158
# [1] 80
# [1] 78

ALAS_W_C_lowcovDMRs <- getData(myDifftiles10p)
write.table(ALAS_W_C_lowcovDMRs, "Methylkit_ALAS_W_C_10plowcovDMRs.txt", sep = "\t", quote = FALSE, row.names = FALSE)

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

# overdisp and chisq
# [1] 20
# [1] 9
# [1] 11


ALAS_W_C_lowcovDMRs <- getData(myDifftiles25p)
write.table(ALAS_W_C_lowcovDMRs, "Methylkit_ALAS_W_C_25plowcovDMRs.txt", sep = "\t", quote = FALSE, row.names = FALSE)

#############################################
# With over dispersion

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

# overdip and chisq
# [1] 158
# [1] 80
# [1] 78

ALAS_W_C_lowcovDMRs <- getData(myDifftiles10p)
write.table(ALAS_W_C_lowcovDMRs, "Methylkit_ALAS_W_C_10_overdisp_DMRs.txt", sep = "\t", quote = FALSE, row.names = FALSE)

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

# overdisp and chisq
# [1] 20
# [1] 9
# [1] 11


ALAS_W_C_lowcovDMRs <- getData(myDifftiles25p)
write.table(ALAS_W_C_lowcovDMRs, "Methylkit_ALAS_W_C_25_overdisp_DMRs.txt", sep = "\t", quote = FALSE, row.names = FALSE)

########################################
# DMRs for Svalbard

library(methylKit)

# Set the working directory to the folder containing the text files
setwd("/lustre04/scratch/celphin/Dryas_large_folders/CpG/Wild_stranded_CpG_report")
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

myDifftiles=calculateDiffMeth(methtiles, overdispersion="MN", test="Chisq")
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

# less than random
# [1] 92
# [1] 37
# [1] 55

SVAL_W_C_lowcovDMRs <- getData(myDifftiles10p)
write.table(SVAL_W_C_lowcovDMRs, "Methylkit_SVAL_W_C_10plowcovDMRs.txt", sep = "\t", quote = FALSE, row.names = FALSE)


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

# less than random
# [1] 11
# [1] 6
# [1] 5


SVAL_W_C_lowcovDMRs <- getData(myDifftiles25p)
write.table(SVAL_W_C_lowcovDMRs, "Methylkit_SVAL_W_C_25plowcovDMRs.txt", sep = "\t", quote = FALSE, row.names = FALSE)

####------------------------------------

# calculate overdisp DMRs

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

# less than random
# [1] 92
# [1] 37
# [1] 55

SVAL_W_C_lowcovDMRs <- getData(myDifftiles10p)
write.table(SVAL_W_C_lowcovDMRs, "Methylkit_SVAL_W_C_10_overdisp_DMRs.txt", sep = "\t", quote = FALSE, row.names = FALSE)


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

# less than random
# [1] 11
# [1] 6
# [1] 5


SVAL_W_C_lowcovDMRs <- getData(myDifftiles25p)
write.table(SVAL_W_C_lowcovDMRs, "Methylkit_SVAL_W_C_25_overdisp_DMRs.txt", sep = "\t", quote = FALSE, row.names = FALSE)


############################################
# DMRs for Sweden

library(methylKit)

# Set the working directory to the folder containing the text files
setwd("/lustre04/scratch/celphin/Dryas_large_folders/CpG/Wild_stranded_CpG_report")
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
#print(LAT_myobj)
save(LAT_myobj, file = "LAT_MethylKit_lowcov_data.RData")
load("LAT_MethylKit_lowcov_data.RData")

meth=unite(LAT_myobj, destrand=FALSE)
save(meth, file = "LAT_united.RData")

#------------------
# plotting
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

myDifftiles=calculateDiffMeth(methtiles, overdispersion="MN", test="Chisq")
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

# overdisp and chisq
# [1] 964
# [1] 325
# [1] 639

LAT_W_C_lowcovDMRs <- getData(myDifftiles10p)
write.table(LAT_W_C_lowcovDMRs, "Methylkit_LAT_W_C_10plowcovDMRs.txt", sep = "\t", quote = FALSE, row.names = FALSE)

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

# overdisp and chisq
# [1] 114
# [1] 49
# [1] 65


LAT_W_C_lowcovDMRs <- getData(myDifftiles25p)
write.table(LAT_W_C_lowcovDMRs, "Methylkit_LAT_W_C_25plowcovDMRs.txt", sep = "\t", quote = FALSE, row.names = FALSE)
###------------------------------

# calculate overdisp DMRs

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

# overdisp and chisq
# [1] 964
# [1] 325
# [1] 639

LAT_W_C_lowcovDMRs <- getData(myDifftiles10p)
write.table(LAT_W_C_lowcovDMRs, "Methylkit_LAT_W_C_10_overdisp_DMRs.txt", sep = "\t", quote = FALSE, row.names = FALSE)

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

# overdisp and chisq
# [1] 114
# [1] 49
# [1] 65


LAT_W_C_lowcovDMRs <- getData(myDifftiles25p)
write.table(LAT_W_C_lowcovDMRs, "Methylkit_LAT_W_C_25_overdisp_DMRs.txt", sep = "\t", quote = FALSE, row.names = FALSE)

########################################
# DMRs for Nunavut

salloc -c1 --time 2:00:00 --mem 120000m --account def-cronk

R 

library(methylKit)

# Set the working directory to the folder containing the text files
setwd("/lustre04/scratch/celphin/Dryas_large_folders/CpG/Wild_stranded_CpG_report")

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
# plotting
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

myDifftiles=calculateDiffMeth(methtiles, covariates=Alex_d.f$covariates, overdispersion="MN",test="Chisq")
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

# 0 with overdispersion, Chisq and covariates

Alex_W_C_lowcovDMRs <- getData(myDifftiles10p)
write.table(Alex_W_C_lowcovDMRs, "Methylkit_Alex_W_C_10plowcovDMRs.txt", sep = "\t", quote = FALSE, row.names = FALSE)

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

Alex_W_C_lowcovDMRs <- getData(myDifftiles25p)
write.table(Alex_W_C_lowcovDMRs, "Methylkit_Alex_W_C_25plowcovDMRs.txt", sep = "\t", quote = FALSE, row.names = FALSE)


#################################

# calculate overdisp DMRs

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

# 0 with overdispersion, Chisq and covariates

Alex_W_C_lowcovDMRs <- getData(myDifftiles10p)
write.table(Alex_W_C_lowcovDMRs, "Methylkit_Alex_W_C_10_overdisp_DMRs.txt", sep = "\t", quote = FALSE, row.names = FALSE)

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

Alex_W_C_lowcovDMRs <- getData(myDifftiles25p)
write.table(Alex_W_C_lowcovDMRs, "Methylkit_Alex_W_C_25_overdisp_DMRs.txt", sep = "\t", quote = FALSE, row.names = FALSE)

##################################
# Compare subsites methylation to eachother


salloc -c1 --time 2:00:00 --mem 120000m --account def-cronk

R 

library(methylKit)

# Set the working directory to the folder containing the text files
setwd("/lustre04/scratch/celphin/Dryas_large_folders/CpG/Wild_stranded_CpG_report")

# List all text files in the folder with a specific pattern 
files <- list.files(pattern = "*CpG_report.txt$")

# Define sample ids (replace with your actual sample names)
sample_ids <- sub("_R1_val_1_bismark_bt2_pe\\.deduplicated\\.CpG_report\\.txt$", "", files)

covariate <- sub("^[^_]*_([A-Za-z]{3}).*", "\\1", sample_ids)

# Define treatment conditions (0 = control, 1 = treated)
treatment <- c(rep(0, length(files)/2), rep(1, length(files)/2))

# Replace Dry sites with 0 and Wet sites with 1
treatment <- gsub("CAS", "0", covariate)
treatment <- gsub("Dry", "0", treatment)
treatment <- gsub("WIL", "0", treatment)
treatment <- gsub("FER", "1", treatment)
treatment <- gsub("MEA", "1", treatment)

# Coverage header: <chromosome> <start position> <end position> <methylation percentage> <count methylated> <count unmethylated>
# https://rdrr.io/bioc/methylKit/man/methRead-methods.html

d.f <- as.data.frame(cbind(sample_ids, treatment, covariate))

Alex_d.f <- d.f[c(which(d.f$covariate=="CAS"), which(d.f$covariate=="WIL"), which(d.f$covariate=="Dry"), which(d.f$covariate=="MEA"), which(d.f$covariate=="FER")),]

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
#print(Alex_myobj)
save(Alex_myobj, file = "AlexWet-Dry_MethylKit_lowcov_data.RData")
load("AlexWet-Dry_MethylKit_lowcov_data.RData")

meth=unite(Alex_myobj, destrand=FALSE)
save(meth, file = "AlexWet-Dry_united.RData")

#-------------------------------
# plotting
load("Alex_united.RData")
pdf("AlexWet-Dry_PCA_PC1_PC2.pdf", width = 10, height = 8)
PCASamples(meth)
dev.off()
pdf("AlexWet-Dry_PCA_PC3_PC4.pdf", width = 10, height = 8)
PCASamples(meth, comp=c(3,4))
dev.off()

pdf("AlexWet-Dry_Dendrogram_plot.pdf", width = 10, height = 8)
clusterSamples(meth, dist="correlation", method="ward.D2", plot=TRUE)
dev.off()

#------------------------------
# make nicer PCA

allDataPCA <- PCASamples(meth, obj.return = TRUE) #Run a PCA analysis on percent methylation for all samples. methylKit uses prcomp to create the PCA matrix
summary(allDataPCA)

library(vegan)
pdf("AlexWet-Dry_PCA_PC1_PC2_ordiplot.pdf", width = 10, height = 8)
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
save(lowtiles, file = "AlexWet-Dry_MethylKit_tiles_lowcov.RData")
load("AlexWet-Dry_MethylKit_tiles_lowcov.RData")

#------------------------
methtiles=unite(lowtiles, destrand=FALSE)
save(methtiles, file = "AlexWet-Dry_united_lowtiles.RData")
load("AlexWet-Dry_united_lowtiles.RData")

#---------------------------
# calculate DMRs

myDifftiles=calculateDiffMeth(methtiles, covariates=Alex_d.f$covariates, overdispersion="MN",test="Chisq")
save(myDifftiles, file = "AlexWet-Dry_MethDiffRegions.RData")
load("AlexWet-Dry_MethDiffRegions.RData")

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

#1
#1
#0

Alex_W_C_lowcovDMRs <- getData(myDifftiles10p)
write.table(Alex_W_C_lowcovDMRs, "Methylkit_AlexWet-Dry_10plowcovDMRs.txt", sep = "\t", quote = FALSE, row.names = FALSE)

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



Alex_W_C_lowcovDMRs <- getData(myDifftiles25p)
write.table(Alex_W_C_lowcovDMRs, "Methylkit_AlexWet-Dry_25plowcovDMRs.txt", sep = "\t", quote = FALSE, row.names = FALSE)


##################################
# Compare High and Low Arctic


salloc -c1 --time 2:00:00 --mem 120000m --account def-cronk

R 

library(methylKit)

# Set the working directory to the folder containing the text files
setwd("/lustre04/scratch/celphin/Dryas_large_folders/CpG/Wild_stranded_CpG_report")

# List all text files in the folder with a specific pattern 
files <- list.files(pattern = "*CpG_report.txt$")

# Define sample ids (replace with your actual sample names)
sample_ids <- sub("_R1_val_1_bismark_bt2_pe\\.deduplicated\\.CpG_report\\.txt$", "", files)

covariate <- sub("^[^_]*_([A-Za-z]{3}).*", "\\1", sample_ids)

# Define treatment conditions (0 = control, 1 = treated)
treatment <- c(rep(0, length(files)/2), rep(1, length(files)/2))

# Replace Dry sites with 0 and Wet sites with 1
treatment <- gsub("CAS", "0", covariate)
treatment <- gsub("DRY", "0", treatment)
treatment <- gsub("WIL", "0", treatment)
treatment <- gsub("FER", "0", treatment)
treatment <- gsub("MEA", "0", treatment)
treatment <- gsub("SVA", "0", treatment)

treatment <- gsub("ALA", "1", treatment)
treatment <- gsub("LAT", "1", treatment)

# Coverage header: <chromosome> <start position> <end position> <methylation percentage> <count methylated> <count unmethylated>
# https://rdrr.io/bioc/methylKit/man/methRead-methods.html

d.f <- as.data.frame(cbind(sample_ids, treatment, covariate))

HL_d.f <- d.f[c(which(d.f$covariate=="CAS"), 
which(d.f$covariate=="WIL"), 
which(d.f$covariate=="DRY"), 
which(d.f$covariate=="MEA"), 
which(d.f$covariate=="FER"),
which(d.f$covariate=="SVA"),
which(d.f$covariate=="ALA"),
which(d.f$covariate=="LAT")),]

load("MethylKit_lowcov_data.RData")

HL_myobj <- reorganize(
  myobj_lowCov,
  as.list(HL_d.f$sample_ids),
  as.numeric(HL_d.f$treatment),
  chunk.size = 1e+06,
  save.db = TRUE,
  suffix = "HL",
  dbtype = "tabix"
)

# Print the object to check the result
#print(HL_myobj)
save(HL_myobj, file = "HL_MethylKit_lowcov_data.RData")
load("HL_MethylKit_lowcov_data.RData")

#------------------------------
# make regions
lowtiles = tileMethylCounts(HL_myobj,win.size=1000,step.size=1000,cov.bases = 10)
head(lowtiles[[1]],3)
save(lowtiles, file = "HL_MethylKit_tiles_lowcov.RData")
load("HL_MethylKit_tiles_lowcov.RData")

#------------------------
methtiles=unite(lowtiles, destrand=FALSE)
save(methtiles, file = "HL_united_lowtiles.RData")
load("HL_united_lowtiles.RData")

#---------------------------
# calculate DMRs

myDifftiles=calculateDiffMeth(methtiles, covariates=HL_d.f$covariates, overdispersion="MN",test="Chisq")
save(myDifftiles, file = "HL_MethDiffRegions.RData")
load("HL_MethDiffRegions.RData")

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

# [1] 12822
# [1] 6943
# [1] 5879


HL_lowcovDMRs <- getData(myDifftiles10p)
write.table(HL_lowcovDMRs, "Methylkit_HL_10plowcovDMRs.txt", sep = "\t", quote = FALSE, row.names = FALSE)

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

# [1] 1464
# [1] 845
# [1] 619


HL_lowcovDMRs <- getData(myDifftiles25p)
write.table(HL_lowcovDMRs, "Methylkit_HL_25plowcovDMRs.txt", sep = "\t", quote = FALSE, row.names = FALSE)



######################################
# Cassiope/Dryas/Willow/Fert/Meadow sites

cd /lustre04/scratch/celphin/Dryas_large_folders/CpG/Wild_stranded_CpG_report

tmux new-session -s methwarm2
tmux attach-session -t methwarm2

salloc -c1 --time 2:00:00 --mem 120000m --account def-cronk

module load StdEnv/2023
module load r/4.4.0

R

library(methylKit)

# Set the working directory to the folder containing the text files
setwd("/lustre04/scratch/celphin/Dryas_large_folders/CpG/Wild_stranded_CpG_report")

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

CASS_d.f <- d.f[which(d.f$covariate=="CAS"),]
WILL_d.f <- d.f[which(d.f$covariate=="WIL"),]
DRY_d.f <- d.f[which(d.f$covariate=="DRY"),]
MEAD_d.f <- d.f[which(d.f$covariate=="MEA"),]
FERT_d.f <- d.f[which(d.f$covariate=="FER"),]

load("MethylKit_lowcov_data.RData")

####################################
# Cassiope

CASS_myobj <- reorganize(
  myobj_lowCov,
  as.list(CASS_d.f$sample_ids),
  as.numeric(CASS_d.f$treatment),
  chunk.size = 1e+06,
  save.db = TRUE,
  suffix = "Alex",
  dbtype = "tabix"
)

# Print the object to check the result
print(CASS_myobj)
save(CASS_myobj, file = "CASS_MethylKit_lowcov_data.RData")
load("CASS_MethylKit_lowcov_data.RData")

meth=unite(CASS_myobj, destrand=FALSE)
save(meth, file = "CASS_united.RData")

#-------------------------------

load("CASS_united.RData")
pdf("CASS_PCA_PC1_PC2.pdf", width = 10, height = 8)
PCASamples(meth)
dev.off()
pdf("CASS_PCA_PC3_PC4.pdf", width = 10, height = 8)
PCASamples(meth, comp=c(3,4))
dev.off()

pdf("CASS_Dendrogram_plot.pdf", width = 10, height = 8)
clusterSamples(meth, dist="correlation", method="ward.D2", plot=TRUE)
dev.off()

#------------------------------
# make nicer PCA

allDataPCA <- PCASamples(meth, obj.return = TRUE) #Run a PCA analysis on percent methylation for all samples. methylKit uses prcomp to create the PCA matrix
summary(allDataPCA)

library(vegan)
pdf("CASS_PCA_PC1_PC2_ordiplot.pdf", width = 10, height = 8)
fig.allDataPCA <- ordiplot(allDataPCA, choices = c(1, 2), type = "none", display = "sites", cex = 0.5, xlim = c(-400, 200), xlab = "", ylab = "", xaxt = "n", yaxt = "n") #Use ordiplot to create base biplot. Do not add any points
points(fig.allDataPCA, "sites", col = c(rep("blue", times = 25), rep("red", times = 25)), pch = c(rep(16, times = 25), rep(17, times = 25)), cex = 3) #Add each sample. Darker samples are ambient, lighter samples are elevated pCO2

ordiellipse(allDataPCA, CASS_d.f$treatment, show.groups = "1", col = "red") #Add confidence ellipse around the samples in elevated pCO2
ordiellipse(allDataPCA, CASS_d.f$treatment, show.groups = "0", col = "blue") #Add confidence ellipse around the samples in ambient pCO2

axis(side =  1, at = seq(-400, 200, 200), col = "grey80", cex.axis = 1.7) #Add x-axis
mtext(side = 1, text = "PC 1", line = 3, cex = 1.5) #Add x-axis label

axis(side =  2, labels = TRUE, col = "grey80", cex.axis = 1.7) #Add y-axis
mtext(side = 2, text = "PC 2", line = 3, cex = 1.5) #Add y-axis label

mtext(side = 3, line = -5, adj = c(-100,0), text = "All CpG Loci", cex = 1.5)

dev.off()


#------------------------------
# make regions
lowtiles = tileMethylCounts(CASS_myobj,win.size=1000,step.size=1000,cov.bases = 10)
head(lowtiles[[1]],3)
save(lowtiles, file = "CASS_MethylKit_tiles_lowcov.RData")
load("CASS_MethylKit_tiles_lowcov.RData")

#------------------------
methtiles=unite(lowtiles, destrand=FALSE)
save(methtiles, file = "CASS_united_lowtiles.RData")
load("CASS_united_lowtiles.RData")

#---------------------------
# calculate DMRs

myDifftiles=calculateDiffMeth(methtiles, overdispersion="MN", test="Chisq")
save(myDifftiles, file = "CASS_MethDiffRegions.RData")
load("CASS_MethDiffRegions.RData")

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

# [1] 4962
# [1] 2739
# [1] 2223

# correct for overdispersion with F-test get 0
# correct for overdispersion with Chisq get:
# [1] 93
# [1] 58
# [1] 35


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

# [1] 187
# [1] 101
# [1] 86

# correct for overdispersion with F-test get 0
# correct for overdispersion with Chisq get:
# [1] 11
# [1] 7
# [1] 4

CASS_W_C_lowcovDMRs <- getData(myDifftiles10p)
write.table(CASS_W_C_lowcovDMRs, "Methylkit_CASS_W_C_10plowcovDMRs.txt", sep = "\t", quote = FALSE, row.names = FALSE)

CASS_W_C_lowcovDMRs <- getData(myDifftiles25p)
write.table(CASS_W_C_lowcovDMRs, "Methylkit_CASS_W_C_25plowcovDMRs.txt", sep = "\t", quote = FALSE, row.names = FALSE)

#####################################
# DRY

files0 <- list.files(pattern = "*CpG_report.txt$")
files <- files0[grep("_DRY", files0)]

# Read all the files into methylRawList using methRead
Dry_myobj0 <- methRead(
  location = as.list(files), 
  sample.id = as.list(DRY_d.f$sample_ids), 
  assembly = "Dryas",
  treatment = as.numeric(DRY_d.f$treatment),
  dbdir = getwd(),
  context = "CpG", 
  dbtype = "tabix",
  resolution = "base", 
  mincov = 3,
  pipeline = "bismarkCytosineReport" 
)

#------------------
# Print the object to check the result
print(Dry_myobj0)
save(Dry_myobj0, file = "Dryas_site_MethylKit_data.RData")

load("Dryas_site_MethylKit_data.RData")

DRY_myobj <- reorganize(
  Dry_myobj0,
  as.list(DRY_d.f$sample_ids),
  as.numeric(DRY_d.f$treatment),
  chunk.size = 1e+06,
  save.db = TRUE,
  suffix = "Alex",
  dbtype = "tabix"
)

# Print the object to check the result
print(DRY_myobj)
save(DRY_myobj, file = "DRY_MethylKit_lowcov_data.RData")
load("DRY_MethylKit_lowcov_data.RData")

meth=unite(DRY_myobj, destrand=FALSE)
save(meth, file = "DRY_united.RData")

#-------------------------------

load("DRY_united.RData")
pdf("DRY_PCA_PC1_PC2.pdf", width = 10, height = 8)
PCASamples(meth)
dev.off()
pdf("DRY_PCA_PC3_PC4.pdf", width = 10, height = 8)
PCASamples(meth, comp=c(3,4))
dev.off()

pdf("DRY_Dendrogram_plot.pdf", width = 10, height = 8)
clusterSamples(meth, dist="correlation", method="ward.D2", plot=TRUE)
dev.off()

#------------------------------
# make nicer PCA

allDataPCA <- PCASamples(meth, obj.return = TRUE) #Run a PCA analysis on percent methylation for all samples. methylKit uses prcomp to create the PCA matrix
summary(allDataPCA)

library(vegan)
pdf("DRY_PCA_PC1_PC2_ordiplot.pdf", width = 10, height = 8)
fig.allDataPCA <- ordiplot(allDataPCA, choices = c(1, 2), type = "none", display = "sites", cex = 0.5, xlim = c(-400, 200), xlab = "", ylab = "", xaxt = "n", yaxt = "n") #Use ordiplot to create base biplot. Do not add any points
points(fig.allDataPCA, "sites", col = c(rep("blue", times = 25), rep("red", times = 25)), pch = c(rep(16, times = 25), rep(17, times = 25)), cex = 3) #Add each sample. Darker samples are ambient, lighter samples are elevated pCO2

ordiellipse(allDataPCA, DRY_d.f$treatment, show.groups = "1", col = "red") #Add confidence ellipse around the samples in elevated pCO2
ordiellipse(allDataPCA, DRY_d.f$treatment, show.groups = "0", col = "blue") #Add confidence ellipse around the samples in ambient pCO2

axis(side =  1, at = seq(-400, 200, 200), col = "grey80", cex.axis = 1.7) #Add x-axis
mtext(side = 1, text = "PC 1", line = 3, cex = 1.5) #Add x-axis label

axis(side =  2, labels = TRUE, col = "grey80", cex.axis = 1.7) #Add y-axis
mtext(side = 2, text = "PC 2", line = 3, cex = 1.5) #Add y-axis label

mtext(side = 3, line = -5, adj = c(-100,0), text = "    a. All CpG Loci", cex = 1.5)

dev.off()


#------------------------------
# make regions
lowtiles = tileMethylCounts(DRY_myobj,win.size=1000,step.size=1000,cov.bases = 10)
head(lowtiles[[1]],3)
save(lowtiles, file = "DRY_MethylKit_tiles_lowcov.RData")
load("DRY_MethylKit_tiles_lowcov.RData")

#------------------------
methtiles=unite(lowtiles, destrand=FALSE)
save(methtiles, file = "DRY_united_lowtiles.RData")
load("DRY_united_lowtiles.RData")

#---------------------------
# calculate DMRs

myDifftiles=calculateDiffMeth(methtiles, overdispersion="MN", test="Chisq")
save(myDifftiles, file = "DRY_MethDiffRegions.RData")
load("DRY_MethDiffRegions.RData")

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

# [1] 5936
# [1] 3049
# [1] 2887

# correct for overdispersion with Chisq get:
# [1] 18
# [1] 9
# [1] 9


DRY_W_C_lowcovDMRs <- getData(myDifftiles10p)
write.table(DRY_W_C_lowcovDMRs, "Methylkit_DRY_W_C_10plowcovDMRs.txt", sep = "\t", quote = FALSE, row.names = FALSE)

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

# [1] 231
# [1] 111
# [1] 120
# correct for overdispersion with Chisq get:


DRY_W_C_lowcovDMRs <- getData(myDifftiles25p)
write.table(DRY_W_C_lowcovDMRs, "Methylkit_DRY_W_C_25plowcovDMRs.txt", sep = "\t", quote = FALSE, row.names = FALSE)

#####################################
# Meadow

MEAD_myobj <- reorganize(
  myobj_lowCov,
  as.list(MEAD_d.f$sample_ids),
  as.numeric(MEAD_d.f$treatment),
  chunk.size = 1e+06,
  save.db = TRUE,
  suffix = "Alex",
  dbtype = "tabix"
)

# Print the object to check the result
#print(MEAD_myobj)
save(MEAD_myobj, file = "MEAD_MethylKit_lowcov_data.RData")
load("MEAD_MethylKit_lowcov_data.RData")

meth=unite(MEAD_myobj, destrand=FALSE)
save(meth, file = "MEAD_united.RData")

#-------------------------------

load("MEAD_united.RData")
pdf("MEAD_PCA_PC1_PC2.pdf", width = 10, height = 8)
PCASamples(meth)
dev.off()
pdf("MEAD_PCA_PC3_PC4.pdf", width = 10, height = 8)
PCASamples(meth, comp=c(3,4))
dev.off()

pdf("MEAD_Dendrogram_plot.pdf", width = 10, height = 8)
clusterSamples(meth, dist="correlation", method="ward.D2", plot=TRUE)
dev.off()

#------------------------------
# make nicer PCA

allDataPCA <- PCASamples(meth, obj.return = TRUE) #Run a PCA analysis on percent methylation for all samples. methylKit uses prcomp to create the PCA matrix
summary(allDataPCA)

library(vegan)
pdf("MEAD_PCA_PC1_PC2_ordiplot.pdf", width = 10, height = 8)
fig.allDataPCA <- ordiplot(allDataPCA, choices = c(1, 2), type = "none", display = "sites", cex = 0.5, xlim = c(-400, 200), xlab = "", ylab = "", xaxt = "n", yaxt = "n") #Use ordiplot to create base biplot. Do not add any points
points(fig.allDataPCA, "sites", col = c(rep("blue", times = 25), rep("red", times = 25)), pch = c(rep(16, times = 25), rep(17, times = 25)), cex = 3) #Add each sample. Darker samples are ambient, lighter samples are elevated pCO2

ordiellipse(allDataPCA, MEAD_d.f$treatment, show.groups = "1", col = "red") #Add confidence ellipse around the samples in elevated pCO2
ordiellipse(allDataPCA, MEAD_d.f$treatment, show.groups = "0", col = "blue") #Add confidence ellipse around the samples in ambient pCO2

axis(side =  1, at = seq(-400, 200, 200), col = "grey80", cex.axis = 1.7) #Add x-axis
mtext(side = 1, text = "PC 1", line = 3, cex = 1.5) #Add x-axis label

axis(side =  2, labels = TRUE, col = "grey80", cex.axis = 1.7) #Add y-axis
mtext(side = 2, text = "PC 2", line = 3, cex = 1.5) #Add y-axis label

mtext(side = 3, line = -5, adj = c(-100,0), text = "    a. All CpG Loci", cex = 1.5)

dev.off()


#------------------------------
# make regions
lowtiles = tileMethylCounts(MEAD_myobj,win.size=1000,step.size=1000,cov.bases = 10)
head(lowtiles[[1]],3)
save(lowtiles, file = "MEAD_MethylKit_tiles_lowcov.RData")
load("MEAD_MethylKit_tiles_lowcov.RData")

#------------------------
methtiles=unite(lowtiles, destrand=FALSE)
save(methtiles, file = "MEAD_united_lowtiles.RData")
load("MEAD_united_lowtiles.RData")

#---------------------------
# calculate DMRs

myDifftiles=calculateDiffMeth(methtiles, overdispersion="MN", test="Chisq")
save(myDifftiles, file = "MEAD_MethDiffRegions.RData")
load("MEAD_MethDiffRegions.RData")

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

# [1] 8526
# [1] 3413
# [1] 5113

# correct for overdispersion with Chisq get:



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

# [1] 477
# [1] 197
# [1] 280

# correct for overdispersion with Chisq get:





MEAD_W_C_lowcovDMRs <- getData(myDifftiles10p)
write.table(MEAD_W_C_lowcovDMRs, "Methylkit_MEAD_W_C_10plowcovDMRs.txt", sep = "\t", quote = FALSE, row.names = FALSE)

MEAD_W_C_lowcovDMRs <- getData(myDifftiles25p)
write.table(MEAD_W_C_lowcovDMRs, "Methylkit_MEAD_W_C_25plowcovDMRs.txt", sep = "\t", quote = FALSE, row.names = FALSE)

#####################################
# Willow

WILL_myobj <- reorganize(
  myobj_lowCov,
  as.list(WILL_d.f$sample_ids),
  as.numeric(WILL_d.f$treatment),
  chunk.size = 1e+06,
  save.db = TRUE,
  suffix = "Alex",
  dbtype = "tabix"
)

# Print the object to check the result
print(WILL_myobj)
save(WILL_myobj, file = "WILL_MethylKit_lowcov_data.RData")
load("WILL_MethylKit_lowcov_data.RData")

meth=unite(WILL_myobj, destrand=FALSE)
save(meth, file = "WILL_united.RData")

#-------------------------------

load("WILL_united.RData")
pdf("WILL_PCA_PC1_PC2.pdf", width = 10, height = 8)
PCASamples(meth)
dev.off()
pdf("WILL_PCA_PC3_PC4.pdf", width = 10, height = 8)
PCASamples(meth, comp=c(3,4))
dev.off()

pdf("WILL_Dendrogram_plot.pdf", width = 10, height = 8)
clusterSamples(meth, dist="correlation", method="ward.D2", plot=TRUE)
dev.off()

#------------------------------
# make nicer PCA

allDataPCA <- PCASamples(meth, obj.return = TRUE) #Run a PCA analysis on percent methylation for all samples. methylKit uses prcomp to create the PCA matrix
summary(allDataPCA)

library(vegan)
pdf("WILL_PCA_PC1_PC2_ordiplot.pdf", width = 10, height = 8)
fig.allDataPCA <- ordiplot(allDataPCA, choices = c(1, 2), type = "none", display = "sites", cex = 0.5, xlim = c(-400, 200), xlab = "", ylab = "", xaxt = "n", yaxt = "n") #Use ordiplot to create base biplot. Do not add any points
points(fig.allDataPCA, "sites", col = c(rep("blue", times = 25), rep("red", times = 25)), pch = c(rep(16, times = 25), rep(17, times = 25)), cex = 3) #Add each sample. Darker samples are ambient, lighter samples are elevated pCO2

ordiellipse(allDataPCA, WILL_d.f$treatment, show.groups = "1", col = "red") #Add confidence ellipse around the samples in elevated pCO2
ordiellipse(allDataPCA, WILL_d.f$treatment, show.groups = "0", col = "blue") #Add confidence ellipse around the samples in ambient pCO2

axis(side =  1, at = seq(-400, 200, 200), col = "grey80", cex.axis = 1.7) #Add x-axis
mtext(side = 1, text = "PC 1", line = 3, cex = 1.5) #Add x-axis label

axis(side =  2, labels = TRUE, col = "grey80", cex.axis = 1.7) #Add y-axis
mtext(side = 2, text = "PC 2", line = 3, cex = 1.5) #Add y-axis label

mtext(side = 3, line = -5, adj = c(-100,0), text = "    a. All CpG Loci", cex = 1.5)

dev.off()


#------------------------------
# make regions
lowtiles = tileMethylCounts(WILL_myobj,win.size=1000,step.size=1000,cov.bases = 10)
head(lowtiles[[1]],3)
save(lowtiles, file = "WILL_MethylKit_tiles_lowcov.RData")
load("WILL_MethylKit_tiles_lowcov.RData")

#------------------------
methtiles=unite(lowtiles, destrand=FALSE)
save(methtiles, file = "WILL_united_lowtiles.RData")
load("WILL_united_lowtiles.RData")

#---------------------------
# calculate DMRs

myDifftiles=calculateDiffMeth(methtiles, overdispersion="MN", test="Chisq")
save(myDifftiles, file = "WILL_MethDiffRegions.RData")
load("WILL_MethDiffRegions.RData")

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

# [1] 7852
# [1] 3714
# [1] 4138

# correct for overdispersion with Chisq get:




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

# [1] 394
# [1] 184
# [1] 210

# correct for overdispersion with Chisq get:




WILL_W_C_lowcovDMRs <- getData(myDifftiles10p)
write.table(WILL_W_C_lowcovDMRs, "Methylkit_WILL_W_C_10plowcovDMRs.txt", sep = "\t", quote = FALSE, row.names = FALSE)

WILL_W_C_lowcovDMRs <- getData(myDifftiles25p)
write.table(WILL_W_C_lowcovDMRs, "Methylkit_WILL_W_C_25plowcovDMRs.txt", sep = "\t", quote = FALSE, row.names = FALSE)

#####################################
# Fertilizer

FERT_myobj <- reorganize(
  myobj_lowCov,
  as.list(FERT_d.f$sample_ids),
  as.numeric(FERT_d.f$treatment),
  chunk.size = 1e+06,
  save.db = TRUE,
  suffix = "Alex",
  dbtype = "tabix"
)

# Print the object to check the result
print(FERT_myobj)
save(FERT_myobj, file = "FERT_MethylKit_lowcov_data.RData")
load("FERT_MethylKit_lowcov_data.RData")

meth=unite(FERT_myobj, destrand=FALSE)
save(meth, file = "FERT_united.RData")

#-------------------------------

load("FERT_united.RData")
pdf("FERT_PCA_PC1_PC2.pdf", width = 10, height = 8)
PCASamples(meth)
dev.off()
pdf("FERT_PCA_PC3_PC4.pdf", width = 10, height = 8)
PCASamples(meth, comp=c(3,4))
dev.off()

pdf("FERT_Dendrogram_plot.pdf", width = 10, height = 8)
clusterSamples(meth, dist="correlation", method="ward.D2", plot=TRUE)
dev.off()

#------------------------------
# make nicer PCA

allDataPCA <- PCASamples(meth, obj.return = TRUE) #Run a PCA analysis on percent methylation for all samples. methylKit uses prcomp to create the PCA matrix
summary(allDataPCA)

library(vegan)
pdf("FERT_PCA_PC1_PC2_ordiplot.pdf", width = 10, height = 8)
fig.allDataPCA <- ordiplot(allDataPCA, choices = c(1, 2), type = "none", display = "sites", cex = 0.5, xlim = c(-400, 200), xlab = "", ylab = "", xaxt = "n", yaxt = "n") #Use ordiplot to create base biplot. Do not add any points
points(fig.allDataPCA, "sites", col = c(rep("blue", times = 25), rep("red", times = 25)), pch = c(rep(16, times = 25), rep(17, times = 25)), cex = 3) #Add each sample. Darker samples are ambient, lighter samples are elevated pCO2

ordiellipse(allDataPCA, FERT_d.f$treatment, show.groups = "1", col = "red") #Add confidence ellipse around the samples in elevated pCO2
ordiellipse(allDataPCA, FERT_d.f$treatment, show.groups = "0", col = "blue") #Add confidence ellipse around the samples in ambient pCO2

axis(side =  1, at = seq(-400, 200, 200), col = "grey80", cex.axis = 1.7) #Add x-axis
mtext(side = 1, text = "PC 1", line = 3, cex = 1.5) #Add x-axis label

axis(side =  2, labels = TRUE, col = "grey80", cex.axis = 1.7) #Add y-axis
mtext(side = 2, text = "PC 2", line = 3, cex = 1.5) #Add y-axis label

mtext(side = 3, line = -5, adj = c(-100,0), text = "    a. All CpG Loci", cex = 1.5)

dev.off()


#------------------------------
# make regions
lowtiles = tileMethylCounts(FERT_myobj,win.size=1000,step.size=1000,cov.bases = 10)
head(lowtiles[[1]],3)
save(lowtiles, file = "FERT_MethylKit_tiles_lowcov.RData")
load("FERT_MethylKit_tiles_lowcov.RData")

#------------------------
methtiles=unite(lowtiles, destrand=FALSE)
save(methtiles, file = "FERT_united_lowtiles.RData")
load("FERT_united_lowtiles.RData")

#---------------------------
# calculate DMRs

myDifftiles=calculateDiffMeth(methtiles, overdispersion="MN", test="Chisq")
save(myDifftiles, file = "FERT_MethDiffRegions.RData")
load("FERT_MethDiffRegions.RData")

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

# [1] 9436
# [1] 4346
# [1] 5090

# correct for overdispersion with Chisq get:



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

# [1] 549
# [1] 256
# [1] 293

# correct for overdispersion with Chisq get:





FERT_W_C_lowcovDMRs <- getData(myDifftiles10p)
write.table(FERT_W_C_lowcovDMRs, "Methylkit_FERT_W_C_10plowcovDMRs.txt", sep = "\t", quote = FALSE, row.names = FALSE)

FERT_W_C_lowcovDMRs <- getData(myDifftiles25p)
write.table(FERT_W_C_lowcovDMRs, "Methylkit_FERT_W_C_25plowcovDMRs.txt", sep = "\t", quote = FALSE, row.names = FALSE)

#####################################
# back in bash compare with metilene and overlay with genome annotations
# see bedGraph intersection notes