###########################
# methylKit on Phenology Seedling data
# https://github.com/al2na/methylKit
# Jan 2025
##############################

# data

cd /lustre04/scratch/celphin/Dryas_large_folders/CpG

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

mkdir Seedling
mv SE* Seedling/

mkdir Phenology
mv MatFl* Phenology

# copy over the other cytosine reports for Sen
cd /lustre04/scratch/celphin/Dryas_large_folders/CpG/Wild_stranded_CpG_report

cp W1.F06.W1e5_CASS10W_544_60_R1_val_1_bismark_bt2_pe.deduplicated.CpG_report.txt ../stranded_CpG_report/Phenology/Sen_Cass_10W_60_544.CpG_report.txt
cp C1.B05.C1d1_CASS4C_524_4_R1_val_1_bismark_bt2_pe.deduplicated.CpG_report.txt ../stranded_CpG_report/Phenology/Sen_Cass_4C_4_524.CpG_report.txt
cp W1.A09.W1f10_CASS5W_525_130_R1_val_1_bismark_bt2_pe.deduplicated.CpG_report.txt ../stranded_CpG_report/Phenology/Sen_Cass_5W_130_525.CpG_report.txt
cp C1.G07.C1h7_FERT5C_1F_97_R1_val_1_bismark_bt2_pe.deduplicated.CpG_report.txt ../stranded_CpG_report/Phenology/Sen_Fert_5C_97_1F.CpG_report.txt
cp W1.B08.W1f8_FERT6W_3F_110_R1_val_1_bismark_bt2_pe.deduplicated.CpG_report.txt ../stranded_CpG_report/Phenology/Sen_Fert_6W_110_3F.CpG_report.txt
cp C1.H05.C1f3_MEAD1C_446_33_R1_val_1_bismark_bt2_pe.deduplicated.CpG_report.txt ../stranded_CpG_report/Phenology/Sen_Mead_1C_33_446.CpG_report.txt
cp W1.E08.W1c9_MEAD1W_444_116_R1_val_1_bismark_bt2_pe.deduplicated.CpG_report.txt ../stranded_CpG_report/Phenology/Sen_Mead_1W_116_444.CpG_report.txt
cp C1.H07.C1b8_WILL3C_414_100_R1_val_1_bismark_bt2_pe.deduplicated.CpG_report.txt ../stranded_CpG_report/Phenology/Sen_Will_3C_100_414.CpG_report.txt
cp W1.C05.W1g1_WILL4W_417_13_R1_val_1_bismark_bt2_pe.deduplicated.CpG_report.txt ../stranded_CpG_report/Phenology/Sen_Will_4W_13_417.CpG_report.txt

cd ../stranded_CpG_report/Phenology/
rename _R1_val_1_bismark_bt2_pe.deduplicated. . *

cd /lustre04/scratch/celphin/Dryas_large_folders/CpG/stranded_CpG_report/Seedling
rename _R1_val_1_bismark_bt2_pe.deduplicated. . *


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


###############################
# load Seedling data

library(methylKit)

# Set the working directory to the folder containing the text files
setwd("/lustre04/scratch/celphin/Dryas_large_folders/CpG/stranded_CpG_report/Seedling")

# List all text files in the folder with a specific pattern 
files <- list.files(pattern = "*CpG_report.txt$")

# Define sample ids (replace with your actual sample names)
sample_ids <- sub("\\.CpG_report\\.txt$", "", files)

covariate  <- substr(sample_ids, 4, 4)

# Define treatment conditions (0 = control, 1 = treated)
treatment <- c(rep(0, length(files)/2), rep(1, length(files)/2))

# Extract first letter (C or W)
first_letter <- substr(sample_ids, 6, 6)

# Replace C with 0 and W with 1
treatment <- gsub("C", "0", first_letter)
treatment <- gsub("W", "1", treatment)

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
  mincov = 3,
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
save(myobj, file = "SE_MethylKit_data.RData")

#--------------------------------------------------
# DMRs for Seedlings

library(methylKit)

load("SE_MethylKit_data.RData")

# subset for just Sweden "L" samples
SE_d.f <- d.f[which(d.f$covariate=="L"),]

SE_myobj <- reorganize(
  myobj,
  as.list(SE_d.f$sample_ids),
  as.numeric(SE_d.f$treatment),
  chunk.size = 1e+06,
  save.db = TRUE,
  suffix = "SE",
  dbtype = "tabix"
)

# Print the object to check the result
#print(SE_myobj)
save(SE_myobj, file = "SE_MethylKit_lowcov_data.RData")
load("SE_MethylKit_lowcov_data.RData")

meth=unite(SE_myobj, destrand=FALSE)
save(meth, file = "SE_united.RData")

#-----------------------------
# Plotting
load("SE_united.RData")
pdf("SE_PCA_PC1_PC2.pdf", width = 10, height = 8)
PCASamples(meth)
dev.off()
pdf("SE_Dendrogram_plot.pdf", width = 10, height = 8)
clusterSamples(meth, dist="correlation", method="ward.D2", plot=TRUE)
dev.off()

#------------------------------
# make nicer PCA

allDataPCA <- PCASamples(meth, obj.return = TRUE) #Run a PCA analysis on percent methylation for all samples. methylKit uses prcomp to create the PCA matrix
summary(allDataPCA)

library(vegan)
pdf("SE_PCA_PC1_PC2_ordiplot.pdf", width = 10, height = 8)
fig.allDataPCA <- ordiplot(allDataPCA, choices = c(1, 2), type = "none", display = "sites", cex = 0.5, xlim = c(-400, 200), xlab = "", ylab = "", xaxt = "n", yaxt = "n") #Use ordiplot to create base biplot. Do not add any points
points(fig.allDataPCA, "sites", col = c(rep("blue", times = 10), rep("red", times = 10)), pch = c(rep(16, times = 10), rep(17, times = 10)), cex = 3) #Add each sample. Darker samples are ambient, lighter samples are elevated pCO2

ordiellipse(allDataPCA, SE_d.f$treatment, show.groups = "1", col = "red") #Add confidence ellipse around the samples in elevated pCO2
ordiellipse(allDataPCA, SE_d.f$treatment, show.groups = "0", col = "blue") #Add confidence ellipse around the samples in ambient pCO2

axis(side =  1, at = seq(-400, 200, 200), col = "grey80", cex.axis = 1.7) #Add x-axis
mtext(side = 1, text = "PC 1", line = 3, cex = 1.5) #Add x-axis label

axis(side =  2, labels = TRUE, col = "grey80", cex.axis = 1.7) #Add y-axis
mtext(side = 2, text = "PC 2", line = 3, cex = 1.5) #Add y-axis label

mtext(side = 3, line = -5, adj = c(-100,0), text = "    a. All CpG Loci", cex = 1.5)

dev.off()

#-----------------------
# make regions
lowtiles = tileMethylCounts(SE_myobj,win.size=1000,step.size=1000,cov.bases = 10)
head(lowtiles[[1]],3)
save(lowtiles, file = "SE_MethylKit_tiles_lowcov.RData")
load("SE_MethylKit_tiles_lowcov.RData")

#------------------------
methtiles=unite(lowtiles, destrand=FALSE)
save(methtiles, file = "SE_united_lowtiles.RData")
load("SE_united_lowtiles.RData")

#---------------------------
# calculate DMRs

myDifftiles=calculateDiffMeth(methtiles, overdispersion="MN",  test="Chisq")
save(myDifftiles, file = "SE_MethDiffRegions.RData")
load("SE_MethDiffRegions.RData")

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

# [1] 6829
# [1] 4210
# [1] 2619

# overdisp, chisq
# [1] 421
# [1] 260
# [1] 161

# overdisp and Chisq - Sweden only
# [1] 1516
# [1] 953
# [1] 563

SE_W_C_DMRs <- getData(myDifftiles10p)
write.table(SE_W_C_DMRs, "Methylkit_SE_W_C_10_DMRs.txt", sep = "\t", quote = FALSE, row.names = FALSE)

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

# [1] 642
# [1] 459
# [1] 183

# overdisp, chisq
# [1] 110
# [1] 70
# [1] 40

# overdisp and Chisq - Sweden only
# [1] 369
# [1] 231
# [1] 138


SE_W_C_DMRs <- getData(myDifftiles25p)
write.table(SE_W_C_DMRs, "Methylkit_SE_W_C_25_DMRs.txt", sep = "\t", quote = FALSE, row.names = FALSE)


###############################
# load Seedling data again to look at High Low comparison

library(methylKit)

# Set the working directory to the folder containing the text files
setwd("/lustre04/scratch/celphin/Dryas_large_folders/CpG/stranded_CpG_report/Seedling")

# List all text files in the folder with a specific pattern 
files <- list.files(pattern = "*CpG_report.txt$")

# Define sample ids (replace with your actual sample names)
sample_ids <- sub("\\.CpG_report\\.txt$", "", files)

covariate  <- substr(sample_ids, 4, 4)

# Define treatment conditions (0 = control, 1 = treated)
treatment <- c(rep(0, length(files)/2), rep(1, length(files)/2))

# Extract first letter (C or W)
first_letter <- substr(sample_ids, 8, 8)

# Replace C with 0 and W with 1
treatment <- gsub("H", "0", first_letter)
treatment <- gsub("L", "1", treatment)

# Coverage header: <chromosome> <start position> <end position> <methylation percentage> <count methylated> <count unmethylated>
# https://rdrr.io/bioc/methylKit/man/methRead-methods.html

d.f <- as.data.frame(cbind(sample_ids, treatment, covariate))

load("SE_MethylKit_data.RData")

# subset for just Sweden "L" samples
SE_HL_d.f <- d.f[which(d.f$covariate=="L"),]

SE_HL_myobj <- reorganize(
  myobj,
  as.list(SE_HL_d.f$sample_ids),
  as.numeric(SE_HL_d.f$treatment),
  chunk.size = 1e+06,
  save.db = TRUE,
  suffix = "SE_HL",
  dbtype = "tabix"
)

# Print the object to check the result
#print(SE_HL_myobj)
save(SE_HL_myobj, file = "SE_HL_MethylKit_lowcov_data.RData")
load("SE_HL_MethylKit_lowcov_data.RData")

#-----------------------
# make regions
lowtiles = tileMethylCounts(SE_HL_myobj,win.size=1000,step.size=1000,cov.bases = 10)
head(lowtiles[[1]],3)
save(lowtiles, file = "SE_HL_MethylKit_tiles_lowcov.RData")
load("SE_HL_MethylKit_tiles_lowcov.RData")

#------------------------
methtiles=unite(lowtiles, destrand=FALSE)
save(methtiles, file = "SE_HL_united_lowtiles.RData")
load("SE_HL_united_lowtiles.RData")

#---------------------------
# calculate DMRs

myDifftiles=calculateDiffMeth(methtiles, overdispersion="MN",  test="Chisq")
save(myDifftiles, file = "SE_HL_MethDiffRegions.RData")
load("SE_HL_MethDiffRegions.RData")

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

# [1] 415
# [1] 380
# [1] 35

SE_HL_DMRs <- getData(myDifftiles10p)
write.table(SE_HL_DMRs, "Methylkit_SE_HL_10_DMRs.txt", sep = "\t", quote = FALSE, row.names = FALSE)

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

# [1] 52
# [1] 48
# [1] 4


SE_HL_DMRs <- getData(myDifftiles25p)
write.table(SE_HL_DMRs, "Methylkit_SE_HL_25_DMRs.txt", sep = "\t", quote = FALSE, row.names = FALSE)

########################################
# load Phenology data

tmux new-session -s methylkit1
tmux attach-session -t methylkit1

module load StdEnv/2023
module load r/4.4.0

R

library(methylKit)

# Set the working directory to the folder containing the text files
setwd("/lustre04/scratch/celphin/Dryas_large_folders/CpG/stranded_CpG_report/Phenology")

# List all text files in the folder with a specific pattern 
files <- list.files(pattern = "*CpG_report.txt$")

# Define sample ids (replace with your actual sample names)
sample_ids <- sub("\\.CpG_report\\.txt$", "", files)

# extract the site
covariate  <- substr(sample_ids, 7, 9)

# Extract first letter (S or M)
first_letter <- substr(sample_ids, 1, 1)

# Replace C with 0 and W with 1
treatment <- gsub("S", "0", first_letter)
treatment <- gsub("M", "1", treatment)

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
  mincov = 3,
  pipeline = "bismarkCytosineReport" 
)


#------------------
# Print the object to check the result
#print(myobj)
save(myobj, file = "Pheno_MethylKit_data.RData")

load("Pheno_MethylKit_data.RData")

#--------------------------------------------------
# DMRs for Phenology

library(methylKit)

# Set the working directory to the folder containing the text files
setwd("/lustre04/scratch/celphin/Dryas_large_folders/CpG/stranded_CpG_report/Phenology")

load("Pheno_MethylKit_data.RData")

Pheno_myobj <- reorganize(
  myobj,
  as.list(d.f$sample_ids),
  as.numeric(d.f$treatment),
  chunk.size = 1e+06,
  save.db = TRUE,
  suffix = "Pheno",
  dbtype = "tabix"
)

# Print the object to check the result
print(Pheno_myobj)
save(Pheno_myobj, file = "Pheno_MethylKit_lowcov_data.RData")
load("Pheno_MethylKit_lowcov_data.RData")

meth=unite(Pheno_myobj, destrand=FALSE)
save(meth, file = "Pheno_united.RData")
load("Pheno_united.RData")
pdf("Pheno_PCA_PC1_PC2.pdf", width = 10, height = 8)
PCASamples(meth)
dev.off()
pdf("Pheno_Dendrogram_plot.pdf", width = 10, height = 8)
clusterSamples(meth, dist="correlation", method="ward.D2", plot=TRUE)
dev.off()

#------------------------------
# make nicer PCA

allDataPCA <- PCASamples(meth, obj.return = TRUE) #Run a PCA analysis on percent methylation for all samples. methylKit uses prcomp to create the PCA matrix
summary(allDataPCA)

library(vegan)
pdf("Pheno_PCA_PC1_PC2_ordiplot.pdf", width = 10, height = 8)
fig.allDataPCA <- ordiplot(allDataPCA, choices = c(1, 2), type = "none", display = "sites", cex = 0.5, xlim = c(-400, 200), xlab = "", ylab = "", xaxt = "n", yaxt = "n") #Use ordiplot to create base biplot. Do not add any points
points(fig.allDataPCA, "sites", col = c(rep("blue", times = 10), rep("red", times = 10)), pch = c(rep(16, times = 10), rep(17, times = 10)), cex = 3) #Add each sample. Darker samples are ambient, lighter samples are elevated pCO2

ordiellipse(allDataPCA, d.f$treatment, show.groups = "1", col = "red") #Add confidence ellipse around the samples in elevated pCO2
ordiellipse(allDataPCA, d.f$treatment, show.groups = "0", col = "blue") #Add confidence ellipse around the samples in ambient pCO2

axis(side =  1, at = seq(-400, 200, 200), col = "grey80", cex.axis = 1.7) #Add x-axis
mtext(side = 1, text = "PC 1", line = 3, cex = 1.5) #Add x-axis label

axis(side =  2, labels = TRUE, col = "grey80", cex.axis = 1.7) #Add y-axis
mtext(side = 2, text = "PC 2", line = 3, cex = 1.5) #Add y-axis label

mtext(side = 3, line = -5, adj = c(-100,0), text = "All CpG Loci", cex = 1.5)

dev.off()

#-----------------------
# make regions
lowtiles = tileMethylCounts(Pheno_myobj,win.size=1000,step.size=1000,cov.bases = 10)
head(lowtiles[[1]],3)
save(lowtiles, file = "Pheno_MethylKit_tiles_lowcov.RData")
load("Pheno_MethylKit_tiles_lowcov.RData")

#------------------------
methtiles=unite(lowtiles, destrand=FALSE)
save(methtiles, file = "Pheno_united_lowtiles.RData")
load("Pheno_united_lowtiles.RData")

#---------------------------
# calculate DMRs

myDifftiles=calculateDiffMeth(methtiles, overdispersion="MN", covariates=d.f$covariates, test="Chisq")
save(myDifftiles, file = "Pheno_MethDiffRegions.RData")
load("Pheno_MethDiffRegions.RData")

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

# [1] 17536
# [1] 1231
# [1] 16305

# with overdispersion and covariates and chisq - check actually CpG
# [1] 14383
# [1] 939
# [1] 13444

Pheno_DMRs <- getData(myDifftiles10p)
write.table(Pheno_DMRs, "Methylkit_Pheno_10_DMRs.txt", sep = "\t", quote = FALSE, row.names = FALSE)

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

# [1] 3594
# [1] 66
# [1] 3528

# with overdispersion and covariates and chisq
# [1] 3530
# [1] 63
# [1] 3467


Pheno_DMRs <- getData(myDifftiles25p)
write.table(Pheno_DMRs, "Methylkit_Pheno_25_DMRs.txt", sep = "\t", quote = FALSE, row.names = FALSE)


#####################################
# back in bash compare with metilene and overlay with genome annotations
# see bedGraph intersection notes