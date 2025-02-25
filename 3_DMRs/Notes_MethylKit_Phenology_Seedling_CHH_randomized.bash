###########################
# methylKit on Phenology Seedling data
# https://github.com/al2na/methylKit
# Jan 2025
##############################

# data

tmux new-session -s methpheno
tmux attach-session -t methpheno

tmux new-session -s methSE
tmux attach-session -t methSE

#-----------------------------
cd /lustre04/scratch/celphin/Dryas_large_folders/CHG_CHH/seedling_data/coverage/seedling/CHH_cytosine_report

cd /lustre04/scratch/celphin/Dryas_large_folders/CHG_CHH/seedling_data/coverage/phenology/CHH_cytosine_report

# copy over the other cytosine reports for Sen
cd /lustre04/scratch/celphin/Dryas_large_folders/CHG_CHH/data/coverage/CHH_cytosine_report

cp W1.F06.W1e5.CASS10W.544.60_Non_CpG.report.txt.CX_report.txt /lustre04/scratch/celphin/Dryas_large_folders/CHG_CHH/seedling_data/coverage/phenology/CHH_cytosine_report/SenFl_Cass_10W_60_544_Non_CpG.report.txt.CX_report.txt
cp C1.B05.C1d1.CASS4C.524.4_Non_CpG.report.txt.CX_report.txt /lustre04/scratch/celphin/Dryas_large_folders/CHG_CHH/seedling_data/coverage/phenology/CHH_cytosine_report/SenFl_Cass_4C_4_524_Non_CpG.report.txt.CX_report.txt
cp W1.A09.W1f10.CASS5W.525.130_Non_CpG.report.txt.CX_report.txt /lustre04/scratch/celphin/Dryas_large_folders/CHG_CHH/seedling_data/coverage/phenology/CHH_cytosine_report/SenFl_Cass_5W_130_525_Non_CpG.report.txt.CX_report.txt
cp C1.G07.C1h7.FERT5C.1F.97_Non_CpG.report.txt.CX_report.txt /lustre04/scratch/celphin/Dryas_large_folders/CHG_CHH/seedling_data/coverage/phenology/CHH_cytosine_report/SenFl_Fert_5C_97_1F_Non_CpG.report.txt.CX_report.txt
cp W1.B08.W1f8.FERT6W.3F.110_Non_CpG.report.txt.CX_report.txt /lustre04/scratch/celphin/Dryas_large_folders/CHG_CHH/seedling_data/coverage/phenology/CHH_cytosine_report/SenFl_Fert_6W_110_3F_Non_CpG.report.txt.CX_report.txt
cp C1.H05.C1f3.MEAD1C.446.33_Non_CpG.report.txt.CX_report.txt /lustre04/scratch/celphin/Dryas_large_folders/CHG_CHH/seedling_data/coverage/phenology/CHH_cytosine_report/SenFl_Mead_1C_33_446_Non_CpG.report.txt.CX_report.txt
cp W1.E08.W1c9.MEAD1W.444.116_Non_CpG.report.txt.CX_report.txt /lustre04/scratch/celphin/Dryas_large_folders/CHG_CHH/seedling_data/coverage/phenology/CHH_cytosine_report/SenFl_Mead_1W_116_444_Non_CpG.report.txt.CX_report.txt
cp C1.H07.C1b8.WILL3C.414.100_Non_CpG.report.txt.CX_report.txt /lustre04/scratch/celphin/Dryas_large_folders/CHG_CHH/seedling_data/coverage/phenology/CHH_cytosine_report/SenFl_Will_3C_100_414_Non_CpG.report.txt.CX_report.txt
cp W1.C05.W1g1.WILL4W.417.13_Non_CpG.report.txt.CX_report.txt /lustre04/scratch/celphin/Dryas_large_folders/CHG_CHH/seedling_data/coverage/phenology/CHH_cytosine_report/SenFl_Will_4W_13_417_Non_CpG.report.txt.CX_report.txt

cd /lustre04/scratch/celphin/Dryas_large_folders/CHG_CHH/seedling_data/coverage/phenology/CHH_cytosine_report
rename _Non_CpG.report.txt.CX_report.txt .CX_report.txt *

cd /lustre04/scratch/celphin/Dryas_large_folders/CHG_CHH/seedling_data/coverage/seedling/CHH_cytosine_report
rename _Non_CpG.report.txt.CX_report.txt .CX_report.txt *


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

tmux new-session -s methSE
tmux attach-session -t methSE

salloc -c1 --time 3:00:00 --mem 120000m --account def-henryg

module load StdEnv/2023
module load r/4.4.0

R

library(methylKit)

# Set the working directory to the folder containing the text files
setwd("/lustre04/scratch/celphin/Dryas_large_folders/CHG_CHH/seedling_data/coverage/seedling/CHH_cytosine_report")

# List all text files in the folder with a specific pattern 
files <- list.files(pattern = "*CX_report.txt$")

# Define sample ids (replace with your actual sample names)
sample_ids <- sub("\\.CX_report\\.txt$", "", files)

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

# subset for just Sweden "L" samples
SE_rand_d.f <- d.f[which(d.f$covariate=="L"),]
SE_files <- files[grep("SE.L.",files)]

# randomly rearrange treatment
treatment0 <- sample(SE_rand_d.f$treatment)
SE_rand_d.f$treatment0 <- treatment0 


# Read all the files into methylRawList using methRead
myobj <- methRead(
  location = as.list(SE_files), 
  sample.id = as.list(SE_rand_d.f$sample_ids), 
  assembly = "Dryas",
  treatment = as.numeric(SE_rand_d.f$treatment0),
  dbdir = getwd(),
  context = "CHH", 
  dbtype = "tabix",
  resolution = "base", 
  mincov = 3,
  pipeline = "bismarkCytosineReport" 
)

#------------------
# Print the object to check the result
#print(myobj)
save(myobj, file = "SE_MethylKit_data_CHH_rand.RData")

#--------------------------------------------------
# DMRs for Seedlings

load("SE_MethylKit_data_CHH_rand.RData")

SE_myobj <-  myobj

meth=unite(SE_myobj, destrand=FALSE)
save(meth, file = "SE_united_CHH_rand.RData")

#---------------
# plotting
load("SE_united_CHH_rand.RData")
pdf("SE_PCA_PC1_PC2_CHH_rand.pdf", width = 10, height = 8)
PCASamples(meth)
dev.off()
pdf("SE_Dendrogram_plot_CHH_rand.pdf", width = 10, height = 8)
clusterSamples(meth, dist="correlation", method="ward.D2", plot=TRUE)
dev.off()

#------------------------------
# make nicer PCA

allDataPCA <- PCASamples(meth, obj.return = TRUE) #Run a PCA analysis on percent methylation for all samples. methylKit uses prcomp to create the PCA matrix
summary(allDataPCA)

library(vegan)
pdf("SE_PCA_PC1_PC2_ordiplot_CHH_rand.pdf", width = 10, height = 8)
fig.allDataPCA <- ordiplot(allDataPCA, choices = c(1, 2), type = "none", display = "sites", cex = 0.5, xlim = c(-400, 200), xlab = "", ylab = "", xaxt = "n", yaxt = "n") #Use ordiplot to create base biplot. Do not add any points
points(fig.allDataPCA, "sites", col = c(rep("blue", times = 10), rep("red", times = 10)), pch = c(rep(16, times = 10), rep(17, times = 10)), cex = 3) #Add each sample. Darker samples are ambient, lighter samples are elevated pCO2

ordiellipse(allDataPCA, SE_rand_d.f$treatment, show.groups = "1", col = "red") #Add confidence ellipse around the samples in elevated pCO2
ordiellipse(allDataPCA, SE_rand_d.f$treatment, show.groups = "0", col = "blue") #Add confidence ellipse around the samples in ambient pCO2

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
save(lowtiles, file = "SE_MethylKit_tiles_lowcov_CHH_rand.RData")
load("SE_MethylKit_tiles_lowcov_CHH_rand.RData")

#------------------------
methtiles=unite(lowtiles, destrand = FALSE, min.per.group = NULL, chunk.size = 3e+05)
save(methtiles, file = "SE_united_lowtiles_CHH_rand.RData")
load("SE_united_lowtiles_CHH_rand.RData")

#---------------------------
# calculate DMRs

myDifftiles=calculateDiffMeth(methtiles, overdispersion="MN", test="Chisq")
save(myDifftiles, file = "SE_MethDiffRegions_CHH_rand.RData")
load("SE_MethDiffRegions_CHH_rand.RData")

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


# overdisp and chisq Sweden only - CHH rand
# 9
# 0
# 9


SE_W_C_DMRs <- getData(myDifftiles10p)
write.table(SE_W_C_DMRs, "Methylkit_SE_W_C_10_DMRs_CHH_rand.txt", sep = "\t", quote = FALSE, row.names = FALSE)

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

# none overdisp and chisq Sweden only


SE_W_C_DMRs <- getData(myDifftiles25p)
write.table(SE_W_C_DMRs, "Methylkit_SE_W_C_25_DMRs_CHH_rand.txt", sep = "\t", quote = FALSE, row.names = FALSE)

########################################
# load Phenology data

tmux new-session -s methpheno
tmux attach-session -t methpheno

salloc -c1 --time 5:00:00 --mem 120000m --account def-henryg


module load StdEnv/2023
module load r/4.4.0

R

library(methylKit)

# Set the working directory to the folder containing the text files
setwd("/lustre04/scratch/celphin/Dryas_large_folders/CHG_CHH/seedling_data/coverage/phenology/CHH_cytosine_report")

# List all text files in the folder with a specific pattern 
files <- list.files(pattern = "*CX_report.txt$")

# Define sample ids (replace with your actual sample names)
sample_ids <- sub("\\.CX_report\\.txt$", "", files)

# extract the site
covariate  <- substr(sample_ids, 7, 9)

# Extract first letter (S or M)
first_letter <- substr(sample_ids, 1, 1)

# Replace C with 0 and W with 1
treatment <- gsub("S", "0", first_letter)
treatment <- gsub("M", "1", treatment)

# Coverage header: <chromosome> <start position> <end position> <methylation percentage> <count methylated> <count unmethylated>
# https://rdrr.io/bioc/methylKit/man/methRead-methods.html

rand_d.f <- as.data.frame(cbind(sample_ids, treatment, covariate))

# randomly rearrange treatment
treatment0 <- sample(rand_d.f$treatment)
rand_d.f$treatment0 <- treatment0 


# Read all the files into methylRawList using methRead
myobj <- methRead(
  location = as.list(files), 
  sample.id = as.list(rand_d.f$sample_ids), 
  assembly = "Dryas",
  treatment = as.numeric(rand_d.f$treatment0),
  dbdir = getwd(),
  context ="CHH", 
  dbtype = "tabix",
  resolution = "base", 
  mincov = 3,
  pipeline = "bismarkCytosineReport" 
)


#------------------
# Print the object to check the result
#print(myobj)
save(myobj, file = "Pheno_MethylKit_CHH_data_rand.RData")

load("Pheno_MethylKit_CHH_data_rand.RData")

# Pheno_myobj <- reorganize(
  # myobj,
  # as.list(d.f$sample_ids),
  # as.numeric(d.f$treatment),
  # chunk.size = 1e+06,
  # save.db = TRUE,
  # suffix = "Pheno",
  # dbtype = "tabix"
# )

# Print the object to check the result
#print(Pheno_myobj)
#save(Pheno_myobj, file = "Pheno_MethylKit_lowcov_data_rand.RData")
#load("Pheno_MethylKit_lowcov_data_rand.RData")

Pheno_myobj <-  myobj
meth=unite(Pheno_myobj, destrand=FALSE)
save(meth, file = "Pheno_united_CHH_rand.RData")

#----------
# plotting
load("Pheno_united_CHH_rand.RData")
pdf("Pheno_CHH_PCA_PC1_PC2_rand.pdf", width = 10, height = 8)
PCASamples(meth)
dev.off()
pdf("Pheno_CHH_Dendrogram_plot_rand.pdf", width = 10, height = 8)
clusterSamples(meth, dist="correlation", method="ward.D2", plot=TRUE)
dev.off()

#------------------------------
# make nicer PCA

allDataPCA <- PCASamples(meth, obj.return = TRUE) #Run a PCA analysis on percent methylation for all samples. methylKit uses prcomp to create the PCA matrix
summary(allDataPCA)

library(vegan)
pdf("Pheno_CHH_PCA_PC1_PC2_ordiplot_rand.pdf", width = 10, height = 8)
fig.allDataPCA <- ordiplot(allDataPCA, choices = c(1, 2), type = "none", display = "sites", cex = 0.5, xlim = c(-400, 200), xlab = "", ylab = "", xaxt = "n", yaxt = "n") #Use ordiplot to create base biplot. Do not add any points
points(fig.allDataPCA, "sites", col = c(rep("blue", times = 10), rep("red", times = 10)), pch = c(rep(16, times = 10), rep(17, times = 10)), cex = 3) #Add each sample. Darker samples are ambient, lighter samples are elevated pCO2

ordiellipse(allDataPCA, rand_d.f$treatment, show.groups = "1", col = "red") #Add confidence ellipse around the samples in elevated pCO2
ordiellipse(allDataPCA, rand_d.f$treatment, show.groups = "0", col = "blue") #Add confidence ellipse around the samples in ambient pCO2

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
save(lowtiles, file = "Pheno_MethylKit_tiles_lowcov_CHH_rand.RData")
load("Pheno_MethylKit_tiles_lowcov_CHH_rand.RData")

#------------------------
methtiles=unite(lowtiles, destrand=FALSE)
save(methtiles, file = "Pheno_united_lowtiles_CHH_rand.RData")
load("Pheno_united_lowtiles_CHH_rand.RData")

#---------------------------
# calculate DMRs

myDifftiles=calculateDiffMeth(methtiles, overdispersion="MN", covariates=rand_d.f$covariates, test="Chisq")
save(myDifftiles, file = "Pheno_MethDiffRegions_CHH_rand.RData")
load("Pheno_MethDiffRegions_CHH_rand.RData")

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

# with overdispersion and chisq  - CHH rand
# 4
# 4
# 0

Pheno_DMRs <- getData(myDifftiles10p)
write.table(Pheno_DMRs, "Methylkit_Pheno_10_DMRs_CHH_rand.txt", sep = "\t", quote = FALSE, row.names = FALSE)


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

# with overdispersion and chisq - CHH rand
# 2
# 2
# 0

Pheno_DMRs <- getData(myDifftiles25p)
write.table(Pheno_DMRs, "Methylkit_Pheno_25_DMRs_CHH_rand.txt", sep = "\t", quote = FALSE, row.names = FALSE)


#####################################
# back in bash compare with metilene and overlay with genome annotations
# see bedGraph intersection notes