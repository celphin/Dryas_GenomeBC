###########################
# methylKit
# https://github.com/al2na/methylKit
# Dec 2024
##############################

# data

tmux new-session -s methwarm1
tmux attach-session -t methwarm1

cd /lustre04/scratch/celphin/Dryas_large_folders/CHG_CHH/data/coverage/CHH_cytosine_report

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


#####################################
# try running again for each site independently - set meth diff to 10 or 25%
## ----------------------------------------------------------------
# DMRs for Alaska

cd /lustre04/scratch/celphin/Dryas_large_folders/CHG_CHH/data/coverage/CHH_cytosine_report

tmux new-session -s methwarm1
tmux attach-session -t methwarm1

salloc -c1 --time 3:00:00 --mem 120000m --account def-henryg

module load StdEnv/2023
module load r/4.4.0

R

library(methylKit)

# Set the working directory to the folder containing the text files
setwd("/lustre04/scratch/celphin/Dryas_large_folders/CHG_CHH/data/coverage/CHH_cytosine_report")

# List all text files in the folder with a specific pattern 
files <- list.files(pattern = "*CX_report.txt$")

# Define sample ids (replace with your actual sample names)
sample_ids <- sub("_Non_CpG\\.report\\.txt\\.CX_report\\.txt$", "", files)

covariate <- sub("^([^\\.]*\\.[^\\.]*\\.[^\\.]*\\.([A-Za-z]{3}).*)", "\\2", sample_ids)

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

rand_d.f <- as.data.frame(cbind(sample_ids, treatment, covariate))

# subset for Alaska
ALAS_rand_d.f <- rand_d.f[which(rand_d.f$covariate=="ALA"),]
ALAS_files <- files[grep("ALAS",files)]

# randomly rearrange treatment
treatment0 <- sample(ALAS_rand_d.f$treatment)
ALAS_rand_d.f$treatment0 <- treatment0 

# Read all the files into methylRawList using methRead
ALAS_myobj <- methRead(
  location = as.list(ALAS_files), 
  sample.id = as.list(ALAS_rand_d.f$sample_ids), 
  assembly = "Dryas",
  treatment = as.numeric(ALAS_rand_d.f$treatment0),
  dbdir = getwd(),
  context ="CHH",
  dbtype = "tabix",
  resolution = "base", 
  mincov = 3,
  pipeline = "bismarkCytosineReport" 
)

# Print the object to check the result
#print(ALAS_myobj)
save(ALAS_myobj, file = "ALAS_MethylKit_lowcov_data_CHH_rand.RData")

#-----------------------
load("ALAS_MethylKit_lowcov_data_CHH_rand.RData")

meth=unite(ALAS_myobj, destrand=FALSE)
save(meth, file = "ALAS_united_CHH_rand.RData")

#-----------------------
# make regions
lowtiles = tileMethylCounts(ALAS_myobj,win.size=1000,step.size=1000,cov.bases = 10)
head(lowtiles[[1]],3)
save(lowtiles, file = "ALAS_MethylKit_tiles_lowcov_CHH_rand.RData")
load("ALAS_MethylKit_tiles_lowcov_CHH_rand.RData")

#------------------------
methtiles=unite(lowtiles, destrand=FALSE)
save(methtiles, file = "ALAS_united_lowtiles_CHH_rand.RData")
load("ALAS_united_lowtiles_CHH_rand.RData")

#---------------------------
# calculate DMRs

myDifftiles=calculateDiffMeth(methtiles, overdispersion="MN", test="Chisq")
save(myDifftiles, file = "ALAS_MethDiffRegions_CHH_rand.RData")
load("ALAS_MethDiffRegions_CHH_rand.RData")

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

# overdisp and chisq
2
2
0


ALAS_W_C_lowcovDMRs <- getData(myDifftiles10p)
write.table(ALAS_W_C_lowcovDMRs, "Methylkit_ALAS_W_C_10_CHH_rand_plowcovDMRs.txt", sep = "\t", quote = FALSE, row.names = FALSE)

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

# none overdisp and chisq

ALAS_W_C_lowcovDMRs <- getData(myDifftiles25p)
write.table(ALAS_W_C_lowcovDMRs, "Methylkit_ALAS_W_C_25_CHH_rand_plowcovDMRs.txt", sep = "\t", quote = FALSE, row.names = FALSE)

########################################
# DMRs for Svalbard

cd /lustre04/scratch/celphin/Dryas_large_folders/CHG_CHH/data/coverage/CHH_cytosine_report

tmux new-session -s methwarm1
tmux attach-session -t methwarm1

salloc -c1 --time 3:00:00 --mem 120000m --account def-henryg

module load StdEnv/2023
module load r/4.4.0

R

library(methylKit)

# Set the working directory to the folder containing the text files
setwd("/lustre04/scratch/celphin/Dryas_large_folders/CHG_CHH/data/coverage/CHH_cytosine_report")

# List all text files in the folder with a specific pattern 
files <- list.files(pattern = "*CX_report.txt$")

# Define sample ids (replace with your actual sample names)
sample_ids <- sub("_Non_CpG\\.report\\.txt\\.CX_report\\.txt$", "", files)

covariate <- sub("^([^\\.]*\\.[^\\.]*\\.[^\\.]*\\.([A-Za-z]{3}).*)", "\\2", sample_ids)


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

rand_d.f <- as.data.frame(cbind(sample_ids, treatment, covariate))

# subset for Svalbard
SVAL_rand_d.f <- rand_d.f[which(rand_d.f$covariate=="SVA"),]
SVAL_files <- files[grep("SVA",files)]

# randomly rearrange treatment
treatment0 <- sample(SVAL_rand_d.f$treatment)
SVAL_rand_d.f$treatment0 <- treatment0 

# Read all the files into methylRawList using methRead
SVAL_myobj <- methRead(
  location = as.list(SVAL_files), 
  sample.id = as.list(SVAL_rand_d.f$sample_ids), 
  assembly = "Dryas",
  treatment = as.numeric(SVAL_rand_d.f$treatment0),
  dbdir = getwd(),
  context ="CHH",
  dbtype = "tabix",
  resolution = "base", 
  mincov = 3,
  pipeline = "bismarkCytosineReport" 
)


# Print the object to check the result
#print(SVAL_myobj)
save(SVAL_myobj, file = "SVAL_MethylKit_lowcov_data_CHH_rand.RData")
load("SVAL_MethylKit_lowcov_data_CHH_rand.RData")

meth=unite(SVAL_myobj, destrand=FALSE)
save(meth, file = "SVAL_united_CHH_rand.RData")

#--------------------------------
# make regions
lowtiles = tileMethylCounts(SVAL_myobj,win.size=1000,step.size=1000,cov.bases = 10)
head(lowtiles[[1]],3)
save(lowtiles, file = "SVAL_MethylKit_tiles_lowcov_CHH_rand.RData")
load("SVAL_MethylKit_tiles_lowcov_CHH_rand.RData")

#------------------------
methtiles=unite(lowtiles, destrand=FALSE)
save(methtiles, file = "SVAL_united_lowtiles_CHH_rand.RData")
load("SVAL_united_lowtiles_CHH_rand.RData")

#---------------------------
# calculate DMRs

myDifftiles=calculateDiffMeth(methtiles, overdispersion="MN", test="Chisq")
save(myDifftiles, file = "SVAL_MethDiffRegions_CHH_rand.RData")
load("SVAL_MethDiffRegions_CHH_rand.RData")

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

# overdisp and chisq
3
0
3

SVAL_W_C_lowcovDMRs <- getData(myDifftiles10p)
write.table(SVAL_W_C_lowcovDMRs, "Methylkit_SVAL_W_C_10_CHH_rand_plowcovDMRs.txt", sep = "\t", quote = FALSE, row.names = FALSE)

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

# none overdisp and chisq


SVAL_W_C_lowcovDMRs <- getData(myDifftiles25p)
write.table(SVAL_W_C_lowcovDMRs, "Methylkit_SVAL_W_C_25_CHH_rand_plowcovDMRs.txt", sep = "\t", quote = FALSE, row.names = FALSE)

############################################
# DMRs for Sweden

cd /lustre04/scratch/celphin/Dryas_large_folders/CHG_CHH/data/coverage/CHH_cytosine_report

tmux new-session -s methwarm4
tmux attach-session -t methwarm4

salloc -c1 --time 3:00:00 --mem 120000m --account def-henryg

module load StdEnv/2023
module load r/4.4.0

R

library(methylKit)

# Set the working directory to the folder containing the text files
setwd("/lustre04/scratch/celphin/Dryas_large_folders/CHG_CHH/data/coverage/CHH_cytosine_report")

# List all text files in the folder with a specific pattern 
files <- list.files(pattern = "*CX_report.txt$")

# Define sample ids (replace with your actual sample names)
sample_ids <- sub("_Non_CpG\\.report\\.txt\\.CX_report\\.txt$", "", files)

covariate <- sub("^([^\\.]*\\.[^\\.]*\\.[^\\.]*\\.([A-Za-z]{3}).*)", "\\2", sample_ids)


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

rand_d.f <- as.data.frame(cbind(sample_ids, treatment, covariate))

# subset for Sweden
LAT_rand_d.f <- rand_d.f[which(rand_d.f$covariate=="LAT"),]
LAT_files <- files[grep("LAT",files)]

# randomly rearrange treatment
treatment0 <- sample(LAT_rand_d.f$treatment)
LAT_rand_d.f$treatment0 <- treatment0 

# Read all the files into methylRawList using methRead
LAT_myobj <- methRead(
  location = as.list(LAT_files), 
  sample.id = as.list(LAT_rand_d.f$sample_ids), 
  assembly = "Dryas",
  treatment = as.numeric(LAT_rand_d.f$treatment0),
  dbdir = getwd(),
  context ="CHH",
  dbtype = "tabix",
  resolution = "base", 
  mincov = 3,
  pipeline = "bismarkCytosineReport" 
)

# Print the object to check the result
#print(LAT_myobj)
save(LAT_myobj, file = "LAT_MethylKit_lowcov_data_CHH_rand.RData")
load("LAT_MethylKit_lowcov_data_CHH_rand.RData")

meth=unite(LAT_myobj, destrand=FALSE)
save(meth, file = "LAT_united_CHH_rand.RData")

#---------------------------
# make regions
lowtiles = tileMethylCounts(LAT_myobj,win.size=1000,step.size=1000,cov.bases = 10)
head(lowtiles[[1]],3)
save(lowtiles, file = "LAT_MethylKit_tiles_lowcov_CHH_rand.RData")
load("LAT_MethylKit_tiles_lowcov_CHH_rand.RData")

#------------------------
methtiles=unite(lowtiles, destrand=FALSE)
save(methtiles, file = "LAT_united_lowtiles_CHH_rand.RData")
load("LAT_united_lowtiles_CHH_rand.RData")

#---------------------------
# calculate DMRs

myDifftiles=calculateDiffMeth(methtiles, overdispersion="MN", test="Chisq")
save(myDifftiles, file = "LAT_MethDiffRegions_CHH_rand.RData")
load("LAT_MethDiffRegions_CHH_rand.RData")

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

# none -overdisp and Chisq


LAT_W_C_lowcovDMRs <- getData(myDifftiles10p)
write.table(LAT_W_C_lowcovDMRs, "Methylkit_LAT_W_C_10_CHH_rand_plowcovDMRs.txt", sep = "\t", quote = FALSE, row.names = FALSE)

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


# none overdisp and Chisq


LAT_W_C_lowcovDMRs <- getData(myDifftiles25p)
write.table(LAT_W_C_lowcovDMRs, "Methylkit_LAT_W_C_25_CHH_rand_plowcovDMRs.txt", sep = "\t", quote = FALSE, row.names = FALSE)

########################
# Try 300bp tiles

#---------------------------
# make regions
lowtiles = tileMethylCounts(LAT_myobj,win.size=300,step.size=300,cov.bases = 10)
head(lowtiles[[1]],3)
save(lowtiles, file = "LAT_MethylKit_tiles_lowcov_CHH_rand300.RData")
load("LAT_MethylKit_tiles_lowcov_CHH_rand300.RData")

#------------------------
methtiles=unite(lowtiles, destrand=FALSE)
save(methtiles, file = "LAT_united_lowtiles_CHH_rand300.RData")
load("LAT_united_lowtiles_CHH_rand300.RData")

#---------------------------
# calculate DMRs

myDifftiles=calculateDiffMeth(methtiles, overdispersion="MN", test="Chisq")
save(myDifftiles, file = "LAT_MethDiffRegions_CHH_rand300.RData")
load("LAT_MethDiffRegions_CHH_rand300.RData")

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

# overdisp and Chisq
# 1
# 1
# 0

LAT_W_C_lowcovDMRs <- getData(myDifftiles10p)
write.table(LAT_W_C_lowcovDMRs, "Methylkit_LAT_W_C_10_CHH_rand_plowcovDMRs_300.txt", sep = "\t", quote = FALSE, row.names = FALSE)

########################################
# DMRs for Nunavut

cd /lustre04/scratch/celphin/Dryas_large_folders/CHG_CHH/data/coverage/CHH_cytosine_report

tmux new-session -s methwarm
tmux attach-session -t methwarm

salloc -c1 --time 7:00:00 --mem 120000m --account def-henryg

module load StdEnv/2023
module load r/4.4.0

R

library(methylKit)

# Set the working directory to the folder containing the text files
setwd("/lustre04/scratch/celphin/Dryas_large_folders/CHG_CHH/data/coverage/CHH_cytosine_report")

# List all text files in the folder with a specific pattern 
files <- list.files(pattern = "*CX_report.txt$")

# Define sample ids (replace with your actual sample names)
sample_ids <- sub("_Non_CpG\\.report\\.txt\\.CX_report\\.txt$", "", files)

covariate <- sub("^([^\\.]*\\.[^\\.]*\\.[^\\.]*\\.([A-Za-z]{3}).*)", "\\2", sample_ids)

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

rand_d.f <- as.data.frame(cbind(sample_ids, treatment, covariate))

Alex_rand_d.f <- rand_d.f[c(which(rand_d.f$covariate=="CAS"), which(rand_d.f$covariate=="WIL"), which(rand_d.f$covariate=="Dry"), which(rand_d.f$covariate=="MEA"), which(rand_d.f$covariate=="FER")),]
Alex_files <- files[c(grep("CAS",files), grep("WIL",files), grep("Dry",files), grep("MEA",files), grep("FER",files))]

# randomly rearrange treatment
treatment0 <- sample(Alex_rand_d.f$treatment)
Alex_rand_d.f$treatment0 <- treatment0 

# Read all the files into methylRawList using methRead
Alex_myobj <- methRead(
  location = as.list(Alex_files), 
  sample.id = as.list(Alex_rand_d.f$sample_ids), 
  assembly = "Dryas",
  treatment = as.numeric(Alex_rand_d.f$treatment0),
  dbdir = getwd(),
  context ="CHH",
  dbtype = "tabix",
  resolution = "base", 
  mincov = 3,
  pipeline = "bismarkCytosineReport" 
)

# Print the object to check the result
#print(Alex_myobj)
save(Alex_myobj, file = "Alex_MethylKit_lowcov_data_CHH_rand.RData")
load("Alex_MethylKit_lowcov_data_CHH_rand.RData")

meth=unite(Alex_myobj, destrand=FALSE)
save(meth, file = "Alex_united_CHH_rand.RData")


#------------------------------
# make regions
lowtiles = tileMethylCounts(Alex_myobj,win.size=1000,step.size=1000,cov.bases = 10)
head(lowtiles[[1]],3)
save(lowtiles, file = "Alex_MethylKit_tiles_lowcov_CHH_rand.RData")
load("Alex_MethylKit_tiles_lowcov_CHH_rand.RData")

#------------------------
methtiles=unite(lowtiles, destrand=FALSE)
save(methtiles, file = "Alex_united_lowtiles_CHH_rand.RData")
load("Alex_united_lowtiles_CHH_rand.RData")

#---------------------------
# calculate DMRs

myDifftiles=calculateDiffMeth(methtiles, overdispersion="MN", test="Chisq")
save(myDifftiles, file = "Alex_MethDiffRegions_CHH_rand.RData")
load("Alex_MethDiffRegions_CHH_rand.RData")

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

# none overdisp and chisq

Alex_W_C_lowcovDMRs <- getData(myDifftiles10p)
write.table(Alex_W_C_lowcovDMRs, "Methylkit_Alex_W_C_10_CHH_rand_plowcovDMRs.txt", sep = "\t", quote = FALSE, row.names = FALSE)

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

# none overdisp and chisq


Alex_W_C_lowcovDMRs <- getData(myDifftiles25p)
write.table(Alex_W_C_lowcovDMRs, "Methylkit_Alex_W_C_25_CHH_rand_plowcovDMRs.txt", sep = "\t", quote = FALSE, row.names = FALSE)

######################################
##################################
# Compare subsites methylation to eachother
# not run - need to edit the file names 


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


#####################################
# back in bash compare with metilene and overlay with genome annotations
# see bedGraph intersection notes