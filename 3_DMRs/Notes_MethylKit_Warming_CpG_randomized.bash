###########################
# methylKit
# https://github.com/al2na/methylKit
# Dec 2024
##############################

# Randomize methylkit

cd /lustre04/scratch/celphin/Dryas_large_folders/CpG/

tmux new-session -s methwarm
tmux attach-session -t methwarm

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

rand_d.f <- as.data.frame(cbind(sample_ids, treatment, covariate))

# randomly rearrange treatment
treatment0 <- sample(treatment)

rand_d.f$treatment0 <- treatment0 

load("united_lowtiles.RData")

methtiles_rand <- reorganize(
  methtiles,
  as.list(rand_d.f$sample_ids),
  as.numeric(rand_d.f$treatment0),
  chunk.size = 1e+06,
  save.db = TRUE,
  suffix = "rand",
  dbtype = "tabix"
)

#---------------------------
# calculate DMRs

myDifftiles=calculateDiffMeth(methtiles, covariates=rand_d.f$covariates,overdispersion="MN", test="Chisq")

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

# with covariates, overdisp, chisq  - none for rand or real


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



##########################################
# DMRs for Alaska

library(methylKit)

# Set the working directory to the folder containing the text files
setwd("/lustre04/scratch/celphin/Dryas_large_folders/CpG/Wild_stranded_CpG_report")

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

myDifftiles=calculateDiffMeth(methtiles_rand,overdispersion="MN", test="Chisq")
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

# rand over disp and chisq
# [1] 6
# [1] 2
# [1] 4

# rand1 
# [1] 146
# [1] 89
# [1] 57

# rand2
# 7,4,3 

# rand 3 3 only moved 3 individuals between W/C
# [1] 55
# [1] 31
# [1] 24


# real overdip and chisq
# [1] 158
# [1] 80
# [1] 78

ALAS_rand_W_C_lowcovDMRs <- getData(myDifftiles10p)
write.table(ALAS_rand_W_C_lowcovDMRs, "Methylkit_ALAS_rand3_W_C_10plowcovDMRs.txt", sep = "\t", quote = FALSE, row.names = FALSE)

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

# real overdisp and chisq
# [1] 56
# [1] 19
# [1] 37

# rand over disp and chisq - none

# rand1 
# [1] 13
# [1] 7
# [1] 6

# rand2
# 1,0,1

# rand3
# [1] 9
# [1] 7
# [1] 2


ALAS_rand_W_C_lowcovDMRs <- getData(myDifftiles25p)
write.table(ALAS_rand_W_C_lowcovDMRs, "Methylkit_ALAS_rand3_W_C_25plowcovDMRs.txt", sep = "\t", quote = FALSE, row.names = FALSE)

######################
#---------------------------
# calculate overdispDMRs

myDifftiles=calculateDiffMeth(methtiles_rand)
save(myDifftiles, file = "ALAS_rand_MethDiffRegions_overdisp.RData")
load("ALAS_rand_MethDiffRegions_overdisp.RData")

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

# rand3 overdisp
# [1] 4581
# [1] 2267
# [1] 2314


ALAS_rand_W_C_lowcovDMRs <- getData(myDifftiles10p)
write.table(ALAS_rand_W_C_lowcovDMRs, "Methylkit_ALAS_rand_W_C_10_overdisp_DMRs.txt", sep = "\t", quote = FALSE, row.names = FALSE)

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

#rand3 overdisp
# [1] 119
# [1] 72
# [1] 47


ALAS_rand_W_C_lowcovDMRs <- getData(myDifftiles25p)
write.table(ALAS_rand_W_C_lowcovDMRs, "Methylkit_ALAS_rand_W_C_25_overdisp_DMRs.txt", sep = "\t", quote = FALSE, row.names = FALSE)

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

myDifftiles=calculateDiffMeth(methtiles_rand,overdispersion="MN", test="Chisq")
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

# rand over disp and chisq
# [1] 128
# [1] 59
# [1] 69

# rand1
[1] 80
[1] 40
[1] 40



SVAL_rand_W_C_lowcovDMRs <- getData(myDifftiles10p)
write.table(SVAL_rand_W_C_lowcovDMRs, "Methylkit_SVAL_rand1_W_C_10plowcovDMRs.txt", sep = "\t", quote = FALSE, row.names = FALSE)

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

# rand over disp and chisq
# [1] 14
# [1] 9
# [1] 5

# rand1 
[1] 16
[1] 8
[1] 8


SVAL_rand_W_C_lowcovDMRs <- getData(myDifftiles25p)
write.table(SVAL_rand_W_C_lowcovDMRs, "Methylkit_SVAL_rand1_W_C_25plowcovDMRs.txt", sep = "\t", quote = FALSE, row.names = FALSE)

#---------------------------
# calculate overdisp DMRs

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

# rand1 overdisp

SVAL_rand_W_C_lowcovDMRs <- getData(myDifftiles10p)
write.table(SVAL_rand_W_C_lowcovDMRs, "Methylkit_SVAL_rand_W_C_10_overdisp_DMRs.txt", sep = "\t", quote = FALSE, row.names = FALSE)

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

# rand1 overdisp
[1] 274
[1] 149
[1] 125


SVAL_rand_W_C_lowcovDMRs <- getData(myDifftiles25p)
write.table(SVAL_rand_W_C_lowcovDMRs, "Methylkit_SVAL_rand_W_C_25_overdisp_DMRs.txt", sep = "\t", quote = FALSE, row.names = FALSE)

#######################################
# Sweden

library(methylKit)

# Set the working directory to the folder containing the text files
setwd("/lustre04/scratch/celphin/Dryas_large_folders/CpG/Wild_stranded_CpG_report")

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

myDifftiles=calculateDiffMeth(methtiles_rand,overdispersion="MN", test="Chisq")
save(myDifftiles, file = "LAT_rand_MethDiffRegions.RData")
load("LAT_rand_MethDiffRegions.RData")

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

# [1] 1597
# [1] 585
# [1] 1012

# rand over disp and chisq
# 1,1,0

# rand1
# 6,5,1

# rand 2
0

# rand3
[1] 26
[1] 20
[1] 6


LAT_rand_W_C_lowcovDMRs <- getData(myDifftiles10p)
write.table(LAT_rand_W_C_lowcovDMRs, "Methylkit_LAT_rand3_W_C_10plowcovDMRs.txt", sep = "\t", quote = FALSE, row.names = FALSE)

#--------------------------
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

# [1] 20
# [1] 5
# [1] 15

# rand over disp and chisq - none

# rand1
# 1 1 0 

# rand 2
0

# rand3
# 2,2,0

LAT_rand_W_C_lowcovDMRs <- getData(myDifftiles25p)
write.table(LAT_rand_W_C_lowcovDMRs, "Methylkit_LAT_rand3_W_C_25plowcovDMRs.txt", sep = "\t", quote = FALSE, row.names = FALSE)

##############################

# calculate overdisp DMRs

myDifftiles=calculateDiffMeth(methtiles_rand)
save(myDifftiles, file = "LAT_rand_MethDiffRegions.RData")
load("LAT_rand_MethDiffRegions.RData")

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

# overdisp
[1] 4327
[1] 2043
[1] 2284



LAT_rand_W_C_lowcovDMRs <- getData(myDifftiles10p)
write.table(LAT_rand_W_C_lowcovDMRs, "Methylkit_LAT_rand_W_C_10_overdisp_DMRs.txt", sep = "\t", quote = FALSE, row.names = FALSE)

#--------------------------
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

# overdisp 
[1] 118
[1] 56
[1] 62


LAT_rand_W_C_lowcovDMRs <- getData(myDifftiles25p)
write.table(LAT_rand_W_C_lowcovDMRs, "Methylkit_LAT_rand_W_C_25_overdisp_DMRs.txt", sep = "\t", quote = FALSE, row.names = FALSE)

####################################
#Nunavut

library(methylKit)

# Set the working directory to the folder containing the text files
setwd("/lustre04/scratch/celphin/Dryas_large_folders/CpG/Wild_stranded_CpG_report")

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

myDifftiles=calculateDiffMeth(methtiles_rand,overdispersion="MN", test="Chisq")
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

# overdisp and chisq - none 

Alex_rand_W_C_lowcovDMRs <- getData(myDifftiles10p)
write.table(Alex_rand_W_C_lowcovDMRs, "Methylkit_Alex_rand1_W_C_10plowcovDMRs.txt", sep = "\t", quote = FALSE, row.names = FALSE)

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

# overdisp and chisq - none 
# rand1 - none

Alex_rand_W_C_lowcovDMRs <- getData(myDifftiles25p)
write.table(Alex_rand_W_C_lowcovDMRs, "Methylkit_Alex_rand1_W_C_25plowcovDMRs.txt", sep = "\t", quote = FALSE, row.names = FALSE)

#---------------------------
# calculate overdisp DMRs

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

# rand
# [1] 707
# [1] 448
# [1] 259

# rand1
# [1] 782
# [1] 404
# [1] 378


# overdisp and chisq - none 

Alex_rand_W_C_lowcovDMRs <- getData(myDifftiles10p)
write.table(Alex_rand_W_C_lowcovDMRs, "Methylkit_Alex_rand_W_C_10_overdisp_DMRs.txt", sep = "\t", quote = FALSE, row.names = FALSE)

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

# rand
# 3
# 3
# 0

# rand1
# [1] 4
# [1] 2
# [1] 2



Alex_rand_W_C_lowcovDMRs <- getData(myDifftiles25p)
write.table(Alex_rand_W_C_lowcovDMRs, "Methylkit_Alex_rand_W_C_25_overdisp_DMRs.txt", sep = "\t", quote = FALSE, row.names = FALSE)

##################################
# Compare High and Low Arctic


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

HL_rand_d.f <- d.f[c(which(d.f$covariate=="CAS"), 
which(d.f$covariate=="WIL"), 
which(d.f$covariate=="DRY"), 
which(d.f$covariate=="MEA"), 
which(d.f$covariate=="FER"),
which(d.f$covariate=="SVA"),
which(d.f$covariate=="ALA"),
which(d.f$covariate=="LAT")),]

load("MethylKit_lowcov_data.RData")

# randomly rearrange treatment
treatment0 <- sample(HL_rand_d.f$treatment)
HL_rand_d.f$treatment0 <- treatment0 


HL_rand_myobj <- reorganize(
  myobj_lowCov,
  as.list(HL_rand_d.f$sample_ids),
  as.numeric(HL_rand_d.f$treatment0),
  chunk.size = 1e+06,
  save.db = TRUE,
  suffix = "HL",
  dbtype = "tabix"
)

# Print the object to check the result
#print(HL_rand_myobj)
save(HL_rand_myobj, file = "HL_rand_MethylKit_lowcov_data.RData")
load("HL_rand_MethylKit_lowcov_data.RData")

#------------------------------
# make regions
lowtiles = tileMethylCounts(HL_rand_myobj,win.size=1000,step.size=1000,cov.bases = 10)
head(lowtiles[[1]],3)
save(lowtiles, file = "HL_rand_MethylKit_tiles_lowcov.RData")
load("HL_rand_MethylKit_tiles_lowcov.RData")

#------------------------
methtiles=unite(lowtiles, destrand=FALSE)
save(methtiles, file = "HL_rand_united_lowtiles.RData")
load("HL_rand_united_lowtiles.RData")

#---------------------------
# calculate DMRs

myDifftiles=calculateDiffMeth(methtiles, covariates=HL_rand_d.f$covariates, overdispersion="MN",test="Chisq")
save(myDifftiles, file = "HL_rand_MethDiffRegions.RData")
load("HL_rand_MethDiffRegions.RData")

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

# 1
# 1
# 0


HL_rand_lowcovDMRs <- getData(myDifftiles10p)
write.table(HL_rand_lowcovDMRs, "Methylkit_HL_rand_10plowcovDMRs.txt", sep = "\t", quote = FALSE, row.names = FALSE)

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

# 1
# 1
# 0


HL_rand_lowcovDMRs <- getData(myDifftiles25p)
write.table(HL_rand_lowcovDMRs, "Methylkit_HL_rand_25plowcovDMRs.txt", sep = "\t", quote = FALSE, row.names = FALSE)

######################################
# Cassiope/Dryas/Willow/Fert/Meadow sites separately randomized

tmux new-session -s methylkit3
tmux attach-session -t methylkit3

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

d.f <- as.data.frame(cbind(sample_ids, covariate))

CASS_rand_d.f <- d.f[which(d.f$covariate=="CAS"),]

# Extract first letter (C or W)
first_letter <- substr(CASS_rand_d.f$sample_ids, 1, 1)

# Replace C with 0 and W with 1
treatment <- gsub("C", "0", first_letter)
treatment <- gsub("W", "1", treatment)

# randomly rearrange treatment
treatment0 <- sample(treatment)

CASS_rand_d.f$treatment0 <- treatment0 

load("CASS_united_lowtiles.RData")

methtiles_rand <- reorganize(
  methtiles,
  as.list(CASS_rand_d.f$sample_ids),
  as.numeric(CASS_rand_d.f$treatment0),
  chunk.size = 1e+06,
  save.db = TRUE,
  suffix = "Alex_rand",
  dbtype = "tabix"
)


#---------------------------
# calculate DMRs

myDifftiles=calculateDiffMeth(methtiles_rand, overdispersion="MN", test="Chisq")
save(myDifftiles, file = "CASS_rand_MethDiffRegions.RData")
load("CASS_rand_MethDiffRegions.RData")

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

# rand
# [1] 4850
# [1] 2965
# [1] 1885

# real W/C
# [1] 4962
# [1] 2739
# [1] 2223

# overdispersion and chisq
# [1] 100
# [1] 70
# [1] 30


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

CASS_rand_W_C_lowcovDMRs <- getData(myDifftiles25p)
write.table(CASS_rand_W_C_lowcovDMRs, "Methylkit_CASS_rand_W_C_25plowcovDMRs.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# rand
# [1] 149
# [1] 84
# [1] 65

# real W/C
# [1] 187
# [1] 101
# [1] 86

# overdispersion and chisq
[1] 13
[1] 7
[1] 6


#####################################
# DRY
library(methylKit)

# Set the working directory to the folder containing the text files
setwd("/lustre04/scratch/celphin/Dryas_large_folders/CpG/Wild_stranded_CpG_report")

# List all text files in the folder with a specific pattern 
files <- list.files(pattern = "*CpG_report.txt$")

# Define sample ids (replace with your actual sample names)
sample_ids <- sub("_R1_val_1_bismark_bt2_pe\\.deduplicated\\.CpG_report\\.txt$", "", files)

covariate <- sub("^[^_]*_([A-Za-z]{3}).*", "\\1", sample_ids)

d.f <- as.data.frame(cbind(sample_ids, covariate))

DRY_rand_d.f <- d.f[which(d.f$covariate=="DRY"),]

# Extract first letter (C or W)
first_letter <- substr(DRY_rand_d.f$sample_ids, 1, 1)

# Replace C with 0 and W with 1
treatment <- gsub("C", "0", first_letter)
treatment <- gsub("W", "1", treatment)

# randomly rearrange treatment
treatment0 <- sample(treatment)

DRY_rand_d.f$treatment0 <- treatment0 

load("DRY_united_lowtiles.RData")

methtiles_rand <- reorganize(
  methtiles,
  as.list(DRY_rand_d.f$sample_ids),
  as.numeric(DRY_rand_d.f$treatment0),
  chunk.size = 1e+06,
  save.db = TRUE,
  suffix = "Alex_rand",
  dbtype = "tabix"
)


#---------------------------
# calculate DMRs

myDifftiles=calculateDiffMeth(methtiles_rand)
save(myDifftiles, file = "DRY_rand_MethDiffRegions.RData")
load("DRY_rand_MethDiffRegions.RData")

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

[1] 5932
[1] 2888
[1] 3044

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

[1] 221
[1] 96
[1] 125


DRY_rand_W_C_lowcovDMRs <- getData(myDifftiles25p)
write.table(DRY_rand_W_C_lowcovDMRs, "Methylkit_DRY_rand_W_C_25plowcovDMRs.txt", sep = "\t", quote = FALSE, row.names = FALSE)

#####################################
# MEAD
library(methylKit)

# Set the working directory to the folder containing the text files
setwd("/lustre04/scratch/celphin/Dryas_large_folders/CpG/Wild_stranded_CpG_report")

# List all text files in the folder with a specific pattern 
files <- list.files(pattern = "*CpG_report.txt$")

# Define sample ids (replace with your actual sample names)
sample_ids <- sub("_R1_val_1_bismark_bt2_pe\\.deduplicated\\.CpG_report\\.txt$", "", files)

covariate <- sub("^[^_]*_([A-Za-z]{3}).*", "\\1", sample_ids)

d.f <- as.data.frame(cbind(sample_ids, covariate))

MEAD_rand_d.f <- d.f[which(d.f$covariate=="MEA"),]

# Extract first letter (C or W)
first_letter <- substr(MEAD_rand_d.f$sample_ids, 1, 1)

# Replace C with 0 and W with 1
treatment <- gsub("C", "0", first_letter)
treatment <- gsub("W", "1", treatment)

# randomly rearrange treatment
treatment0 <- sample(treatment)

MEAD_rand_d.f$treatment0 <- treatment0 

load("MEAD_united_lowtiles.RData")

methtiles_rand <- reorganize(
  methtiles,
  as.list(MEAD_rand_d.f$sample_ids),
  as.numeric(MEAD_rand_d.f$treatment0),
  chunk.size = 1e+06,
  save.db = TRUE,
  suffix = "Alex_rand",
  dbtype = "tabix"
)


#---------------------------
# calculate DMRs

myDifftiles=calculateDiffMeth(methtiles_rand)
save(myDifftiles, file = "MEAD_rand_MethDiffRegions.RData")
load("MEAD_rand_MethDiffRegions.RData")

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

[1] 8334
[1] 4572
[1] 3762


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

[1] 423
[1] 243
[1] 180


MEAD_rand_W_C_lowcovDMRs <- getData(myDifftiles25p)
write.table(MEAD_rand_W_C_lowcovDMRs, "Methylkit_MEAD_rand_W_C_25plowcovDMRs.txt", sep = "\t", quote = FALSE, row.names = FALSE)


#####################################
# WILL
library(methylKit)

# Set the working directory to the folder containing the text files
setwd("/lustre04/scratch/celphin/Dryas_large_folders/CpG/Wild_stranded_CpG_report")

# List all text files in the folder with a specific pattern 
files <- list.files(pattern = "*CpG_report.txt$")

# Define sample ids (replace with your actual sample names)
sample_ids <- sub("_R1_val_1_bismark_bt2_pe\\.deduplicated\\.CpG_report\\.txt$", "", files)

covariate <- sub("^[^_]*_([A-Za-z]{3}).*", "\\1", sample_ids)

d.f <- as.data.frame(cbind(sample_ids, covariate))

WILL_rand_d.f <- d.f[which(d.f$covariate=="WIL"),]

# Extract first letter (C or W)
first_letter <- substr(WILL_rand_d.f$sample_ids, 1, 1)

# Replace C with 0 and W with 1
treatment <- gsub("C", "0", first_letter)
treatment <- gsub("W", "1", treatment)

# randomly rearrange treatment
treatment0 <- sample(treatment)

WILL_rand_d.f$treatment0 <- treatment0 

load("WILL_united_lowtiles.RData")

methtiles_rand <- reorganize(
  methtiles,
  as.list(WILL_rand_d.f$sample_ids),
  as.numeric(WILL_rand_d.f$treatment0),
  chunk.size = 1e+06,
  save.db = TRUE,
  suffix = "Alex_rand",
  dbtype = "tabix"
)


#---------------------------
# calculate DMRs

myDifftiles=calculateDiffMeth(methtiles_rand)
save(myDifftiles, file = "WILL_rand_MethDiffRegions.RData")
load("WILL_rand_MethDiffRegions.RData")

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

# [1] 7840
# [1] 3550
# [1] 4290

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

# [1] 405
# [1] 180
# [1] 225


WILL_rand_W_C_lowcovDMRs <- getData(myDifftiles25p)
write.table(WILL_rand_W_C_lowcovDMRs, "Methylkit_WILL_rand_W_C_25plowcovDMRs.txt", sep = "\t", quote = FALSE, row.names = FALSE)


###############################
# FERT
library(methylKit)

# Set the working directory to the folder containing the text files
setwd("/lustre04/scratch/celphin/Dryas_large_folders/CpG/Wild_stranded_CpG_report")

# List all text files in the folder with a specific pattern 
files <- list.files(pattern = "*CpG_report.txt$")

# Define sample ids (replace with your actual sample names)
sample_ids <- sub("_R1_val_1_bismark_bt2_pe\\.deduplicated\\.CpG_report\\.txt$", "", files)

covariate <- sub("^[^_]*_([A-Za-z]{3}).*", "\\1", sample_ids)

d.f <- as.data.frame(cbind(sample_ids, covariate))

FERT_rand_d.f <- d.f[which(d.f$covariate=="FER"),]

# Extract first letter (C or W)
first_letter <- substr(FERT_rand_d.f$sample_ids, 1, 1)

# Replace C with 0 and W with 1
treatment <- gsub("C", "0", first_letter)
treatment <- gsub("W", "1", treatment)

# randomly rearrange treatment
treatment0 <- sample(treatment)

FERT_rand_d.f$treatment0 <- treatment0 

load("FERT_united_lowtiles.RData")

methtiles_rand <- reorganize(
  methtiles,
  as.list(FERT_rand_d.f$sample_ids),
  as.numeric(FERT_rand_d.f$treatment0),
  chunk.size = 1e+06,
  save.db = TRUE,
  suffix = "Alex_rand",
  dbtype = "tabix"
)


#---------------------------
# calculate DMRs

myDifftiles=calculateDiffMeth(methtiles_rand)
save(myDifftiles, file = "FERT_rand_MethDiffRegions.RData")
load("FERT_rand_MethDiffRegions.RData")

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

# [1] 8801
# [1] 3858
# [1] 4943

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

# [1] 447
# [1] 195
# [1] 252

FERT_rand_W_C_lowcovDMRs <- getData(myDifftiles25p)
write.table(FERT_rand_W_C_lowcovDMRs, "Methylkit_FERT_rand_W_C_25plowcovDMRs.txt", sep = "\t", quote = FALSE, row.names = FALSE)

#####################################
# back in bash compare with metilene and overlay with genome annotations
# see bedGraph intersection notes