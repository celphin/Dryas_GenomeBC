########################################
# Dryas differential gene expression data
# Determining differentially expressed genes, plotting results
# Nov 2024
#####################################
# following https://bioconductor.org/packages/release/workflows/vignettes/RnaSeqGeneEdgeRQL/inst/doc/edgeRQL.html
###########################################################
# try differential expression analysis - plotting
# EdgeR

tmux new-session -s RNA
tmux attach-session -t RNA

# Cedar1
cd /home/celphin/projects/rrg-rieseber-ac/rpp-rieseber/celphin/Dryas/RNAseq_analysis/

# copy to Beluga scratch/
cd /home/celphin/scratch/Dryas/RNAseq_analysis/

#---------------------------------
# https://github.com/owensgl/biol525D/tree/master/Topic_6
# https://bioconductor.org/packages/release/bioc/html/edgeR.html
# http://combine-australia.github.io/RNAseq-R/slides/RNASeq_filtering_qc.pdf 
# also see https://rpubs.com/jrgonzalezISGlobal/enrichment
# https://bioinformatics-core-shared-training.github.io/cruk-bioinf-sschool/Day3/rnaSeq_DE.pdf


# get R on server
# open R
module load StdEnv/2020
module load r/4.2.1
module load gdal
module load udunits
module load python
#export R_LIBS_USER=/home/msandler/R/x86_64-pc-linux-gnu-library/4.2.1/
export R_LIBS_USER=/home/celphin/R/x86_64-pc-linux-gnu-library/4.2.1/

R

#install the package edgeR
 if (!require("BiocManager", quietly = TRUE))
     install.packages("BiocManager")
 BiocManager::install(version = "3.15")

 BiocManager::install("edgeR")
 install.packages('gplots')
 install.packages('stringr')
 
#paste it in here (i.e. replace my path with yours):
#setwd ("/scratch/msandler/RNAseq_analysis")
setwd ("/home/celphin/scratch/Dryas/RNAseq_analysis/")

#find your working directory:
getwd()

#load the libraries you will need 
library(edgeR)
library(gplots)
library(stringr)

#-----------------------
#read in the data
mydata <- read.table("./RNA_site_tables/gene_names_expression_table1.txt", header=TRUE)

CDS <- mydata[,1]
exp_data <- mydata[,c(2:length(mydata))]

#make a DGElist object out of the data
dge <- DGEList (counts = exp_data, genes=rownames(exp_data))

#--------------------------
# NORMALIZATION

# calculate the normalization factors to adjust the effective library size relative to 
# other libraries in the dataset (norm.factor that minimizes log fold change among samples)
# https://www.rdocumentation.org/packages/edgeR/versions/3.14.0/topics/calcNormFactors
dge_norm <- calcNormFactors(dge, method="TMM")

# Normalize raw counts to CPM (Counts Per Million)
# https://www.rdocumentation.org/packages/edgeR/versions/3.14.0/topics/cpm
norm_cpm <- cpm(dge_norm, normalized.lib.sizes = TRUE)

#-------------------------
# Setup study design matrix

dge_norm$samples$group <-  as.factor(str_extract(rownames(dge_norm$samples), "^[WC]"))
dge_norm$samples$site <-  as.factor(str_extract(rownames(dge_norm$samples), "(?<=_)[A-Za-z]+"))

designMat <- model.matrix( ~ 0 + group + site, data=dge_norm$samples)

#---------------------------
# FILTERING 

#filter out the genes with less than 1 cpm in 6 or fewer libraries (a somewhat arbitrary choice). 
#Genes are usually dropped if they can't possibly be expressed in all the samples for any of the conditions.
# filtered_dge <- dge.filt[rowSums(norm_cpm > 1) >= 6,]
# cpm.filtered_dge <- cpm.dge[rowSums(cpm.dge > 1) >= 6,]

# https://rdrr.io/bioc/edgeR/man/filterByExpr.html
# filter and see dimensions
keep.exprs <- filterByExpr(dge_norm, design = designMat, min.count = 10, min.total.count = 15, large.n = 10, min.prop = 0.7)
dge.filt <- dge_norm[keep.exprs,]
dim(dge.filt)

#[1] 23222    69
# 8*8+6 (missing one Swed control)

#------------------------
# PLOTTING 

# Before normalization
png("./plots/Before_normalization_counts.png", width = 700, height = 500)
boxplot(exp_data, main = "Raw Counts", las=2)
dev.off()

# After normalization (using CPM)
png("./plots/After_normalization_counts.png", width = 700, height = 500)
boxplot(norm_cpm, main = "Normalized Counts", las=2)
dev.off()

# plot mean variance
png("./plots/voom_mean_variance_trend.png", width = 700, height = 500)
v <- voom(dge.filt, design = designMat, plot = TRUE)
dev.off()

#generate a multi-dimensional scaling plot
col_treat <- as.character (dge.filt$samples$group)
col_treat [col_treat == "C"] <- "blue"
col_treat [col_treat == "W"] <- "red"

col_site <- as.character (dge.filt$samples$site)
col_site [col_treat == "Alex"] <- 15 # square
col_site [col_treat == "Norw"] <- 16 # circle
col_site [col_treat == "Alas"] <- 17 # triangle
col_site [col_treat == "Swed"] <- 18 # diamond

#col_treat [col_treat == "H"] <- "green"
#col_treat [col_treat == "L"] <- "yellow"

#plot the MDS graph
png("./plots/W_C_RNA_MDS.png", width = 700, height = 500)
plotMDS (dge.filt, col = col_treat)
dev.off()

png("./plots/Site_RNA_MDS.png", width = 700, height = 500)
plotMDS (dge.filt, col = col_treat, pch=col_site)
dev.off()

#----------------------------
# Estimating dispersion 
dge.filt <- estimateDisp(dge.filt, design=designMat)

#fit the common + tagwise dispersion models
dge.filt <- estimateGLMCommonDisp(dge.filt, design=designMat)
dge.filt <- estimateGLMTrendedDisp(dge.filt, design=designMat)
dge.filt <- estimateGLMTagwiseDisp(dge.filt, design=designMat)

png("./plots/Biological_coeff_variation.png", width = 700, height = 500)
plotBCV(dge.filt)
dev.off()

#-----------------------
# Fitting the model
# https://seqqc.wordpress.com/2020/11/28/10-tips-tricks-for-complex-model-matrix-designs-in-dge-analysis/
# https://support.bioconductor.org/p/9136560/


#Set the model to use. This one includes the intercept, but other models can be 
# specified that omit the intercept or that have more complex designs. See EdgeR manual for details.
fit_QLF <- glmQLFit(dge.filt, design=designMat)

#fit a GLM to the data using the tagwise dispersion estimate
fit <- glmFit(dge.filt, designMat, dispersion = dge.filt$tagwise.dispersion)

# fit linear model to voom
lmfit <-  lmFit(v, designMat)

#---------------------
# Contrasts

fit$design
            # groupC groupW siteAlex siteNorw siteSeed siteSwed
# C_Alas_16        1      0        0        0        0        0
# C_Alas_1         1      0        0        0        0        0
# C_Alas_20        1      0        0        0        0        0
# C_Alas_5         1      0        0        0        0        0
# C_Alex_11        1      0        1        0        0        0


contrasts <- makeContrasts(groupWvsgroupC = groupW-groupC,
                           levels = designMat)

          # Contrasts
# Levels     groupWvsgroupC
  # groupC               -1
  # groupW                1
  # siteAlex              0
  # siteNorw              0
  # siteSeed              0
  # siteSwed              0

#------------------------------
# Performing differential expression analysis
# https://rdrr.io/bioc/edgeR/man/glmfit.html

qlf <- glmQLFTest(fit_QLF, contrast=contrasts)
lrt <- glmLRT(fit, contrast=contrasts)


fit2 <- contrasts.fit(lmfit, contrasts)
fit2 <- eBayes(fit2)

#--------------------
# Extract the top differentially expressed genes
#?topTags

WC_result_eBayes = topTable(fit2, n=Inf, coef = "groupWvsgroupC", p.value = 0.05)
                                         # genes      logFC    AveExpr         t
# Do1_05_a00001G00678V1.1 Do1_05_a00001G00678V1.1  0.4082366  4.9961193  5.189187
# Do1_04_a00002G00429V1.1 Do1_04_a00002G00429V1.1 -1.6381272 -1.0473707 -5.457280
# Do1_01_a00007G00164V1.1 Do1_01_a00007G00164V1.1  2.5823275 -2.5398183  5.321775
# Do1_03_a00002G01336V1.1 Do1_03_a00002G01336V1.1 -0.9156573  0.7377493 -4.973376
# Do1_03_a00001G00267V1.1 Do1_03_a00001G00267V1.1 -1.6954531  2.8551757 -4.911546
# Do1_02_a00003G01059V1.1 Do1_02_a00003G01059V1.1 -0.3332865  5.7836771 -4.834632
                             # P.Value  adj.P.Val        B
# Do1_05_a00001G00678V1.1 3.089154e-06 0.02339931 4.388781
# Do1_04_a00002G00429V1.1 1.165338e-06 0.02170753 4.201616
# Do1_01_a00007G00164V1.1 1.910538e-06 0.02170753 4.176496
# Do1_03_a00002G01336V1.1 6.702159e-06 0.03795597 3.562598
# Do1_03_a00001G00267V1.1 8.351516e-06 0.03795597 3.476994
# Do1_02_a00003G01059V1.1 1.096676e-05 0.04153479 3.175710


WC_result_LRT <- topTags(lrt)

# Coefficient:  -1*groupC 1*groupW
                                          # genes      logFC   logCPM       LR
# Do1_02_a00003G01208V1.1 Do1_02_a00003G01208V1.1 -1.9316686 3.762947 25.24659
# Do1_05_a00001G00678V1.1 Do1_05_a00001G00678V1.1  0.4122129 5.085566 20.31381
# Do1_01_a00001G01531V1.1 Do1_01_a00001G01531V1.1 -0.7713982 3.073233 18.67167
# Do1_05_a00001G00394V1.1 Do1_05_a00001G00394V1.1  0.5395761 6.089037 18.33425
# Do1_01_a00001G00911V1.1 Do1_01_a00001G00911V1.1  0.6284556 4.263097 18.22913
# Do1_06_a00001G00770V1.1 Do1_06_a00001G00770V1.1 -1.0175513 4.346709 17.44931
# Do1_02_a00003G01059V1.1 Do1_02_a00003G01059V1.1 -0.3331901 5.829629 17.44339
# Do1_02_a00003G01238V1.1 Do1_02_a00003G01238V1.1  1.1042641 2.564388 17.29871
# Do1_05_a00001G01718V1.1 Do1_05_a00001G01718V1.1  0.7143276 2.750032 16.82159
# Do1_03_a00002G01336V1.1 Do1_03_a00002G01336V1.1 -0.8842956 1.093737 16.25318
                              # PValue        FDR
# Do1_02_a00003G01208V1.1 5.044877e-07 0.01146398
# Do1_05_a00001G00678V1.1 6.572432e-06 0.07467597
# Do1_01_a00001G01531V1.1 1.552722e-05 0.08901463
# Do1_05_a00001G00394V1.1 1.853449e-05 0.08901463
# Do1_01_a00001G00911V1.1 1.958604e-05 0.08901463
# Do1_06_a00001G00770V1.1 2.950718e-05 0.09072642
# Do1_02_a00003G01059V1.1 2.959925e-05 0.09072642
# Do1_02_a00003G01238V1.1 3.194030e-05 0.09072642
# Do1_05_a00001G01718V1.1 4.106356e-05 0.10368094
# Do1_03_a00002G01336V1.1 5.541654e-05 0.11382603


WC_result_QLF <- topTags(qlf)
# Coefficient:  -1*groupC 1*groupW
                                          # genes      logFC    logCPM        F
# Do1_05_a00001G00678V1.1 Do1_05_a00001G00678V1.1  0.4122718  5.085566 25.02695
# Do1_02_a00003G01059V1.1 Do1_02_a00003G01059V1.1 -0.3331354  5.829629 23.32719
# Do1_02_a00003G01208V1.1 Do1_02_a00003G01208V1.1 -1.9335646  3.762947 21.84757
# Do1_04_a00001G00224V1.1 Do1_04_a00001G00224V1.1  3.5219559 -2.492429 21.36874
# Do1_01_a00004G01253V1.1 Do1_01_a00004G01253V1.1 -0.3327786  2.967671 20.67348
# Do1_01_a00001G01531V1.1 Do1_01_a00001G01531V1.1 -0.7714805  3.073233 18.80383
# Do1_05_a00001G00394V1.1 Do1_05_a00001G00394V1.1  0.5395776  6.089037 18.52961
# Do1_01_a00001G00911V1.1 Do1_01_a00001G00911V1.1  0.6284348  4.263097 18.45767
# Do1_03_a00001G00904V1.1 Do1_03_a00001G00904V1.1  0.3783682  4.803131 18.07284
# Do1_04_a00001G00673V1.1 Do1_04_a00001G00673V1.1  0.3675828  5.685339 17.86829
                              # PValue       FDR
# Do1_05_a00001G00678V1.1 6.094265e-06 0.1278321
# Do1_02_a00003G01059V1.1 1.125085e-05 0.1278321
# Do1_02_a00003G01208V1.1 1.942058e-05 0.1319802
# Do1_04_a00001G00224V1.1 2.323187e-05 0.1319802
# Do1_01_a00004G01253V1.1 3.020453e-05 0.1372736
# Do1_01_a00001G01531V1.1 6.203528e-05 0.1905427
# Do1_05_a00001G00394V1.1 6.906413e-05 0.1905427
# Do1_01_a00001G00911V1.1 7.104189e-05 0.1905427
# Do1_03_a00001G00904V1.1 8.266954e-05 0.1905427
# Do1_04_a00001G00673V1.1 8.964057e-05 0.1905427


save(WC_result_eBayes, file='edgeR_TopTags_Bayes.RData') #We will need this later
save(WC_result_LRT, file='edgeR_TopTags_LRT.RData') #We will need this later
save(WC_result_QLF, file='edgeR_TopTags_QLF.RData') #We will need this later

#--------------
# Filter based on significance criteria (e.g., adjusted p-value and fold-change)

# eBayes
mask <- WC_result_eBayes$adj.P.Val < 0.05 & abs(WC_result_eBayes$logFC) > 1
deGenes <- rownames(WC_result_eBayes[mask, ])
head(deGenes)

# [1] "Do1_03_a00001G00267V1.1" "Do1_01_a00007G00164V1.1"
# [3] "Do1_04_a00002G00429V1.1" "Do1_02_a00003G01208V1.1"
# [5] "Do1_05_a00001G01169V1.1" "Do1_01_a00007G00166V1.1"

#------------
# LRT
DE_genes <- as.data.frame(WC_result_LRT)  # Extract the result table
DE_genes <- DE_genes[DE_genes$FDR < 0.05 & abs(DE_genes$logFC) > 0, ]

head(DE_genes)
                                         # genes     logFC   logCPM       LR
# Do1_05_a00001G00678V1.1 Do1_05_a00001G00678V1.1 0.4125184 5.014348 26.06703
# Do1_05_a00001G00394V1.1 Do1_05_a00001G00394V1.1 0.5555086 6.033763 23.56654
                              # PValue         FDR
# Do1_05_a00001G00678V1.1 3.297668e-07 0.007657844
# Do1_05_a00001G00394V1.1 1.206677e-06 0.014010730

#--------------
summary(decideTests(lrt, method="hierarchical", adjust.method="fdr", p.value=0.05,lfc=0))

       # -1*groupC 1*groupW
# Down                    0
# NotSig              23220
# Up                      2

summary(decideTests(fit2, method="hierarchical", adjust.method="fdr", p.value=0.05, lfc=0))

       # groupWvsgroupC
# Down               21
# NotSig          23186
# Up                 15

summary(decideTests(fit2, method="hierarchical", adjust.method="fdr", p.value=0.05, lfc=1))
       # groupWvsgroupC
# Down                6
# NotSig          23212
# Up                  4

#------------
# get the DERs and their information 
lrt$table$fdr<-p.adjust(lrt$table$PValue, method='fdr')

dim(lrt$table[lrt$table$fdr<0.05,])
# 2 5
length(lrt$table$PValue[lrt$table$PValue<0.05])
# 1985

DERs_FDR <- lrt$table[lrt$table$fdr<0.05,]

dim(DERs_FDR)
# treat 29

write.table(DERs_FDR, file = "./RNA_DERs_Nov2024_W_C_Total.txt", quote = FALSE, row.names=TRUE, col.names=TRUE, sep="\t")

                           # logFC   logCPM       LR       PValue         fdr
# Do1_05_a00001G00394V1.1 0.5555086 6.033763 23.56654 1.206677e-06 0.014010730
# Do1_05_a00001G00678V1.1 0.4125184 5.014348 26.06703 3.297668e-07 0.007657844

#---------------------
# total gene universe for GO enrichment

geneUniverse <- rownames(lrt$table)
length(geneUniverse)
# 23222

#---------------------
# Finally, we can plot the log-fold changes of all the genes, and the highlight those that are differentially expressed
#?decideTests
deGenes <- decideTestsDGE(lrt, p=0.05)
deGenes <- rownames(lrt)[as.logical(deGenes)]

png("./plots/Plot_smear.png", width = 700, height = 500)
plotSmear(lrt, de.tags=deGenes)
abline(h=c(-1, 1), col=2)
dev.off()

#####################
# Try glmm to include random effects (site)
# https://myles-lewis.github.io/glmmSeq/articles/glmmSeq.html#metadata
# https://myles-lewis.github.io/glmmSeq/

install.packages("glmmSeq")

library(glmmSeq)
set.seed(1234)

tpm <- dge.filt

disp  <- setNames(edgeR::estimateDisp(tpm)$tagwise.dispersion, rownames(tpm))

sizeFactors <- calcNormFactors(counts, method="TMM")

results <- glmmSeq(~ treatment + (1 | site),
                   countdata = tpm,
                   metadata = metadata,
                   dispersion = disp,
                   progress = TRUE)

names(attributes(results))

stats <- summary(results)

results1 <- glmmQvals(results)










##################################
# Run for each site

tmux new-session -s RNA
tmux attach-session -t RNA

# Cedar1
cd /home/celphin/projects/rrg-rieseber-ac/rpp-rieseber/celphin/Dryas/RNAseq_analysis/

# copy to Beluga scratch/
cd /home/celphin/scratch/Dryas/RNAseq_analysis/


# get R on server
# open R
module load StdEnv/2020
module load r/4.2.1
module load gdal
module load udunits
module load python
#export R_LIBS_USER=/home/msandler/R/x86_64-pc-linux-gnu-library/4.2.1/
export R_LIBS_USER=/home/celphin/R/x86_64-pc-linux-gnu-library/4.2.1/

R

setwd ("/home/celphin/projects/rrg-rieseber-ac/rpp-rieseber/celphin/Dryas/RNAseq_analysis/")
setwd ("/home/celphin/scratch/Dryas/RNAseq_analysis/")

#find your working directory:
getwd()

#load the libraries you will need 
library(edgeR)
library(gplots)
library(stringr)
library(statmod)
#-----------------------
#read in the data
mydata <- read.table("./RNA_site_tables/gene_names_expression_table1.txt", header=TRUE)

CDS <- mydata[,1]
exp_data <- mydata[,c(2:length(mydata))]

#make a DGElist object out of the data
dge <- DGEList (counts = exp_data, genes=rownames(exp_data))

# normalize
dge_norm <- calcNormFactors(dge, method="TMM")
norm_cpm <- cpm(dge_norm, normalized.lib.sizes = TRUE)

#Design matrix
treat <- as.factor (sapply (strsplit(colnames(exp_data),split = "_"),"[[",1))
site <- as.factor (sapply (strsplit(colnames(exp_data),split = "_"),"[[",2))
treat_site <- paste0(as.character(site), "_", as.character(treat))
treat_site <- as.factor(treat_site)

designMat <- model.matrix( ~ 0 + treat_site, data=dge_norm$samples)

# Filter
keep.exprs <- filterByExpr(dge_norm, design = designMat)
dge.filt <- dge_norm[keep.exprs,]
dim(dge.filt)
# [1] 24834    69

# Dispersion
dge.filt <- estimateGLMCommonDisp(dge.filt, design=designMat)
dge.filt <- estimateGLMTrendedDisp(dge.filt, design=designMat)
dge.filt <- estimateGLMTagwiseDisp(dge.filt, design=designMat)

# Fit models
fit <- glmFit(dge.filt, designMat,  robust=TRUE)
fit_QLF <- glmQLFit(dge.filt, designMat, robust=TRUE)

#-----
# Alaska

con <- makeContrasts(treat_siteAlas_W - treat_siteAlas_C, levels=designMat)
qlf <- glmQLFTest(fit_QLF, contrast=con)
lrt <- glmLRT(fit, contrast=con)
topTags(lrt)
summary(decideTests(lrt, method="hierarchical", adjust.method="fdr", p.value=0.05,lfc=0))

       # -1*treat_siteAlas_C 1*treat_siteAlas_W
# Down                                       37
# NotSig                                  24762
# Up                                         35


summary(decideTests(qlf, method="hierarchical", adjust.method="fdr", p.value=0.05,lfc=0))

       # -1*treat_siteAlas_C 1*treat_siteAlas_W
# Down                                       34
# NotSig                                  24775
# Up                                         25

Alaska_W_C_DERs <- lrt$table[which(!decideTests(lrt, method="hierarchical", adjust.method="fdr", p.value=0.05,lfc=1)==0),]
write.table(Alaska_W_C_DERs, file = "./RNA_Alaska_W_C_DERs.txt", quote = FALSE, row.names=TRUE, col.names=TRUE, sep="\t")

#----------------
Alaska_W_C_DERs_updown <- as.data.frame(decideTests(lrt, method="hierarchical", adjust.method="fdr", p.value=0.05,lfc=0))
Alaska_W_C_DERs_updown$Gene <- rownames(Alaska_W_C_DERs_updown)
colnames(Alaska_W_C_DERs_updown) <- c("UpDown", "Gene")
Alaska_W_C_DERs_updown0 <- Alaska_W_C_DERs_updown[-which(Alaska_W_C_DERs_updown$UpDown==0),]
write.table(Alaska_W_C_DERs_updown0, file = "./RNA_Alaska_W_C_DERs_updown.txt", quote = FALSE, row.names=FALSE, col.names=TRUE, sep="\t")

#-------------
png("./plots/Plot_MD_Alaska.png", width = 700, height = 500)
plotMD(lrt)
dev.off()

pdf("./plots/Plot_MD_Alaska.pdf")
plotMD(lrt)
dev.off()


#-----
# Norway

con <- makeContrasts(treat_siteNorw_W - treat_siteNorw_C, levels=designMat)
qlf <- glmQLFTest(fit_QLF, contrast=con)
lrt <- glmLRT(fit, contrast=con)
topTags(lrt)
summary(decideTests(lrt, method="hierarchical", adjust.method="fdr", p.value=0.05,lfc=0))

       # -1*treat_siteNorw_C 1*treat_siteNorw_W
# Down                                        9
# NotSig                                  24820
# Up                                          5


summary(decideTests(qlf, method="hierarchical", adjust.method="fdr", p.value=0.05,lfc=0))
       # -1*treat_siteNorw_C 1*treat_siteNorw_W
# Down                                        3
# NotSig                                  24829
# Up                                          2

Norway_W_C_DERs <-  lrt$table[which(!decideTests(lrt, method="hierarchical", adjust.method="fdr", p.value=0.05,lfc=1)==0),]
write.table(Norway_W_C_DERs, file = "./RNA_Norway_W_C_DERs.txt", quote = FALSE, row.names=TRUE, col.names=TRUE, sep="\t")

#----------------
Norway_W_C_DERs_updown <- as.data.frame(decideTests(lrt, method="hierarchical", adjust.method="fdr", p.value=0.05,lfc=0))
Norway_W_C_DERs_updown$Gene <- rownames(Norway_W_C_DERs_updown)
colnames(Norway_W_C_DERs_updown) <- c("UpDown", "Gene")
Norway_W_C_DERs_updown0 <- Norway_W_C_DERs_updown[-which(Norway_W_C_DERs_updown$UpDown==0),]
write.table(Norway_W_C_DERs_updown0, file = "./RNA_Norway_W_C_DERs_updown.txt", quote = FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
#----------------------

png("./plots/Plot_MD_Norway.png", width = 700, height = 500)
plotMD(lrt)
dev.off()

pdf("./plots/Plot_MD_Norway.pdf")
plotMD(lrt)
dev.off()

#-----
# Sweden

con <- makeContrasts(treat_siteSwed_W - treat_siteSwed_C, levels=designMat)
qlf <- glmQLFTest(fit_QLF, contrast=con)
lrt <- glmLRT(fit, contrast=con)
topTags(lrt)
summary(decideTests(lrt, method="hierarchical", adjust.method="fdr", p.value=0.05,lfc=0))

       # -1*treat_siteSwed_C 1*treat_siteSwed_W
# Down                                       25
# NotSig                                  24732
# Up                                         77

summary(decideTests(qlf, method="hierarchical", adjust.method="fdr", p.value=0.05,lfc=0))

       # -1*treat_siteSwed_C 1*treat_siteSwed_W
# Down                                        8
# NotSig                                  24806
# Up                                         20

Sweden_W_C_DERs <-   lrt$table[which(!decideTests(lrt, method="hierarchical", adjust.method="fdr", p.value=0.05,lfc=1)==0),]
write.table(Sweden_W_C_DERs, file = "./RNA_Sweden_W_C_DERs.txt", quote = FALSE, row.names=TRUE, col.names=TRUE, sep="\t")

#----------------
Sweden_W_C_DERs_updown <- as.data.frame(decideTests(lrt, method="hierarchical", adjust.method="fdr", p.value=0.05,lfc=0))
Sweden_W_C_DERs_updown$Gene <- rownames(Sweden_W_C_DERs_updown)
colnames(Sweden_W_C_DERs_updown) <- c("UpDown", "Gene")
Sweden_W_C_DERs_updown0 <- Sweden_W_C_DERs_updown[-which(Sweden_W_C_DERs_updown$UpDown==0),]
write.table(Sweden_W_C_DERs_updown0, file = "./RNA_Sweden_W_C_DERs_updown.txt", quote = FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
#----------------------

png("./plots/Plot_MD_Swed.png", width = 700, height = 500)
plotMD(lrt)
dev.off()

pdf("./plots/Plot_MD_Swed.pdf")
plotMD(lrt)
dev.off()

#-----
# Alex fiord
con <- makeContrasts(treat_siteAlex_W  - treat_siteAlex_C, levels=designMat)
qlf <- glmQLFTest(fit_QLF, contrast=con)
lrt <- glmLRT(fit, contrast=con)
topTags(lrt)
summary(decideTests(lrt, method="hierarchical", adjust.method="fdr", p.value=0.05,lfc=0))

       # -1*treat_siteAlex_C 1*treat_siteAlex_W
# Down                                       33
# NotSig                                  24790
# Up                                         11

summary(decideTests(qlf, method="hierarchical", adjust.method="fdr", p.value=0.05,lfc=0))
       # -1*treat_siteAlex_C 1*treat_siteAlex_W
# Down                                       14
# NotSig                                  24817
# Up                                          3

Alex_W_C_DERs <-  lrt$table[which(!decideTests(lrt, method="hierarchical", adjust.method="fdr", p.value=0.05,lfc=1)==0),]
write.table(Alex_W_C_DERs, file = "./RNA_Alex_W_C_DERs.txt", quote = FALSE, row.names=TRUE, col.names=TRUE, sep="\t")

#----------------
Alex_W_C_DERs_updown <- as.data.frame(decideTests(lrt, method="hierarchical", adjust.method="fdr", p.value=0.05,lfc=0))
Alex_W_C_DERs_updown$Gene <- rownames(Alex_W_C_DERs_updown)
colnames(Alex_W_C_DERs_updown) <- c("UpDown", "Gene")
Alex_W_C_DERs_updown0 <- Alex_W_C_DERs_updown[-which(Alex_W_C_DERs_updown$UpDown==0),]
write.table(Alex_W_C_DERs_updown0, file = "./RNA_Alex_W_C_DERs_updown.txt", quote = FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
#----------------------

png("./plots/Plot_MD_Alex.png", width = 700, height = 500)
plotMD(lrt)
dev.off()

pdf("./plots/Plot_MD_Alex.pdf")
plotMD(lrt)
dev.off()


#-----
# Seedlings

con <- makeContrasts(treat_siteSeed_W - treat_siteSeed_C, levels=designMat)
qlf <- glmQLFTest(fit_QLF, contrast=con)
lrt <- glmLRT(fit, contrast=con)
topTags(lrt)
summary(decideTests(lrt, method="hierarchical", adjust.method="fdr", p.value=0.05,lfc=0))

       # -1*treat_siteSeed_C 1*treat_siteSeed_W
# Down                                       19
# NotSig                                  24767
# Up                                         48


summary(decideTests(qlf, method="hierarchical", adjust.method="fdr", p.value=0.05,lfc=0))
       # -1*treat_siteSeed_C 1*treat_siteSeed_W
# Down                                        5
# NotSig                                  24823
# Up                                          6


Seedling_W_C_DERs <-  lrt$table[which(!decideTests(lrt, method="hierarchical", adjust.method="fdr", p.value=0.05,lfc=1)==0),]
write.table(Seedling_W_C_DERs, file = "./RNA_Seedling_W_C_DERs.txt", quote = FALSE, row.names=TRUE, col.names=TRUE, sep="\t")

#----------------
Seedling_W_C_DERs_updown <- as.data.frame(decideTests(lrt, method="hierarchical", adjust.method="fdr", p.value=0.05,lfc=0))
Seedling_W_C_DERs_updown$Gene <- rownames(Seedling_W_C_DERs_updown)
colnames(Seedling_W_C_DERs_updown) <- c("UpDown", "Gene")
Seedling_W_C_DERs_updown0 <- Seedling_W_C_DERs_updown[-which(Seedling_W_C_DERs_updown$UpDown==0),]
write.table(Seedling_W_C_DERs_updown0, file = "./RNA_Seedling_W_C_DERs_updown.txt", quote = FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
#----------------------

png("./plots/Plot_MD_Seedlings.png", width = 700, height = 500)
plotMD(lrt)
dev.off()

pdf("./plots/Plot_MD_Seedlings.pdf")
plotMD(lrt)
dev.off()


#####################################
# Heatmaps

#make total heatmap by getting the counts per million from each gene and turning them relative proportions (columns add up to 1)
logcpm <- cpm(dge.filt, log=TRUE)

sub1 <- colSums(logcpm)
sub2 <- matrix(rep(sub1,nrow(logcpm)), c(nrow (logcpm),ncol(logcpm)),byrow = TRUE)
sub3 <- logcpm / sub2

#subset this matrix to just get the genes from above
names1 <- row.names (sub3)
length(names1)
#24834

names2 <- c(row.names(Seedling_W_C_DERs), 
            row.names(Sweden_W_C_DERs), 
            row.names(Norway_W_C_DERs), 
            row.names(Alaska_W_C_DERs),
            row.names(Alex_W_C_DERs))

duplicated(names2) 

# remove duplicates
names3 <- unique(names2)

index2 <- names1 %in% names3
heatmap1 <- sub3[index2,]

#------------------------------
# All sites
##     par(mar = c(bottom, left, top, right))
# plot used for GenomeBC report
#play around with the options to make the plot fit what you like for options type ?heatmap.2

# order by site
# http://www.sthda.com/english/wiki/reordering-data-frame-columns-in-r
colnames(heatmap1)
 # [1] "C_Alas_1"    "C_Alas_12"   "C_Alas_14"   "C_Alas_15"   "C_Alas_16"
 # [6] "C_Alas_20"   "C_Alas_5"    "C_Alas_7"    "C_Alex_11"   "C_Alex_2"
# [11] "C_Alex_21"   "C_Alex_22"   "C_Alex_23"   "C_Alex_24"   "C_Alex_6"
# [16] "C_Alex_7"    "C_Norw_1"    "C_Norw_10"   "C_Norw_2_B"  "C_Norw_3_B"
# [21] "C_Norw_5"    "C_Norw_6"    "C_Norw_7"    "C_Norw_9_B"  "C_Seed_1"
# [26] "C_Seed_2"    "C_Seed_3"    "C_Swed_11"   "C_Swed_15"   "C_Swed_16"
# [31] "C_Swed_2"    "C_Swed_3"    "C_Swed_5"    "C_Swed_8"    "W_Alas_1"
# [36] "W_Alas_12"   "W_Alas_13"   "W_Alas_15"   "W_Alas_2"    "W_Alas_3"
# [41] "W_Alas_5"    "W_Alas_9"    "W_Alex_1"    "W_Alex_11"   "W_Alex_21"
# [46] "W_Alex_22"   "W_Alex_23"   "W_Alex_24"   "W_Alex_6"    "W_Alex_7"
# [51] "W_Norw_10_B" "W_Norw_11"   "W_Norw_2"    "W_Norw_3_B"  "W_Norw_4_B"
# [56] "W_Norw_6"    "W_Norw_7"    "W_Norw_8"    "W_Seed_1"    "W_Seed_2"
# [61] "W_Seed_3"    "W_Swed_10"   "W_Swed_13"   "W_Swed_14"   "W_Swed_16"
# [66] "W_Swed_17"   "W_Swed_2"    "W_Swed_5"    "W_Swed_6"

heatmap2 <- heatmap1 #[, col_order]

png("./plots/W_C_RNA_heatmap.png", width = 700, height = 500)
heatmap.2(heatmap2, scale = "row", trace = "none", Colv = FALSE, 
          labCol = colnames(heatmap2), labRow = rownames(heatmap2), 
          col = 'rainbow', dendrogram = "none", colsep = c(7), 
          margins = c(10, 15), sepcolor = "white", sepwidth = 0.1, 
          cexRow = 1, cexCol = 1)
dev.off()

pdf("./plots/W_C_RNA_heatmap.pdf", width = 12, height = 8)
heatmap.2(heatmap2, scale = "row", trace = "none", Colv = FALSE, 
          labCol = colnames(heatmap2), labRow = rownames(heatmap2), 
          col = 'rainbow', dendrogram = "none", colsep = c(7), 
          margins = c(10, 15), sepcolor = "white", sepwidth = 0.1, 
          cexRow = 1, cexCol = 1)
dev.off()

# write out file of genes that differ a lot
write.table(heatmap2, file = "./RNA_DER_May2023_W_C_Total.txt", quote = FALSE, row.names=TRUE, col.names=TRUE, sep="\t")


##################################################
# Alex

names2 <- c(row.names(Alex_W_C_DERs))

duplicated(names2) 

# remove duplicates
names3 <- unique(names2)

index2 <- names1 %in% names3
heatmap1 <- sub3[index2,]

colnames(heatmap1)
 # [1] "C_Alas_1"    "C_Alas_12"   "C_Alas_14"   "C_Alas_15"   "C_Alas_16"
 # [6] "C_Alas_20"   "C_Alas_5"    "C_Alas_7"    "C_Alex_11"   "C_Alex_2"
# [11] "C_Alex_21"   "C_Alex_22"   "C_Alex_23"   "C_Alex_24"   "C_Alex_6"
# [16] "C_Alex_7"    "C_Norw_1"    "C_Norw_10"   "C_Norw_2_B"  "C_Norw_3_B"
# [21] "C_Norw_5"    "C_Norw_6"    "C_Norw_7"    "C_Norw_9_B"  "C_Seed_1"
# [26] "C_Seed_2"    "C_Seed_3"    "C_Swed_11"   "C_Swed_15"   "C_Swed_16"
# [31] "C_Swed_2"    "C_Swed_3"    "C_Swed_5"    "C_Swed_8"    "W_Alas_1"
# [36] "W_Alas_12"   "W_Alas_13"   "W_Alas_15"   "W_Alas_2"    "W_Alas_3"
# [41] "W_Alas_5"    "W_Alas_9"    "W_Alex_1"    "W_Alex_11"   "W_Alex_21"
# [46] "W_Alex_22"   "W_Alex_23"   "W_Alex_24"   "W_Alex_6"    "W_Alex_7"
# [51] "W_Norw_10_B" "W_Norw_11"   "W_Norw_2"    "W_Norw_3_B"  "W_Norw_4_B"
# [56] "W_Norw_6"    "W_Norw_7"    "W_Norw_8"    "W_Seed_1"    "W_Seed_2"
# [61] "W_Seed_3"    "W_Swed_10"   "W_Swed_13"   "W_Swed_14"   "W_Swed_16"
# [66] "W_Swed_17"   "W_Swed_2"    "W_Swed_5"    "W_Swed_6"

heatmap2 <- heatmap1[,c(9:16,43:50)]

heatmap2_ordered_by_row_sum <- heatmap2[order(rowSums(heatmap2[,c(1:8)]), decreasing = TRUE), ]

png("./plots/Alex_W_C_RNA_heatmap2.png", width = 700, height = 500)
heatmap.2(heatmap2, scale = "row", trace = "none", Colv = FALSE, 
          labCol = colnames(heatmap2), labRow = rownames(heatmap2), 
          col = 'rainbow', dendrogram = "none", colsep = c(8), 
          margins = c(10, 15), sepcolor = "white", sepwidth = 0.1, 
          cexRow = 1, cexCol = 1)
dev.off()

pdf("./plots/Alex_W_C_RNA_heatmap2.pdf", width = 12, height = 8)
heatmap.2(heatmap2, scale = "row", trace = "none", Colv = FALSE, 
          labCol = colnames(heatmap2), labRow = rownames(heatmap2), 
          col = 'rainbow', dendrogram = "none", colsep = c(8), 
          margins = c(10, 15), sepcolor = "white", sepwidth = 0.1, 
          cexRow = 1, cexCol = 1)
dev.off()

# write out file of genes that differ a lot
write.table(heatmap2, file = "./RNA_DER_W_C_Alex.txt", quote = FALSE, row.names=TRUE, col.names=TRUE, sep="\t")


#------------------------------
# Norway

names2 <- c(row.names(Norway_W_C_DERs))

duplicated(names2) 

# remove duplicates
names3 <- unique(names2)

index2 <- names1 %in% names3
heatmap1 <- sub3[index2,]

colnames(heatmap1)
 # [1] "C_Alas_1"    "C_Alas_12"   "C_Alas_14"   "C_Alas_15"   "C_Alas_16"
 # [6] "C_Alas_20"   "C_Alas_5"    "C_Alas_7"    "C_Alex_11"   "C_Alex_2"
# [11] "C_Alex_21"   "C_Alex_22"   "C_Alex_23"   "C_Alex_24"   "C_Alex_6"
# [16] "C_Alex_7"    "C_Norw_1"    "C_Norw_10"   "C_Norw_2_B"  "C_Norw_3_B"
# [21] "C_Norw_5"    "C_Norw_6"    "C_Norw_7"    "C_Norw_9_B"  "C_Seed_1"
# [26] "C_Seed_2"    "C_Seed_3"    "C_Swed_11"   "C_Swed_15"   "C_Swed_16"
# [31] "C_Swed_2"    "C_Swed_3"    "C_Swed_5"    "C_Swed_8"    "W_Alas_1"
# [36] "W_Alas_12"   "W_Alas_13"   "W_Alas_15"   "W_Alas_2"    "W_Alas_3"
# [41] "W_Alas_5"    "W_Alas_9"    "W_Alex_1"    "W_Alex_11"   "W_Alex_21"
# [46] "W_Alex_22"   "W_Alex_23"   "W_Alex_24"   "W_Alex_6"    "W_Alex_7"
# [51] "W_Norw_10_B" "W_Norw_11"   "W_Norw_2"    "W_Norw_3_B"  "W_Norw_4_B"
# [56] "W_Norw_6"    "W_Norw_7"    "W_Norw_8"    "W_Seed_1"    "W_Seed_2"
# [61] "W_Seed_3"    "W_Swed_10"   "W_Swed_13"   "W_Swed_14"   "W_Swed_16"
# [66] "W_Swed_17"   "W_Swed_2"    "W_Swed_5"    "W_Swed_6"

heatmap2 <- heatmap1[,c(17:24, 51:58)]
colnames(heatmap2)

heatmap2_ordered_by_row_sum <- heatmap2[order(rowSums(heatmap2[,c(1:8)]), decreasing = TRUE), ]

png("./plots/Norway_W_C_RNA_heatmap2.png", width = 700, height = 500)
heatmap.2(heatmap2, scale = "row", trace = "none", Colv = FALSE, 
          labCol = colnames(heatmap2), labRow = rownames(heatmap2), 
          col = 'rainbow', dendrogram = "none", colsep = c(8), 
          margins = c(10, 15), sepcolor = "white", sepwidth = 0.1, 
          cexRow = 1, cexCol = 1)
dev.off()

pdf("./plots/Norway_W_C_RNA_heatmap2.pdf", width = 12, height = 8)
heatmap.2(heatmap2, scale = "row", trace = "none", Colv = FALSE, 
          labCol = colnames(heatmap2), labRow = rownames(heatmap2), 
          col = 'rainbow', dendrogram = "none", colsep = c(8), 
          margins = c(10, 15), sepcolor = "white", sepwidth = 0.1, 
          cexRow = 1, cexCol = 1)
dev.off()

# write out file of genes that differ a lot
write.table(heatmap2, file = "./RNA_DER_W_C_Norway.txt", quote = FALSE, row.names=TRUE, col.names=TRUE, sep="\t")

#------------------------------
# Alaska

names2 <- c(row.names(Alaska_W_C_DERs))

duplicated(names2) 

# remove duplicates
names3 <- unique(names2)

index2 <- names1 %in% names3
heatmap1 <- sub3[index2,]

colnames(heatmap1)
 # [1] "C_Alas_1"    "C_Alas_12"   "C_Alas_14"   "C_Alas_15"   "C_Alas_16"
 # [6] "C_Alas_20"   "C_Alas_5"    "C_Alas_7"    "C_Alex_11"   "C_Alex_2"
# [11] "C_Alex_21"   "C_Alex_22"   "C_Alex_23"   "C_Alex_24"   "C_Alex_6"
# [16] "C_Alex_7"    "C_Norw_1"    "C_Norw_10"   "C_Norw_2_B"  "C_Norw_3_B"
# [21] "C_Norw_5"    "C_Norw_6"    "C_Norw_7"    "C_Norw_9_B"  "C_Seed_1"
# [26] "C_Seed_2"    "C_Seed_3"    "C_Swed_11"   "C_Swed_15"   "C_Swed_16"
# [31] "C_Swed_2"    "C_Swed_3"    "C_Swed_5"    "C_Swed_8"    "W_Alas_1"
# [36] "W_Alas_12"   "W_Alas_13"   "W_Alas_15"   "W_Alas_2"    "W_Alas_3"
# [41] "W_Alas_5"    "W_Alas_9"    "W_Alex_1"    "W_Alex_11"   "W_Alex_21"
# [46] "W_Alex_22"   "W_Alex_23"   "W_Alex_24"   "W_Alex_6"    "W_Alex_7"
# [51] "W_Norw_10_B" "W_Norw_11"   "W_Norw_2"    "W_Norw_3_B"  "W_Norw_4_B"
# [56] "W_Norw_6"    "W_Norw_7"    "W_Norw_8"    "W_Seed_1"    "W_Seed_2"
# [61] "W_Seed_3"    "W_Swed_10"   "W_Swed_13"   "W_Swed_14"   "W_Swed_16"
# [66] "W_Swed_17"   "W_Swed_2"    "W_Swed_5"    "W_Swed_6"

heatmap2 <- heatmap1[,c(1:8, 35:42)]
colnames(heatmap2)

heatmap2_ordered_by_row_sum <- heatmap2[order(rowSums(heatmap2[,c(1:8)]), decreasing = TRUE), ]

png("./plots/Alaska_W_C_RNA_heatmap2.png", width = 700, height = 500)
heatmap.2(heatmap2, scale = "row", trace = "none", Colv = FALSE, 
          labCol = colnames(heatmap2), labRow = rownames(heatmap2), 
          col = 'rainbow', dendrogram = "none", colsep = c(8), 
          margins = c(10, 15), sepcolor = "white", sepwidth = 0.1, 
          cexRow = 1, cexCol = 1)
dev.off()

pdf("./plots/Alaska_W_C_RNA_heatmap2.pdf", width = 12, height = 8)
heatmap.2(heatmap2, scale = "row", trace = "none", Colv = FALSE, 
          labCol = colnames(heatmap2), labRow = rownames(heatmap2), 
          col = 'rainbow', dendrogram = "none", colsep = c(8), 
          margins = c(10, 15), sepcolor = "white", sepwidth = 0.1, 
          cexRow = 1, cexCol = 1)
dev.off()

# write out file of genes that differ a lot
write.table(heatmap2, file = "./RNA_DER_W_C_Alaska.txt", quote = FALSE, row.names=TRUE, col.names=TRUE, sep="\t")

#----------------------------------
# Sweden

names2 <- c(row.names(Sweden_W_C_DERs))

duplicated(names2) 

# remove duplicates
names3 <- unique(names2)

index2 <- names1 %in% names3
heatmap1 <- sub3[index2,]

colnames(heatmap1)
 # [1] "C_Alas_1"    "C_Alas_12"   "C_Alas_14"   "C_Alas_15"   "C_Alas_16"
 # [6] "C_Alas_20"   "C_Alas_5"    "C_Alas_7"    "C_Alex_11"   "C_Alex_2"
# [11] "C_Alex_21"   "C_Alex_22"   "C_Alex_23"   "C_Alex_24"   "C_Alex_6"
# [16] "C_Alex_7"    "C_Norw_1"    "C_Norw_10"   "C_Norw_2_B"  "C_Norw_3_B"
# [21] "C_Norw_5"    "C_Norw_6"    "C_Norw_7"    "C_Norw_9_B"  "C_Seed_1"
# [26] "C_Seed_2"    "C_Seed_3"    "C_Swed_11"   "C_Swed_15"   "C_Swed_16"
# [31] "C_Swed_2"    "C_Swed_3"    "C_Swed_5"    "C_Swed_8"    "W_Alas_1"
# [36] "W_Alas_12"   "W_Alas_13"   "W_Alas_15"   "W_Alas_2"    "W_Alas_3"
# [41] "W_Alas_5"    "W_Alas_9"    "W_Alex_1"    "W_Alex_11"   "W_Alex_21"
# [46] "W_Alex_22"   "W_Alex_23"   "W_Alex_24"   "W_Alex_6"    "W_Alex_7"
# [51] "W_Norw_10_B" "W_Norw_11"   "W_Norw_2"    "W_Norw_3_B"  "W_Norw_4_B"
# [56] "W_Norw_6"    "W_Norw_7"    "W_Norw_8"    "W_Seed_1"    "W_Seed_2"
# [61] "W_Seed_3"    "W_Swed_10"   "W_Swed_13"   "W_Swed_14"   "W_Swed_16"
# [66] "W_Swed_17"   "W_Swed_2"    "W_Swed_5"    "W_Swed_6"

heatmap2 <- heatmap1[,c(28:34, 62:69)]
colnames(heatmap2)

heatmap2_ordered_by_row_sum <- heatmap2[order(rowSums(heatmap2[,c(1:7)]), decreasing = TRUE), ]

png("./plots/Sweden_W_C_RNA_heatmap2.png", width = 700, height = 500)
heatmap.2(heatmap2, scale = "row", trace = "none", Colv = FALSE, 
          labCol = colnames(heatmap2), labRow = rownames(heatmap2), 
          col = 'rainbow', dendrogram = "none", colsep = c(7), 
          margins = c(10, 15), sepcolor = "white", sepwidth = 0.1, 
          cexRow = 1, cexCol = 1)
dev.off()

pdf("./plots/Sweden_W_C_RNA_heatmap2.pdf", width = 12, height = 8)
heatmap.2(heatmap2, scale = "row", trace = "none", Colv = FALSE, 
          labCol = colnames(heatmap2), labRow = rownames(heatmap2), 
          col = 'rainbow', dendrogram = "none", colsep = c(7), 
          margins = c(10, 15), sepcolor = "white", sepwidth = 0.1, 
          cexRow = 1, cexCol = 1)
dev.off()

# write out file of genes that differ a lot
write.table(heatmap2, file = "./RNA_DER_W_C_Sweden.txt", quote = FALSE, row.names=TRUE, col.names=TRUE, sep="\t")

#----------------------------------
# Seedlings

names2 <- c(row.names(Seedling_W_C_DERs))

duplicated(names2) 

# remove duplicates
names3 <- unique(names2)

index2 <- names1 %in% names3
heatmap1 <- sub3[index2,]

colnames(heatmap1)
 # [1] "C_Alas_1"    "C_Alas_12"   "C_Alas_14"   "C_Alas_15"   "C_Alas_16"
 # [6] "C_Alas_20"   "C_Alas_5"    "C_Alas_7"    "C_Alex_11"   "C_Alex_2"
# [11] "C_Alex_21"   "C_Alex_22"   "C_Alex_23"   "C_Alex_24"   "C_Alex_6"
# [16] "C_Alex_7"    "C_Norw_1"    "C_Norw_10"   "C_Norw_2_B"  "C_Norw_3_B"
# [21] "C_Norw_5"    "C_Norw_6"    "C_Norw_7"    "C_Norw_9_B"  "C_Seed_1"
# [26] "C_Seed_2"    "C_Seed_3"    "C_Swed_11"   "C_Swed_15"   "C_Swed_16"
# [31] "C_Swed_2"    "C_Swed_3"    "C_Swed_5"    "C_Swed_8"    "W_Alas_1"
# [36] "W_Alas_12"   "W_Alas_13"   "W_Alas_15"   "W_Alas_2"    "W_Alas_3"
# [41] "W_Alas_5"    "W_Alas_9"    "W_Alex_1"    "W_Alex_11"   "W_Alex_21"
# [46] "W_Alex_22"   "W_Alex_23"   "W_Alex_24"   "W_Alex_6"    "W_Alex_7"
# [51] "W_Norw_10_B" "W_Norw_11"   "W_Norw_2"    "W_Norw_3_B"  "W_Norw_4_B"
# [56] "W_Norw_6"    "W_Norw_7"    "W_Norw_8"    "W_Seed_1"    "W_Seed_2"
# [61] "W_Seed_3"    "W_Swed_10"   "W_Swed_13"   "W_Swed_14"   "W_Swed_16"
# [66] "W_Swed_17"   "W_Swed_2"    "W_Swed_5"    "W_Swed_6"

heatmap2 <- heatmap1[,c(25:27, 59:61)]
colnames(heatmap2)

heatmap2_ordered_by_row_sum <- heatmap2[order(rowSums(heatmap2[,c(1:3)]), decreasing = TRUE), ]

png("./plots/Seedling_W_C_RNA_heatmap2.png", width = 700, height = 500)
heatmap.2(heatmap2, scale = "row", trace = "none", Colv = FALSE, 
          labCol = colnames(heatmap2), labRow = rownames(heatmap2), 
          col = 'rainbow', dendrogram = "none", colsep = c(3), 
          margins = c(10, 15), sepcolor = "white", sepwidth = 0.1, 
          cexRow = 1, cexCol = 1)
dev.off()

pdf("./plots/Seedling_W_C_RNA_heatmap2.pdf", width = 12, height = 8)
heatmap.2(heatmap2, scale = "row", trace = "none", Colv = FALSE, 
          labCol = colnames(heatmap2), labRow = rownames(heatmap2), 
          col = 'rainbow', dendrogram = "none", colsep = c(3), 
          margins = c(10, 15), sepcolor = "white", sepwidth = 0.1, 
          cexRow = 1, cexCol = 1)
dev.off()

# write out file of genes that differ a lot
write.table(heatmap2, file = "./RNA_DER_W_C_Seedling.txt", quote = FALSE, row.names=TRUE, col.names=TRUE, sep="\t")

###################################
# Sweden differential expression - specific genes plotting 
library dplyr)
library(ggplot2)

head(logcpm)

#plot the individual expression values from a single gene:

logcpm_names <- colnames(logcpm)
values <- as.data.frame( t(rbind( logcpm_names, logcpm)))

data <- values %>%
  mutate(
    site = str_extract(logcpm_names, "_[A-Za-z]+"),  # Extract site (e.g., "Norw", "Seed", "Swed", etc.)  
    treatment_type = ifelse(grepl("^W_", logcpm_names), "W", "C")  # Identify W vs C (W is for "W_" and C for "C_")
  )

#------------------------------
# Do1_05_a00001G00678V1.1 Do1_05_a00001G00678V1.1  0.4082366  4.9961193  5.189187
# Do1_04_a00002G00429V1.1 Do1_04_a00002G00429V1.1 -1.6381272 -1.0473707 -5.457280
# Do1_01_a00007G00164V1.1 Do1_01_a00007G00164V1.1  2.5823275 -2.5398183  5.321775
# Do1_03_a00002G01336V1.1 Do1_03_a00002G01336V1.1 -0.9156573  0.7377493 -4.973376
# Do1_03_a00001G00267V1.1 Do1_03_a00001G00267V1.1 -1.6954531  2.8551757 -4.911546
# Do1_02_a00003G01059V1.1 Do1_02_a00003G01059V1.1 -0.3332865  5.7836771 -4.834632

#-------------------------
# Do1_05_a00001G00678V1.1
# Do1_05_a00001G00678V1.1,TCP family transcription factor,1,1,1,1,1,AT3G47620,LotjaGi4g1v0248700,Ro05_G03227,Medtr4g108370,FvH4_5g15150
# https://www.arabidopsis.org/locus?name=AT3G47620

subdat1<- data[,c("site","treatment_type", "Do1_05_a00001G00678V1.1")]
data$Do1_05_a00001G00678V1.1 <- as.numeric(data$Do1_05_a00001G00678V1.1)
png("./plots/W_C_RNA_Do1_05_a00001G00678V1.1.png", width = 1000, height = 700)
print(
  ggplot(data, aes(x = site, y = Do1_05_a00001G00678V1.1, fill = treatment_type)) +
    geom_boxplot(position = "dodge") +  # Use geom_boxplot instead of geom_bar
    labs(title = "Treatment Values by Site and Treatment Type", 
         x = "Site", y = "Treatment Value") +
    scale_fill_manual(values = c("W" = "red", "C" = "skyblue")) +  # Customize colors for W and C
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability
)
dev.off()


#---------------------------
#Do1_04_a00002G00429V1.1

subdat1<- data[,c("site","treatment_type", "Do1_04_a00002G00429V1.1")]
data$Do1_04_a00002G00429V1.1 <- as.numeric(data$Do1_04_a00002G00429V1.1)
png("./plots/W_C_RNA_Do1_04_a00002G00429V1.1.png", width = 1000, height = 700)
print(
  ggplot(data, aes(x = site, y = Do1_04_a00002G00429V1.1, fill = treatment_type)) +
    geom_boxplot(position = "dodge") +  # Use geom_boxplot instead of geom_bar
    labs(title = "Treatment Values by Site and Treatment Type", 
         x = "Site", y = "Treatment Value") +
    scale_fill_manual(values = c("W" = "red", "C" = "skyblue")) +  # Customize colors for W and C
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability
)
dev.off()

#---------------------------
#Do1_01_a00007G00164V1.1

subdat1<- data[,c("site","treatment_type", "Do1_01_a00007G00164V1.1")]
data$Do1_01_a00007G00164V1.1 <- as.numeric(data$Do1_01_a00007G00164V1.1)
png("./plots/W_C_RNA_Do1_01_a00007G00164V1.1.png", width = 1000, height = 700)
print(
  ggplot(data, aes(x = site, y = Do1_01_a00007G00164V1.1, fill = treatment_type)) +
    geom_boxplot(position = "dodge") +  # Use geom_boxplot instead of geom_bar
    labs(title = "Treatment Values by Site and Treatment Type", 
         x = "Site", y = "Treatment Value") +
    scale_fill_manual(values = c("W" = "red", "C" = "skyblue")) +  # Customize colors for W and C
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability
)
dev.off()

#---------------------------
#Do1_03_a00002G01336V1.1

subdat1<- data[,c("site","treatment_type", "Do1_03_a00002G01336V1.1")]
data$Do1_03_a00002G01336V1.1 <- as.numeric(data$Do1_03_a00002G01336V1.1)
png("./plots/W_C_RNA_Do1_03_a00002G01336V1.1.png", width = 1000, height = 700)
print(
  ggplot(data, aes(x = site, y = Do1_03_a00002G01336V1.1, fill = treatment_type)) +
    geom_point(size=3)+
    geom_boxplot(position = "dodge") +  # Use geom_boxplot instead of geom_bar
    labs(title = "Treatment Values by Site and Treatment Type", 
         x = "Site", y = "Treatment Value") +
    scale_fill_manual(values = c("W" = "red", "C" = "skyblue")) +  # Customize colors for W and C
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability
)
dev.off()

#---------------------------
#Do1_03_a00001G00267V1.1

subdat1<- data[,c("site","treatment_type", "Do1_03_a00001G00267V1.1")]
data$Do1_03_a00001G00267V1.1 <- as.numeric(data$Do1_03_a00001G00267V1.1)
png("./plots/W_C_RNA_Do1_03_a00001G00267V1.1.png", width = 1000, height = 700)
print(
  ggplot(data, aes(x = site, y = Do1_03_a00001G00267V1.1, fill = treatment_type)) +
    geom_boxplot(position = "dodge") +  # Use geom_boxplot instead of geom_bar
    labs(title = "Treatment Values by Site and Treatment Type", 
         x = "Site", y = "Treatment Value") +
    scale_fill_manual(values = c("W" = "red", "C" = "skyblue")) +  # Customize colors for W and C
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability
)
dev.off()

#---------------------------
#Do1_02_a00003G01059V1.1

subdat1<- data[,c("site","treatment_type", "Do1_02_a00003G01059V1.1")]
data$Do1_02_a00003G01059V1.1 <- as.numeric(data$Do1_02_a00003G01059V1.1)
png("./plots/W_C_RNA_Do1_02_a00003G01059V1.1.png", width = 1000, height = 700)
print(
  ggplot(data, aes(x = site, y = Do1_02_a00003G01059V1.1, fill = treatment_type)) +
    geom_boxplot(position = "dodge") +  # Use geom_boxplot instead of geom_bar
    labs(title = "Treatment Values by Site and Treatment Type", 
         x = "Site", y = "Treatment Value") +
    scale_fill_manual(values = c("W" = "red", "C" = "skyblue")) +  # Customize colors for W and C
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability
)
dev.off()


#######################################
# Do1_01_a00001G00358V1.1

subdat1<- data[,c("site","treatment_type", "Do1_01_a00001G00358V1.1")]
data$Do1_01_a00001G00358V1.1 <- as.numeric(data$Do1_01_a00001G00358V1.1)
png("./plots/W_C_RNA_Do1_01_a00001G00358V1.1.png", width = 1000, height = 700)
print(
  ggplot(data, aes(x = site, y = Do1_01_a00001G00358V1.1, fill = treatment_type)) +
    geom_boxplot(position = "dodge") +  # Use geom_boxplot instead of geom_bar
    labs(title = "Treatment Values by Site and Treatment Type", 
         x = "Site", y = "Treatment Value") +
    scale_fill_manual(values = c("W" = "red", "C" = "skyblue")) +  # Customize colors for W and C
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability
)
dev.off()


#----------------------------------
# Do1_05_a00001G00328V1.1

subdat1<- data[,c("site","treatment_type", "Do1_05_a00001G00328V1.1")]
data$Do1_05_a00001G00328V1.1 <- as.numeric(data$Do1_05_a00001G00328V1.1)
png("./plots/W_C_RNA_Do1_05_a00001G00328V1.1.png", width = 1000, height = 700)
print(
  ggplot(data, aes(x = site, y = Do1_05_a00001G00328V1.1, fill = treatment_type)) +
    geom_boxplot(position = "dodge") +  # Use geom_boxplot instead of geom_bar
    labs(title = "Treatment Values by Site and Treatment Type", 
         x = "Site", y = "Treatment Value") +
    scale_fill_manual(values = c("W" = "red", "C" = "skyblue")) +  # Customize colors for W and C
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability
)
dev.off()


#----------------------------------
#Do1_02_a00003G00097V1.1
grep Do1_02_a00003G00097V1.1 Dryas_octopetala_H1.csv
# Do1_02_a00003G00097V1.1,"glutathione S-transferase, amine-terminal domain protein;phi class glutathione S-transferase",7,1,0,4,2,"AT2G47730,AT1G02920,AT1G02930,AT2G02930,AT4G02520,AT1G02940,AT1G02950",LotjaGi5g1v0252100,,"Medtr1g088825,Medtr1g088840,Medtr1g088845,Medtr1g492670","FvH4_2g25200,FvH4_2g25210"
# https://www.arabidopsis.org/locus?key=28997
# EARLY RESPONSIVE TO DEHYDRATION 11

subdat1<- data[,c("site","treatment_type", "Do1_02_a00003G00097V1.1")]
data$Do1_02_a00003G00097V1.1 <- as.numeric(data$Do1_02_a00003G00097V1.1)
png("./plots/W_C_RNA_Do1_02_a00003G00097V1.1.png", width = 1000, height = 700)
print(
  ggplot(data, aes(x = site, y = Do1_02_a00003G00097V1.1, fill = treatment_type)) +
    geom_boxplot(position = "dodge") +  # Use geom_boxplot instead of geom_bar
    labs(title = "Treatment Values by Site and Treatment Type", 
         x = "Site", y = "Treatment Value") +
    scale_fill_manual(values = c("W" = "red", "C" = "skyblue")) +  # Customize colors for W and C
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability
)
dev.off()
#----------------------------------
#Do1_06_a00002G00149V1.1

subdat1<- data[,c("site","treatment_type", "Do1_06_a00002G00149V1.1")]
data$Do1_06_a00002G00149V1.1 <- as.numeric(data$Do1_06_a00002G00149V1.1)
png("./plots/W_C_RNA_Do1_06_a00002G00149V1.1.png", width = 1000, height = 700)
print(
  ggplot(data, aes(x = site, y = Do1_06_a00002G00149V1.1, fill = treatment_type)) +
    geom_boxplot(position = "dodge") +  # Use geom_boxplot instead of geom_bar
    labs(title = "Treatment Values by Site and Treatment Type", 
         x = "Site", y = "Treatment Value") +
    scale_fill_manual(values = c("W" = "red", "C" = "skyblue")) +  # Customize colors for W and C
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability
)
dev.off()

#----------------------------------
#Do1_01_a00001G02208V1.1

subdat1<- data[,c("site","treatment_type", "Do1_01_a00001G02208V1.1")]
data$Do1_01_a00001G02208V1.1 <- as.numeric(data$Do1_01_a00001G02208V1.1)
png("./plots/W_C_RNA_Do1_01_a00001G02208V1.1.png", width = 1000, height = 700)
print(
  ggplot(data, aes(x = site, y = Do1_01_a00001G02208V1.1, fill = treatment_type)) +
    geom_boxplot(position = "dodge") +  # Use geom_boxplot instead of geom_bar
    labs(title = "Treatment Values by Site and Treatment Type", 
         x = "Site", y = "Treatment Value") +
    scale_fill_manual(values = c("W" = "red", "C" = "skyblue")) +  # Customize colors for W and C
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability
)
dev.off()

##################################
# Add interproscan data

cd /home/celphin/projects/rpp-rieseber/celphin/Dryas/Reference_genomes/annotation

grep Do1_05_a00001G00678V1.1 Dryas_octopetala_H1.csv
#Do1_05_a00001G00678V1.1,TCP family transcription factor,1,1,1,1,1,AT3G47620,LotjaGi4g1v0248700,Ro05_G03227,Medtr4g108370,FvH4_5g15150

grep Do1_02_a00003G00097V1.1 Dryas_octopetala_H1.csv
#Do1_02_a00003G00097V1.1,"glutathione S-transferase, amine-terminal domain protein;phi class glutathione S-transferase",7,1,0,4,2,"AT2G47730,AT1G02920,AT1G02930,AT2G02930,AT4G02520,AT1G02940,AT1G02950",LotjaGi5g1v0252100,,"Medtr1g088825,Medtr1g088840,Medtr1g088845,Medtr1g492670","FvH4_2g25200,FvH4_2g25210"



grep Do1_04_a00002G00429V1.1 Dryas_octopetala_H1.csv # plant development and defense
grep Do1_01_a00007G00164V1.1 Dryas_octopetala_H1.csv
grep Do1_03_a00002G01336V1.1 Dryas_octopetala_H1.csv # linked to stress-induced responses and programmed cell death
grep Do1_03_a00001G00267V1.1 Dryas_octopetala_H1.csv
grep Do1_02_a00003G01059V1.1 Dryas_octopetala_H1.csv # veg growth

Do1_04_a00002G00429V1.1,germin family protein,19,5,8,13,8,"AT4G14630,AT5G38940,AT5G39130,AT5G39160,AT5G39180,AT3G04180,AT3G04190,AT5G38930,AT5G39190,AT3G04200,AT5G39150,AT3G05950,AT5G38910,AT5G38960,AT3G04170,AT3G04150,AT5G39120,AT5G39100,AT5G39110","LotjaGi1g1v0620000,LotjaGi1g1v0620100,LotjaGi2g1v0242900,LotjaGi2g1v0243100,LotjaGi2g1v0243300","Ro05_G31491,Ro05_G31492,Ro05_G31493,Ro05_G31494,Ro05_G31495,Ro06_G28962,Ro06_G28965,Ro06_G28966","Medtr1g079490,Medtr2g019250,Medtr4g017030,Medtr4g017040,Medtr4g017050,Medtr5g046410,Medtr5g046430,Medtr6g005310,Medtr6g005330,Medtr6g005340,Medtr6g005350,Medtr6g005360,Medtr6g005380","FvH4_2g10250,FvH4_3g20400,FvH4_5g18050,FvH4_5g18070,FvH4_5g18440,FvH4_5g18450,FvH4_5g18460,FvH4_5g18470"
Do1_01_a00007G00164V1.1,,0,0,0,0,0,,,,,
Do1_03_a00002G01336V1.1,pirin-like plant protein,1,1,1,1,1,AT1G50590,LotjaGi1g1v0021700,Ro03_G24134,Medtr5g094250,FvH4_3g23750
Do1_03_a00001G00267V1.1,,0,1,1,0,1,,LotjaGi1g1v0015300,Ro03_G14687,,FvH4_3g42750
Do1_02_a00003G01059V1.1,high-affinity nickel-transport family protein,2,2,1,1,1,"AT2G16800,AT4G35080","LotjaGi1g1v0558300,LotjaGi3g1v0519700",Ro05_G25796,Medtr4g073010,FvH4_5g30360
# https://bar.utoronto.ca/thalemine/report.do?id=5448674



###################################
# on local machine
cd /home/Owner/Desktop
scp -v celphin@cedar.computecanada.ca:/home/celphin/projects/rpp-rieseber/celphin/Dryas/RNAseq_analysis/plots/*png .
scp -v celphin@cedar.computecanada.ca:/home/celphin/projects/rpp-rieseber/celphin/Dryas/RNAseq_analysis/output/*pdf .
scp -v celphin@cedar.computecanada.ca:/home/celphin/projects/rpp-rieseber/celphin/Dryas/RNAseq_analysis/RNA_heatmap_ordered.png .
scp -v celphin@cedar.computecanada.ca:/home/celphin/projects/rpp-rieseber/celphin/Dryas/Reference_genomes/annotation/Dryas_octopetala_H1.xlsx .
scp -v celphin@cedar.computecanada.ca:/home/celphin/projects/rpp-rieseber/celphin/Dryas/RNAseq_analysis/Dryas_octopetala_H1.supercontigs.fa .
scp -v celphin@cedar.computecanada.ca:/home/celphin/projects/rpp-rieseber/celphin/Dryas/RNAseq_analysis/Dryas_octopetala_H1.gff3 .

mkdir LAT_DMRs
scp -v celphin@cedar.computecanada.ca:/home/celphin/projects/rpp-rieseber/celphin/Dryas/DMRs/June2022_Metilene_DMRs/LAT_DMRs* .
scp -v celphin@cedar.computecanada.ca:/home/celphin/projects/rpp-rieseber/celphin/Dryas/RNAseq_analysis/*RNA_Do1_01_a00001G00358V1.1.png .
scp -v celphin@cedar.computecanada.ca:/home/celphin/projects/rpp-rieseber/celphin/Dryas/RNAseq_analysis/*RNA_Do1_05_a00001G00328V1.1.png .
scp -v celphin@cedar.computecanada.ca:/home/celphin/projects/rpp-rieseber/celphin/Dryas/RNAseq_analysis/W_C_*png .
scp -v celphin@cedar.computecanada.ca:/home/celphin/projects/rpp-rieseber/celphin/Other/biol525d.tar.gz

###############################
# transcript count in transcriptome
grep -c ">" Dryas_octopetala_H1.transcript.fa
41181

# 28218 genes had 10 or more reads map in Lat
# 16716 genes had 200 or more reads map in LAT

##########################################
# Count the number of DEGs that overlap between sites and wtih seedlings

# Extract the first column from both files
cut -f1 RNA_Alaska_W_C_DERs.txt | sort > Alaska_first_column_sorted.txt
cut -f1 RNA_Sweden_W_C_DERs.txt | sort > Swed_first_column_sorted.txt
cut -f1 RNA_Seedling_W_C_DERs.txt | sort > SE_first_column_sorted.txt
cut -f1 RNA_Alex_W_C_DERs.txt | sort > Alex_first_column_sorted.txt
cut -f1 RNA_Norway_W_C_DERs.txt | sort > Norway_first_column_sorted.txt

# Use 'comm' to find the common values and count them
common_count=$(comm -12 Alaska_first_column_sorted.txt Swed_first_column_sorted.txt | wc -l)
echo $common_count

common_count=$(comm -12 SE_first_column_sorted.txt Swed_first_column_sorted.txt | wc -l)
echo $common_count

common_count=$(comm -12 Alex_first_column_sorted.txt Swed_first_column_sorted.txt | wc -l)
echo $common_count

common_count=$(comm -12 Norway_first_column_sorted.txt Swed_first_column_sorted.txt | wc -l)
echo $common_count

common_count=$(comm -12 SE_first_column_sorted.txt Alaska_first_column_sorted.txt | wc -l)
echo $common_count

common_count=$(comm -12 Alex_first_column_sorted.txt Alaska_first_column_sorted.txt | wc -l)
echo $common_count

common_count=$(comm -12 Norway_first_column_sorted.txt Alaska_first_column_sorted.txt | wc -l)
echo $common_count


# Swed/Seedling 2
# Swed/Alaska 4
# Alex Swed 1
# Norway/Swed 3

# SE/Alaska 1
# Alex/Alas 2
# Norw/Alas 1

common_count=$(comm -12 Alaska_first_column_sorted.txt SE_first_column_sorted.txt | wc -l)
echo $common_count

common_count=$(comm -12 Swed_first_column_sorted.txt SE_first_column_sorted.txt | wc -l)
echo $common_count

common_count=$(comm -12 Alex_first_column_sorted.txt SE_first_column_sorted.txt | wc -l)
echo $common_count

common_count=$(comm -12 Norway_first_column_sorted.txt SE_first_column_sorted.txt | wc -l)
echo $common_count

1
2 Do1_05_a00001G00111V1.1 NADH OXIDOREDUCTASE
1
1

cat Alaska_first_column_sorted.txt Swed_first_column_sorted.txt SE_first_column_sorted.txt \
Alex_first_column_sorted.txt Norway_first_column_sorted.txt > All_DEGs.txt

sort All_DEGs.txt | uniq -c | awk '$1 > 1'

sort All_DEGs.txt | uniq -c | awk '$1 > 1'
      2 Do1_00148G00016V1.1
      2 Do1_00148G00027V1.1
      2 Do1_00148G00037V1.1
      2 Do1_01_a00001G01424V1.1
      2 Do1_04_a00002G00133V1.1
      2 Do1_04_a00002G00429V1.1
      2 Do1_05_a00001G00111V1.1
      2 Do1_05_a00003G00085V1.1
      5 logFC


##################################
# Intersect the DMR genes for CHH and CpG and DEGs for each site

cd /home/celphin/scratch/Dryas/RNAseq_analysis
mkdir DMR_DEG_comparison
cd DMR_DEG_comparison

# copy over CHH DMR genes
cp /home/celphin/scratch/Dryas/snpEff/Dryas_DMRs/CHH/*_CHH_W_C.bed.out .

# copy CpG DMR genes
cp /home/celphin/scratch/Dryas/snpEff/out/*_W_C.bed.out .

# copy over the DEG significant lists
cp ../RNA_Alaska_W_C_DERs.txt .
cp ../RNA_Svalbard_W_C_DERs.txt .
cp ../RNA_Sweden_W_C_DERs.txt .
cp ../RNA_Nunavut_W_C_DERs.txt .

cd  /home/celphin/scratch/Dryas/RNAseq_analysis/DMR_DEG_comparison

#----------------------------
# Nunavut
awk '{print $1}' RNA_Alex_W_C_DERs.txt |wc -l
# 42

for value in $(awk 'NR > 1 {print $1}' RNA_Alex_W_C_DERs.txt ); do
    echo "Searching for $value in *.out files..."
    # Search for the Nunavut_DEGs_overlap_DMRs value in all *.out files
    grep -H "$value" Nunavut_*.out | awk -v term="$value" '{print term ": " $0}' >> Nunavut_DEGs_overlap_DMRs.txt
done
wc -l Nunavut_DEGs_overlap_DMRs.txt
#4 

awk '{print $1}' Nunavut_DEGs_overlap_DMRs.txt | sort -u | wc -l
# 2

cat Nunavut_DEGs_overlap_DMRs.txt
# Do1_01_a00002G00362V1.1: Nunavut_CHH_W_C.bed.out:Do1_01_a00002  2172057 2172229 line_131;Intergenic:Do1_01_a00002G00362-Do1_01_a00002G00363;Upstream:2590|Transcript:Do1_01_a00002G00362V1.1:protein_coding|Gene:Do1_01_a00002G00362:protein_coding
# Do1_01_a00002G00362V1.1: Nunavut_W_C.bed.out:Do1_01_a00002      2171563 2172008 line_27;Intergenic:Do1_01_a00002G00362-Do1_01_a00002G00363;Upstream:2096|Transcript:Do1_01_a00002G00362V1.1:protein_coding|Gene:Do1_01_a00002G00362:protein_coding
# Do1_01_a00002G00362V1.1: Nunavut_W_C.bed.out:Do1_01_a00002      2172008 2172207 line_28;Upstream:2541|Transcript:Do1_01_a00002G00362V1.1:protein_coding|Gene:Do1_01_a00002G00362:protein_coding;Intergenic:Do1_01_a00002G00362-Do1_01_a00002G00363
# Do1_04_a00001G01911V1.1: Nunavut_CHH_W_C.bed.out:Do1_04_a00001  11058419        11058547        line_635;Upstream:1876|Transcript:Do1_04_a00001G01910V1.1:protein_coding|Gene:Do1_04_a00001G01910:protein_coding;Upstream:324|Transcript:Do1_04_a00001G01911V1.1:protein_coding|Gene:Do1_04_a00001G01911:protein_coding;Intergenic:Do1_04_a00001G01910-Do1_04_a00001G01911

# Search GO terms
sed -i 's/V1.1://g' *_DEGs_overlap_DMRs.txt

for value in $(awk 'NR > 1 {print $1}' Nunavut_DEGs_overlap_DMRs.txt | sort -u ); do
    echo "Searching for $value in GO term file..."
    grep -H "$value"  interproscan_dryas_full3.tsv >> GO_terms_Nunavut_DEGs_overlap_DMRs.txt
done



#------------------------
# Sweden
awk '{print $1}'  RNA_Sweden_W_C_DERs.txt |wc -l
# 99
for value in $(awk 'NR > 1 {print $1}' RNA_Sweden_W_C_DERs.txt ); do
    echo "Searching for $value in *.out files..."
    # Search for the Sweden_DEGs_overlap_DMRslue in all *.out files
    grep -H "$value" Sweden_*.out | awk -v term="$value" '{print term ": " $0}' >> Sweden_DEGs_overlap_DMRs.txt
done
wc -l Sweden_DEGs_overlap_DMRs.txt
23
awk '{print $1}' Sweden_DEGs_overlap_DMRs.txt| sort -u | wc -l
21

Do1_00161G00004V1.1: Sweden_W_C.bed.out:Do1_00161       21292   21584   line_12;Upstream:3020|Transcript:Do1_00161G00004V1.1:protein_coding|Gene:Do1_00161G00004:protein_coding;Exon:2:2:RETAINED|Transcript:Do1_00161G00005V1.1:protein_coding|Gene:Do1_00161G00005:protein_coding;Intergenic:Do1_00161G00004-Do1_00161G00005;Upstream:1804|Transcript:Do1_00161G00006V1.1:protein_coding|Gene:Do1_00161G00006:protein_coding;Downstream:176|Transcript:Do1_00161G00005V1.1:protein_coding|Gene:Do1_00161G00005:protein_coding
Do1_01_a00001G00214V1.1: Sweden_W_C.bed.out:Do1_01_a00001       1059564 1059798 line_38;Intron:7:7:RETAINED-RETAINED|Transcript:Do1_01_a00001G00216V1.1:protein_coding|Gene:Do1_01_a00001G00216:protein_coding;Downstream:3463|Transcript:Do1_01_a00001G00214V1.1:protein_coding|Gene:Do1_01_a00001G00214:protein_coding;Downstream:2692|Transcript:Do1_01_a00001G00215V1.1:protein_coding|Gene:Do1_01_a00001G00215:protein_coding
Do1_01_a00001G00214V1.1: Sweden_W_C.bed.out:Do1_01_a00001       1051991 1052562 line_39;Downstream:486|Transcript:Do1_01_a00001G00213V1.1:protein_coding|Gene:Do1_01_a00001G00213:protein_coding;Downstream:5227|Transcript:Do1_01_a00001G00216V1.1:protein_coding|Gene:Do1_01_a00001G00216:protein_coding;Upstream:443|Transcript:Do1_01_a00001G00214V1.1:protein_coding|Gene:Do1_01_a00001G00214:protein_coding;Upstream:4557|Transcript:Do1_01_a00001G00215V1.1:protein_coding|Gene:Do1_01_a00001G00215:protein_coding;Intergenic:Do1_01_a00001G00213-Do1_01_a00001G00214;Exon:1:5:RETAINED|Transcript:Do1_01_a00001G00214V1.1:protein_coding|Gene:Do1_01_a00001G00214:protein_coding
Do1_01_a00001G02184V1.1: Sweden_CHH_W_C.bed.out:Do1_01_a00001   12235130        12235536        line_141;Downstream:4542|Transcript:Do1_01_a00001G02184V1.1:protein_coding|Gene:Do1_01_a00001G02184:protein_coding;Downstream:1414|Transcript:Do1_01_a00001G02185V1.1:protein_coding|Gene:Do1_01_a00001G02185:protein_coding;Intergenic:Do1_01_a00001G02185-Do1_01_a00001G02186
Do1_01_a00003G00798V1.1: Sweden_W_C.bed.out:Do1_01_a00003       4988127 4988265 line_191;Downstream:3389|Transcript:Do1_01_a00003G00798V1.1:protein_coding|Gene:Do1_01_a00003G00798:protein_coding;Intergenic:Do1_01_a00003G00799-Do1_01_a00003G00800;Downstream:279|Transcript:Do1_01_a00003G00799V1.1:protein_coding|Gene:Do1_01_a00003G00799:protein_coding;Exon:1:4:RETAINED|Transcript:Do1_01_a00003G00800V1.1:protein_coding|Gene:Do1_01_a00003G00800:protein_coding;Upstream:121|Transcript:Do1_01_a00003G00800V1.1:protein_coding|Gene:Do1_01_a00003G00800:protein_coding
Do1_02_a00001G00643V1.1: Sweden_CHH_W_C.bed.out:Do1_02_a00001   4204748 4204988 line_431;Intergenic:Do1_02_a00001G00642-Do1_02_a00001G00643;Downstream:1578|Transcript:Do1_02_a00001G00643V1.1:protein_coding|Gene:Do1_02_a00001G00643:protein_coding
Do1_02_a00001G00643V1.1: Sweden_W_C.bed.out:Do1_02_a00001       4213499 4213622 line_319;Exon:1:29:RETAINED|Transcript:Do1_02_a00001G00643V1.1:protein_coding|Gene:Do1_02_a00001G00643:protein_coding
Do1_02_a00002G00086V1.1: Sweden_CHH_W_C.bed.out:Do1_02_a00002   399374  399652  line_492;Downstream:2737|Transcript:Do1_02_a00002G00087V1.1:protein_coding|Gene:Do1_02_a00002G00087:protein_coding;Upstream:1957|Transcript:Do1_02_a00002G00086V1.1:protein_coding|Gene:Do1_02_a00002G00086:protein_coding;Downstream:3858|Transcript:Do1_02_a00002G00088V1.1:protein_coding|Gene:Do1_02_a00002G00088:protein_coding;Intergenic:Do1_02_a00002G00086-Do1_02_a00002G00087
Do1_02_a00003G01189V1.1: Sweden_CHH_W_C.bed.out:Do1_02_a00003   5957427 5957522 line_538;Intergenic:Do1_02_a00003G01188-Do1_02_a00003G01189;Downstream:675|Transcript:Do1_02_a00003G01188V1.1:protein_coding|Gene:Do1_02_a00003G01188:protein_coding;Upstream:2341|Transcript:Do1_02_a00003G01189V1.1:protein_coding|Gene:Do1_02_a00003G01189:protein_coding
Do1_02_a00003G01269V1.1: Sweden_CHH_W_C.bed.out:Do1_02_a00003   6330996 6331243 line_544;Intergenic:Do1_02_a00003G01268-Do1_02_a00003G01269;Upstream:4022|Transcript:Do1_02_a00003G01269V1.1:protein_coding|Gene:Do1_02_a00003G01269:protein_coding
Do1_03_a00002G01515V1.1: Sweden_W_C.bed.out:Do1_03_a00002       8261778 8261888 line_631;Intergenic:Do1_03_a00002G01515-Do1_03_a00002G01516;Downstream:959|Transcript:Do1_03_a00002G01515V1.1:protein_coding|Gene:Do1_03_a00002G01515:protein_coding
Do1_05_a00001G00111V1.1: Sweden_CHH_W_C.bed.out:Do1_05_a00001   558703  558914  line_1286;Downstream:979|Transcript:Do1_05_a00001G00111V1.1:protein_coding|Gene:Do1_05_a00001G00111:protein_coding;Intergenic:Do1_05_a00001G00110-Do1_05_a00001G00111
Do1_05_a00001G01731V1.1: Sweden_W_C.bed.out:Do1_05_a00001       8858012 8858374 line_1036;Upstream:5065|Transcript:Do1_05_a00001G01732V1.1:protein_coding|Gene:Do1_05_a00001G01732:protein_coding;Downstream:3525|Transcript:Do1_05_a00001G01730V1.1:protein_coding|Gene:Do1_05_a00001G01730:protein_coding;Exon:1:1:RETAINED|Transcript:Do1_05_a00001G01731V1.1:protein_coding|Gene:Do1_05_a00001G01731:protein_coding
Do1_05_a00001G01744V1.1: Sweden_W_C.bed.out:Do1_05_a00001       8941886 8942095 line_1037;Upstream:2072|Transcript:Do1_05_a00001G01742V1.1:protein_coding|Gene:Do1_05_a00001G01742:protein_coding;Intergenic:Do1_05_a00001G01743-Do1_05_a00001G01744;Upstream:263|Transcript:Do1_05_a00001G01743V1.1:protein_coding|Gene:Do1_05_a00001G01743:protein_coding;Upstream:3998|Transcript:Do1_05_a00001G01744V1.1:protein_coding|Gene:Do1_05_a00001G01744:protein_coding
Do1_05_a00002G00016V1.1: Sweden_W_C.bed.out:Do1_05_a00002       107390  107696  line_1078;Exon:7:7:RETAINED|Transcript:Do1_05_a00002G00014V1.1:protein_coding|Gene:Do1_05_a00002G00014:protein_coding;Upstream:1998|Transcript:Do1_05_a00002G00015V1.1:protein_coding|Gene:Do1_05_a00002G00015:protein_coding;Upstream:4612|Transcript:Do1_05_a00002G00016V1.1:protein_coding|Gene:Do1_05_a00002G00016:protein_coding
Do1_05_a00003G00079V1.1: Sweden_CHH_W_C.bed.out:Do1_05_a00003   402294  402497  line_1458;Intergenic:Do1_05_a00003G00080-Do1_05_a00003G00081;Downstream:722|Transcript:Do1_05_a00003G00080V1.1:protein_coding|Gene:Do1_05_a00003G00080:protein_coding;Upstream:3541|Transcript:Do1_05_a00003G00079V1.1:protein_coding|Gene:Do1_05_a00003G00079:protein_coding;Downstream:5199|Transcript:Do1_05_a00003G00082V1.1:protein_coding|Gene:Do1_05_a00003G00082:protein_coding;Downstream:1756|Transcript:Do1_05_a00003G00081V1.1:protein_coding|Gene:Do1_05_a00003G00081:protein_coding
Do1_05_a00003G00082V1.1: Sweden_CHH_W_C.bed.out:Do1_05_a00003   402294  402497  line_1458;Intergenic:Do1_05_a00003G00080-Do1_05_a00003G00081;Downstream:722|Transcript:Do1_05_a00003G00080V1.1:protein_coding|Gene:Do1_05_a00003G00080:protein_coding;Upstream:3541|Transcript:Do1_05_a00003G00079V1.1:protein_coding|Gene:Do1_05_a00003G00079:protein_coding;Downstream:5199|Transcript:Do1_05_a00003G00082V1.1:protein_coding|Gene:Do1_05_a00003G00082:protein_coding;Downstream:1756|Transcript:Do1_05_a00003G00081V1.1:protein_coding|Gene:Do1_05_a00003G00081:protein_coding
Do1_05_a00003G00085V1.1: Sweden_CHH_W_C.bed.out:Do1_05_a00003   421747  421972  line_1459;Intergenic:Do1_05_a00003G00084-Do1_05_a00003G00085;Upstream:1944|Transcript:Do1_05_a00003G00084V1.1:protein_coding|Gene:Do1_05_a00003G00084:protein_coding;Downstream:4995|Transcript:Do1_05_a00003G00085V1.1:protein_coding|Gene:Do1_05_a00003G00085:protein_coding
Do1_05_a00003G00319V1.1: Sweden_W_C.bed.out:Do1_05_a00003       1787043 1787199 line_1147;Exon:4:8:RETAINED|Transcript:Do1_05_a00003G00320V1.1:protein_coding|Gene:Do1_05_a00003G00320:protein_coding;Upstream:3699|Transcript:Do1_05_a00003G00321V1.1:protein_coding|Gene:Do1_05_a00003G00321:protein_coding;Downstream:4788|Transcript:Do1_05_a00003G00319V1.1:protein_coding|Gene:Do1_05_a00003G00319:protein_coding
Do1_06_a00001G01083V1.1: Sweden_CHH_W_C.bed.out:Do1_06_a00001   6339632 6339724 line_1607;Downstream:4780|Transcript:Do1_06_a00001G01083V1.1:protein_coding|Gene:Do1_06_a00001G01083:protein_coding;Intergenic:Do1_06_a00001G01083-Do1_06_a00001G01084
Do1_06_a00001G01665V1.1: Sweden_CHH_W_C.bed.out:Do1_06_a00001   9827492 9827601 line_1638;Upstream:4644|Transcript:Do1_06_a00001G01665V1.1:protein_coding|Gene:Do1_06_a00001G01665:protein_coding;Intergenic:Do1_06_a00001G01664-Do1_06_a00001G01665;Upstream:395|Transcript:Do1_06_a00001G01664V1.1:protein_coding|Gene:Do1_06_a00001G01664:protein_coding;Downstream:1328|Transcript:Do1_06_a00001G01663V1.1:protein_coding|Gene:Do1_06_a00001G01663:protein_coding
Do1_07_a00002G01553V1.1: Sweden_W_C.bed.out:Do1_07_a00002       7853137 7853522 line_1515;Exon:1:1:RETAINED|Transcript:Do1_07_a00002G01552V1.1:protein_coding|Gene:Do1_07_a00002G01552:protein_coding;Upstream:3802|Transcript:Do1_07_a00002G01550V1.1:protein_coding|Gene:Do1_07_a00002G01550:protein_coding;Upstream:1230|Transcript:Do1_07_a00002G01551V1.1:protein_coding|Gene:Do1_07_a00002G01551:protein_coding;Upstream:4671|Transcript:Do1_07_a00002G01549V1.1:protein_coding|Gene:Do1_07_a00002G01549:protein_coding;Upstream:0|Transcript:Do1_07_a00002G01552V1.1:protein_coding|Gene:Do1_07_a00002G01552:protein_coding;Downstream:2042|Transcript:Do1_07_a00002G01553V1.1:protein_coding|Gene:Do1_07_a00002G01553:protein_coding;Intergenic:Do1_07_a00002G01552-Do1_07_a00002G01553
Do1_07_a00004G00048V1.1: Sweden_CHH_W_C.bed.out:Do1_07_a00004   287837  287900  line_1991;Upstream:1650|Transcript:Do1_07_a00004G00048V1.1:protein_coding|Gene:Do1_07_a00004G00048:protein_coding;Intergenic:Do1_07_a00004G00047-Do1_07_a00004G00048

# Search GO terms

for value in $(awk 'NR > 1 {print $1}' Sweden_DEGs_overlap_DMRs.txt | sort -u ); do
    echo "Searching for $value in GO term file..."
    grep -H "$value"  interproscan_dryas_full3.tsv >> GO_terms_Sweden_DEGs_overlap_DMRs.txt
done


#----------------------------
# Alaska
awk '{print $1}' RNA_Alaska_W_C_DERs.txt |wc -l
# 65

for value in $(awk 'NR > 1 {print $1}' RNA_Alaska_W_C_DERs.txt ); do
    echo "Searching for $value in *.out files..."
    # Search for the Alaska_DEGs_overlap_DMRslue in all *.out files
    grep -H "$value" Alaska_*.out | awk -v term="$value" '{print term ": " $0}'  >> Alaska_DEGs_overlap_DMRs.txt
done
wc -l Alaska_DEGs_overlap_DMRs.txt
# 7
awk '{print $1}' Alaska_DEGs_overlap_DMRs.txt| sort -u | wc -l
# 7
Do1_01_a00001G00762V1.1: Alaska_W_C.bed.out:Do1_01_a00001       3955238 3955782 line_21;Gene:Do1_01_a00001G00761:protein_coding;Upstream:1141|Transcript:Do1_01_a00001G00760V1.1:protein_coding|Gene:Do1_01_a00001G00760:protein_coding;Downstream:3187|Transcript:Do1_01_a00001G00762V1.1:protein_coding|Gene:Do1_01_a00001G00762:protein_coding
Do1_01_a00001G01048V1.1: Alaska_W_C.bed.out:Do1_01_a00001       5607706 5607853 line_32;Downstream:3699|Transcript:Do1_01_a00001G01050V1.1:protein_coding|Gene:Do1_01_a00001G01050:protein_coding;Downstream:2575|Transcript:Do1_01_a00001G01047V1.1:protein_coding|Gene:Do1_01_a00001G01047:protein_coding;Upstream:2761|Transcript:Do1_01_a00001G01049V1.1:protein_coding|Gene:Do1_01_a00001G01049:protein_coding;Intron:4:4:RETAINED-RETAINED|Transcript:Do1_01_a00001G01048V1.1:protein_coding|Gene:Do1_01_a00001G01048:protein_coding
Do1_02_a00001G00594V1.1: Alaska_CHH_W_C.bed.out:Do1_02_a00001   3950478 3950584 line_409;Upstream:1688|Transcript:Do1_02_a00001G00595V1.1:protein_coding|Gene:Do1_02_a00001G00595:protein_coding;Upstream:714|Transcript:Do1_02_a00001G00596V1.1:protein_coding|Gene:Do1_02_a00001G00596:protein_coding;Upstream:4505|Transcript:Do1_02_a00001G00594V1.1:protein_coding|Gene:Do1_02_a00001G00594:protein_coding;Intergenic:Do1_02_a00001G00595-Do1_02_a00001G00596
Do1_02_a00004G01021V1.1: Alaska_CHH_W_C.bed.out:Do1_02_a00004   6382379 6382442 line_628;Downstream:1622|Transcript:Do1_02_a00004G01020V1.1:protein_coding|Gene:Do1_02_a00004G01020:protein_coding;Intergenic:Do1_02_a00004G01019-Do1_02_a00004G01020;Upstream:2653|Transcript:Do1_02_a00004G01021V1.1:protein_coding|Gene:Do1_02_a00004G01021:protein_coding
Do1_03_a00002G00493V1.1: Alaska_W_C.bed.out:Do1_03_a00002       2787693 2787975 line_324;Downstream:0|Transcript:Do1_03_a00002G00495V1.1:protein_coding|Gene:Do1_03_a00002G00495:protein_coding;Upstream:3201|Transcript:Do1_03_a00002G00494V1.1:protein_coding|Gene:Do1_03_a00002G00494:protein_coding;Intergenic:Do1_03_a00002G00494-Do1_03_a00002G00495;Exon:1:1:RETAINED|Transcript:Do1_03_a00002G00495V1.1:protein_coding|Gene:Do1_03_a00002G00495:protein_coding;Downstream:941|Transcript:Do1_03_a00002G00496V1.1:protein_coding|Gene:Do1_03_a00002G00496:protein_coding;Intergenic:Do1_03_a00002G00495-Do1_03_a00002G00496;Upstream:11|Transcript:Do1_03_a00002G00495V1.1:protein_coding|Gene:Do1_03_a00002G00495:protein_coding;Upstream:3800|Transcript:Do1_03_a00002G00493V1.1:protein_coding|Gene:Do1_03_a00002G00493:protein_coding
Do1_03_a00002G01930V1.1: Alaska_CHH_W_C.bed.out:Do1_03_a00002   10532006        10532093        line_798;Downstream:4828|Transcript:Do1_03_a00002G01928V1.1:protein_coding|Gene:Do1_03_a00002G01928:protein_coding;Downstream:4205|Transcript:Do1_03_a00002G01929V1.1:protein_coding|Gene:Do1_03_a00002G01929:protein_coding;Upstream:1979|Transcript:Do1_03_a00002G01930V1.1:protein_coding|Gene:Do1_03_a00002G01930:protein_coding;Intergenic:Do1_03_a00002G01930-Do1_03_a00002G01931
Do1_04_a00003G00270V1.1: Alaska_CHH_W_C.bed.out:Do1_04_a00003   1581541 1581932 line_1096;Upstream:3377|Transcript:Do1_04_a00003G00271V1.1:protein_coding|Gene:Do1_04_a00003G00271:protein_coding;Downstream:3839|Transcript:Do1_04_a00003G00269V1.1:protein_coding|Gene:Do1_04_a00003G00269:protein_coding;Downstream:3237|Transcript:Do1_04_a00003G00270V1.1:protein_coding|Gene:Do1_04_a00003G00270:protein_coding;Intergenic:Do1_04_a00003G00270-Do1_04_a00003G00271

# Search GO terms

for value in $(awk 'NR > 1 {print $1}' Alaska_DEGs_overlap_DMRs.txt | sort -u ); do
    echo "Searching for $value in GO term file..."
    grep -H "$value"  interproscan_dryas_full3.tsv >> GO_terms_Alaska_DEGs_overlap_DMRs.txt
done

#---------------------
# Svalbard
awk '{print $1}' RNA_Norway_W_C_DERs.txt |wc -l
# 15

for value in $(awk 'NR > 1 {print $1}' RNA_Norway_W_C_DERs.txt ); do
    echo "Searching for $value in *.out files..."
    # Search for the Norway_DEGs_overlap_DMRslue in all *.out files
    grep -H "$value" *.out | awk -v term="$value" '{print term ": " $0}' >> Norway_DEGs_overlap_DMRs.txt
done
wc -l Norway_DEGs_overlap_DMRs.txt
10
awk '{print $1}' Norway_DEGs_overlap_DMRs.txt| sort -u | wc -l 
8

Do1_00237G00005V1.1: Nunavut_CHH_W_C.bed.out:Do1_00237  33052   33158   line_15;Intergenic:Do1_00237G00004-Do1_00237G00005;Downstream:1428|Transcript:Do1_00237G00005V1.1:protein_coding|Gene:Do1_00237G00005:protein_coding
Do1_01_a00001G01424V1.1: Svalbard_CHH_W_C.bed.out:Do1_01_a00001 7631941 7632076 line_51;Exon:2:4:RETAINED|Transcript:Do1_01_a00001G01423V1.1:protein_coding|Gene:Do1_01_a00001G01423:protein_coding;Downstream:3704|Transcript:Do1_01_a00001G01424V1.1:protein_coding|Gene:Do1_01_a00001G01424:protein_coding;Upstream:4908|Transcript:Do1_01_a00001G01422V1.1:protein_coding|Gene:Do1_01_a00001G01422:protein_coding;Exon:1:4:RETAINED|Transcript:Do1_01_a00001G01423V1.1:protein_coding|Gene:Do1_01_a00001G01423:protein_coding
Do1_03_a00001G00132V1.1: Alaska_CHH_W_C.bed.out:Do1_03_a00001   666244  666398  line_665;Upstream:1576|Transcript:Do1_03_a00001G00131V1.1:protein_coding|Gene:Do1_03_a00001G00131:protein_coding;Intergenic:Do1_03_a00001G00131-Do1_03_a00001G00132;Downstream:3506|Transcript:Do1_03_a00001G00132V1.1:protein_coding|Gene:Do1_03_a00001G00132:protein_coding
Do1_03_a00001G00132V1.1: Svalbard_CHH_W_C.bed.out:Do1_03_a00001 666208  666405  line_489;Downstream:3542|Transcript:Do1_03_a00001G00132V1.1:protein_coding|Gene:Do1_03_a00001G00132:protein_coding;Intergenic:Do1_03_a00001G00131-Do1_03_a00001G00132;Upstream:1540|Transcript:Do1_03_a00001G00131V1.1:protein_coding|Gene:Do1_03_a00001G00131:protein_coding
Do1_04_a00001G01663V1.1: Alaska_CHH_W_C.bed.out:Do1_04_a00001   9749297 9749386 line_989;Intron:17:17:RETAINED-RETAINED|Transcript:Do1_04_a00001G01664V1.1:protein_coding|Gene:Do1_04_a00001G01664:protein_coding;Upstream:2764|Transcript:Do1_04_a00001G01662V1.1:protein_coding|Gene:Do1_04_a00001G01662:protein_coding;Upstream:1228|Transcript:Do1_04_a00001G01663V1.1:protein_coding|Gene:Do1_04_a00001G01663:protein_coding
Do1_04_a00001G02249V1.1: Sweden_CHH_W_C.bed.out:Do1_04_a00001   12747905        12748144        line_1135;Downstream:949|Transcript:Do1_04_a00001G02250V1.1:protein_coding|Gene:Do1_04_a00001G02250:protein_coding;Upstream:1863|Transcript:Do1_04_a00001G02249V1.1:protein_coding|Gene:Do1_04_a00001G02249:protein_coding;Intergenic:Do1_04_a00001G02249-Do1_04_a00001G02250;Upstream:1612|Transcript:Do1_04_a00001G02251V1.1:protein_coding|Gene:Do1_04_a00001G02251:protein_coding
Do1_05_a00001G00493V1.1: Sweden_W_C.bed.out:Do1_05_a00001       2371094 2371354 line_991;Exon:1:1:RETAINED|Transcript:Do1_05_a00001G00492V1.1:protein_coding|Gene:Do1_05_a00001G00492:protein_coding;Upstream:1873|Transcript:Do1_05_a00001G00493V1.1:protein_coding|Gene:Do1_05_a00001G00493:protein_coding;Upstream:1137|Transcript:Do1_05_a00001G00491V1.1:protein_coding|Gene:Do1_05_a00001G00491:protein_coding;Upstream:101|Transcript:Do1_05_a00001G00492V1.1:protein_coding|Gene:Do1_05_a00001G00492:protein_coding;Intergenic:Do1_05_a00001G00491-Do1_05_a00001G00492
Do1_05_a00003G00085V1.1: Svalbard_CHH_W_C.bed.out:Do1_05_a00003 429829  429883  line_958;Upstream:3930|Transcript:Do1_05_a00003G00086V1.1:protein_coding|Gene:Do1_05_a00003G00086:protein_coding;Intergenic:Do1_05_a00003G00085-Do1_05_a00003G00086;Upstream:2855|Transcript:Do1_05_a00003G00085V1.1:protein_coding|Gene:Do1_05_a00003G00085:protein_coding
Do1_05_a00003G00085V1.1: Sweden_CHH_W_C.bed.out:Do1_05_a00003   421747  421972  line_1459;Intergenic:Do1_05_a00003G00084-Do1_05_a00003G00085;Upstream:1944|Transcript:Do1_05_a00003G00084V1.1:protein_coding|Gene:Do1_05_a00003G00084:protein_coding;Downstream:4995|Transcript:Do1_05_a00003G00085V1.1:protein_coding|Gene:Do1_05_a00003G00085:protein_coding
Do1_06_a00001G02298V1.1: SE_W_C.bed.out:Do1_06_a00001   12932542        12932826        line_666;Intergenic:Do1_06_a00001G02298-Do1_06_a00001G02299;Downstream:3839|Transcript:Do1_06_a00001G02300V1.1:protein_coding|Gene:Do1_06_a00001G02300:protein_coding;Downstream:2969|Transcript:Do1_06_a00001G02299V1.1:protein_coding|Gene:Do1_06_a00001G02299:protein_coding;Downstream:0|Transcript:Do1_06_a00001G02298V1.1:protein_coding|Gene:Do1_06_a00001G02298:protein_coding;Upstream:1311|Transcript:Do1_06_a00001G02297V1.1:protein_coding|Gene:Do1_06_a00001G02297:protein_coding;Exon:2:2:RETAINED|Transcript:Do1_06_a00001G02298V1.1:protein_coding|Gene:Do1_06_a00001G02298:protein_coding

# Search GO terms

for value in $(awk 'NR > 1 {print $1}' Norway_DEGs_overlap_DMRs.txt | sort -u ); do
    echo "Searching for $value in GO term file..."
    grep -H "$value"  interproscan_dryas_full3.tsv >> GO_terms_Norway_DEGs_overlap_DMRs.txt
done

#-----------------------
# total

awk '{print $1}' RNA_DERs_Nov2024_W_C_Total.txt |wc -l
# 3

for value in $(awk 'NR > 1 {print $1}' RNA_DERs_Nov2024_W_C_Total.txt ); do
    echo "Searching for $value in *.out files..."
    # Search for the Norway_DEGs_overlap_DMRslue in all *.out files
    grep -H "$value" *.out | awk -v term="$value" '{print term ": " $0}' >> Total_DEGs_overlap_DMRs.txt
done
wc -l Total_DEGs_overlap_DMRs.txt
0


######################
# Check GO terms

cat GO_terms_* > overall_GO_terms.txt
awk '{print $NF}' overall_GO_terms.txt | sed 's/"//g'

GO:0006281
GO:0020037
GO:0008234
GO:0016301
GO:0006629
GO:0005515
GO:0015074
GO:0005515
GO:0050660
GO:0006508
GO:0042744
GO:0016491
GO:0006364
GO:0008757
GO:0007166
GO:0030247
GO:0016788
GO:0051087

	GO:0006281	DNA repair		NaN	6.065	2.861%	0.460	0.646
	GO:0006364	rRNA processing		NaN	5.790	1.518%	0.622	0.107
	GO:0006508	proteolysis		NaN	6.335	5.329%	0.752	0.206
	GO:0006629	lipid metabolic process		NaN	6.386	5.980%	0.841	0.050
	GO:0007166	cell surface receptor signaling pathway		NaN	5.850	1.743%	0.846	0.014
	GO:0015074	DNA integration		NaN	5.485	0.752%	0.600	0.384
	GO:0042744	hydrogen peroxide catabolic process		NaN	4.887	0.190%	0.900	0.000
	
	GO:0005515	protein binding		NaN	6.579	8.271%	0.948	0.066
	GO:0008234	cysteine-type peptidase activity		NaN	5.444	0.606%	0.904	0.000
	GO:0008757	S-adenosylmethionine-dependent methyltransferase activity		NaN	5.788	1.340%	0.889	0.042
	GO:0016301	kinase activity		NaN	6.426	5.816%	0.881	0.354
	GO:0016491	oxidoreductase activity		NaN	6.705	11.075%	0.940	0.061
	GO:0016788	hydrolase activity, acting on ester bonds		NaN	6.425	5.810%	0.892	0.296
	GO:0020037	heme binding		NaN	5.829	1.474%	0.937	0.044
	GO:0030247	polysaccharide binding		NaN	4.837	0.150%	0.963	0.000
	GO:0050660	flavin adenine dinucleotide binding		NaN	5.877	1.646%	0.936	0.138
	GO:0051087	protein-folding chaperone binding		NaN	5.097	0.273%	0.961	0.036

#########################
# Get counts of increased/decreased expression and Hyper/hypo methylation

# CpG DMRs
cd /home/celphin/scratch/Dryas/CrossMap/file_for_converting/BedGraphs_From_Metilene/

for file in $(ls *.bedGraph); do
echo ${file}
grep "-" ${file} | wc -l 
grep -v "-" ${file} | wc -l 
wc -l ${file}
done

Alaska_W_C.bedGraph
471
395
866 Alaska_W_C.bedGraph
Mat_Sen.bedGraph
26
4963
4989 Mat_Sen.bedGraph
Nunavut_W_C.bedGraph
169
153
322 Nunavut_W_C.bedGraph
Parent_W_C.bedGraph
181
268
449 Parent_W_C.bedGraph
SE_L_H.bedGraph
277
149
426 SE_L_H.bedGraph
SE_W_C.bedGraph
435
453
888 SE_W_C.bedGraph
Svalbard_W_C.bedGraph
158
200
358 Svalbard_W_C.bedGraph
Sweden_W_C.bedGraph
958
779
1737 Sweden_W_C.bedGraph
Wild_Lat_L_H.bedGraph
6992
7994
14986 Wild_Lat_L_H.bedGraph
Wild_W_C.bedGraph
203
182
385 Wild_W_C.bedGraph

#-----------------
cd /home/celphin/scratch/Dryas/CrossMap/file_for_converting/Bedgraphs_Intersected_Subtracted

for file in $(ls *.bedGraph); do
echo ${file}
grep "-" ${file} | wc -l 
grep -v "-" ${file} | wc -l 
wc -l ${file}
done

intersect_SE_L_H_Wild_L_H.bedGraph
215
123
338 intersect_SE_L_H_Wild_L_H.bedGraph
intersect_SE_W_C_P_W_C.bedGraph
70
82
152 intersect_SE_W_C_P_W_C.bedGraph
intersect_SE_W_C_SE_L_H.bedGraph
50
58
108 intersect_SE_W_C_SE_L_H.bedGraph
intersect_SE_W_C_Wild_W_C.bedGraph
44
40
84 intersect_SE_W_C_Wild_W_C.bedGraph
intersect_Wild_W_C_Mat_Sen.bedGraph
29
40
69 intersect_Wild_W_C_Mat_Sen.bedGraph
total_subtract_SE_W_C_P_W_C.bedGraph
365
371
736 total_subtract_SE_W_C_P_W_C.bedGraph
total_subtract_SE_W_C_SE_L_H.bedGraph
385
395
780 total_subtract_SE_W_C_SE_L_H.bedGraph
total_subtract_W_C_Mat_Sen.bedGraph
174
142
316 total_subtract_W_C_Mat_Sen.bedGraph

for file in $(ls *.bedgraph); do
echo ${file}
grep "-" ${file} | wc -l 
grep -v "-" ${file} | wc -l 
wc -l ${file}
done

#-------------------
# CHH DMRs
cd /home/celphin/scratch/Dryas/CrossMap/file_for_converting/CHH

for file in $(ls *.bedGraph); do
echo ${file}
grep "-" ${file} | wc -l 
grep -v "-" ${file} | wc -l 
wc -l ${file}
done

Alaska_CHH_W_C_CHH.bedGraph
436
1501
1937 Alaska_CHH_W_C_CHH.bedGraph
ALEX_SVAL_ALAS_Sites_intersect_DMRs_CHH.bedGraph
8
12
20 ALEX_SVAL_ALAS_Sites_intersect_DMRs_CHH.bedGraph
ALEX_SVAL_SWED_Sites_intersect_DMRs_CHH.bedGraph
13
13
26 ALEX_SVAL_SWED_Sites_intersect_DMRs_CHH.bedGraph
ALEX_SWED_ALAS_Sites_intersect_DMRs_CHH.bedGraph
5
14
19 ALEX_SWED_ALAS_Sites_intersect_DMRs_CHH.bedGraph
ALL_Sites_intersect_DMRs_CHH.bedGraph
2
4
6 ALL_Sites_intersect_DMRs_CHH.bedGraph
Nunavut_CHH_W_C_CHH.bedGraph
320
866
1186 Nunavut_CHH_W_C_CHH.bedGraph
Svalbard_CHH_W_C_CHH.bedGraph
725
672
1397 Svalbard_CHH_W_C_CHH.bedGraph
SVAL_SWED_ALAS_Sites_intersect_DMRs_CHH.bedGraph
22
15
37 SVAL_SWED_ALAS_Sites_intersect_DMRs_CHH.bedGraph
Sweden_CHH_W_C_CHH.bedGraph
1138
1068
2206 Sweden_CHH_W_C_CHH.bedGraph


################################
# DEGs
# need to look in r

cd /home/celphin/scratch/Dryas/RNAseq_analysis

# counts
grep -v "-" RNA_Alaska_W_C_DERs_updown.txt | wc -l 
grep -v "-" RNA_Norway_W_C_DERs_updown.txt | wc -l 
grep -v "-" RNA_Sweden_W_C_DERs_updown.txt | wc -l 
grep -v "-" RNA_Alex_W_C_DERs_updown.txt | wc -l 

grep  "-" RNA_Alaska_W_C_DERs_updown.txt | wc -l 
grep  "-" RNA_Norway_W_C_DERs_updown.txt | wc -l 
grep  "-" RNA_Sweden_W_C_DERs_updown.txt | wc -l 
grep  "-" RNA_Alex_W_C_DERs_updown.txt | wc -l 

# make FILES
cd /home/celphin/scratch/Dryas/RNAseq_analysis/snpEff

grep -v "-" ../RNA_Alaska_W_C_DERs_updown.txt |awk '{print $2}' > RNA_Alaska_W_C_DERs_up
grep -v "-" ../RNA_Norway_W_C_DERs_updown.txt |awk '{print $2}' > RNA_Norway_W_C_DERs_up
grep -v "-" ../RNA_Sweden_W_C_DERs_updown.txt |awk '{print $2}' > RNA_Sweden_W_C_DERs_up
grep -v "-" ../RNA_Alex_W_C_DERs_updown.txt |awk '{print $2}' > RNA_Alex_W_C_DERs_up

grep "-" ../RNA_Alaska_W_C_DERs_updown.txt |awk '{print $2}' > RNA_Alaska_W_C_DERs_down
grep "-" ../RNA_Norway_W_C_DERs_updown.txt |awk '{print $2}' > RNA_Norway_W_C_DERs_down
grep "-" ../RNA_Sweden_W_C_DERs_updown.txt |awk '{print $2}' > RNA_Sweden_W_C_DERs_down
grep "-" ../RNA_Alex_W_C_DERs_updown.txt |awk '{print $2}' > RNA_Alex_W_C_DERs_down

# remove the V1.1

sed -i 's/V1\.1//g' RNA*_W_C_DERs_*

#######################################
# compare to DMRs for matching genes
cd /home/celphin/scratch/Dryas/snpEff/

cp /home/celphin/scratch/Dryas/snpEff/out/*.out /home/celphin/scratch/Dryas/RNAseq_analysis/snpEff
cp /home/celphin/scratch/Dryas/snpEff/Dryas_DMRs/*.bedGraph /home/celphin/scratch/Dryas/RNAseq_analysis/snpEff

cp /home/celphin/scratch/Dryas/snpEff/Dryas_DMRs/CHH/*.out /home/celphin/scratch/Dryas/RNAseq_analysis/snpEff/CHH
cp /home/celphin/scratch/Dryas/snpEff/Dryas_DMRs/CHH/*.bedGraph /home/celphin/scratch/Dryas/RNAseq_analysis/snpEff/CHH

#--------------------------
cd /home/celphin/scratch/Dryas/RNAseq_analysis/snpEff

# remove the header lines of the snpEff file
for file in *.out; do
  sed -i.bak '/^#/d' "$file" 
done

# join with hyper and hypo info
paste Alaska_W_C.bedGraph Alaska_W_C.bed.out > Alaska_W_C_combined.out
paste Nunavut_W_C.bedGraph Nunavut_W_C.bed.out > Nunavut_W_C_combined.out
paste Parent_W_C.bedGraph Parent_W_C.bed.out > Parent_W_C_combined.out
paste SE_W_C.bedGraph SE_W_C.bed.out > SE_W_C_combined.out
paste Svalbard_W_C.bedGraph Svalbard_W_C.bed.out > Svalbard_W_C_combined.out
paste Sweden_W_C.bedGraph Sweden_W_C.bed.out > Sweden_W_C_combined.out

# split DMRs by hyper and hyp
grep -v "-" Alaska_W_C_combined.out  > Alaska_W_C_up.out
grep "-" Alaska_W_C_combined.out  > Alaska_W_C_down.out

grep -v "-" Nunavut_W_C_combined.out  > Nunavut_W_C_up.out
grep "-" Nunavut_W_C_combined.out > Nunavut_W_C_down.out

grep -v "-" Sweden_W_C_combined.out > Sweden_W_C_up.out
grep "-" Sweden_W_C_combined.out  > Sweden_W_C_down.out

grep -v "-" Svalbard_W_C_combined.out  > Svalbard_W_C_up.out
grep "-" Svalbard_W_C_combined.out  > Svalbard_W_C_down.out

grep -v "-" SE_W_C_combined.out  > SE_W_C_up.out
grep "-" SE_W_C_combined.out  > SE_W_C_down.out

grep -v "-" Parent_W_C_combined.out  > Parent_W_C_up.out
grep "-" Parent_W_C_combined.out  > Parent_W_C_down.out

#----------------------------------
# compare hyper and hypo methylation to up and down regulated genes

# Loop through all .up files
for gene_file in *up; do
    if [[ -f "$gene_file" ]]; then
        # Output file for results related to .up genes
        output_up_hyper="up_hyper_results.txt"
        output_up_hypo="up_hypo_results.txt"

        # Read each gene from the current .up file
        while read -r gene; do
            # Search for the gene in all .out files and append the gene name as a column to the results
            find . -type f -name "*_up.out" -exec grep -H "$gene" {} \; | awk -v gene="$gene" '{print $0 "\tGene: " gene }' | awk -v gene_file="$gene_file"  '{print $0 "\tFile: " gene_file }' >> "$output_up_hyper"
        done < "$gene_file"
        # Read each gene from the current .up file
        while read -r gene; do
            # Search for the gene in all .out files and append the gene name as a column to the results
            find . -type f -name "*_down.out" -exec grep -H "$gene" {} \; | awk -v gene="$gene"  '{print $0 "\tGene: " gene }' | awk -v gene_file="$gene_file"  '{print $0 "\tFile: " gene_file }' >> "$output_up_hypo"
        done < "$gene_file"
    fi
done

for gene_file in *down; do
    if [[ -f "$gene_file" ]]; then
        # Output file for results related to .up genes
        output_down_hyper="down_hyper_results.txt"
        output_down_hypo="down_hypo_results.txt"

        # Read each gene from the current .up file
        while read -r gene; do
            # Search for the gene in all .out files and append the gene name as a column to the results
            find . -type f -name "*_up.out" -exec grep -H "$gene" {} \; | awk -v gene="$gene" '{print $0 "\tGene: " gene }' | awk -v gene_file="$gene_file"  '{print $0 "\tFile: " gene_file }'  >> "$output_down_hyper"
        done < "$gene_file"
        # Read each gene from the current .up file
        while read -r gene; do
            # Search for the gene in all .out files and append the gene name as a column to the results
            find . -type f -name "*_down.out" -exec grep -H "$gene" {} \; | awk -v gene="$gene"  '{print $0 "\tGene: " gene }' | awk -v gene_file="$gene_file"  '{print $0 "\tFile: " gene_file }' >> "$output_down_hypo"
        done < "$gene_file"
    fi
done

wc -l up_hyper_results.txt
wc -l up_hypo_results.txt
wc -l down_hyper_results.txt
wc -l down_hypo_results.txt


8 up_hyper_results.txt
41 up_hypo_results.txt

6 down_hyper_results.txt
24 down_hypo_results.txt


#-------------------------
# Extract the gene name and types

FILES=$(ls *.txt)
echo ${FILES}

for file in ${FILES}; do
awk -F'[;:]' '{print $3 $7}' "$file" > "$file"1
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


down_hyper_results.txt 
Upstream 5
Exon 1

up_hyper_results.txt 
Upstream 2
Downstream 1
Gene 1
Exon 4

down_hypo_results.txt 
Upstream 7
Downstream 7
Intergenic 6
Exon 4

up_hypo_results.txt 
Upstream 14
Downstream 8
Gene 3
Intergenic 9
Intron 1
Exon 6

##########################
# CHH

cd /home/celphin/scratch/Dryas/RNAseq_analysis/snpEff/CHH

# remove the header lines of the snpEff file
for file in *.out; do
  sed -i.bak '/^#/d' "$file" 
done

# join with hyper and hypo info
paste Alaska_CHH_W_C.bedGraph Alaska_CHH_W_C.bed.out > Alaska_CHH_W_C_combined.out
paste Nunavut_CHH_W_C.bedGraph Nunavut_CHH_W_C.bed.out > Nunavut_CHH_W_C_combined.out
paste Svalbard_CHH_W_C.bedGraph Svalbard_CHH_W_C.bed.out > Svalbard_CHH_W_C_combined.out
paste Sweden_CHH_W_C.bedGraph Sweden_CHH_W_C.bed.out > Sweden_CHH_W_C_combined.out

# split DMRs by hyper and hyp
grep -v "-" Alaska_CHH_W_C_combined.out  > Alaska_CHH_W_C_up.out
grep "-" Alaska_CHH_W_C_combined.out  > Alaska_CHH_W_C_down.out

grep -v "-" Nunavut_CHH_W_C_combined.out  > Nunavut_CHH_W_C_up.out
grep "-" Nunavut_CHH_W_C_combined.out > Nunavut_CHH_W_C_down.out

grep -v "-" Sweden_CHH_W_C_combined.out > Sweden_CHH_W_C_up.out
grep "-" Sweden_CHH_W_C_combined.out  > Sweden_CHH_W_C_down.out

grep -v "-" Svalbard_CHH_W_C_combined.out  > Svalbard_CHH_W_C_up.out
grep "-" Svalbard_CHH_W_C_combined.out  > Svalbard_CHH_W_C_down.out

cp ../RNA_* .

# compare hyper and hypo methylation to up and down regulated genes

# Loop through all .up files
for gene_file in *up; do
    if [[ -f "$gene_file" ]]; then
        # Output file for results related to .up genes
        output_up_hyper="up_hyper_results.txt"
        output_up_hypo="up_hypo_results.txt"

        # Read each gene from the current .up file
        while read -r gene; do
            # Search for the gene in all .out files and append the gene name as a column to the results
            find . -type f -name "*_up.out" -exec grep -H "$gene" {} \; | awk -v gene="$gene" '{print $0 "\tGene: " gene }' | awk -v gene_file="$gene_file"  '{print $0 "\tFile: " gene_file }' >> "$output_up_hyper"
        done < "$gene_file"
        # Read each gene from the current .up file
        while read -r gene; do
            # Search for the gene in all .out files and append the gene name as a column to the results
            find . -type f -name "*_down.out" -exec grep -H "$gene" {} \; | awk -v gene="$gene"  '{print $0 "\tGene: " gene }' | awk -v gene_file="$gene_file"  '{print $0 "\tFile: " gene_file }' >> "$output_up_hypo"
        done < "$gene_file"
    fi
done

for gene_file in *down; do
    if [[ -f "$gene_file" ]]; then
        # Output file for results related to .up genes
        output_down_hyper="down_hyper_results.txt"
        output_down_hypo="down_hypo_results.txt"

        # Read each gene from the current .up file
        while read -r gene; do
            # Search for the gene in all .out files and append the gene name as a column to the results
            find . -type f -name "*_up.out" -exec grep -H "$gene" {} \; | awk -v gene="$gene" '{print $0 "\tGene: " gene }' | awk -v gene_file="$gene_file"  '{print $0 "\tFile: " gene_file }'  >> "$output_down_hyper"
        done < "$gene_file"
        # Read each gene from the current .up file
        while read -r gene; do
            # Search for the gene in all .out files and append the gene name as a column to the results
            find . -type f -name "*_down.out" -exec grep -H "$gene" {} \; | awk -v gene="$gene"  '{print $0 "\tGene: " gene }' | awk -v gene_file="$gene_file"  '{print $0 "\tFile: " gene_file }' >> "$output_down_hypo"
        done < "$gene_file"
    fi
done

wc -l up_hyper_results.txt
wc -l up_hypo_results.txt
wc -l down_hyper_results.txt
wc -l down_hypo_results.txt


1 up_hyper_results.txt

45 up_hypo_results.txt

0 down_hyper_results.txt

57 down_hypo_results.txt


#-------------------------
# Extract the gene name and types

FILES=$(ls *.txt)
echo ${FILES}

for file in ${FILES}; do
awk -F'[;:]' '{print $3 $7}' "$file" > "$file"1
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

down_hypo_results.txt Upstream
18
down_hypo_results.txt Downstream
17
down_hypo_results.txt Gene
0
down_hypo_results.txt Intergenic
18
down_hypo_results.txt Intron
1
down_hypo_results.txt Exon
3

up_hyper_results.txt Upstream
1
up_hyper_results.txt Downstream
0
up_hyper_results.txt Gene
0
up_hyper_results.txt Intergenic
0
up_hyper_results.txt Intron
0
up_hyper_results.txt Exon
0

up_hypo_results.txt Upstream
17
up_hypo_results.txt Downstream
13
up_hypo_results.txt Gene
0
up_hypo_results.txt Intergenic
15
up_hypo_results.txt Intron
0
up_hypo_results.txt Exon
0

#---------------------------
# check if percent on TEs is higher for hyper or hypo meth

module load StdEnv/2023 bedtools/2.31.0

cd /home/celphin/scratch/Dryas/RNAseq_analysis/snpEff/CHH

for file in *.bedGraph; do
    # Extract rows with negative values in the fourth column and save to a file
    awk -F'\t' '$4 < 0 {print $1, $2, $3}' "$file" > "${file%.bedGraph}_hypo.bed"
    
    # Extract rows with positive values in the fourth column and save to a file
    awk -F'\t' '$4 > 0 {print $1, $2, $3}' "$file" > "${file%.bedGraph}_hyper.bed"
done

for file in *.bed; do
    # Replace spaces with tabs in the file and overwrite the file
    sed -i 's/ /\t/g' "$file"
done

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

1154 Alaska_CHH_W_C_hyper.bed-All-TE.bed
289 Alaska_CHH_W_C_hypo.bed-All-TE.bed
708 Nunavut_CHH_W_C_hyper.bed-All-TE.bed
263 Nunavut_CHH_W_C_hypo.bed-All-TE.bed
532 Svalbard_CHH_W_C_hyper.bed-All-TE.bed
575 Svalbard_CHH_W_C_hypo.bed-All-TE.bed
839 Sweden_CHH_W_C_hyper.bed-All-TE.bed
781 Sweden_CHH_W_C_hypo.bed-All-TE.bed

#---------------------
cd /home/celphin/scratch/Dryas/RNAseq_analysis/snpEff/

for file in *.bedGraph; do
    # Extract rows with negative values in the fourth column and save to a file
    awk -F'\t' '$4 < 0 {print $1, $2, $3}' "$file" > "${file%.bedGraph}_hypo.bed"
    
    # Extract rows with positive values in the fourth column and save to a file
    awk -F'\t' '$4 > 0 {print $1, $2, $3}' "$file" > "${file%.bedGraph}_hyper.bed"
done

for file in *.bed; do
    # Replace spaces with tabs in the file and overwrite the file
    sed -i 's/ /\t/g' "$file"
done

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

153 Alaska_W_C_hyper.bed-All-TE.bed
177 Alaska_W_C_hypo.bed-All-TE.bed
65 Nunavut_W_C_hyper.bed-All-TE.bed
74 Nunavut_W_C_hypo.bed-All-TE.bed
105 Parent_W_C_hyper.bed-All-TE.bed
71 Parent_W_C_hypo.bed-All-TE.bed
89 SE_W_C_hyper.bed-All-TE.bed
136 SE_W_C_hypo.bed-All-TE.bed
72 Svalbard_W_C_hyper.bed-All-TE.bed
64 Svalbard_W_C_hypo.bed-All-TE.bed
292 Sweden_W_C_hyper.bed-All-TE.bed
289 Sweden_W_C_hypo.bed-All-TE.bed
78 Wild_W_C_hyper.bed-All-TE.bed
95 Wild_W_C_hypo.bed-All-TE.bed

#---------------------
# methylkit
cd /home/celphin/scratch/Dryas/CpG/stranded_CpG_report/DMRs

for file in *.bedGraph; do
    # Extract rows with negative values in the fourth column and save to a file
    awk -F'\t' '$4 < 0 {print $1, $2, $3}' "$file" > "${file%.bedGraph}_hypo.bed"
    
    # Extract rows with positive values in the fourth column and save to a file
    awk -F'\t' '$4 > 0 {print $1, $2, $3}' "$file" > "${file%.bedGraph}_hyper.bed"
done

for file in *.bed; do
    # Replace spaces with tabs in the file and overwrite the file
    sed -i 's/ /\t/g' "$file"
done

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

2 intersect_ALAS_LAT_W_C_25plowcovDMRs_hyper.bed-All-TE.bed
0 intersect_ALAS_LAT_W_C_25plowcovDMRs_hypo.bed-All-TE.bed
1 intersect_Alex_LAT_SVAL_plowcovDMRs_hyper.bed-All-TE.bed
0 intersect_Alex_LAT_SVAL_plowcovDMRs_hypo.bed-All-TE.bed
4 intersect_Alex_LAT_W_C_25plowcovDMRs_hyper.bed-All-TE.bed
4 intersect_Alex_LAT_W_C_25plowcovDMRs_hypo.bed-All-TE.bed
0 intersect_ALL_plowcovDMRs_hyper.bed-All-TE.bed
0 intersect_ALL_plowcovDMRs_hypo.bed-All-TE.bed
7 intersect_SVAL_Alex_W_C_25plowcovDMRs_hyper.bed-All-TE.bed
4 intersect_SVAL_Alex_W_C_25plowcovDMRs_hypo.bed-All-TE.bed
6 intersect_SVAL_LAT_W_C_25plowcovDMRs_hyper.bed-All-TE.bed
1 intersect_SVAL_LAT_W_C_25plowcovDMRs_hypo.bed-All-TE.bed
32 Methylkit_ALAS_W_C_25plowcovDMRs_hyper.bed-All-TE.bed
37 Methylkit_ALAS_W_C_25plowcovDMRs_hypo.bed-All-TE.bed
1247 Methylkit_ALAS_W_C_lowcovDMRs_hyper.bed-All-TE.bed
1133 Methylkit_ALAS_W_C_lowcovDMRs_hypo.bed-All-TE.bed
203 Methylkit_Alex_W_C_10plowcovDMRs_hyper.bed-All-TE.bed
154 Methylkit_Alex_W_C_10plowcovDMRs_hypo.bed-All-TE.bed
112 Methylkit_LAT_W_C_25plowcovDMRs_hyper.bed-All-TE.bed
83 Methylkit_LAT_W_C_25plowcovDMRs_hypo.bed-All-TE.bed
77 Methylkit_SVAL_W_C_25plowcovDMRs_hyper.bed-All-TE.bed
59 Methylkit_SVAL_W_C_25plowcovDMRs_hypo.bed-All-TE.bed
10 Methylkit_Wild_W_C_DMRs_hyper.bed-All-TE.bed
7 Methylkit_Wild_W_C_DMRs_hypo.bed-All-TE.bed
26 Methylkit_Wild_W_C_lowcovDMRs_hyper.bed-All-TE.bed
26 Methylkit_Wild_W_C_lowcovDMRs_hypo.bed-All-TE.bed


wc -l *hyper.bed
     3 intersect_ALAS_LAT_W_C_25plowcovDMRs_hyper.bed
     2 intersect_Alex_LAT_SVAL_plowcovDMRs_hyper.bed
     9 intersect_Alex_LAT_W_C_25plowcovDMRs_hyper.bed
     9 intersect_SVAL_Alex_W_C_25plowcovDMRs_hyper.bed
     7 intersect_SVAL_LAT_W_C_25plowcovDMRs_hyper.bed
    67 Methylkit_ALAS_W_C_25plowcovDMRs_hyper.bed
  2476 Methylkit_ALAS_W_C_lowcovDMRs_hyper.bed
   376 Methylkit_Alex_W_C_10plowcovDMRs_hyper.bed
   183 Methylkit_LAT_W_C_25plowcovDMRs_hyper.bed
   140 Methylkit_SVAL_W_C_25plowcovDMRs_hyper.bed
    21 Methylkit_Wild_W_C_DMRs_hyper.bed
   116 Methylkit_Wild_W_C_lowcovDMRs_hyper.bed
   
wc -l *hypo.bed
     1 intersect_ALAS_LAT_W_C_25plowcovDMRs_hypo.bed
     0 intersect_Alex_LAT_SVAL_plowcovDMRs_hypo.bed
     4 intersect_Alex_LAT_W_C_25plowcovDMRs_hypo.bed
     4 intersect_SVAL_Alex_W_C_25plowcovDMRs_hypo.bed
     2 intersect_SVAL_LAT_W_C_25plowcovDMRs_hypo.bed
    72 Methylkit_ALAS_W_C_25plowcovDMRs_hypo.bed
  2610 Methylkit_ALAS_W_C_lowcovDMRs_hypo.bed
   302 Methylkit_Alex_W_C_10plowcovDMRs_hypo.bed
   145 Methylkit_LAT_W_C_25plowcovDMRs_hypo.bed
   115 Methylkit_SVAL_W_C_25plowcovDMRs_hypo.bed
    11 Methylkit_Wild_W_C_DMRs_hypo.bed
    96 Methylkit_Wild_W_C_lowcovDMRs_hypo.bed

# hypermethylation of TEs with warming chambers is slightly higher

###############################
# Sort overlap files by first column

sort -k1,1 up_hypo_results.txt > up_hypo_results_sorted.txt

# really driven by overlap in Sweden

more up_hypo_results_sorted.txt
./Alaska_W_C_down.out:Do1_01_a00001     10414157        10414286        24.791812       Do1_01_a00001   10
414157  10414286        line_59;Downstream:3137|Transcript:Do1_01_a00001G01909V1.1:protein_coding|Gene:Do1
_01_a00001G01909:protein_coding;Intergenic:Do1_01_a00001G01908-Do1_01_a00001G01909;Downstream:2141|Transcr
ipt:Do1_01_a00001G01908V1.1:protein_coding|Gene:Do1_01_a00001G01908:protein_coding      Gene: Do1_01_a0000
1G01909 File: RNA_Sweden_W_C_DERs_up
./Alaska_W_C_down.out:Do1_01_a00001     3955238 3955782 -6.108593       Do1_01_a00001   3955238 3955782 li
ne_21;Gene:Do1_01_a00001G00761:protein_coding;Upstream:1141|Transcript:Do1_01_a00001G00760V1.1:protein_cod
ing|Gene:Do1_01_a00001G00760:protein_coding;Downstream:3187|Transcript:Do1_01_a00001G00762V1.1:protein_cod
ing|Gene:Do1_01_a00001G00762:protein_coding     Gene: Do1_01_a00001G00762       File: RNA_Alaska_W_C_DERs_
up
./Alaska_W_C_down.out:Do1_01_a00007     793375  793496  47.069723       Do1_01_a00007   793375  793496  li
ne_167;Upstream:3476|Transcript:Do1_01_a00007G00167V1.1:protein_coding|Gene:Do1_01_a00007G00167:protein_co
ding;Upstream:4305|Transcript:Do1_01_a00007G00164V1.1:protein_coding|Gene:Do1_01_a00007G00164:protein_codi
ng;Downstream:3800|Transcript:Do1_01_a00007G00165V1.1:protein_coding|Gene:Do1_01_a00007G00165:protein_codi
ng;Intergenic:Do1_01_a00007G00166-Do1_01_a00007G00167;Downstream:2546|Transcript:Do1_01_a00007G00166V1.1:p
rotein_coding|Gene:Do1_01_a00007G00166:protein_coding   Gene: Do1_01_a00007G00164       File: RNA_Sweden_W
_C_DERs_up
./Alaska_W_C_down.out:Do1_03_a00002     2787693 2787975 -16.751389      Do1_03_a00002   2787693 2787975 li
ne_324;Downstream:0|Transcript:Do1_03_a00002G00495V1.1:protein_coding|Gene:Do1_03_a00002G00495:protein_cod
ing;Upstream:3201|Transcript:Do1_03_a00002G00494V1.1:protein_coding|Gene:Do1_03_a00002G00494:protein_codin
g;Intergenic:Do1_03_a00002G00494-Do1_03_a00002G00495;Exon:1:1:RETAINED|Transcript:Do1_03_a00002G00495V1.1:
protein_coding|Gene:Do1_03_a00002G00495:protein_coding;Downstream:941|Transcript:Do1_03_a00002G00496V1.1:p
rotein_coding|Gene:Do1_03_a00002G00496:protein_coding;Intergenic:Do1_03_a00002G00495-Do1_03_a00002G00496;U
pstream:11|Transcript:Do1_03_a00002G00495V1.1:protein_coding|Gene:Do1_03_a00002G00495:protein_coding;Upstr
eam:3800|Transcript:Do1_03_a00002G00493V1.1:protein_coding|Gene:Do1_03_a00002G00493:protein_coding      Ge
ne: Do1_03_a00002G00493 File: RNA_Alaska_W_C_DERs_up
./Alaska_W_C_down.out:Do1_06_a00001     7563000 7563197 -33.187374      Do1_06_a00001   7563000 7563197 li
ne_646;Intergenic:Do1_06_a00001G01290-Do1_06_a00001G01291;Upstream:2430|Transcript:Do1_06_a00001G01290V1.1
:protein_coding|Gene:Do1_06_a00001G01290:protein_coding;Downstream:3861|Transcript:Do1_06_a00001G01291V1.1
:protein_coding|Gene:Do1_06_a00001G01291:protein_coding Gene: Do1_06_a00001G01291       File: RNA_Sweden_W
_C_DERs_up
./Alaska_W_C_down.out:Do1_06_a00002     4486112 4486383 21.591812       Do1_06_a00002   4486112 4486383 li
ne_687;Downstream:1636|Transcript:Do1_06_a00002G00928V1.1:protein_coding|Gene:Do1_06_a00002G00928:protein_
coding;Upstream:672|Transcript:Do1_06_a00002G00927V1.1:protein_coding|Gene:Do1_06_a00002G00927:protein_cod
ing;Intergenic:Do1_06_a00002G00927-Do1_06_a00002G00928  Gene: Do1_06_a00002G00928       File: RNA_Sweden_W
_C_DERs_up
./Alaska_W_C_down.out:Do1_06_a00002     587971  588152  34.489606       Do1_06_a00002   587971  588152  li
ne_670;Upstream:2445|Transcript:Do1_06_a00002G00118V1.1:protein_coding|Gene:Do1_06_a00002G00118:protein_co
ding;Upstream:1324|Transcript:Do1_06_a00002G00117V1.1:protein_coding|Gene:Do1_06_a00002G00117:protein_codi
ng;Intergenic:Do1_06_a00002G00115-Do1_06_a00002G00116;Downstream:3096|Transcript:Do1_06_a00002G00119V1.1:p
rotein_coding|Gene:Do1_06_a00002G00119:protein_coding;Downstream:2975|Transcript:Do1_06_a00002G00115V1.1:p
rotein_coding|Gene:Do1_06_a00002G00115:protein_coding;Downstream:386|Transcript:Do1_06_a00002G00116V1.1:pr
otein_coding|Gene:Do1_06_a00002G00116:protein_coding    Gene: Do1_06_a00002G00119       File: RNA_Sweden_W
_C_DERs_up
./Alaska_W_C_down.out:Do1_06_a00002     588236  588420  36.740881       Do1_06_a00002   588236  588420  li
ne_671;Exon:1:1:RETAINED|Transcript:Do1_06_a00002G00116V1.1:protein_coding|Gene:Do1_06_a00002G00116:protei
n_coding;Upstream:2180|Transcript:Do1_06_a00002G00118V1.1:protein_coding|Gene:Do1_06_a00002G00118:protein_
coding;Intergenic:Do1_06_a00002G00115-Do1_06_a00002G00116;Downstream:121|Transcript:Do1_06_a00002G00116V1.
1:protein_coding|Gene:Do1_06_a00002G00116:protein_coding;Upstream:1059|Transcript:Do1_06_a00002G00117V1.1:
protein_coding|Gene:Do1_06_a00002G00117:protein_coding;Downstream:3240|Transcript:Do1_06_a00002G00115V1.1:
protein_coding|Gene:Do1_06_a00002G00115:protein_coding;Downstream:2831|Transcript:Do1_06_a00002G00119V1.1:
protein_coding|Gene:Do1_06_a00002G00119:protein_coding  Gene: Do1_06_a00002G00119       File: RNA_Sweden_W
_C_DERs_up
#------------------
./Nunavut_W_C_down.out:Do1_01_a00002    2973631 2973746 -12.314129      Do1_01_a00002   2973631 2973746 li
ne_30;Upstream:135|Transcript:Do1_01_a00002G00497V1.1:protein_coding|Gene:Do1_01_a00002G00497:protein_codi
ng;Gene:Do1_01_a00002G00498:protein_coding;Downstream:891|Transcript:Do1_01_a00002G00496V1.1:protein_codin
g|Gene:Do1_01_a00002G00496:protein_coding;Upstream:3600|Transcript:Do1_01_a00002G00493V1.1:protein_coding|
Gene:Do1_01_a00002G00493:protein_coding;Upstream:5006|Transcript:Do1_01_a00002G00499V1.1:protein_coding|Ge
ne:Do1_01_a00002G00499:protein_coding;Downstream:2534|Transcript:Do1_01_a00002G00494V1.1:protein_coding|Ge
ne:Do1_01_a00002G00494:protein_coding;Downstream:2035|Transcript:Do1_01_a00002G00495V1.1:protein_coding|Ge
ne:Do1_01_a00002G00495:protein_coding   Gene: Do1_01_a00002G00498       File: RNA_Alaska_W_C_DERs_up
./Nunavut_W_C_down.out:Do1_01_a00002    441249  441476  19.829772       Do1_01_a00002   441249  441476  li
ne_25;Intergenic:Do1_01_a00002G00067-Do1_01_a00002G00068        Gene: Do1_01_a00002G00067       File: RNA_
Sweden_W_C_DERs_up
./Nunavut_W_C_down.out:Do1_04_a00005    194843  195113  -11.450868      Do1_04_a00005   194843  195113  li
ne_175;Upstream:3323|Transcript:Do1_04_a00005G00025V1.1:protein_coding|Gene:Do1_04_a00005G00025:protein_co
ding;Intergenic:Do1_04_a00005G00024-Do1_04_a00005G00025 Gene: Do1_04_a00005G00024       File: RNA_Sweden_W
_C_DERs_up

#---------------------
./SE_W_C_down.out:Do1_01_a00001 2735512 2735769 -39.384367      Do1_01_a00001   2735512 2735769 line_13;Up
stream:3859|Transcript:Do1_01_a00001G00556V1.1:protein_coding|Gene:Do1_01_a00001G00556:protein_coding;Inte
rgenic:Do1_01_a00001G00555-Do1_01_a00001G00556;Downstream:1379|Transcript:Do1_01_a00001G00555V1.1:protein_
coding|Gene:Do1_01_a00001G00555:protein_coding  Gene: Do1_01_a00001G00556       File: RNA_Sweden_W_C_DERs_
up
./SE_W_C_down.out:Do1_01_a00007 788623  789090  27.228386       Do1_01_a00007   788623  789090  line_111;I
ntergenic:Do1_01_a00007G00163-Do1_01_a00007G00164;Upstream:1621|Transcript:Do1_01_a00007G00163V1.1:protein
_coding|Gene:Do1_01_a00007G00163:protein_coding;Upstream:0|Transcript:Do1_01_a00007G00164V1.1:protein_codi
ng|Gene:Do1_01_a00007G00164:protein_coding;Intergenic:Do1_01_a00007G00164-Do1_01_a00007G00165;Upstream:125
0|Transcript:Do1_01_a00007G00166V1.1:protein_coding|Gene:Do1_01_a00007G00166:protein_coding;Exon:1:1:RETAI
NED|Transcript:Do1_01_a00007G00164V1.1:protein_coding|Gene:Do1_01_a00007G00164:protein_coding;Downstream:1
50|Transcript:Do1_01_a00007G00164V1.1:protein_coding|Gene:Do1_01_a00007G00164:protein_coding;Upstream:685|
Transcript:Do1_01_a00007G00165V1.1:protein_coding|Gene:Do1_01_a00007G00165:protein_coding       Gene: Do1_
01_a00007G00164 File: RNA_Sweden_W_C_DERs_up
./SE_W_C_down.out:Do1_02_a00003 6756922 6757188 -50.739069      Do1_02_a00003   6756922 6757188 line_192;U
pstream:4084|Transcript:Do1_02_a00003G01352V1.1:protein_coding|Gene:Do1_02_a00003G01352:protein_coding;Dow
nstream:1076|Transcript:Do1_02_a00003G01350V1.1:protein_coding|Gene:Do1_02_a00003G01350:protein_coding;Exo
n:1:1:RETAINED|Transcript:Do1_02_a00003G01351V1.1:protein_coding|Gene:Do1_02_a00003G01351:protein_coding -Gene: Do1_02_a00003G01351       File: RNA_Alaska_W_C_DERs_up
./SE_W_C_down.out:Do1_03_a00002 2787713 2787943 -28.135316      Do1_03_a00002   2787713 2787943 line_288;D
ownstream:0|Transcript:Do1_03_a00002G00495V1.1:protein_coding|Gene:Do1_03_a00002G00495:protein_coding;Down
stream:921|Transcript:Do1_03_a00002G00496V1.1:protein_coding|Gene:Do1_03_a00002G00496:protein_coding;Exon:
1:1:RETAINED|Transcript:Do1_03_a00002G00495V1.1:protein_coding|Gene:Do1_03_a00002G00495:protein_coding;Ups
tream:3820|Transcript:Do1_03_a00002G00493V1.1:protein_coding|Gene:Do1_03_a00002G00493:protein_coding;Inter
genic:Do1_03_a00002G00495-Do1_03_a00002G00496;Upstream:3221|Transcript:Do1_03_a00002G00494V1.1:protein_cod
ing|Gene:Do1_03_a00002G00494:protein_coding     Gene: Do1_03_a00002G00493       File: RNA_Alaska_W_C_DERs_
up
./SE_W_C_down.out:Do1_06_a00001 12932542        12932826        26.974964       Do1_06_a00001   12932542 -12932826        line_666;Intergenic:Do1_06_a00001G02298-Do1_06_a00001G02299;Downstream:3839|Transcript:Do1
_06_a00001G02300V1.1:protein_coding|Gene:Do1_06_a00001G02300:protein_coding;Downstream:2969|Transcript:Do1
_06_a00001G02299V1.1:protein_coding|Gene:Do1_06_a00001G02299:protein_coding;Downstream:0|Transcript:Do1_06
_a00001G02298V1.1:protein_coding|Gene:Do1_06_a00001G02298:protein_coding;Upstream:1311|Transcript:Do1_06_a
00001G02297V1.1:protein_coding|Gene:Do1_06_a00001G02297:protein_coding;Exon:2:2:RETAINED|Transcript:Do1_06
_a00001G02298V1.1:protein_coding|Gene:Do1_06_a00001G02298:protein_coding        Gene: Do1_06_a00001G02298-File: RNA_Norway_W_C_DERs_up
./SE_W_C_down.out:Do1_06_a00002 3939270 3939419 -34.189996      Do1_06_a00002   3939270 3939419 line_694;U
pstream:1920|Transcript:Do1_06_a00002G00826V1.1:protein_coding|Gene:Do1_06_a00002G00826:protein_coding;Dow
nstream:3951|Transcript:Do1_06_a00002G00824V1.1:protein_coding|Gene:Do1_06_a00002G00824:protein_coding;Dow
nstream:4581|Transcript:Do1_06_a00002G00830V1.1:protein_coding|Gene:Do1_06_a00002G00830:protein_coding;Ups
tream:68|Transcript:Do1_06_a00002G00827V1.1:protein_coding|Gene:Do1_06_a00002G00827:protein_coding;Downstr
eam:1450|Transcript:Do1_06_a00002G00829V1.1:protein_coding|Gene:Do1_06_a00002G00829:protein_coding;Exon:1:
1:RETAINED|Transcript:Do1_06_a00002G00827V1.1:protein_coding|Gene:Do1_06_a00002G00827:protein_coding;Upstr
eam:2283|Transcript:Do1_06_a00002G00825V1.1:protein_coding|Gene:Do1_06_a00002G00825:protein_coding;Upstrea
m:445|Transcript:Do1_06_a00002G00828V1.1:protein_coding|Gene:Do1_06_a00002G00828:protein_coding;Intergenic
:Do1_06_a00002G00826-Do1_06_a00002G00827        Gene: Do1_06_a00002G00830       File: RNA_Sweden_W_C_DERs_
up
#------------------------------------
./Svalbard_W_C_down.out:Do1_01_a00001   3955454 3955755 -11.929379      Do1_01_a00001   3955454 3955755 li
ne_9;Gene:Do1_01_a00001G00761:protein_coding;Downstream:2971|Transcript:Do1_01_a00001G00762V1.1:protein_co
ding|Gene:Do1_01_a00001G00762:protein_coding;Upstream:1357|Transcript:Do1_01_a00001G00760V1.1:protein_codi
ng|Gene:Do1_01_a00001G00760:protein_coding      Gene: Do1_01_a00001G00762       File: RNA_Alaska_W_C_DERs_
up
./Svalbard_W_C_down.out:Do1_01_a00002   938512  939081  -22.306690      Do1_01_a00002   938512  939081  li
ne_25;Upstream:1312|Transcript:Do1_01_a00002G00147V1.1:protein_coding|Gene:Do1_01_a00002G00147:protein_cod
ing;Upstream:2658|Transcript:Do1_01_a00002G00148V1.1:protein_coding|Gene:Do1_01_a00002G00148:protein_codin
g;Intergenic:Do1_01_a00002G00147-Do1_01_a00002G00148    Gene: Do1_01_a00002G00148       File: RNA_Sweden_W
_C_DERs_up

#------------------------------
./Sweden_W_C_down.out:Do1_00161 21292   21584   -22.554756      Do1_00161       21292   21584   line_12;Up
stream:3020|Transcript:Do1_00161G00004V1.1:protein_coding|Gene:Do1_00161G00004:protein_coding;Exon:2:2:RET
AINED|Transcript:Do1_00161G00005V1.1:protein_coding|Gene:Do1_00161G00005:protein_coding;Intergenic:Do1_001
61G00004-Do1_00161G00005;Upstream:1804|Transcript:Do1_00161G00006V1.1:protein_coding|Gene:Do1_00161G00006:
protein_coding;Downstream:176|Transcript:Do1_00161G00005V1.1:protein_coding|Gene:Do1_00161G00005:protein_c
oding   Gene: Do1_00161G00004   File: RNA_Sweden_W_C_DERs_up
./Sweden_W_C_down.out:Do1_01_a00001     1051991 1052562 -15.741293      Do1_01_a00001   1051991 1052562 li
ne_39;Downstream:486|Transcript:Do1_01_a00001G00213V1.1:protein_coding|Gene:Do1_01_a00001G00213:protein_co
ding;Downstream:5227|Transcript:Do1_01_a00001G00216V1.1:protein_coding|Gene:Do1_01_a00001G00216:protein_co
ding;Upstream:443|Transcript:Do1_01_a00001G00214V1.1:protein_coding|Gene:Do1_01_a00001G00214:protein_codin
g;Upstream:4557|Transcript:Do1_01_a00001G00215V1.1:protein_coding|Gene:Do1_01_a00001G00215:protein_coding;
Intergenic:Do1_01_a00001G00213-Do1_01_a00001G00214;Exon:1:5:RETAINED|Transcript:Do1_01_a00001G00214V1.1:pr
otein_coding|Gene:Do1_01_a00001G00214:protein_coding    Gene: Do1_01_a00001G00214       File: RNA_Sweden_W
_C_DERs_up
./Sweden_W_C_down.out:Do1_01_a00001     1052681 1052804 -13.843138      Do1_01_a00001   1052681 1052804 li
ne_37;Gene:Do1_01_a00001G00214:protein_coding;Downstream:4537|Transcript:Do1_01_a00001G00216V1.1:protein_c
oding|Gene:Do1_01_a00001G00216:protein_coding;Downstream:1176|Transcript:Do1_01_a00001G00213V1.1:protein_c
oding|Gene:Do1_01_a00001G00213:protein_coding;Upstream:3867|Transcript:Do1_01_a00001G00215V1.1:protein_cod
ing|Gene:Do1_01_a00001G00215:protein_coding     Gene: Do1_01_a00001G00214       File: RNA_Sweden_W_C_DERs_
up
./Sweden_W_C_down.out:Do1_01_a00001     1059564 1059798 -24.216919      Do1_01_a00001   1059564 1059798 li
ne_38;Intron:7:7:RETAINED-RETAINED|Transcript:Do1_01_a00001G00216V1.1:protein_coding|Gene:Do1_01_a00001G00
216:protein_coding;Downstream:3463|Transcript:Do1_01_a00001G00214V1.1:protein_coding|Gene:Do1_01_a00001G00
214:protein_coding;Downstream:2692|Transcript:Do1_01_a00001G00215V1.1:protein_coding|Gene:Do1_01_a00001G00
215:protein_coding      Gene: Do1_01_a00001G00214       File: RNA_Sweden_W_C_DERs_up
./Sweden_W_C_down.out:Do1_01_a00001     11134580        11134766        10.475884       Do1_01_a00001   11
134580  11134766        line_114;Downstream:1068|Transcript:Do1_01_a00001G02027V1.1:protein_coding|Gene:Do
1_01_a00001G02027:protein_coding;Intergenic:Do1_01_a00001G02027-Do1_01_a00001G02028     Gene: Do1_01_a0000
1G02027 File: RNA_Alex_W_C_DERs_up
./Sweden_W_C_down.out:Do1_01_a00003     4988127 4988265 40.687564       Do1_01_a00003   4988127 4988265 li
ne_191;Downstream:3389|Transcript:Do1_01_a00003G00798V1.1:protein_coding|Gene:Do1_01_a00003G00798:protein_
coding;Intergenic:Do1_01_a00003G00799-Do1_01_a00003G00800;Downstream:279|Transcript:Do1_01_a00003G00799V1.
1:protein_coding|Gene:Do1_01_a00003G00799:protein_coding;Exon:1:4:RETAINED|Transcript:Do1_01_a00003G00800V
1.1:protein_coding|Gene:Do1_01_a00003G00800:protein_coding;Upstream:121|Transcript:Do1_01_a00003G00800V1.1
:protein_coding|Gene:Do1_01_a00003G00800:protein_coding Gene: Do1_01_a00003G00798       File: RNA_Sweden_W
_C_DERs_up

./Sweden_W_C_down.out:Do1_02_a00001     4213499 4213622 -8.961747       Do1_02_a00001   4213499 4213622 li
ne_319;Exon:1:29:RETAINED|Transcript:Do1_02_a00001G00643V1.1:protein_coding|Gene:Do1_02_a00001G00643:prote
in_coding       Gene: Do1_02_a00001G00643       File: RNA_Sweden_W_C_DERs_up
./Sweden_W_C_down.out:Do1_02_a00003     6757018 6757189 -26.791575      Do1_02_a00003   6757018 6757189 li
ne_387;Exon:1:1:RETAINED|Transcript:Do1_02_a00003G01351V1.1:protein_coding|Gene:Do1_02_a00003G01351:protei
n_coding;Downstream:1172|Transcript:Do1_02_a00003G01350V1.1:protein_coding|Gene:Do1_02_a00003G01350:protei
n_coding;Upstream:3988|Transcript:Do1_02_a00003G01352V1.1:protein_coding|Gene:Do1_02_a00003G01352:protein_
coding  Gene: Do1_02_a00003G01351       File: RNA_Alaska_W_C_DERs_up
./Sweden_W_C_down.out:Do1_02_a00003     9786846 9786959 -50.486828      Do1_02_a00003   9786846 9786959 li
ne_412;Upstream:655|Transcript:Do1_02_a00003G01881V1.1:protein_coding|Gene:Do1_02_a00003G01881:protein_cod
ing;Downstream:4720|Transcript:Do1_02_a00003G01880V1.1:protein_coding|Gene:Do1_02_a00003G01880:protein_cod
ing;Downstream:2715|Transcript:Do1_02_a00003G01882V1.1:protein_coding|Gene:Do1_02_a00003G01882:protein_cod
ing;Intergenic:Do1_02_a00003G01880-Do1_02_a00003G01881  Gene: Do1_02_a00003G01881       File: RNA_Alaska_W
_C_DERs_up
./Sweden_W_C_down.out:Do1_03_a00002     10528980        10529136        -24.067710      Do1_03_a00002   10
528980  10529136        line_647;Exon:2:5:RETAINED|Transcript:Do1_03_a00002G01930V1.1:protein_coding|Gene:
Do1_03_a00002G01930:protein_coding;Downstream:1802|Transcript:Do1_03_a00002G01928V1.1:protein_coding|Gene:
Do1_03_a00002G01928:protein_coding;Upstream:2169|Transcript:Do1_03_a00002G01927V1.1:protein_coding|Gene:Do
1_03_a00002G01927:protein_coding;Downstream:1179|Transcript:Do1_03_a00002G01929V1.1:protein_coding|Gene:Do
1_03_a00002G01929:protein_coding        Gene: Do1_03_a00002G01930       File: RNA_Alaska_W_C_DERs_up
./Sweden_W_C_down.out:Do1_03_a00002     11938686        11938752        -34.326561      Do1_03_a00002   11
938686  11938752        line_666;Upstream:2887|Transcript:Do1_03_a00002G02135V1.1:protein_coding|Gene:Do1_
03_a00002G02135:protein_coding;Intergenic:Do1_03_a00002G02134-Do1_03_a00002G02135;Upstream:433|Transcript:
Do1_03_a00002G02134V1.1:protein_coding|Gene:Do1_03_a00002G02134:protein_coding  Gene: Do1_03_a00002G02134-File: RNA_Alaska_W_C_DERs_up
./Sweden_W_C_down.out:Do1_03_a00002     8261778 8261888 -30.088264      Do1_03_a00002   8261778 8261888 li
ne_631;Intergenic:Do1_03_a00002G01515-Do1_03_a00002G01516;Downstream:959|Transcript:Do1_03_a00002G01515V1.
1:protein_coding|Gene:Do1_03_a00002G01515:protein_coding        Gene: Do1_03_a00002G01515       File: RNA_
Sweden_W_C_DERs_up
./Sweden_W_C_down.out:Do1_04_a00005     175106  175216  20.636618       Do1_04_a00005   175106  175216  li
ne_970;Intergenic:Do1_04_a00005G00023-Do1_04_a00005G00024       Gene: Do1_04_a00005G00024       File: RNA_
Sweden_W_C_DERs_up
./Sweden_W_C_down.out:Do1_05_a00001     5188839 5188896 -15.401527      Do1_05_a00001   5188839 5188896 li
ne_1006;Intergenic:Do1_05_a00001G01056-Do1_05_a00001G01057      Gene: Do1_05_a00001G01057       File: RNA_
Sweden_W_C_DERs_up
./Sweden_W_C_down.out:Do1_05_a00001     8941886 8942095 -33.537105      Do1_05_a00001   8941886 8942095 li
ne_1037;Upstream:2072|Transcript:Do1_05_a00001G01742V1.1:protein_coding|Gene:Do1_05_a00001G01742:protein_c
oding;Intergenic:Do1_05_a00001G01743-Do1_05_a00001G01744;Upstream:263|Transcript:Do1_05_a00001G01743V1.1:p
rotein_coding|Gene:Do1_05_a00001G01743:protein_coding;Upstream:3998|Transcript:Do1_05_a00001G01744V1.1:pro
tein_coding|Gene:Do1_05_a00001G01744:protein_coding     Gene: Do1_05_a00001G01744       File: RNA_Sweden_W
_C_DERs_up
./Sweden_W_C_down.out:Do1_05_a00002     107390  107696  -23.376064      Do1_05_a00002   107390  107696  li
ne_1078;Exon:7:7:RETAINED|Transcript:Do1_05_a00002G00014V1.1:protein_coding|Gene:Do1_05_a00002G00014:prote
in_coding;Upstream:1998|Transcript:Do1_05_a00002G00015V1.1:protein_coding|Gene:Do1_05_a00002G00015:protein
_coding;Upstream:4612|Transcript:Do1_05_a00002G00016V1.1:protein_coding|Gene:Do1_05_a00002G00016:protein_c
oding   Gene: Do1_05_a00002G00016       File: RNA_Sweden_W_C_DERs_up
./Sweden_W_C_down.out:Do1_05_a00004     582199  582458  8.291380        Do1_05_a00004   582199  582458  li
ne_1165;Downstream:1730|Transcript:Do1_05_a00004G00096V1.1:protein_coding|Gene:Do1_05_a00004G00096:protein
_coding;Intergenic:Do1_05_a00004G00096-Do1_05_a00004G00097      Gene: Do1_05_a00004G00097       File: RNA_
Sweden_W_C_DERs_up
./Sweden_W_C_down.out:Do1_06_a00002     12074391        12074558        -23.236069      Do1_06_a00002   12
074391  12074558        line_1393;Upstream:45|Transcript:Do1_06_a00002G02285V1.1:protein_coding|Gene:Do1_0
6_a00002G02285:protein_coding;Upstream:4675|Transcript:Do1_06_a00002G02283V1.1:protein_coding|Gene:Do1_06_
a00002G02283:protein_coding;Downstream:3505|Transcript:Do1_06_a00002G02284V1.1:protein_coding|Gene:Do1_06_
a00002G02284:protein_coding;Intergenic:Do1_06_a00002G02285-Do1_06_a00002G02286  Gene: Do1_06_a00002G02285-File: RNA_Alaska_W_C_DERs_up
./Sweden_W_C_down.out:Do1_06_a00002     5342908 5343304 21.881194       Do1_06_a00002   5342908 5343304 li
ne_1334;Intergenic:Do1_06_a00002G01112-Do1_06_a00002G01113      Gene: Do1_06_a00002G01112       File: RNA_
Sweden_W_C_DERs_up
./Sweden_W_C_down.out:Do1_07_a00002     7853137 7853522 -15.313078      Do1_07_a00002   7853137 7853522 li
ne_1515;Exon:1:1:RETAINED|Transcript:Do1_07_a00002G01552V1.1:protein_coding|Gene:Do1_07_a00002G01552:prote
in_coding;Upstream:3802|Transcript:Do1_07_a00002G01550V1.1:protein_coding|Gene:Do1_07_a00002G01550:protein
_coding;Upstream:1230|Transcript:Do1_07_a00002G01551V1.1:protein_coding|Gene:Do1_07_a00002G01551:protein_c
oding;Upstream:4671|Transcript:Do1_07_a00002G01549V1.1:protein_coding|Gene:Do1_07_a00002G01549:protein_cod
ing;Upstream:0|Transcript:Do1_07_a00002G01552V1.1:protein_coding|Gene:Do1_07_a00002G01552:protein_coding;D
ownstream:2042|Transcript:Do1_07_a00002G01553V1.1:protein_coding|Gene:Do1_07_a00002G01553:protein_coding;I
ntergenic:Do1_07_a00002G01552-Do1_07_a00002G01553       Gene: Do1_07_a00002G01553       File: RNA_Sweden_W
_C_DERs_up
./Sweden_W_C_down.out:Do1_07_a00005     1435801 1435956 7.980132        Do1_07_a00005   1435801 1435956 li
ne_1672;Upstream:1278|Transcript:Do1_07_a00005G00229V1.1:protein_coding|Gene:Do1_07_a00005G00229:protein_c
oding;Intergenic:Do1_07_a00005G00228-Do1_07_a00005G00229;Upstream:3503|Transcript:Do1_07_a00005G00230V1.1:
protein_coding|Gene:Do1_07_a00005G00230:protein_coding;Upstream:3608|Transcript:Do1_07_a00005G00228V1.1:pr
otein_coding|Gene:Do1_07_a00005G00228:protein_coding    Gene: Do1_07_a00005G00229       File: RNA_Alaska_W
_C_DERs_up

