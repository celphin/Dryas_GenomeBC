########################################
# Dryas differential gene expression data
# Determining differentially expressed genes, plotting results
# Sept 2023 - Sweden, Feb 2023 other sites
#####################################

###########################################################
# try differential expression analysis - plotting
# EdgeR

tmux new-session -s RNA
tmux attach-session -t RNA

# copy data to shared folder for Marie
#---------------------------------
# https://github.com/owensgl/biol525D/tree/master/Topic_6
# https://bioconductor.org/packages/release/bioc/html/edgeR.html
# http://combine-australia.github.io/RNAseq-R/slides/RNASeq_filtering_qc.pdf 

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
# if (!require("BiocManager", quietly = TRUE))
    # install.packages("BiocManager")
# BiocManager::install(version = "3.15")

# BiocManager::install("edgeR")
# install.packages('gplots')

#paste it in here (i.e. replace my path with yours):
#setwd ("/scratch/msandler/RNAseq_analysis")
setwd ("/home/celphin/projects/rrg-rieseber-ac/rpp-rieseber/celphin/Dryas/RNAseq_analysis/")

#find your working directory:
getwd()

#load the libraries you will need 
library(edgeR)
library(gplots)
library(stringr)

#read in the data
mydata <- read.table("./RNA_site_tables/gene_names_expression_table1.txt", header=TRUE)

CDS <- mydata[,1]
exp_data <- mydata[,c(2:length(mydata))]

#extract and turn the column names into a factor
treat <- as.factor (sapply (strsplit(colnames(exp_data),split = "_"),"[[",1))
site <- as.factor (sapply (strsplit(colnames(exp_data),split = "_"),"[[",2))
mixed<- paste0(treat, "_", site)

treat_site <- paste0(as.character(site), as.character(treat))
treat_site <- as.factor(treat_site)

#make a DGElist object out of the data
dge <- DGEList (counts = exp_data, genes=rownames(exp_data))

#calculate the normalization factors to adjust the effective library size relative to other libraries in the dataset (norm.factor that minimizes log fold change among samples)
dge <- calcNormFactors (dge)

#----------------------
# also see https://rpubs.com/jrgonzalezISGlobal/enrichment
# https://bioinformatics-core-shared-training.github.io/cruk-bioinf-sschool/Day3/rnaSeq_DE.pdf

# filter and see dimensions
keep.exprs <- filterByExpr(dge)
dge.filt <- dge[keep.exprs,]
dim(dge.filt)
# 15856    58

dge.filt$samples$group <-  as.factor(str_extract(rownames(dge.filt$samples), "^[WC]"))
dge.filt$samples$site <-  as.factor(str_extract(rownames(dge.filt$samples), "(?<=_)[A-Za-z]+"))

# plot
designMat <- model.matrix( ~ group+site, data=dge.filt$samples)

png("./plots/voom_mean_variance_trend.png", width = 700, height = 500)
v <- voom(dge.filt, design = designMat, plot = TRUE)
dev.off()

fit <- lmFit(v, designMat)
fit <- eBayes(fit)
topTable(fit)

                            # logFC  AveExpr         t      P.Value adj.P.Val
# Do1_02_a00003G01059V1.1 -0.3736246 5.783677 -4.585718 2.418313e-05 0.2272247
# Do1_05_a00001G00678V1.1  0.4868721 4.996119  4.132682 1.152688e-04 0.5415327
# Do1_05_a00001G01393V1.1 -0.2364523 5.309064 -4.022222 1.668231e-04 0.6269878
# Do1_01_a00004G01253V1.1 -0.3654033 2.920314 -4.754991 1.326829e-05 0.2272247
# Do1_01_a00001G00911V1.1  0.6073389 4.125502  4.154200 1.072021e-04 0.5415327
# Do1_03_a00002G01021V1.1 -0.3666808 6.127259 -3.616250 6.219597e-04 0.8243400
# Do1_04_a00001G00236V1.1 -0.6364836 4.384837 -3.744945 4.129610e-04 0.8243400
# Do1_04_a00001G02342V1.1  0.2748636 4.832356  3.649093 5.606347e-04 0.8243400
# Do1_01_a00001G00364V1.1 -0.3931061 5.075336 -3.503560 8.847180e-04 0.8243400
# Do1_04_a00001G01745V1.1 -0.2929932 6.026556 -3.402461 1.207479e-03 0.8243400
                                 # B
# Do1_02_a00003G01059V1.1 -0.1606955
# Do1_05_a00001G00678V1.1 -1.0265966
# Do1_05_a00001G01393V1.1 -1.1399468
# Do1_01_a00004G01253V1.1 -1.2243403
# Do1_01_a00001G00911V1.1 -1.2939314
# Do1_03_a00002G01021V1.1 -1.7415305
# Do1_04_a00001G00236V1.1 -1.7825713
# Do1_04_a00001G02342V1.1 -1.8051183
# Do1_01_a00001G00364V1.1 -1.9829141
# Do1_04_a00001G01745V1.1 -2.0717539

tt <- topTable(fit, n=Inf)
mask <- tt$adj.P.Val < 0.05 & abs(tt$logFC) > log2(2)
deGenes <- rownames(tt[mask, ])
head(deGenes)
# character(0)

length(deGenes)
# 0

geneUniverse <- rownames(tt)
length(geneUniverse)
# 18792

#---------------------
dge <- estimateGLMCommonDisp(dge, design=designMat)
dge <- estimateGLMTrendedDisp(dge, design=designMat)
dge <- estimateGLMTagwiseDisp(dge, design=designMat)

png("./plots/Biological_coeff_variation.png", width = 700, height = 500)
plotBCV(dge)
dev.off()

# We can now find our differentially expressed genes. After fitting the model, we can use the topTags() 
# function to explore the results, and set theresholds to identify subsets of differentially expressed genes
fit <- glmFit(dge, designMat)
lrt <- glmLRT(fit, coef=4)
edgeR_result <- topTags(lrt)
#?topTags
toptags_Warming <- topTags(lrt,n=15000)
save(toptags_Warming, file='edgeR_TopTags.RData') #We will need this later

# Finally, we can plot the log-fold changes of all the genes, and the highlight those that are differentially expressed
#?decideTests
deGenes <- decideTestsDGE(lrt, p=0.001)
deGenes <- rownames(lrt)[as.logical(deGenes)]

png("./plots/Plot_smear.png", width = 700, height = 500)
plotSmear(lrt, de.tags=deGenes)
abline(h=c(-1, 1), col=2)
dev.off()

#----------------
# plot deGenes 

DE_genes <- edgeR_result$table  # Extract the result table
DE_genes <- DE_genes[DE_genes$FDR < 0.05 & abs(DE_genes$logFC) > 1, ]


[which(dge.filt$genes$genes==DE_genes#genes)]

#generate a multi-dimensional scaling plot
col_treat <- as.character (treat)
col_treat [col_treat == "C"] <- "blue"
col_treat [col_treat == "W"] <- "red"

col_site <- as.character (site)
col_site [col_treat == "Alex"] <- 15 # square
col_site [col_treat == "Norw"] <- 16 # circle
col_site [col_treat == "Alas"] <- 17 # triangle
col_site [col_treat == "Swed"] <- 18 # diamond

#col_treat [col_treat == "H"] <- "green"
#col_treat [col_treat == "L"] <- "yellow"

#plot the MDS graph
png("./plots/W_C_RNA_MDS.png", width = 700, height = 500)
plotMDS (filtered_dge, col = col_treat)
dev.off()

png("./plots/Site_RNA_MDS.png", width = 700, height = 500)
plotMDS (filtered_dge, col = col_treat, pch=col_site)
dev.off()


##################################
# https://bioinformatics-core-shared-training.github.io/cruk-bioinf-sschool/Day3/rnaSeq_DE.pdf

#calculate the counts per million
cpm.dge <- cpm(dge)

#filter out the genes with less than 1 cpm in 6 or fewer libraries (a somewhat arbitrary choice). 
#Genes are usually dropped if they can't possibly be expressed in all the samples for any of the conditions.
filtered_dge <- dge[rowSums(cpm.dge > 1) >= 6,]
cpm.filtered_dge <- cpm.dge[rowSums(cpm.dge > 1) >= 6,]

#generate a multi-dimensional scaling plot
col_treat <- as.character (treat)
col_treat [col_treat == "C"] <- "blue"
col_treat [col_treat == "W"] <- "red"

col_site <- as.character (site)
col_site [col_treat == "Alex"] <- 15 # square
col_site [col_treat == "Norw"] <- 16 # circle
col_site [col_treat == "Alas"] <- 17 # triangle
col_site [col_treat == "Swed"] <- 18 # diamond

#col_treat [col_treat == "H"] <- "green"
#col_treat [col_treat == "L"] <- "yellow"

#plot the MDS graph
png("./plots/W_C_RNA_MDS.png", width = 700, height = 500)
plotMDS (filtered_dge, col = col_treat)
dev.off()

png("./plots/Site_RNA_MDS.png", width = 700, height = 500)
plotMDS (filtered_dge, col = col_treat, pch=col_site)
dev.off()

#Set the model to use. This one includes the intercept, but other models can be specified that omit the intercept or that have more complex designs. See EdgeR manual for details.
design <- model.matrix (~treat_site)

#fit the common + tagwise dispersion models
filtered_dge <- estimateGLMCommonDisp (filtered_dge, design)
filtered_dge <- estimateGLMTagwiseDisp (filtered_dge, design)

#fit a GLM to the data using the tagwise dispersion estimate
glm.filtered_dge <- glmFit (filtered_dge, design, dispersion = filtered_dge$tagwise.dispersion)
lrt.filtered_dge <- glmLRT (glm.filtered_dge)

#get the topTags out of the model
#top <- topTags (lrt.filtered_dge, n = 1000)$table
#top <- topTags (lrt.filtered_dge, n = 200)$table
top <- topTags (lrt.filtered_dge, n = 2000)$table # p-value all < 5e-02 

#------------------
# gets the DERs and their information here

fdr<-p.adjust(lrt.filtered_dge$table$PValue, method='fdr')

dim(lrt.filtered_dge$table[fdr<0.05,])
length(lrt.filtered_dge$table$PValue[lrt.filtered_dge$table$PValue<0.05])

lrt.list3 <- cbind(lrt.filtered_dge$table, fdr)

DERs_FDR <- lrt.list3[lrt.list3$fdr<0.05,]

dim(DERs_FDR)
# treat*site 21
# treat 24
# treat+site 6058 
# treat_site 2600

write.table(DERs_FDR, file = "./RNA_DERs_Oct2023_W_C_Total.txt", quote = FALSE, row.names=TRUE, col.names=TRUE, sep="\t")


                            logFC     logCPM       LR       PValue         fdr
# Do1_01_a00001G00030V1.1 -1.555119  0.6418423 18.68095 1.545183e-05 0.018544255
# Do1_01_a00001G00085V1.1 -2.313605  2.3305637 28.13965 1.128694e-07 0.001221876
# Do1_01_a00001G00274V1.1 -1.558471  3.0357361 19.20198 1.175912e-05 0.016283671
# Do1_01_a00001G01424V1.1 -1.166022  4.6918350 20.21781 6.910660e-06 0.016283671
# Do1_01_a00003G00772V1.1 -2.098630  0.2479064 19.35295 1.086511e-05 0.016283671
# Do1_02_a00003G01208V1.1 -1.947013  3.7626088 22.40433 2.208758e-06 0.006627010
# Do1_03_a00001G00983V1.1 -1.592588  2.9575249 16.33764 5.300064e-05 0.041483369
# Do1_03_a00001G01365V1.1  3.317710  1.2698965 27.01589 2.017897e-07 0.001221876
# Do1_03_a00001G01606V1.1  2.201592  0.9893801 19.46556 1.024302e-05 0.016283671
# Do1_03_a00002G00567V1.1 -1.476117  5.8557118 16.41551 5.086744e-05 0.041483369
# Do1_03_a00002G00615V1.1  2.974555  1.6162456 18.23439 1.953201e-05 0.021975958
# Do1_04_a00001G01481V1.1 -2.287186  4.5534372 25.11579 5.398887e-07 0.002429769
# Do1_04_a00002G00429V1.1 -1.825682  0.2126465 20.02648 7.637724e-06 0.016283671
# Do1_04_a00004G00131V1.1  3.032924  3.0171184 18.96812 1.329212e-05 0.017091765
# Do1_05_a00003G00079V1.1 -2.978305 -0.3997939 26.99841 2.036233e-07 0.001221876
# Do1_05_a00003G00319V1.1 -2.590086  2.5616164 15.95157 6.498387e-05 0.048743320
# Do1_06_a00001G01291V1.1  2.247983 -0.4754747 16.91118 3.917021e-05 0.035257104
# Do1_06_a00001G02298V1.1  1.402241  0.3019848 19.65192 9.290992e-06 0.016283671
# Do1_06_a00001G02407V1.1 -1.969638  2.1677893 17.49011 2.888057e-05 0.027363585
# Do1_06_a00001G02519V1.1 -2.410007  2.4129034 16.53544 4.774933e-05 0.040932543
# Do1_06_a00001G02551V1.1  2.743289  0.2693729 22.99003 1.628441e-06 0.005863041
# Do1_06_a00002G00928V1.1  1.558237  1.7748365 17.72458 2.553061e-05 0.027035414
# Do1_07_a00002G00350V1.1 -1.835707  2.9102510 19.38836 1.066553e-05 0.016283671
# Do1_07_a00002G02243V1.1 -1.235973  3.6158298 17.54269 2.809291e-05 0.027363585

#-----------------------------------
# try another method
install.packages("statmod")
library(statmod)
design <- model.matrix (~0+treat_site)
fit <- glmQLFit(filtered_dge, design, robust=TRUE)
head(fit$coefficients)

#-----
con <- makeContrasts(treat_siteAlasW - treat_siteAlasC, levels=design)
qlf <- glmQLFTest(fit, contrast=con)
topTags(qlf)
summary(decideTests(qlf))

       # -1*treat_siteAlasC 1*treat_siteAlasW
# Down                                      9
# NotSig                                17989
# Up                                        4

Alaska_W_C_DERs <- qlf$table[which(!decideTests(qlf)==0),]
write.table(Alaska_W_C_DERs, file = "./RNA_Alaska_W_C_DERs.txt", quote = FALSE, row.names=TRUE, col.names=TRUE, sep="\t")

#-----
con <- makeContrasts(treat_siteNorwW - treat_siteNorwC, levels=design)
qlf <- glmQLFTest(fit, contrast=con)
topTags(qlf)
summary(decideTests(qlf))

       # -1*treat_siteNorwC 1*treat_siteNorwW
# Down                                      0
# NotSig                                18002
# Up                                        0

Norway_W_C_DERs <- qlf$table[which(!decideTests(qlf)==0),]
write.table(Norway_W_C_DERs, file = "./RNA_Norway_W_C_DERs.txt", quote = FALSE, row.names=TRUE, col.names=TRUE, sep="\t")

#-----
con <- makeContrasts(treat_siteSwedW - treat_siteSwedC, levels=design)
qlf <- glmQLFTest(fit, contrast=con)
topTags(qlf)
summary(decideTests(qlf))


       # -1*treat_siteSwedC 1*treat_siteSwedW
# Down                                      2
# NotSig                                17998
# Up                                        2

Sweden_W_C_DERs <-  qlf$table[which(!decideTests(qlf)==0),]
write.table(Sweden_W_C_DERs, file = "./RNA_Sweden_W_C_DERs.txt", quote = FALSE, row.names=TRUE, col.names=TRUE, sep="\t")

#-----
con <- makeContrasts(treat_siteAlexW  - treat_siteAlexC, levels=design)
qlf <- glmQLFTest(fit, contrast=con)
topTags(qlf)
summary(decideTests(qlf))

       # -1*treat_siteAlexC 1*treat_siteAlexW
# Down                                    129
# NotSig                                17782
# Up                                       91

Alex_W_C_DERs <- qlf$table[which(!decideTests(qlf)==0),]
write.table(Alex_W_C_DERs, file = "./RNA_Alex_W_C_DERs.txt", quote = FALSE, row.names=TRUE, col.names=TRUE, sep="\t")

#-----
con <- makeContrasts(treat_siteSeedW - treat_siteSeedC, levels=design)
qlf <- glmQLFTest(fit, contrast=con)
topTags(qlf)
summary(decideTests(qlf))
       # -1*treat_siteSeedC 1*treat_siteSeedW
# Down                                      7
# NotSig                                17988
# Up                                        7


Seedling_W_C_DERs <- qlf$table[which(!decideTests(qlf)==0),]
write.table(Seedling_W_C_DERs, file = "./RNA_Seedling_W_C_DERs.txt", quote = FALSE, row.names=TRUE, col.names=TRUE, sep="\t")


#---------------------------
plotMD(qlf)
tr <- glmTreat(fit, contrast=con, lfc=log2(1.2))
topTags(tr)
summary(decideTests(tr))
plotMD(tr)

#----------------------------------
#make a heatmap by getting the counts per million from each gene and turning them relative proportions (columns add up to 1)
sub1 <- colSums (cpm.filtered_dge)
sub2 <- matrix (rep(sub1,nrow (cpm.filtered_dge)), c (nrow (cpm.filtered_dge),ncol(cpm.filtered_dge)),byrow = TRUE)
sub3 <- cpm.filtered_dge / sub2

#subset this matrix to just get the genes from above
names1 <- row.names (sub3)
names2 <- row.names (top)
index2 <- names1 %in% names2
heatmap1 <- sub3[index2,]

#------------------------------
# All sites
##     par(mar = c(bottom, left, top, right))
# plot used for GenomeBC report
#play around with the options to make the plot fit what you like for options type ?heatmap.2

# order by site
# http://www.sthda.com/english/wiki/reordering-data-frame-columns-in-r
colnames(heatmap1)
# ("C_Alas_16",   "C_Alas_1",    "C_Alas_20",   "C_Alas_5",    "C_Alex_11",
# "C_Alex_21",   "C_Alex_22",   "C_Alex_23",   "C_Alex_24",   "C_Alex_2",
# "C_Alex_6",    "C_Alex_7",    "C_Norw_10",   "C_Norw_1",    "C_Norw_3_B",
# "C_Norw_5",    "C_Norw_6",    "C_Norw_7",    "C_Seed_1",    "C_Seed_2",
# "C_Seed_3",    "C_Swed_11",   "C_Swed_15",   "C_Swed_16",   "C_Swed_2",
# "C_Swed_3",    "C_Swed_5",    "C_Swed_8",    "W_Alas_12",   "W_Alas_13",
# "W_Alas_15",   "W_Alas_1",    "W_Alas_2",    "W_Alas_3",    "W_Alas_5",
# "W_Alas_9",    "W_Alex_11",   "W_Alex_1",    "W_Alex_21",   "W_Alex_22",
# "W_Alex_6",    "W_Norw_10_B" "W_Norw_2",    "W_Norw_3_B"  "W_Norw_4_B",
# "W_Norw_6",    "W_Norw_7",    "W_Seed_1",    "W_Seed_2",    "W_Seed_3",
# "W_Swed_10",   "W_Swed_13",   "W_Swed_14",   "W_Swed_16",   "W_Swed_17",
# "W_Swed_2",    "W_Swed_5",    "W_Swed_6")

col_order <- c("C_Alas_16",   "C_Alas_1",    "C_Alas_20",   "C_Alas_5",  "W_Alas_15",   "W_Alas_1",    "W_Alas_2",    "W_Alas_3",    "W_Alas_5",  "W_Alas_9",  "W_Alas_12",   "W_Alas_13",
"C_Alex_11",  "C_Alex_21",   "C_Alex_22",   "C_Alex_23",   "C_Alex_24",   "C_Alex_2",  "C_Alex_6",    "C_Alex_7",    "W_Alex_11",   "W_Alex_1",    "W_Alex_21",   "W_Alex_22",  "W_Alex_6",   
"C_Norw_10",   "C_Norw_1",    "C_Norw_3_B",  "C_Norw_5",    "C_Norw_6",    "C_Norw_7",    "W_Norw_10_B", "W_Norw_2",    "W_Norw_3_B",  "W_Norw_4_B", "W_Norw_6",    "W_Norw_7",    
"C_Swed_15",   "C_Swed_16",   "C_Swed_2", "C_Swed_11",   "C_Swed_3",    "C_Swed_5",    "C_Swed_8",   "W_Swed_10",   "W_Swed_13",   "W_Swed_14",   "W_Swed_16",   "W_Swed_17",  "W_Swed_2",    "W_Swed_5",    "W_Swed_6", 
"C_Seed_1",    "C_Seed_2", "C_Seed_3", "W_Seed_1",    "W_Seed_2",    "W_Seed_3")

heatmap2 <- heatmap1[, col_order]

png("./plots/W_C_RNA_heatmap.png", width = 3200, height = 1000)
heatmap.2 (heatmap2, scale = "row", trace = "none", Colv = FALSE, labCol = colnames(heatmap2), labRow = rownames(heatmap2), col='rainbow', dendrogram = "none", colsep =c(5), margins =c(20,70), sepcolor = "white", sepwidth = 0.1, cexRow=2, cexCol=2)
dev.off()

png("./plots/W_C_RNA_heatmap_legend.png", width = 1000, height = 800)
heatmap.2 (heatmap2, scale = "row", trace = "none", notecex=3, notecol="black", tracecol="black", Colv = FALSE, labCol = colnames(heatmap2), labRow = rownames(heatmap2), col='rainbow', dendrogram = "none", colsep =c(5), margins =c(0.2, 0.2), sepcolor = "white", sepwidth = 0.1, cexRow=6, cexCol=6)
dev.off()

# write out file of genes that differ a lot

write.table(heatmap2, file = "./RNA_DER_May2023_W_C_Total.txt", quote = FALSE, row.names=TRUE, col.names=TRUE, sep="\t")

#------------------------------
# Alex

#read in the data
#mydata <- read.table("./RNA_site_tables/gene_names_expression_table1.txt", header=TRUE)

CDS <- mydata[,1]
exp_data <- mydata[,c(6:13,38:42)]

 # [1] "C_Alas_16"   "C_Alas_1"    "C_Alas_20"   "C_Alas_5"    "C_Alex_11"
 # [6] "C_Alex_21"   "C_Alex_22"   "C_Alex_23"   "C_Alex_24"   "C_Alex_2"
# [11] "C_Alex_6"    "C_Alex_7"    "C_Norw_10"   "C_Norw_1"    "C_Norw_3_B"
# [16] "C_Norw_5"    "C_Norw_6"    "C_Norw_7"    "C_Seed_1"    "C_Seed_2"
# [21] "C_Seed_3"    "C_Swed_11"   "C_Swed_15"   "C_Swed_16"   "C_Swed_2"
# [26] "C_Swed_3"    "C_Swed_5"    "C_Swed_8"    "W_Alas_12"   "W_Alas_13"
# [31] "W_Alas_15"   "W_Alas_1"    "W_Alas_2"    "W_Alas_3"    "W_Alas_5"
# [36] "W_Alas_9"    "W_Alex_11"   "W_Alex_1"    "W_Alex_21"   "W_Alex_22"
# [41] "W_Alex_6"    "W_Norw_10_B" "W_Norw_2"    "W_Norw_3_B"  "W_Norw_4_B"
# [46] "W_Norw_6"    "W_Norw_7"    "W_Seed_1"    "W_Seed_2"    "W_Seed_3"
# [51] "W_Swed_10"   "W_Swed_13"   "W_Swed_14"   "W_Swed_16"   "W_Swed_17"
# [56] "W_Swed_2"    "W_Swed_5"    "W_Swed_6"


#extract and turn the column names into a factor
treat <- as.factor (sapply (strsplit(colnames(exp_data),split = "_"),"[[",1))
site <- as.factor (sapply (strsplit(colnames(exp_data),split = "_"),"[[",2))

treat_site <- paste0(as.character(site), as.character(treat))
treat_site <- as.factor(treat_site)

#make a DGElist object out of the data
dge <- DGEList (counts = exp_data, group = treat)

#calculate the normalization factors to adjust the effective library size relative to other libraries in the dataset (norm.factor that minimizes log fold change among samples)
dge <- calcNormFactors (dge)

#calculate the counts per million
cpm.dge <- cpm(dge)

#filter out the genes with less than 1 cpm in 6 or fewer libraries (a somewhat arbitrary choice). 
#Genes are usually dropped if they can't possibly be expressed in all the samples for any of the conditions.
filtered_dge <- dge[rowSums(cpm.dge > 1) >= 6,]
cpm.filtered_dge <- cpm.dge[rowSums(cpm.dge > 1) >= 6,]

#generate a multi-dimensional scaling plot
col_treat <- as.character (treat)
col_treat [col_treat == "C"] <- "blue"
col_treat [col_treat == "W"] <- "red"

col_site <- as.character (site)
col_site [col_treat == "Alex"] <- 15 # square
col_site [col_treat == "Norw"] <- 16 # circle
col_site [col_treat == "Alas"] <- 17 # triangle
col_site [col_treat == "Swed"] <- 18 # diamond

#col_treat [col_treat == "H"] <- "green"
#col_treat [col_treat == "L"] <- "yellow"

#plot the MDS graph
png("./plots/Alex_W_C_RNA_MDS.png", width = 700, height = 500)
plotMDS (filtered_dge, col = col_treat)
dev.off()

#Set the model to use. This one includes the intercept, but other models can be specified that omit the intercept or that have more complex designs. See EdgeR manual for details.
design <- model.matrix (~treat)

#fit the common + tagwise dispersion models
filtered_dge <- estimateGLMCommonDisp (filtered_dge, design)
filtered_dge <- estimateGLMTagwiseDisp (filtered_dge, design)

#fit a GLM to the data using the tagwise dispersion estimate
glm.filtered_dge <- glmFit (filtered_dge, design, dispersion = filtered_dge$tagwise.dispersion)
lrt.filtered_dge <- glmLRT (glm.filtered_dge)

#get the topTags out of the model
#top <- topTags (lrt.filtered_dge, n = 1000)$table
#top <- topTags (lrt.filtered_dge, n = 200)$table
top <- topTags (lrt.filtered_dge, n = 500)$table # p-value all < 5e-02, v2

fdr<-p.adjust(lrt.filtered_dge$table$PValue, method='fdr')

dim(lrt.filtered_dge$table[fdr<0.05,])
length(lrt.filtered_dge$table$PValue[lrt.filtered_dge$table$PValue<0.05])

#make a heatmap by getting the counts per million from each gene and turning them relative proportions (columns add up to 1)
sub1 <- colSums (cpm.filtered_dge)
sub2 <- matrix (rep(sub1,nrow (cpm.filtered_dge)), c (nrow (cpm.filtered_dge),ncol(cpm.filtered_dge)),byrow = TRUE)
sub3 <- cpm.filtered_dge / sub2

#subset this matrix to just get the genes from above
names1 <- row.names (sub3)
names2 <- row.names (top)
index2 <- names1 %in% names2
heatmap1 <- sub3[index2,]

##     par(mar = c(bottom, left, top, right))
# plot used for GenomeBC report
#play around with the options to make the plot fit what you like for options type ?heatmap.2

# order by site
# http://www.sthda.com/english/wiki/reordering-data-frame-columns-in-r
col_order <- c("C_Alex_11",  "C_Alex_21",   "C_Alex_22",   "C_Alex_23",   "C_Alex_24",   "C_Alex_2",  "C_Alex_6",    "C_Alex_7",    "W_Alex_11",   "W_Alex_1",    "W_Alex_21",   "W_Alex_22",  "W_Alex_6")
heatmap2 <- heatmap1[, col_order]

png("./plots/Alex_W_C_RNA_heatmap2.png", width = 3200, height = 2500)
heatmap.2 (heatmap1, scale = "row", trace = "none", Colv = FALSE, labCol = colnames(heatmap1), labRow = rownames(heatmap1), col='rainbow', dendrogram = "none", colsep =c(8), margins =c(20,70), sepcolor = "white", sepwidth = 0.1, cexRow=3, cexCol=6)
dev.off()

png("./plots/Alex_W_C_RNA_heatmap_legend2.png", width = 1000, height = 800)
heatmap.2 (heatmap1, scale = "row", trace = "none", notecex=3, notecol="black", tracecol="black", Colv = FALSE, labCol = colnames(heatmap1), labRow = rownames(heatmap1), col='rainbow', dendrogram = "none", colsep =c(5), margins =c(0.2, 0.2), sepcolor = "white", sepwidth = 0.1, cexRow=6, cexCol=6)
dev.off()

# write out file of genes that differ a lot

write.table(heatmap2, file = "./RNA_DER_May2023_W_C_Alex.txt", quote = FALSE, row.names=TRUE, col.names=TRUE, sep="\t")

#------------------------------
# Norway

#read in the data
#mydata <- read.table("./RNA_site_tables/gene_names_expression_table1.txt", header=TRUE)

CDS <- mydata[,1]
exp_data <- mydata[,c(14:19, 43:48)]

 # [1] "C_Alas_16"   "C_Alas_1"    "C_Alas_20"   "C_Alas_5"    "C_Alex_11"
 # [6] "C_Alex_21"   "C_Alex_22"   "C_Alex_23"   "C_Alex_24"   "C_Alex_2"
# [11] "C_Alex_6"    "C_Alex_7"    "C_Norw_10"   "C_Norw_1"    "C_Norw_3_B"
# [16] "C_Norw_5"    "C_Norw_6"    "C_Norw_7"    "C_Seed_1"    "C_Seed_2"
# [21] "C_Seed_3"    "C_Swed_11"   "C_Swed_15"   "C_Swed_16"   "C_Swed_2"
# [26] "C_Swed_3"    "C_Swed_5"    "C_Swed_8"    "W_Alas_12"   "W_Alas_13"
# [31] "W_Alas_15"   "W_Alas_1"    "W_Alas_2"    "W_Alas_3"    "W_Alas_5"
# [36] "W_Alas_9"    "W_Alex_11"   "W_Alex_1"    "W_Alex_21"   "W_Alex_22"
# [41] "W_Alex_6"    "W_Norw_10_B" "W_Norw_2"    "W_Norw_3_B"  "W_Norw_4_B"
# [46] "W_Norw_6"    "W_Norw_7"    "W_Seed_1"    "W_Seed_2"    "W_Seed_3"
# [51] "W_Swed_10"   "W_Swed_13"   "W_Swed_14"   "W_Swed_16"   "W_Swed_17"
# [56] "W_Swed_2"    "W_Swed_5"    "W_Swed_6"


#extract and turn the column names into a factor
treat <- as.factor (sapply (strsplit(colnames(exp_data),split = "_"),"[[",1))
site <- as.factor (sapply (strsplit(colnames(exp_data),split = "_"),"[[",2))

treat_site <- paste0(as.character(site), as.character(treat))
treat_site <- as.factor(treat_site)

#make a DGElist object out of the data
dge <- DGEList (counts = exp_data, group = treat)

#calculate the normalization factors to adjust the effective library size relative to other libraries in the dataset (norm.factor that minimizes log fold change among samples)
dge <- calcNormFactors (dge)

#calculate the counts per million
cpm.dge <- cpm(dge)

#filter out the genes with less than 1 cpm in 6 or fewer libraries (a somewhat arbitrary choice). 
#Genes are usually dropped if they can't possibly be expressed in all the samples for any of the conditions.
filtered_dge <- dge[rowSums(cpm.dge > 1) >= 6,]
cpm.filtered_dge <- cpm.dge[rowSums(cpm.dge > 1) >= 6,]

#generate a multi-dimensional scaling plot
col_treat <- as.character (treat)
col_treat [col_treat == "C"] <- "blue"
col_treat [col_treat == "W"] <- "red"

col_site <- as.character (site)
col_site [col_treat == "Alex"] <- 15 # square
col_site [col_treat == "Norw"] <- 16 # circle
col_site [col_treat == "Alas"] <- 17 # triangle
col_site [col_treat == "Swed"] <- 18 # diamond

#col_treat [col_treat == "H"] <- "green"
#col_treat [col_treat == "L"] <- "yellow"

#plot the MDS graph
png("./plots/Norway_W_C_RNA_MDS.png", width = 700, height = 500)
plotMDS (filtered_dge, col = col_treat)
dev.off()


#Set the model to use. This one includes the intercept, but other models can be specified that omit the intercept or that have more complex designs. See EdgeR manual for details.
design <- model.matrix (~treat)

#fit the common + tagwise dispersion models
filtered_dge <- estimateGLMCommonDisp (filtered_dge, design)
filtered_dge <- estimateGLMTagwiseDisp (filtered_dge, design)

#fit a GLM to the data using the tagwise dispersion estimate
glm.filtered_dge <- glmFit (filtered_dge, design, dispersion = filtered_dge$tagwise.dispersion)
lrt.filtered_dge <- glmLRT (glm.filtered_dge)

#get the topTags out of the model
#top <- topTags (lrt.filtered_dge, n = 1000)$table
#top <- topTags (lrt.filtered_dge, n = 200)$table
top <- topTags (lrt.filtered_dge, n = 500)$table # p-value all < 5e-02, v2

fdr<-p.adjust(lrt.filtered_dge$table$PValue, method='fdr')

dim(lrt.filtered_dge$table[fdr<0.05,])
length(lrt.filtered_dge$table$PValue[lrt.filtered_dge$table$PValue<0.05])

#make a heatmap by getting the counts per million from each gene and turning them relative proportions (columns add up to 1)
sub1 <- colSums (cpm.filtered_dge)
sub2 <- matrix (rep(sub1,nrow (cpm.filtered_dge)), c (nrow (cpm.filtered_dge),ncol(cpm.filtered_dge)),byrow = TRUE)
sub3 <- cpm.filtered_dge / sub2

#subset this matrix to just get the genes from above
names1 <- row.names (sub3)
names2 <- row.names (top)
index2 <- names1 %in% names2
heatmap1 <- sub3[index2,]

##     par(mar = c(bottom, left, top, right))
# plot used for GenomeBC report
#play around with the options to make the plot fit what you like for options type ?heatmap.2

# order by site
# http://www.sthda.com/english/wiki/reordering-data-frame-columns-in-r
col_order <- c("C_Norw_10",   "C_Norw_1",    "C_Norw_3_B",  "C_Norw_5",    "C_Norw_6",    "C_Norw_7",    "W_Norw_10_B", "W_Norw_2",    "W_Norw_3_B",  "W_Norw_4_B", "W_Norw_6",    "W_Norw_7")
heatmap2 <- heatmap1[, col_order]

png("./plots/Norway_W_C_RNA_heatmap2.png", width = 3200, height = 2500)
heatmap.2 (heatmap1, scale = "row", trace = "none", Colv = FALSE, labCol = colnames(heatmap1), labRow = rownames(heatmap1), col='rainbow', dendrogram = "none", colsep =c(6), margins =c(20,70), sepcolor = "white", sepwidth = 0.1, cexRow=3, cexCol=6)
dev.off()

png("./plots/Norway_W_C_RNA_heatmap_legend2.png", width = 1000, height = 800)
heatmap.2 (heatmap1, scale = "row", trace = "none", notecex=3, notecol="black", tracecol="black", Colv = FALSE, labCol = colnames(heatmap1), labRow = rownames(heatmap1), col='rainbow', dendrogram = "none", colsep =c(5), margins =c(0.2, 0.2), sepcolor = "white", sepwidth = 0.1, cexRow=6, cexCol=6)
dev.off()

# write out file of genes that differ a lot

write.table(heatmap2, file = "./RNA_DER_May2023_W_C_Norway.txt", quote = FALSE, row.names=TRUE, col.names=TRUE, sep="\t")

#------------------------------
# Alaska

#read in the data
#mydata <- read.table("./RNA_site_tables/gene_names_expression_table1.txt", header=TRUE)

CDS <- mydata[,1]

exp_data <- mydata[,c(2:5, 30:37)]

 # [1] "C_Alas_16"   "C_Alas_1"    "C_Alas_20"   "C_Alas_5"    "C_Alex_11"
 # [6] "C_Alex_21"   "C_Alex_22"   "C_Alex_23"   "C_Alex_24"   "C_Alex_2"
# [11] "C_Alex_6"    "C_Alex_7"    "C_Norw_10"   "C_Norw_1"    "C_Norw_3_B"
# [16] "C_Norw_5"    "C_Norw_6"    "C_Norw_7"    "C_Seed_1"    "C_Seed_2"
# [21] "C_Seed_3"    "C_Swed_11"   "C_Swed_15"   "C_Swed_16"   "C_Swed_2"
# [26] "C_Swed_3"    "C_Swed_5"    "C_Swed_8"    "W_Alas_12"   "W_Alas_13"
# [31] "W_Alas_15"   "W_Alas_1"    "W_Alas_2"    "W_Alas_3"    "W_Alas_5"
# [36] "W_Alas_9"    "W_Alex_11"   "W_Alex_1"    "W_Alex_21"   "W_Alex_22"
# [41] "W_Alex_6"    "W_Norw_10_B" "W_Norw_2"    "W_Norw_3_B"  "W_Norw_4_B"
# [46] "W_Norw_6"    "W_Norw_7"    "W_Seed_1"    "W_Seed_2"    "W_Seed_3"
# [51] "W_Swed_10"   "W_Swed_13"   "W_Swed_14"   "W_Swed_16"   "W_Swed_17"
# [56] "W_Swed_2"    "W_Swed_5"    "W_Swed_6"


#extract and turn the column names into a factor
treat <- as.factor (sapply (strsplit(colnames(exp_data),split = "_"),"[[",1))
site <- as.factor (sapply (strsplit(colnames(exp_data),split = "_"),"[[",2))

treat_site <- paste0(as.character(site), as.character(treat))
treat_site <- as.factor(treat_site)

#make a DGElist object out of the data
dge <- DGEList (counts = exp_data, group = treat)

#calculate the normalization factors to adjust the effective library size relative to other libraries in the dataset (norm.factor that minimizes log fold change among samples)
dge <- calcNormFactors (dge)

#calculate the counts per million
cpm.dge <- cpm(dge)

#filter out the genes with less than 1 cpm in 6 or fewer libraries (a somewhat arbitrary choice). 
#Genes are usually dropped if they can't possibly be expressed in all the samples for any of the conditions.
filtered_dge <- dge[rowSums(cpm.dge > 1) >= 6,]
cpm.filtered_dge <- cpm.dge[rowSums(cpm.dge > 1) >= 6,]

#generate a multi-dimensional scaling plot
col_treat <- as.character (treat)
col_treat [col_treat == "C"] <- "blue"
col_treat [col_treat == "W"] <- "red"

col_site <- as.character (site)
col_site [col_treat == "Alex"] <- 15 # square
col_site [col_treat == "Norw"] <- 16 # circle
col_site [col_treat == "Alas"] <- 17 # triangle
col_site [col_treat == "Swed"] <- 18 # diamond

#col_treat [col_treat == "H"] <- "green"
#col_treat [col_treat == "L"] <- "yellow"

#plot the MDS graph
png("./plots/Alaska_W_C_RNA_MDS.png", width = 700, height = 500)
plotMDS (filtered_dge, col = col_treat)
dev.off()


#Set the model to use. This one includes the intercept, but other models can be specified that omit the intercept or that have more complex designs. See EdgeR manual for details.
design <- model.matrix (~treat)

#fit the common + tagwise dispersion models
filtered_dge <- estimateGLMCommonDisp (filtered_dge, design)
filtered_dge <- estimateGLMTagwiseDisp (filtered_dge, design)

#fit a GLM to the data using the tagwise dispersion estimate
glm.filtered_dge <- glmFit (filtered_dge, design, dispersion = filtered_dge$tagwise.dispersion)
lrt.filtered_dge <- glmLRT (glm.filtered_dge)

#get the topTags out of the model
#top <- topTags (lrt.filtered_dge, n = 1000)$table
#top <- topTags (lrt.filtered_dge, n = 200)$table
top <- topTags (lrt.filtered_dge, n = 500)$table # p-value all < 5e-02, v2

fdr<-p.adjust(lrt.filtered_dge$table$PValue, method='fdr')

dim(lrt.filtered_dge$table[fdr<0.05,])
length(lrt.filtered_dge$table$PValue[lrt.filtered_dge$table$PValue<0.05])

#make a heatmap by getting the counts per million from each gene and turning them relative proportions (columns add up to 1)
sub1 <- colSums (cpm.filtered_dge)
sub2 <- matrix (rep(sub1,nrow (cpm.filtered_dge)), c (nrow (cpm.filtered_dge),ncol(cpm.filtered_dge)),byrow = TRUE)
sub3 <- cpm.filtered_dge / sub2

#subset this matrix to just get the genes from above
names1 <- row.names (sub3)
names2 <- row.names (top)
index2 <- names1 %in% names2
heatmap1 <- sub3[index2,]

##     par(mar = c(bottom, left, top, right))
# plot used for GenomeBC report
#play around with the options to make the plot fit what you like for options type ?heatmap.2

# order by site
# http://www.sthda.com/english/wiki/reordering-data-frame-columns-in-r
col_order <- c("C_Alas_16",   "C_Alas_1",    "C_Alas_20",   "C_Alas_5",  "W_Alas_15",   "W_Alas_1",    "W_Alas_2",    "W_Alas_3",    "W_Alas_5",  "W_Alas_9",  "W_Alas_12",   "W_Alas_13")
heatmap2 <- heatmap1[, col_order]

png("./plots/Alaska_W_C_RNA_heatmap2.png", width = 3200, height = 2500)
heatmap.2 (heatmap1, scale = "row", trace = "none", Colv = FALSE, labCol = colnames(heatmap1), labRow = rownames(heatmap1), col='rainbow', dendrogram = "none", colsep =c(4), margins =c(20,70), sepcolor = "white", sepwidth = 0.1, cexRow=3, cexCol=6)
dev.off()

png("./plots/Alaska_W_C_RNA_heatmap_legend2.png", width = 1000, height = 800)
heatmap.2 (heatmap1, scale = "row", trace = "none", notecex=3, notecol="black", tracecol="black", Colv = FALSE, labCol = colnames(heatmap1), labRow = rownames(heatmap1), col='rainbow', dendrogram = "none", colsep =c(5), margins =c(0.2, 0.2), sepcolor = "white", sepwidth = 0.1, cexRow=6, cexCol=6)
dev.off()

# write out file of genes that differ a lot

write.table(heatmap2, file = "./RNA_DER_May2023_W_C_Alaska.txt", quote = FALSE, row.names=TRUE, col.names=TRUE, sep="\t")

#----------------------------------
# Sweden

#read in the data
#mydata <- read.table("./RNA_site_tables/gene_names_expression_table1.txt", header=TRUE)

CDS <- mydata[,1]
exp_data <- mydata[,c(23:29,52:59)]

 # [1] "C_Alas_16"   "C_Alas_1"    "C_Alas_20"   "C_Alas_5"    "C_Alex_11"
 # [6] "C_Alex_21"   "C_Alex_22"   "C_Alex_23"   "C_Alex_24"   "C_Alex_2"
# [11] "C_Alex_6"    "C_Alex_7"    "C_Norw_10"   "C_Norw_1"    "C_Norw_3_B"
# [16] "C_Norw_5"    "C_Norw_6"    "C_Norw_7"    "C_Seed_1"    "C_Seed_2"
# [21] "C_Seed_3"    "C_Swed_11"   "C_Swed_15"   "C_Swed_16"   "C_Swed_2"
# [26] "C_Swed_3"    "C_Swed_5"    "C_Swed_8"    "W_Alas_12"   "W_Alas_13"
# [31] "W_Alas_15"   "W_Alas_1"    "W_Alas_2"    "W_Alas_3"    "W_Alas_5"
# [36] "W_Alas_9"    "W_Alex_11"   "W_Alex_1"    "W_Alex_21"   "W_Alex_22"
# [41] "W_Alex_6"    "W_Norw_10_B" "W_Norw_2"    "W_Norw_3_B"  "W_Norw_4_B"
# [46] "W_Norw_6"    "W_Norw_7"    "W_Seed_1"    "W_Seed_2"    "W_Seed_3"
# [51] "W_Swed_10"   "W_Swed_13"   "W_Swed_14"   "W_Swed_16"   "W_Swed_17"
# [56] "W_Swed_2"    "W_Swed_5"    "W_Swed_6"

#extract and turn the column names into a factor
treat <- as.factor (sapply (strsplit(colnames(exp_data),split = "_"),"[[",1))
site <- as.factor (sapply (strsplit(colnames(exp_data),split = "_"),"[[",2))

treat_site <- paste0(as.character(site), as.character(treat))
treat_site <- as.factor(treat_site)

#make a DGElist object out of the data
dge <- DGEList (counts = exp_data, group = treat)

#calculate the normalization factors to adjust the effective library size relative to other libraries in the dataset (norm.factor that minimizes log fold change among samples)
dge <- calcNormFactors (dge)

#calculate the counts per million
cpm.dge <- cpm(dge)

#filter out the genes with less than 1 cpm in 6 or fewer libraries (a somewhat arbitrary choice). 
#Genes are usually dropped if they can't possibly be expressed in all the samples for any of the conditions.
filtered_dge <- dge[rowSums(cpm.dge > 1) >= 6,]
cpm.filtered_dge <- cpm.dge[rowSums(cpm.dge > 1) >= 6,]

#generate a multi-dimensional scaling plot
col_treat <- as.character (treat)
col_treat [col_treat == "C"] <- "blue"
col_treat [col_treat == "W"] <- "red"

col_site <- as.character (site)
col_site [col_treat == "Alex"] <- 15 # square
col_site [col_treat == "Norw"] <- 16 # circle
col_site [col_treat == "Alas"] <- 17 # triangle
col_site [col_treat == "Swed"] <- 18 # diamond

#col_treat [col_treat == "H"] <- "green"
#col_treat [col_treat == "L"] <- "yellow"

#plot the MDS graph
png("./plots/Sweden_W_C_RNA_MDS.png", width = 700, height = 500)
plotMDS (filtered_dge, col = col_treat)
dev.off()

png("./plots/Sweden_Site_RNA_MDS.png", width = 700, height = 500)
plotMDS (filtered_dge, col = col_treat, pch=col_site)
dev.off()

#Set the model to use. This one includes the intercept, but other models can be specified that omit the intercept or that have more complex designs. See EdgeR manual for details.
design <- model.matrix (~treat)

#fit the common + tagwise dispersion models
filtered_dge <- estimateGLMCommonDisp (filtered_dge, design)
filtered_dge <- estimateGLMTagwiseDisp (filtered_dge, design)

#fit a GLM to the data using the tagwise dispersion estimate
glm.filtered_dge <- glmFit (filtered_dge, design, dispersion = filtered_dge$tagwise.dispersion)
lrt.filtered_dge <- glmLRT (glm.filtered_dge)

#get the topTags out of the model
#top <- topTags (lrt.filtered_dge, n = 1000)$table
#top <- topTags (lrt.filtered_dge, n = 200)$table
top <- topTags (lrt.filtered_dge, n = 500)$table # p-value all < 5e-02, v2

fdr<-p.adjust(lrt.filtered_dge$table$PValue, method='fdr')

dim(lrt.filtered_dge$table[fdr<0.05,])
length(lrt.filtered_dge$table$PValue[lrt.filtered_dge$table$PValue<0.05])

#make a heatmap by getting the counts per million from each gene and turning them relative proportions (columns add up to 1)
sub1 <- colSums (cpm.filtered_dge)
sub2 <- matrix (rep(sub1,nrow (cpm.filtered_dge)), c (nrow (cpm.filtered_dge),ncol(cpm.filtered_dge)),byrow = TRUE)
sub3 <- cpm.filtered_dge / sub2

#subset this matrix to just get the genes from above
names1 <- row.names (sub3)
names2 <- row.names (top)
index2 <- names1 %in% names2
heatmap1 <- sub3[index2,]

##     par(mar = c(bottom, left, top, right))
# plot used for GenomeBC report
#play around with the options to make the plot fit what you like for options type ?heatmap.2

# order by site
# http://www.sthda.com/english/wiki/reordering-data-frame-columns-in-r
col_order <- c("C_Swed_15",   "C_Swed_16",   "C_Swed_2", "C_Swed_11",   "C_Swed_3",    "C_Swed_5",    "C_Swed_8",   "W_Swed_10",   "W_Swed_13",   "W_Swed_14",   "W_Swed_16",   "W_Swed_17",  "W_Swed_2",    "W_Swed_5",    "W_Swed_6", "C_Seed_1",    "C_Seed_2", "C_Seed_3", "W_Seed_1",    "W_Seed_2",    "W_Seed_3")
heatmap2 <- heatmap1[, col_order]

png("./plots/Sweden_W_C_RNA_heatmap2.png", width = 3200, height = 2500)
heatmap.2 (heatmap1, scale = "row", trace = "none", Colv = FALSE, labCol = colnames(heatmap1), labRow = rownames(heatmap1), col='rainbow', dendrogram = "none", colsep =c(7), margins =c(20,70), sepcolor = "white", sepwidth = 0.1, cexRow=3, cexCol=6)
dev.off()

png("./plots/Sweden_W_C_RNA_heatmap_legend2.png", width = 1000, height = 800)
heatmap.2 (heatmap1, scale = "row", trace = "none", notecex=3, notecol="black", tracecol="black", Colv = FALSE, labCol = colnames(heatmap1), labRow = rownames(heatmap1), col='rainbow', dendrogram = "none", colsep =c(5), margins =c(0.2, 0.2), sepcolor = "white", sepwidth = 0.1, cexRow=6, cexCol=6)
dev.off()

subdat1<- exp_data["Do1_01_a00001G00358V1.1",]

png("W_C_Sweden_RNA_Do1_01_a00001G00358V1.1.png", width = 1000, height = 700)
stripchart(log10(as.numeric (subdat1)) ~ treat, vertical = T, xlim = c (0.5,2.5), method = "jitter", ylab = "", main="Parent RNAseq Do1_01_a00001G00358V1", pch=16, cex=4, cex.main=1, cex.axis=2, cex.lab=2)
dev.off()

# write out file of genes that differ a lot

write.table(heatmap2, file = "./RNA_DER_May2023_W_C_Sweden.txt", quote = FALSE, row.names=TRUE, col.names=TRUE, sep="\t")

#----------------------------------
# Seedlings

#read in the data
#mydata <- read.table("./RNA_site_tables/gene_names_expression_table1.txt", header=TRUE)

CDS <- mydata[,1]
exp_data <- mydata[,c(20:22,49:51)]

 # [1] "C_Alas_16"   "C_Alas_1"    "C_Alas_20"   "C_Alas_5"    "C_Alex_11"
 # [6] "C_Alex_21"   "C_Alex_22"   "C_Alex_23"   "C_Alex_24"   "C_Alex_2"
# [11] "C_Alex_6"    "C_Alex_7"    "C_Norw_10"   "C_Norw_1"    "C_Norw_3_B"
# [16] "C_Norw_5"    "C_Norw_6"    "C_Norw_7"    "C_Seed_1"    "C_Seed_2"
# [21] "C_Seed_3"    "C_Swed_11"   "C_Swed_15"   "C_Swed_16"   "C_Swed_2"
# [26] "C_Swed_3"    "C_Swed_5"    "C_Swed_8"    "W_Alas_12"   "W_Alas_13"
# [31] "W_Alas_15"   "W_Alas_1"    "W_Alas_2"    "W_Alas_3"    "W_Alas_5"
# [36] "W_Alas_9"    "W_Alex_11"   "W_Alex_1"    "W_Alex_21"   "W_Alex_22"
# [41] "W_Alex_6"    "W_Norw_10_B" "W_Norw_2"    "W_Norw_3_B"  "W_Norw_4_B"
# [46] "W_Norw_6"    "W_Norw_7"    "W_Seed_1"    "W_Seed_2"    "W_Seed_3"
# [51] "W_Swed_10"   "W_Swed_13"   "W_Swed_14"   "W_Swed_16"   "W_Swed_17"
# [56] "W_Swed_2"    "W_Swed_5"    "W_Swed_6"

#extract and turn the column names into a factor
treat <- as.factor (sapply (strsplit(colnames(exp_data),split = "_"),"[[",1))
site <- as.factor (sapply (strsplit(colnames(exp_data),split = "_"),"[[",2))

treat_site <- paste0(as.character(site), as.character(treat))
treat_site <- as.factor(treat_site)

#make a DGElist object out of the data
dge <- DGEList (counts = exp_data, group = treat)

#calculate the normalization factors to adjust the effective library size relative to other libraries in the dataset (norm.factor that minimizes log fold change among samples)
dge <- calcNormFactors (dge)

#calculate the counts per million
cpm.dge <- cpm(dge)

#filter out the genes with less than 1 cpm in 6 or fewer libraries (a somewhat arbitrary choice). 
#Genes are usually dropped if they can't possibly be expressed in all the samples for any of the conditions.
filtered_dge <- dge[rowSums(cpm.dge > 1) >= 6,]
cpm.filtered_dge <- cpm.dge[rowSums(cpm.dge > 1) >= 6,]

#generate a multi-dimensional scaling plot
col_treat <- as.character (treat)
col_treat [col_treat == "C"] <- "blue"
col_treat [col_treat == "W"] <- "red"

col_site <- as.character (site)
col_site [col_treat == "Alex"] <- 15 # square
col_site [col_treat == "Norw"] <- 16 # circle
col_site [col_treat == "Alas"] <- 17 # triangle
col_site [col_treat == "Swed"] <- 18 # diamond

#col_treat [col_treat == "H"] <- "green"
#col_treat [col_treat == "L"] <- "yellow"

#plot the MDS graph
png("./plots/Seed_W_C_RNA_MDS.png", width = 700, height = 500)
plotMDS (filtered_dge, col = col_treat)
dev.off()

#Set the model to use. This one includes the intercept, but other models can be specified that omit the intercept or that have more complex designs. See EdgeR manual for details.
design <- model.matrix (~treat)

#fit the common + tagwise dispersion models
filtered_dge <- estimateGLMCommonDisp (filtered_dge, design)
filtered_dge <- estimateGLMTagwiseDisp (filtered_dge, design)

#fit a GLM to the data using the tagwise dispersion estimate
glm.filtered_dge <- glmFit (filtered_dge, design, dispersion = filtered_dge$tagwise.dispersion)
lrt.filtered_dge <- glmLRT (glm.filtered_dge)

#get the topTags out of the model
#top <- topTags (lrt.filtered_dge, n = 1000)$table
#top <- topTags (lrt.filtered_dge, n = 200)$table
top <- topTags (lrt.filtered_dge, n = 500)$table # p-value all < 5e-02, v2

fdr<-p.adjust(lrt.filtered_dge$table$PValue, method='fdr')

dim(lrt.filtered_dge$table[fdr<0.05,])
length(lrt.filtered_dge$table$PValue[lrt.filtered_dge$table$PValue<0.05])

#make a heatmap by getting the counts per million from each gene and turning them relative proportions (columns add up to 1)
sub1 <- colSums (cpm.filtered_dge)
sub2 <- matrix (rep(sub1,nrow (cpm.filtered_dge)), c (nrow (cpm.filtered_dge),ncol(cpm.filtered_dge)),byrow = TRUE)
sub3 <- cpm.filtered_dge / sub2

#subset this matrix to just get the genes from above
names1 <- row.names (sub3)
names2 <- row.names (top)
index2 <- names1 %in% names2
heatmap1 <- sub3[index2,]

##     par(mar = c(bottom, left, top, right))
# plot used for GenomeBC report
#play around with the options to make the plot fit what you like for options type ?heatmap.2

# order by site
# http://www.sthda.com/english/wiki/reordering-data-frame-columns-in-r
col_order <- c("C_Swed_15",   "C_Swed_16",   "C_Swed_2", "C_Swed_11",   "C_Swed_3",    "C_Swed_5",    "C_Swed_8",   "W_Swed_10",   "W_Swed_13",   "W_Swed_14",   "W_Swed_16",   "W_Swed_17",  "W_Swed_2",    "W_Swed_5",    "W_Swed_6", "C_Seed_1",    "C_Seed_2", "C_Seed_3", "W_Seed_1",    "W_Seed_2",    "W_Seed_3")
heatmap2 <- heatmap1[, col_order]

png("./plots/Seed_W_C_RNA_heatmap2.png", width = 3200, height = 2500)
heatmap.2 (heatmap1, scale = "row", trace = "none", Colv = FALSE, labCol = colnames(heatmap1), labRow = rownames(heatmap1), col='rainbow', dendrogram = "none", colsep =c(5), margins =c(20,70), sepcolor = "white", sepwidth = 0.1, cexRow=2, cexCol=2)
dev.off()

png("./plots/Seed_W_C_RNA_heatmap_legend2.png", width = 1000, height = 800)
heatmap.2 (heatmap1, scale = "row", trace = "none", notecex=3, notecol="black", tracecol="black", Colv = FALSE, labCol = colnames(heatmap1), labRow = rownames(heatmap1), col='rainbow', dendrogram = "none", colsep =c(5), margins =c(0.2, 0.2), sepcolor = "white", sepwidth = 0.1, cexRow=6, cexCol=6)
dev.off()

subdat1seed <- exp_data["Do1_01_a00001G00358V1.1",]

png("W_C_RNA_Seed_Do1_01_a00001G00358V1.1.png", width = 1000, height = 700)
stripchart (log10(as.numeric (subdat1seed)) ~ treat, vertical = T, xlim = c (0.5,2.5), method = "jitter", ylab = "",  main="Seedling RNAseq Do1_01_a00001G00358V1", pch=16, cex=4, cex.main=1, cex.axis=2, cex.lab=2)
dev.off()


# write out file of genes that differ a lot

write.table(heatmap2, file = "./RNA_DER_May2023_W_C_Seedlings.txt", quote = FALSE, row.names=TRUE, col.names=TRUE, sep="\t")


#-------------------------
#order genes based on most to least expressed in first individual
#heatmap1_df <- as.data.frame(heatmap1)
#heatmap1_df$C_L11 <- as.numeric(as.character(heatmap1_df$C_L11))
#heatmap_ordered_df <- heatmap1_df[order(heatmap1_df$C_L11),]
#heatmap_ordered <- as.matrix(heatmap_ordered_df)

#png("W_C_RNA_heatmap_ordered.png", width = 1500, height = 700)
#heatmap.2 (heatmap_ordered, trace = "none", scale = "row", Rowv = FALSE, Colv = FALSE, labCol = colnames(heatmap_ordered), cexCol = 1, dendrogram = "none", labRow = rownames(heatmap_ordered), colsep =c(5), sepcolor = "white", sepwidth = 0.03)
#dev.off()

###################################
# Sweden differential expression - specific genes plotting 
#Do1_01_a00001G00358V1.1

#plot the individual expression values from a single gene:

subdat1<- exp_data["Do1_01_a00001G00358V1.1",]
png("W_C_Sweden_RNA_Do1_01_a00001G00358V1.1.png", width = 1000, height = 700)
stripchart(log10(as.numeric (subdat1)) ~ treat, vertical = T, xlim = c (0.5,2.5), method = "jitter", ylab = "", main="Parent RNAseq Do1_01_a00001G00358V1", pch=16, cex=4, cex.main=1, cex.axis=2, cex.lab=2)
dev.off()

# Parents
                       # C_Swed_11 C_Swed_15 C_Swed_16 C_Swed_2 C_Swed_3
# Do1_01_a00001G00358V1.1       800       830       791      712     1313
                        # C_Swed_5 C_Swed_8 W_Swed_10 W_Swed_13 W_Swed_14
# Do1_01_a00001G00358V1.1      783      415      1240       204       153
                        # W_Swed_16 W_Swed_17 W_Swed_2 W_Swed_5 W_Swed_6
# Do1_01_a00001G00358V1.1       172       398      105      368     2143

#-----------
subdat1seed <- exp_data["Do1_01_a00001G00358V1.1",]

png("W_C_RNA_Seed_Do1_01_a00001G00358V1.1.png", width = 1000, height = 700)
stripchart (log10(as.numeric (subdat1seed)) ~ treat, vertical = T, xlim = c (0.5,2.5), method = "jitter", ylab = "",  main="Seedling RNAseq Do1_01_a00001G00358V1", pch=16, cex=4, cex.main=1, cex.axis=2, cex.lab=2)
dev.off()


#Seedlings
                        # C_Seed_1 C_Seed_2 C_Seed_3 W_Seed_1 W_Seed_2 W_Seed_3
# Do1_01_a00001G00358V1.1      284      151       45      762     2112     1527


#----------------------------------
# Do1_05_a00001G00328V1.1

#plot the individual expression values from a single gene:
subdat1 <- exp_data["Do1_05_a00001G00328V1.1",]

png("W_C_RNA_Do1_05_a00001G00328V1.1.png", width = 700, height = 500)
stripchart (log10(as.numeric (subdat1)) ~ treat, vertical = T, xlim = c (0.5,2.5), method = "jitter", ylab = "Log10 Expression count", pch=16)
dev.off()
#----------------------------------
#Do1_02_a00003G00097V1.1
#plot the individual expression values from a single gene:
subdat1 <- exp_data["Do1_02_a00003G00097V1.1",]

png("W_C_RNA_Do1_02_a00003G00097V1.1.png", width = 700, height = 500)
stripchart (log10(as.numeric (subdat1)) ~ treat, vertical = T, xlim = c (0.5,2.5), method = "jitter", ylab = "Log10 Expression count", pch=16)
dev.off()

#----------------------------------
#Do1_06_a00002G00149V1.1
#plot the individual expression values from a single gene:
subdat1 <- exp_data["Do1_06_a00002G00149V1.1",]

png("W_C_RNA_Do1_06_a00002G00149V1.1.png", width = 700, height = 500)
stripchart (log10(as.numeric (subdat1)) ~ treat, vertical = T, xlim = c (0.5,2.5), method = "jitter", ylab = "Log10 Expression count", pch=16)
dev.off()
#----------------------------------
#Do1_01_a00001G02208V1.1
#plot the individual expression values from a single gene:
subdat1 <- exp_data["Do1_01_a00001G02208V1.1",]

png("W_C_RNA_Do1_01_a00001G02208V1.1.png", width = 700, height = 500)
stripchart (log10(as.numeric (subdat1)) ~ treat, vertical = T, xlim = c (0.5,2.5), method = "jitter", ylab = "Log10 Expression count", pch=16)
dev.off()

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


