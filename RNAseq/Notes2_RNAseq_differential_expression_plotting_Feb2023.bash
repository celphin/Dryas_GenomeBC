########################################
# Dryas differential gene expression data
# Determining differentially expressed genes, plotting results
# Sept 2023 - Sweden, Feb 2023 other sites
#####################################

###########################################################
# try differential expression analysis - plotting
# EdgeR

# https://github.com/owensgl/biol525D/tree/master/Topic_6
# https://bioconductor.org/packages/release/bioc/html/edgeR.html
# http://combine-australia.github.io/RNAseq-R/slides/RNASeq_filtering_qc.pdf 

cd /home/celphin/projects/rpp-rieseber/celphin/Dryas/RNAseq_analysis/

# get R on server
# open R
module load StdEnv/2020
module load r/4.2.1
module load gdal
module load udunits
module load python
export R_LIBS_USER=/home/celphin/R/x86_64-pc-linux-gnu-library/4.2.1/

R

#install the package edgeR
# if (!require("BiocManager", quietly = TRUE))
    # install.packages("BiocManager")
# BiocManager::install(version = "3.15")

# BiocManager::install("edgeR")
# install.packages('gplots')

#paste it in here (i.e. replace my path with yours):
setwd ("/home/celphin/projects/rpp-rieseber/celphin/Dryas/RNAseq_analysis/")

#find your working directory:
getwd()

#load the libraries you will need 
library ("edgeR")
library ("gplots")

#read in the data
mydata <- read.table("./RNA_site_tables/gene_names_expression_table1.txt", header=TRUE)

CDS <- mydata[,1]
exp_data <- mydata[,c(2:length(mydata))]

#extract and turn the column names into a factor
treat <- as.factor (sapply (strsplit(colnames(exp_data),split = "_"),"[[",1))
site <- as.factor (sapply (strsplit(colnames(exp_data),split = "_"),"[[",2))

treat_site <- paste0(as.character(site), as.character(treat))
treat_site <- as.factor(treat_site)

#make a DGElist object out of the data
list1 <- DGEList (counts = exp_data, group = treat)

#calculate the normalization factors to adjust the effective library size relative to other libraries in the dataset (norm.factor that minimizes log fold change among samples)
list1 <- calcNormFactors (list1)

#calculate the counts per million
cpm.list1 <- cpm(list1)

#filter out the genes with less than 1 cpm in 6 or fewer libraries (a somewhat arbitrary choice). 
#Genes are usually dropped if they can't possibly be expressed in all the samples for any of the conditions.
list2 <- list1[rowSums(cpm.list1 > 1) >= 6,]
cpm.list2 <- cpm.list1[rowSums(cpm.list1 > 1) >= 6,]

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
jpeg("./plots/W_C_RNA_MDS.jpg", width = 700, height = 500)
plotMDS (list2, col = col_treat)
dev.off()

jpeg("./plots/Site_RNA_MDS.jpg", width = 700, height = 500)
plotMDS (list2, col = col_treat, pch=col_site)
dev.off()

#Set the model to use. This one includes the intercept, but other models can be specified that omit the intercept or that have more complex designs. See EdgeR manual for details.
design <- model.matrix (~treat)

#fit the common + tagwise dispersion models
list2 <- estimateGLMCommonDisp (list2, design)
list2 <- estimateGLMTagwiseDisp (list2, design)

#fit a GLM to the data using the tagwise dispersion estimate
glm.list2 <- glmFit (list2, design, dispersion = list2$tagwise.dispersion)
lrt.list2 <- glmLRT (glm.list2)

#get the topTags out of the model
#top <- topTags (lrt.list2, n = 1000)$table
#top <- topTags (lrt.list2, n = 200)$table
top <- topTags (lrt.list2, n = 40)$table # p-value all < 5e-02 

fdr<-p.adjust(lrt.list2$table$PValue, method='fdr')

dim(lrt.list2$table[fdr<0.05,])
length(lrt.list2$table$PValue[lrt.list2$table$PValue<0.05])

#make a heatmap by getting the counts per million from each gene and turning them relative proportions (columns add up to 1)
sub1 <- colSums (cpm.list2)
sub2 <- matrix (rep(sub1,nrow (cpm.list2)), c (nrow (cpm.list2),ncol(cpm.list2)),byrow = TRUE)
sub3 <- cpm.list2 / sub2

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

jpeg("./plots/W_C_RNA_heatmap.jpg", width = 3200, height = 1000)
heatmap.2 (heatmap2, scale = "row", trace = "none", Colv = FALSE, labCol = colnames(heatmap2), labRow = rownames(heatmap2), col='rainbow', dendrogram = "none", colsep =c(5), margins =c(20,70), sepcolor = "white", sepwidth = 0.1, cexRow=2, cexCol=2)
dev.off()

jpeg("./plots/W_C_RNA_heatmap_legend.jpg", width = 1000, height = 800)
heatmap.2 (heatmap2, scale = "row", trace = "none", notecex=3, notecol="black", tracecol="black", Colv = FALSE, labCol = colnames(heatmap2), labRow = rownames(heatmap2), col='rainbow', dendrogram = "none", colsep =c(5), margins =c(0.2, 0.2), sepcolor = "white", sepwidth = 0.1, cexRow=6, cexCol=6)
dev.off()

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
list1 <- DGEList (counts = exp_data, group = treat)

#calculate the normalization factors to adjust the effective library size relative to other libraries in the dataset (norm.factor that minimizes log fold change among samples)
list1 <- calcNormFactors (list1)

#calculate the counts per million
cpm.list1 <- cpm(list1)

#filter out the genes with less than 1 cpm in 6 or fewer libraries (a somewhat arbitrary choice). 
#Genes are usually dropped if they can't possibly be expressed in all the samples for any of the conditions.
list2 <- list1[rowSums(cpm.list1 > 1) >= 6,]
cpm.list2 <- cpm.list1[rowSums(cpm.list1 > 1) >= 6,]

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
jpeg("./plots/Alex_W_C_RNA_MDS.jpg", width = 700, height = 500)
plotMDS (list2, col = col_treat)
dev.off()

#Set the model to use. This one includes the intercept, but other models can be specified that omit the intercept or that have more complex designs. See EdgeR manual for details.
design <- model.matrix (~treat)

#fit the common + tagwise dispersion models
list2 <- estimateGLMCommonDisp (list2, design)
list2 <- estimateGLMTagwiseDisp (list2, design)

#fit a GLM to the data using the tagwise dispersion estimate
glm.list2 <- glmFit (list2, design, dispersion = list2$tagwise.dispersion)
lrt.list2 <- glmLRT (glm.list2)

#get the topTags out of the model
#top <- topTags (lrt.list2, n = 1000)$table
#top <- topTags (lrt.list2, n = 200)$table
top <- topTags (lrt.list2, n = 40)$table # p-value all < 5e-02, v2

fdr<-p.adjust(lrt.list2$table$PValue, method='fdr')

dim(lrt.list2$table[fdr<0.05,])
length(lrt.list2$table$PValue[lrt.list2$table$PValue<0.05])

#make a heatmap by getting the counts per million from each gene and turning them relative proportions (columns add up to 1)
sub1 <- colSums (cpm.list2)
sub2 <- matrix (rep(sub1,nrow (cpm.list2)), c (nrow (cpm.list2),ncol(cpm.list2)),byrow = TRUE)
sub3 <- cpm.list2 / sub2

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
#col_order <- c("C_Alex_11",  "C_Alex_21",   "C_Alex_22",   "C_Alex_23",   "C_Alex_24",   "C_Alex_2",  "C_Alex_6",    "C_Alex_7",    "W_Alex_11",   "W_Alex_1",    "W_Alex_21",   "W_Alex_22",  "W_Alex_6")
#heatmap2 <- heatmap1[, col_order]

jpeg("./plots/Alex_W_C_RNA_heatmap2.jpg", width = 3200, height = 2500)
heatmap.2 (heatmap1, scale = "row", trace = "none", Colv = FALSE, labCol = colnames(heatmap1), labRow = rownames(heatmap1), col='rainbow', dendrogram = "none", colsep =c(8), margins =c(20,70), sepcolor = "white", sepwidth = 0.1, cexRow=3, cexCol=6)
dev.off()

jpeg("./plots/Alex_W_C_RNA_heatmap_legend2.jpg", width = 1000, height = 800)
heatmap.2 (heatmap1, scale = "row", trace = "none", notecex=3, notecol="black", tracecol="black", Colv = FALSE, labCol = colnames(heatmap1), labRow = rownames(heatmap1), col='rainbow', dendrogram = "none", colsep =c(5), margins =c(0.2, 0.2), sepcolor = "white", sepwidth = 0.1, cexRow=6, cexCol=6)
dev.off()
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
list1 <- DGEList (counts = exp_data, group = treat)

#calculate the normalization factors to adjust the effective library size relative to other libraries in the dataset (norm.factor that minimizes log fold change among samples)
list1 <- calcNormFactors (list1)

#calculate the counts per million
cpm.list1 <- cpm(list1)

#filter out the genes with less than 1 cpm in 6 or fewer libraries (a somewhat arbitrary choice). 
#Genes are usually dropped if they can't possibly be expressed in all the samples for any of the conditions.
list2 <- list1[rowSums(cpm.list1 > 1) >= 6,]
cpm.list2 <- cpm.list1[rowSums(cpm.list1 > 1) >= 6,]

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
jpeg("./plots/Norway_W_C_RNA_MDS.jpg", width = 700, height = 500)
plotMDS (list2, col = col_treat)
dev.off()


#Set the model to use. This one includes the intercept, but other models can be specified that omit the intercept or that have more complex designs. See EdgeR manual for details.
design <- model.matrix (~treat)

#fit the common + tagwise dispersion models
list2 <- estimateGLMCommonDisp (list2, design)
list2 <- estimateGLMTagwiseDisp (list2, design)

#fit a GLM to the data using the tagwise dispersion estimate
glm.list2 <- glmFit (list2, design, dispersion = list2$tagwise.dispersion)
lrt.list2 <- glmLRT (glm.list2)

#get the topTags out of the model
#top <- topTags (lrt.list2, n = 1000)$table
#top <- topTags (lrt.list2, n = 200)$table
top <- topTags (lrt.list2, n = 40)$table # p-value all < 5e-02, v2

fdr<-p.adjust(lrt.list2$table$PValue, method='fdr')

dim(lrt.list2$table[fdr<0.05,])
length(lrt.list2$table$PValue[lrt.list2$table$PValue<0.05])

#make a heatmap by getting the counts per million from each gene and turning them relative proportions (columns add up to 1)
sub1 <- colSums (cpm.list2)
sub2 <- matrix (rep(sub1,nrow (cpm.list2)), c (nrow (cpm.list2),ncol(cpm.list2)),byrow = TRUE)
sub3 <- cpm.list2 / sub2

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
#col_order <- c("C_Norw_10",   "C_Norw_1",    "C_Norw_3_B",  "C_Norw_5",    "C_Norw_6",    "C_Norw_7",    "W_Norw_10_B", "W_Norw_2",    "W_Norw_3_B",  "W_Norw_4_B", "W_Norw_6",    "W_Norw_7")
#heatmap2 <- heatmap1[, col_order]

jpeg("./plots/Norway_W_C_RNA_heatmap2.jpg", width = 3200, height = 2500)
heatmap.2 (heatmap1, scale = "row", trace = "none", Colv = FALSE, labCol = colnames(heatmap1), labRow = rownames(heatmap1), col='rainbow', dendrogram = "none", colsep =c(6), margins =c(20,70), sepcolor = "white", sepwidth = 0.1, cexRow=3, cexCol=6)
dev.off()

jpeg("./plots/Norway_W_C_RNA_heatmap_legend2.jpg", width = 1000, height = 800)
heatmap.2 (heatmap1, scale = "row", trace = "none", notecex=3, notecol="black", tracecol="black", Colv = FALSE, labCol = colnames(heatmap1), labRow = rownames(heatmap1), col='rainbow', dendrogram = "none", colsep =c(5), margins =c(0.2, 0.2), sepcolor = "white", sepwidth = 0.1, cexRow=6, cexCol=6)
dev.off()
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
list1 <- DGEList (counts = exp_data, group = treat)

#calculate the normalization factors to adjust the effective library size relative to other libraries in the dataset (norm.factor that minimizes log fold change among samples)
list1 <- calcNormFactors (list1)

#calculate the counts per million
cpm.list1 <- cpm(list1)

#filter out the genes with less than 1 cpm in 6 or fewer libraries (a somewhat arbitrary choice). 
#Genes are usually dropped if they can't possibly be expressed in all the samples for any of the conditions.
list2 <- list1[rowSums(cpm.list1 > 1) >= 6,]
cpm.list2 <- cpm.list1[rowSums(cpm.list1 > 1) >= 6,]

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
jpeg("./plots/Alaska_W_C_RNA_MDS.jpg", width = 700, height = 500)
plotMDS (list2, col = col_treat)
dev.off()


#Set the model to use. This one includes the intercept, but other models can be specified that omit the intercept or that have more complex designs. See EdgeR manual for details.
design <- model.matrix (~treat)

#fit the common + tagwise dispersion models
list2 <- estimateGLMCommonDisp (list2, design)
list2 <- estimateGLMTagwiseDisp (list2, design)

#fit a GLM to the data using the tagwise dispersion estimate
glm.list2 <- glmFit (list2, design, dispersion = list2$tagwise.dispersion)
lrt.list2 <- glmLRT (glm.list2)

#get the topTags out of the model
#top <- topTags (lrt.list2, n = 1000)$table
#top <- topTags (lrt.list2, n = 200)$table
top <- topTags (lrt.list2, n = 40)$table # p-value all < 5e-02, v2

fdr<-p.adjust(lrt.list2$table$PValue, method='fdr')

dim(lrt.list2$table[fdr<0.05,])
length(lrt.list2$table$PValue[lrt.list2$table$PValue<0.05])

#make a heatmap by getting the counts per million from each gene and turning them relative proportions (columns add up to 1)
sub1 <- colSums (cpm.list2)
sub2 <- matrix (rep(sub1,nrow (cpm.list2)), c (nrow (cpm.list2),ncol(cpm.list2)),byrow = TRUE)
sub3 <- cpm.list2 / sub2

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
#col_order <- c("C_Alas_16",   "C_Alas_1",    "C_Alas_20",   "C_Alas_5",  "W_Alas_15",   "W_Alas_1",    "W_Alas_2",    "W_Alas_3",    "W_Alas_5",  "W_Alas_9",  "W_Alas_12",   "W_Alas_13")
#heatmap2 <- heatmap1[, col_order]

jpeg("./plots/Alaska_W_C_RNA_heatmap2.jpg", width = 3200, height = 2500)
heatmap.2 (heatmap1, scale = "row", trace = "none", Colv = FALSE, labCol = colnames(heatmap1), labRow = rownames(heatmap1), col='rainbow', dendrogram = "none", colsep =c(4), margins =c(20,70), sepcolor = "white", sepwidth = 0.1, cexRow=3, cexCol=6)
dev.off()

jpeg("./plots/Alaska_W_C_RNA_heatmap_legend2.jpg", width = 1000, height = 800)
heatmap.2 (heatmap1, scale = "row", trace = "none", notecex=3, notecol="black", tracecol="black", Colv = FALSE, labCol = colnames(heatmap1), labRow = rownames(heatmap1), col='rainbow', dendrogram = "none", colsep =c(5), margins =c(0.2, 0.2), sepcolor = "white", sepwidth = 0.1, cexRow=6, cexCol=6)
dev.off()
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
list1 <- DGEList (counts = exp_data, group = treat)

#calculate the normalization factors to adjust the effective library size relative to other libraries in the dataset (norm.factor that minimizes log fold change among samples)
list1 <- calcNormFactors (list1)

#calculate the counts per million
cpm.list1 <- cpm(list1)

#filter out the genes with less than 1 cpm in 6 or fewer libraries (a somewhat arbitrary choice). 
#Genes are usually dropped if they can't possibly be expressed in all the samples for any of the conditions.
list2 <- list1[rowSums(cpm.list1 > 1) >= 6,]
cpm.list2 <- cpm.list1[rowSums(cpm.list1 > 1) >= 6,]

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
jpeg("./plots/Sweden_W_C_RNA_MDS.jpg", width = 700, height = 500)
plotMDS (list2, col = col_treat)
dev.off()

jpeg("./plots/Sweden_Site_RNA_MDS.jpg", width = 700, height = 500)
plotMDS (list2, col = col_treat, pch=col_site)
dev.off()

#Set the model to use. This one includes the intercept, but other models can be specified that omit the intercept or that have more complex designs. See EdgeR manual for details.
design <- model.matrix (~treat)

#fit the common + tagwise dispersion models
list2 <- estimateGLMCommonDisp (list2, design)
list2 <- estimateGLMTagwiseDisp (list2, design)

#fit a GLM to the data using the tagwise dispersion estimate
glm.list2 <- glmFit (list2, design, dispersion = list2$tagwise.dispersion)
lrt.list2 <- glmLRT (glm.list2)

#get the topTags out of the model
#top <- topTags (lrt.list2, n = 1000)$table
#top <- topTags (lrt.list2, n = 200)$table
top <- topTags (lrt.list2, n = 40)$table # p-value all < 5e-02, v2

fdr<-p.adjust(lrt.list2$table$PValue, method='fdr')

dim(lrt.list2$table[fdr<0.05,])
length(lrt.list2$table$PValue[lrt.list2$table$PValue<0.05])

#make a heatmap by getting the counts per million from each gene and turning them relative proportions (columns add up to 1)
sub1 <- colSums (cpm.list2)
sub2 <- matrix (rep(sub1,nrow (cpm.list2)), c (nrow (cpm.list2),ncol(cpm.list2)),byrow = TRUE)
sub3 <- cpm.list2 / sub2

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
#col_order <- c("C_Swed_15",   "C_Swed_16",   "C_Swed_2", "C_Swed_11",   "C_Swed_3",    "C_Swed_5",    "C_Swed_8",   "W_Swed_10",   "W_Swed_13",   "W_Swed_14",   "W_Swed_16",   "W_Swed_17",  "W_Swed_2",    "W_Swed_5",    "W_Swed_6", 
#"C_Seed_1",    "C_Seed_2", "C_Seed_3", "W_Seed_1",    "W_Seed_2",    "W_Seed_3")
#heatmap2 <- heatmap1[, col_order]

jpeg("./plots/Sweden_W_C_RNA_heatmap2.jpg", width = 3200, height = 2500)
heatmap.2 (heatmap1, scale = "row", trace = "none", Colv = FALSE, labCol = colnames(heatmap1), labRow = rownames(heatmap1), col='rainbow', dendrogram = "none", colsep =c(7), margins =c(20,70), sepcolor = "white", sepwidth = 0.1, cexRow=3, cexCol=6)
dev.off()

jpeg("./plots/Sweden_W_C_RNA_heatmap_legend2.jpg", width = 1000, height = 800)
heatmap.2 (heatmap1, scale = "row", trace = "none", notecex=3, notecol="black", tracecol="black", Colv = FALSE, labCol = colnames(heatmap1), labRow = rownames(heatmap1), col='rainbow', dendrogram = "none", colsep =c(5), margins =c(0.2, 0.2), sepcolor = "white", sepwidth = 0.1, cexRow=6, cexCol=6)
dev.off()

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
list1 <- DGEList (counts = exp_data, group = treat)

#calculate the normalization factors to adjust the effective library size relative to other libraries in the dataset (norm.factor that minimizes log fold change among samples)
list1 <- calcNormFactors (list1)

#calculate the counts per million
cpm.list1 <- cpm(list1)

#filter out the genes with less than 1 cpm in 6 or fewer libraries (a somewhat arbitrary choice). 
#Genes are usually dropped if they can't possibly be expressed in all the samples for any of the conditions.
list2 <- list1[rowSums(cpm.list1 > 1) >= 6,]
cpm.list2 <- cpm.list1[rowSums(cpm.list1 > 1) >= 6,]

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
jpeg("./plots/Seed_W_C_RNA_MDS.jpg", width = 700, height = 500)
plotMDS (list2, col = col_treat)
dev.off()

#Set the model to use. This one includes the intercept, but other models can be specified that omit the intercept or that have more complex designs. See EdgeR manual for details.
design <- model.matrix (~treat)

#fit the common + tagwise dispersion models
list2 <- estimateGLMCommonDisp (list2, design)
list2 <- estimateGLMTagwiseDisp (list2, design)

#fit a GLM to the data using the tagwise dispersion estimate
glm.list2 <- glmFit (list2, design, dispersion = list2$tagwise.dispersion)
lrt.list2 <- glmLRT (glm.list2)

#get the topTags out of the model
#top <- topTags (lrt.list2, n = 1000)$table
#top <- topTags (lrt.list2, n = 200)$table
top <- topTags (lrt.list2, n = 40)$table # p-value all < 5e-02, v2

fdr<-p.adjust(lrt.list2$table$PValue, method='fdr')

dim(lrt.list2$table[fdr<0.05,])
length(lrt.list2$table$PValue[lrt.list2$table$PValue<0.05])

#make a heatmap by getting the counts per million from each gene and turning them relative proportions (columns add up to 1)
sub1 <- colSums (cpm.list2)
sub2 <- matrix (rep(sub1,nrow (cpm.list2)), c (nrow (cpm.list2),ncol(cpm.list2)),byrow = TRUE)
sub3 <- cpm.list2 / sub2

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
#col_order <- c("C_Swed_15",   "C_Swed_16",   "C_Swed_2", "C_Swed_11",   "C_Swed_3",    "C_Swed_5",    "C_Swed_8",   "W_Swed_10",   "W_Swed_13",   "W_Swed_14",   "W_Swed_16",   "W_Swed_17",  "W_Swed_2",    "W_Swed_5",    "W_Swed_6", 
#"C_Seed_1",    "C_Seed_2", "C_Seed_3", "W_Seed_1",    "W_Seed_2",    "W_Seed_3")
#heatmap2 <- heatmap1[, col_order]

jpeg("./plots/Seed_W_C_RNA_heatmap2.jpg", width = 3200, height = 2500)
heatmap.2 (heatmap1, scale = "row", trace = "none", Colv = FALSE, labCol = colnames(heatmap1), labRow = rownames(heatmap1), col='rainbow', dendrogram = "none", colsep =c(5), margins =c(20,70), sepcolor = "white", sepwidth = 0.1, cexRow=2, cexCol=2)
dev.off()

jpeg("./plots/Seed_W_C_RNA_heatmap_legend2.jpg", width = 1000, height = 800)
heatmap.2 (heatmap1, scale = "row", trace = "none", notecex=3, notecol="black", tracecol="black", Colv = FALSE, labCol = colnames(heatmap1), labRow = rownames(heatmap1), col='rainbow', dendrogram = "none", colsep =c(5), margins =c(0.2, 0.2), sepcolor = "white", sepwidth = 0.1, cexRow=6, cexCol=6)
dev.off()

#-------------------------
#order genes based on most to least expressed in first individual
#heatmap1_df <- as.data.frame(heatmap1)
#heatmap1_df$C_L11 <- as.numeric(as.character(heatmap1_df$C_L11))
#heatmap_ordered_df <- heatmap1_df[order(heatmap1_df$C_L11),]
#heatmap_ordered <- as.matrix(heatmap_ordered_df)

#jpeg("W_C_RNA_heatmap_ordered.jpg", width = 1500, height = 700)
#heatmap.2 (heatmap_ordered, trace = "none", scale = "row", Rowv = FALSE, Colv = FALSE, labCol = colnames(heatmap_ordered), cexCol = 1, dendrogram = "none", labRow = rownames(heatmap_ordered), colsep =c(5), sepcolor = "white", sepwidth = 0.03)
#dev.off()

###################################
# Sweden differential expression - specific genes plotting 
#Do1_01_a00001G00358V1.1

#plot the individual expression values from a single gene:
subdat1 <- exp_data["Do1_01_a00001G00358V1.1",]

jpeg("W_C_RNA_Do1_01_a00001G00358V1.1.jpg", width = 700, height = 500)
stripchart (log10(as.numeric (subdat1)) ~ treat, vertical = T, xlim = c (0.5,2.5), method = "jitter", ylab = "Log10 Expression count", pch=16)
dev.off()

#----------------------------------
# Do1_05_a00001G00328V1.1

#plot the individual expression values from a single gene:
subdat1 <- exp_data["Do1_05_a00001G00328V1.1",]

jpeg("W_C_RNA_Do1_05_a00001G00328V1.1.jpg", width = 700, height = 500)
stripchart (log10(as.numeric (subdat1)) ~ treat, vertical = T, xlim = c (0.5,2.5), method = "jitter", ylab = "Log10 Expression count", pch=16)
dev.off()
#----------------------------------
#Do1_02_a00003G00097V1.1
#plot the individual expression values from a single gene:
subdat1 <- exp_data["Do1_02_a00003G00097V1.1",]

jpeg("W_C_RNA_Do1_02_a00003G00097V1.1.jpg", width = 700, height = 500)
stripchart (log10(as.numeric (subdat1)) ~ treat, vertical = T, xlim = c (0.5,2.5), method = "jitter", ylab = "Log10 Expression count", pch=16)
dev.off()

#----------------------------------
#Do1_06_a00002G00149V1.1
#plot the individual expression values from a single gene:
subdat1 <- exp_data["Do1_06_a00002G00149V1.1",]

jpeg("W_C_RNA_Do1_06_a00002G00149V1.1.jpg", width = 700, height = 500)
stripchart (log10(as.numeric (subdat1)) ~ treat, vertical = T, xlim = c (0.5,2.5), method = "jitter", ylab = "Log10 Expression count", pch=16)
dev.off()
#----------------------------------
#Do1_01_a00001G02208V1.1
#plot the individual expression values from a single gene:
subdat1 <- exp_data["Do1_01_a00001G02208V1.1",]

jpeg("W_C_RNA_Do1_01_a00001G02208V1.1.jpg", width = 700, height = 500)
stripchart (log10(as.numeric (subdat1)) ~ treat, vertical = T, xlim = c (0.5,2.5), method = "jitter", ylab = "Log10 Expression count", pch=16)
dev.off()

###################################
# on local machine
cd /home/Owner/Desktop
scp -v celphin@cedar.computecanada.ca:/home/celphin/projects/rpp-rieseber/celphin/Dryas/RNAseq_analysis/plots/*jpg .
scp -v celphin@cedar.computecanada.ca:/home/celphin/projects/rpp-rieseber/celphin/Dryas/RNAseq_analysis/output/*pdf .
scp -v celphin@cedar.computecanada.ca:/home/celphin/projects/rpp-rieseber/celphin/Dryas/RNAseq_analysis/RNA_heatmap_ordered.jpg .
scp -v celphin@cedar.computecanada.ca:/home/celphin/projects/rpp-rieseber/celphin/Dryas/Reference_genomes/annotation/Dryas_octopetala_H1.xlsx .
scp -v celphin@cedar.computecanada.ca:/home/celphin/projects/rpp-rieseber/celphin/Dryas/RNAseq_analysis/Dryas_octopetala_H1.supercontigs.fa .
scp -v celphin@cedar.computecanada.ca:/home/celphin/projects/rpp-rieseber/celphin/Dryas/RNAseq_analysis/Dryas_octopetala_H1.gff3 .

mkdir LAT_DMRs
scp -v celphin@cedar.computecanada.ca:/home/celphin/projects/rpp-rieseber/celphin/Dryas/DMRs/June2022_Metilene_DMRs/LAT_DMRs* .
scp -v celphin@cedar.computecanada.ca:/home/celphin/projects/rpp-rieseber/celphin/Dryas/RNAseq_analysis/*RNA_Do1_01_a00001G00358V1.1.jpg .
scp -v celphin@cedar.computecanada.ca:/home/celphin/projects/rpp-rieseber/celphin/Dryas/RNAseq_analysis/*RNA_Do1_05_a00001G00328V1.1.jpg .
scp -v celphin@cedar.computecanada.ca:/home/celphin/projects/rpp-rieseber/celphin/Dryas/RNAseq_analysis/W_C_*jpg .
scp -v celphin@cedar.computecanada.ca:/home/celphin/projects/rpp-rieseber/celphin/Other/biol525d.tar.gz

###############################
# transcript count in transcriptome
grep -c ">" Dryas_octopetala_H1.transcript.fa
41181

# 28218 genes had 10 or more reads map in Lat
# 16716 genes had 200 or more reads map in LAT

##########################################
# plot DMRs as heatmaps


tmux new-session -s R
tmux attach-session -t R

cd /home/celphin/projects/rpp-rieseber/celphin/Dryas/RNAseq_analysis/

# get R on server
# open R
module load StdEnv/2020
module load r/4.2.1
module load gdal
module load udunits
module load python
export R_LIBS_USER=/home/celphin/R/x86_64-pc-linux-gnu-library/4.2.1/

R

#paste it in here (i.e. replace my path with yours):
setwd ("/home/celphin/projects/rpp-rieseber/celphin/Dryas/RNAseq_analysis/")

#find your working directory:
getwd()

#load the libraries you will need 
library ("edgeR")
library ("gplots")

#read in the data
mydata <- read.table("gene_names_expression_table.txt", header=TRUE)

CDS <- mydata[,1]
exp_data <- mydata[,c(2,3)]

#extract and turn the column names into a factor
treat <- as.factor (sapply (strsplit(colnames(exp_data),split = ""),"[[",1))

#make a DGElist object out of the data
list1 <- DGEList (counts = exp_data, group = treat)

#calculate the normalization factors to adjust the effective library size relative to other libraries in the dataset (norm.factor that minimizes log fold change among samples)
list1 <- calcNormFactors (list1)

#calculate the counts per million
cpm.list1 <- cpm(list1)

#filter out the genes with less than 1 cpm in 6 or fewer libraries (a somewhat arbitrary choice). 
#Genes are usually dropped if they can't possibly be expressed in all the samples for any of the conditions.
list2 <- list1[rowSums(cpm.list1 > 1) >= 6,]
cpm.list2 <- cpm.list1[rowSums(cpm.list1 > 1) >= 6,]

#generate a multi-dimensional scaling plot
col_treat <- as.character (treat)
col_treat [col_treat == "C"] <- "blue"
col_treat [col_treat == "W"] <- "red"

#plot the MDS graph
jpeg("RNA_MDS.jpg", width = 700, height = 500)
plotMDS (list2, col = col_treat)
#Only 2 columns of data: need at least 3
dev.off()

#Set the model to use. This one includes the intercept, but other models can be specified that omit the intercept or that have more complex designs. See EdgeR manual for details.
design <- model.matrix (~treat)

#fit the common + tagwise dispersion models
list2 <- estimateGLMCommonDisp (list2, design)
list2 <- estimateGLMTagwiseDisp (list2, design)

#fit a GLM to the data using the tagwise dispersion estimate
glm.list2 <- glmFit (list2, design, dispersion = list2$tagwise.dispersion)
lrt.list2 <- glmLRT (glm.list2)

#get the topTags out of the model
top <- topTags (lrt.list2, n = 1000)$table

fdr<-p.adjust(lrt.list2$table$PValue, method='fdr')

dim(lrt.list2$table[fdr<0.05,])
length(lrt.list2$table$PValue[lrt.list2$table$PValue<0.05])

#make a heatmap by getting the counts per million from each gene and turning them relative proportions (columns add up to 1)
sub1 <- colSums (cpm.list2)
sub2 <- matrix (rep(sub1,nrow (cpm.list2)), c (nrow (cpm.list2),ncol(cpm.list2)),byrow = TRUE)
sub3 <- cpm.list2 / sub2

#subset this matrix to just get the genes from above
names1 <- row.names (sub3)
names2 <- row.names (top)
index2 <- names1 %in% names2
heatmap1 <- sub3[index2,]

#play around with the options to make the plot fit what you like for options type ?heatmap.2
jpeg("RNA_heatmap.jpg", width = 700, height = 500)
heatmap.2 (heatmap1, trace = "none", scale = "row", Colv = FALSE, labCol = treat,cexCol = 1, dendrogram = "none", labRow = FALSE, colsep =c(6), sepcolor = "white", sepwidth = 0.03)
dev.off()

#plot the individual expression values from a single gene:
subdat1 <- data["Do1_00107G00003V1.1",]

jpeg("RNA_Do1_00107G00003.jpg", width = 700, height = 500)
stripchart (log10(as.numeric (subdat1)) ~ treat, vertical = T, xlim = c (0.5,2.5), method = "jitter", ylab = "Log10 Expression count")
dev.off()


#############################################
# Compare with DMR positions
# Find closest gene to each DMR
# find any union between DMRs and RNAseq data

# first which share gene contig numbers?




###############################################
# use gff3 and gene annotations to determine information about genes with patterns - what do they do??






################################
# Sept 2022 Results for Sweden and seedlings
#------------------------
> top - only 6 samples 5 sweden and one seedling
                            logFC    logCPM       LR       PValue          FDR
Do1_01_a00001G01943V1.1 -3.716785  3.761415 84.59324 3.665160e-20 4.618835e-16  Do1_01_a00001G01943V1.1	cytochrome P450 family 81 protein;cytochrome P450 family monooxygenase;cytochrome P450 family protein	14	5	6	8	3	AT2G23190,AT4G37370,AT2G23220,AT1G66540,AT4G37400,AT4G37340,AT4G37320,AT4G37330,AT5G36220,AT5G57220,AT4G37360,AT3G28740,AT4G37410,AT4G37430	LotjaGi2g1v0405300,LotjaGi2g1v0405500,LotjaGi4g1v0019000,LotjaGi6g1v0245700,LotjaGi6g1v0245800	Ro01_G29762,Ro01_G29766,Ro01_G30289,Ro01_G30291,Ro01_G30292,Ro06_G28278	Medtr2g437840,Medtr2g437880,Medtr4g035330,Medtr4g094772,Medtr4g094775,Medtr4g095050,Medtr5g016410,Medtr5g016440	FvH4_1g16870,FvH4_1g16930,FvH4_1g23710
Seedling - Do1_04_a00004G00178V1.1  2.509768  4.504091 44.20357 2.959407e-11 1.864722e-07  AT4G11373  Serine/threonine-specific protein kinase
Do1_02_a00004G00438V1.1 -2.469737  3.814916 41.09731 1.448353e-10 6.084047e-07
Do1_03_a00001G00439V1.1 -2.445216  5.856561 37.27391 1.026486e-09 3.233943e-06  Do1_03_a00001G00439V1.1	nodulin-like/MFS transporter	1	4	1	2	2	AT2G30300	LotjaGi1g1v0110700,LotjaGi2g1v0229900_LC,LotjaGi6g1v0088600,LotjaGi6g1v0106400	Ro03_G10642	Medtr3g466040,Medtr3g466110	FvH4_2g11110,FvH4_3g40010
Do1_02_a00003G02101V1.1 -2.542581  5.257169 36.40598 1.602112e-09 4.037962e-06
Do1_01_a00001G02578V1.1  2.973467 10.216749 33.43870 7.354576e-09 1.544706e-05
Do1_05_a00001G02101V1.1  2.084820  3.025305 28.94044 7.463816e-08 1.343700e-04
Do1_02_a00004G00448V1.1 -2.550266  5.092723 28.53135 9.219373e-08 1.452282e-04
Do1_06_a00001G02309V1.1  2.648086  3.954880 28.26877 1.055855e-07 1.478432e-04
Do1_07_a00001G00314V1.1 -2.351876  2.978853 26.98296 2.052570e-07 2.586648e-04
Do1_02_a00001G00593V1.1  2.017230  4.041317 25.66742 4.056172e-07 4.646898e-04   Do1_02_a00001G00593V1.1	glutamate receptor 3.2	9	4	8	2	11	AT5G11210,AT2G29100,AT5G27100,AT2G29110,AT2G24710,AT5G11180,AT2G29120,AT4G31710,AT2G24720	LotjaGi1g1v0373900_LC,LotjaGi1g1v0374200,LotjaGi1g1v0374400_LC,LotjaGi1g1v0376200	Ro02_G10462,Ro02_G16376,Ro02_G21240,Ro02_G21241,Ro02_G25047,Ro02_G25780,Ro02_G27245,Ro03_G26229	Medtr3g105595,Medtr3g105610	FvH4_2g11280,FvH4_2g11340,FvH4_2g11350,FvH4_2g37130,FvH4_2g37140,FvH4_2g37150,FvH4_2g37160,FvH4_3g22790,FvH4_3g22810,FvH4_3g22830,FvH4_7g09130
Do1_03_a00001G00084V1.1  2.108269  2.366549 25.29033 4.931747e-07 4.981012e-04
Do1_03_a00001G00121V1.1 -2.678973  5.520926 25.21119 5.138324e-07 4.981012e-04
Do1_07_a00002G00673V1.1 -1.691772  4.042921 23.87594 1.027478e-06 9.132174e-04
Do1_04_a00001G00835V1.1  3.093479  3.593327 23.76756 1.086991e-06 9.132174e-04
Do1_05_a00001G00461V1.1  1.815305  4.290923 23.17340 1.480309e-06 1.157897e-03
Do1_01_a00005G00038V1.1 -2.021675  1.972812 23.07012 1.561994e-06 1.157897e-03
Do1_03_a00002G00921V1.1 -1.984726  3.986131 21.83972 2.964009e-06 2.075135e-03
Do1_03_a00002G00043V1.1 -2.265313  2.121081 21.36723 3.791969e-06 2.515073e-03
Seedling - Do1_03_a00001G01245V1.1  2.617007  3.270904 21.10359 4.351128e-06 2.741646e-03   H2O and Nitrogen stress response https://www.arabidopsis.org/servlets/TairObject?type=locus&name=AT1G29290
Do1_04_a00001G02314V1.1 -2.111120  3.256554 20.22990 6.867128e-06 4.120931e-03
Do1_06_a00003G00036V1.1  2.887388  6.672914 19.75629 8.797148e-06 5.039166e-03
Do1_06_a00002G01334V1.1  1.721644  8.558247 19.23468 1.155940e-05 6.188454e-03
Do1_07_a00006G00097V1.1 -2.226138  4.307729 19.19768 1.178566e-05 6.188454e-03
Do1_02_a00001G00094V1.1 -1.540599  4.251726 18.87986 1.392145e-05 6.987094e-03
Do1_06_a00001G02951V1.1  1.749372  3.713850 18.77126 1.473713e-05 6.987094e-03
Do1_04_a00002G00114V1.1  1.699485  9.843537 18.74136 1.496997e-05 6.987094e-03
Do1_a00045G00229V1.1    -2.313917  4.880837 18.62185 1.593836e-05 7.173399e-03
Do1_03_a00002G00403V1.1 -2.011494  2.702186 18.28615 1.900838e-05 8.260126e-03

cd /home/celphin/projects/rpp-rieseber/celphin/Dryas/Reference_genomes/

#-------------------------------
# top genes diff exp when including seed and lat samples
> top
                            logFC    logCPM       LR       PValue          FDR
Do1_02_a00004G00437V1.1 -3.061705 2.4188967 40.15780 2.342542e-10 3.564177e-06
Do1_03_a00001G00267V1.1 -4.248927 3.2511249 27.72753 1.396611e-07 1.062472e-03
Do1_02_a00004G00438V1.1 -2.084952 3.6033711 26.46941 2.677452e-07 1.357914e-03
Do1_04_a00001G01021V1.1 -2.403769 1.5073000 21.73363 3.132513e-06 1.191529e-02
Do1_05_a00001G00108V1.1  1.406956 1.9842271 18.34410 1.843890e-05 4.238594e-02
Do1_05_a00001G00111V1.1 -2.445288 2.2379062 18.24666 1.940657e-05 4.238594e-02
Do1_01_a00001G00326V1.1  4.591127 3.5994562 18.23746 1.950060e-05 4.238594e-02
Do1_02_a00001G00017V1.1  4.956106 5.2588043 17.71497 2.565990e-05 4.880192e-02
Do1_02_a00003G00374V1.1  1.366918 4.8938583 17.20089 3.362780e-05 5.684966e-02
Do1_05_a00004G00097V1.1  3.258051 2.2709375 16.91922 3.900455e-05 5.791244e-02
Do1_01_a00004G00182V1.1 -1.462268 0.3829396 16.78472 4.186900e-05 5.791244e-02
Do1_01_a00001G01816V1.1 -1.765247 0.7743929 16.26757 5.499731e-05 6.677046e-02
Do1_07_a00002G01981V1.1 -3.840328 3.5903958 16.19814 5.705001e-05 6.677046e-02
Do1_06_a00001G02495V1.1  3.529003 3.7907165 15.82915 6.932603e-05 6.873449e-02
Do1_06_a00001G02407V1.1 -3.131659 2.8204730 15.58299 7.896180e-05 6.873449e-02
Do1_07_a00002G00394V1.1  2.438521 1.8965830 15.56177 7.985296e-05 6.873449e-02
Do1_02_a00001G00624V1.1  1.539024 6.4207157 15.55123 8.029936e-05 6.873449e-02
Do1_00161G00004V1.1      3.568757 5.1007585 15.49000 8.294295e-05 6.873449e-02
Do1_03_a00002G01018V1.1 -1.178729 2.7992128 15.42526 8.583340e-05 6.873449e-02
Do1_02_a00004G00448V1.1 -2.137437 4.8225760 15.18262 9.759774e-05 6.965044e-02
Do1_06_a00002G01517V1.1  2.728417 3.8550647 15.14940 9.933006e-05 6.965044e-02
Do1_05_a00001G02101V1.1  1.501972 2.9948903 15.12334 1.007105e-04 6.965044e-02
Do1_07_a00006G00097V1.1 -1.589055 4.0907945 14.83034 1.176276e-04 7.781322e-02
Do1_06_a00001G02951V1.1  1.429018 3.8957041 14.57723 1.345304e-04 8.528669e-02
Do1_03_a00002G02560V1.1  1.524894 4.4835134 14.13888 1.697981e-04 1.033391e-01
Do1_02_a00003G01138V1.1  2.596465 6.7901841 14.00359 1.824620e-04 1.060888e-01
Do1_06_a00002G00928V1.1  2.976371 3.0136546 13.94476 1.882614e-04 1.060888e-01
Do1_06_a00001G01079V1.1 -1.603874 3.7612594 13.63564 2.219322e-04 1.153583e-01


#---------------------------------------------
# top 10 for sweden 
                            logFC    logCPM       LR       PValue          FDR
Do1_01_a00004G00316V1.1 -4.937726 3.3503933 56.36843 6.008750e-14 8.648994e-10   Do1_01_a00004G00316V1.1	glutathione S-transferase, amine-terminal domain protein	2	1	1	1	1	AT3G03190,AT5G17220	LotjaGi6g1v0052500	Ro01_G19766	Medtr3g064700	FvH4_1g27460
Do1_02_a00004G00438V1.1 -2.609394 3.8628151 54.94937 1.236752e-13 8.900905e-10   Do1_02_a00004G00438V1.1		0	0	5	0	4			Ro01_G12265,Ro02_G21447,Ro02_G27047,Ro02_G27328,Ro02_G35470		FvH4_1g20010,FvH4_2g12570,FvH4_2g12590,FvH4_2g12610
Do1_02_a00004G00448V1.1 -2.576089 5.2025815 45.48762 1.536065e-11 7.370038e-08   Do1_02_a00004G00451V1.1	60S ribosomal L8-like protein	3	4	2	2	2	AT2G18020,AT3G51190,AT4G36130	LotjaGi1g1v0474100,LotjaGi2g1v0005300,LotjaGi2g1v0382300,LotjaGi4g1v0127000	Ro03_G17469,Ro03_G33007	Medtr5g021800,Medtr7g076940	FvH4_3g11370,FvH4_3g19850
Do1_02_a00004G00437V1.1 -3.582977 2.4061586 43.88584 3.480984e-11 1.252632e-07   Do1_02_a00004G00437V1.1		0	0	5	0	4			Ro01_G12265,Ro02_G21447,Ro02_G27047,Ro02_G27328,Ro02_G35470		FvH4_1g20010,FvH4_2g12570,FvH4_2g12590,FvH4_2g12610
Do1_02_a00003G00097V1.1 -2.718762 2.8171720 39.80851 2.801217e-10 8.064142e-07   Do1_02_a00003G00097V1.1	glutathione S-transferase, amine-terminal domain protein;phi class glutathione S-transferase	7	1	0	4	2	AT2G47730,AT1G02920,AT1G02930,AT2G02930,AT4G02520,AT1G02940,AT1G02950	LotjaGi5g1v0252100		Medtr1g088825,Medtr1g088840,Medtr1g088845,Medtr1g492670	FvH4_2g25200,FvH4_2g25210
Do1_01_a00001G00200V1.1 -5.672076 5.2291364 37.24621 1.041171e-09 2.497770e-06   Do1_01_a00001G00200V1.1	glycoside hydrolase family 18 protein	9	2	1	3	0	AT4G19810,AT4G19750,AT4G19820,AT4G19800,AT4G19770,AT4G19730,AT4G19720,AT4G19760,AT4G19740	LotjaGi4g1v0284700,LotjaGi6g1v0206300	Ro07_G12717	Medtr3g110280,Medtr3g110320,Medtr4g116920
Do1_06_a00002G00149V1.1 -3.894506 2.4260804 35.51684 2.528549e-09 5.199418e-06   Do1_06_a00002G00149V1.1	Chitinase / Hevein / PR-4 / Wheatwin2;Nodule Cysteine-Rich (NCR) secreted peptide	1	1	1	3	3	AT3G04720	LotjaGi1g1v0782200	Ro06_G08266	Medtr1g080780,Medtr6g083200,Medtr6g083310	FvH4_3g05950,FvH4_3g26010,FvH4_6g02840
Do1_01_a00001G00197V1.1 -5.634628 4.1237503 34.03102 5.424035e-09 9.747487e-06

# top 100 for Sweden
                            logFC    logCPM       LR       PValue          FDR
Do1_01_a00004G00316V1.1 -4.937726 3.3503933 56.36843 6.008750e-14 8.648994e-10
Do1_02_a00004G00438V1.1 -2.609394 3.8628151 54.94937 1.236752e-13 8.900905e-10
Do1_02_a00004G00448V1.1 -2.576089 5.2025815 45.48762 1.536065e-11 7.370038e-08
Do1_02_a00004G00437V1.1 -3.582977 2.4061586 43.88584 3.480984e-11 1.252632e-07
***Do1_02_a00003G00097V1.1 -2.718762 2.8171720 39.80851 2.801217e-10 8.064142e-07 # found in seedlings too
Do1_01_a00001G00200V1.1 -5.672076 5.2291364 37.24621 1.041171e-09 2.497770e-06
***Do1_06_a00002G00149V1.1 -3.894506 2.4260804 35.51684 2.528549e-09 5.199418e-06  # found in seedlings too
Do1_01_a00001G00197V1.1 -5.634628 4.1237503 34.03102 5.424035e-09 9.747487e-06
Do1_04_a00001G01834V1.1 -2.894693 2.8670213 33.80420 6.094719e-09 9.747487e-06
Do1_03_a00002G02091V1.1 -2.646768 2.4397700 31.38149 2.119920e-08 3.051412e-05
Do1_03_a00004G00250V1.1 -2.612075 1.3290812 29.60662 5.292379e-08 6.925319e-05
Do1_04_a00004G00178V1.1  3.265860 4.6756189 28.79568 8.043038e-08 9.647624e-05
Do1_03_a00001G00547V1.1 -2.235252 1.6944055 27.94382 1.248894e-07 1.346425e-04
Do1_03_a00001G01348V1.1 -2.257250 1.6906384 27.85202 1.309570e-07 1.346425e-04
Do1_06_a00001G01641V1.1 -2.328598 1.7763909 26.57152 2.539599e-07 2.376102e-04
Do1_02_a00001G01081V1.1 -2.557534 3.1987546 26.49573 2.641214e-07 2.376102e-04
Do1_05_a00001G01780V1.1 -2.384909 3.8935027 26.26532 2.975849e-07 2.519669e-04
Do1_05_a00003G00294V1.1 -1.493102 3.9630156 26.04701 3.332051e-07 2.664530e-04
Do1_01_a00001G01205V1.1  2.101949 3.2381505 24.90712 6.015980e-07 4.557579e-04
Do1_05_a00001G00730V1.1 -2.550101 3.4788031 24.77787 6.433179e-07 4.629959e-04
Do1_07_a00004G00042V1.1 -3.020202 1.6528234 23.91449 1.007109e-06 6.845356e-04
Do1_05_a00001G01862V1.1 -2.292063 3.8695890 23.84108 1.046254e-06 6.845356e-04
Do1_07_a00006G00097V1.1 -1.979547 4.3475373 23.53632 1.225779e-06 7.671247e-04                                                               
Do1_07_a00002G00376V1.1  2.714913 2.4064678 21.91694 2.847102e-06 1.707550e-03
Do1_02_a00001G01347V1.1 -2.064196 3.0290462 21.44015 3.650472e-06 2.101796e-03
Do1_03_a00002G00043V1.1 -2.208495 2.1173633 21.01065 4.567368e-06 2.528565e-03
Do1_04_a00001G03040V1.1 -2.933031 2.3156119 20.83864 4.996467e-06 2.663672e-03
Do1_06_a00002G01303V1.1  1.512377 8.1303000 20.35476 6.433290e-06 3.307171e-03
Do1_04_a00001G03042V1.1 -3.912075 4.9950145 20.23516 6.848273e-06 3.399105e-03
***Do1_01_a00001G02208V1.1 -2.431277 3.5826875 20.03725 7.594807e-06 3.643988e-03 # found in seedlings too
Do1_01_a00001G01851V1.1 -2.574457 3.4205898 19.30295 1.115339e-05 4.753073e-03
***Do1_01_a00001G00358V1.1 -1.822147 4.1516709 19.27778 1.130139e-05 4.753073e-03 # found in seedlings too - reversed
Do1_06_a00002G01434V1.1 -1.667573 4.3889936 19.24057 1.152382e-05 4.753073e-03
Do1_00161G00004V1.1      3.924668 5.4232814 19.22930 1.159205e-05 4.753073e-03
Do1_03_a00001G00439V1.1 -1.505353 5.8552443 19.14341 1.212555e-05 4.753073e-03
Do1_05_a00001G01513V1.1 -1.524542 2.5981144 19.11538 1.230496e-05 4.753073e-03
*** Do1_05_a00001G00328V1.1 -1.777468 2.3282554 19.05713 1.268629e-05 4.753073e-03 # found in seedlings too
Do1_06_a00001G00062V1.1 -1.833667 1.1502824 18.99262 1.312250e-05 4.753073e-03
Do1_03_a00002G01018V1.1 -1.365641 2.9412341 18.98065 1.320510e-05 4.753073e-03
Do1_06_a00003G00036V1.1  2.655676 6.3976855 18.98016 1.320848e-05 4.753073e-03
Do1_04_a00001G00990V1.1 -4.598518 4.4665685 18.84660 1.416634e-05 4.973422e-03
Do1_01_a00001G00356V1.1 -1.643048 1.5057245 18.38641 1.803396e-05 6.180496e-03
Do1_06_a00001G00091V1.1 -2.766267 4.0821241 18.20944 1.978957e-05 6.624442e-03
Do1_01_a00004G01379V1.1 -2.147850 0.9833846 17.42561 2.987726e-05 9.293107e-03
Do1_03_a00002G00280V1.1 -2.054114 5.2362859 17.41415 3.005793e-05 9.293107e-03
Do1_04_a00001G01827V1.1 -2.064466 1.4288917 17.41126 3.010380e-05 9.293107e-03
Do1_04_a00001G01620V1.1 -1.586619 1.8153516 17.36227 3.088976e-05 9.293107e-03
Do1_03_a00002G01336V1.1 -1.420071 0.8888886 17.32566 3.149062e-05 9.293107e-03
Do1_01_a00004G00182V1.1 -1.676674 0.4366857 17.31693 3.163556e-05 9.293107e-03
Do1_04_a00001G01981V1.1 -1.844841 3.2228330 17.14751 3.458620e-05 9.956674e-03
Do1_02_a00003G01072V1.1 -1.383869 0.9506148 16.92808 3.882306e-05 1.095724e-02
Do1_06_a00002G02412V1.1  1.878878 2.8631565 16.86689 4.009493e-05 1.103700e-02
Do1_04_a00001G01021V1.1 -2.364177 1.7274232 16.84130 4.063921e-05 1.103700e-02
Do1_04_a00001G02846V1.1 -1.272319 4.6369178 16.77971 4.197970e-05 1.104564e-02
Do1_01_a00004G01382V1.1 -3.974093 1.6585795 16.76952 4.220578e-05 1.104564e-02
Do1_07_a00002G01981V1.1 -4.262134 4.0277386 16.54464 4.751808e-05 1.221384e-02
Do1_06_a00002G00905V1.1 -1.546116 2.2856647 16.27776 5.470226e-05 1.381376e-02
Do1_02_a00002G00363V1.1 -1.652936 1.2832137 16.14074 5.880536e-05 1.454349e-02
Do1_03_a00002G02460V1.1  2.799788 2.4572234 16.11491 5.961276e-05 1.454349e-02
Do1_05_a00001G00111V1.1 -2.712357 2.3834375 16.04986 6.169625e-05 1.480093e-02
Do1_05_a00001G00596V1.1 -1.512903 2.0407681 15.94153 6.532942e-05 1.541560e-02
Do1_05_a00001G00108V1.1  1.406473 1.8279896 15.81620 6.980220e-05 1.601148e-02
Do1_07_a00004G00465V1.1 -1.529871 0.8691872 15.80870 7.007942e-05 1.601148e-02
Do1_06_a00001G02407V1.1 -3.527966 3.2256476 15.75924 7.193575e-05 1.617880e-02
Do1_04_a00001G01197V1.1  1.342146 3.2913900 15.46409 8.408753e-05 1.862086e-02
Do1_07_a00002G01266V1.1 -1.782649 0.6953799 15.33013 9.026552e-05 1.968609e-02
Do1_06_a00002G01541V1.1 -2.077214 1.9856818 15.20124 9.664012e-05 2.076176e-02
Do1_06_a00002G01112V1.1  1.826609 1.7495469 15.13523 1.000779e-04 2.118414e-02
Do1_04_a00004G00161V1.1  2.473197 4.1847507 15.09459 1.022557e-04 2.133143e-02
Do1_01_a00001G02078V1.1 -1.526306 1.7987827 15.02982 1.058255e-04 2.176075e-02
Do1_01_a00002G00102V1.1 -2.937303 1.4908883 14.93658 1.111863e-04 2.254107e-02
Do1_03_a00002G02348V1.1 -2.715939 2.0066346 14.86411 1.155406e-04 2.309849e-02
Do1_04_a00001G01481V1.1 -1.569561 2.7614658 14.78392 1.205592e-04 2.362476e-02
Do1_07_a00001G00240V1.1 -1.649805 2.0511994 14.76995 1.214556e-04 2.362476e-02
Do1_07_a00002G02494V1.1 -2.573611 2.4949976 14.62535 1.311388e-04 2.489595e-02
Do1_07_a00006G00157V1.1 -1.892021 1.5484233 14.62088 1.314501e-04 2.489595e-02
Do1_05_a00001G00393V1.1  1.411951 1.4012667 14.43072 1.454109e-04 2.718240e-02
Do1_06_a00002G00118V1.1 -1.700487 0.6845367 14.33092 1.533254e-04 2.783812e-02
Do1_01_a00001G00911V1.1  1.281405 4.1916495 14.32662 1.536766e-04 2.783812e-02
Do1_a00045G00281V1.1     1.803696 1.9615394 14.31387 1.547207e-04 2.783812e-02
Do1_07_a00002G01011V1.1 -3.615584 3.1373984 14.28605 1.570241e-04 2.790377e-02
Do1_05_a00001G01791V1.1 -2.850870 4.6495228 14.23812 1.610743e-04 2.827444e-02
Do1_01_a00001G02361V1.1  1.440814 7.1364535 14.21505 1.630613e-04 2.827837e-02
***Do1_03_a00001G00809V1.1 -1.456689 2.3343819 14.16103 1.678101e-04 2.851475e-02 *** also changes in offspring with H and L
Do1_04_a00001G02626V1.1 -1.346572 5.3766293 14.15458 1.683864e-04 2.851475e-02
Do1_01_a00001G02662V1.1 -2.628102 1.4904016 14.10998 1.724269e-04 2.885946e-02
Do1_01_a00002G00188V1.1 -1.265952 3.5788784 13.91347 1.914218e-04 3.167040e-02
Do1_04_a00005G00185V1.1 -2.862011 4.7953014 13.88741 1.940948e-04 3.174773e-02
Do1_06_a00002G00928V1.1  3.663478 3.3019270 13.84028 1.990239e-04 3.218820e-02
Do1_03_a00001G00530V1.1  1.466239 1.8384588 13.80415 2.028876e-04 3.244850e-02
Do1_03_a00001G00377V1.1  1.476212 2.6528397 13.77415 2.061532e-04 3.260845e-02
Do1_05_a00001G02839V1.1  1.158489 2.0568761 13.73332 2.106838e-04 3.296285e-02
Do1_04_a00001G00232V1.1 -1.928271 2.6669384 13.56960 2.298789e-04 3.557933e-02
Do1_04_a00001G01172V1.1 -1.040471 7.4656205 13.52968 2.348197e-04 3.595739e-02
Do1_07_a00002G01798V1.1  1.557759 0.6296865 13.32095 2.624571e-04 3.976640e-02
Do1_04_a00001G02314V1.1 -1.544218 2.9946892 13.27302 2.692529e-04 4.037111e-02
Do1_06_a00001G02410V1.1 -1.509296 3.9386858 13.13601 2.896735e-04 4.263921e-02
Do1_02_a00004G00186V1.1 -1.338257 4.1863529 13.13194 2.903045e-04 4.263921e-02
Do1_07_a00002G01793V1.1  1.208946 3.6795129 13.05976 3.017074e-04 4.375931e-02
Do1_05_a00001G02219V1.1 -2.224671 2.3220940 13.04552 3.040108e-04 4.375931e-02

#-----------------------------
# total all sites DMRs
Do1_01_a00004 4252688 4254212  ID: Do1_01_a00004G00692
Do1_03_a00002 5353037 5355212  ID: Do1_03_a00002G00984
Do1_03_a00002 12037975 12039326 ID: Do1_03_a00002G02151
Do1_03_a00002 12806950 12808303 ID: Do1_03_a00002G02281
Do1_03_a00004 1829129 1830473   ID: Do1_03_a00004G00312
Do1_05_a00001 17164475 17165996 ID: Do1_05_a00001G03065
Do1_05_a00003 1152846 1153987   ID: Do1_05_a00003G00206
Do1_05_a00003 1480241 1481901  ID: Do1_05_a00003G00261 - no clear gene
Do1_05_a00004 324260 325804    ID: Do1_05_a00004G00060 and 61
Do1_06_a00002 2718020 2719459  ID: Do1_06_a00002G00568
Do1_07_a00002 11005173 11006813 ID: Do1_07_a00002G02083
Do1_07_a00002 12442551 12444127  ID: Do1_07_a00002G02343
Do1_07_a00004 2484955 2486400   ID: Do1_07_a00004G00413
Do1_a00028 195651 196795 ID: Do1_a00028G00040

#----------------------------------
# seedling differential expression (warming vs control parents)

> top
                            logFC    logCPM       LR       PValue          FDR
Do1_07_a00002G02277V1.1 -4.127065  3.515385 54.91575 1.258088e-13 1.589720e-09
Do1_a00045G00108V1.1     4.806968  4.929651 42.08255 8.750033e-11 5.528271e-07
Do1_07_a00001G00368V1.1  2.923005  2.960172 36.52403 1.507954e-09 5.081562e-06
***Do1_01_a00001G00358V1.1  3.830079  5.023377 36.39810 1.608598e-09 5.081562e-06 # found in parents and seedlings!! reversed
***Do1_05_a00001G00328V1.1  2.989188  2.816902 34.24767 4.852570e-09 1.100638e-05 # found in parents and seedlings!!
Do1_04_a00001G02757V1.1 -2.494756  3.437855 34.10332 5.226201e-09 1.100638e-05
Do1_05_a00001G02550V1.1 -2.637157  3.375944 32.00690 1.536260e-08 2.486711e-05
***Do1_01_a00001G02208V1.1  4.600832  5.900678 31.95930 1.574366e-08 2.486711e-05 # found in parents and seedlings!!
Do1_07_a00002G00277V1.1 -2.924173  4.490836 29.38813 5.923925e-08 8.317190e-05
Do1_01_a00004G01687V1.1  5.305166  7.801042 28.27123 1.054518e-07 1.332488e-04
Do1_05_a00001G00945V1.1  3.105283  5.395090 28.04557 1.184919e-07 1.361148e-04
***Do1_06_a00002G00149V1.1  4.055430  4.412215 26.92856 2.111158e-07 2.223050e-04 # found in parents too!
Do1_04_a00001G02266V1.1  2.329759  3.066702 26.26270 2.979886e-07 2.896450e-04
Do1_01_a00001G02221V1.1  2.189273  6.166206 24.86787 6.139736e-07 5.541550e-04
***Do1_02_a00003G00097V1.1  4.756650  5.719489 24.54903 7.244275e-07 6.102578e-04 # found in parents too!
Do1_00230G00005V1.1     -2.556451  3.043149 23.55050 1.216779e-06 9.609512e-04
Do1_07_a00005G00229V1.1 -4.260806  5.382966 23.00737 1.613818e-06 1.199542e-03
Do1_01_a00001G01228V1.1 -2.613710  2.828726 21.24120 4.049654e-06 2.842857e-03
***Do1_04_a00001G02267V1.1  3.951548  4.189929 20.99570 4.603143e-06 3.061332e-03 # also seen in offspring with H and L
Do1_02_a00004G01142V1.1  2.627041  5.364756 20.52283 5.892437e-06 3.722842e-03
Do1_05_a00003G00219V1.1 -2.011393  2.580611 20.41000 6.250241e-06 3.760859e-03
Do1_04_a00001G01716V1.1  2.897706  3.814696 19.27316 1.132875e-05 6.223918e-03


Do1_04_a00001G01717V1.1  2.897706  3.814696 19.27316 1.132875e-05 6.223918e-03
Do1_03_a00002G00774V1.1 -1.911874  2.817522 19.12639 1.223416e-05 6.441284e-03
Do1_06_a00002G01422V1.1  3.769847 10.928787 18.29637 1.890669e-05 9.556198e-03
Do1_05_a00001G00119V1.1 -2.219207  1.669764 18.17875 2.011105e-05 9.590702e-03
Do1_06_a00002G01118V1.1 -1.813189  6.548455 18.14292 2.049295e-05 9.590702e-03
Do1_05_a00001G00258V1.1  1.687598  3.774643 17.60462 2.719265e-05 1.227165e-02
Do1_04_a00001G01392V1.1  2.272308  3.741957 16.74081 4.284929e-05 1.757273e-02
Do1_a00045G00218V1.1     3.409980  3.155601 16.69474 4.390256e-05 1.757273e-02 # found in H and L
Do1_03_a00002G02422V1.1 -1.755947  3.126458 16.64463 4.507781e-05 1.757273e-02
Do1_04_a00001G00703V1.1 -1.845622  4.022794 16.64263 4.512530e-05 1.757273e-02
Do1_00107G00006V1.1     -2.989527  5.227925 16.61065 4.589270e-05 1.757273e-02
Do1_05_a00001G02341V1.1 -2.701115  5.258763 16.25991 5.521988e-05 2.052231e-02
Do1_06_a00001G02659V1.1  2.973678  2.992329 16.08864 6.044550e-05 2.182255e-02
Do1_07_a00002G00155V1.1  1.814681  2.888628 15.97567 6.416186e-05 2.252081e-02
Do1_04_a00001G01370V1.1 -1.943306  3.955606 15.86828 6.790726e-05 2.319125e-02
Do1_01_a00001G00050V1.1  1.858584  7.497361 15.78926 7.080336e-05 2.354398e-02
Do1_02_a00001G01951V1.1  4.249625  5.693181 15.68216 7.492765e-05 2.427656e-02 # found in H and L
Do1_06_a00002G01474V1.1  1.663296  6.556697 15.61319 7.771082e-05 2.454885e-02
Do1_02_a00001G00160V1.1  2.182521  4.180311 14.98251 1.085124e-04 3.344299e-02
Do1_06_a00001G00167V1.1  3.188600  5.714918 14.84088 1.169721e-04 3.519188e-02
Do1_02_a00003G01840V1.1  2.191786  2.859280 14.23897 1.610016e-04 4.731201e-02
Do1_06_a00001G01193V1.1  2.113457  6.525282 14.14351 1.693805e-04 4.864301e-02

44 significant differentially expressed genes

#-------------------------------
# seedling low vs high environments

> top                                                                                                                                                 [2/1931]
                            logFC   logCPM       LR       PValue          FDR
**Do1_03_a00001G00809V1.1 -6.121451 6.362631 31.09422 2.458019e-08 0.0003105953 # seen in parents
Do1_07_a00001G00557V1.1  3.028724 3.382305 27.15209 1.880614e-07 0.0011881719
Do1_06_a00002G00856V1.1  2.316840 3.724584 25.98154 3.446972e-07 0.0013699127
Do1_03_a00001G01064V1.1  2.302416 8.209176 25.53845 4.336539e-07 0.0013699127
Do1_02_a00003G00882V1.1 -2.729326 5.112801 23.68195 1.136435e-06 0.0028719982
**Do1_04_a00001G02267V1.1 -3.932401 4.189619 20.95099 4.711839e-06 0.0092820331 # seen in W and C origins seedling
Do1_02_a00003G00284V1.1 -2.398292 2.538923 20.78366 5.141994e-06 0.0092820331
Do1_01_a00001G01844V1.1  3.421113 4.866874 20.03315 7.611113e-06 0.0111151420
Do1_a00045G00218V1.1    -3.570593 3.155009 19.95786 7.916768e-06 0.0111151420 # seen in W and C origins seedling
Do1_01_a00001G01541V1.1  2.477973 3.625364 19.64500 9.324677e-06 0.0117826614
Do1_06_a00002G00043V1.1 -4.978372 6.665097 19.45072 1.032289e-05 0.0118581834
Do1_06_a00001G01852V1.1 -3.899395 5.607365 19.17177 1.194673e-05 0.0125799071
Do1_06_a00001G01952V1.1 -3.741381 5.355062 18.88939 1.385212e-05 0.0133523743
Do1_03_a00002G00790V1.1 -1.768679 5.217219 18.76395 1.479370e-05 0.0133523743
Do1_02_a00001G01951V1.1 -4.474077 5.693061 18.51764 1.683388e-05 0.0141808595 # seen in W and C origins seedling
Do1_04_a00001G03073V1.1 -3.079991 5.494363 18.16970 2.020677e-05 0.0159582991
Do1_07_a00002G01878V1.1 -4.029855 6.265620 17.56894 2.770772e-05 0.0205949817
Do1_02_a00003G01138V1.1 -3.850941 8.891820 17.22102 3.327343e-05 0.0214667413
Do1_06_a00003G00075V1.1 -2.881696 3.399148 17.18632 3.388686e-05 0.0214667413
Do1_04_a00001G00479V1.1 -1.853462 3.289071 17.17144 3.415327e-05 0.0214667413
Do1_04_a00001G01725V1.1 -2.075464 2.682200 17.08859 3.567597e-05 0.0214667413
Do1_03_a00001G00649V1.1  1.967929 2.744378 16.74312 4.279722e-05 0.0231389392
Do1_00278G00003V1.1     -2.340140 6.398020 16.70944 4.356371e-05 0.0231389392
Do1_06_a00001G01609V1.1 -3.654998 6.006032 16.65818 4.475689e-05 0.0231389392
Do1_07_a00004G00348V1.1  2.389688 4.149870 16.61532 4.577979e-05 0.0231389392
Do1_03_a00003G00131V1.1 -2.140674 2.511495 16.28297 5.455205e-05 0.0265122980
Do1_06_a00002G01284V1.1 -2.053918 5.494367 16.01835 6.273137e-05 0.0293582813
Do1_04_a00001G00539V1.1 -1.902693 2.608443 15.39644 8.715243e-05 0.0393306467
Do1_07_a00002G00942V1.1 -2.867124 5.234822 14.91923 1.122132e-04 0.0488939872
Do1_01_a00001G00753V1.1 -3.914782 3.889935 14.84578 1.166687e-04 0.0491408424


#-----------------------------
# Sweden DMRs
 cd /home/celphin/projects/rpp-rieseber/celphin/Dryas/DMRs/June2022_Metilene_DMRs/LAT_DMRs* .

# put in IGV to view where matches genes with differential expression match mathylation - none clearly
cd /home/celphin/projects/rpp-rieseber/celphin/Dryas/RNAseq_analysis/output/

grep "Do1_01_a00004G00692" C_L*.genes.results
C_L11.genes.results:Do1_01_a00004G00692V1.1     Do1_01_a00004G00692V1.1 229.00  46.63   0.00    0.00    0.00
C_L15.genes.results:Do1_01_a00004G00692V1.1     Do1_01_a00004G00692V1.1 229.00  44.27   3.00    1.60    2.37
C_L16.genes.results:Do1_01_a00004G00692V1.1     Do1_01_a00004G00692V1.1 229.00  41.19   0.00    0.00    0.00
C_L2.genes.results:Do1_01_a00004G00692V1.1      Do1_01_a00004G00692V1.1 229.00  44.43   0.00    0.00    0.00
C_L3.genes.results:Do1_01_a00004G00692V1.1      Do1_01_a00004G00692V1.1 229.00  44.95   1.00    0.48    0.73
grep "Do1_01_a00004G00692" W_L*.genes.results
W_L13.genes.results:Do1_01_a00004G00692V1.1     Do1_01_a00004G00692V1.1 229.00  48.03   0.00    0.00    0.00
W_L14.genes.results:Do1_01_a00004G00692V1.1     Do1_01_a00004G00692V1.1 229.00  46.37   0.00    0.00    0.00
W_L16.genes.results:Do1_01_a00004G00692V1.1     Do1_01_a00004G00692V1.1 229.00  42.60   0.00    0.00    0.00
W_L17.genes.results:Do1_01_a00004G00692V1.1     Do1_01_a00004G00692V1.1 229.00  46.41   0.00    0.00    0.00
W_L2.genes.results:Do1_01_a00004G00692V1.1      Do1_01_a00004G00692V1.1 229.00  46.41   0.00    0.00    0.00

grep "Do1_03_a00002G00984" C_L*.genes.results
grep "Do1_03_a00002G00984" W_L*.genes.results


