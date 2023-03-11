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


