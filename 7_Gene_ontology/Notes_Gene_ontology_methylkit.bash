###########################
# Gene ontology analysis
# Nov 2024 
# ErmineJ
##########################

# copy over relevant files

tmux new-session -s GO1
tmux attach-session -t GO1

mkdir /lustre04/scratch/celphin/Dryas/GO_enrichment/
cd /lustre04/scratch/celphin/Dryas/GO_enrichment/


cd /home/celphin/scratch/Dryas/MS_Dryas_Merged_Data/original_data
cp interproscan_dryas_full.tsv /lustre04/scratch/celphin/Dryas/GO_enrichment/
cp Dryas_octopetala_H1.gff3 /lustre04/scratch/celphin/Dryas/GO_enrichment/
cd ..

# get merged data from metilene
#cp home/celphin/scratch/Dryas/MS_Dryas_Merged_Data/original_data/Gene_DMR_RNA_GO_Merged_table.tsv /lustre04/scratch/celphin/Dryas/GO_enrichment/metilene

# get merged data from methylkit
cp /lustre04/scratch/celphin/Dryas/methylkit_merged_data/genes_RNA_MethylkitDMR_merged_data.tsv /lustre04/scratch/celphin/Dryas/GO_enrichment

cd /lustre04/scratch/celphin/Dryas/GO_enrichment/
# count of total genes
grep mRNA Dryas_octopetala_H1.gff3 | wc -l
# 41181


##########################
# GENE IDs for DMRs and DEGs

module load StdEnv/2023
module load r/4.4.0

R

library(dplyr)
library(tidyr)

path="/lustre04/scratch/celphin/Dryas/GO_enrichment"
gene_ont <- read.delim(paste0(path,"/interproscan_dryas_full3.tsv"), header = TRUE, sep = "\t", na.strings = "-")

# load data 
DMR_DEG <- read.delim(paste0(path,"/genes_RNA_MethylkitDMR_merged_data.tsv"), header = TRUE, sep = "\t")


#---------------------------
# Get Origins of DEGs and DMRs from filenames
head(DMR_DEG)
colnames(DMR_DEG)

 # [1] "Scaffold"    "Gene_Start"  "Gene_End"    "Gene"        "gene"
 # [6] "INTPRO"      "descrip1"    "descrip2"    "GOterm"      "logFC"
# [11] "logCPM"      "LR"          "PValue"      "RNAsite"     "UpDown"
# [16] "chr"         "start"       "end"         "site"        "context"
# [21] "random"      "perdiff"     "program"     "strand"      "pvalue"
# [26] "qvalue.x"    "meth.diff.x" "qvalue.y"    "meth.diff.y" "NumbCpG"
# [31] "Mean1"       "Mean2"

#------------------------
#Get a list of the specific genes for each group 

unique(DMR_DEG$site)

 # [1] NA         "Pheno"    "SE_HL"    "LAT_W_C"  "HL"       "ALAS_W_C"
 # [7] "Alex_W_C" "MEAD_W_C" "WILL_W_C" "FERT_W_C" "SVAL_W_C" "SE_W_C"
# [13] "DRY_W_C"  "Total"    "Parent"   "CASS_W_C"


# "ALAS_W_C"
DMR_Alaska_W_C <- DMR_DEG[which(DMR_DEG$site=="ALAS_W_C" & DMR_DEG$perdiff== 10 & DMR_DEG$random=="non-rand" & DMR_DEG$context=="CpG" & DMR_DEG$program=="non-overdisp"),]
Alaska_W_C <- as.data.frame(cbind(DMR_Alaska_W_C$Gene, DMR_Alaska_W_C$qvalue.x))
Alaska_W_C <- distinct(Alaska_W_C)
write.table(Alaska_W_C, "Alaska_W_C_genes.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
#142

# "LAT_W_C"
DMR_Sweden_W_C <- DMR_DEG[which(DMR_DEG$site=="LAT_W_C" & DMR_DEG$perdiff== 10 & DMR_DEG$random=="non-rand" & DMR_DEG$context=="CpG" & DMR_DEG$program=="non-overdisp"),]
Sweden_W_C <- as.data.frame(cbind(DMR_Sweden_W_C$Gene, DMR_Sweden_W_C$qvalue.x))
Sweden_W_C <- distinct(Sweden_W_C)
write.table(Sweden_W_C, "Sweden_W_C_genes.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
# 808

# "CASS_W_C"
DMR_Nunavut_W_C <- DMR_DEG[which(DMR_DEG$site=="CASS_W_C" & DMR_DEG$perdiff== 10 & DMR_DEG$random=="non-rand" & DMR_DEG$context=="CpG" & DMR_DEG$program=="non-overdisp"),]
Nunavut_W_C <- as.data.frame(cbind(DMR_Nunavut_W_C$Gene, DMR_Nunavut_W_C$qvalue.x))
Nunavut_W_C <- distinct(Nunavut_W_C)
write.table(Nunavut_W_C, "Nunavut_W_C_genes.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
#78

# "SVAL_W_C"
DMR_Svalbard_W_C <- DMR_DEG[which(DMR_DEG$site=="SVAL_W_C" & DMR_DEG$perdiff== 10 & DMR_DEG$random=="non-rand" & DMR_DEG$context=="CpG" & DMR_DEG$program=="non-overdisp"),]
Svalbard_W_C <- as.data.frame(cbind(DMR_Svalbard_W_C$Gene, DMR_Svalbard_W_C$qvalue.x))
Svalbard_W_C <- distinct(Svalbard_W_C)
write.table(Svalbard_W_C, "Svalbard_W_C_genes.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
#103

# "Pheno"
DMR_Mat_Sen <- DMR_DEG[which(DMR_DEG$site=="Pheno" & DMR_DEG$perdiff== 10 & DMR_DEG$random=="non-rand" & DMR_DEG$context=="CpG" & DMR_DEG$program=="non-overdisp"),]
Mat_Sen <- as.data.frame(cbind(DMR_Mat_Sen$Gene, DMR_Mat_Sen$qvalue.x))
Mat_Sen <- distinct(Mat_Sen)
write.table(Mat_Sen, "Mat_Sen_genes.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
# 8584

# "HL"
DMR_Wild_Lat_L_H <- DMR_DEG[which(DMR_DEG$site=="HL" & DMR_DEG$perdiff== 10 & DMR_DEG$random=="non-rand" & DMR_DEG$context=="CpG" & DMR_DEG$program=="non-overdisp"),]
Wild_Lat_L_H <- as.data.frame(cbind(DMR_Wild_Lat_L_H$Gene, DMR_Wild_Lat_L_H$qvalue.x))
Wild_Lat_L_H <- distinct(Wild_Lat_L_H)
write.table(Wild_Lat_L_H, "Wild_Lat_L_H_genes.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
#7621

# "SE_W_C"
DMR_SE_W_C <- DMR_DEG[which(DMR_DEG$site=="SE_W_C" & DMR_DEG$perdiff== 10 & DMR_DEG$random=="non-rand" & DMR_DEG$context=="CpG" & DMR_DEG$program=="non-overdisp"),]
SE_W_C <- as.data.frame(cbind(DMR_SE_W_C$Gene, DMR_SE_W_C$qvalue.x))
SE_W_C <- distinct(SE_W_C)
write.table(SE_W_C, "SE_W_C_genes.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
#1434

# "SE_HL"
DMR_SE_L_H <- DMR_DEG[which(DMR_DEG$site=="SE_HL" & DMR_DEG$perdiff== 10 & DMR_DEG$random=="non-rand" & DMR_DEG$context=="CpG" & DMR_DEG$program=="non-overdisp"),]
SE_L_H <- as.data.frame(cbind(DMR_SE_L_H$Gene, DMR_SE_L_H$qvalue.x))
SE_L_H <- distinct(SE_L_H)
write.table(SE_L_H, "SE_L_H_genes.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
# 388

#-------------------
# Extract RNA info separately

unique(DMR_DEG$RNAsite)
# [1] NA         "SE_W_C"   "LAT_W_C"  "NORW_W_C" "ALAS_W_C" "ALEX_W_C"

# "LAT_W_C"
LAT_RNA <- DMR_DEG[which(DMR_DEG$RNAsite=="LAT_W_C"),]
LAT_RNA_gene <- as.data.frame(cbind(LAT_RNA$Gene, LAT_RNA$PValue))
LAT_RNA_gene <- distinct(LAT_RNA_gene)
write.table(LAT_RNA_gene, "Sweden_RNA_genes.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
#59

# "ALAS_W_C"
ALAS_RNA <- DMR_DEG[which(DMR_DEG$RNAsite=="ALAS_W_C"),]
ALAS_RNA_gene <- as.data.frame(cbind(ALAS_RNA$Gene, ALAS_RNA$PValue))
ALAS_RNA_gene <- distinct(ALAS_RNA_gene)
write.table(ALAS_RNA_gene , "Alaska_RNA_genes.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
#19

# "ALEX_W_C"
ALEX_RNA <- DMR_DEG[which(DMR_DEG$RNAsite=="ALEX_W_C"),]
ALEX_RNA_gene <- as.data.frame(cbind(ALEX_RNA$Gene, ALEX_RNA$PValue))
ALEX_RNA_gene <- distinct(ALEX_RNA_gene)
write.table(ALEX_RNA_gene , "Nunavut_RNA_genes.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
#21

# "NORW_W_C"
NORW_RNA <- DMR_DEG[which(DMR_DEG$RNAsite=="NORW_W_C"),]
NORW_RNA_gene <- as.data.frame(cbind(NORW_RNA$Gene, NORW_RNA$PValue))
NORW_RNA_gene <- distinct(NORW_RNA_gene)
write.table(NORW_RNA_gene , "Norway_RNA_genes.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
#6

# "SE_W_C"
SE_RNA <- DMR_DEG[which(DMR_DEG$RNAsite=="SE_W_C"),]
SE_RNA_gene <- as.data.frame(cbind(SE_RNA$Gene, SE_RNA$PValue))
SE_RNA_gene <- distinct(SE_RNA_gene)
write.table(SE_RNA_gene , "Seedling_RNA_genes.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
#45

#--------------------------
# Extract for random as well

# "ALAS_W_C"
DMR_Alaska_W_C <- DMR_DEG[which(DMR_DEG$site=="ALAS_W_C" & DMR_DEG$perdiff== 10 & DMR_DEG$random=="rand" & DMR_DEG$context=="CpG"),]
Alaska_W_C <- as.data.frame(cbind(DMR_Alaska_W_C$Gene, DMR_Alaska_W_C$qvalue.x))
Alaska_W_C <- distinct(Alaska_W_C)
write.table(Alaska_W_C, "Alaska_W_C_randgenes.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
#13

# "LAT_W_C"
DMR_Sweden_W_C <- DMR_DEG[which(DMR_DEG$site=="LAT_W_C" & DMR_DEG$perdiff== 10 & DMR_DEG$random=="rand" & DMR_DEG$context=="CpG"),]
Sweden_W_C <- as.data.frame(cbind(DMR_Sweden_W_C$Gene, DMR_Sweden_W_C$qvalue.x))
Sweden_W_C <- distinct(Sweden_W_C)
write.table(Sweden_W_C, "Sweden_W_C_randgenes.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
#0

# "CASS_W_C"
DMR_Nunavut_W_C <- DMR_DEG[which(DMR_DEG$site=="CASS_W_C" & DMR_DEG$perdiff== 10 & DMR_DEG$random=="rand" & DMR_DEG$context=="CpG"),]
Nunavut_W_C <- as.data.frame(cbind(DMR_Nunavut_W_C$Gene, DMR_Nunavut_W_C$qvalue.x))
Nunavut_W_C <- distinct(Nunavut_W_C)
write.table(Nunavut_W_C, "Nunavut_W_C_randgenes.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
#0

# "SVAL_W_C"
DMR_Svalbard_W_C <- DMR_DEG[which(DMR_DEG$site=="SVAL_W_C" & DMR_DEG$perdiff== 10 & DMR_DEG$random=="rand" & DMR_DEG$context=="CpG"),]
Svalbard_W_C <- as.data.frame(cbind(DMR_Svalbard_W_C$Gene, DMR_Svalbard_W_C$qvalue.x))
Svalbard_W_C <- distinct(Svalbard_W_C)
write.table(Svalbard_W_C, "Svalbard_W_C_randgenes.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
#123

# "Pheno"
DMR_Mat_Sen <- DMR_DEG[which(DMR_DEG$site=="Pheno" & DMR_DEG$perdiff== 10 & DMR_DEG$random=="rand" & DMR_DEG$context=="CpG"),]
Mat_Sen <- as.data.frame(cbind(DMR_Mat_Sen$Gene, DMR_Mat_Sen$qvalue.x))
Mat_Sen <- distinct(Mat_Sen)
write.table(Mat_Sen, "Mat_Sen_randgenes.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
#0

# "HL"
DMR_Wild_Lat_L_H <- DMR_DEG[which(DMR_DEG$site=="HL" & DMR_DEG$perdiff== 10 & DMR_DEG$random=="rand" & DMR_DEG$context=="CpG"),]
Wild_Lat_L_H <- as.data.frame(cbind(DMR_Wild_Lat_L_H$Gene, DMR_Wild_Lat_L_H$qvalue.x))
Wild_Lat_L_H <- distinct(Wild_Lat_L_H)
write.table(Wild_Lat_L_H, "Wild_Lat_L_H_randgenes.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
#0

# "SE_W_C"
DMR_SE_W_C <- DMR_DEG[which(DMR_DEG$site=="SE_W_C" & DMR_DEG$perdiff== 10 & DMR_DEG$random=="rand" & DMR_DEG$context=="CpG"),]
SE_W_C <- as.data.frame(cbind(DMR_SE_W_C$Gene, DMR_SE_W_C$qvalue.x))
SE_W_C <- distinct(SE_W_C)
write.table(SE_W_C, "SE_W_C_randgenes.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
#0

# "SE_HL"
DMR_SE_L_H <- DMR_DEG[which(DMR_DEG$site=="SE_HL" & DMR_DEG$perdiff== 10 & DMR_DEG$random=="rand" & DMR_DEG$context=="CpG"),]
SE_L_H <- as.data.frame(cbind(DMR_SE_L_H$Gene, DMR_SE_L_H$qvalue.x))
SE_L_H <- distinct(SE_L_H)
write.table(SE_L_H, "SE_L_H_randgenes.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
#8

#-----------------------------------
# Extract for metilene

# "ALAS_W_C"
DMR_Alaska_W_C <- DMR_DEG[which(DMR_DEG$site=="ALAS_W_C" & DMR_DEG$random=="non-rand" & DMR_DEG$context=="CpG" & DMR_DEG$program=="metilene"),]
Alaska_W_C <- as.data.frame(cbind(DMR_Alaska_W_C$Gene, DMR_Alaska_W_C$qvalue.y))
Alaska_W_C <- distinct(Alaska_W_C)
write.table(Alaska_W_C, "Alaska_W_C_metilene_genes.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)


# "LAT_W_C"
DMR_Sweden_W_C <- DMR_DEG[which(DMR_DEG$site=="LAT_W_C" & DMR_DEG$random=="non-rand" & DMR_DEG$context=="CpG" & DMR_DEG$program=="metilene"),]
Sweden_W_C <- as.data.frame(cbind(DMR_Sweden_W_C$Gene, DMR_Sweden_W_C$qvalue.y))
Sweden_W_C <- distinct(Sweden_W_C)
write.table(Sweden_W_C, "Sweden_W_C_metilene_genes.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
# 808

# "CASS_W_C"
DMR_Nunavut_W_C <- DMR_DEG[which(DMR_DEG$site=="CASS_W_C" & DMR_DEG$random=="non-rand" & DMR_DEG$context=="CpG" & DMR_DEG$program=="metilene"),]
Nunavut_W_C <- as.data.frame(cbind(DMR_Nunavut_W_C$Gene, DMR_Nunavut_W_C$qvalue.y))
Nunavut_W_C <- distinct(Nunavut_W_C)
write.table(Nunavut_W_C, "Nunavut_W_C_metilene_genes.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
#78

# "SVAL_W_C"
DMR_Svalbard_W_C <- DMR_DEG[which(DMR_DEG$site=="SVAL_W_C" & DMR_DEG$random=="non-rand" & DMR_DEG$context=="CpG" & DMR_DEG$program=="metilene"),]
Svalbard_W_C <- as.data.frame(cbind(DMR_Svalbard_W_C$Gene, DMR_Svalbard_W_C$qvalue.y))
Svalbard_W_C <- distinct(Svalbard_W_C)
write.table(Svalbard_W_C, "Svalbard_W_C_metilene_genes.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
#103

# "Pheno"
DMR_Mat_Sen <- DMR_DEG[which(DMR_DEG$site=="Pheno" & DMR_DEG$random=="non-rand" & DMR_DEG$context=="CpG" & DMR_DEG$program=="metilene"),]
Mat_Sen <- as.data.frame(cbind(DMR_Mat_Sen$Gene, DMR_Mat_Sen$qvalue.y))
Mat_Sen <- distinct(Mat_Sen)
write.table(Mat_Sen, "Mat_Sen_metilene_genes.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
# 8584

# "HL"
DMR_Wild_Lat_L_H <- DMR_DEG[which(DMR_DEG$site=="HL" & DMR_DEG$random=="non-rand" & DMR_DEG$context=="CpG" & DMR_DEG$program=="metilene"),]
Wild_Lat_L_H <- as.data.frame(cbind(DMR_Wild_Lat_L_H$Gene, DMR_Wild_Lat_L_H$qvalue.y))
Wild_Lat_L_H <- distinct(Wild_Lat_L_H)
write.table(Wild_Lat_L_H, "Wild_Lat_L_H_metilene_genes.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
#7621

# "SE_W_C"
DMR_SE_W_C <- DMR_DEG[which(DMR_DEG$site=="SE_W_C" & DMR_DEG$random=="non-rand" & DMR_DEG$context=="CpG" & DMR_DEG$program=="metilene"),]
SE_W_C <- as.data.frame(cbind(DMR_SE_W_C$Gene, DMR_SE_W_C$qvalue.y))
SE_W_C <- distinct(SE_W_C)
write.table(SE_W_C, "SE_W_C_metilene_genes.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
#1434

# "SE_HL"
DMR_SE_L_H <- DMR_DEG[which(DMR_DEG$site=="SE_HL" & DMR_DEG$random=="non-rand" & DMR_DEG$context=="CpG" & DMR_DEG$program=="metilene"),]
SE_L_H <- as.data.frame(cbind(DMR_SE_L_H$Gene, DMR_SE_L_H$qvalue.y))
SE_L_H <- distinct(SE_L_H)
write.table(SE_L_H, "SE_L_H_metilene_genes.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
# 388


#-----------------------------
# Extract for overdispersion

# "ALAS_W_C"
DMR_Alaska_W_C <- DMR_DEG[which(DMR_DEG$site=="ALAS_W_C" & DMR_DEG$perdiff== 10 & DMR_DEG$random=="non-rand" & DMR_DEG$context=="CpG" & DMR_DEG$program=="overdisp"),]
Alaska_W_C <- as.data.frame(cbind(DMR_Alaska_W_C$Gene, DMR_Alaska_W_C$qvalue.x))
Alaska_W_C <- distinct(Alaska_W_C)
write.table(Alaska_W_C, "Alaska_W_C_overdisp_genes.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
#142

# "LAT_W_C"
DMR_Sweden_W_C <- DMR_DEG[which(DMR_DEG$site=="LAT_W_C" & DMR_DEG$perdiff== 10 & DMR_DEG$random=="non-rand" & DMR_DEG$context=="CpG" & DMR_DEG$program=="overdisp"),]
Sweden_W_C <- as.data.frame(cbind(DMR_Sweden_W_C$Gene, DMR_Sweden_W_C$qvalue.x))
Sweden_W_C <- distinct(Sweden_W_C)
write.table(Sweden_W_C, "Sweden_W_C_overdisp_genes.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
# 808

# "CASS_W_C"
DMR_Nunavut_W_C <- DMR_DEG[which(DMR_DEG$site=="CASS_W_C" & DMR_DEG$perdiff== 10 & DMR_DEG$random=="non-rand" & DMR_DEG$context=="CpG" & DMR_DEG$program=="overdisp"),]
Nunavut_W_C <- as.data.frame(cbind(DMR_Nunavut_W_C$Gene, DMR_Nunavut_W_C$qvalue.x))
Nunavut_W_C <- distinct(Nunavut_W_C)
write.table(Nunavut_W_C, "Nunavut_W_C_overdisp_genes.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
#78

# "SVAL_W_C"
DMR_Svalbard_W_C <- DMR_DEG[which(DMR_DEG$site=="SVAL_W_C" & DMR_DEG$perdiff== 10 & DMR_DEG$random=="non-rand" & DMR_DEG$context=="CpG" & DMR_DEG$program=="overdisp"),]
Svalbard_W_C <- as.data.frame(cbind(DMR_Svalbard_W_C$Gene, DMR_Svalbard_W_C$qvalue.x))
Svalbard_W_C <- distinct(Svalbard_W_C)
write.table(Svalbard_W_C, "Svalbard_W_C_overdisp_genes.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
#103

# "Pheno"
DMR_Mat_Sen <- DMR_DEG[which(DMR_DEG$site=="Pheno" & DMR_DEG$perdiff== 10 & DMR_DEG$random=="non-rand" & DMR_DEG$context=="CpG" & DMR_DEG$program=="overdisp"),]
Mat_Sen <- as.data.frame(cbind(DMR_Mat_Sen$Gene, DMR_Mat_Sen$qvalue.x))
Mat_Sen <- distinct(Mat_Sen)
write.table(Mat_Sen, "Mat_Sen_overdisp_genes.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
# 8584

# "HL"
DMR_Wild_Lat_L_H <- DMR_DEG[which(DMR_DEG$site=="HL" & DMR_DEG$perdiff== 10 & DMR_DEG$random=="non-rand" & DMR_DEG$context=="CpG" & DMR_DEG$program=="overdisp"),]
Wild_Lat_L_H <- as.data.frame(cbind(DMR_Wild_Lat_L_H$Gene, DMR_Wild_Lat_L_H$qvalue.x))
Wild_Lat_L_H <- distinct(Wild_Lat_L_H)
write.table(Wild_Lat_L_H, "Wild_Lat_L_H_overdisp_genes.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
#7621

# "SE_W_C"
DMR_SE_W_C <- DMR_DEG[which(DMR_DEG$site=="SE_W_C" & DMR_DEG$perdiff== 10 & DMR_DEG$random=="non-rand" & DMR_DEG$context=="CpG" & DMR_DEG$program=="overdisp"),]
SE_W_C <- as.data.frame(cbind(DMR_SE_W_C$Gene, DMR_SE_W_C$qvalue.x))
SE_W_C <- distinct(SE_W_C)
write.table(SE_W_C, "SE_W_C_overdisp_genes.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
#1434

# "SE_HL"
DMR_SE_L_H <- DMR_DEG[which(DMR_DEG$site=="SE_HL" & DMR_DEG$perdiff== 10 & DMR_DEG$random=="non-rand" & DMR_DEG$context=="CpG" & DMR_DEG$program=="overdisp"),]
SE_L_H <- as.data.frame(cbind(DMR_SE_L_H$Gene, DMR_SE_L_H$qvalue.x))
SE_L_H <- distinct(SE_L_H)
write.table(SE_L_H, "SE_L_H_overdisp_genes.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
# 388




##############################
# Get summary counts for number of genes in various categories

DMR_DEG_summary <- DMR_DEG %>%
  group_by(site, perdiff,context, RNAsite, random) %>%
  summarize(count = n())

print(DMR_DEG_summary, n=137)

write.table(DMR_DEG_summary, "DMR_DEG_summary_counts.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

####################################
# GO ont data

nrow(gene_ont)
# 14963

length(unique(gene_ont$INTPRO))
# 5759

#-----------------------------
# format for ermineJ
# https://erminej.msl.ubc.ca/help/input-files/
# https://erminej.msl.ubc.ca/help/input-files/gene-annotations/

gene_mappings <- gene_ont %>%
  mutate(gene2 = gene) %>%
  select(gene, gene2, descrip2, GOterm)
# write output for each spp
write.table(gene_mappings, "GO_gene_mappings.ermineJ.txt", sep = "\t", quote = FALSE, row.names = FALSE)

#-------------------
# Prepare gene score files
# https://erminej.msl.ubc.ca/help/input-files/gene-scores/

# Manually read the GFF3 file while ignoring lines starting with '#'
gff_data <- read.table("Dryas_octopetala_H1.gff3", header = FALSE, comment.char = "#", sep = "\t")

gene_list <- gff_data[which(gff_data$V3=="gene"), ]
total_gene_list <-  gsub("ID=", "", unique(gene_list$V9))

write.table(total_gene_list, "total_gene_list.txt", sep = "\t", quote = FALSE, row.names = FALSE)

q()
n

#################################
# in bash
# Prepare gene score files - for non-overdisp, random, RNA
# https://erminej.msl.ubc.ca/help/input-files/gene-scores/

cd /lustre04/scratch/celphin/Dryas/GO_enrichment
rename randgenes.txt rand_genes.txt *

awk 'BEGIN{FS="\t"}{print $1,"1"}' total_gene_list.txt | sort -u > Dryas_blank_geneset

# Alaska_W_C Sweden_W_C Nunavut_W_C Svalbard_W_C Mat_Sen Wild_Lat_L_H SE_W_C SE_L_H Sweden_RNA Alaska_RNA Nunavut_RNA Norway_RNA Seedling_RNA Alaska_W_C_rand Sweden_W_C_rand Nunavut_W_C_rand Svalbard_W_C_rand Mat_Sen_rand Wild_Lat_L_H_rand SE_W_C_rand SE_L_H_rand

for taxon in Alaska_W_C Sweden_W_C Nunavut_W_C Svalbard_W_C Mat_Sen Wild_Lat_L_H SE_W_C SE_L_H \
Sweden_RNA Alaska_RNA Nunavut_RNA Norway_RNA Seedling_RNA \
Alaska_W_C_rand Sweden_W_C_rand Nunavut_W_C_rand Svalbard_W_C_rand \
Mat_Sen_rand Wild_Lat_L_H_rand SE_W_C_rand SE_L_H_rand ; \
do \
echo "$taxon"
cp Dryas_blank_geneset "$taxon"_geneset

# Loop through the gene list for each taxon, extract the qvalue and replace 1 with the qvalue in the gene set file
sort "$taxon"_genes.txt | \
while read gene qvalue ; do \
    # Replace 1 with the qvalue in the gene set file
    sed -i "s/$gene 1/$gene $qvalue/g" "$taxon"_geneset ; \
done ; done

for taxon in Alaska_W_C Sweden_W_C Nunavut_W_C Svalbard_W_C Mat_Sen Wild_Lat_L_H SE_W_C SE_L_H \
Sweden_RNA Alaska_RNA Nunavut_RNA Norway_RNA Seedling_RNA \
Alaska_W_C_rand Sweden_W_C_rand Nunavut_W_C_rand Svalbard_W_C_rand \
Mat_Sen_rand Wild_Lat_L_H_rand SE_W_C_rand SE_L_H_rand ; \
do sed -i 's/ /\t/g' "$taxon"_geneset; done

# remove V1.1 in GO_gene_mappings.ermineJ.txt
sed 's/V1\.1//g' GO_gene_mappings.ermineJ.txt > GO_gene_mappings1.ermineJ.txt

################################
# Prepare gene score files - for non-overdisp, random, RNA
# https://erminej.msl.ubc.ca/help/input-files/gene-scores/

cd /lustre04/scratch/celphin/Dryas/GO_enrichment

# Alaska_W_C_overdisp Sweden_W_C_overdisp Nunavut_W_C_overdisp Svalbard_W_C_overdisp Mat_Sen_overdisp Wild_Lat_L_H_overdisp SE_W_C_overdisp SE_L_H_overdisp Alaska_W_C_metilene Sweden_W_C_metilene Nunavut_W_C_metilene Svalbard_W_C_metilene Mat_Sen_metilene Wild_Lat_L_H_metilene SE_W_C_metilene SE_L_H_metilene

for taxon in Alaska_W_C_overdisp Sweden_W_C_overdisp Nunavut_W_C_overdisp \
Svalbard_W_C_overdisp Mat_Sen_overdisp Wild_Lat_L_H_overdisp SE_W_C_overdisp \
SE_L_H_overdisp Alaska_W_C_metilene Sweden_W_C_metilene Nunavut_W_C_metilene \
Svalbard_W_C_metilene Mat_Sen_metilene Wild_Lat_L_H_metilene SE_W_C_metilene SE_L_H_metilene ; \
do \
echo "$taxon"
cp Dryas_blank_geneset "$taxon"_geneset

# Loop through the gene list for each taxon, extract the qvalue and replace 1 with the qvalue in the gene set file
sort "$taxon"_genes.txt | \
while read gene qvalue ; do \
    # Replace 1 with the qvalue in the gene set file
    sed -i "s/$gene 1/$gene $qvalue/g" "$taxon"_geneset ; \
done ; done

for taxon in Alaska_W_C_overdisp Sweden_W_C_overdisp Nunavut_W_C_overdisp \
Svalbard_W_C_overdisp Mat_Sen_overdisp Wild_Lat_L_H_overdisp SE_W_C_overdisp \
SE_L_H_overdisp Alaska_W_C_metilene Sweden_W_C_metilene Nunavut_W_C_metilene \
Svalbard_W_C_metilene Mat_Sen_metilene Wild_Lat_L_H_metilene SE_W_C_metilene SE_L_H_metilene ; \
do sed -i 's/ /\t/g' "$taxon"_geneset; done

# remove V1.1 in GO_gene_mappings.ermineJ.txt
sed 's/V1\.1//g' GO_gene_mappings.ermineJ.txt > GO_gene_mappings1.ermineJ.txt


############################################
# get counts
cd /lustre04/scratch/celphin/Dryas/GO_enrichment
wc -l *_genes.txt

   # 1456 Alaska_W_C_genes.txt
   # 9339 Mat_Sen_genes.txt
    # 555 Nunavut_W_C_genes.txt
    # 753 Parent_W_C_genes.txt
   # 2599 RNA_genes.txt
    # 653 SE_L_H_genes.txt
   # 1282 SE_W_C_genes.txt
    # 531 Svalbard_W_C_genes.txt
   # 2902 Sweden_W_C_genes.txt
  # 31161 Wild_Lat_L_H_genes.txt
    # 573 Wild_W_C_genes.txt

     # 19 Alaska_RNA_genes.txt
    # 150 Alaska_W_C_genes.txt
     # 13 Alaska_W_C_rand_genes.txt
  # 14278 Mat_Sen_genes.txt
      # 0 Mat_Sen_rand_genes.txt
      # 6 Norway_RNA_genes.txt
     # 21 Nunavut_RNA_genes.txt
     # 84 Nunavut_W_C_genes.txt
      # 0 Nunavut_W_C_rand_genes.txt
     # 45 Seedling_RNA_genes.txt
    # 405 SE_L_H_genes.txt
      # 8 SE_L_H_rand_genes.txt
   # 1631 SE_W_C_genes.txt
      # 0 SE_W_C_rand_genes.txt
    # 103 Svalbard_W_C_genes.txt
    # 123 Svalbard_W_C_rand_genes.txt
     # 59 Sweden_RNA_genes.txt
    # 890 Sweden_W_C_genes.txt
      # 0 Sweden_W_C_rand_genes.txt
  # 12853 Wild_Lat_L_H_genes.txt
      # 0 Wild_Lat_L_H_rand_genes.txt


##############################################
# ermineJ

# https://erminej.msl.ubc.ca/help/tutorials/running-an-analysis-ora/

# As of ErmineJ 3, when using the ‘ORA’ method you have the option to use a simple “hit list” of genes,
# rather than preparing a score file yourself (a “quick list”). Caution: If you use this feature, 
# the “non-hits” will be all the rest of the genes listed in your annotation file. That might not 
# be appropriate if the annotation file includes genes that were not assayed in your experiment. 
# This is most likely to be a problem if your annotation file is a list of all the genes in the genome

#In home, install ermineJ:
cd ~
wget https://home.pavlab.msl.ubc.ca/ermineJ/distributions/ermineJ-3.2-generic-bundle.zip
unzip ermineJ-3.2-generic-bundle.zip
module load StdEnv/2020  python/3.10.2 scipy-stack/2021a diamond/2.1.7 fastme/2.1.6.2 gcc/9.3.0 mcl/14.137 java/13.0.2 r/4.2.2 glpk/5.0
cd ermineJ-3.2/bin/
chmod +x ermineJ.sh

#Get gene ontology file (gene set file)
mkdir ~/ermineJ.data
cd ~/ermineJ.data
wget https://release.geneontology.org/2024-04-24/ontology/go.obo

#######################################################################################################################################

tmux new-session -s GO1
tmux attach-session -t GO1 #methwarm

salloc -c1 --time 3:00:00 --mem 120000m --account def-henryg

cd /lustre04/scratch/celphin/Dryas/GO_enrichment

ERMINEJ_HOME=/home/celphin/ermineJ-3.2
export JAVA_HOME=/cvmfs/soft.computecanada.ca/easybuild/software/2020/Core/java/13.0.2/

module load StdEnv/2020 
module load java/13.0.2

for taxon in Alaska_W_C Sweden_W_C Nunavut_W_C Svalbard_W_C Mat_Sen Wild_Lat_L_H SE_W_C SE_L_H \
Sweden_RNA Alaska_RNA Nunavut_RNA Norway_RNA Seedling_RNA \
Alaska_W_C_rand Sweden_W_C_rand Nunavut_W_C_rand Svalbard_W_C_rand \
Mat_Sen_rand Wild_Lat_L_H_rand SE_W_C_rand SE_L_H_rand ; \
do \
echo $taxon 
$ERMINEJ_HOME/bin/ermineJ.sh \
-a GO_gene_mappings1.ermineJ.txt \
-s "$taxon"_geneset \
-c /home/celphin/ermineJ.data/go.obo \
-n ORA -t 0.01 \
--genesOut -aspects BCM \
-o "$taxon".ermine.results -y 5 ; done


# Starting analysis ...
# Reading scores from /lustre04/scratch/celphin/Dryas/GO_enrichment/SE_L_H_rand_geneset
# Reading gene scores from column 2 ...
# 14960 (36.33%) of the scores were usable (others may not have genes in the annotations?)
# 3 elements in your gene score file had no gene sets and were ignored.
# Usable scores for 14960 distinct genes found (99.98%)
# Creating a subsetted annotation set for 14960/14963 elements) ...
# INFO: Subclone: 1507ms
# INFO: Subclone annotations: 1514ms
# Starting ORA analysis
# Hit list (8 genes) enrichment for multifunctionality: P = 0.495
# 300 gene sets analyzed ...
# 600 gene sets analyzed ...
# 900 gene sets analyzed ...
# 1200 gene sets analyzed ...
# 1500 gene sets analyzed ...
# 'Hits' are not significantly multifunctionality-biased, no multifunctionality correction needed
# Finished with ORA computations: 8 elements passed your threshold.

#---------------------
# Try combining some of the genesets
cat Alaska_W_C_geneset Sweden_W_C_geneset | sort -u >  Alaska_Sweden_W_C_geneset
cat Svalbard_W_C_geneset Nunavut_W_C_geneset | sort -u >  Nunavut_Svalbard_W_C_geneset
cat Alaska_RNA_geneset Nunavut_RNA_geneset Norway_RNA_geneset | sort -u >  Total_RNA_geneset
cat Alaska_W_C_rand_geneset Sweden_W_C_rand_geneset Nunavut_W_C_rand_geneset Svalbard_W_C_rand_geneset \
Mat_Sen_rand_geneset Wild_Lat_L_H_rand_geneset SE_W_C_rand_geneset SE_L_H_rand_geneset | sort -u > Total_random_geneset

# remove duplicates and select the lowest value

for taxon in Alaska_Sweden_W_C Nunavut_Svalbard_W_C Total_RNA Total_random ; \
do \
echo $taxon 
$ERMINEJ_HOME/bin/ermineJ.sh \
-a GO_gene_mappings1.ermineJ.txt \
-s "$taxon"_geneset \
-c /home/celphin/ermineJ.data/go.obo \
-n ORA -t 0.01 \
--genesOut -aspects BCM \
-o "$taxon".ermine.results -y 5 ; done

# Picked up JAVA_TOOL_OPTIONS: -Xmx2g
# INFO: Data directory is /home/celphin/ermineJ.data
# DEBUG: Custom gene sets directory is /home/celphin/ermineJ.data/genesets
# INFO: Gene symbols for each term will be output
# DEBUG: gene score threshold set to 0.01
# Reading GO descriptions from /home/celphin/ermineJ.data/go.obo ...
# Reading gene annotations from /lustre04/scratch/celphin/Dryas/GO_enrichment/GO_gene_mappings1.e
# rmineJ.txt ...
# Read 2500 elements ...
# Read 5000 elements ...
# Read 7500 elements ...
# Read 10000 elements ...
# Read 12500 elements ...
# No gene sets found in /home/celphin/ermineJ.data/genesets
# Inferring annotations in graph ...
# 3000 genes examined for term parents ... ...
# 6000 genes examined for term parents ... ...
# 9000 genes examined for term parents ... ...
# 12000 genes examined for term parents ... ...
# Added 194927 inferred annotations (affected 14836/14963 genes) ...
# Pruning: 921/3985 sets removed: obsolete (0), too small (921) or too big (0) terms pruned. ...
# There are 3064 gene sets in the annotations, checking for redundancy ... ...
# 2000 sets checked for redundancy, 757 found ... ...
# 1455/3064 gene sets are redundant with at least one other. ...
# INFO: Redundancy check: 161ms
# INFO: Multifunctionality computation: 1143ms
# INFO: Total annotation setup: 1823ms
# Initializing gene class mapping ...
# Done with setup ...
# Ready.

# Starting analysis ...
# Reading scores from /lustre04/scratch/celphin/Dryas/GO_enrichment/Total_random_geneset
# Reading gene scores from column 2 ...
# Repeated identifier: Do1_01_a00001G00168, keeping original value.
# 14960 (36.33%) of the scores were usable (others may not have genes in the annotations?)
# 3 elements in your gene score file had no gene sets and were ignored.
# 144 identifiers in your gene score file were repeats. Only the first occurrence encountered was
 # kept in each case.
# Usable scores for 14960 distinct genes found (99.98%)
# Creating a subsetted annotation set for 14960/14963 elements) ...
# INFO: Multifunctionality computation: 1185ms
# INFO: Subclone: 1539ms
# INFO: Subclone annotations: 1544ms
# Starting ORA analysis
# Hit list (117 genes) enrichment for multifunctionality: P = 0.00448
# 300 gene sets analyzed ...
# 600 gene sets analyzed ...
# 900 gene sets analyzed ...
# 1200 gene sets analyzed ...
# 1500 gene sets analyzed ...
# 0 top groups will be monitored for multifunctionality sensitivity
# Insufficient enrichment found, skipping multifunctionality correction
# Finished with ORA computations: 117 elements passed your threshold.
# Multiple test correction for 1667 scored sets.
# Multifunctionality correlation is 0.02 for 14960 values
# Done!

# Alaska_Sweden_W_C
# Finished with ORA computations: 815 elements passed your threshold.

# Nunavut_Svalbard_W_C
# Finished with ORA computations: 152 elements passed your threshold.

# Total_RNA
# Finished with ORA computations: 2 elements passed your threshold.

#Total_random
# Finished with ORA computations: 117 elements passed your threshold.

#---------------------------------
# Run for metilene and overdisp data


for taxon in Alaska_W_C_overdisp Sweden_W_C_overdisp Nunavut_W_C_overdisp \
Svalbard_W_C_overdisp Mat_Sen_overdisp Wild_Lat_L_H_overdisp SE_W_C_overdisp \
SE_L_H_overdisp Alaska_W_C_metilene Sweden_W_C_metilene Nunavut_W_C_metilene \
Svalbard_W_C_metilene Mat_Sen_metilene Wild_Lat_L_H_metilene SE_W_C_metilene SE_L_H_metilene ; \
do \
echo $taxon 
$ERMINEJ_HOME/bin/ermineJ.sh \
-a GO_gene_mappings1.ermineJ.txt \
-s "$taxon"_geneset \
-c /home/celphin/ermineJ.data/go.obo \
-n ORA -t 0.01 \
--genesOut -aspects BCM \
-o "$taxon".ermine.results -y 5 ; done





############################
# Results exploration

# Look at the files for each set combined
# Subset those with a p-value < 0.05
for taxon in Alaska_W_C Sweden_W_C Nunavut_W_C Svalbard_W_C Mat_Sen Wild_Lat_L_H SE_W_C SE_L_H \
Sweden_RNA Alaska_RNA Nunavut_RNA Norway_RNA Seedling_RNA \
Alaska_W_C_rand Sweden_W_C_rand Nunavut_W_C_rand Svalbard_W_C_rand \
Mat_Sen_rand Wild_Lat_L_H_rand SE_W_C_rand SE_L_H_rand ; \
do grep "!" "$taxon".ermine.results  | awk -F '\t' '$7 <0.5 { print }'| awk  -F $'\t' '{print $2, $3, $7, $1}' |sort  > "$taxon"_sig.ermine.results
done

# Extract just GO terms and count duplicates
cat *_sig.ermine.results | awk -F '\t' '$7 <0.5 { print }' | awk  -F $'\t' '{print $3, $2}' | sort | uniq -c  
#  4948

mkdir sig.ermine.results
mv *sig.ermine.results sig.ermine.results

mkdir ermine.results
mv *ermine.results ermine.results

mkdir Revigio
mv Revigio* Revigio
#####################
# look at all the data combined

cd /lustre04/scratch/celphin/Dryas/GO_enrichment/ermine.results


for file in *.ermine.results; do
  tail -n +27 "$file" | awk -F'\t' -v fname="$file" '$6 > 0 {print $3 "\t" $2 "\t" fname "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 }'
done | sort -t $'\t' -k1,1 > combined_sorted_file.txt

awk -F'\t' '{print $1 "\t" $2}' combined_sorted_file.txt | sort | uniq -c | sort -n > value_counts.txt

grep "trichome morphogenesis" combined_sorted_file.txt

grep "response to cold" combined_sorted_file.txt

grep "response to water" combined_sorted_file.txt

grep "oxidoreductase activity" combined_sorted_file.txt

grep "shoot system development" combined_sorted_file.txt

grep "system development" combined_sorted_file.txt

grep "defense response to other organism" combined_sorted_file.txt

grep "light" combined_sorted_file.txt

grep "positive regulation" combined_sorted_file.txt


GO:0098542      defense response to other organism      Alaska_Sweden_W_C.ermine.results        100     100      4.000   0.80117054      1.000
GO:0098542      defense response to other organism      Alaska_W_C.ermine.results               100     100      2.000    0.24291749      1.000
GO:0098542      defense response to other organism      Alaska_W_C_metilene.ermine.results      100     100      9.000   0.01870974      1.000
GO:0098542      defense response to other organism      Alaska_W_C_overdisp.ermine.results      93      93       22.000  0.36690591      1.000
GO:0098542      defense response to other organism      Mat_Sen.ermine.results                  91      91      53.000      0.23333652      0.67203156
GO:0098542      defense response to other organism      Mat_Sen_metilene.ermine.results          93      93     5.000    0.7328553       1.000
GO:0098542      defense response to other organism      Nunavut_RNA.ermine.results             100     100      1.000    0.13145871      1.000
GO:0098542      defense response to other organism      Nunavut_Svalbard_W_C.erm ine.results     100    100      1.000   0.64107842      1.000
GO:0098542      defense response to other organism      SE_L_H.ermine.results                  100     100     2.000       0.73517512      1.000
GO:0098542      defense response to other organism      SE_L_H_metilene.ermine.results         100     100    1.000    0.90646637      1.000
GO:0098542      defense response to other organism      SE_W_C.ermine.results                 100     100     5.000       0.96810934      1.000
GO:0098542      defense response to other organism      SE_W_C_metilene.ermine.results         100     100      5.000    0.4779764       1.000
GO:0098542      defense response to other organism      Svalbard_W_C.ermine.results             100     100     1.000    0.5000204       1.000
GO:0098542      defense response to other organism      Svalbard_W_C_metilene.ermine.results    100     100      1.000   0.81457363      1.000
GO:0098542      defense response to other organism      Svalbard_W_C_overdisp.ermine.results    93      93       30.000  0.17914744      0.82941163
GO:0098542      defense response to other organism      Sweden_W_C.ermine.results              100     100      3.000    0.91171961      1.000
GO:0098542      defense response to other organism      Sweden_W_C_metilene.ermine.results      99     99       11.000  0.09357311      1.000
GO:0098542      defense response to other organism      Sweden_W_C_overdisp.ermine.results      92     92       32.000  0.09128341      0.82021804
GO:0098542      defense response to other organism      Wild_Lat_L_H.ermine.results            93      93       35.000   0.96103226      1.000
GO:0098542      defense response to other organism      Wild_Lat_L_H_metilene.ermine.results    88     88       36.000  0.73912346      1.000



###################################
# summarize with : http://revigo.irb.hr/

for taxon in Alaska_W_C Sweden_W_C Nunavut_W_C Svalbard_W_C Mat_Sen Wild_Lat_L_H SE_W_C SE_L_H \
Sweden_RNA Alaska_RNA Nunavut_RNA Norway_RNA Seedling_RNA \
Alaska_W_C_rand Sweden_W_C_rand Nunavut_W_C_rand Svalbard_W_C_rand \
Mat_Sen_rand Wild_Lat_L_H_rand SE_W_C_rand SE_L_H_rand ; \
do cat "$taxon".ermine.results | awk -F '\t' '$7 <0.05 { print }' | awk  -F $'\t' '{print $3}' | \
sort  | uniq -c | awk  -F " " '{print $2, $1}' > Revigio_"$taxon".txt 
done

###################################
# Look at all and what species they are from
cd /lustre04/scratch/celphin/Dryas/GO_enrichment

grep "!" *sig.ermine.results  | awk -F '\t' '$7 <1 { print }'| awk  -F $'\t' '{print $1, $2, $3, $7}' |sort  > List_all_GO.txt


############################
# Defense genes

# Alaska 
response to fungus 12
Do1_01_a00001G01727|Do1_01_a00001G01728|Do1_01_a00001G01733|Do1_01_a00001G01734|Do1_01_a00001G01735|Do1_01_a00001G01741|Do1_01_a00001G01742|Do1_01_a00001G01743|Do1_01_a00001G01744|Do1_01_a00001G01746|Do1_01_a00001G01754|Do1_01_a00001G01755|Do1_01_a00001G01757|Do1_01_a00001G01760|Do1_01_a00001G01761|Do1_01_a00001G01765|Do1_01_a00001G01770|Do1_01_a00001G02267|Do1_01_a00001G02274|Do1_01_a00003G00001|Do1_02_a00003G02060|Do1_02_a00003G02106|Do1_04_a00001G00056|Do1_04_a00004G00044|Do1_04_a00004G00088|Do1_05_a00001G00574|Do1_05_a00001G01355|Do1_05_a00001G01358|Do1_05_a00001G01360|Do1_05_a00001G01361|Do1_05_a00001G01362|Do1_05_a00001G01363|Do1_06_a00001G00217|Do1_06_a00002G00149|Do1_06_a00002G00150|Do1_06_a00002G01992|Do1_06_a00002G02295|Do1_06_a00002G02298|Do1_06_a00002G02312|Do1_07_a00001G00124|Do1_07_a00001G00135|Do1_07_a00001G00268|Do1_07_a00001G00512|Do1_07_a00001G00596|Do1_07_a00002G02049|Do1_07_a00004G00205|Do1_07_a00004G01058|

defense response to fungus 12
Do1_01_a00001G01727|Do1_01_a00001G01728|Do1_01_a00001G01733|Do1_01_a00001G01734|Do1_01_a00001G01735|Do1_01_a00001G01741|Do1_01_a00001G01742|Do1_01_a00001G01743|Do1_01_a00001G01744|Do1_01_a00001G01746|Do1_01_a00001G01754|Do1_01_a00001G01755|Do1_01_a00001G01757|Do1_01_a00001G01760|Do1_01_a00001G01761|Do1_01_a00001G01765|Do1_01_a00001G01770|Do1_01_a00001G02267|Do1_01_a00001G02274|Do1_01_a00003G00001|Do1_02_a00003G02060|Do1_02_a00003G02106|Do1_04_a00001G00056|Do1_04_a00004G00044|Do1_04_a00004G00088|Do1_05_a00001G00574|Do1_05_a00001G01355|Do1_05_a00001G01358|Do1_05_a00001G01360|Do1_05_a00001G01361|Do1_05_a00001G01362|Do1_05_a00001G01363|Do1_06_a00001G00217|Do1_06_a00002G00149|Do1_06_a00002G00150|Do1_06_a00002G01992|Do1_06_a00002G02295|Do1_06_a00002G02298|Do1_06_a00002G02312|Do1_07_a00001G00124|Do1_07_a00001G00135|Do1_07_a00001G00268|Do1_07_a00001G00512|Do1_07_a00001G00596|Do1_07_a00002G02049|Do1_07_a00004G00205|Do1_07_a00004G01058|

response to external biotic stimulus 10
Do1_01_a00001G01727|Do1_01_a00001G01728|Do1_01_a00001G01733|Do1_01_a00001G01734|Do1_01_a00001G01735|Do1_01_a00001G01741|Do1_01_a00001G01742|Do1_01_a00001G01743|Do1_01_a00001G01744|Do1_01_a00001G01746|Do1_01_a00001G01754|Do1_01_a00001G01755|Do1_01_a00001G01757|Do1_01_a00001G01760|Do1_01_a00001G01761|Do1_01_a00001G01765|Do1_01_a00001G01770|Do1_01_a00001G01935|Do1_01_a00001G02267|Do1_01_a00001G02274|Do1_01_a00001G02636|Do1_01_a00003G00001|Do1_02_a00001G00907|Do1_02_a00003G00790|Do1_02_a00003G01423|Do1_02_a00003G01427|Do1_02_a00003G01718|Do1_02_a00003G02060|Do1_02_a00003G02106|Do1_02_a00004G00421|Do1_02_a00004G00422|Do1_02_a00004G01171|Do1_03_a00002G01111|Do1_03_a00002G01327|Do1_03_a00002G01395|Do1_03_a00002G01466|Do1_03_a00002G01589|Do1_03_a00002G01590|Do1_03_a00002G01591|Do1_03_a00002G02055|Do1_03_a00002G02474|Do1_03_a00003G00250|Do1_04_a00001G00056|Do1_04_a00001G00217|Do1_04_a00001G00303|Do1_04_a00001G00314|Do1_04_a00001G00315|Do1_04_a00001G00316|Do1_04_a00001G00668|Do1_04_a00001G00672|Do1_04_a00001G00716|Do1_04_a00001G00992|Do1_04_a00001G02322|Do1_04_a00001G03227|Do1_04_a00004G00044|Do1_04_a00004G00088|Do1_05_a00001G00574|Do1_05_a00001G01108|Do1_05_a00001G01355|Do1_05_a00001G01358|Do1_05_a00001G01360|Do1_05_a00001G01361|Do1_05_a00001G01362|Do1_05_a00001G01363|Do1_05_a00003G00011|Do1_05_a00003G00201|Do1_06_a00001G00217|Do1_06_a00001G01034|Do1_06_a00001G02027|Do1_06_a00001G02700|Do1_06_a00002G00004|Do1_06_a00002G00149|Do1_06_a00002G00150|Do1_06_a00002G00187|Do1_06_a00002G00731|Do1_06_a00002G00735|Do1_06_a00002G00996|Do1_06_a00002G01629|Do1_06_a00002G01635|Do1_06_a00002G01637|Do1_06_a00002G01651|Do1_06_a00002G01784|Do1_06_a00002G01992|Do1_06_a00002G02295|Do1_06_a00002G02298|Do1_06_a00002G02312|Do1_07_a00001G00124|Do1_07_a00001G00135|Do1_07_a00001G00268|Do1_07_a00001G00512|Do1_07_a00001G00529|Do1_07_a00001G00596|Do1_07_a00002G00835|Do1_07_a00002G00836|Do1_07_a00002G02049|Do1_07_a00004G00205|Do1_07_a00004G01044|Do1_07_a00004G01046|Do1_07_a00004G01058|Do1_07_a00006G00136|

response to other organism 10
Do1_01_a00001G01727|Do1_01_a00001G01728|Do1_01_a00001G01733|Do1_01_a00001G01734|Do1_01_a00001G01735|Do1_01_a00001G01741|Do1_01_a00001G01742|Do1_01_a00001G01743|Do1_01_a00001G01744|Do1_01_a00001G01746|Do1_01_a00001G01754|Do1_01_a00001G01755|Do1_01_a00001G01757|Do1_01_a00001G01760|Do1_01_a00001G01761|Do1_01_a00001G01765|Do1_01_a00001G01770|Do1_01_a00001G01935|Do1_01_a00001G02267|Do1_01_a00001G02274|Do1_01_a00001G02636|Do1_01_a00003G00001|Do1_02_a00001G00907|Do1_02_a00003G00790|Do1_02_a00003G01423|Do1_02_a00003G01427|Do1_02_a00003G01718|Do1_02_a00003G02060|Do1_02_a00003G02106|Do1_02_a00004G00421|Do1_02_a00004G00422|Do1_02_a00004G01171|Do1_03_a00002G01111|Do1_03_a00002G01327|Do1_03_a00002G01395|Do1_03_a00002G01466|Do1_03_a00002G01589|Do1_03_a00002G01590|Do1_03_a00002G01591|Do1_03_a00002G02055|Do1_03_a00002G02474|Do1_03_a00003G00250|Do1_04_a00001G00056|Do1_04_a00001G00217|Do1_04_a00001G00303|Do1_04_a00001G00314|Do1_04_a00001G00315|Do1_04_a00001G00316|Do1_04_a00001G00668|Do1_04_a00001G00672|Do1_04_a00001G00716|Do1_04_a00001G00992|Do1_04_a00001G02322|Do1_04_a00001G03227|Do1_04_a00004G00044|Do1_04_a00004G00088|Do1_05_a00001G00574|Do1_05_a00001G01108|Do1_05_a00001G01355|Do1_05_a00001G01358|Do1_05_a00001G01360|Do1_05_a00001G01361|Do1_05_a00001G01362|Do1_05_a00001G01363|Do1_05_a00003G00011|Do1_05_a00003G00201|Do1_06_a00001G00217|Do1_06_a00001G01034|Do1_06_a00001G02027|Do1_06_a00001G02700|Do1_06_a00002G00004|Do1_06_a00002G00149|Do1_06_a00002G00150|Do1_06_a00002G00187|Do1_06_a00002G00731|Do1_06_a00002G00735|Do1_06_a00002G00996|Do1_06_a00002G01629|Do1_06_a00002G01635|Do1_06_a00002G01637|Do1_06_a00002G01651|Do1_06_a00002G01784|Do1_06_a00002G01992|Do1_06_a00002G02295|Do1_06_a00002G02298|Do1_06_a00002G02312|Do1_07_a00001G00124|Do1_07_a00001G00135|Do1_07_a00001G00268|Do1_07_a00001G00512|Do1_07_a00001G00529|Do1_07_a00001G00596|Do1_07_a00002G00835|Do1_07_a00002G00836|Do1_07_a00002G02049|Do1_07_a00004G00205|Do1_07_a00004G01044|Do1_07_a00004G01046|Do1_07_a00004G01058|Do1_07_a00006G00136|

defense response to other organism 12
Do1_01_a00001G01727|Do1_01_a00001G01728|Do1_01_a00001G01733|Do1_01_a00001G01734|Do1_01_a00001G01735|Do1_01_a00001G01741|Do1_01_a00001G01742|Do1_01_a00001G01743|Do1_01_a00001G01744|Do1_01_a00001G01746|Do1_01_a00001G01754|Do1_01_a00001G01755|Do1_01_a00001G01757|Do1_01_a00001G01760|Do1_01_a00001G01761|Do1_01_a00001G01765|Do1_01_a00001G01770|Do1_01_a00001G01935|Do1_01_a00001G02267|Do1_01_a00001G02274|Do1_01_a00001G02636|Do1_01_a00003G00001|Do1_02_a00001G00907|Do1_02_a00003G00790|Do1_02_a00003G01423|Do1_02_a00003G01427|Do1_02_a00003G01718|Do1_02_a00003G02060|Do1_02_a00003G02106|Do1_02_a00004G00421|Do1_02_a00004G00422|Do1_02_a00004G01171|Do1_03_a00002G01111|Do1_03_a00002G01327|Do1_03_a00002G01395|Do1_03_a00002G01466|Do1_03_a00002G01589|Do1_03_a00002G01590|Do1_03_a00002G01591|Do1_03_a00002G02055|Do1_03_a00002G02474|Do1_03_a00003G00250|Do1_04_a00001G00056|Do1_04_a00001G00217|Do1_04_a00001G00303|Do1_04_a00001G00314|Do1_04_a00001G00315|Do1_04_a00001G00316|Do1_04_a00001G00668|Do1_04_a00001G00672|Do1_04_a00001G00716|Do1_04_a00001G00992|Do1_04_a00001G02322|Do1_04_a00001G03227|Do1_04_a00004G00044|Do1_04_a00004G00088|Do1_05_a00001G00574|Do1_05_a00001G01108|Do1_05_a00001G01355|Do1_05_a00001G01358|Do1_05_a00001G01360|Do1_05_a00001G01361|Do1_05_a00001G01362|Do1_05_a00001G01363|Do1_05_a00003G00011|Do1_05_a00003G00201|Do1_06_a00001G00217|Do1_06_a00001G01034|Do1_06_a00001G02027|Do1_06_a00001G02700|Do1_06_a00002G00004|Do1_06_a00002G00149|Do1_06_a00002G00150|Do1_06_a00002G00187|Do1_06_a00002G00731|Do1_06_a00002G00735|Do1_06_a00002G00996|Do1_06_a00002G01629|Do1_06_a00002G01635|Do1_06_a00002G01637|Do1_06_a00002G01651|Do1_06_a00002G01784|Do1_06_a00002G01992|Do1_06_a00002G02295|Do1_06_a00002G02298|Do1_06_a00002G02312|Do1_07_a00001G00124|Do1_07_a00001G00135|Do1_07_a00001G00268|Do1_07_a00001G00512|Do1_07_a00001G00529|Do1_07_a00001G00596|Do1_07_a00002G00835|Do1_07_a00002G00836|Do1_07_a00002G02049|Do1_07_a00004G00205|Do1_07_a00004G01044|Do1_07_a00004G01046|Do1_07_a00004G01058|Do1_07_a00006G00136|

#------------
#Seedling

response to external biotic stimulus 12
Do1_01_a00001G01727|Do1_01_a00001G01728|Do1_01_a00001G01733|Do1_01_a00001G01734|Do1_01_a00001G01735|Do1_01_a00001G01741|Do1_01_a00001G01742|Do1_01_a00001G01743|Do1_01_a00001G01744|Do1_01_a00001G01746|Do1_01_a00001G01754|Do1_01_a00001G01755|Do1_01_a00001G01757|Do1_01_a00001G01760|Do1_01_a00001G01761|Do1_01_a00001G01765|Do1_01_a00001G01770|Do1_01_a00001G01935|Do1_01_a00001G02267|Do1_01_a00001G02274|Do1_01_a00001G02636|Do1_01_a00003G00001|Do1_02_a00001G00907|Do1_02_a00003G00790|Do1_02_a00003G01423|Do1_02_a00003G01427|Do1_02_a00003G01718|Do1_02_a00003G02060|Do1_02_a00003G02106|Do1_02_a00004G00421|Do1_02_a00004G00422|Do1_02_a00004G01171|Do1_03_a00002G01111|Do1_03_a00002G01327|Do1_03_a00002G01395|Do1_03_a00002G01466|Do1_03_a00002G01589|Do1_03_a00002G01590|Do1_03_a00002G01591|Do1_03_a00002G02055|Do1_03_a00002G02474|Do1_03_a00003G00250|Do1_04_a00001G00056|Do1_04_a00001G00217|Do1_04_a00001G00303|Do1_04_a00001G00314|Do1_04_a00001G00315|Do1_04_a00001G00316|Do1_04_a00001G00668|Do1_04_a00001G00672|Do1_04_a00001G00716|Do1_04_a00001G00992|Do1_04_a00001G02322|Do1_04_a00001G03227|Do1_04_a00004G00044|Do1_04_a00004G00088|Do1_05_a00001G00574|Do1_05_a00001G01108|Do1_05_a00001G01355|Do1_05_a00001G01358|Do1_05_a00001G01360|Do1_05_a00001G01361|Do1_05_a00001G01362|Do1_05_a00001G01363|Do1_05_a00003G00011|Do1_05_a00003G00201|Do1_06_a00001G00217|Do1_06_a00001G01034|Do1_06_a00001G02027|Do1_06_a00001G02700|Do1_06_a00002G00004|Do1_06_a00002G00149|Do1_06_a00002G00150|Do1_06_a00002G00187|Do1_06_a00002G00731|Do1_06_a00002G00735|Do1_06_a00002G00996|Do1_06_a00002G01629|Do1_06_a00002G01635|Do1_06_a00002G01637|Do1_06_a00002G01651|Do1_06_a00002G01784|Do1_06_a00002G01992|Do1_06_a00002G02295|Do1_06_a00002G02298|Do1_06_a00002G02312|Do1_07_a00001G00124|Do1_07_a00001G00135|Do1_07_a00001G00268|Do1_07_a00001G00512|Do1_07_a00001G00529|Do1_07_a00001G00596|Do1_07_a00002G00835|Do1_07_a00002G00836|Do1_07_a00002G02049|Do1_07_a00004G00205|Do1_07_a00004G01044|Do1_07_a00004G01046|Do1_07_a00004G01058|Do1_07_a00006G00136|

response to other organism 12
Do1_01_a00001G01727|Do1_01_a00001G01728|Do1_01_a00001G01733|Do1_01_a00001G01734|Do1_01_a00001G01735|Do1_01_a00001G01741|Do1_01_a00001G01742|Do1_01_a00001G01743|Do1_01_a00001G01744|Do1_01_a00001G01746|Do1_01_a00001G01754|Do1_01_a00001G01755|Do1_01_a00001G01757|Do1_01_a00001G01760|Do1_01_a00001G01761|Do1_01_a00001G01765|Do1_01_a00001G01770|Do1_01_a00001G01935|Do1_01_a00001G02267|Do1_01_a00001G02274|Do1_01_a00001G02636|Do1_01_a00003G00001|Do1_02_a00001G00907|Do1_02_a00003G00790|Do1_02_a00003G01423|Do1_02_a00003G01427|Do1_02_a00003G01718|Do1_02_a00003G02060|Do1_02_a00003G02106|Do1_02_a00004G00421|Do1_02_a00004G00422|Do1_02_a00004G01171|Do1_03_a00002G01111|Do1_03_a00002G01327|Do1_03_a00002G01395|Do1_03_a00002G01466|Do1_03_a00002G01589|Do1_03_a00002G01590|Do1_03_a00002G01591|Do1_03_a00002G02055|Do1_03_a00002G02474|Do1_03_a00003G00250|Do1_04_a00001G00056|Do1_04_a00001G00217|Do1_04_a00001G00303|Do1_04_a00001G00314|Do1_04_a00001G00315|Do1_04_a00001G00316|Do1_04_a00001G00668|Do1_04_a00001G00672|Do1_04_a00001G00716|Do1_04_a00001G00992|Do1_04_a00001G02322|Do1_04_a00001G03227|Do1_04_a00004G00044|Do1_04_a00004G00088|Do1_05_a00001G00574|Do1_05_a00001G01108|Do1_05_a00001G01355|Do1_05_a00001G01358|Do1_05_a00001G01360|Do1_05_a00001G01361|Do1_05_a00001G01362|Do1_05_a00001G01363|Do1_05_a00003G00011|Do1_05_a00003G00201|Do1_06_a00001G00217|Do1_06_a00001G01034|Do1_06_a00001G02027|Do1_06_a00001G02700|Do1_06_a00002G00004|Do1_06_a00002G00149|Do1_06_a00002G00150|Do1_06_a00002G00187|Do1_06_a00002G00731|Do1_06_a00002G00735|Do1_06_a00002G00996|Do1_06_a00002G01629|Do1_06_a00002G01635|Do1_06_a00002G01637|Do1_06_a00002G01651|Do1_06_a00002G01784|Do1_06_a00002G01992|Do1_06_a00002G02295|Do1_06_a00002G02298|Do1_06_a00002G02312|Do1_07_a00001G00124|Do1_07_a00001G00135|Do1_07_a00001G00268|Do1_07_a00001G00512|Do1_07_a00001G00529|Do1_07_a00001G00596|Do1_07_a00002G00835|Do1_07_a00002G00836|Do1_07_a00002G02049|Do1_07_a00004G00205|Do1_07_a00004G01044|Do1_07_a00004G01046|Do1_07_a00004G01058|Do1_07_a00006G00136|

defense response to other organism 12
Do1_01_a00001G01727|Do1_01_a00001G01728|Do1_01_a00001G01733|Do1_01_a00001G01734|Do1_01_a00001G01735|Do1_01_a00001G01741|Do1_01_a00001G01742|Do1_01_a00001G01743|Do1_01_a00001G01744|Do1_01_a00001G01746|Do1_01_a00001G01754|Do1_01_a00001G01755|Do1_01_a00001G01757|Do1_01_a00001G01760|Do1_01_a00001G01761|Do1_01_a00001G01765|Do1_01_a00001G01770|Do1_01_a00001G01935|Do1_01_a00001G02267|Do1_01_a00001G02274|Do1_01_a00001G02636|Do1_01_a00003G00001|Do1_02_a00001G00907|Do1_02_a00003G00790|Do1_02_a00003G01423|Do1_02_a00003G01427|Do1_02_a00003G01718|Do1_02_a00003G02060|Do1_02_a00003G02106|Do1_02_a00004G00421|Do1_02_a00004G00422|Do1_02_a00004G01171|Do1_03_a00002G01111|Do1_03_a00002G01327|Do1_03_a00002G01395|Do1_03_a00002G01466|Do1_03_a00002G01589|Do1_03_a00002G01590|Do1_03_a00002G01591|Do1_03_a00002G02055|Do1_03_a00002G02474|Do1_03_a00003G00250|Do1_04_a00001G00056|Do1_04_a00001G00217|Do1_04_a00001G00303|Do1_04_a00001G00314|Do1_04_a00001G00315|Do1_04_a00001G00316|Do1_04_a00001G00668|Do1_04_a00001G00672|Do1_04_a00001G00716|Do1_04_a00001G00992|Do1_04_a00001G02322|Do1_04_a00001G03227|Do1_04_a00004G00044|Do1_04_a00004G00088|Do1_05_a00001G00574|Do1_05_a00001G01108|Do1_05_a00001G01355|Do1_05_a00001G01358|Do1_05_a00001G01360|Do1_05_a00001G01361|Do1_05_a00001G01362|Do1_05_a00001G01363|Do1_05_a00003G00011|Do1_05_a00003G00201|Do1_06_a00001G00217|Do1_06_a00001G01034|Do1_06_a00001G02027|Do1_06_a00001G02700|Do1_06_a00002G00004|Do1_06_a00002G00149|Do1_06_a00002G00150|Do1_06_a00002G00187|Do1_06_a00002G00731|Do1_06_a00002G00735|Do1_06_a00002G00996|Do1_06_a00002G01629|Do1_06_a00002G01635|Do1_06_a00002G01637|Do1_06_a00002G01651|Do1_06_a00002G01784|Do1_06_a00002G01992|Do1_06_a00002G02295|Do1_06_a00002G02298|Do1_06_a00002G02312|Do1_07_a00001G00124|Do1_07_a00001G00135|Do1_07_a00001G00268|Do1_07_a00001G00512|Do1_07_a00001G00529|Do1_07_a00001G00596|Do1_07_a00002G00835|Do1_07_a00002G00836|Do1_07_a00002G02049|Do1_07_a00004G00205|Do1_07_a00004G01044|Do1_07_a00004G01046|Do1_07_a00004G01058|Do1_07_a00006G00136|

response to fungus 8
Do1_01_a00001G01727|Do1_01_a00001G01728|Do1_01_a00001G01733|Do1_01_a00001G01734|Do1_01_a00001G01735|Do1_01_a00001G01741|Do1_01_a00001G01742|Do1_01_a00001G01743|Do1_01_a00001G01744|Do1_01_a00001G01746|Do1_01_a00001G01754|Do1_01_a00001G01755|Do1_01_a00001G01757|Do1_01_a00001G01760|Do1_01_a00001G01761|Do1_01_a00001G01765|Do1_01_a00001G01770|Do1_01_a00001G02267|Do1_01_a00001G02274|Do1_01_a00003G00001|Do1_02_a00003G02060|Do1_02_a00003G02106|Do1_04_a00001G00056|Do1_04_a00004G00044|Do1_04_a00004G00088|Do1_05_a00001G00574|Do1_05_a00001G01355|Do1_05_a00001G01358|Do1_05_a00001G01360|Do1_05_a00001G01361|Do1_05_a00001G01362|Do1_05_a00001G01363|Do1_06_a00001G00217|Do1_06_a00002G00149|Do1_06_a00002G00150|Do1_06_a00002G01992|Do1_06_a00002G02295|Do1_06_a00002G02298|Do1_06_a00002G02312|Do1_07_a00001G00124|Do1_07_a00001G00135|Do1_07_a00001G00268|Do1_07_a00001G00512|Do1_07_a00001G00596|Do1_07_a00002G02049|Do1_07_a00004G00205|Do1_07_a00004G01058|

defense response to fungus 8
Do1_01_a00001G01727|Do1_01_a00001G01728|Do1_01_a00001G01733|Do1_01_a00001G01734|Do1_01_a00001G01735|Do1_01_a00001G01741|Do1_01_a00001G01742|Do1_01_a00001G01743|Do1_01_a00001G01744|Do1_01_a00001G01746|Do1_01_a00001G01754|Do1_01_a00001G01755|Do1_01_a00001G01757|Do1_01_a00001G01760|Do1_01_a00001G01761|Do1_01_a00001G01765|Do1_01_a00001G01770|Do1_01_a00001G02267|Do1_01_a00001G02274|Do1_01_a00003G00001|Do1_02_a00003G02060|Do1_02_a00003G02106|Do1_04_a00001G00056|Do1_04_a00004G00044|Do1_04_a00004G00088|Do1_05_a00001G00574|Do1_05_a00001G01355|Do1_05_a00001G01358|Do1_05_a00001G01360|Do1_05_a00001G01361|Do1_05_a00001G01362|Do1_05_a00001G01363|Do1_06_a00001G00217|Do1_06_a00002G00149|Do1_06_a00002G00150|Do1_06_a00002G01992|Do1_06_a00002G02295|Do1_06_a00002G02298|Do1_06_a00002G02312|Do1_07_a00001G00124|Do1_07_a00001G00135|Do1_07_a00001G00268|Do1_07_a00001G00512|Do1_07_a00001G00596|Do1_07_a00002G02049|Do1_07_a00004G00205|Do1_07_a00004G01058|

polysaccharide binding 8
Do1_00288G00001|Do1_01_a00001G02379|Do1_01_a00003G00528|Do1_01_a00004G01626|Do1_01_a00004G01634|Do1_01_a00004G01635|Do1_01_a00004G01640|Do1_01_a00004G01663|Do1_01_a00004G01664|Do1_01_a00004G01666|Do1_02_a00001G01061|Do1_02_a00001G01062|Do1_02_a00004G01357|Do1_02_a00004G01366|Do1_03_a00001G00512|Do1_03_a00001G00568|Do1_03_a00001G01525|Do1_03_a00001G01536|Do1_03_a00001G01539|Do1_03_a00001G01547|Do1_03_a00001G01549|Do1_03_a00001G01554|Do1_03_a00002G01509|Do1_04_a00001G01858|Do1_04_a00001G03447|Do1_04_a00001G03448|Do1_04_a00005G00107|Do1_04_a00005G00108|Do1_04_a00005G00110|Do1_04_a00005G00111|Do1_04_a00005G00130|Do1_05_a00001G00301|Do1_05_a00001G00306|Do1_05_a00001G01433|Do1_05_a00001G01441|Do1_05_a00001G02117|Do1_05_a00001G02371|Do1_05_a00001G02372|Do1_05_a00001G02393|Do1_05_a00003G00031|Do1_05_a00003G00082|Do1_05_a00003G00276|Do1_06_a00001G01824|Do1_06_a00002G00471|Do1_06_a00002G00904|Do1_06_a00002G01038|Do1_07_a00001G00051|Do1_07_a00002G00060|Do1_07_a00002G00061|Do1_07_a00002G00064|Do1_07_a00002G01316|Do1_07_a00002G01465|Do1_07_a00004G00063|Do1_07_a00004G00068|Do1_07_a00004G00075|Do1_07_a00004G00417|Do1_07_a00004G00764|Do1_07_a00004G00780|Do1_07_a00004G00782|Do1_07_a00004G00784|Do1_07_a00004G00796|Do1_07_a00004G00812|Do1_07_a00004G01082|Do1_a00045G00108|

#--------------------
# Sweden_W_C

polysaccharide binding 8
Do1_00288G00001|Do1_01_a00001G02379|Do1_01_a00003G00528|Do1_01_a00004G01626|Do1_01_a00004G01635|Do1_01_a00004G01663|Do1_01_a00004G01664|Do1_01_a00004G01666|Do1_02_a00001G01061|Do1_02_a00001G01062|Do1_02_a00004G01357|Do1_02_a00004G01366|Do1_03_a00001G00512|Do1_03_a00001G00568|Do1_03_a00001G01525|Do1_03_a00001G01527|Do1_03_a00001G01536|Do1_03_a00001G01539|Do1_03_a00001G01545|Do1_03_a00001G01547|Do1_03_a00001G01549|Do1_03_a00001G01554|Do1_03_a00002G01509|Do1_04_a00001G01858|Do1_04_a00001G03447|Do1_04_a00001G03448|Do1_04_a00005G00107|Do1_04_a00005G00108|Do1_04_a00005G00110|Do1_04_a00005G00111|Do1_04_a00005G00130|Do1_05_a00001G00301|Do1_05_a00001G00306|Do1_05_a00001G01433|Do1_05_a00001G01441|Do1_05_a00001G02117|Do1_05_a00001G02371|Do1_05_a00001G02372|Do1_05_a00001G02393|Do1_05_a00003G00031|Do1_05_a00003G00082|Do1_05_a00003G00276|Do1_06_a00001G01824|Do1_06_a00002G00471|Do1_06_a00002G00904|Do1_06_a00002G01038|Do1_07_a00001G00051|Do1_07_a00002G00060|Do1_07_a00002G00061|Do1_07_a00002G00064|Do1_07_a00002G01316|Do1_07_a00002G01465|Do1_07_a00004G00063|Do1_07_a00004G00068|Do1_07_a00004G00075|Do1_07_a00004G00417|Do1_a00045G00108|

response to fungus 5


defense response to fungus 5

response to biotic stimulus 8
biological process involved in interspecies interaction between organisms
response to external biotic stimulus
response to other organism
defense response to other organism

Do1_01_a00001G01727|Do1_01_a00001G01728|Do1_01_a00001G01733|Do1_01_a00001G01734|Do1_01_a00001G01735|Do1_01_a00001G01741|Do1_01_a00001G01742|Do1_01_a00001G01743|Do1_01_a00001G01744|Do1_01_a00001G01746|Do1_01_a00001G01754|Do1_01_a00001G01755|Do1_01_a00001G01757|Do1_01_a00001G01760|Do1_01_a00001G01761|Do1_01_a00001G01765|Do1_01_a00001G01770|Do1_01_a00001G01935|Do1_01_a00001G02267|Do1_01_a00001G02274|Do1_01_a00001G02636|Do1_01_a00003G00001|Do1_02_a00001G00907|Do1_02_a00003G00790|Do1_02_a00003G01423|Do1_02_a00003G01427|Do1_02_a00003G01718|Do1_02_a00003G02060|Do1_02_a00003G02106|Do1_02_a00004G00421|Do1_02_a00004G00422|Do1_02_a00004G01171|Do1_03_a00002G01111|Do1_03_a00002G01327|Do1_03_a00002G01395|Do1_03_a00002G01466|Do1_03_a00002G01589|Do1_03_a00002G01590|Do1_03_a00002G01591|Do1_03_a00002G02055|Do1_03_a00002G02474|Do1_03_a00003G00250|Do1_04_a00001G00056|Do1_04_a00001G00217|Do1_04_a00001G00303|Do1_04_a00001G00314|Do1_04_a00001G00315|Do1_04_a00001G00316|Do1_04_a00001G00668|Do1_04_a00001G00672|Do1_04_a00001G00716|Do1_04_a00001G00992|Do1_04_a00001G02322|Do1_04_a00001G03227|Do1_04_a00004G00044|Do1_04_a00004G00088|Do1_05_a00001G00510|Do1_05_a00001G00574|Do1_05_a00001G01108|Do1_05_a00001G01355|Do1_05_a00001G01358|Do1_05_a00001G01360|Do1_05_a00001G01361|Do1_05_a00001G01362|Do1_05_a00001G01363|Do1_05_a00003G00011|Do1_05_a00003G00201|Do1_06_a00001G00217|Do1_06_a00001G01034|Do1_06_a00001G02027|Do1_06_a00001G02700|Do1_06_a00002G00004|Do1_06_a00002G00149|Do1_06_a00002G00150|Do1_06_a00002G00187|Do1_06_a00002G00731|Do1_06_a00002G00735|Do1_06_a00002G00996|Do1_06_a00002G01629|Do1_06_a00002G01635|Do1_06_a00002G01637|Do1_06_a00002G01651|Do1_06_a00002G01784|Do1_06_a00002G01992|Do1_06_a00002G02295|Do1_06_a00002G02312|Do1_07_a00001G00124|Do1_07_a00001G00135|Do1_07_a00001G00268|Do1_07_a00001G00512|Do1_07_a00001G00529|Do1_07_a00001G00596|Do1_07_a00002G00835|Do1_07_a00002G00836|Do1_07_a00002G02049|Do1_07_a00004G00205|Do1_07_a00004G01044|Do1_07_a00004G01046|Do1_07_a00004G01058|Do1_07_a00006G00136|


#-----------------------
# genes
DoctH0-5_RagTag	Liftoff	gene	11999864	12000933	.	-	.	ID=Do1_01_a00001G01728
DoctH0-5_RagTag	Liftoff	gene	11968082	11970716	.	-	.	ID=Do1_01_a00001G01733


DoctH0-8_RagTag	Liftoff	gene	3915855	3916178	.	+	.	ID=Do1_02_a00003G00790

Do1_02_a00003G00810










