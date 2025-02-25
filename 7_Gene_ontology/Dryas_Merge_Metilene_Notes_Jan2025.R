##############################
# Merge DMR genes and RNAseq files
# October 2023
##############################
#mkdir /lustre04/scratch/celphin/Dryas/methylkit_merged_data/

# Get Methylkit overdisp DMRs
cd /lustre04/scratch/celphin/Dryas/methylkit_merged_data/
cp /lustre04/scratch/celphin/Dryas_large_folders/intersections/*.txt /lustre04/scratch/celphin/Dryas/methylkit_merged_data/
cp /lustre04/scratch/celphin/Dryas_large_folders/intersections/txt/* /lustre04/scratch/celphin/Dryas/methylkit_merged_data/

# Get Methylkit overdisp DMR genes
cd /lustre04/scratch/celphin/Dryas/methylkit_merged_data/
cp  /lustre04/scratch/celphin/Dryas/snpEff/methylkit/*.uniq_genes /lustre04/scratch/celphin/Dryas/methylkit_merged_data/

# get metilene DMRs
cp /lustre04/scratch/celphin/Dryas/DMR_bedGraph_files/*.bedGraph /lustre04/scratch/celphin/Dryas/methylkit_merged_data/


# annotation info
cd /lustre04/scratch/celphin/Dryas/MS_Dryas_Merged_Data/original_data
cp Dryas_octopetala_H1.gff3 /lustre04/scratch/celphin/Dryas/methylkit_merged_data/

# Interproscan
cd /lustre04/scratch/celphin/Dryas/methylkit_merged_data/
cp /lustre04/scratch/celphin/Dryas/GO_enrichment/interproscan_dryas_full3.tsv .
cp /home/celphin/scratch/Dryas/MS_Dryas_Merged_Data/original_data/interproscan_dryas_full.tsv .


##################################################
# Process Interproscan file

# Interproscan output
#Do1_04_a00001G01549V1.1 a660bcfabce7c9b57fe024301da0870e        560     SUPERFAMILY     SSF48452        TPR-like       138     295     1.23E-12        T       29-09-2023      IPR011990       Tetratricopeptide-like helicaldomain superfamily      GO:0005515

####################################
# Edit Interproscan file to remove duplicates
# format 
awk -v FS="\t" '{print $1 "\t" $4 "\t" $6 "\t" $12 "\t" $13 "\t" $14}' interproscan_dryas_full0.tsv | wc -l 
# 216 861

awk -v FS="\t" '{print $1 "\t" $4 "\t" $6 "\t" $12 "\t" $13 "\t" $14}' interproscan_dryas_full0.tsv | sort | uniq |wc -l
# 157 061

awk -v FS="\t" '{print $1 "\t" $4 "\t" $6 "\t" $12 "\t" $13 "\t" $14}' interproscan_dryas_full0.tsv | sort | uniq > interproscan_dryas_full.tsv

grep -v $'\t''-'$'\t''-'$'\t' interproscan_dryas_full.tsv > interproscan_dryas_full1.tsv

sed 's/|/,/g' interproscan_dryas_full1.tsv | sort -u > interproscan_dryas_full2.tsv

wc -l interproscan_dryas_full2.tsv
# 95 869

##############################################
tmux new-session -s GO
tmux attach-session -t GO

cd /lustre04/scratch/celphin/Dryas/methylkit_merged_data/

module load StdEnv/2023
module load r/4.4.0

R

# install.packages("dplyr", "tidyr")

library(dplyr)
library(tidyr)
library(stringr)

#-------------------------
# load GO ont data

path="/lustre04/scratch/celphin/Dryas/methylkit_merged_data"
Gene_ont_file <- "interproscan_dryas_full2.tsv"
gene_ont <- read.delim(paste0(path,"/", Gene_ont_file), header = FALSE, sep = "\t", na.strings = "-", colClasses = c("character", "character", "character", "character"))

colnames(gene_ont) <- c( "gene", "Pfam", "descrip1", "INTPRO", "descrip2", "GOterm")
length(unique(gene_ont$INTPRO))
# 8644

nrow(gene_ont)
#  95869 # getting read in properly, half of the rows missing before

#-----------------------------
# formatting Interproscan to have no duplicates of genes - one row per gene

# collapse GO terms
collapsed_go_terms1 <- gene_ont %>%
  group_by(gene) %>%
  summarize(
    descrip1 = paste(unique(descrip1), collapse = ","),  # Collapse unique descriptions
    descrip2 = paste(unique(descrip2), collapse = ","),  # Collapse unique descriptions
    GOterm = paste(unique(GOterm), collapse = ","),    # Collapse unique GO terms
    INTPRO = paste(unique(INTPRO), collapse = ","), .groups="keep"      # Collapse unique IPR terms
  )

# remove duplicate values in a list
cleaned_tibble <- collapsed_go_terms1 %>%
  separate_rows(GOterm, sep = ",") %>%  # Split the GOterm string into multiple rows
  separate_rows(INTPRO, sep = ",") %>%  # Split the INTPRO string into multiple rows
  separate_rows(descrip1, sep = ",") %>%  # Split the descrip string into multiple rows
  separate_rows(descrip2, sep = ",") %>%  # Split the descrip string into multiple rows
  filter(GOterm != "NA") %>%              # Remove rows with '-'
  distinct(gene, GOterm, INTPRO,descrip1,descrip2,.keep_all = TRUE ) #%>% # Keep unique terms with gene info

# recollapse 
collapsed_go_terms2 <- cleaned_tibble %>%
  group_by(gene) %>%
  summarize(
    INTPRO = paste(sort(unique(INTPRO)), collapse = ","),           # Collapse IPR terms with unique values
    descrip1 = paste(sort(unique(descrip1)), collapse = ","),        # Collapse descriptions with unique values
    descrip2 = paste(sort(unique(descrip2)), collapse = ","),        # Collapse descriptions with unique values
    GOterm = paste(sort(unique(GOterm)), collapse = ",") , .groups="keep"         # Collapse GO terms with unique values
  )

# format for ermineJ
collapsed_go_terms <- collapsed_go_terms2 %>%
  #mutate(gene2 = gene) %>%
  select(gene, everything())

collapsed_go_terms_df <- as.data.frame(collapsed_go_terms)

gene_ont <- collapsed_go_terms_df

nrow(gene_ont)
# 14963

# write out file
utils::write.table(x=gene_ont , file=paste0(path,"/interproscan_dryas_full3.tsv"), append = FALSE, sep = "\t", dec = ".", row.names = FALSE, col.names = TRUE)


#####################################################################
#Process gene annotation gff3 file:

gff3_file <- read.table("Dryas_octopetala_H1.gff3", header=FALSE)
gff3_scaffold_start_end <- gff3_file[, c(1, 3, 4, 5, 9)]
dryas_genes <- gff3_scaffold_start_end[gff3_scaffold_start_end[, 2] == "gene", ]
dryas_genes <- dryas_genes[, c(1, 3, 4, 5)]
gene_col_names <- c("Scaffold", "Gene_Start", "Gene_End", "Gene")
colnames(dryas_genes) <- gene_col_names
#remove ID prefix:
dryas_genes$Gene <- sub("^ID=", "", dryas_genes$Gene)

print(head(dryas_genes))
#Should have start, end, gene

#################################
# Read in edited functional annotation
interproscan <- read.table("/lustre04/scratch/celphin/Dryas/methylkit_merged_data/interproscan_dryas_full3.tsv", header=TRUE)
print(head(interproscan))

interproscan$Gene <- sub("V1.1", "", interproscan$gene)

#---------------
# join to the total annotation
dryas_genes_interpro <- left_join(dryas_genes, interproscan, by="Gene")

####################################
# Get file lists

file_path <- "/lustre04/scratch/celphin/Dryas/methylkit_merged_data"
file_list <- list.files(file_path, full.names=TRUE)

methylkitDMRs_list <- file_list[grep("Methylkit", file_list)]
RNA_list <- file_list[grep("RNA", file_list)]

methylkitDMRs_info_list <- methylkitDMRs_list[grep(".txt", methylkitDMRs_list)]
methylkitDMRs_genenames_list <- methylkitDMRs_list[grep(".uniq_genes", methylkitDMRs_list)]

RNA_info_list <- RNA_list[grep("_DERs.txt", RNA_list)]
RNA_updown_list <- RNA_list[grep("updown", RNA_list)]


#####################################################################
#Read in the methylkit DMRs, add filename column, join for all data
# *ethylkit*.txt

# Read the data into a list of data frames, treating the first row as headers
list_of_data <- lapply(methylkitDMRs_info_list, function(file_path) {
  read.table(file_path, header = TRUE)  # Read each file with header
})

# Extract basenames from the file paths
file_basenames <- basename(methylkitDMRs_info_list)

# Add the file basename as a new column to each data frame in the list
list_of_data_with_basename <- lapply(seq_along(list_of_data), function(i) {
  data <- list_of_data[[i]]
  data$file_basename <- file_basenames[i]  # Add the basename as a new column
  return(data)
})

# Combine all data frames into one using rbind
methylkitDMRs_info_combined_data <- do.call(rbind, list_of_data_with_basename)

head(methylkitDMRs_info_combined_data)

colnames(methylkitDMRs_info_combined_data)

# Extract DMR site/origin info
methylkitDMRs_info_combined_data$origin <- methylkitDMRs_info_combined_data$file_basename
methylkitDMRs_info_combined_data$origin <- sub("Rand_", "", methylkitDMRs_info_combined_data$origin)
methylkitDMRs_info_combined_data$origin <- sub("Methylkit_", "", methylkitDMRs_info_combined_data$origin)
methylkitDMRs_info_combined_data$origin <- sub("_DMRs.txt", "", methylkitDMRs_info_combined_data$origin)

methylkitDMRs_info_combined_data$site <- methylkitDMRs_info_combined_data$origin
methylkitDMRs_info_combined_data$context <- methylkitDMRs_info_combined_data$origin
methylkitDMRs_info_combined_data$random <- methylkitDMRs_info_combined_data$origin
methylkitDMRs_info_combined_data$perdiff <- methylkitDMRs_info_combined_data$origin

# find sites
methylkitDMRs_info_combined_data$site[grep("SE_W_C", methylkitDMRs_info_combined_data$origin)]<- "SE_W_C"
methylkitDMRs_info_combined_data$site[grep("Pheno", methylkitDMRs_info_combined_data$origin)]<- "Pheno"
methylkitDMRs_info_combined_data$site[grep("HL", methylkitDMRs_info_combined_data$origin)]<- "HL"
methylkitDMRs_info_combined_data$site[grep("SE_HL", methylkitDMRs_info_combined_data$origin)]<- "SE_HL"

methylkitDMRs_info_combined_data$site[grep("MEAD", methylkitDMRs_info_combined_data$origin)]<- "MEAD_W_C"
methylkitDMRs_info_combined_data$site[grep("WILL", methylkitDMRs_info_combined_data$origin)]<- "WILL_W_C"
methylkitDMRs_info_combined_data$site[grep("DRY", methylkitDMRs_info_combined_data$origin)]<- "DRY_W_C"
methylkitDMRs_info_combined_data$site[grep("CASS", methylkitDMRs_info_combined_data$origin)]<- "CASS_W_C"
methylkitDMRs_info_combined_data$site[grep("FERT", methylkitDMRs_info_combined_data$origin)]<- "FERT_W_C"

methylkitDMRs_info_combined_data$site[grep("LAT", methylkitDMRs_info_combined_data$origin)]<- "LAT_W_C"
methylkitDMRs_info_combined_data$site[grep("ALAS", methylkitDMRs_info_combined_data$origin)]<- "ALAS_W_C"
methylkitDMRs_info_combined_data$site[grep("SVAL", methylkitDMRs_info_combined_data$origin)]<- "SVAL_W_C"

# find contexts
methylkitDMRs_info_combined_data$context[grep("CHH", methylkitDMRs_info_combined_data$origin)]<- "CHH"
methylkitDMRs_info_combined_data$context[-grep("CHH", methylkitDMRs_info_combined_data$origin)]<- "CpG"

# random
methylkitDMRs_info_combined_data$random[grep("rand", methylkitDMRs_info_combined_data$origin)]<- "rand"
methylkitDMRs_info_combined_data$random[-grep("rand", methylkitDMRs_info_combined_data$origin)]<- "non-rand"
methylkitDMRs_info_combined_data$random[grep("Rand", methylkitDMRs_info_combined_data$origin)]<- "rand"

# percent DMR diff
methylkitDMRs_info_combined_data$perdiff[grep("10", methylkitDMRs_info_combined_data$origin)]<- "10"
methylkitDMRs_info_combined_data$perdiff[grep("25", methylkitDMRs_info_combined_data$origin)]<- "25"


unique(as.factor(methylkitDMRs_info_combined_data$site))
unique(as.factor(methylkitDMRs_info_combined_data$context))
unique(as.factor(methylkitDMRs_info_combined_data$random))
unique(as.factor(methylkitDMRs_info_combined_data$perdiff))

methylkitDMRs_info_combined_data0 <- methylkitDMRs_info_combined_data
methylkitDMRs_info_combined_data = subset(methylkitDMRs_info_combined_data0, select = -c(file_basename, origin) )

#####################################################################
#Read in the methylkit DMR nearest gene names, add filename column, join for all data
#*.uniq_genes

# Read the data into a list of data frames, treating the first row as headers
list_of_data <- lapply(methylkitDMRs_genenames_list, function(file_path) {
  read.table(file_path, header = FALSE)  # Read each file with header
})

# Extract basenames from the file paths
file_basenames <- basename(methylkitDMRs_genenames_list)

# Add the file basename as a new column to each data frame in the list
list_of_data_with_basename <- lapply(seq_along(list_of_data), function(i) {
  data <- list_of_data[[i]]
  data$file_basename <- file_basenames[i]  # Add the basename as a new column
  return(data)
})

# Combine all data frames into one using rbind
methylkitDMRs_genenames_combined_data <- do.call(rbind, list_of_data_with_basename)
colnames(methylkitDMRs_genenames_combined_data)<- c("chr", "start", "end", "Gene", "file_basename")
methylkitDMRs_genenames_combined_data$gene <- methylkitDMRs_genenames_combined_data$Gene
methylkitDMRs_genenames_combined_data$Gene <- sub("V1.1", "", methylkitDMRs_genenames_combined_data$gene)

head(methylkitDMRs_genenames_combined_data)

#----------------------
# Extract DMR site/origin info
methylkitDMRs_genenames_combined_data$origin <- methylkitDMRs_genenames_combined_data$file_basename
methylkitDMRs_genenames_combined_data$origin <- sub("Rand_", "", methylkitDMRs_genenames_combined_data$origin)
methylkitDMRs_genenames_combined_data$origin <- sub("Methylkit_", "", methylkitDMRs_genenames_combined_data$origin)
methylkitDMRs_genenames_combined_data$origin <- sub("_DMRs.txt", "", methylkitDMRs_genenames_combined_data$origin)


methylkitDMRs_genenames_combined_data$site <- methylkitDMRs_genenames_combined_data$origin
methylkitDMRs_genenames_combined_data$context <- methylkitDMRs_genenames_combined_data$origin
methylkitDMRs_genenames_combined_data$random <- methylkitDMRs_genenames_combined_data$origin
methylkitDMRs_genenames_combined_data$perdiff <- methylkitDMRs_genenames_combined_data$origin

# find sites
methylkitDMRs_genenames_combined_data$site[grep("SE_W_C", methylkitDMRs_genenames_combined_data$origin)]<- "SE_W_C"
methylkitDMRs_genenames_combined_data$site[grep("Pheno", methylkitDMRs_genenames_combined_data$origin)]<- "Pheno"
methylkitDMRs_genenames_combined_data$site[grep("HL", methylkitDMRs_genenames_combined_data$origin)]<- "HL"
methylkitDMRs_genenames_combined_data$site[grep("SE_HL", methylkitDMRs_genenames_combined_data$origin)]<- "SE_HL"

methylkitDMRs_genenames_combined_data$site[grep("MEAD", methylkitDMRs_genenames_combined_data$origin)]<- "MEAD_W_C"
methylkitDMRs_genenames_combined_data$site[grep("WILL", methylkitDMRs_genenames_combined_data$origin)]<- "WILL_W_C"
methylkitDMRs_genenames_combined_data$site[grep("DRY", methylkitDMRs_genenames_combined_data$origin)]<- "DRY_W_C"
methylkitDMRs_genenames_combined_data$site[grep("CASS", methylkitDMRs_genenames_combined_data$origin)]<- "CASS_W_C"
methylkitDMRs_genenames_combined_data$site[grep("FERT", methylkitDMRs_genenames_combined_data$origin)]<- "FERT_W_C"

methylkitDMRs_genenames_combined_data$site[grep("LAT", methylkitDMRs_genenames_combined_data$origin)]<- "LAT_W_C"
methylkitDMRs_genenames_combined_data$site[grep("ALAS", methylkitDMRs_genenames_combined_data$origin)]<- "ALAS_W_C"
methylkitDMRs_genenames_combined_data$site[grep("SVAL", methylkitDMRs_genenames_combined_data$origin)]<- "SVAL_W_C"

# find contexts
methylkitDMRs_genenames_combined_data$context[grep("CHH", methylkitDMRs_genenames_combined_data$origin)]<- "CHH"
methylkitDMRs_genenames_combined_data$context[-grep("CHH", methylkitDMRs_genenames_combined_data$origin)]<- "CpG"

# random
methylkitDMRs_genenames_combined_data$random[grep("rand", methylkitDMRs_genenames_combined_data$origin)]<- "rand"
methylkitDMRs_genenames_combined_data$random[-grep("rand", methylkitDMRs_genenames_combined_data$origin)]<- "non-rand"
methylkitDMRs_genenames_combined_data$random[grep("Rand", methylkitDMRs_genenames_combined_data$origin)]<- "rand"

# percent DMR diff
methylkitDMRs_genenames_combined_data$perdiff[grep("10", methylkitDMRs_genenames_combined_data$origin)]<- "10"
methylkitDMRs_genenames_combined_data$perdiff[grep("25", methylkitDMRs_genenames_combined_data$origin)]<- "25"


unique(as.factor(methylkitDMRs_genenames_combined_data$site))
unique(as.factor(methylkitDMRs_genenames_combined_data$context))
unique(as.factor(methylkitDMRs_genenames_combined_data$random))
unique(as.factor(methylkitDMRs_genenames_combined_data$perdiff))

methylkitDMRs_genenames_combined_data0 <- methylkitDMRs_genenames_combined_data
methylkitDMRs_genenames_combined_data = subset(methylkitDMRs_genenames_combined_data0, select = -c(file_basename, origin) )

#####################################################################
#Read in the DEG gene names, add filename column, join for all data
# RNA*_DERs.txt

# Read the data into a list of data frames, treating the first row as headers
list_of_data <- lapply(RNA_info_list, function(file_path) {
  read.table(file_path, header = TRUE)  # Read each file with header
})

# Extract basenames from the file paths
file_basenames <- basename(RNA_info_list)

# Add the file basename as a new column to each data frame in the list
list_of_data_with_basename <- lapply(seq_along(list_of_data), function(i) {
  data <- list_of_data[[i]]
  data$file_basename <- file_basenames[i]  # Add the basename as a new column
  return(data)
})

# Combine all data frames into one using rbind
RNA_info_combined_data <- do.call(rbind, list_of_data_with_basename)

head(RNA_info_combined_data)

colnames(RNA_info_combined_data)
RNA_info_combined_data$gene <- rownames(RNA_info_combined_data)
RNA_info_combined_data$Gene <- sub("V1.1", "", RNA_info_combined_data$gene)


# Extract DEG site/origin info and remove filename
RNA_info_combined_data$origin <- RNA_info_combined_data$file_basename
RNA_info_combined_data$origin <- sub("RNA_", "", RNA_info_combined_data$origin)
RNA_info_combined_data$origin <- sub("_DERs.txt", "", RNA_info_combined_data$origin)

RNA_info_combined_data$RNAsite <- RNA_info_combined_data$origin

# find sites
RNA_info_combined_data$RNAsite[grep("Seedling_W_C", RNA_info_combined_data$origin)]<- "SE_W_C"
RNA_info_combined_data$RNAsite[grep("Alex_W_C", RNA_info_combined_data$origin)]<- "ALEX_W_C"
RNA_info_combined_data$RNAsite[grep("Sweden_W_C", RNA_info_combined_data$origin)]<- "LAT_W_C"
RNA_info_combined_data$RNAsite[grep("Alaska_W_C", RNA_info_combined_data$origin)]<- "ALAS_W_C"
RNA_info_combined_data$RNAsite[grep("Norway_W_C", RNA_info_combined_data$origin)]<- "NORW_W_C"

unique(as.factor(RNA_info_combined_data$RNAsite))

RNA_info_combined_data0 <- RNA_info_combined_data
RNA_info_combined_data = subset(RNA_info_combined_data0, select = -c(file_basename, origin) )


#####################################################################
#Read in the DEG up or down regulation info, add filename column, join for all data
# RNA*_updown.txt

# Read the data into a list of data frames, treating the first row as headers
list_of_data <- lapply(RNA_updown_list, function(file_path) {
  read.table(file_path, header = TRUE)  # Read each file with header
})

# Extract basenames from the file paths
file_basenames <- basename(RNA_updown_list)

# Add the file basename as a new column to each data frame in the list
list_of_data_with_basename <- lapply(seq_along(list_of_data), function(i) {
  data <- list_of_data[[i]]
  data$file_basename <- file_basenames[i]  # Add the basename as a new column
  return(data)
})

# Combine all data frames into one using rbind
RNA_updown_combined_data <- do.call(rbind, list_of_data_with_basename)

RNA_updown_combined_data$gene <- RNA_updown_combined_data$Gene
RNA_updown_combined_data$Gene <- sub("V1.1", "", RNA_updown_combined_data$gene)

head(RNA_updown_combined_data)

colnames(RNA_updown_combined_data)

# Extract DEG site/origin info and remove filename
RNA_updown_combined_data$origin <- RNA_updown_combined_data$file_basename
RNA_updown_combined_data$origin <- sub("RNA_", "", RNA_updown_combined_data$origin)
RNA_updown_combined_data$origin <- sub("_DERs_updown.txt", "", RNA_updown_combined_data$origin)

RNA_updown_combined_data$RNAsite <- RNA_updown_combined_data$origin

# find sites
RNA_updown_combined_data$RNAsite[grep("Seedling_W_C", RNA_updown_combined_data$origin)]<- "SE_W_C"
RNA_updown_combined_data$RNAsite[grep("Alex_W_C", RNA_updown_combined_data$origin)]<- "ALEX_W_C"
RNA_updown_combined_data$RNAsite[grep("Sweden_W_C", RNA_updown_combined_data$origin)]<- "LAT_W_C"
RNA_updown_combined_data$RNAsite[grep("Alaska_W_C", RNA_updown_combined_data$origin)]<- "ALAS_W_C"
RNA_updown_combined_data$RNAsite[grep("Norway_W_C", RNA_updown_combined_data$origin)]<- "NORW_W_C"

unique(as.factor(RNA_updown_combined_data$RNAsite))

RNA_updown_combined_data0 <- RNA_updown_combined_data
RNA_updown_combined_data = subset(RNA_updown_combined_data0, select = -c(file_basename, origin) )



#################################################################
# join all by the gene names

DMR_total_data <- left_join(methylkitDMRs_genenames_combined_data, methylkitDMRs_info_combined_data, by = join_by(chr, start, end, site, context, random, perdiff))
head(DMR_total_data)

#--------------------
# join DEGs with up and down regulation by gene name
RNA_total_data <- left_join(RNA_info_combined_data, RNA_updown_combined_data, by=join_by(Gene, gene, RNAsite))
head(RNA_total_data)

#---------------------
# join annotation and interproscan data with DMRdata and DEG data

genes_RNA_merged_data <- left_join(dryas_genes_interpro, RNA_total_data, by = join_by(Gene, gene))
genes_RNA_DMR_merged_data <- left_join(genes_RNA_merged_data, DMR_total_data, by = join_by(Gene, gene))

head(genes_RNA_DMR_merged_data)

write.table(genes_RNA_DMR_merged_data, "genes_RNA_MethylkitDMR_merged_data.tsv", sep = "\t", quote = FALSE, row.names = FALSE)














########################
#OLDER notes
######################



#------------------------
# make wide by file type??

# Remove duplicate rows based on all columns
data_no_duplicates <- data %>% distinct()

# Reshape the data to wide format
# Assuming you want to spread by the "Gene_ID" and have the "DMR_Type" as columns

wide_data <- data_no_duplicates %>%
  pivot_wider(names_from = DMR_Type, values_from = Score, values_fill = NA)

# View the result
head(wide_data)




#--------------------
# extract site etc. info from filenames

# DMRs - SITE, Conext, %diff

# DEGs - SITE



#-----------------
# remove filename columns 


# try removing duplicates 


#-----
sd_chr_list1 <- stringr::str_split(file_basenames, "_", simplify =TRUE)
sd_chr_list2 <- sd_chr_list1[,3]
names(lsd) <- sd_chr_list2
RepeatAbundance_total <- dplyr::bind_rows(lsd, .id = 'chromosome')
colnames(RepeatAbundance_total) <- c("Chromosome", "Genome_position", "RepeatAbundance")


####################################
    #Turn into: Scaffold, Ext_Start, Ext_End, Gene, Origin
read_blast <- function(file_name) {
    title_string <- gsub("cleaned_blast_ref_(.*?)\\.out", "\\1", file_name)
    blast_table <- read.table(file_name, header=FALSE)
    subset_blast_table <- blast_table[, c(1, 2)]
    subset_blast_table <- subset_blast_table %>% mutate(Origin = title_string)
    colnames(subset_blast_table) <- c("DMR", "Gene", "Origin")
    subset_blast_table <- separate(subset_blast_table, DMR, into = c("Scaffold", "Ext_Start", "Ext_End"), sep = "[:,\\-]")
   # print(head(subset_blast_table, n=10))
    return(data.frame(subset_blast_table))
}

#Read blasts:
Phenology_blast <- read_blast("original_data/cleaned_blast_ref_Mat_Sen.out")
# needs compute module

#Seedlings:
Seedling_W_C_blast <- read_blast("original_data/cleaned_blast_ref_SE_W_C.out")
Seedling_L_H_blast <- read_blast("original_data/cleaned_blast_ref_SE_L_H.out")

#Site Specific
Alaska_W_C_blast <- read_blast("original_data/cleaned_blast_ref_Alaska_W_C.out") 
Svalbard_W_C_blast <- read_blast("original_data/cleaned_blast_ref_Svalbard_W_C.out")
Nunavut_W_C_blast <- read_blast("original_data/cleaned_blast_ref_Nunavut_W_C.out")
Sweden_W_C_blast <- read_blast("original_data/cleaned_blast_ref_Sweden_W_C.out")

#Wild
Wild_W_C_blast <- read_blast("original_data/cleaned_blast_ref_Wild_W_C.out")
Parent_W_C_blast <- read_blast("original_data/cleaned_blast_ref_Parent_W_C.out")

#Wild_DO_DI_blast <- read.blast("original_data/cleaned_blast_ref_Wild_Species_DO_DI.out") #Missing blast output
Wild_L_H_blast <-read_blast("original_data/cleaned_blast_ref_Wild_Lat_L_H.out")
print("All blast files loaded")


###########################################################################
#Merge methylkit:
#Merge with methylkit data:
    #Read methylkit bedGraph (post filter)
    #Merge by scaffold and location
###In progress
merge_methylkit <- function(blast_table, methylkit_file_name) {
     methylkit_table <- read.table(methylkit_file_name, header=FALSE) 
    colnames(methylkit_table) <- methylkit_colnames
    print(head(methylkit_table))    
    merged <- merge(blast_table, methylkit_table, by = "Scaffold")
    merged <- merged %>% filter(DMR_Start >= Ext_Start, DMR_End <= Ext_End)
    merged <- merged[c("Origin", "Gene", "Scaffold", "DMR_Start", "DMR_End", "q-value", "meanmethyl", "CpG", "meanW", "meanC")]
    merged$Gene <- sub("V1.1$", "", merged$Gene)
    return(merged)
}

#Attach methylkit info:
Phenology <- merge_methylkit(Phenology, "original_data/Mat_Sen_70_5_4_0.9_qval.0.001.out")

#site specific places
Alaska_W_C <- merge_methylkit(Alaska_W_C_blast, "original_data/Alaska_W_C_70_5_4_0.9_qval.0.001.out") 
Svalbard_W_C <- merge_methylkit(Svalbard_W_C_blast, "original_data/Svalbard_W_C_70_5_4_0.9_qval.0.001.out")
Nunavut_W_C <- merge_methylkit(Nunavut_W_C_blast, "original_data/Nunavut_W_C_70_5_4_0.9_qval.0.001.out")
Sweden_W_C <- merge_methylkit(Sweden_W_C_blast, "original_data/Sweden_W_C_70_5_4_0.9_qval.0.001.out")
#Wild
Wild_W_C <- merge_methylkit(Wild_W_C_blast, "original_data/Wild_W_C_70_5_4_0.9_qval.0.001.out")
Wild_L_H <- merge_methylkit(Wild_L_H_blast, "original_data/Wild_Lat_L_H_70_5_4_0.9_qval.0.001.out")
Parent_W_C <- merge_methylkit(Parent_W_C_blast, "original_data/P_W_C_70_5_4_0.9_qval.0.001.out")
#Seedlings:
SE_W_C <- merge_methylkit(Seedling_W_C_blast, "original_data/SE_W_C_70_5_4_0.9_qval.0.001.out")
SE_L_H <- merge_methylkit(Seedling_L_H_blast, "original_data/SE_L_H_70_5_4_0.9_qval.0.001.out")
print("Metilene loaded")


############################################################################



############################################################################
#process_go terms

process_go_terms <- function(goterm_file) {
    goterm_file_lines <- readLines(goterm_file)
    pattern <- c("GO:")
    filtered_goterm_lines <- goterm_file_lines[grep(pattern, goterm_file_lines)]
    filtered_goterms <- data.frame(filtered_goterm_lines)
    colnames(filtered_goterms) <- c("Lines")
    split_file <- data.frame(separate(filtered_goterms, Lines, into = c("IPRscan", "GO_Name", "GO_Terms"), sep = ">|;"))
    goterms <- split_file[,c("GO_Name", "GO_Terms")]
    goterms$GO_Terms <- gsub("\\s", "", goterms$GO_Terms)
    no_dub_goterms <- distinct(goterms, GO_Terms, GO_Name)
    return(data.frame(no_dub_goterms))
    #print(head(split_file))
}
go_terms <- process_go_terms(gotermsfile)
print(head(go_terms))

############################################################################
dmr_table <- rbind(Phenology, Alaska_W_C, Svalbard_W_C, Sweden_W_C, Nunavut_W_C, Wild_W_C, Wild_L_H, Parent_W_C, SE_W_C, SE_L_H)
############################################################################
#Attatch go_Term_names to genes
merge_go_names <- function(interproscan_terms, go_terms){
    go_term_names <- merge(interproscan_table, go_terms, by = "GO_Terms")
    go_term_names$Gene <- sub("V1.1$", "", go_term_names$Gene)
    collapsed_go_terms <- go_term_names %>% group_by(Gene) %>% summarize(GO_Terms = paste(GO_Terms, collapse = ","), GO_Name = paste(GO_Name, collapse = ","))
    print(head(collapsed_go_terms))
    return(data.frame(collapsed_go_terms))

}

collapsed_go_terms <- merge_go_names(interproscan_table, go_terms)

gene_dmr_table <- merge(dmr_table, dryas_genes, by = "Gene")
#Modify DMR table
gene_dmr_table$Origin <- sub("^original_data/", "", gene_dmr_table$Origin)
#write.table(gene_dmr_table, "Gene_DMR_Total_Merged_table2.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

##########################################################################
#Save dmrs with go table
##########################################################################

gene_dmr_go_table <- merge(gene_dmr_table, collapsed_go_terms, by = "Gene", all.x = TRUE)
gene_dmr_go_table <- distinct(gene_dmr_go_table)

# Gene Origin  Scaffold DMR_Start DMR_End q-value meanmethyl CpG  meanW  meanC Gene_Start Gene_End GO_Terms Go_Names

#write.table(gene_dmr_go_table, "Gene_DMR_Total_GO_Merged_table2.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

##########################################################################
#Read RNA 
read_rna <- function(rna_file) {
    rna_cols <- c("Gene", "logFC", "logCPM", "LR", "PValue", "fdr")
    rna_table <- read.table(rna_file, header=TRUE)
    colnames(rna_table) <- rna_cols
    rna_table$Gene <- sub("V1.1$", "", rna_table$Gene)
    print(head(rna_table))
    return(rna_table)
}

rna_table <- read_rna("original_data/RNA_DERS_Oct2023_W_C_Total_No_Headers.txt")

#######################################################################
#Attach GO_Terms, and names
#!!! Start end did not work!!!
rna_table <- merge(rna_table, dryas_genes, by = "Gene")
rna_go_table <- merge(rna_table, collapsed_go_terms, by = "Gene", all.x = TRUE)
print(head(rna_go_table, n=20))

#col names : Gene     logFC     logCPM        LR       PValue          fdr  Gene_Start Gene_End GO_Terms GO_Name

#######################################################################
#Edit RNA :
rna_go_table <- rna_go_table %>% mutate(Origin = "RNA")
rna_go_table <- rna_go_table %>% mutate(Scaffold = str_extract(Gene, "Do1_[^G]+"))
print(head(rna_go_table))
#write.table(rna_go_table, "Gene_RNA_Total_GO_Merged_table2.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

######################################################################
#Combine:
total_table <- bind_rows(gene_dmr_go_table, rna_go_table)
print(head(total_table))

#write.table(total_table, "Gene_DMR_RNA_GO_Merged_table.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

##############################################################
# Separate out the RNA data from the wild warming DMRs
# left join by gene IDs

print(head(rna_go_table))
print(head(gene_dmr_go_table))
DMRs_RNA <- left_join(gene_dmr_go_table, rna_go_table, by="Gene")
print(head(DMRs_RNA))

write.table(DMRs_RNA, "Gene_DMR_RNA_GO_Left_join_table.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

#-------------------------------------
# subset for only the DER and DMR genes/regions

DMRs_DER <- DMRs_RNA[which(!is.na(DMRs_RNA$Gene_Start.y)),]
print(head(DMRs_DER))

write.table(DMRs_DER, "DMR_DERs_Left_join_table.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

#------------------
unique(DMRs_DER$Origin.x)
 # [1] "Mat_Sen"      "SE_W_C"       "Alaska_W_C"   "Wild_Lat_L_H" "Nunavut_W_C"
 # [6] "Svalbard_W_C" "SE_L_H"       "Wild_W_C"     "Sweden_W_C"   "Parent_W_C"

DMRs_DER_Wild_W_C <- DMRs_DER[which(DMRs_DER$Origin.x=="Wild_W_C"),]
nrow(DMRs_DER_Wild_W_C)
# 44

DMRs_DER_SE_W_C <- DMRs_DER[which(DMRs_DER$Origin.x=="SE_W_C"),]
nrow(DMRs_DER_SE_W_C)
#130

DMRs_DER_Alaska_W_C <- DMRs_DER[which(DMRs_DER$Origin.x=="Alaska_W_C"),]
nrow(DMRs_DER_Alaska_W_C)
#120

DMRs_DER_Lat_L_H <- DMRs_DER[which(DMRs_DER$Origin.x=="Wild_Lat_L_H"),]
nrow(DMRs_DER_Lat_L_H)
#2676

DMRs_DER_Nunavut_W_C <- DMRs_DER[which(DMRs_DER$Origin.x=="Nunavut_W_C"),]
nrow(DMRs_DER_Nunavut_W_C)
#39

#-----------------
# difference between gene and DMR positions

print(head(DMRs_DER))
DMRs_DER$DMR_DER_diffstart <- abs(DMRs_DER$DMR_Start - DMRs_DER$Gene_Start.x)

DMRs_DER_closesubset <- DMRs_DER[which(DMRs_DER$DMR_DER_diffstart>20000),]

nrow(DMRs_DER_closesubset)
#1928

DMRs_DER_W_C_close <- DMRs_DER_closesubset[-which(DMRs_DER_closesubset$Origin.x=="Wild_Lat_L_H"|DMRs_DER_closesubset$Origin.x=="SE_L_H"|DMRs_DER_closesubset$Origin.x=="Mat_Sen"),]
nrow(DMRs_DER_W_C_close)
# 254

write.table(DMRs_DER_W_C_close, "DMR_DERs_Left_join_W_C_closesubset.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

#----------------------------------
# number of genes that fit close subset

DMRs_DER_Wild_W_C <- DMRs_DER_closesubset[which(DMRs_DER_closesubset$Origin.x=="Wild_W_C"),]
nrow(DMRs_DER_Wild_W_C)
# 6

DMRs_DER_SE_W_C <- DMRs_DER_closesubset[which(DMRs_DER_closesubset$Origin.x=="SE_W_C"),]
nrow(DMRs_DER_SE_W_C)
#25

DMRs_DER_Alaska_W_C <- DMRs_DER_closesubset[which(DMRs_DER_closesubset$Origin.x=="Alaska_W_C"),]
nrow(DMRs_DER_Alaska_W_C)
#53

DMRs_DER_Lat_L_H <- DMRs_DER_closesubset[which(DMRs_DER_closesubset$Origin.x=="Wild_Lat_L_H"),]
nrow(DMRs_DER_Lat_L_H)
#1403

DMRs_DER_Nunavut_W_C <- DMRs_DER_closesubset[which(DMRs_DER_closesubset$Origin.x=="Nunavut_W_C"),]
nrow(DMRs_DER_Nunavut_W_C)
#14

#---------------
#bash

grep "Wild_W_C" DMR_DERs_Left_join_table.tsv | sort -u  | uniq | wc -l 
# 44

###################################################################
#Get Phenogram with origin as phenotype
#total_table <- read.delim("Gene_DMR_RNA_GO_Merged_table.tsv", na.)
pheno_table <- total_table[, c("Scaffold", "Gene_Start", "Origin")]
print(head(pheno_table))
pheno_table <- pheno_table[!duplicated(pheno_table), ]
#result <- any(grepl("RNA", pheno_table$Origin))
#print(result)

#Extract true 
rows_to_keep <- grepl("^Do1_[[:alnum:]]+_.*$", pheno_table$Scaffold)
filtered_pheno_table <- pheno_table[rows_to_keep, ]
print(head(filtered_pheno_table))
colnames(filtered_pheno_table) <- c("CHR", "POS", "PHENOTYPE")
result <- any(grepl("RNA", filtered_pheno_table$PHENOTYPE))
print(result)
#write.table(filtered_pheno_table, "Dryas_Total_Origin_Phenogram.txt", sep = "\t", quote = FALSE, row.names = FALSE)

filtered_pheno_table$CHR <- gsub(".*0(\\d+)_a0*(\\d+)", "\\1\\2", filtered_pheno_table$CHR)
write.table(filtered_pheno_table, "Dryas_Total_Origin_Phenogram_Chr.txt", sep = "\t", quote = FALSE, row.names = FALSE)
###################################################################
#Isolate repeating values
phenogram_table <- read.delim(Dryas_Total_Origin_Phenogram_Chr.txt, sep = "\t")
sorted_phenogram_table <- phenogram_table[order(phenogram_table$CHR, phenogram_table$POS), ]

#duplicates <- duplicated(sorted_phenogram_table[, c("CHR", "POS")])
intersecting_phenogram_table <- sorted_phenogram_table[duplicated(sorted_phenogram_table[, c("CHR", "POS")]) | duplicated(sorted_phenogram_table[, c("CHR", "POS")], fromLast = TRUE), ]
print(head(intersecting_phenogram_table))
# Filter the data frame to keep only duplicated rows
write.table(intersecting_phenogram_table, "Dryas_Intersecting_Phenogram_Chr.txt", sep = "\t", quote = FALSE, row.names = FALSE)
# Display the duplicated data frame

sorted_W_C_pheno <- sorted_phenogram_table[grepl("SE_W_C|Wild_W_C|Parent_W_C|RNA", sorted_phenogram_table$PHENOTYPE), ]
intersecting_W_C_table <- sorted_W_C_pheno[duplicated(sorted_W_C_pheno[, c("CHR", "POS")]) | duplicated(sorted_W_C_pheno[, c("CHR", "POS")], fromLast = TRUE), ]
print(head(intersecting_phenogram_table))

write.table(intersecting_W_C_table, "Dryas_Intersecting_W_C_Phenogram_Chr.txt", sep = "\t", quote = FALSE, row.names = FALSE)