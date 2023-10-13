#module load r/4.2.1
#export R_LIBS_USER=/home/msandler/R/x86_64-pc-linux-gnu-library/4.2.1/
#tmux session cedar5, LeftMerge
#In R:

library(dplyr)
library(tidyr)
#Get Alaska file:


read_blast <- function(file_name) {
    title_string <- gsub("cleaned_blast_ref_(.*?)\\.out", "\\1", file_name)
    blast_table <- read.table("cleaned_blast_ref_Alaska_W_C.out", header=FALSE)
    subset_blast_table <- blast_table[, c(1, 2)]
    subset_blast_table <- subset_blast_table %>% mutate(Origin = title_string)
    colnames(subset_blast_table) <- c("DMR", "Gene", "Origin")
    separate(subset_blast_table, DMR, into = c("Scaffold", "Ext_Start", "Ext_End"), sep = "[:,\\-]")
   # print(head(subset_blast_table, n=10))
    return(data.frame(subset_blast_table))
}


Alaska_W_C <- read_blast("cleaned_blast_ref_Alaska_W_C.out") 
Svalbard_W_C <- read_blast("cleaned_blast_ref_Svalbard_W_C.out")
Nunavut_W_C <- read_blast("cleaned_blast_ref_Nunavut_W_C.out")
Sweden_W_C <- read_blast("cleaned_blast_ref_Sweden_W_C.out")
Wild_W_C <- read_blast("cleaned_blast_ref_Wild_W_C.out")
Parent_W_C <- read_blast("cleaned_blast_ref_Parent_W_C.out")
SE_W_C <- read_blast("cleaned_blast_ref_SE_W_C.out")

#Simple version of merge, not collapsed by origin, might delete later, just GeneID, DMR, 
combined_table <- rbind(Alaska_W_C, Svalbard_W_C, Sweden_W_C, Wild_W_C, Parent_W_C, SE_W_C)
write.table(combined_table, "Total_DMR_Gene_Origin_Not_Collapse_W_C.csv", sep = "\t", quote = FALSE, row.names = FALSE)

#Load gff3 into nice format
gff3_file <- read.table("Dryas_octopetala_H1.gff3", header=FALSE)
gff3_scaffold_start_end <- gff3_file[, c(1, 3, 4, 5, 9)]
dryas_genes <- gff3_scaffold_start_end[gff3_scaffold_start_end[, 2] == "gene", ]
dryas_genes <- dryas_genes[, c(3, 4, 5)]
gene_col_names <- c("Gene_Start", "Gene_End", "Gene")
colnames(dryas_genes) <- gene_col_names
#remove ID prefix:
dryas_genes$Gene <- sub("^ID=", "", dryas_genes$Gene)
write.table(dryas_genes, "Dryas_Genes.csv", sep = "\t", quote = FALSE, row.names = FALSE)


#For Alaska:
#attach metilene info

metilene_colnames <- colnames("Scaffold", "DMR_Start", "DMR_End", "metilene_val")


###Alaska example set
Alaska_Metilene <- read.table("original_data/Alaska_W_C.bedGraph", header=FALSE) 
colnames(Alaska_Metilene) <- metilene_colnames

library(tidyr)
#Add this to function above :
Alaska_W_C <- separate(Alaska_W_C, DMR, into = c("Scaffold", "Ext_Start", "Ext_End"), sep = "[:,\\-]")

Alaska <- merge(Alaska_W_C, Alaska_Metilene, by = "Scaffold")
Alaska <- Alaska %>% filter(DMR_Start >= Ext_Start, DMR_End <= Ext_End)
Alaska <- Alaska[c("Origin", "Gene", "Scaffold", "DMR_Start", "DMR_End", "metilene_val")]
Alaska$Gene <- sub("V1.1$", "", Alaska$Gene)
#Next account for the extension (see that DMR_Start < start, and end > DMR_End in start-end)

#Do this later with whole table:

AlaskaGeneID <- merge(Alaska, dryas_genes, by = "Gene")
write.table(AlaskaGeneID, "Alaska_Total_Merged.tsv", sep = "\t", quote = FALSE, row.names = TRUE)



###In progress
merge_metilene <- function(blast_table, metilene_file_name) {
    metilene_table <- read.table(metilene_file_name, header=FALSE) 
    colnames(metilene_table) <- metilene_colnames


    #Add this to function above :
    
    merged <- merge(blast_table, metilene_table, by = "Scaffold")
    merged <- merged %>% filter(DMR_Start >= Ext_Start, DMR_End <= Ext_End)
    merged <- merged[c("Origin", "Gene", "Scaffold", "DMR_Start", "DMR_End", "metilene_val")]
    merged$Gene <- sub("V1.1$", "", merged$Gene)
    return(merged)
}



##REDO COMBINE ALL THE TABLES BUT WITH DMR INFO: 
combined_table <- rbind(Alaska_W_C, Svalbard_W_C, Sweden_W_C, Wild_W_C, Parent_W_C, SE_W_C)
write.table(combined_table, "Total_Genes", sep = "\t", quote = FALSE, row.names = FALSE)



#Next: Merge by Gene Id for associated gene start/end