
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

metilene_colnames <- c("Scaffold", "DMR_Start", "DMR_End", "metilene_val")


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
    print(head(metilene_table))
    print(head(blast_table))

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
#####################################################################

#Cleanly:

read_blast <- function(file_name) {
    title_string <- gsub("cleaned_blast_ref_(.*?)\\.out", "\\1", file_name)
    blast_table <- read.table("cleaned_blast_ref_Alaska_W_C.out", header=FALSE)
    subset_blast_table <- blast_table[, c(1, 2)]
    subset_blast_table <- subset_blast_table %>% mutate(Origin = title_string)
    colnames(subset_blast_table) <- c("DMR", "Gene", "Origin")
    subset_blast_table <- separate(subset_blast_table, DMR, into = c("Scaffold", "Ext_Start", "Ext_End"), sep = "[:,\\-]")
    print(head(subset_blast_table, n=10))
    return(data.frame(subset_blast_table))
}


merge_metilene <- function(blast_table, metilene_file_name) {
    metilene_table <- read.table(metilene_file_name, header=FALSE) 
    colnames(metilene_table) <- metilene_colnames
    print(head(metilene_table))

    #Add this to function above :
    
    merged <- merge(blast_table, metilene_table, by = "Scaffold")
    merged <- merged %>% filter(DMR_Start >= Ext_Start, DMR_End <= Ext_End)
    merged <- merged[c("Origin", "Gene", "Scaffold", "DMR_Start", "DMR_End", "metilene_val")]
    merged$Gene <- sub("V1.1$", "", merged$Gene)
    return(merged)
}

metilene_colnames <- c("Scaffold", "DMR_Start", "DMR_End", "metilene_val")


Alaska_W_C <- read_blast("cleaned_blast_ref_Alaska_W_C.out") 
Svalbard_W_C <- read_blast("cleaned_blast_ref_Svalbard_W_C.out")
Nunavut_W_C <- read_blast("cleaned_blast_ref_Nunavut_W_C.out")
Sweden_W_C <- read_blast("cleaned_blast_ref_Sweden_W_C.out")
Wild_W_C <- read_blast("cleaned_blast_ref_Wild_W_C.out")
Parent_W_C <- read_blast("cleaned_blast_ref_Parent_W_C.out")
SE_W_C <- read_blast("cleaned_blast_ref_SE_W_C.out")

Alaska <- merge_metilene(Alaska_W_C, "original_data/Alaska_W_C.bedGraph") 
Svalbard <- merge_metilene(Svalbard_W_C, "original_data/Svalbard_W_C.bedGraph")
Nunavut <- merge_metilene(Nunavut_W_C, "original_data/Nunavut_W_C.bedGraph")
Sweden <- merge_metilene(Sweden_W_C, "original_data/Sweden_W_C.bedGraph")
Wild <- merge_metilene(Wild_W_C, "original_data/Wild_W_C.bedGraph")
Parent <- merge_metilene(Parent_W_C, "original_data/Parent_W_C.bedGraph")
SE <- merge_metilene(SE_W_C, "original_data/SE_W_C.bedGraph")

dmr_table <- rbind(Alaska, Svalbard, Sweden, Nunavut, Wild, Parent, SE)
gene_dmr_table <- merge(dmr_table, dryas_genes, by = "Gene")
write.table(gene_dmr_table, "Gene_DMR_Merged_table.tsv", sep = "\t", quote = FALSE, row.names = FALSE)



process_interproscan_file <- function(interpro_file){
    interpro_table <- read.table(interpro_file, header = FALSE)
    colnames(interpro_table) <- c("GeneName", "GO_Terms")
    interpro_table <- separate_rows(interpro_table, GO_Terms, sep = "\\|")
    no_dub_interpro_table <- distinct(interpro_table, GeneName, GO_Terms)
    no_dub_interpro_table <- data.frame(no_dub_interpro_table)
    return(no_dub_interpro_table)
}


process_go_terms <- function(goterm_file) {
    goterm_file_lines <- readLines(goterm_file)
    pattern <- c("GO:")
    filtered_goterm_lines <- goterm_file_lines[grep(pattern, goterm_file_lines)]
    filtered_goterms <- data.frame(filtered_goterm_lines)
    colnames(filtered_goterms) <- c("Lines")
    split_file <- data.frame(separate(filtered_goterms, Lines, into = c("IPRscan", "Name", "GO_Terms"), sep = ">|;"))
    goterms <- split_file[,c("Name", "GO_Terms")]
    goterms$GO_Terms <- gsub("\\s", "", goterms$GO_Terms)
    no_dub_goterms <- distinct(goterms, GO_Terms, Name)
    return(data.frame(no_dub_goterms))
    #print(head(split_file))
}

interproscanfile <- "original_data/dryas_goterm_file.tsv"
gotermsfile <- "original_data/interpro2go"

interproscan_table <- process_interproscan_file(interproscanfile)
go_terms <- process_go_terms(gotermsfile)

go_term_names <- merge(interproscan_table, go_terms, by = "GO_Terms")
go_term_names <- go_term_names[, c("GeneName", "Name")]
colnames(go_term_names) <- c("Gene", "GO_Terms")
go_term_names$Gene <- sub("V1.1$", "", go_term_names$Gene)
collapsed_go_terms <- go_term_names %>% group_by(Gene) %>% summarize(GO_Terms = paste(GO_Terms, collapse = ","))


gene_dmr_go_table <- merge(gene_dmr_table, collapsed_go_terms, by = "Gene", all.x = TRUE)
gene_dmr_go_table <- distinct(gene_dmr_go_table)
write.table(gene_dmr_go_table, "Gene_DMR_GO_Merged_table.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

####################################################
#Testing full metilene:

read_blast <- function(file_name) {
    title_string <- gsub("cleaned_blast_ref_(.*?)\\.out", "\\1", file_name)
    blast_table <- read.table("cleaned_blast_ref_Alaska_W_C.out", header=FALSE)
    subset_blast_table <- blast_table[, c(1, 2)]
    subset_blast_table <- subset_blast_table %>% mutate(Origin = title_string)
    colnames(subset_blast_table) <- c("DMR", "Gene", "Origin")
    subset_blast_table <- separate(subset_blast_table, DMR, into = c("Scaffold", "Ext_Start", "Ext_End"), sep = "[:,\\-]")
   # print(head(subset_blast_table, n=10))
    return(data.frame(subset_blast_table))
}

#Site specific
Alaska_W_C_blast <- read_blast("original_data/cleaned_blast_ref_Alaska_W_C.out ") 
Svalbard_W_C_blast <- read_blast("original_data/cleaned_blast_ref_Svalbard_W_C.out")
Nunavut_W_C_blast <- read_blast("original_data/cleaned_blast_ref_Nunavut_W_C.out")
Sweden_W_C_blast <- read_blast("original_data/cleaned_blast_ref_Sweden_W_C.out")
#Wild
Wild_W_C_blast <- read_blast("original_data/cleaned_blast_ref_Wild_W_C.out")
Parent_W_C_blast <- read_blast("original_data/cleaned_blast_ref_Parent_W_C.out")
Wild_DO_DI_blast <- read.blast("original_data/cleaned_blast_ref_Wild_Species_DO_DI.out")
Wild_L_H_blast <-read_blast("original_data/cleaned_blast_ref_Wild_Lat_Low_Hi.out")


SE_W_C <- read_blast("original_data/cleaned_blast_ref_SE_W_C.out")

Phenology <- merge_metilene(Pheno_Mat_Sen, "original_data/ Mat_Sen_70_5_4_0.9_qval.0.001.out")

Alaska <- merge_metilene(Alaska_W_C, "original_data/Alaska_W_C_70_5_4_0.9_qval.0.001.out") 
Svalbard <- merge_metilene(Svalbard_W_C, "original_data/Svalbard_W_C.bedGraph")
Nunavut <- merge_metilene(Nunavut_W_C, "original_data/Nunavut_W_C.bedGraph")
Sweden <- merge_metilene(Sweden_W_C, "original_data/Sweden_W_C.bedGraph")
Wild <- merge_metilene(Wild_W_C, "original_data/Wild_W_C.bedGraph")
Parent <- merge_metilene(Parent_W_C, "original_data/Parent_W_C.bedGraph")
SE <- merge_metilene(SE_W_C, "original_data/SE_W_C.bedGraph"



 Mat_Sen_70_5_4_0.9_qval.0.001.out      SE_W_C_70_5_4_0.9_qval.0.001.out        Wild_Species_DO_DI_70_5_4_0.9_qval.0.001.out
dryas_goterm_file.tsv                 Nunavut_W_C_70_5_4_0.9_qval.0.001.out  Svalbard_W_C_70_5_4_0.9_qval.0.001.out  Wild_W_C_70_5_4_0.9_qval.0.001.out
interpro2go                           P_W_C_70_5_4_0.9_qval.0.001.out        Sweden_W_C_70_5_4_0.9_qval.0.001.out
interproscan_dryas_full.tsv           SE_L_H_70_5_4_0.9_qval.0.001.out       Wild_Lat_L_H_70_5_4_0.9_qval.0.001.out


#Note these values may be wrong
metilene_colnames <- c("Scaffold", "DMR_Start", "DMR_End", "q-value", "meanmethyl", "CpG", "meanW", "meanC") #Percent?

metilene_table <- read.table("original_data/Mat_Sen_70_5_4_0.9_qval.0.001.out", header=FALSE) 
colnames(metilene_table) <- metilene_colnames
 print(head(metilene_table))



#########################

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

#Phenology
Phenology_blast <- read_blast("original_data/cleaned_blast_ref_Mat_Sen.out")
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

metilene_colnames <- c("Scaffold", "DMR_Start", "DMR_End", "q-value", "meanmethyl", "CpG", "meanW", "meanC") 

merge_metilene <- function(blast_table, metilene_file_name) {
    metilene_table <- read.table(metilene_file_name, header=FALSE) 
    colnames(metilene_table) <- metilene_colnames
    print(head(metilene_table))

    #Add this to function above :
    
    merged <- merge(blast_table, metilene_table, by = "Scaffold")
    merged <- merged %>% filter(DMR_Start >= Ext_Start, DMR_End <= Ext_End)
    merged <- merged[c("Origin", "Gene", "Scaffold", "DMR_Start", "DMR_End", "q-value", "meanmethyl", "CpG", "meanW", "meanC")]
    merged$Gene <- sub("V1.1$", "", merged$Gene)
    return(merged)
}


Phenology <- merge_metilene(Phenology, "original_data/Mat_Sen_70_5_4_0.9_qval.0.001.out")

#site specific places
Alaska_W_C <- merge_metilene(Alaska_W_C_blast, "original_data/Alaska_W_C_70_5_4_0.9_qval.0.001.out") 
Svalbard_W_C <- merge_metilene(Svalbard_W_C_blast, "original_data/Svalbard_W_C_70_5_4_0.9_qval.0.001.out")
Nunavut_W_C <- merge_metilene(Nunavut_W_C_blast, "original_data/Nunavut_W_C_70_5_4_0.9_qval.0.001.out")
Sweden_W_C <- merge_metilene(Sweden_W_C_blast, "original_data/Sweden_W_C_70_5_4_0.9_qval.0.001.out")
#Wild
Wild_W_C <- merge_metilene(Wild_W_C_blast, "original_data/Wild_W_C_70_5_4_0.9_qval.0.001.out")
Wild_L_H <- merge_metilene(Wild_L_H_blast, "original_data/Wild_Lat_L_H_70_5_4_0.9_qval.0.001.out")
Parent_W_C <- merge_metilene(Parent_W_C_blast, "original_data/P_W_C_70_5_4_0.9_qval.0.001.out")
#Seedlings:
SE_W_C <- merge_metilene(Seedling_W_C_blast, "original_data/SE_W_C_70_5_4_0.9_qval.0.001.out")
SE_L_H <- merge_metilene(Seedling_L_H_blast, "original_data/SE_L_H_70_5_4_0.9_qval.0.001.out")

dmr_table <- rbind(Phenology, Alaska_W_C, Svalbard_W_C, Sweden_W_C, Nunavut_W_C, Wild_W_C, Wild_L_H, Parent_W_C, SE_W_C, SE_L_H)

process_gff3 <- function(file_name) {
    gff3_file <- read.table("Dryas_octopetala_H1.gff3", header=FALSE)
    gff3_scaffold_start_end <- gff3_file[, c(1, 3, 4, 5, 9)]
    dryas_genes <- gff3_scaffold_start_end[gff3_scaffold_start_end[, 2] == "gene", ]
    dryas_genes <- dryas_genes[, c(3, 4, 5)]
    gene_col_names <- c("Gene_Start", "Gene_End", "Gene")
    colnames(dryas_genes) <- gene_col_names
    #remove ID prefix:
    dryas_genes$Gene <- sub("^ID=", "", dryas_genes$Gene)
    write.table(dryas_genes, "Dryas_Genes.csv", sep = "\t", quote = FALSE, row.names = FALSE)
    return(dryas_genes)
}

process_interproscan_file <- function(interpro_file){
    interpro_table <- read.table(interpro_file, header = FALSE)
    colnames(interpro_table) <- c("GeneName", "GO_Terms")
    interpro_table <- separate_rows(interpro_table, GO_Terms, sep = "\\|")
    no_dub_interpro_table <- distinct(interpro_table, GeneName, GO_Terms)
    no_dub_interpro_table <- data.frame(no_dub_interpro_table)
    return(no_dub_interpro_table)
}



process_go_terms <- function(goterm_file) {
    goterm_file_lines <- readLines(goterm_file)
    pattern <- c("GO:")
    filtered_goterm_lines <- goterm_file_lines[grep(pattern, goterm_file_lines)]
    filtered_goterms <- data.frame(filtered_goterm_lines)
    colnames(filtered_goterms) <- c("Lines")
    split_file <- data.frame(separate(filtered_goterms, Lines, into = c("IPRscan", "Name", "GO_Terms"), sep = ">|;"))
    goterms <- split_file[,c("Name", "GO_Terms")]
    goterms$GO_Terms <- gsub("\\s", "", goterms$GO_Terms)
    no_dub_goterms <- distinct(goterms, GO_Terms, Name)
    return(data.frame(no_dub_goterms))
    #print(head(split_file))
}

merge_go_names <- function(interproscan_terms, go_terms){
    go_term_names <- merge(interproscan_table, go_terms, by = "GO_Terms")
    go_term_names <- go_term_names[, c("GeneName", "Name")]
    colnames(go_term_names) <- c("Gene", "GO_Terms")
    go_term_names$Gene <- sub("V1.1$", "", go_term_names$Gene)
    collapsed_go_terms <- go_term_names %>% group_by(Gene) %>% summarize(GO_Terms = paste(GO_Terms, collapse = ","))
    return(collapsed_go_terms)

}



interproscanfile <- "original_data/dryas_goterm_file.tsv"
gotermsfile <- "original_data/interpro2go"
dryas_genes <- process_gff3("Dryas_octopetala_H1.gff3")

read_rna <- function(rna_file) {
    rna_cols <- c("Gene", "logFC", "logCPM", "LR", "PValue", "fdr")
    rna_table <- read.table(rna_file, header=TRUE)
    colnames(rna_table) <- rna_cols
    print(head(rna_table))
    return(rna_table)
}

rna_table <- read_rna("original_data/RNA_DERS_Oct2023_W_C_Total_No_Headers.txt")

gene_dmr_rna_table <- merge(gene_dmr_rna_table, rna_table, by = "Gene")