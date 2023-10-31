#October 2023
#module load r/4.2.1
#export R_LIBS_USER=/home/msandler/R/x86_64-pc-linux-gnu-library/4.2.1/
#Need some memory for Wild_Lat data:
#salloc -c1 --time 4:00:00 --mem 120000m --account def-rieseber
#In R
#To do:    
    #Add instructions on how to preprocess files 
#Can be run as script or just inputted into R

library(dplyr)
library(tidyr)
#Get Alaska file:

#####################################################################
#Constants:
metilene_colnames <- c("Scaffold", "DMR_Start", "DMR_End", "q-value", "meanmethyl", "CpG", "meanW", "meanC")
interproscanfile <- "original_data/dryas_goterm_file.tsv"
gotermsfile <- "original_data/interpro2go"

#####################################################################
#Process gff3 file:
process_gff3 <- function(file_name) {
    gff3_file <- read.table(file_name, header=FALSE)
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
dryas_genes <- process_gff3("Dryas_octopetala_H1.gff3")
print(head(dryas_genes))
#Should have start, end, gene
#####################################################################
#Read cleaned blast file:
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
#Merge metilene:
#Merge with metilene data:
    #Read metilene bedGraph (post filter)
    #Merge by scaffold and location
###In progress
merge_metilene <- function(blast_table, metilene_file_name) {
     metilene_table <- read.table(metilene_file_name, header=FALSE) 
    colnames(metilene_table) <- metilene_colnames
    print(head(metilene_table))    
    merged <- merge(blast_table, metilene_table, by = "Scaffold")
    merged <- merged %>% filter(DMR_Start >= Ext_Start, DMR_End <= Ext_End)
    merged <- merged[c("Origin", "Gene", "Scaffold", "DMR_Start", "DMR_End", "q-value", "meanmethyl", "CpG", "meanW", "meanC")]
    merged$Gene <- sub("V1.1$", "", merged$Gene)
    return(merged)
}

#Attach metilene info:
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
print("Metilene loaded")
############################################################################
#Works until here
###########################################################################

#Interproscan file processing
process_interproscan_file <- function(interpro_file){
    interpro_table <- read.table(interpro_file, header = FALSE)
    colnames(interpro_table) <- c("GeneName", "GO_Terms")
    interpro_table <- separate_rows(interpro_table, GO_Terms, sep = "\\|")
    no_dub_interpro_table <- distinct(interpro_table, GeneName, GO_Terms)
    no_dub_interpro_table <- data.frame(no_dub_interpro_table)
    print(head(no_dub_interpro_table))
    return(no_dub_interpro_table)
}

interproscan_table <- process_interproscan_file(interproscanfile)


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

##########################################################################
#Modify DMR table
gene_dmr_go_table$Origin <- sub("^original_data/", "", gene_dmr_go_table$Origin)
##########################################################################
gene_dmr_table <- merge(dmr_table, dryas_genes, by = "Gene")
write.table(gene_dmr_table, "Gene_DMR_Total_Merged_table.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
##########################################################################
#Save dmrs with go table

gene_dmr_go_table <- merge(gene_dmr_table, collapsed_go_terms, by = "Gene", all.x = TRUE)
gene_dmr_go_table <- distinct(gene_dmr_go_table)

# Gene Origin  Scaffold DMR_Start DMR_End q-value meanmethyl CpG  meanW  meanC Gene_Start Gene_End GO_Terms Go_Names


write.table(gene_dmr_go_table, "Gene_DMR_Total_GO_Merged_table.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
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
rna_table <- merge(rna_table, dryas_genes, by = "Gene")
rna_go_table <- merge(rna_table, collapsed_go_terms, by = "Gene", all.x = TRUE)
print(head(rna_go_table))
write.table(rna_go_table, "Gene_RNA_Total_GO_Merged_table.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

#col names : Gene     logFC     logCPM        LR       PValue          fdr  Gene_Start Gene_End GO_Terms GO_Name

#######################################################################