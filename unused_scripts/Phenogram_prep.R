##############################################################
#Load necessary libraries
library(stringr)
library(tidyverse)
##############################################################
#Set Global Variables
annotationfile <- "FINAL_Oxy_dig.AED_1.gff3"
interproscanfile <- "Oxy_goterm_file.tsv"
gotermsfile <- "goslim_plant.obo"
##############################################################
#Define functions:

#Process_annotation_file:
    #Arguments: path name to file with: start, end positions, and gene name
    #Returned value: Table of the form Start, End, GeneName
process_annotation_file <- function(ann_file) {
    annotation_table <- read.table(ann_file, header = FALSE)
    ncol(annotation_table) #9 columns
    genepos_table <- annotation_table[,c(9, 4, 5)]
    colnames(genepos_table) <- c("GeneName", "Start", "End")
    genenames <- str_extract(genepos_table$GeneName, "(?<=ID=)[^;-]+")
    genepos_table <- cbind(genepos_table, genenames)
    genepos_table <- subset(genepos_table, select = -GeneName)
    colnames(genepos_table) <- c("Start", "End", "GeneName")
    #print('The first 5 lines of your positions file are:')
    #print(head(genepos_table, n = 5))
    return(genepos_table)
}

#Make_interproscan_file
    #Argument: File of format Genename, PFamID, GO_TERM
    #Returns: table with interproscan values:
process_interproscan_file <- function(interpro_file){
    interpro_table <- read.table(interpro_file, header = FALSE)
    colnames(interpro_table) <- c("GeneName", "GO_Terms")
    interpro_table <- separate_rows(interpro_table, GO_Terms, sep = "\\|")
    no_dub_interpro_table <- distinct(interpro_table, GeneName, GO_Terms)
    no_dub_interpro_table <- data.frame(no_dub_interpro_table)
    #print(no_dub_interpro_table)
    #print('The first 5 lines of your interproscan table are:')
    #print(head(interpro_table, n = 5))
    return(no_dub_interpro_table)
}

split_go_terms <- function(goterm_file) {
    #load file
    goterm_file_lines <- readLines(goterm_file)
    #get only namespaces + go terms
    pattern <- "^(id: GO|namespace|name:)"
    rm_pattern <- c("external", "name: regulates", "name: occurs in", "part of", "during", "negatively regulates", "positively regulates", "has part", "name: term tracker item")
    filtered_goterms <- goterm_file_lines[grep(pattern, goterm_file_lines)]
    #print(filtered_goterms)
    for (l in filtered_goterms) {
        if (!any(grepl(paste(rm_pattern, collapse = "|"), l))) {
            filtered_goterms <- c(filtered_goterms, l)
        }
    }
    go_table <- data.frame(matrix(filtered_goterms, ncol = 3, byrow = TRUE))
    colnames(go_table) <- c("GO_Terms", "Name", "Domain")
    domains <- go_table$Domain

    bio_process <- data.frame(make_domain_table(go_table, domains, "biological_process"))
    cell_comp <- data.frame(make_domain_table(go_table, domains, "cellular_component"))
    mol_fxn <- data.frame(make_domain_table(go_table, domains, "molecular_function"))

    return(list(bio_process, cell_comp, mol_fxn))
}

#Helper fxn: Splits domain in obo file by keyword
make_domain_table <- function (go_table, domains, detect_string) {
    dv <- c()
    dv <-go_table %>% filter(str_detect(domains, detect_string))
    #print(paste0("The first 5 lines of the ",detect_string," table are: "))
    #print(head(dv, n = 5))
    goterms <- dv$GO_Terms
    dv$GO_Terms <- gsub("id: ", "", goterms)
    names <- dv$Name
    dv$Name <- gsub("name: ", "", names)
    return(dv)

}

#Subset interprofile:
    #Args interproscan table, list of ids to subset by
    #returns subsetted table
subset_interproscan_table <- function(interpro_table, gotable){
    ids <- gotable$GO_Terms
    #print("First 5 GO Terms ids being filtered for are: ")
    #print(head(ids))
    interpro_subset <- subset(interpro_table, interpro_table[, 3] %in% ids)
    print(nrow(interpro_subset))
    interpro_subset_name <- merge(interpro_subset, gotable, by = "GO_Terms", all.x=FALSE, all.y=FALSE)
    ##print(head(interpro_subset_name))
    ##print(interpro_subset_name$Name)
    return(interpro_subset_name)
}

#Make_pheno_table
    #Args: I
        #Interproscan table of format 
    #returns: phenogram table
make_pheno_table <- function(interproscan_subset, genepos_table) {
    #Subset table:
    #Merges Subsetted interproscan with Gene position table
    total_table <- merge(interproscan_subset, genepos_table, by = "GeneName", all.x = FALSE, all.y = FALSE)
    #Isolate Chromosome numbers
    total_table$Chr_Number <- str_extract(total_table$GeneName, "(?<=Chr0)[0-9]")
    no_dub_table <- distinct(total_table, GeneName, Name, .keep_all = TRUE)
    no_dub_table <- data.frame(no_dub_table)
    print(no_dub_table)
    #Removes any duplicates of Chromosome name + GO term
    
    #no_doubles_table <- total_table[!doubles_rows, ]
    #Creates correctly formatted table for Phenogram input
    ##print(no_doubles_table$Name)
    phenogram_pos_table <- no_dub_table[,c("Chr_Number", "Start", "Name")]
    colnames(phenogram_pos_table) <- c("CHR", "POS", "PHENOTYPE")
    return(phenogram_pos_table)
}

#make_pheno_file:
    #Args: 
        #A table compatible with Phenogram Input
        #filename
make_pheno_file <- function(pheno_table, file_name) {
    write.table(pheno_table, file = file_name, sep = "\t", quote = FALSE, row.names = FALSE)
}

#Subsets a ready pheno table by name:
    #Args: Ready Pheno_table
    #String to look for 
subset_by_name <- function(pheno_table, str) {
    filtered_table <- subset(pheno_table, grepl(str, PHENOTYPE))
    return(filtered_table)
}

subset_goterms_by_name <- function(goterms_list, str_list) {
    filtered_table <- data.frame()
    goterms_table <- data.frame(goterms_list)
    for (name in str_list) {
        per_term_rows <- goterms_table[grep(name, goterms_table$Name), ]
        filtered_table <- rbind(filtered_table, per_term_rows)
    }   
    return(filtered_table)

}

subset_interproscan_table <- function(interpro_table, gotable){
    ids <- gotable$GO_Terms
    #print("First 5 GO Terms ids being filtered for are: ")
    #print(head(ids))
    interpro_subset <- subset(interpro_table, interpro_table[, 2] %in% ids)
    print(nrow(interpro_subset))
    interpro_subset_name <- merge(interpro_subset, gotable, by = "GO_Terms", all.x=FALSE, all.y=FALSE)
    #print("The first 5 lines of your interproscan subset table are")
    #print(head(interpro_subset_name))
    ##print(interpro_subset_name$Name)
    return(interpro_subset_name)
}



##############################################################
#Make input file:
#Creates table of Start, end, Gene Name
#Creates table of Gene Name, PfamId, GO_Term

#  return(list(mean = mean_val, median = median_val, maximum = max_val))
#Process files:
genepos_table <- process_annotation_file(annotationfile)
interproscan_table <- process_interproscan_file(interproscanfile)
#Splits into domains:
split_domain <- split_go_terms(gotermsfile)

bio_process <- split_domain[[1]]
cell_comp <- split_domain[[2]]
mol_fxn <- split_domain[[3]]



#Call Split mol_fxn

#mol_fxn$Go_Term <- mol_fxn[mol_fxn != "GO:0005515"]

#Subset Interproscan by domain
interproscan_bioproc <- subset_interproscan_table(interproscan_table, bio_process)
interproscan_cellcomp <- subset_interproscan_table(interproscan_table, cell_comp)
interproscan_molfxn <- subset_interproscan_table(interproscan_table, mol_fxn)

#Make tables
bioproc_table <- make_pheno_table(interproscan_bioproc, genepos_table)
cellcomp_table <- make_pheno_table(interproscan_cellcomp, genepos_table)
molfxn_table <- make_pheno_table(interproscan_molfxn, genepos_table)
#print(paste0("Bioprocess has: ",nrow(bioproc_table)))
#print(paste0("Cellcomp has: ",nrow(cellcomp_table)))
#print(paste0("Molfxn has: ",nrow(molfxn_table)))
#sum <- nrow(bioproc_table)+nrow(cellcomp_table)+nrow(molfxn_table)
##print(paste0("For a total of: ", sum))

molfxn_activity <- subset_by_name(molfxn_table, "activity") 
molfxn_binding <- subset_by_name(molfxn_table, "binding") 


#Write files
make_pheno_file(bioproc_table, "Bio_Pro_Phenogram_input.txt")
make_pheno_file(cellcomp_table, "Cell_Comp_Phenogram_input.txt")
make_pheno_file(molfxn_activity, "Mol_Fxn_Activity_Phenogram_input.txt")
make_pheno_file(molfxn_binding, "Mol_Fxn_Binding_Phenogram_input.txt")

goterms_of_interest <- c("response to biotic stimulus", "response to stress")
stress_go_terms <- subset_goterms_by_name(bio_process, goterms_of_interest)

interproscan_stress <- subset_interproscan_table(interproscan_table, stress_go_terms)
stress_table <- make_pheno_table(interproscan_stress, genepos_table)
make_pheno_file(stress_table, "Stress_Stimuli_Phenogram_Oxyria_input.txt")