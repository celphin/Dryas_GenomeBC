###########################
# Gene ontology analysis
# Nov 2024 
# ErmineJ
##########################

# copy over relevant files
# beluga3
tmux new-session -s GO
tmux attach-session -t GO

mkdir /home/celphin/scratch/Dryas/GO_enrichment

cd /home/celphin/scratch/Dryas/MS_Dryas_Merged_Data/original_data
cp interproscan_dryas_full.tsv /home/celphin/scratch/Dryas/GO_enrichment/
cp Dryas_octopetala_H1.gff3 /home/celphin/scratch/Dryas/GO_enrichment/
cd ..
cp Gene_DMR_RNA_GO_Merged_table.tsv /home/celphin/scratch/Dryas/GO_enrichment/
cd /home/celphin/scratch/Dryas/GO_enrichment

grep mRNA Dryas_octopetala_H1.gff3 | wc -l
# 41181

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

#-------------------------
# load GO ont data

module load  StdEnv/2020 r/4.2.2
export R_LIBS_USER=/home/celphin/R/x86_64-pc-linux-gnu-library/4.1.0/

R

path="/home/celphin/scratch/Dryas/GO_enrichment"
Gene_ont_file <- "interproscan_dryas_full2.tsv"
gene_ont <- read.delim(paste0(path,"/", Gene_ont_file), header = FALSE, sep = "\t", na.strings = "-", colClasses = c("character", "character", "character", "character"))

colnames(gene_ont) <- c( "gene", "Pfam", "descrip1", "INTPRO", "descrip2", "GOterm")
length(unique(gene_ont$INTPRO))
# 8644

nrow(gene_ont)
#  95869 # getting read in properly, half of the rows missing before

#-----------------------------
# formatting Interproscan to have no duplicates of genes - one row per gene
library(dplyr)
library(tidyr)

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

##########################
# GENE IDs for DMRs and DEGs

cd /home/celphin/scratch/Dryas/GO_enrichment
wc -l Gene_DMR_RNA_GO_Merged_table.tsv
# 51805 Gene_DMR_RNA_GO_Merged_table.tsv

#----------------------------------
tmux new-session -s GO
tmux attach-session -t GO

module load  StdEnv/2020 r/4.2.2
export R_LIBS_USER=/home/celphin/R/x86_64-pc-linux-gnu-library/4.1.0/

R
library(dplyr)
library(tidyr)

path="/home/celphin/scratch/Dryas/GO_enrichment"
gene_ont <- read.delim(paste0(path,"/interproscan_dryas_full3.tsv"), header = TRUE, sep = "\t", na.strings = "-")
DMR_DEG <- read.delim(paste0(path,"/Gene_DMR_RNA_GO_Merged_table.tsv"), header = TRUE, sep = "\t")

#Get a list of the specific genes for each group 
head(DMR_DEG)

             # Gene       Origin  Scaffold DMR_Start DMR_End    q.value
# 1 Do1_00107G00001      Mat_Sen Do1_00107     43844   44183 1.2495e-17
# 2 Do1_00107G00001      Mat_Sen Do1_00107    457781  457843 1.9232e-04
# 3 Do1_00107G00001 Wild_Lat_L_H Do1_00107      5214    5389 1.9126e-10
# 4 Do1_00107G00001      Mat_Sen Do1_00107    309861  309911 1.3805e-05
# 5 Do1_00107G00001 Wild_Lat_L_H Do1_00107    442512  442804 9.1255e-44
# 6 Do1_00107G00001 Wild_Lat_L_H Do1_00107    457808  457939 9.3147e-04
  # meanmethyl CpG  meanW  meanC Gene_Start Gene_End GO_Terms GO_Name logFC
# 1   4.311270  41 94.651 90.340       3059     4636     <NA>    <NA>    NA
# 2  -8.347944  12 81.937 90.285       3059     4636     <NA>    <NA>    NA
# 3 -13.033270  22 77.568 90.601       3059     4636     <NA>    <NA>    NA
# 4  10.383184  17 92.144 81.760       3059     4636     <NA>    <NA>    NA
# 5 -28.801401  20 36.889 65.690       3059     4636     <NA>    <NA>    NA
# 6  -5.939124  14 84.956 90.895       3059     4636     <NA>    <NA>    NA
  # logCPM LR PValue fdr
# 1     NA NA     NA  NA
# 2     NA NA     NA  NA
# 3     NA NA     NA  NA
# 4     NA NA     NA  NA
# 5     NA NA     NA  NA
# 6     NA NA     NA  NA

unique(DMR_DEG$Origin)

# "Wild_W_C"
DMR_Wild_W_C <- DMR_DEG[which(DMR_DEG$Origin=="Wild_W_C"),]
Wild_W_C <- cbind(DMR_Wild_W_C$Gene, DMR_Wild_W_C$q.value)
write.table(Wild_W_C, "Wild_W_C_genes.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# "Mat_Sen"
DMR_Mat_Sen <- DMR_DEG[which(DMR_DEG$Origin=="Mat_Sen"),]
Mat_Sen <- cbind(DMR_Mat_Sen$Gene, DMR_Mat_Sen$q.value)
write.table(Mat_Sen, "Mat_Sen_genes.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# "Wild_Lat_L_H"
DMR_Wild_Lat_L_H <- DMR_DEG[which(DMR_DEG$Origin=="Wild_Lat_L_H"),]
Wild_Lat_L_H <- cbind(DMR_Wild_Lat_L_H$Gene, DMR_Wild_Lat_L_H$q.value)
write.table(Wild_Lat_L_H, "Wild_Lat_L_H_genes.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# "Sweden_W_C"
DMR_Sweden_W_C <- DMR_DEG[which(DMR_DEG$Origin=="Sweden_W_C"),]
Sweden_W_C <- cbind(DMR_Sweden_W_C$Gene, DMR_Sweden_W_C$q.value)
write.table(Sweden_W_C, "Sweden_W_C_genes.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# "Alaska_W_C"
DMR_Alaska_W_C <- DMR_DEG[which(DMR_DEG$Origin=="Alaska_W_C"),]
Alaska_W_C <- cbind(DMR_Alaska_W_C$Gene, DMR_Alaska_W_C$q.value)
write.table(Alaska_W_C, "Alaska_W_C_genes.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# "Nunavut_W_C"
DMR_Nunavut_W_C <- DMR_DEG[which(DMR_DEG$Origin=="Nunavut_W_C"),]
Nunavut_W_C <- cbind(DMR_Nunavut_W_C$Gene, DMR_Nunavut_W_C$q.value)
write.table(Nunavut_W_C, "Nunavut_W_C_genes.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# "Svalbard_W_C"
DMR_Svalbard_W_C <- DMR_DEG[which(DMR_DEG$Origin=="Svalbard_W_C"),]
Svalbard_W_C <- cbind(DMR_Svalbard_W_C$Gene, DMR_Svalbard_W_C$q.value)
write.table(Svalbard_W_C, "Svalbard_W_C_genes.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)


# "Parent_W_C"
DMR_Parent_W_C <- DMR_DEG[which(DMR_DEG$Origin=="Parent_W_C"),]
Parent_W_C <- cbind(DMR_Parent_W_C$Gene, DMR_Parent_W_C$q.value)
write.table(Parent_W_C, "Parent_W_C_genes.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# "SE_W_C"
DMR_SE_W_C <- DMR_DEG[which(DMR_DEG$Origin=="SE_W_C"),]
SE_W_C <- cbind(DMR_SE_W_C$Gene, DMR_SE_W_C$q.value)
write.table(SE_W_C, "SE_W_C_genes.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# "SE_L_H"
DMR_SE_L_H <- DMR_DEG[which(DMR_DEG$Origin=="SE_L_H"),]
SE_L_H <- cbind(DMR_SE_L_H$Gene, DMR_SE_L_H$q.value)
write.table(SE_L_H, "SE_L_H_genes.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# "RNA"
DMR_RNA <- DMR_DEG[which(DMR_DEG$Origin=="RNA"),]
RNA <- cbind(DMR_RNA$Gene, DMR_RNA$fdr)
write.table(RNA, "RNA_genes.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)


#-------------------------
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

#--------------------------------
# in bash
# Prepare gene score file 
# https://erminej.msl.ubc.ca/help/input-files/gene-scores/

awk 'BEGIN{FS="\t"}{print $1,"1"}' total_gene_list.txt | sort -u > Dryas_blank_geneset

for taxon in Wild_W_C Mat_Sen Wild_Lat_L_H Sweden_W_C Alaska_W_C Nunavut_W_C Svalbard_W_C RNA Parent_W_C SE_W_C SE_L_H ; \
do \
echo "$taxon"
cp Dryas_blank_geneset "$taxon"_geneset

# Loop through the gene list for each taxon, extract the qvalue and replace 1 with the qvalue in the gene set file
sort "$taxon"_genes.txt | \
while read gene qvalue ; do \
    # Replace 1 with the qvalue in the gene set file
    sed -i "s/$gene 1/$gene $qvalue/g" "$taxon"_geneset ; \
done ; done

for taxon in Wild_W_C Mat_Sen Wild_Lat_L_H Sweden_W_C Alaska_W_C Nunavut_W_C Svalbard_W_C RNA Parent_W_C SE_W_C SE_L_H ; \
do sed -i 's/ /\t/g' "$taxon"_geneset; done

# remove V1.1 in GO_gene_mappings.ermineJ.txt
sed 's/V1\.1//g' GO_gene_mappings.ermineJ.txt > GO_gene_mappings1.ermineJ.txt

############################################
# get counts
cd /home/celphin/scratch/Dryas/GO_enrichment
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


##############################################
# ermineJ

# https://erminej.msl.ubc.ca/help/tutorials/running-an-analysis-ora/

# As of ErmineJ 3, when using the ‘ORA’ method you have the option to use a simple “hit list” of genes,
# rather than preparing a score file yourself (a “quick list”). Caution: If you use this feature, 
# the “non-hits” will be all the rest of the genes listed in your annotation file. That might not 
# be appropriate if the annotation file includes genes that were not assayed in your experiment. 
# This is most likely to be a problem if your annotation file is a list of all the genes in the genome

tmux new-session -s GO
tmux attach-session -t GO

salloc -c1 --time 2:00:00 --mem 120000m --account def-cronk

cd /home/celphin/scratch/Dryas/GO_enrichment

ERMINEJ_HOME=/home/celphin/ermineJ-3.2
export JAVA_HOME=/cvmfs/soft.computecanada.ca/easybuild/software/2020/Core/java/13.0.2/

module load java/13.0.2

for taxon in Wild_W_C Mat_Sen Wild_Lat_L_H Sweden_W_C Alaska_W_C Nunavut_W_C Svalbard_W_C RNA Parent_W_C SE_W_C SE_L_H ; \
do $ERMINEJ_HOME/bin/ermineJ.sh \
-a GO_gene_mappings1.ermineJ.txt \
-s "$taxon"_geneset \
-c /home/celphin/ermineJ.data/go.obo \
-n ORA -t 0.01 \
--genesOut -aspects BCM \
-o "$taxon".ermine.results -y 5 ; done

# Starting ORA analysis
# Hit list (182 genes) enrichment for multifunctionality: P = 0.998
# 300 gene sets analyzed ...
# 600 gene sets analyzed ...
# 900 gene sets analyzed ...
# 1200 gene sets analyzed ...
# 1500 gene sets analyzed ...
# 'Hits' are not significantly multifunctionality-biased, no multifunctionality correction needed
# Finished with ORA computations: 182 elements passed your threshold.
# Multiple test correction for 1667 scored sets.
# Multifunctionality correlation is 0.01 for 14954 values
# Done!
# [Gemma 2024-12-01 23:06:53,276] INFO [main] ubic.erminej.ResultsPrinter.getDestination(177) | Writing results t
# o Wild_W_C.ermine.results

############################
# Results exploration

# Look at the files for each set combined
# Subset those with a p-value < 0.05
for taxon in Wild_W_C Mat_Sen Wild_Lat_L_H Sweden_W_C Alaska_W_C Nunavut_W_C Svalbard_W_C RNA Parent_W_C SE_W_C SE_L_H ; \
do grep "!" "$taxon".ermine.results  | awk -F '\t' '$7 <0.5 { print }'| awk  -F $'\t' '{print $2, $3, $7, $1}' |sort  > "$taxon"_sig.ermine.results
done

# Extract just GO terms and count duplicates
cat *_sig.ermine.results | awk -F '\t' '$7 <0.5 { print }' | awk  -F $'\t' '{print $3, $2}' | sort | uniq -c  
# 3404

###################################
# summarize with : http://revigo.irb.hr/

for taxon in Wild_W_C Mat_Sen Wild_Lat_L_H Sweden_W_C Alaska_W_C Nunavut_W_C Svalbard_W_C RNA Parent_W_C SE_W_C SE_L_H ; \
do cat "$taxon".ermine.results | awk -F '\t' '$7 <0.05 { print }' | awk  -F $'\t' '{print $3}' | \
sort  | uniq -c | awk  -F " " '{print $2, $1}' > Revigio_"$taxon".txt 
done

###################################
# Look at all and what species they are from
cd /home/celphin/scratch/Dryas/GO_enrichment

grep "!" *totalfam*_genesets.ermine.results  | awk -F '\t' '$7 <1 { print }'| awk  -F $'\t' '{print $1, $2, $3, $7}' |sort  > List_all_GO.txt


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










