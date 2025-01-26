###########################
# snpEff to get gene features and TE overlap
# https://github.com/al2na/methylKit
# Jan 2025
##############################

# get feature overlap with SNPEff and TEs

tmux new-session -s snpEff1
tmux attach-session -t snpEff1

mkdir /lustre04/scratch/celphin/Dryas/snpEff/methylkit
cp /lustre04/scratch/celphin/Dryas_large_folders/intersections/Methylkit_*bedGraph /lustre04/scratch/celphin/Dryas/snpEff/methylkit
cp /lustre04/scratch/celphin/Dryas_large_folders/intersections/Rand_*bedGraph /lustre04/scratch/celphin/Dryas/snpEff/methylkit

cp /lustre04/scratch/celphin/Dryas_large_folders/intersections/Methylkit_Pheno_10_DMRs_CHH.bedGraph /lustre04/scratch/celphin/Dryas/snpEff/methylkit
cp /lustre04/scratch/celphin/Dryas_large_folders/intersections/Methylkit_Pheno_25_DMRs_CHH.bedGraph /lustre04/scratch/celphin/Dryas/snpEff/methylkit

cd /lustre04/scratch/celphin/Dryas/snpEff/methylkit/

FILES=$(ls /lustre04/scratch/celphin/Dryas/snpEff/methylkit/*.bedGraph)
echo ${FILES}


for file in ${FILES}; do
    # Run the command on the file and create the .bed file
    grep -v track "$file" | awk '{print $1 "\t" $2 "\t" $3}' > "${file%.bedGraph}.bed"
    echo "Processed $file"
done

# run SNPEff
module load StdEnv/2023 java/21.0.1
mkdir gene_features
cd /lustre04/scratch/celphin/Dryas/snpEff/methylkit/gene_features/

FILES=$(ls /lustre04/scratch/celphin/Dryas/snpEff/methylkit/*.bed)
echo ${FILES}

cd /lustre04/scratch/celphin/Dryas/snpEff

for file in ${FILES}; do
java -Xmx8g -jar snpEff.jar -i bed OldDoct "$file" > "$file".out
echo "Processed $file"
done


# extract the immediate feature types and count them
cd /lustre04/scratch/celphin/Dryas/snpEff/methylkit/
FILES=$(ls *.out)
echo ${FILES}

for file in ${FILES}; do
awk -F'[;:]' '{print $1 $2}' "$file" > "$file"1
echo 
echo "$file Upstream"
grep "Upstream" "$file"1 | wc -l 
echo "$file Downstream"
grep "Downstream" "$file"1 | wc -l 
echo "$file Gene"
grep "Gene" "$file"1 | wc -l 
echo "$file Intergenic"
grep "Intergenic" "$file"1 | wc -l 
echo "$file Intron"
grep "Intron" "$file"1 | wc -l 
echo "$file Exon"
grep "Exon" "$file"1 | wc -l 
done > methylkit_gene_features.txt

# Methylkit_Pheno_10_DMRs_CHH.bed.out Upstream
# 5145
# Methylkit_Pheno_10_DMRs_CHH.bed.out Downstream
# 4020
# Methylkit_Pheno_10_DMRs_CHH.bed.out Gene
# 167
# Methylkit_Pheno_10_DMRs_CHH.bed.out Intergenic
# 5055
# Methylkit_Pheno_10_DMRs_CHH.bed.out Intron
# 209
# Methylkit_Pheno_10_DMRs_CHH.bed.out Exon
# 967

# Methylkit_Pheno_25_DMRs_CHH.bed.out Upstream
# 1292
# Methylkit_Pheno_25_DMRs_CHH.bed.out Downstream
# 980
# Methylkit_Pheno_25_DMRs_CHH.bed.out Gene
# 39
# Methylkit_Pheno_25_DMRs_CHH.bed.out Intergenic
# 1402
# Methylkit_Pheno_25_DMRs_CHH.bed.out Intron
# 52
# Methylkit_Pheno_25_DMRs_CHH.bed.out Exon
# 197

#-----------------------------
# make into table automatically - not working

# module load StdEnv/2023
# module load r/4.4.0

# R

# # Load necessary libraries
# library(tidyr)
# library(dplyr)
# library(zoo)

# # Read the input data as a character vector (line by line)
# data_raw <- readLines("methylkit_gene_features.txt")

# # Initialize empty vectors to store the parsed data
# file_names <- c()
# categories <- c()
# values <- c()

# # Loop through the lines and extract the data
# for (i in seq(1, length(data_raw), by = 2)) {
  # # Read the file name (every odd line)
  # file_names <- c(file_names, data_raw[i])
  
  # # Read the corresponding category (every even line)
  # categories <- c(categories, data_raw[i + 1])
  
  # # Read the corresponding value (next odd line, ensuring it's numeric)
  # values <- c(values, as.numeric(data_raw[i + 2]))
# }

# # Combine into a data frame
# data <- data.frame(Methylkit_File = file_names, Category = categories, Value = values, stringsAsFactors = FALSE)

# # Reshape the data into wide format using pivot_wider
# data_wide <- data %>%
  # pivot_wider(names_from = Category, values_from = Value, values_fill = list(Value = 0))

# # Print the reshaped data to the console (optional)
# print(data_wide)

# # Save the reshaped data to a new file
# write.table(data_wide, "output_file.txt", sep = "\t", quote = FALSE, row.names = FALSE)



############################

# Look at transposons
# https://github.com/RobertsLab/project-gigas-oa-meth/blob/master/code/10-Genomic-Location-of-DML.ipynb 

module load StdEnv/2023 bedtools/2.31.0

cd /lustre04/scratch/celphin/Dryas/snpEff/methylkit/
FILES=$(ls *.bed)
echo ${FILES}

for file in ${FILES}; do
intersectBed \
-u \
-a ${file} \
-b /lustre04/scratch/celphin/Dryas/snpEff/data/OldDoct/OldDoct.genome.fa.mod.EDTA.TEanno.gff3 \
> ${file}-All-TE.bed

#head ${file}-All-TE.bed
wc -l ${file}-All-TE.bed

done

wc -l *-All-TE.bed


#-------------------------
# could subet the TEs by type and rerun for each type

awk '{print $3}' /lustre04/scratch/celphin/Dryas/snpEff/data/OldDoct/OldDoct.genome.fa.mod.EDTA.TEanno.gff3 | sort | uniq

# CACTA_TIR_transposon
# Copia_LTR_retrotransposon
# Gypsy_LTR_retrotransposon
# hAT_TIR_transposon
# helitron
# identity
# long_terminal_repeat
# LTR_retrotransposon
# Mutator_TIR_transposon
# Nov
# PIF_Harbinger_TIR_transposon
# repeat_region
# sequence_ontology
# site
# target_site_duplication
# Tc1_Mariner_TIR_transposon

cd /lustre04/scratch/celphin/Dryas/snpEff/methylkit/
# Extract unique values from the third column (or any other column)
unique_values=$(awk '{print $3}' /lustre04/scratch/celphin/Dryas/snpEff/data/OldDoct/OldDoct.genome.fa.mod.EDTA.TEanno.gff3 | sort | uniq)

# Loop through each unique value and extract lines to a new file
for value in $unique_values; do
    grep "$value" /lustre04/scratch/celphin/Dryas/snpEff/data/OldDoct/OldDoct.genome.fa.mod.EDTA.TEanno.gff3 > "${value}.gff3"
    echo "Extracted lines for $value into ${value}.gff3"
done

mkdir gff3
mv *.gff3 gff3/


#-----------------
module load StdEnv/2023 bedtools/2.31.0
unique_values=$(awk '{print $3}' /lustre04/scratch/celphin/Dryas/snpEff/data/OldDoct/OldDoct.genome.fa.mod.EDTA.TEanno.gff3 | sort | uniq)

cd /lustre04/scratch/celphin/Dryas/snpEff/methylkit
FILES=$(ls ./*-All-TE.bed)
echo ${FILES}

for value in $unique_values; do
for file in ${FILES}; do
intersectBed \
-u \
-a ${file} \
-b ./gff3/${value}.gff3 \
> ${file}-${value}.bed

wc -l ${file}-${value}.bed

done
done > DMR_TE_types.txt

#----------------------
# make into a table

module load StdEnv/2023
module load r/4.4.0

R

# Load necessary libraries
library(tidyr)
library(dplyr)
library(zoo)

# read in data
data_lines <- readLines("DMR_TE_types.txt")    # Read the lines from the file

# Step 2: Split the lines into a data frame (count and file path)
data <- do.call(rbind, strsplit(data_lines, " ", fixed = TRUE))
counts <- as.numeric(data[, 1])  # Extract the counts (first column)
files <- data[, 2]              # Extract the file paths (second column)

# Step 3: Define the methylkit and TE files
methylkit_files <- unique(sub("-All-TE.bed.*", "", files)) # Extract methylkit file names
te_files <- unique(sub(".*-All-TE.bed-", "", files))      # Extract TE file names

# Step 4: Create a matrix to store the counts
matrix_data <- matrix(0, nrow = length(methylkit_files), ncol = length(te_files),
                      dimnames = list(methylkit_files, te_files))

# Step 5: Populate the matrix with the counts
for (i in 1:length(counts)) {
  methylkit_file <- sub("-All-TE.bed.*", "", files[i])  # Extract methylkit file
  te_file <- sub(".*-All-TE.bed-", "", files[i])        # Extract TE file
  matrix_data[methylkit_file, te_file] <- counts[i]     # Assign the count to the matrix
}

# Step 6: Print the matrix (or save it to a file, if needed)
print(matrix_data)

write.table(matrix_data, "DMR_TE_types_table.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names=TRUE)


#############################
# Extract the gene names for each DMR
# try with just one file
cd /lustre04/scratch/celphin/Dryas/snpEff/methylkit/

awk 'BEGIN {OFS="\t"} {
    if ($0 !~ /^#/) {
        chromo = $1
        start = $2
        end = $3
        genes = $4
        # Use a regular expression to match only the gene names starting with "Do1_0" and ending with "V1.1"
        while (match(genes, /Do1_0[^\s]*V1.1/)) {
            gene_name = substr(genes, RSTART, RLENGTH)
            # Print chromo, start, end, and the gene name, trimming any unwanted characters
            print chromo, start, end, gene_name
            # Remove the matched gene name from the genes string and continue searching
            genes = substr(genes, RSTART + RLENGTH)
        }
    }
}'  Methylkit_WILL_rand_W_C_25plowcovDMRs.bed.out > Methylkit_WILL_rand_W_C_25plowcovDMRs.bed.out2

head Methylkit_WILL_rand_W_C_25plowcovDMRs.bed.out2
sort Methylkit_WILL_rand_W_C_25plowcovDMRs.bed.out2 | uniq > Methylkit_WILL_rand_W_C_25plowcovDMRs.bed.uniq_genes.out

head Methylkit_WILL_rand_W_C_25plowcovDMRs.bed.uniq_genes.out

# Do1_00174       48001   49000   Do1_00174G00007V1.1
# Do1_00174       48001   49000   Do1_00174G00008V1.1
# Do1_01_a00001   10352001        10353000        Do1_01_a00001G01901V1.1
# Do1_01_a00001   10352001        10353000        Do1_01_a00001G01902V1.1
# Do1_01_a00001   10947001        10948000        Do1_01_a00001G01999V1.1
# Do1_01_a00001   11095001        11096000        Do1_01_a00001G02021V1.1
# Do1_01_a00001   11095001        11096000        Do1_01_a00001G02022V1.1
# Do1_01_a00001   11095001        11096000        Do1_01_a00001G02023V1.1
# Do1_01_a00001   11095001        11096000        Do1_01_a00001G02024V1.1
# Do1_01_a00001   11318001        11319000        Do1_01_a00001G02056V1.1

wc -l Methylkit_WILL_rand_W_C_25plowcovDMRs.bed.out
#408 Methylkit_WILL_rand_W_C_25plowcovDMRs.bed.out

wc -l Methylkit_WILL_rand_W_C_25plowcovDMRs.bed.uniq_genes.out
#899 Methylkit_WILL_rand_W_C_25plowcovDMRs.bed.uniq_genes.out

#---------------------------
# extract also the information about position
# not working
# awk 'BEGIN {OFS="\t"} {
    # if ($0 !~ /^#/) {
        # chromo = $1
        # start = $2
        # end = $3
        # genes = $4

        # # Extract the gene name (the first part after Gene:), without the :protein_coding part
        # gene_name = gensub(/.*Gene:([^:;|]*)/, "\\1", "g", genes)

        # # Split the fields by semicolon
        # split(genes, fields, ";")

        # # Loop through the fields to find the relevant annotations
        # first_gene_found = 0
        # for (i in fields) {
            # if (fields[i] ~ /Gene:/ && first_gene_found == 0) {
                # first_gene_found = 1
                # print chromo, start, end, gene_name
            # }
            # # If it's an annotation (e.g., Transcript, Downstream, Exon, etc.)
            # if (fields[i] ~ /Transcript:/ || fields[i] ~ /Downstream:/ || fields[i] ~ /Exon:/ || fields[i] ~ /Upstream:/) {
                # print chromo, start, end, gene_name, fields[i]
            # }
        # }
    # }
# }'  Methylkit_WILL_rand_W_C_25plowcovDMRs.bed.out > Methylkit_WILL_rand_W_C_25plowcovDMRs.bed.out3

# head Methylkit_WILL_rand_W_C_25plowcovDMRs.bed.out3

#-----------------------
# to loop
# extract the immediate feature types and count them
cd /lustre04/scratch/celphin/Dryas/snpEff/methylkit/
FILES=$(ls *.out)
echo ${FILES}

FILES=$(ls Methylkit_Pheno_*_DMRs_CHH.bed.out)
echo ${FILES}


for file in ${FILES}; do

awk 'BEGIN {OFS="\t"} {
    if ($0 !~ /^#/) {
        chromo = $1
        start = $2
        end = $3
        genes = $4
        # Use a regular expression to match only the gene names starting with "Do1_0" and ending with "V1.1"
        while (match(genes, /Do1_0[^\s]*V1.1/)) {
            gene_name = substr(genes, RSTART, RLENGTH)
            # Print chromo, start, end, and the gene name, trimming any unwanted characters
            print chromo, start, end, gene_name
            # Remove the matched gene name from the genes string and continue searching
            genes = substr(genes, RSTART + RLENGTH)
        }
    }
}'  ${file} > ${file}2

sort ${file}2 | uniq > ${file}.uniq_genes

echo ${file} 

done

cp  /lustre04/scratch/celphin/Dryas/snpEff/methylkit/*.uniq_genes /lustre04/scratch/celphin/Dryas/methylkit_merged_data/

###################################
# join with Interproscan data and GO terms

# See Merging and then Gene ontology notes file
