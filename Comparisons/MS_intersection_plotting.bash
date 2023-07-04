####################################################
#Experimenting with left joining Heat map plots for intersection of RNA seq and Phenology DMRS 
#Step 1: Follow Cassandra's Notes3b to make heat map of WC - MatSen
####################################################
#Require:
    #interproscan_goterms_overlap_rna_subtract_W_C_Mat_Sen.tsv 
    #metilene_total_subtract_W_C_Mat_Sen.bedGraph

##################################################

tmux new-session -s DMR
tmux attach-session -t DMR

salloc -c1 --time 3:00:00 --mem 120000m --account def-rieseber

#-----------------------------------------------
cd scratch
cp ~/projects/def-rieseber/Dryas_shared_data/MS_blast_input_bedgraphs/total_subtract_W_C_Mat_Sen.bedGraph .
cp ../Wild_W_C_input_files/Wild_metilene_W_C.input .


infile="total_subtract_W_C_Mat_Sen.bedGraph"
wc -l $infile
#316
#Count columns 
awk -F'\t' '{print NF}' "Wild_metilene_W_C.input" | sort -nu | tail -n 1 
# 104 (102samples)

#Sets chrom to be first field, start and end to be 2nd and 3rd end and prints
#This loop creates a lot of txt files using the metilene file
for X in {1..316}
do
chrom=$(awk -v X=${X} 'NR == X {print $1}' $infile)
start=$(awk -v X=${X} 'NR == X {print $2}' $infile)
end=$(awk -v X=${X} 'NR == X {print $3}' $infile)
echo "${chrom}_${start}_${end}"
awk -v chrom="$chrom" -v start="$start" -v end="$end" -F "\t" '{ if(($1 == chrom) && ($2 <= end && $2 >= start)) { print } }' Wild_metilene_W_C.input  > "${chrom}_${start}_${end}.txt"

#For each DMR of interest (in subtracted domain)
#Creates by file values for DMRS
for i in {1..104}
do
  awk -v i="$i" -F "\t" '{Total=Total+$i} {print Total}' "${chrom}_${start}_${end}.txt" >  TmpLine
  tail -n 1 TmpLine > Tmpvalue
  cat Tmpvalue >> "Line_${chrom}_${start}_${end}.txt"
done

done

paste Line* > ${infile}_DMRs_summed.txt
ls Line* > ${infile}_DMRs_list.txt

mkdir Line_files
mv Line* Line_files/
#########################################################################################
#Exit compute node


module load StdEnv/2020
module load r/4.2.1
module load gdal
module load udunits
module load python
export R_LIBS_USER=/home/msandler/R/x86_64-pc-linux-gnu-library/4.2.1/

R

#install the package edgeR
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.15")

BiocManager::install("edgeR")
install.packages('gplots')

#paste it in here (i.e. replace my path with yours):
setwd ("/home/msandler/scratch/Wild_Metilene/Plotting")

#find your working directory:
getwd()

#load the libraries you will need 
library ("edgeR")
library ("gplots")
#library ("Equitable")
library(RColorBrewer)


filename="total_subtract_W_C_Mat_Sen.bedGraph_DMRs_summed.txt"
filename_list="total_subtract_W_C_Mat_Sen.bedGraph_DMRs_list.txt"

#read in the data
mydata <- read.table(filename, header=FALSE)
colnam <- read.table(filename_list, header=FALSE)

nrow(colnam)
#316
ncol(mydata)
#316
colnam =colnam[1:316,]
nrow(mydata)
#104

rownam <- c("Chrom", "Pos", "W_CASS_04W_519_NA",  "W_WILL_7W_448_107",  "W_FERT_14W_6F_126",  "W_LATJ_04W_NA_188",  "W_LATC_1W_12_219",  
            "W_ALAS_0W_14_249","W_LATC_9W_11_216",  "W_DRY_6W_31_147",  "W_SVAL_6W_6_274",  "W_SVAL_18W_18_271",  
            "W_LATD_4W_8_211",  "W_MEAD_2W_450_70",     "W_ALAS_00W_00_228",  "W_DRY_9W_50_185",  "W_LATD_2W_5_206",  "W_ALAS_0W_8_265", 
            "W_WILL_5W_421_154",  "W_LATJ_02W_00_193","W_ALAS_0W_7_263",  "W_MEAD_6W_466_163",  "W_WILL_4W_417_13",  "W_MEAD_1W_444_116", 
            "W_MEAD_08W_4R0_NA",  "W_LATC_3W_16_220",    "W_LATD_2W_4_212",  "W_WILL_1W_403_67",  "W_SVAL_16W_16_277",  "W_DRY_1W_3_39",  
            "W_DRY_3W_15_69",  "W_CASS_9W_539_128", "W_ALAS_0W_18_248",  "W_SVAL_0W_00_267",  "W_FERT_22W_12F_111",  "W_CASS_17W_574_137", 
            "W_CASS_10W_544_60",  "W_ALAS_0W_15_242","W_ALAS_0W_16_239",  "W_FERT_30W_14F_40",  "W_LATD_4W_9_207",  "W_MEAD_7W_470_173", 
            "W_ALAS_0W_17_235",  "W_SVAL_0W_NA_270", "W_CASS_7W_600_19",  "W_SVAL_8W_8_272",  "W_ALAS_00W_00_232",  "W_DRY_8W_45_155", 
            "W_FERT_6W_3F_110",  "W_LATC_5W_18_190",     "W_CASS_5W_525_130",   "W_ALAS_0W_3_236",  "W_WILL_10W_437_84", "C_WILL_5C_422_31","C_CASS_10C_548_144", 
			"C_CASS_4C_524_4",  "C_FERT_39C_20F_71",  "C_FERT_5C_1F_97",  "C_ALAS_00C_00_227",  "C_SVAL_16C_16_276",  
            "C_CASS_8C_535_54",  "C_ALAS_0C_10_246",    "C_DRY_9C_53_149",  "C_CASS_09C_00_541",  "C_MEAD_7C_473_95",  "C_MEAD_03C_1R0_00",  "C_MEAD_2C_451_76", 
            "C_ALAS_0C_13_254", "C_LATD_5C_2_201", "C_LATD_2C_7_209",  "C_WILL_7C_445_125",  "C_LATJ_00C_00_187",  "C_MEAD_1C_446_33",  "C_ALAS_0C_18_229", 
            "C_DRY_5C_28_92", "C_DRY_2C_10_52",  "C_LATD_2C_6_198",  "C_LATD_1C_4_223",  "C_SVAL_0C_00_268",  "C_SVAL_0C_00_269",  "C_WILL_10C_440_16", 
            "C_LATD_2C_1_203","C_LATD_5C_5_191",  "C_SVAL_12C_12_273",  "C_FERT_13C_7F_112",  "C_LATJ_02C_00_194",  "C_ALAS_0C_19_261", 
            "C_CASS_17C_576_175","C_DRY_10C_60_41",  "C_MEAD_6C_468_22",  "C_WILL_3C_414_100",  "C_SVAL_49C_49_278",  "C_FERT_31C_15F_170", 
            "C_ALAS_0C_00_231",  "C_SVAL_8C_8_275","C_ALAS_0C_5_238",  "C_ALAS_0C_12_256",  "C_LATD_4C_3_196",  "C_ALAS_0C_4_240", "C_CASS_5C_529_159", 
            "C_DRY_4C_23_82",  "C_LATD_5C_20_199",  "C_WILL_1C_406_152",  "C_ALAS_0C_3_258")

colnames(mydata) <- colnam
rownames(mydata) <- rownam

##############################
# total data plotting

#Keeps all but first two columns
mydata0 <-  mydata[-c(1,2),]

#Sort
mydata1 <- mydata0[order(row.names(mydata0)), ]

#transposes matrix
mydata2 <-  t(as.matrix(mydata1))

# https://r-graph-gallery.com/215-the-heatmap-function.html

jpeg(paste0(filename,"Total_Subtract_W_C_Phenology_DMR_heatmap_rainbow_order.jpg"), width = 3000, height = 1000)
heatmap.2 (mydata2, scale = "row", trace = "none", Colv = FALSE, labCol = colnames(mydata2), labRow = rownames(mydata2), col='rainbow', dendrogram = "none", colsep =c(51), margins =c(20,70), sepcolor = "white", sepwidth = 0.1, cexRow=2, cexCol=2)
#dev.off() (still do this it just breaks bash for some reaosn)


jpeg(paste0(filename,"Total_Subtract_W_C_Phenology_DMR_heatmap_red_order.jpg"), width = 3000, height = 1000)
heatmap.2 (mydata2, scale = "row", trace = "none", Colv = FALSE, labCol = colnames(mydata2), labRow = rownames(mydata2),  dendrogram = "none", colsep =c(51), margins =c(20,70), sepcolor = "white", sepwidth = 0.1, cexRow=2, cexCol=2)
#dev.off()

##############################

#On my local computer:
scp -v msandler@cedar.computecanada.ca:/home/msandler/scratch/Wild_Metilene/Plotting/*.jpg .

mkdir Subtract_W_C_Phenology_Plotting
mv * Subtract_W_C_Phenology_Plotting/

#############################################
#Plotting RNAseq DMRS:
#Step 1: get RNA DMR intersection file
#Step 2: Get list of genes
#Step 3: Get list of associated DMRS from blast file
#Step 4: Select those genes from bedgraph file to make subset
#Step 5: Plot them, similarly to this heatmap

###########################################
#Note this is where left join would be nice but haven't figured it out yet
#Messy work for pseudo left join
cp /home/msandler/projects/def-rieseber/msandler/Scratch_interproscan/interproscan_goterms_overlap_rna_subtract_W_C_Mat_Sen.tsv .
#Get list of gene ID's
awk  '!a[$1]++ {print $1}' interproscan_goterms_overlap_rna_subtract_W_C_Mat_Sen.tsv > gene_ids.txt

#copy reference blast here
cp /home/msandler/projects/def-rieseber/msandler/MS_blast_output/Blast_ref_output/blast_ref_total_subtract_W_C_Mat_Sen.out .

#Subset -> all intersection of first and 2nd column
#FNR iterates through second file
#NR through first file 
#records all first field of first file in id
#If 2nd field of 2nd file in id then prints
#Sets seperator to 
awk -F '\t' 'FNR==NR { ids[$1]; next } $2 in ids' gene_ids.txt blast_ref_total_subtract_W_C_Mat_Sen.out > subset.out 

#Removes any repeated sections (even if for different genes)
awk '!seen[$1]++' subset.out  > RNA_Warming_Minus_Phenology_blast.out

#Get intersection from metilene bedgraph DMRS with these DMRSut


#Cross over to the bedgraph

cp ~/projects/def-rieseber/Dryas_shared_data/MS_blast_input_bedgraphs/total_subtract_W_C_Mat_Sen.bedGraph .



#Maybe work sequentially:
    #Pipe all the awks
        #Matching first field
        #any interlap at all keep DMRS?? than 2nd field
        #

module load bedtools/2.30.0

#splits first field based on :,- as field seprators into the array a
awk '{split($1, a, /[:-]/); print a[1] "\t" a[2] "\t" a[3] "\t" $1}' RNA_Warming_Minus_Phenology_blast.out > Split_fields_genes.bedGraph

#WORKS 
bedtools intersect -wa -a total_subtract_W_C_Mat_Sen.bedGraph -b Split_fields_genes.bedGraph > RNA_total_subtract_W_C_Mat_Sen.bedGraph


RNA_total_subtract_W_C_Mat_Sen.bedGraph

#####################################################################
tmux new-session -s DMR
tmux attach-session -t DMR

salloc -c1 --time 3:00:00 --mem 120000m --account def-rieseber


infile="RNA_total_subtract_W_C_Mat_Sen.bedGraph"
wc -l $infile
#11
#Count columns 
cp ../Wild_W_C_input_files/Wild_metilene_W_C.input .
awk -F'\t' '{print NF}' "Wild_metilene_W_C.input" | sort -nu | tail -n 1 
# 104 (102samples)

#Sets chrom to be first field, start and end to be 2nd and 3rd end and prints
#This loop creates a lot of txt files using the metilene file
for X in {1..11}
do
chrom=$(awk -v X=${X} 'NR == X {print $1}' $infile)
start=$(awk -v X=${X} 'NR == X {print $2}' $infile)
end=$(awk -v X=${X} 'NR == X {print $3}' $infile)
echo "${chrom}_${start}_${end}"
awk -v chrom="$chrom" -v start="$start" -v end="$end" -F "\t" '{ if(($1 == chrom) && ($2 <= end && $2 >= start)) { print } }' Wild_metilene_W_C.input  > "${chrom}_${start}_${end}.txt"

#For each DMR of interest (in subtracted domain)
#Creates by file values for DMRS
for i in {1..104}
do
  awk -v i="$i" -F "\t" '{Total=Total+$i} {print Total}' "${chrom}_${start}_${end}.txt" >  TmpLine
  tail -n 1 TmpLine > Tmpvalue
  cat Tmpvalue >> "Line_${chrom}_${start}_${end}.txt"
done

done

paste Line* > ${infile}_DMRs_summed.txt
ls Line* > ${infile}_DMRs_list.txt

mkdir Line_files
mv Line* Line_files/

#########################################################################################
#Exit compute node


module load StdEnv/2020
module load r/4.2.1
module load gdal
module load udunits
module load python
export R_LIBS_USER=/home/msandler/R/x86_64-pc-linux-gnu-library/4.2.1/

R

#install the package edgeR
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.15")

BiocManager::install("edgeR")
install.packages('gplots')

#paste it in here (i.e. replace my path with yours):
setwd ("/home/msandler/scratch/Wild_Metilene/Plotting")

#find your working directory:
getwd()

#load the libraries you will need 
library ("edgeR")
library ("gplots")
#library ("Equitable")
library(RColorBrewer)


filename="RNA_total_subtract_W_C_Mat_Sen.bedGraph_DMRs_summed.txt"
filename_list="RNA_total_subtract_W_C_Mat_Sen.bedGraph_DMRs_list.txt"

#read in the data
mydata <- read.table(filename, header=FALSE)
colnam <- read.table(filename_list, header=FALSE)

nrow(colnam)
#11
ncol(mydata)
#11
colnam =colnam[1:11,]
nrow(mydata)
#104

rownam <- c("Chrom", "Pos", "W_CASS_04W_519_NA",  "W_WILL_7W_448_107",  "W_FERT_14W_6F_126",  "W_LATJ_04W_NA_188",  "W_LATC_1W_12_219",  
            "W_ALAS_0W_14_249","W_LATC_9W_11_216",  "W_DRY_6W_31_147",  "W_SVAL_6W_6_274",  "W_SVAL_18W_18_271",  
            "W_LATD_4W_8_211",  "W_MEAD_2W_450_70",     "W_ALAS_00W_00_228",  "W_DRY_9W_50_185",  "W_LATD_2W_5_206",  "W_ALAS_0W_8_265", 
            "W_WILL_5W_421_154",  "W_LATJ_02W_00_193","W_ALAS_0W_7_263",  "W_MEAD_6W_466_163",  "W_WILL_4W_417_13",  "W_MEAD_1W_444_116", 
            "W_MEAD_08W_4R0_NA",  "W_LATC_3W_16_220",    "W_LATD_2W_4_212",  "W_WILL_1W_403_67",  "W_SVAL_16W_16_277",  "W_DRY_1W_3_39",  
            "W_DRY_3W_15_69",  "W_CASS_9W_539_128", "W_ALAS_0W_18_248",  "W_SVAL_0W_00_267",  "W_FERT_22W_12F_111",  "W_CASS_17W_574_137", 
            "W_CASS_10W_544_60",  "W_ALAS_0W_15_242","W_ALAS_0W_16_239",  "W_FERT_30W_14F_40",  "W_LATD_4W_9_207",  "W_MEAD_7W_470_173", 
            "W_ALAS_0W_17_235",  "W_SVAL_0W_NA_270", "W_CASS_7W_600_19",  "W_SVAL_8W_8_272",  "W_ALAS_00W_00_232",  "W_DRY_8W_45_155", 
            "W_FERT_6W_3F_110",  "W_LATC_5W_18_190",     "W_CASS_5W_525_130",   "W_ALAS_0W_3_236",  "W_WILL_10W_437_84", "C_WILL_5C_422_31","C_CASS_10C_548_144", 
			"C_CASS_4C_524_4",  "C_FERT_39C_20F_71",  "C_FERT_5C_1F_97",  "C_ALAS_00C_00_227",  "C_SVAL_16C_16_276",  
            "C_CASS_8C_535_54",  "C_ALAS_0C_10_246",    "C_DRY_9C_53_149",  "C_CASS_09C_00_541",  "C_MEAD_7C_473_95",  "C_MEAD_03C_1R0_00",  "C_MEAD_2C_451_76", 
            "C_ALAS_0C_13_254", "C_LATD_5C_2_201", "C_LATD_2C_7_209",  "C_WILL_7C_445_125",  "C_LATJ_00C_00_187",  "C_MEAD_1C_446_33",  "C_ALAS_0C_18_229", 
            "C_DRY_5C_28_92", "C_DRY_2C_10_52",  "C_LATD_2C_6_198",  "C_LATD_1C_4_223",  "C_SVAL_0C_00_268",  "C_SVAL_0C_00_269",  "C_WILL_10C_440_16", 
            "C_LATD_2C_1_203","C_LATD_5C_5_191",  "C_SVAL_12C_12_273",  "C_FERT_13C_7F_112",  "C_LATJ_02C_00_194",  "C_ALAS_0C_19_261", 
            "C_CASS_17C_576_175","C_DRY_10C_60_41",  "C_MEAD_6C_468_22",  "C_WILL_3C_414_100",  "C_SVAL_49C_49_278",  "C_FERT_31C_15F_170", 
            "C_ALAS_0C_00_231",  "C_SVAL_8C_8_275","C_ALAS_0C_5_238",  "C_ALAS_0C_12_256",  "C_LATD_4C_3_196",  "C_ALAS_0C_4_240", "C_CASS_5C_529_159", 
            "C_DRY_4C_23_82",  "C_LATD_5C_20_199",  "C_WILL_1C_406_152",  "C_ALAS_0C_3_258")

colnames(mydata) <- colnam
rownames(mydata) <- rownam

##############################
# total data plotting

#Keeps all but first two columns
mydata0 <-  mydata[-c(1,2),]

#Sort
mydata1 <- mydata0[order(row.names(mydata0)), ]

#transposes matrix
mydata2 <-  t(as.matrix(mydata1))

# https://r-graph-gallery.com/215-the-heatmap-function.html

jpeg(paste0(filename,"Subtract_W_C_Phenology_DMR_heatmap_rainbow_order.jpg"), width = 3000, height = 1000)
heatmap.2 (mydata2, scale = "row", trace = "none", Colv = FALSE, labCol = colnames(mydata2), labRow = rownames(mydata2), col='rainbow', dendrogram = "none", colsep =c(51), margins =c(20,70), sepcolor = "white", sepwidth = 0.1, cexRow=2, cexCol=2)
#dev.off() #(still do this it just breaks bash for some reaosn)


jpeg(paste0(filename,"Subtract_W_C_Phenology_DMR_heatmap_red_order.jpg"), width = 3000, height = 1000)
heatmap.2 (mydata2, scale = "row", trace = "none", Colv = FALSE, labCol = colnames(mydata2), labRow = rownames(mydata2),  dendrogram = "none", colsep =c(51), margins =c(20,70), sepcolor = "white", sepwidth = 0.1, cexRow=2, cexCol=2)
#dev.off()

mkdir Methylation_Heatmap_W_C_Minus_Phenology








#######################
#Next step:
#Make RNA expression heatmap for 8 genes of interest (original)

#1 get interproscanRNA intersect file
#Get list of genes
#Intersect back to RNA seq file
#Plot heat map of RNA seq


cp /home/msandler/projects/def-rieseber/Dryas_shared_data/CE_RNAseq_DERs/RNA_DER_May2023_W_C_Total.txt .
cp /home/msandler/projects/def-rieseber/Dryas_shared_data/CE_RNAseq_raw_data/gene_names_expression_table1.txt .
cp /home/msandler/projects/def-rieseber/msandler/Scratch_interproscan/interproscan_goterms_overlap_rna_subtract_W_C_Mat_Sen.tsv .

awk '!a[$1]++ {print $1}' interproscan_goterms_overlap_rna_subtract_W_C_Mat_Sen.tsv > gene_ids.txt

awk 'FNR==NR{values[$0]; next} $1 in values' gene_ids.txt gene_names_expression_table1.txt > rna_w_c_gene_expressions.txt


#13 genes
#Need gene_names expressions table
#Add this to beggining of rna_w_c_gene_expression.txt

sed -i '1s/^/sequence C_Alas_16 C_Alas_1 C_Alas_20 C_Alas_5 C_Alex_11 C_Alex_21 C_Alex_22 C_Alex_23 C_Alex_24 C_Alex_2 C_
Alex_6 C_Alex_7 C_Norw_10 C_Norw_1 C_Norw_3_B C_Norw_5 C_Norw_6 C_Norw_7 C_Seed_1 C_Seed_2 C_Seed_3 C_Swed_1
1 C_Swed_15 C_Swed_16 C_Swed_2 C_Swed_3 C_Swed_5 C_Swed_8 W_Alas_12 W_Alas_13 W_Alas_15 W_Alas_1 W_Alas_2 W_
Alas_3 W_Alas_5 W_Alas_9 W_Alex_11 W_Alex_1 W_Alex_21 W_Alex_22 W_Alex_6 W_Norw_10_B W_Norw_2 W_Norw_3_B W_N
orw_4_B W_Norw_6 W_Norw_7 W_Seed_1 W_Seed_2 W_Seed_3 W_Swed_10 W_Swed_13 W_Swed_14 W_Swed_16 W_Swed_17 W_Swe
d_2 W_Swed_5 W_Swed_6\n/' file.txt






module load StdEnv/2020
module load r/4.2.1
module load gdal
module load udunits
module load python
export R_LIBS_USER=/home/msandler/R/x86_64-pc-linux-gnu-library/4.2.1/



R


#install the package edgeR
# if (!require("BiocManager", quietly = TRUE))
    # install.packages("BiocManager")
# BiocManager::install(version = "3.15")

# BiocManager::install("edgeR")
# install.packages('gplots')

#paste it in here (i.e. replace my path with yours):
setwd ("/scratch/msandler/Wild_Metilene/Plotting")

#find your working directory:
getwd()

#load the libraries you will need 
library ("edgeR")
library ("gplots")

#read in the data
mydata <- read.table("./rna_w_c_gene_expressions.txt", header=TRUE)

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
jpeg("./W_C_Min_Mat_Sen_RNA_MDS.jpg", width = 700, height = 500)
plotMDS (list2, col = col_treat)
dev.off()

jpeg("./Site_RNA_W_C_Min_Mat_MDS.jpg", width = 700, height = 500)
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
top <- topTags (lrt.list2, n = 2000)$table # p-value all < 5e-02 

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

jpeg("./W_C_min_Phen_RNA_heatmap.jpg", width = 3200, height = 1000)
heatmap.2 (heatmap2, scale = "row", trace = "none", Colv = FALSE, labCol = colnames(heatmap2), labRow = rownames(heatmap2), col='rainbow', dendrogram = "none", colsep =c(5), margins =c(20,70), sepcolor = "white", sepwidth = 0.1, cexRow=2, cexCol=2)
dev.off()

jpeg("./W_C_min_Phen_RNA_heatmap_legend.jpg", width = 1000, height = 800)
heatmap.2 (heatmap2, scale = "row", trace = "none", notecex=3, notecol="black", tracecol="black", Colv = FALSE, labCol = colnames(heatmap2), labRow = rownames(heatmap2), col='rainbow', dendrogram = "none", colsep =c(5), margins =c(0.2, 0.2), sepcolor = "white", sepwidth = 0.1, cexRow=6, cexCol=6)
dev.off()

# write out file of genes that differ a lot

write.table(heatmap2, file = "./RNA_DER_May2023_W_C_Min_Phen_Total.txt", quote = FALSE, row.names=TRUE, col.names=TRUE, sep="\t")


