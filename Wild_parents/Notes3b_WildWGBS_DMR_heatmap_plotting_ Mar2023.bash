###########################
# Plotting all the Dryas DMRs at each site - and overall?
# March 2023
############################

# Goals:
# One heat map per site of all the DMRs
# One total heat map of DMRs

#########################
cd /home/celphin/projects/rpp-rieseber/celphin/Dryas/DMRs/June2022_Metilene_DMRs/data/

# upzip all the bedgraph files
# gunzip -v /home/celphin/projects/rpp-rieseber/celphin/Dryas/DMRs/June2022_Metilene_DMRs/data/*.gz


#########################
#test loop
cd /home/celphin/projects/rpp-rieseber/celphin/Dryas/DMRs/June2022_Metilene_DMRs/

more ./All_sites/ALL_Sites_intersect_DMRs_qval.0.05.txt
Do1_01_a00004   10489852        10490314        -9.776923

chrom="Do1_01_a00004"
start=10489852
end=10490314

DMR="Do1_01_a00004   10489852        10490314        -9.776923"
chrom=$(echo $line| awk -F\ '{print $1}' ALL_Sites_intersect_DMRs_qval.0.05.txt)
start=$(echo $line| awk -F\ '{print $2}' ${DMR})
end=$(echo $line| awk -F\ '{print $3}' ${DMR})

# https://www.tim-dennis.com/data/tech/2016/08/09/using-awk-filter-rows.html
# https://linuxhint.com/awk_command_variables/
# filter only rows in range of DMR
awk -v chrom="$chrom" -v start="$start" -v end="$end" -F "\t" '{ if(($1 == chrom) && ($2 <= end && $2 >= start)) { print } }' ./data/metilene_W_C.input  > "${chrom}_${start}_${end}.txt"

# counts columns
awk -F'\t' '{print NF}' "${chrom}_${start}_${end}.txt" | sort -nu | tail -n 1 
# 104 columns
# 102 samples

# https://www.baeldung.com/linux/add-column-of-numbers
# https://www.cyberciti.biz/faq/unix-linux-iterate-over-a-variable-range-of-numbers-in-bash/
# Sums each column from 3 to 104 for each indiv 
for i in {1..104}
do
  awk -v i="$i" -F "\t" '{Total=Total+$i} {print Total}' "${chrom}_${start}_${end}.txt" >  TmpLine
  tail -n 1 TmpLine > Tmpvalue
  cat Tmpvalue >> "Line_${chrom}_${start}_${end}.txt"
done

######################################
# run loop
tmux new-session -s DMR
tmux attach-session -t DMR

salloc -c1 --time 3:00:00 --mem 120000m --account def-rieseber

# loop through filtered metiliene output for each site

cd /home/celphin/projects/rpp-rieseber/celphin/Dryas/DMRs/June2022_Metilene_DMRs/

#infile="./All_sites/ALL_Sites_intersect_DMRs_qval.0.05.txt"
#infile="./TOTAL_DMRs/Total_DMRs_100-10-10_July2022_qval.1e-20.out"  # note cut short at just 8 genes
infile="./LAT_site/LAT_metilene_W_C_Sept2022_70_5_4_5_qval.1e-5.out"
# should try:
#infile="./LAT_site/LAT_DMRs_300-30-15_qval.1e-10.out"

# for X in wc -l $infile
wc -l $infile

for X in {1..212}
do
chrom=$(awk -v X=${X} 'NR == X {print $1}' $infile)
start=$(awk -v X=${X} 'NR == X {print $2}' $infile)
end=$(awk -v X=${X} 'NR == X {print $3}' $infile)
echo "${chrom}_${start}_${end}"

# https://www.tim-dennis.com/data/tech/2016/08/09/using-awk-filter-rows.html
# https://linuxhint.com/awk_command_variables/
# filter only rows in range of DMR
awk -v chrom="$chrom" -v start="$start" -v end="$end" -F "\t" '{ if(($1 == chrom) && ($2 <= end && $2 >= start)) { print } }' ./data/metilene_W_C.input  > "${chrom}_${start}_${end}.txt"

# counts columns # should be 104 always - or number of samples
# awk -F'\t' '{print NF}' "${chrom}_${start}_${end}.txt" | sort -nu | tail -n 1 
# 104 columns
# 102 samples

# https://www.baeldung.com/linux/add-column-of-numbers
# https://www.cyberciti.biz/faq/unix-linux-iterate-over-a-variable-range-of-numbers-in-bash/
# Sums each column from 3 to 104 for each indiv 
for i in {1..104}
do
  awk -v i="$i" -F "\t" '{Total=Total+$i} {print Total}' "${chrom}_${start}_${end}.txt" >  TmpLine
  tail -n 1 TmpLine > Tmpvalue
  cat Tmpvalue >> "Line_${chrom}_${start}_${end}.txt"
done

done

paste Line* > ${infile}_DMRs_summed.txt
ls Line* > ${infile}_DMRs_list.txt

mv Line* Line_files/


###########################################################
# try differential methylation analysis - plotting
# and EdgeR plotting for DMRs

# https://github.com/owensgl/biol525D/tree/master/Topic_6
# https://bioconductor.org/packages/release/bioc/html/edgeR.html
# http://combine-australia.github.io/RNAseq-R/slides/RNASeq_filtering_qc.pdf 

cd /home/celphin/projects/rpp-rieseber/celphin/Dryas/DMRs/June2022_Metilene_DMRs

# get R on server
# open R
module load StdEnv/2020
module load r/4.2.1
module load gdal
module load udunits
module load python
export R_LIBS_USER=/home/celphin/R/x86_64-pc-linux-gnu-library/4.2.1/

R

#--------------------------
# EdgeR
# DMR plotting

#install the package edgeR
# if (!require("BiocManager", quietly = TRUE))
    # install.packages("BiocManager")
# BiocManager::install(version = "3.15")

# BiocManager::install("edgeR")
# install.packages('gplots')

#paste it in here (i.e. replace my path with yours):
setwd ("/home/celphin/projects/rpp-rieseber/celphin/Dryas/DMRs/June2022_Metilene_DMRs/")

#find your working directory:
getwd()

#load the libraries you will need 
library ("edgeR")
library ("gplots")
#library ("Equitable")
library(RColorBrewer)

#filename="./All_sites/ALL_Sites_intersect_DMRs_qval.0.05.txt_DMRs_summed.txt"
#filename_list="./All_sites/ALL_Sites_intersect_DMRs_qval.0.05.txt_DMRs_list.txt"

filename="./LAT_site/LAT_metilene_W_C_Sept2022_70_5_4_5_qval.1e-5.out_DMRs_summed.txt"
filename_list="./LAT_site/LAT_metilene_W_C_Sept2022_70_5_4_5_qval.1e-5.out_DMRs_list.txt"

#read in the data
mydata <- read.table(filename, header=FALSE)
colnam <- read.table(filename_list, header=FALSE)

nrow(colnam)
ncol(mydata)
colnam =colnam[1:212,]
nrow(mydata)

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

mydata0 <-  mydata[-c(1,2),]

mydata1 <- mydata0[order(row.names(mydata0)), ]

mydata2 <-  t(as.matrix(mydata1))

# https://r-graph-gallery.com/215-the-heatmap-function.html


jpeg(paste0(filename,"_Total_W_C_DMR_heatmap_rainbow_order.jpg"), width = 3000, height = 1000)
heatmap.2 (mydata2, scale = "row", trace = "none", Colv = FALSE, labCol = colnames(mydata2), labRow = rownames(mydata2), col='rainbow', dendrogram = "none", colsep =c(51), margins =c(20,70), sepcolor = "white", sepwidth = 0.1, cexRow=2, cexCol=2)
dev.off()

jpeg(paste0(filename,"_Total_W_C_DMR_heatmap_red_order.jpg"), width = 3000, height = 1000)
heatmap.2 (mydata2, scale = "row", trace = "none", Colv = FALSE, labCol = colnames(mydata2), labRow = rownames(mydata2),  dendrogram = "none", colsep =c(51), margins =c(20,70), sepcolor = "white", sepwidth = 0.1, cexRow=2, cexCol=2)
dev.off()

##############################
# Plot heatmap of site specific 

mydatas <- mydata[order(row.names(mydata)), ]

mydatas$ID <- rownames(mydatas)

q <- as.data.frame(t(as.matrix(as.data.frame((strsplit(as.character(mydatas$ID), "_"))))))
mydata <- cbind(mydatas, q)

colnames(mydata)[which(colnames(mydata)=="V1")] <- "Treat" 
colnames(mydata)[which(colnames(mydata)=="V2")] <- "SubSite"
colnames(mydata)[which(colnames(mydata)=="V3")] <- "Plot"
colnames(mydata)[which(colnames(mydata)=="V4")] <- "FieldID"
colnames(mydata)[which(colnames(mydata)=="V5")] <- "PlantID"

#makes sites
unique(mydata$SubSite)

mydata$Site <- sub('LATD', "LATJ", mydata$SubSite, ignore.case = FALSE)
mydata$Site <- sub('LATC', "LATJ", mydata$Site, ignore.case = FALSE)

mydata$Site <- sub('WILL', "Alex", mydata$Site, ignore.case = FALSE)
mydata$Site <- sub('CASS', "Alex", mydata$Site, ignore.case = FALSE)
mydata$Site <- sub('DRY', "Alex", mydata$Site, ignore.case = FALSE)
mydata$Site <- sub('FERT', "Alex", mydata$Site, ignore.case = FALSE)
mydata$Site <- sub('MEAD', "Alex", mydata$Site, ignore.case = FALSE)

unique(mydata$Site)

mydata <- mydata[-which(mydata$Treat=="Chrom"|mydata$Treat=="Pos"),]

colnames(mydata)

library(plyr)
library(dplyr)
library(tidyr)

mydata$Site <- as.factor(mydata$Site)

#------
# Sweden
LATJ_data  <- filter(mydata, Site=="LATJ")
# remove columns
LATJ_data1 <- subset(LATJ_data, select = -c(ID, Treat, SubSite, Plot, FieldID, PlantID, Site))
XX <- LATJ_data1

XX <-  as.matrix(XX)

jpeg(paste0(filename,"_LATJ_W_C_DMR_heatmap_order.jpg"), width = 6000, height = 2000)
heatmap.2 (XX, scale = "column", trace = "none", Rowv=NA, Colv = TRUE, labCol = colnames(XX), labRow = rownames(XX), col='rainbow', dendrogram = "none", rowsep =c(10), margins =c(20,70), sepcolor = "white", sepwidth = 0.1, cexRow=2, cexCol=2)
dev.off()

#------
# Alex
Alex_data  <- filter(mydata, Site=="Alex")
# remove columns
Alex_data1 <- subset(Alex_data, select = -c(ID, Treat, SubSite, Plot, FieldID, PlantID, Site))
XX <- Alex_data1

XX <-  t(as.matrix(XX))

jpeg(paste0(filename,"_Alex_W_C_DMR_heatmap_red_order.jpg"), width = 3000, height = 1000)
heatmap.2 (XX, scale = "row", trace = "none", Colv = FALSE, labCol = colnames(XX), labRow = rownames(XX),  dendrogram = "none", colsep =c(25), margins =c(20,70), sepcolor = "white", sepwidth = 0.1, cexRow=2, cexCol=2)
dev.off()

#------
# Alaska
ALAS_data  <- filter(mydata, Site=="ALAS")
# remove columns
ALAS_data1 <- subset(ALAS_data, select = -c(ID, Treat, SubSite, Plot, FieldID, PlantID, Site))
XX <- ALAS_data1

XX <-  t(as.matrix(XX))

jpeg(paste0(filename,"_ALAS_W_C_DMR_heatmap_red_order.jpg"), width = 3000, height = 1000)
heatmap.2 (XX, scale = "row", trace = "none", Colv = FALSE, labCol = colnames(XX), labRow = rownames(XX),  dendrogram = "none", colsep =c(10), margins =c(20,70), sepcolor = "white", sepwidth = 0.1, cexRow=2, cexCol=2)
dev.off()

#------
# Svalbard
SVAL_data  <- filter(mydata, Site=="SVAL")
# remove columns
SVAL_data1 <- subset(SVAL_data, select = -c(ID, Treat, SubSite, Plot, FieldID, PlantID, Site))
XX <- SVAL_data1

XX <-  t(as.matrix(XX))

jpeg(paste0(filename,"_SVAL_W_C_DMR_heatmap_red_order.jpg"), width = 3000, height = 1000)
heatmap.2 (XX, scale = "row", trace = "none", Colv = FALSE, labCol = colnames(XX), labRow = rownames(XX),  dendrogram = "none", colsep =c(10), margins =c(20,70), sepcolor = "white", sepwidth = 0.1, cexRow=2, cexCol=2)
dev.off()




################
# order genes by expression amount


#-----------------------
# plot barplots of total results

# take from wide to long

gather(mydata2)

# for each gene and individual - amount of meth - split by site/ treatment

jpeg(paste0("DMRs_qval.0.05_W_C_DMR_total_barplot.jpg"), width = 3000, height = 1000)

dev.off()



##################################
# plot specific DMR in bar plot

#install.packages("ggplot2")
#install.packages("ggpubr")

R

#libraries
library(ggplot2)
library(ggpubr)

setwd ("/home/celphin/projects/rpp-rieseber/celphin/Dryas/DMRs/June2022_Metilene_DMRs/")

#import data
#filename = "./Do1_01_a00004_10489852_10490314.txt" # DONE
#filename = "./Do1_05_a00003_1480241_1480901.txt" 
#filename = "./Do1_05_a00003_1480550_1480822.txt" 

filename = "./Do1_05_a00003_1152846_1152987.txt" # this is the amidase

#														chr	         start	end		q-value	        meanmethyl		CpG										meanW	meanC   
# LAT_metilene_W_C_Sept2022_70_5_4_5_qval.1e-5.out:  Do1_05_a00003  1152846 1152987 1.1428e-25      16.729108       24     								    63.045  46.316
# metilene_W_C_Sept2022_150_5_4:                     Do1_05_a00003  1152846 1152987 1.3592e-25      16.729108       24      3.9524e-14      9.0799e-32      63.045  46.316

# need to fix filenames and run
#filename = "./Do1_04_a00001_11637435_11638165.txt"
#filename = "./Do1_05_a00002_2523208_2523547.txt"


data <- read.table(filename, sep="\t", dec=".", header = TRUE)

# replace "W_(....)([0-9]{1,2}W_) with "W_\1_\2

# from metilene_W_C.input
colnames(data) <- c("Chrom", "Pos", "W_CASS_04W_519_NA",  "W_WILL_7W_448_107",  "W_FERT_14W_6F_126",  "W_LATJ_04W_NA_188",  "W_LATC_1W_12_219",  
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

################################	
			
datat <-  as.data.frame(t(as.matrix(data)))

# section to remove position specific information and sum for whole DMR for each indiv
datat1 <- datat[-c(1,2),]

#factors to numeric
X <- length(datat1)
factorToNumeric <- function(f) as.numeric(as.character(f))
datat1[1:X] <- lapply(datat[1:X], factorToNumeric)

# sum rows
datat1$total <- rowSums(datat1[1:X], na.rm=TRUE)
datat1$ID <- rownames(datat1)

mydata1 <- as.data.frame(cbind(datat1$ID ,datat1$total))
colnames(mydata1) <- c("ID", "TotalMeth")

# split ID

q <- as.data.frame(t(as.matrix(as.data.frame((strsplit(as.character(mydata1$ID), "_"))))))
mydata1 <- cbind(mydata1, q)

colnames(mydata1) <- c("ID", "TotalMeth","Treat", "SubSite", "Plot", "Plant")

mydata1$TotalMeth <- as.numeric(as.character(mydata1$TotalMeth))

#makes sites
unique(mydata1$SubSite)

mydata1$Site <- sub('LATD', "LATJ", mydata1$SubSite, ignore.case = FALSE)
mydata1$Site <- sub('LATC', "LATJ", mydata1$Site, ignore.case = FALSE)

mydata1$Site <- sub('WILL', "Alex", mydata1$Site, ignore.case = FALSE)
mydata1$Site <- sub('CASS', "Alex", mydata1$Site, ignore.case = FALSE)
mydata1$Site <- sub('DRY', "Alex", mydata1$Site, ignore.case = FALSE)
mydata1$Site <- sub('FERT', "Alex", mydata1$Site, ignore.case = FALSE)
mydata1$Site <- sub('MEAD', "Alex", mydata1$Site, ignore.case = FALSE)

mydata1$TotalMeth[which(mydata1$TotalMeth> (15000))] <- NA

jpeg(paste0(filename, "boxplot_subsite_warm_methylation.jpg"), width = 1000, height = 707)
ggboxplot(mydata1, x = "SubSite", y = "TotalMeth",  color = "Treat", add = "jitter", shape = "Treat")
dev.off()

jpeg(paste0(filename, "boxplot_site_warm_methylation.jpg"), width = 1000, height = 707)
ggboxplot(mydata1, x = "Site", y = "TotalMeth",  color = "Treat", add = "jitter", shape = "Treat")
dev.off()

write.csv(mydata1, file = paste0(filename,"DMR.csv"), quote = FALSE, row.names=FALSE)

jpeg(paste0(filename, "splitboxplot_site_warm_methylation.jpg"), width = 1000, height = 707)
ggplot(data=mydata1, aes(x=Treat, y=TotalMeth))+
  geom_boxplot(aes(x=Treat, y=TotalMeth, fill=Treat))+
  facet_wrap(~Site, scales = "free")+ theme_classic()+
  scale_fill_manual(values=c( "#89C5DA", "#DA5724"))+
  labs(fill="Treatment")
dev.off()

#######################################
# try making averaged bedGraph of each region W and C for a given DMR

# split ID

data$DMR_Pos <- paste0(data$Chrom, "_", data$Pos)
datat<-as.data.frame(t(as.matrix(data)))

colnames(datat) <- datat[105,]
datat1 <- datat[-c(1,2,105),]

Plant_IDs <- rownames(datat1)
q <- as.data.frame(t(as.matrix(as.data.frame((strsplit(Plant_IDs, "_"))))))
mydata <- cbind(datat1, q)

colnames(mydata)[which(colnames(mydata)=="V1")] <- "Treat" 
colnames(mydata)[which(colnames(mydata)=="V2")] <- "SubSite"
colnames(mydata)[which(colnames(mydata)=="V3")] <- "Plot"
colnames(mydata)[which(colnames(mydata)=="V4")] <- "FieldID"
colnames(mydata)[which(colnames(mydata)=="V5")] <- "PlantID"

#makes sites
unique(mydata$SubSite)

mydata$Site <- sub('LATD', "LATJ", mydata$SubSite, ignore.case = FALSE)
mydata$Site <- sub('LATC', "LATJ", mydata$Site, ignore.case = FALSE)

mydata$Site <- sub('WILL', "Alex", mydata$Site, ignore.case = FALSE)
mydata$Site <- sub('CASS', "Alex", mydata$Site, ignore.case = FALSE)
mydata$Site <- sub('DRY', "Alex", mydata$Site, ignore.case = FALSE)
mydata$Site <- sub('FERT', "Alex", mydata$Site, ignore.case = FALSE)
mydata$Site <- sub('MEAD', "Alex", mydata$Site, ignore.case = FALSE)

unique(mydata$Site)

# group by Site and Treatment
library(plyr)
library(dplyr)
library(tidyr)

colnames(mydata)
mydata1 <- gather(mydata, Pos, Methylation, Do1_05_a00003_1152848:Do1_05_a00003_1152987)

mydata1$Methylation <- as.numeric(mydata1$Methylation )

# average each column of the DMR
DMR_avg <- mydata1  %>%
  group_by(Site, Treat, Pos) %>%
  dplyr::summarise(avgMeth = mean(Methylation, na.rm=TRUE)) 

# split Pos
Pos_split <- as.data.frame(t(as.matrix(as.data.frame((strsplit(DMR_avg$Pos, "_"))))))

Start <- as.numeric(Pos_split[,4])
One <- as.numeric(rep(1, length(Start)))
End <-  Start + One

Chrom <- paste0(Pos_split[,1], "_", Pos_split[,2], "_", Pos_split[,3])

DMR_avg1 <- cbind(Chrom, Start, End, DMR_avg)

colnames(DMR_avg1) <- c("Chrom",    "Start",    "End",    "Site",    "Treat",   "Pos",     "avgMeth")

#DMR_avg2 <- DMR_avg1[-which(is.na(DMR_avg1$avgMeth)),]


# filter by Site and Treat
# total
tmpdata <- filter(DMR_avg1, Treat=="C")
Total_C_data <- tmpdata[,c("Chrom", "Start", "End","avgMeth")]
write.table(Total_C_data, file =  paste0(filename, "_total_C.bedGraph"), sep="\t", quote = FALSE, row.names=FALSE, col.names=FALSE)

tmpdata <- filter(DMR_avg1, Treat=="W")
Total_W_data <- tmpdata[,c("Chrom", "Start", "End","avgMeth")]
write.table(Total_W_data, file =  paste0(filename, "_total_W.bedGraph"), sep="\t", quote = FALSE, row.names=FALSE, col.names=FALSE)

# Alex
tmpdata <- filter(DMR_avg1, Site=="Alex" & Treat=="C")
Alex_C_data <- tmpdata[,c("Chrom", "Start", "End","avgMeth")]
write.table(Alex_C_data, file =  paste0(filename, "_Alex_C.bedGraph"), sep="\t",quote = FALSE, row.names=FALSE, col.names=FALSE)

tmpdata <- filter(DMR_avg1, Site=="Alex" & Treat=="W")
Alex_W_data <- tmpdata[,c("Chrom", "Start", "End","avgMeth")]
write.table(Alex_W_data, file =  paste0(filename, "_Alex_W.bedGraph"), sep="\t", quote = FALSE, row.names=FALSE, col.names=FALSE)

#SVAL
tmpdata <- filter(DMR_avg1, Site=="SVAL" & Treat=="C")
SVAL_C_data <- tmpdata[,c("Chrom", "Start", "End","avgMeth")]
write.table(SVAL_C_data, file =  paste0(filename, "_SVAL_C.bedGraph"),sep="\t", quote = FALSE, row.names=FALSE, col.names=FALSE)

tmpdata <- filter(DMR_avg1, Site=="SVAL" & Treat=="W")
SVAL_W_data <- tmpdata[,c("Chrom", "Start", "End","avgMeth")]
write.table(SVAL_W_data, file =  paste0(filename, "_SVAL_W.bedGraph"), sep="\t", quote = FALSE, row.names=FALSE, col.names=FALSE)

#ALAS
tmpdata <- filter(DMR_avg1, Site=="ALAS" & Treat=="C")
ALAS_C_data <- tmpdata[,c("Chrom", "Start", "End","avgMeth")]
write.table(ALAS_C_data, file =  paste0(filename, "_ALAS_C.bedGraph"), sep="\t", quote = FALSE, row.names=FALSE, col.names=FALSE)

tmpdata <- filter(DMR_avg1, Site=="ALAS" & Treat=="W")
ALAS_W_data <- tmpdata[,c("Chrom", "Start", "End","avgMeth")]
write.table(ALAS_W_data, file =  paste0(filename, "_ALAS_W.bedGraph"), sep="\t", quote = FALSE, row.names=FALSE, col.names=FALSE)

#LATJ
tmpdata <- filter(DMR_avg1, Site=="LATJ" & Treat=="C")
LATJ_C_data <- tmpdata[,c("Chrom", "Start", "End","avgMeth")]
write.table(LATJ_C_data, file = paste0(filename, "_LATJ_C.bedGraph"), sep="\t", quote = FALSE, row.names=FALSE, col.names=FALSE)

tmpdata <- filter(DMR_avg1, Site=="LATJ" & Treat=="W")
LATJ_W_data <- tmpdata[,c("Chrom", "Start", "End","avgMeth")]
write.table(LATJ_W_data, file =  paste0(filename, "_LATJ_W.bedGraph"), sep="\t", quote = FALSE, row.names=FALSE, col.names=FALSE)

# load into IGV 