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

######################################
tmux new-session -s DMR
tmux attach-session -t DMR

salloc -c1 --time 3:00:00 --mem 120000m --account def-rieseber

# loop through filtered metiliene output for each site

cd /home/celphin/projects/rpp-rieseber/celphin/Dryas/DMRs/June2022_Metilene_DMRs/

for DMR in metilene_output_June2022_maxd100_minC20_mindiff20_filtered.txt
do
chrom=$DMRcolumn1
start=$DMRcolumn2
end=$DMRcolumn3

# goes throuhg metilene input and subsets for DMR positions - should get one file per DMR
grep "^${chrom}" metilene_W_C.input | awk '^"${chrom}"[ \t].*[ \t]' "${start}<=$1 && $1>=${end}" {print} > "${chrom}_${start}_${end}.txt"

done

#########################
# test loop above
#https://www.geeksforgeeks.org/awk-command-unixlinux-examples/
#awk

#test
#chrom=QANW01012871.1	
#start=644636
#end=645586
grep "^${chrom}" metilene_W_C.input | awk '${start}<=$1 && $1>=${end} > 2 {print $0}' > "${chrom}_${start}_${end}.txt"


###########################################################
# try differential methylation analysis - plotting
# and EdgeR plotting for DMRs

# https://github.com/owensgl/biol525D/tree/master/Topic_6
# https://bioconductor.org/packages/release/bioc/html/edgeR.html
# http://combine-australia.github.io/RNAseq-R/slides/RNASeq_filtering_qc.pdf 

cd /home/celphin/projects/rpp-rieseber/celphin/Dryas/RNAseq_analysis/

# get R on server
# open R
module load StdEnv/2020
module load r/4.2.1
module load gdal
module load udunits
module load python
export R_LIBS_USER=/home/celphin/R/x86_64-pc-linux-gnu-library/4.2.1/

R

#---------------------------
# plot specific DMR

#install.packages("ggplot2")
#install.packages("ggpubr")

#libraries
library(ggplot2)
library(ggpubr)

#import data
data <- read.table("C:/Users/gaiaa/Documents/Cassandra/PhD/GenomeBC_Dryas/Wild_Dryas_WGBS/Analysis/DMRs/QANW01009334_91.txt", header = TRUE)

colnames(data) <- c("chrom", "pos", "W_ALAS0W_3_236", "W_ALAS0W_7_263", "W_ALAS0W_15_242", "W_WILL4W_417_13", "W_LATD2W_5_206", "W_LATJ_02W_193", "W_SVAL18W_18_271", "W_ALAS0W_14_249", "W_LATD2W_4_212", "W_LATC9W_11_216", "W_SVAL8W_8_272", "W_DRY9W_50_185", "W_ALAS_00W_232", "W_CASS9W_539_128", "W_FERT22W_12F_111", "W_LATD4W_9_207", "W_ALAS0W_18_248", "W_MEAD6W_466_163", "W_LATC3W_16_220", "W_WILL10W_437_84", "W_LATJ_04W_188", "W_MEAD_08W_4R0", "W_ALAS0W_16_239", "W_DRY3W_15_69", "W_WILL5W_421_154", "W_ALAS0W_17_235", "W_FERT30W_14F_40", "W_MEAD7W_470_173", "W_FERT14W_6F_126", "W_DRY6W_31_147", "W_FERT6W_3F_110", "W_SVAL16W_16_277", "W_MEAD1W_444_116", "W_DRY8W_45_155", "W_CASS17W_574_137", "W_DRY1W_3_39", "W_WILL1W_403_67", "W_LATD4W_8_211", "W_SVAL6W_6_274", "W_MEAD2W_450_70", "W_WILL7W_448_107", "W_LATC5W_18_190", "W_LATC1W_12_219", "W_CASS7W_600_19", "W_CASS10W_544_60", "W_CASS5W_525_130", "W_ALAS_00W_228", "W_SVAL_0W_267", "W_CASS_04W_519", "W_ALAS0W_8_265", "C_ALAS_00C_227", "C_LATJ_02C_194", "C_LATD1C_4_223", "C_LATD2C_7_209", "C_WILL1C_406_152", "C_DRY4C_23_82", "C_SVAL_0C_268", "C_WILL3C_414_100", "C_MEAD2C_451_76", "C_SVAL_0C_269", "C_SVAL16C_16_276", "C_LATD5C_2_201", "C_CASS17C_576_175", "C_LATD2C_1_203", "C_DRY10C_60_41", "C_SVAL12C_12_273", "C_SVAL49C_49_278", "C_LATD5C_20_199", "C_SVAL8C_8_275", "C_MEAD6C_468_22", "C_FERT5C_1F_97", "C_WILL10C_440_16", "C_ALAS0C_10_246", "C_LATD4C_3_196", "C_FERT39C_20F_71", "C_WILL7C_445_125", "C_DRY2C_10_52", "C_MEAD7C_473_95", "C_LATD5C_5_191", "C_ALAS0C_19_261", "C_WILL5C_422_31", "C_FERT31C_15F_170", "C_CASS4C_524_4", "C_FERT13C_7F_112", "C_CASS_09C_541", "C_CASS5C_529_159", "C_DRY5C_28_92", "C_LATD2C_6_198", "C_ALAS0C_4_240", "C_CASS10C_548_144", "C_ALAS0C_18_229", "C_MEAD_03C_1R0", "C_DRY9C_53_149", "C_ALAS0C_3_258", "C_ALAS0C_13_254", "C_LATJ_00C_187", "C_ALAS0C_5_238", "C_ALAS_00C_231", "C_MEAD1C_446_33", "C_CASS8C_535_54", "C_ALAS0C_12_256")

datat <-  as.data.frame(t(as.matrix(data)))

datat$ID <- c("chrom", "pos", "W_ALAS0W_3_236", "W_ALAS0W_7_263", "W_ALAS0W_15_242", "W_WILL4W_417_13", "W_LATD2W_5_206", "W_LATJ_02W_193", "W_SVAL18W_18_271", "W_ALAS0W_14_249", "W_LATD2W_4_212", "W_LATC9W_11_216", "W_SVAL8W_8_272", "W_DRY9W_50_185", "W_ALAS_00W_232", "W_CASS9W_539_128", "W_FERT22W_12F_111", "W_LATD4W_9_207", "W_ALAS0W_18_248", "W_MEAD6W_466_163", "W_LATC3W_16_220", "W_WILL10W_437_84", "W_LATJ_04W_188", "W_MEAD_08W_4R0", "W_ALAS0W_16_239", "W_DRY3W_15_69", "W_WILL5W_421_154", "W_ALAS0W_17_235", "W_FERT30W_14F_40", "W_MEAD7W_470_173", "W_FERT14W_6F_126", "W_DRY6W_31_147", "W_FERT6W_3F_110", "W_SVAL16W_16_277", "W_MEAD1W_444_116", "W_DRY8W_45_155", "W_CASS17W_574_137", "W_DRY1W_3_39", "W_WILL1W_403_67", "W_LATD4W_8_211", "W_SVAL6W_6_274", "W_MEAD2W_450_70", "W_WILL7W_448_107", "W_LATC5W_18_190", "W_LATC1W_12_219", "W_CASS7W_600_19", "W_CASS10W_544_60", "W_CASS5W_525_130", "W_ALAS_00W_228", "W_SVAL_0W_267", "W_CASS_04W_519", "W_ALAS0W_8_265", "C_ALAS_00C_227", "C_LATJ_02C_194", "C_LATD1C_4_223", "C_LATD2C_7_209", "C_WILL1C_406_152", "C_DRY4C_23_82", "C_SVAL_0C_268", "C_WILL3C_414_100", "C_MEAD2C_451_76", "C_SVAL_0C_269", "C_SVAL16C_16_276", "C_LATD5C_2_201", "C_CASS17C_576_175", "C_LATD2C_1_203", "C_DRY10C_60_41", "C_SVAL12C_12_273", "C_SVAL49C_49_278", "C_LATD5C_20_199", "C_SVAL8C_8_275", "C_MEAD6C_468_22", "C_FERT5C_1F_97", "C_WILL10C_440_16", "C_ALAS0C_10_246", "C_LATD4C_3_196", "C_FERT39C_20F_71", "C_WILL7C_445_125", "C_DRY2C_10_52", "C_MEAD7C_473_95", "C_LATD5C_5_191", "C_ALAS0C_19_261", "C_WILL5C_422_31", "C_FERT31C_15F_170", "C_CASS4C_524_4", "C_FERT13C_7F_112", "C_CASS_09C_541", "C_CASS5C_529_159", "C_DRY5C_28_92", "C_LATD2C_6_198", "C_ALAS0C_4_240", "C_CASS10C_548_144", "C_ALAS0C_18_229", "C_MEAD_03C_1R0", "C_DRY9C_53_149", "C_ALAS0C_3_258", "C_ALAS0C_13_254", "C_LATJ_00C_187", "C_ALAS0C_5_238", "C_ALAS_00C_231", "C_MEAD1C_446_33", "C_CASS8C_535_54", "C_ALAS0C_12_256")

datat <- datat[-c(1,2),]

#factors to numeric
X <- length(datat)
factorToNumeric <- function(f) as.numeric(as.character(f))
datat[1:X] <- lapply(datat[1:X], factorToNumeric)


# sum rows
datat$total <- rowSums(datat[1:X], na.rm=TRUE)
datat$ID <- rownames(datat)

mydata <- as.data.frame(cbind(datat$ID ,datat$total))

colnames(mydata) <- c("ID", "TotalMeth")

# split ID

q <- as.data.frame(t(as.matrix(as.data.frame((strsplit(as.character(mydata$ID), "_"))))))
mydata <- cbind(mydata, q)

colnames(mydata) <- c("ID", "TotalMeth","Treat", "SitePlot", "Plot", "Plant")

#makes sites
unique(mydata$SitePlot)


mydata$Site <- sub('ALAS.*', "ALAS", mydata$SitePlot, ignore.case = FALSE)
mydata$Site <- sub('LAT.*', "LATJ", mydata$Site, ignore.case = FALSE)
mydata$Site <- sub('SVAL.*', "SVAL", mydata$Site, ignore.case = FALSE)
mydata$Site <- sub('CASS.*', "CASS", mydata$Site, ignore.case = FALSE)
mydata$Site <- sub('MEAD.*', "MEAD", mydata$Site, ignore.case = FALSE)
mydata$Site <- sub('DRY.*', "DRY", mydata$Site, ignore.case = FALSE)
mydata$Site <- sub('FERT.*', "FERT", mydata$Site, ignore.case = FALSE)
mydata$Site <- sub('WILL.*', "WILL", mydata$Site, ignore.case = FALSE)

unique(mydata$Site)

mydata$TotalMeth <- as.numeric(as.character(mydata$TotalMeth))

jpeg("C:/Users/gaiaa/Documents/Cassandra/PhD/GenomeBC_Dryas/Wild_Dryas_WGBS/Analysis/DMRs/treat_methylation_Q9334region.jpg", width = 1000, height = 707)
ggboxplot(mydata, x = "Site", y = "TotalMeth",  color = "Treat", add = "jitter", shape = "Treat")
dev.off()

write.csv(mydata, file = "C:/Users/gaiaa/Documents/Cassandra/PhD/GenomeBC_Dryas/Wild_Dryas_WGBS/Analysis/DMRs/9334_senescence.csv", quote = FALSE, row.names=FALSE)

##########################
# for each DMR - take total methylation sum/total positions  - make one file of all DMRs 
# rows DMRs and columns individuals 

#import data
data <- read.table("C:/Users/gaiaa/Documents/Cassandra/PhD/GenomeBC_Dryas/Wild_Dryas_WGBS/Analysis/DMRs/*.txt", header = TRUE)

colnames(data) <- c("chrom", "pos", "W_ALAS0W_3_236", "W_ALAS0W_7_263", "W_ALAS0W_15_242", "W_WILL4W_417_13", "W_LATD2W_5_206", "W_LATJ_02W_193", "W_SVAL18W_18_271", "W_ALAS0W_14_249", "W_LATD2W_4_212", "W_LATC9W_11_216", "W_SVAL8W_8_272", "W_DRY9W_50_185", "W_ALAS_00W_232", "W_CASS9W_539_128", "W_FERT22W_12F_111", "W_LATD4W_9_207", "W_ALAS0W_18_248", "W_MEAD6W_466_163", "W_LATC3W_16_220", "W_WILL10W_437_84", "W_LATJ_04W_188", "W_MEAD_08W_4R0", "W_ALAS0W_16_239", "W_DRY3W_15_69", "W_WILL5W_421_154", "W_ALAS0W_17_235", "W_FERT30W_14F_40", "W_MEAD7W_470_173", "W_FERT14W_6F_126", "W_DRY6W_31_147", "W_FERT6W_3F_110", "W_SVAL16W_16_277", "W_MEAD1W_444_116", "W_DRY8W_45_155", "W_CASS17W_574_137", "W_DRY1W_3_39", "W_WILL1W_403_67", "W_LATD4W_8_211", "W_SVAL6W_6_274", "W_MEAD2W_450_70", "W_WILL7W_448_107", "W_LATC5W_18_190", "W_LATC1W_12_219", "W_CASS7W_600_19", "W_CASS10W_544_60", "W_CASS5W_525_130", "W_ALAS_00W_228", "W_SVAL_0W_267", "W_CASS_04W_519", "W_ALAS0W_8_265", "C_ALAS_00C_227", "C_LATJ_02C_194", "C_LATD1C_4_223", "C_LATD2C_7_209", "C_WILL1C_406_152", "C_DRY4C_23_82", "C_SVAL_0C_268", "C_WILL3C_414_100", "C_MEAD2C_451_76", "C_SVAL_0C_269", "C_SVAL16C_16_276", "C_LATD5C_2_201", "C_CASS17C_576_175", "C_LATD2C_1_203", "C_DRY10C_60_41", "C_SVAL12C_12_273", "C_SVAL49C_49_278", "C_LATD5C_20_199", "C_SVAL8C_8_275", "C_MEAD6C_468_22", "C_FERT5C_1F_97", "C_WILL10C_440_16", "C_ALAS0C_10_246", "C_LATD4C_3_196", "C_FERT39C_20F_71", "C_WILL7C_445_125", "C_DRY2C_10_52", "C_MEAD7C_473_95", "C_LATD5C_5_191", "C_ALAS0C_19_261", "C_WILL5C_422_31", "C_FERT31C_15F_170", "C_CASS4C_524_4", "C_FERT13C_7F_112", "C_CASS_09C_541", "C_CASS5C_529_159", "C_DRY5C_28_92", "C_LATD2C_6_198", "C_ALAS0C_4_240", "C_CASS10C_548_144", "C_ALAS0C_18_229", "C_MEAD_03C_1R0", "C_DRY9C_53_149", "C_ALAS0C_3_258", "C_ALAS0C_13_254", "C_LATJ_00C_187", "C_ALAS0C_5_238", "C_ALAS_00C_231", "C_MEAD1C_446_33", "C_CASS8C_535_54", "C_ALAS0C_12_256")
datat <-  as.data.frame(t(as.matrix(data)))
datat$ID <- c("chrom", "pos", "W_ALAS0W_3_236", "W_ALAS0W_7_263", "W_ALAS0W_15_242", "W_WILL4W_417_13", "W_LATD2W_5_206", "W_LATJ_02W_193", "W_SVAL18W_18_271", "W_ALAS0W_14_249", "W_LATD2W_4_212", "W_LATC9W_11_216", "W_SVAL8W_8_272", "W_DRY9W_50_185", "W_ALAS_00W_232", "W_CASS9W_539_128", "W_FERT22W_12F_111", "W_LATD4W_9_207", "W_ALAS0W_18_248", "W_MEAD6W_466_163", "W_LATC3W_16_220", "W_WILL10W_437_84", "W_LATJ_04W_188", "W_MEAD_08W_4R0", "W_ALAS0W_16_239", "W_DRY3W_15_69", "W_WILL5W_421_154", "W_ALAS0W_17_235", "W_FERT30W_14F_40", "W_MEAD7W_470_173", "W_FERT14W_6F_126", "W_DRY6W_31_147", "W_FERT6W_3F_110", "W_SVAL16W_16_277", "W_MEAD1W_444_116", "W_DRY8W_45_155", "W_CASS17W_574_137", "W_DRY1W_3_39", "W_WILL1W_403_67", "W_LATD4W_8_211", "W_SVAL6W_6_274", "W_MEAD2W_450_70", "W_WILL7W_448_107", "W_LATC5W_18_190", "W_LATC1W_12_219", "W_CASS7W_600_19", "W_CASS10W_544_60", "W_CASS5W_525_130", "W_ALAS_00W_228", "W_SVAL_0W_267", "W_CASS_04W_519", "W_ALAS0W_8_265", "C_ALAS_00C_227", "C_LATJ_02C_194", "C_LATD1C_4_223", "C_LATD2C_7_209", "C_WILL1C_406_152", "C_DRY4C_23_82", "C_SVAL_0C_268", "C_WILL3C_414_100", "C_MEAD2C_451_76", "C_SVAL_0C_269", "C_SVAL16C_16_276", "C_LATD5C_2_201", "C_CASS17C_576_175", "C_LATD2C_1_203", "C_DRY10C_60_41", "C_SVAL12C_12_273", "C_SVAL49C_49_278", "C_LATD5C_20_199", "C_SVAL8C_8_275", "C_MEAD6C_468_22", "C_FERT5C_1F_97", "C_WILL10C_440_16", "C_ALAS0C_10_246", "C_LATD4C_3_196", "C_FERT39C_20F_71", "C_WILL7C_445_125", "C_DRY2C_10_52", "C_MEAD7C_473_95", "C_LATD5C_5_191", "C_ALAS0C_19_261", "C_WILL5C_422_31", "C_FERT31C_15F_170", "C_CASS4C_524_4", "C_FERT13C_7F_112", "C_CASS_09C_541", "C_CASS5C_529_159", "C_DRY5C_28_92", "C_LATD2C_6_198", "C_ALAS0C_4_240", "C_CASS10C_548_144", "C_ALAS0C_18_229", "C_MEAD_03C_1R0", "C_DRY9C_53_149", "C_ALAS0C_3_258", "C_ALAS0C_13_254", "C_LATJ_00C_187", "C_ALAS0C_5_238", "C_ALAS_00C_231", "C_MEAD1C_446_33", "C_CASS8C_535_54", "C_ALAS0C_12_256")


library(tidyverse)
dir_path <- '~/path/to/data/directory/'
file_pattern <- 'Df\\.[0-9]\\.csv' # regex pattern to match the file name format

read_dir <- function(dir_path, file_name){
  read_csv(paste0(dir_path, file_name)) %>% 
    mutate(file_name = file_name) %>%                # add the file name as a column              
    gather(variable, value, A:B) %>%                 # convert the data from wide to long
    group_by(file_name, variable) %>% 
    summarize(sum = sum(value, na.rm = TRUE),
              min = min(value, na.rm = TRUE),
              mean = mean(value, na.rm = TRUE),
              median = median(value, na.rm = TRUE),
              max = max(value, na.rm = TRUE))
  }

df_summary <- 
  list.files(dir_path, pattern = file_pattern) %>% 
  map_df(~ read_dir(dir_path, .))


##########################
# EdgeR
# DMR plotting

#install the package edgeR
# if (!require("BiocManager", quietly = TRUE))
    # install.packages("BiocManager")
# BiocManager::install(version = "3.15")

# BiocManager::install("edgeR")
# install.packages('gplots')

#paste it in here (i.e. replace my path with yours):
setwd ("/home/celphin/projects/rpp-rieseber/celphin/Dryas/RNAseq_analysis/")

#find your working directory:
getwd()

#load the libraries you will need 
library ("edgeR")
library ("gplots")

#read in the data
mydata <- read.table("./RNA_site_tables/gene_names_expression_table1.txt", header=TRUE)

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
jpeg("./plots/W_C_RNA_MDS.jpg", width = 700, height = 500)
plotMDS (list2, col = col_treat)
dev.off()

jpeg("./plots/Site_RNA_MDS.jpg", width = 700, height = 500)
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
top <- topTags (lrt.list2, n = 40)$table # p-value all < 5e-02 

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

jpeg("./plots/W_C_RNA_heatmap.jpg", width = 3200, height = 1000)
heatmap.2 (heatmap2, scale = "row", trace = "none", Colv = FALSE, labCol = colnames(heatmap2), labRow = rownames(heatmap2), col='rainbow', dendrogram = "none", colsep =c(5), margins =c(20,70), sepcolor = "white", sepwidth = 0.1, cexRow=2, cexCol=2)
dev.off()

jpeg("./plots/W_C_RNA_heatmap_legend.jpg", width = 1000, height = 800)
heatmap.2 (heatmap2, scale = "row", trace = "none", notecex=3, notecol="black", tracecol="black", Colv = FALSE, labCol = colnames(heatmap2), labRow = rownames(heatmap2), col='rainbow', dendrogram = "none", colsep =c(5), margins =c(0.2, 0.2), sepcolor = "white", sepwidth = 0.1, cexRow=6, cexCol=6)
dev.off()

