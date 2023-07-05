###############################################################
#Notes on making box plot for desired DMR:
    #Messy for now, script describes how to do
###############################################################

#W_List
sed 's/^W..........._/"W_/g' ./W_list.txt | \
sed 's/^W.........._/"W_/g' | \
sed 's/^W........_/"W_/g' | \
sed 's/^W......._/"_/g' | \
sed 's/_R1_val_1_bismark_bt2_pe.deduplicated.bedGraph_sorted_tab.bedGraph/", /g' | \
sed ':a;N;$!ba;s/\n/ /g'

"W_MEAD2W_450_70",  "W_ALAS0W_14_249",  "W_FERT30W_14F_40",  "W_DRY8W_45_155",  "W_LATD4W_8_211",  "W_ALAS0W_8_265",  "W_LATD2W_4_212",  "W_LATJ_02W_193",  "W_SVAL18W_18_271",  "W_LATC1W_12_219",  "W_ALAS0W_18_248",  "W_CASS_04W_519",  "W_WILL5W_421_154",  "_ALAS_00W_228",  "W_ALAS0W_15_242",  "W_WILL1W_403_67",  "W_LATC5W_18_190",  "W_CASS17W_574_137",  "W_DRY3W_15_69",  "W_MEAD6W_466_163",  "W_DRY1W_3_39",  "W_LATD2W_5_206",  "W_MEAD7W_470_173",  "W_SVAL8W_8_272",  "W_ALAS0W_17_235",  "W_CASS10W_544_60",  "W_SVAL6W_6_274",  "W_FERT14W_6F_126",  "W_FERT22W_12F_111",  "W_LATC9W_11_216",  "W_DRY9W_50_185",  "W_CASS5W_525_130",  "W_WILL4W_417_13",  "_ALAS_00W_232",  "W_ALAS0W_7_263",  "W_WILL7W_448_107",  "W_CASS9W_539_128",  "W_LATC3W_16_220",  "W_SVAL_0W_267",  "W_SVAL16W_16_277",  "W_FERT6W_3F_110",  "W_LATD4W_9_207",  "W_DRY6W_31_147",  "W_CASS7W_600_19",  "W_MEAD_08W_4R0",  "W_MEAD1W_444_116",  "W_ALAS0W_3_236",  "W_ALAS0W_16_239",  "W_WILL10W_437_84",  "_LATJ_04W_188",  "W_SVAL_0W_270",
"W_MEAD2W_450_70",  "W_ALAS0W_14_249",  "W_FERT30W_14F_40",  "W_DRY8W_45_155",  "W_LATD4W_8_211",  "W_ALAS0W_8_265",  "W_LATD2W_4_212",  "W_LATJ_02W_193",  "W_SVAL18W_18_271",  "W_LATC1W_12_219",  "W_ALAS0W_18_248",  "W_CASS_04W_519",  "W_WILL5W_421_154",  "W_ALAS_00W_228",  "W_ALAS0W_15_242",  "W_WILL1W_403_67",  "W_LATC5W_18_190",  "W_CASS17W_574_137",  "W_DRY3W_15_69",  "W_MEAD6W_466_163",  "W_DRY1W_3_39",  "W_LATD2W_5_206",  "W_MEAD7W_470_173",  "W_SVAL8W_8_272",  "W_ALAS0W_17_235",  "W_CASS10W_544_60",  "W_SVAL6W_6_274",  "W_FERT14W_6F_126",  "W_FERT22W_12F_111",  "W_LATC9W_11_216",  "W_DRY9W_50_185",  "W_CASS5W_525_130",  "W_WILL4W_417_13",  "W_ALAS_00W_232",  "W_ALAS0W_7_263",  "W_WILL7W_448_107",  "W_CASS9W_539_128",  "W_LATC3W_16_220",  "W_SVAL_0W_267",  "W_SVAL16W_16_277",  "W_FERT6W_3F_110",  "W_LATD4W_9_207",  "W_DRY6W_31_147",  "W_CASS7W_600_19",  "W_MEAD_08W_4R0",  "W_MEAD1W_444_116",  "W_ALAS0W_3_236",  "W_ALAS0W_16_239",  "W_WILL10W_437_84",  "W_LATJ_04W_188",  "W_SVAL_0W_270",
#C_List
sed 's/^C..........._/"C_/g' ./C_list.txt | \
sed 's/^C.........._/"C_/g' | \
sed 's/^C........_/"C_/g' | \
sed 's/^C......._/"C_/g' | \
sed 's/_R1_val_1_bismark_bt2_pe.deduplicated.bedGraph_sorted_tab.bedGraph/", /g' | \
sed ':a;N;$!ba;s/\n/ /g'

"C_WILL5C_422_31",  "C_ALAS0C_19_261",  "C_CASS8C_535_54",  "C_FERT31C_15F_170",  "C_SVAL8C_8_275",  "C_ALAS0C_13_254",  "C_MEAD6C_468_22",  "C_ALAS0C_12_256",  "C_ALAS_00C_227",  "C_SVAL_0C_268",  "C_WILL7C_445_125",  "C_ALAS0C_10_246",  "C_CASS4C_524_4",  "C_LATD2C_6_198",  "C_SVAL16C_16_276",  "C_DRY9C_53_149",  "C_LATD5C_20_199",  "C_FERT13C_7F_112",  "C_MEAD2C_451_76",  "C_LATD5C_2_201",  "C_DRY5C_28_92",  "C_LATD5C_5_191",  "C_LATJ_02C_194",  "C_WILL10C_440_16",  "C_LATD2C_1_203",  "C_DRY10C_60_41",  "C_CASS_09C_541",  "C_SVAL49C_49_278",  "C_ALAS0C_3_258",  "C_CASS10C_548_144",  "C_WILL3C_414_100",  "C_FERT39C_20F_71",  "C_MEAD1C_446_33",  "C_CASS5C_529_159",  "C_WILL1C_406_152",  "C_FERT5C_1F_97",  "C_LATD1C_4_223",  "C_DRY4C_23_82",  "C_SVAL12C_12_273",  "C_LATJ_00C_187",  "C_LATD2C_7_209",  "C_MEAD7C_473_95",  "C_MEAD_03C_1R0",  "C_CASS17C_576_175",  "C_ALAS0C_18_229",  "C_DRY2C_10_52",  "C_ALAS_00C_231",  "C_ALAS0C_5_238",  "C_SVAL_0C_269",  "C_LATD4C_3_196",  "C_ALAS0C_4_240",

chrom="Do1_01_a00001"
start=1104484
end=1104878
cd "${directory}/data/"

grep "^${chrom}" Wild_metilene_W_C.input > "${chrom}.txt"
awk -F "\t" '{ if(($2 >= 1104484) && ($2 <= 1104878)) { print } }' "${chrom}.txt"  > "${chrom}_${start}_${end}.txt"

#Manually adjust above species for different labels
#With:
    #chrom="Do1_01_a00001"
    #start=1104484
    #end=1104878 
#run:
Rscript Box_Script_Do1_01_a00001.R

Rscript Box_Script.R "Do1_01_a00001" 1104484 1104878

#########
#Loop to get all the plots:


while read dmr
do
chrom=$(awk '{print $1}' <<< $dmr)
start=$(awk '{print $2}' <<< $dmr)
end=$(awk '{print $3}' <<< $dmr)


grep "^${chrom}" data/Wild_metilene_W_C.input > "data/${chrom}.txt"
awk -F "\t" -v start="$start" -v end="$end" '{ if(($2 >= start) && ($2 <= end)) { print } }' "data/${chrom}.txt"  > "data/${chrom}_${start}_${end}.txt"


Rscript Box_Script.R ${chrom} $start $end

done < RNA_total_subtract_W_C_Mat_Sen.bedGraph


#check plots
On local machine
scp -v msandler@cedar.computecanada.ca:/home/msandler/scratch/Wild_Metilene/Plotting/DMR_BoxPlots/plots/* .




Rscript Box_Script.R ${chrom} $start $end





############




#Messy from here below



directory=
#Use file:

#Try with:
#Do1_01_a00001   1104484 



module load StdEnv/2020
module load r/4.2.1
module load gdal
module load udunits
module load python
export R_LIBS_USER=/home/msandler/R/x86_64-pc-linux-gnu-library/4.2.1/

#Double check everything after running

R

#outputname=metilene_$h1_$h2_Sept2022_${maxdist}_${mincpgs}_${mindiff}
in_metilene="Wild_metilene_W_C.input"
h1='W'
h2='C'
maxdist=50
mincpgs=5
mindiff=0.9
threads=32
qval=0.001
chrom="Do1_01_a00001"
start=1104484
end=1104878
directory="/home/msandler/scratch/Wild_Metilene/Plotting/Boxplots_DMRs"

setwd (directory)

#install.packages('ggplot2')
#install.packages('ggpubr')

library(ggplot2)
library(ggpubr)

#import data
data <- read.table(paste0(directory,"/data/",chrom,"_",start,"_",end,".txt"), header = FALSE)

# parents
colnames(data) <- c("chrom", "pos", "W_MEAD2W_450_70",  "W_ALAS0W_14_249",  "W_FERT30W_14F_40",  "W_DRY8W_45_155",  "W_LATD4W_8_211",  "W_ALAS0W_8_265",  "W_LATD2W_4_212",  "W_LATJ_02W_193",  "W_SVAL18W_18_271",  "W_LATC1W_12_219",  "W_ALAS0W_18_248",  "W_CASS_04W_519",  "W_WILL5W_421_154",  "W_ALAS_00W_228",  "W_ALAS0W_15_242",  "W_WILL1W_403_67",  "W_LATC5W_18_190",  "W_CASS17W_574_137",  "W_DRY3W_15_69",  "W_MEAD6W_466_163",  "W_DRY1W_3_39",  "W_LATD2W_5_206",  "W_MEAD7W_470_173",  "W_SVAL8W_8_272",  "W_ALAS0W_17_235",  "W_CASS10W_544_60",  "W_SVAL6W_6_274",  "W_FERT14W_6F_126",  "W_FERT22W_12F_111",  "W_LATC9W_11_216",  "W_DRY9W_50_185",  "W_CASS5W_525_130",  "W_WILL4W_417_13",  "W_ALAS_00W_232",  "W_ALAS0W_7_263",  "W_WILL7W_448_107",  "W_CASS9W_539_128",  "W_LATC3W_16_220",  "W_SVAL_0W_267",  "W_SVAL16W_16_277",  "W_FERT6W_3F_110",  "W_LATD4W_9_207",  "W_DRY6W_31_147",  "W_CASS7W_600_19",  "W_MEAD_08W_4R0",  "W_MEAD1W_444_116",  "W_ALAS0W_3_236",  "W_ALAS0W_16_239",  "W_WILL10W_437_84",  "W_LATJ_04W_188",  "W_SVAL_0W_270", "C_WILL5C_422_31",  "C_ALAS0C_19_261",  "C_CASS8C_535_54",  "C_FERT31C_15F_170",  "C_SVAL8C_8_275",  "C_ALAS0C_13_254",  "C_MEAD6C_468_22",  "C_ALAS0C_12_256",  "C_ALAS_00C_227",  "C_SVAL_0C_268",  "C_WILL7C_445_125",  "C_ALAS0C_10_246",  "C_CASS4C_524_4",  "C_LATD2C_6_198",  "C_SVAL16C_16_276",  "C_DRY9C_53_149",  "C_LATD5C_20_199",  "C_FERT13C_7F_112",  "C_MEAD2C_451_76",  "C_LATD5C_2_201",  "C_DRY5C_28_92",  "C_LATD5C_5_191",  "C_LATJ_02C_194",  "C_WILL10C_440_16",  "C_LATD2C_1_203",  "C_DRY10C_60_41",  "C_CASS_09C_541",  "C_SVAL49C_49_278",  "C_ALAS0C_3_258",  "C_CASS10C_548_144",  "C_WILL3C_414_100",  "C_FERT39C_20F_71",  "C_MEAD1C_446_33",  "C_CASS5C_529_159",  "C_WILL1C_406_152",  "C_FERT5C_1F_97",  "C_LATD1C_4_223",  "C_DRY4C_23_82",  "C_SVAL12C_12_273",  "C_LATJ_00C_187",  "C_LATD2C_7_209",  "C_MEAD7C_473_95",  "C_MEAD_03C_1R0",  "C_CASS17C_576_175",  "C_ALAS0C_18_229",  "C_DRY2C_10_52",  "C_ALAS_00C_231",  "C_ALAS0C_5_238",  "C_SVAL_0C_269",  "C_LATD4C_3_196",  "C_ALAS0C_4_240")
datat <-  as.data.frame(t(as.matrix(data)))
datat$ID <- c("chrom", "pos", "W_MEAD2W_450_70",  "W_ALAS0W_14_249",  "W_FERT30W_14F_40",  "W_DRY8W_45_155",  "W_LATD4W_8_211",  "W_ALAS0W_8_265",  "W_LATD2W_4_212",  "W_LATJ_02W_193",  "W_SVAL18W_18_271",  "W_LATC1W_12_219",  "W_ALAS0W_18_248",  "W_CASS_04W_519",  "W_WILL5W_421_154",  "W_ALAS_00W_228",  "W_ALAS0W_15_242",  "W_WILL1W_403_67",  "W_LATC5W_18_190",  "W_CASS17W_574_137",  "W_DRY3W_15_69",  "W_MEAD6W_466_163",  "W_DRY1W_3_39",  "W_LATD2W_5_206",  "W_MEAD7W_470_173",  "W_SVAL8W_8_272",  "W_ALAS0W_17_235",  "W_CASS10W_544_60",  "W_SVAL6W_6_274",  "W_FERT14W_6F_126",  "W_FERT22W_12F_111",  "W_LATC9W_11_216",  "W_DRY9W_50_185",  "W_CASS5W_525_130",  "W_WILL4W_417_13",  "W_ALAS_00W_232",  "W_ALAS0W_7_263",  "W_WILL7W_448_107",  "W_CASS9W_539_128",  "W_LATC3W_16_220",  "W_SVAL_0W_267",  "W_SVAL16W_16_277",  "W_FERT6W_3F_110",  "W_LATD4W_9_207",  "W_DRY6W_31_147",  "W_CASS7W_600_19",  "W_MEAD_08W_4R0",  "W_MEAD1W_444_116",  "W_ALAS0W_3_236",  "W_ALAS0W_16_239",  "W_WILL10W_437_84",  "W_LATJ_04W_188",  "W_SVAL_0W_270", "C_WILL5C_422_31",  "C_ALAS0C_19_261",  "C_CASS8C_535_54",  "C_FERT31C_15F_170",  "C_SVAL8C_8_275",  "C_ALAS0C_13_254",  "C_MEAD6C_468_22",  "C_ALAS0C_12_256",  "C_ALAS_00C_227",  "C_SVAL_0C_268",  "C_WILL7C_445_125",  "C_ALAS0C_10_246",  "C_CASS4C_524_4",  "C_LATD2C_6_198",  "C_SVAL16C_16_276",  "C_DRY9C_53_149",  "C_LATD5C_20_199",  "C_FERT13C_7F_112",  "C_MEAD2C_451_76",  "C_LATD5C_2_201",  "C_DRY5C_28_92",  "C_LATD5C_5_191",  "C_LATJ_02C_194",  "C_WILL10C_440_16",  "C_LATD2C_1_203",  "C_DRY10C_60_41",  "C_CASS_09C_541",  "C_SVAL49C_49_278",  "C_ALAS0C_3_258",  "C_CASS10C_548_144",  "C_WILL3C_414_100",  "C_FERT39C_20F_71",  "C_MEAD1C_446_33",  "C_CASS5C_529_159",  "C_WILL1C_406_152",  "C_FERT5C_1F_97",  "C_LATD1C_4_223",  "C_DRY4C_23_82",  "C_SVAL12C_12_273",  "C_LATJ_00C_187",  "C_LATD2C_7_209",  "C_MEAD7C_473_95",  "C_MEAD_03C_1R0",  "C_CASS17C_576_175",  "C_ALAS0C_18_229",  "C_DRY2C_10_52",  "C_ALAS_00C_231",  "C_ALAS0C_5_238",  "C_SVAL_0C_269",  "C_LATD4C_3_196",  "C_ALAS0C_4_240")
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

# parents
colnames(mydata) <- c("ID", "TotalMeth","Treat", "SitePlot", "Plot", "Plant")



#makes sites - parents
unique(mydata$SitePlot) # parent

mydata$Site <- sub('ALAS.*', "ALASKA", mydata$SitePlot, ignore.case = FALSE)
mydata$Site <- sub('LAT.*', "SWEDEN", mydata$Site, ignore.case = FALSE)
mydata$Site <- sub('SVAL.*', "SVALBARD", mydata$Site, ignore.case = FALSE)

mydata$Site <- sub('CASS.*', "NUNAVUT", mydata$Site, ignore.case = FALSE)
mydata$Site <- sub('MEAD.*', "NUNAVUT", mydata$Site, ignore.case = FALSE)
mydata$Site <- sub('DRY.*', "NUNAVUT", mydata$Site, ignore.case = FALSE)
mydata$Site <- sub('FERT.*', "NUNAVUT", mydata$Site, ignore.case = FALSE)
mydata$Site <- sub('WILL.*', "NUNAVUT", mydata$Site, ignore.case = FALSE)


mydata$sSite <- sub('ALAS.*', "ALASKA", mydata$SitePlot, ignore.case = FALSE)
mydata$sSite <- sub('LAT.*', "SWEDEN", mydata$Site, ignore.case = FALSE)
mydata$sSite <- sub('SVAL.*', "SVALBARD", mydata$Site, ignore.case = FALSE)

mydata$sSite <- sub('CASS.*', "CASS", mydata$Site, ignore.case = FALSE)
mydata$sSite <- sub('MEAD.*', "MEAD", mydata$Site, ignore.case = FALSE)
mydata$sSite <- sub('DRY.*', "DRY", mydata$Site, ignore.case = FALSE)
mydata$sSite <- sub('FERT.*', "FERT", mydata$Site, ignore.case = FALSE)
mydata$sSite <- sub('WILL.*', "WILL", mydata$Site, ignore.case = FALSE)




# plots
mydata$TotalMeth <- as.numeric(as.character(mydata$TotalMeth))

jpeg(paste0(directory,"_",chrom,"_",start,"_",end,"_manuallyadjusted.jpg"), width = 700, height = 500)
ggboxplot(mydata, x = "Site", y = "TotalMeth",  color = "Treat", add = "jitter", shape = "Treat")
dev.off()

####
On local machine
scp -v msandler@cedar.computecanada.ca:/home/msandler/scratch/Wild_Metilene/Plotting/*.jpg .




###################
#checking items
chrom="Do1_01_a00001"
start=1104484
end=1104878
cd "${directory}/data/"

grep "^${chrom}" Wild_metilene_W_C.input > "${chrom}.txt"
awk -F "\t" '{ if(($2 >= 1104484) && ($2 <= 1104878)) { print } }' "${chrom}.txt"  > "${chrom}_${start}_${end}.txt"

R
in_metilene="Wild_metilene_W_C.input"
h1='W'
h2='C'
maxdist=50
mincpgs=5
mindiff=0.9
threads=32
qval=0.001
chrom="Do1_01_a00001"
start=1104484
end=1104878
directory="/home/msandler/scratch/Wild_Metilene/Plotting/Boxplots_DMRs"

setwd ("/home/msandler/scratch/Wild_Metilene/Plotting/Boxplots_DMRs")


library(ggplot2)
library(ggpubr)

#import data
data <- read.table(paste0(directory,"/data/",chrom,"_",start,"_",end,".txt"), header = FALSE)

#paremts
colnames(data) <- c("chrom", "pos", "W_MEAD2W_450_70",  "W_ALAS0W_14_249",  "W_FERT30W_14F_40",  "W_DRY8W_45_155",  "W_LATD4W_8_211",  "W_ALAS0W_8_265",  "W_LATD2W_4_212",  "W_LATJ_02W_193",  "W_SVAL18W_18_271",  "W_LATC1W_12_219",  "W_ALAS0W_18_248",  "W_CASS_04W_519",  "W_WILL5W_421_154",  "W_ALAS_00W_228",  "W_ALAS0W_15_242",  "W_WILL1W_403_67",  "W_LATC5W_18_190",  "W_CASS17W_574_137",  "W_DRY3W_15_69",  "W_MEAD6W_466_163",  "W_DRY1W_3_39",  "W_LATD2W_5_206",  "W_MEAD7W_470_173",  "W_SVAL8W_8_272",  "W_ALAS0W_17_235",  "W_CASS10W_544_60",  "W_SVAL6W_6_274",  "W_FERT14W_6F_126",  "W_FERT22W_12F_111",  "W_LATC9W_11_216",  "W_DRY9W_50_185",  "W_CASS5W_525_130",  "W_WILL4W_417_13",  "W_ALAS_00W_232",  "W_ALAS0W_7_263",  "W_WILL7W_448_107",  "W_CASS9W_539_128",  "W_LATC3W_16_220",  "W_SVAL_0W_267",  "W_SVAL16W_16_277",  "W_FERT6W_3F_110",  "W_LATD4W_9_207",  "W_DRY6W_31_147",  "W_CASS7W_600_19",  "W_MEAD_08W_4R0",  "W_MEAD1W_444_116",  "W_ALAS0W_3_236",  "W_ALAS0W_16_239",  "W_WILL10W_437_84",  "W_LATJ_04W_188",  "W_SVAL_0W_270", "C_WILL5C_422_31",  "C_ALAS0C_19_261",  "C_CASS8C_535_54",  "C_FERT31C_15F_170",  "C_SVAL8C_8_275",  "C_ALAS0C_13_254",  "C_MEAD6C_468_22",  "C_ALAS0C_12_256",  "C_ALAS_00C_227",  "C_SVAL_0C_268",  "C_WILL7C_445_125",  "C_ALAS0C_10_246",  "C_CASS4C_524_4",  "C_LATD2C_6_198",  "C_SVAL16C_16_276",  "C_DRY9C_53_149",  "C_LATD5C_20_199",  "C_FERT13C_7F_112",  "C_MEAD2C_451_76",  "C_LATD5C_2_201",  "C_DRY5C_28_92",  "C_LATD5C_5_191",  "C_LATJ_02C_194",  "C_WILL10C_440_16",  "C_LATD2C_1_203",  "C_DRY10C_60_41",  "C_CASS_09C_541",  "C_SVAL49C_49_278",  "C_ALAS0C_3_258",  "C_CASS10C_548_144",  "C_WILL3C_414_100",  "C_FERT39C_20F_71",  "C_MEAD1C_446_33",  "C_CASS5C_529_159",  "C_WILL1C_406_152",  "C_FERT5C_1F_97",  "C_LATD1C_4_223",  "C_DRY4C_23_82",  "C_SVAL12C_12_273",  "C_LATJ_00C_187",  "C_LATD2C_7_209",  "C_MEAD7C_473_95",  "C_MEAD_03C_1R0",  "C_CASS17C_576_175",  "C_ALAS0C_18_229",  "C_DRY2C_10_52",  "C_ALAS_00C_231",  "C_ALAS0C_5_238",  "C_SVAL_0C_269",  "C_LATD4C_3_196",  "C_ALAS0C_4_240")

datat <-  as.data.frame(t(as.matrix(data)))

datat$ID <- c("chrom", "pos", "W_MEAD2W_450_70",  "W_ALAS0W_14_249",  "W_FERT30W_14F_40",  "W_DRY8W_45_155",  "W_LATD4W_8_211",  "W_ALAS0W_8_265",  "W_LATD2W_4_212",  "W_LATJ_02W_193",  "W_SVAL18W_18_271",  "W_LATC1W_12_219",  "W_ALAS0W_18_248",  "W_CASS_04W_519",  "W_WILL5W_421_154",  "W_ALAS_00W_228",  "W_ALAS0W_15_242",  "W_WILL1W_403_67",  "W_LATC5W_18_190",  "W_CASS17W_574_137",  "W_DRY3W_15_69",  "W_MEAD6W_466_163",  "W_DRY1W_3_39",  "W_LATD2W_5_206",  "W_MEAD7W_470_173",  "W_SVAL8W_8_272",  "W_ALAS0W_17_235",  "W_CASS10W_544_60",  "W_SVAL6W_6_274",  "W_FERT14W_6F_126",  "W_FERT22W_12F_111",  "W_LATC9W_11_216",  "W_DRY9W_50_185",  "W_CASS5W_525_130",  "W_WILL4W_417_13",  "W_ALAS_00W_232",  "W_ALAS0W_7_263",  "W_WILL7W_448_107",  "W_CASS9W_539_128",  "W_LATC3W_16_220",  "W_SVAL_0W_267",  "W_SVAL16W_16_277",  "W_FERT6W_3F_110",  "W_LATD4W_9_207",  "W_DRY6W_31_147",  "W_CASS7W_600_19",  "W_MEAD_08W_4R0",  "W_MEAD1W_444_116",  "W_ALAS0W_3_236",  "W_ALAS0W_16_239",  "W_WILL10W_437_84",  "W_LATJ_04W_188",  "W_SVAL_0W_270", "C_WILL5C_422_31",  "C_ALAS0C_19_261",  "C_CASS8C_535_54",  "C_FERT31C_15F_170",  "C_SVAL8C_8_275",  "C_ALAS0C_13_254",  "C_MEAD6C_468_22",  "C_ALAS0C_12_256",  "C_ALAS_00C_227",  "C_SVAL_0C_268",  "C_WILL7C_445_125",  "C_ALAS0C_10_246",  "C_CASS4C_524_4",  "C_LATD2C_6_198",  "C_SVAL16C_16_276",  "C_DRY9C_53_149",  "C_LATD5C_20_199",  "C_FERT13C_7F_112",  "C_MEAD2C_451_76",  "C_LATD5C_2_201",  "C_DRY5C_28_92",  "C_LATD5C_5_191",  "C_LATJ_02C_194",  "C_WILL10C_440_16",  "C_LATD2C_1_203",  "C_DRY10C_60_41",  "C_CASS_09C_541",  "C_SVAL49C_49_278",  "C_ALAS0C_3_258",  "C_CASS10C_548_144",  "C_WILL3C_414_100",  "C_FERT39C_20F_71",  "C_MEAD1C_446_33",  "C_CASS5C_529_159",  "C_WILL1C_406_152",  "C_FERT5C_1F_97",  "C_LATD1C_4_223",  "C_DRY4C_23_82",  "C_SVAL12C_12_273",  "C_LATJ_00C_187",  "C_LATD2C_7_209",  "C_MEAD7C_473_95",  "C_MEAD_03C_1R0",  "C_CASS17C_576_175",  "C_ALAS0C_18_229",  "C_DRY2C_10_52",  "C_ALAS_00C_231",  "C_ALAS0C_5_238",  "C_SVAL_0C_269",  "C_LATD4C_3_196",  "C_ALAS0C_4_240")
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

# parents
colnames(mydata) <- c("ID", "TotalMeth","Treat", "SitePlot", "Plot", "Plant")

#makes sites - parents
unique(mydata$SitePlot) # parent

mydata$Site <- sub('ALAS.*', "ALASKA", mydata$SitePlot, ignore.case = FALSE)
mydata$Site <- sub('LAT.*', "SWEDEN", mydata$Site, ignore.case = FALSE)
mydata$Site <- sub('SVAL.*', "SVALBARD", mydata$Site, ignore.case = FALSE)

mydata$Site <- sub('CASS.*', "NUNAVUT", mydata$Site, ignore.case = FALSE)
mydata$Site <- sub('MEAD.*', "NUNAVUT", mydata$Site, ignore.case = FALSE)
mydata$Site <- sub('DRY.*', "NUNAVUT", mydata$Site, ignore.case = FALSE)
mydata$Site <- sub('FERT.*', "NUNAVUT", mydata$Site, ignore.case = FALSE)
mydata$Site <- sub('WILL.*', "NUNAVUT", mydata$Site, ignore.case = FALSE)


mydata$sSite <- sub('ALAS.*', "ALASKA", mydata$SitePlot, ignore.case = FALSE)
mydata$sSite <- sub('LAT.*', "SWEDEN", mydata$Site, ignore.case = FALSE)
mydata$sSite <- sub('SVAL.*', "SVALBARD", mydata$Site, ignore.case = FALSE)

mydata$sSite <- sub('CASS.*', "CASS", mydata$Site, ignore.case = FALSE)
mydata$sSite <- sub('MEAD.*', "MEAD", mydata$Site, ignore.case = FALSE)
mydata$sSite <- sub('DRY.*', "DRY", mydata$Site, ignore.case = FALSE)
mydata$sSite <- sub('FERT.*', "FERT", mydata$Site, ignore.case = FALSE)
mydata$sSite <- sub('WILL.*', "WILL", mydata$Site, ignore.case = FALSE)




# plots
mydata$TotalMeth <- as.numeric(as.character(mydata$TotalMeth))

jpeg(paste0(directory,"_",chrom,"_",start,"_",end,"_manuallyadjusted.jpg"), width = 700, height = 500)
ggboxplot(mydata, x = "Site", y = "TotalMeth",  color = "Treat", add = "jitter", shape = "Treat")
dev.off()



#Try
Rscript Box_Script.R


#make Plot unique to Site
mydata[(ncol(mydata)+1)] <- paste(mydata$Site,mydata$Treat, sep=".")
names(mydata)[(ncol(mydata))] <- "MixPlot"
unique(mydata$MixPlot)
length(unique(mydata$MixPlot))

#make Plot unique to Region
mydata[(ncol(mydata)+1)] <- paste(mydata$Region,mydata$Treat, sep=".")
names(mydata)[(ncol(mydata))] <- "RegPlot"
unique(mydata$RegPlot)
length(unique(mydata$RegPlot))

mydata$MixPlot <- as.factor(mydata$MixPlot)
mydata$RegPlot <- as.factor(mydata$RegPlot)
mydata$Treat <- as.factor(mydata$Treat)

#colour
colour <- c("blue", "red","blue", "red","blue", "red","blue", "red","blue", "red","blue", "red","blue", "red","blue", "red")

colour_reg <- c("blue", "red","blue", "red","blue", "red","blue", "red")


#     par(mar = c(bottom, left, top, right))
# default: par(mar = c(5.1, 4.1, 4.1, 2.1))

jpeg(paste0(directory,"/MixPlot_",chrom,"_",start,"_",end,".jpg"), width = 1700, height = 1500)
par(mgp=c(3,3,0))
par(mar = c(35, 20, 4.1, 2.1))
plot(TotalMeth~MixPlot, data = mydata, cex=1.5, ylim =c(min(mydata$TotalMeth, na.rm=TRUE),max(mydata$TotalMeth, na.rm=TRUE)), xlab = "", ylab="", col=colour, las=2, cex.lab=5, cex.axis=5)
stripchart(TotalMeth~MixPlot, vertical = TRUE, method = "jitter", pch = 16, cex = 3,
           col = "black", data = mydata, add=TRUE)
mtext(side=1, text="", line=30, cex=5)
mtext(side=2, text=paste0("Do1_01_a00001   1104484  "," DNA Methylation"), line=15, cex=5)
dev.off()

jpeg(paste0(directory,"/RegPlot_",chrom,"_",start,"_",end,".jpg"), width = 1700, height = 1500)
par(mgp=c(3,3,0))
par(mar = c(15, 20, 4.1, 2.1))
plot(TotalMeth~RegPlot, data = mydata, cex=1.5, ylim =c(min(mydata$TotalMeth, na.rm=TRUE),max(mydata$TotalMeth, na.rm=TRUE)), xlab = "", ylab="", col=colour_reg, las=2, cex.lab=5, cex.axis=5)
stripchart(TotalMeth~RegPlot, vertical = TRUE, method = "jitter", pch = 16, cex = 3,
           col = "black", data = mydata, add=TRUE)
mtext(side=1, text="Site", line=10, cex=5)
mtext(side=2, text=paste0(chrom," DNA Methylation"), line=15, cex=5)
dev.off()

jpeg(paste0(directory,"/TreatPlot_",chrom,"_",start,"_",end,".jpg"), width = 1700, height = 1500)
par(mgp=c(3,3,0))
par(mar = c(15, 20, 4.1, 2.1))
plot(TotalMeth~Treat, data = mydata, cex=1.5, ylim =c(min(mydata$TotalMeth, na.rm=TRUE),max(mydata$TotalMeth, na.rm=TRUE)), xlab = "", ylab="", col=colour_reg, las=2, cex.lab=5, cex.axis=5)
stripchart(TotalMeth~Treat, vertical = TRUE, method = "jitter", pch = 16, cex = 3,
           col = "black", data = mydata, add=TRUE)
mtext(side=1, text="Warming Treatment", line=10, cex=5)
mtext(side=2, text=paste0("Do1_01_a00001   1104484 "," DNA Methylation"), line=15, cex=5)
dev.off()

#################
#Try
plot(TotalMeth~MixPlot, data = mydata, cex=1.2, ylim =c(min(mydata$TotalMeth, na.rm=TRUE),max(mydata$TotalMeth, na.rm=TRUE)), xlab = "Warming MixPlot", ylab="DNA Methylation", col=colour)
stripchart(TotalMeth~MixPlot, vertical = TRUE, method = "jitter", pch = 16,cex.lab=1.5,cex.axis=1.2, col = "black", data = mydata, add=TRUE)


plot(TotalMeth~RegPlot, data = mydata, cex=1.2, ylim =c(min(mydata$TotalMeth, na.rm=TRUE),max(mydata$TotalMeth, na.rm=TRUE)), xlab = "Warming RegPlot", ylab="DNA Methylation", col=colour_reg)
stripchart(TotalMeth~MixPlot, vertical = TRUE, method = "jitter", pch = 16,cex.lab=1.5,cex.axis=1.2, col = "black", data = mydata, add=TRUE)

##################
jpeg(paste0(directory,"/MixPlot_",chrom,"_",start,"_",end,".jpg"), width = 1700, height = 1500)
par(mgp=c(3,3,0))
par(mar = c(35, 20, 4.1, 2.1))
plot(TotalMeth~MixPlot, data = mydata, cex=1.5, ylim =c(min(mydata$TotalMeth, na.rm=TRUE),max(mydata$TotalMeth, na.rm=TRUE)), xlab = "", ylab="", col=colour, las=2, cex.lab=5, cex.axis=5)
stripchart(TotalMeth~MixPlot, vertical = TRUE, method = "jitter", pch = 16, cex = 3,
           col = "black", data = mydata, add=TRUE)
mtext(side=1, text="", line=30, cex=5)
mtext(side=2, text=paste0("DMR:", chrom,"_",start,"_",end,"_DNA Methylation"), line=15, cex=5)
dev.off()


colour <- c("blue", "red","blue", "red","blue", "red","blue", "red")