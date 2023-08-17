

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

# plots
library(ggplot2)
library(ggpubr)

mydata$TotalMeth <- as.numeric(as.character(mydata$TotalMeth))

jpeg("C:/Users/gaiaa/Documents/Cassandra/PhD/GenomeBC_Dryas/Wild_Dryas_WGBS/Analysis/DMRs/treat_methylation_Q9334region.jpg", width = 1000, height = 707)
ggboxplot(mydata, x = "Site", y = "TotalMeth",  color = "Treat", add = "jitter", shape = "Treat")
dev.off()

write.csv(mydata, file = "C:/Users/gaiaa/Documents/Cassandra/PhD/GenomeBC_Dryas/Wild_Dryas_WGBS/Analysis/DMRs/9334_senescence.csv", quote = FALSE, row.names=FALSE)

##########################
########################
#######################
#import data
data <- read.table("C:/Users/gaiaa/Documents/Cassandra/PhD/GenomeBC_Dryas/Wild_Dryas_WGBS/Analysis/DMRs/QANW01001165_34.txt", header = TRUE)

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

# plots
library(ggplot2)
library(ggpubr)

mydata$TotalMeth <- as.numeric(as.character(mydata$TotalMeth))

jpeg("C:/Users/gaiaa/Documents/Cassandra/PhD/GenomeBC_Dryas/Wild_Dryas_WGBS/Analysis/DMRs/treat_methylation_1165region.jpg", width = 1000, height = 707)
ggboxplot(mydata, x = "Site", y = "TotalMeth",  color = "Treat", add = "jitter", shape = "Treat")
dev.off()




##########################
########################
#######################
#import data
data <- read.table("C:/Users/gaiaa/Documents/Cassandra/PhD/GenomeBC_Dryas/Wild_Dryas_WGBS/Analysis/DMRs/QANW01005326_165.txt", header = TRUE)

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

# plots
library(ggplot2)
library(ggpubr)

mydata$TotalMeth <- as.numeric(as.character(mydata$TotalMeth))

jpeg("C:/Users/gaiaa/Documents/Cassandra/PhD/GenomeBC_Dryas/Wild_Dryas_WGBS/Analysis/DMRs/treat_methylation_5326region.jpg", width = 1000, height = 707)
ggboxplot(mydata, x = "Site", y = "TotalMeth",  color = "Treat", add = "jitter", shape = "Treat")
dev.off()

##########################
########################
#######################
#import data
data <- read.table("C:/Users/gaiaa/Documents/Cassandra/PhD/GenomeBC_Dryas/Wild_Dryas_WGBS/Analysis/DMRs/QANW01004277_64.txt", header = TRUE)

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

# plots
library(ggplot2)
library(ggpubr)

mydata$TotalMeth <- as.numeric(as.character(mydata$TotalMeth))

jpeg("C:/Users/gaiaa/Documents/Cassandra/PhD/GenomeBC_Dryas/Wild_Dryas_WGBS/Analysis/DMRs/treat_methylation_4277region.jpg", width = 1000, height = 707)
ggboxplot(mydata, x = "Site", y = "TotalMeth",  color = "Treat", add = "jitter", shape = "Treat")
dev.off()

##########################
########################
#######################
#import data
data <- read.table("C:/Users/gaiaa/Documents/Cassandra/PhD/GenomeBC_Dryas/Wild_Dryas_WGBS/Analysis/DMRs/QANW0100494_427.txt", header = TRUE)

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

# plots
library(ggplot2)
library(ggpubr)

mydata$TotalMeth <- as.numeric(as.character(mydata$TotalMeth))

jpeg("C:/Users/gaiaa/Documents/Cassandra/PhD/GenomeBC_Dryas/Wild_Dryas_WGBS/Analysis/DMRs/treat_methylation_494region.jpg", width = 1000, height = 707)
ggboxplot(mydata, x = "Site", y = "TotalMeth",  color = "Treat", add = "jitter", shape = "Treat")
dev.off()

##########################
########################
#######################
#import data
data <- read.table("C:/Users/gaiaa/Documents/Cassandra/PhD/GenomeBC_Dryas/Wild_Dryas_WGBS/Analysis/DMRs/QANW01004202_170.txt", header = TRUE)

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

# plots
library(ggplot2)
library(ggpubr)

mydata$TotalMeth <- as.numeric(as.character(mydata$TotalMeth))

jpeg("C:/Users/gaiaa/Documents/Cassandra/PhD/GenomeBC_Dryas/Wild_Dryas_WGBS/Analysis/DMRs/treat_methylation_4202region.jpg", width = 1000, height = 707)
ggboxplot(mydata, x = "Site", y = "TotalMeth",  color = "Treat", add = "jitter", shape = "Treat")
dev.off()









##########################
########################
#######################
#import data
data <- read.table("C:/Users/gaiaa/Documents/Cassandra/PhD/GenomeBC_Dryas/Wild_Dryas_WGBS/Analysis/DMRs/QANW01009069_364.txt", header = TRUE)

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

# plots
library(ggplot2)
library(ggpubr)

mydata$TotalMeth <- as.numeric(as.character(mydata$TotalMeth))

jpeg("C:/Users/gaiaa/Documents/Cassandra/PhD/GenomeBC_Dryas/Wild_Dryas_WGBS/Analysis/DMRs/treat_methylation_9069region.jpg", width = 1000, height = 707)
ggboxplot(mydata, x = "Site", y = "TotalMeth",  color = "Treat", add = "jitter", shape = "Treat")
dev.off()

##########################
########################
#######################
#import data
data <- read.table("C:/Users/gaiaa/Documents/Cassandra/PhD/GenomeBC_Dryas/Wild_Dryas_WGBS/Analysis/DMRs/QANW01003660_10.txt", header = TRUE)

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


mydata$Region <- sub('CASS', "Alex", mydata$Site, ignore.case = FALSE)
mydata$Region <- sub('MEAD', "Alex", mydata$Region, ignore.case = FALSE)
mydata$Region <- sub('DRY', "Alex", mydata$Region, ignore.case = FALSE)
mydata$Region <- sub('FERT', "Alex", mydata$Region, ignore.case = FALSE)
mydata$Region <- sub('WILL', "Alex", mydata$Region, ignore.case = FALSE)

unique(mydata$Site)

# plots
library(ggplot2)
library(ggpubr)

mydata$TotalMeth <- as.numeric(as.character(mydata$TotalMeth))

jpeg("C:/Users/gaiaa/Documents/Cassandra/PhD/GenomeBC_Dryas/Wild_Dryas_WGBS/Analysis/DMRs/treat_methylation_3660region.jpg", width = 1000, height = 707)
ggboxplot(mydata, x = "Site", y = "TotalMeth",  color = "Treat", add = "jitter", shape = "Treat")
dev.off()

##########################
########################
#######################
#import data
data <- read.table("C:/Users/gaiaa/Documents/Cassandra/PhD/GenomeBC_Dryas/Wild_Dryas_WGBS/Analysis/DMRs/QANW01004202_170.txt", header = TRUE)

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

mydata$Region <- sub('CASS', "Alex", mydata$Site, ignore.case = FALSE)
mydata$Region <- sub('MEAD', "Alex", mydata$Region, ignore.case = FALSE)
mydata$Region <- sub('DRY', "Alex", mydata$Region, ignore.case = FALSE)
mydata$Region <- sub('FERT', "Alex", mydata$Region, ignore.case = FALSE)
mydata$Region <- sub('WILL', "Alex", mydata$Region, ignore.case = FALSE)

unique(mydata$Site)

# plots
library(ggplot2)
library(ggpubr)

mydata$TotalMeth <- as.numeric(as.character(mydata$TotalMeth))

jpeg("C:/Users/gaiaa/Documents/Cassandra/PhD/GenomeBC_Dryas/Wild_Dryas_WGBS/Analysis/DMRs/treat_methylation_4202region.jpg", width = 1000, height = 707)
ggboxplot(mydata, x = "Site", y = "TotalMeth",  color = "Treat", add = "jitter", shape = "Treat")
dev.off()


###################################################
#linear mixed effect models (lme)

library(nlme) #lme
library(lme4) #lmer
library(lmerTest) 

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

#Variables for assignment
unique(mydata$Site) #random factor - Site name
unique(mydata$Treat) #fixed effect - warming chamber(OTC) or control
unique(mydata$Plot) #random factor - Plot number
unique(mydata$Plant) # sample ID from Plot

#log transform TotalMeth
mydata$log.TotalMeth  <- log(mydata$TotalMeth)

#how does TotalMeth vary with Treat
plot(TotalMeth~Treat, data = mydata, cex=1.2,ylab="Total Methylation")
stripchart(TotalMeth~Treat, vertical = TRUE, method = "jitter", pch = 3,cex.lab=1.5,cex.axis=1.2,
           col = "blue", data = mydata, add=TRUE)

mydata$MixPlot <- as.factor(mydata$MixPlot)
mydata$RegPlot <- as.factor(mydata$RegPlot)

#how does TotalMeth vary with MixPlot
plot(TotalMeth~MixPlot, data = mydata, cex=1.2,ylab="Total Methylation")
stripchart(TotalMeth~MixPlot, vertical = TRUE, method = "jitter", pch = 3,cex.lab=1.5,cex.axis=1.2,
           col = "blue", data = mydata, add=TRUE)

#how does log(TotalMeth) vary with Treat
plot(log.TotalMeth~Treat, data = mydata, cex=1.2, xlab = "Warming Treat", ylab="log(Specific Leaf Area (mm^2/mg)")
stripchart(log.TotalMeth~Treat, vertical = TRUE, method = "jitter", pch = 3,cex.lab=1.5,cex.axis=1.2,
           col = "blue", data = mydata, add=TRUE)

#variables to include in linear mixed model
fix1.warming <- mydata$MixPlot
response.TotalMethlog <- mydata$log.TotalMeth
response.TotalMeth <- mydata$TotalMeth
rand1.Site <- mydata$Site

#linear mixed models
z <- lme(response.TotalMeth ~ fix1.warming, random = ~ 1|rand1.Site , data = mydata, na.action = na.exclude)

#with log transformed SLA values
#z1 <- lme(response.TotalMethlog ~ fix1.warming, random = ~ 1|rand1.Site , data = mydata, na.action = na.exclude)

#not including random effects - for interest
#z2 <- lm(response.TotalMethlog ~ fix1.warming, data = mydata, na.action = na.exclude)

# parameter estimates
summary(z)
#summary(z1)
#summary(z2)

# conf. intervals for fixed effects and variances
intervals(z)   
#intervals(z1) 

# Plot of residuals against predicted values
plot(z)          
#plot(z1)

#are residuals normally distributed?
hist(resid(z))
#hist(resid(z1))

# variance components for random effects  
VarCorr(z)            
#VarCorr(z1)               

# Tests of fixed effects (fitted sequentially if  more than 1 fixed effect)
anova(z)            
#anova(z1)

library(visreg)
visreg(z, xvar = "fix1.warming", xlab="Warming Treat", ylab="log(Specific Leaf Area (mm^2/mg)")

plot(fitted(z)~Treat, data = mydata, cex=1.2, ylim =c(min(mydata$TotalMeth, na.rm=TRUE),max(mydata$TotalMeth, na.rm=TRUE)), xlab = "Warming Treat", ylab="DNA Methylation")
stripchart(TotalMeth~Treat, vertical = TRUE, method = "jitter", pch = 3,cex.lab=1.5,cex.axis=1.2,
           col = "blue", data = mydata, add=TRUE)

#colour
colour <- c("blue", "red","blue", "red","blue", "red","blue", "red","blue", "red","blue", "red","blue", "red","blue", "red")

colour_reg <- c("blue", "red","blue", "red","blue", "red","blue", "red")


plot(TotalMeth~MixPlot, data = mydata, cex=1.2, ylim =c(min(mydata$TotalMeth, na.rm=TRUE),max(mydata$TotalMeth, na.rm=TRUE)), xlab = "Warming MixPlot", ylab="DNA Methylation", col=colour)
stripchart(TotalMeth~MixPlot, vertical = TRUE, method = "jitter", pch = 16,cex.lab=1.5,cex.axis=1.2,
           col = "black", data = mydata, add=TRUE)


plot(TotalMeth~RegPlot, data = mydata, cex=1.2, ylim =c(min(mydata$TotalMeth, na.rm=TRUE),max(mydata$TotalMeth, na.rm=TRUE)), xlab = "Warming RegPlot", ylab="DNA Methylation", col=colour_reg)
stripchart(TotalMeth~MixPlot, vertical = TRUE, method = "jitter", pch = 16,cex.lab=1.5,cex.axis=1.2,
           col = "black", data = mydata, add=TRUE)
