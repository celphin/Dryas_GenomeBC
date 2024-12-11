############################
# Plotting

# PCA analysis

# Join all the data into population file and sample file

# Admixture
# - Determine Admixture groups order
# - Set groups colours
# - Plot PCA images with colours
# - Plot map of sites with dominant colours
# - Make Admixture barplot (total, Ancient, Ellemere)
# - Make Admixture maps (samples and population scales)
# - Kluane correlation and barplot

# PopStats
# - Icetime, Diversity and Temperature correlations
# - Plots of Pi and Tajima's D per population
# - Maps of all population avgeraged stats

# Table of PopStats for paper

# Site specific Admixture plots

#####################################
# download files from Cedar to computer

# Admixture
# .4.Q

# # PopStats
# Fst_mean.txt
# Fst_weighted.txt
# Fst_count.txt
# Fst_sites.txt
# population_TajimaD.txt
# summed_site_Pi.txt
# summed_wind_Pi.txt
# Het_data_by_pop.txt

# # PCA GDS file
# Cassiope_noMER_r10i.recode.gds


######################################

# # install packages
# install.packages("devtools")
# devtools::install_github("MikkoVihtakari/PlotSvalbard", upgrade = "never")
# devtools::install_github("celphin/Equitable")
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("SNPRelate")
#
# install.packages(c("tidyverse","dplyr","vctrs", "tidyr","reshape","geosphere","ade4","ape","scatterpie","data.table","maps","readr","stringr","ggplot2"))

#----------------------------------
# libraries

library(Equitable)
#library(PlotSvalbard)
library(tidyverse)
library(dplyr)
library(tidyr)
library(reshape)
library(geosphere)
library(ade4)
library(ape)
library(scatterpie)
library(data.table)
library(maps)
library(readr)
library(stringr)
library(SNPRelate)
library(ggplot2)

##############################
# Load in the data

setwd("~/GitHub/Dryas_GenomeBC/9_SNP_calling")

#Load the gds file - for PCA
genofile <- snpgdsOpen("./Data_rmbiased/Dryas_filtered_biasSNPs.gds")

# load list of samples in filtered vcf
samples_imiss <- read.table("./Data_rmbiased/Dryas_filtered_biasSNPs.imiss", header = TRUE)

# load in Admixture data - code setup for 5 populations
Admix_tbl=read.table("./Data_rmbiased/Dryas_filtered.4.Q")

# load detailed information about all 371 individuals
#sample_information  <- read.csv("./Data_rmbiased/M_caridnalis_genomics_meta_data.csv", header = TRUE)

# Population lat long
#pop_lat_long <- read.csv("./Data_rmbiased/56_pop_lat.csv", header = TRUE)

# PopStats
#Het_data <- read.table("./Data_filtered/PopStats/Het_data_by_pop.txt", header = TRUE)
#Pi_wind_data <- read.table("./Data_filtered/PopStats/summed_wind_Pi.txt", header=TRUE, sep=" ")
#Pi_site_data <- read.table("./Data_filtered/PopStats/summed_site_Pi.txt", header=TRUE, sep=" ")
#Tajima_data <- read.table("./Data_filtered/PopStats/population_TajimaD.txt", header=TRUE, sep=" ")


#################################
# PCA
# https://owensgl.github.io/biol525D/Topic_8-9/pca.html

# create GDS file on server
# snpgdsVCF2GDS("/scratch/celphin/GBS_Cassiope/Mar2020_SNPFiltering_PopStats/R_plots/FiltXXg9mac5minq30r60i_chrom_rLD.vcf.gz",
# "/scratch/celphin/GBS_Cassiope/Mar2020_SNPFiltering_PopStats/R_plots/FiltXXg9mac5minq30r60i_chrom_rLD.gds",
# method="biallelic.only")

#Prune for linkage
snpset_pruned <- snpgdsLDpruning(genofile, autosome.only=F)

#Make a list of sites we're keeping.
snpset.id <- unlist(snpset_pruned)

#Run the PCA
pca0 <- snpgdsPCA(genofile, num.thread = 1, eigen.cnt = 16, snp.id = snpset.id, missing.rate = 0.1, autosome.only = F)
pca5 <- snpgdsPCA(genofile, num.thread = 1, eigen.cnt = 16, snp.id = snpset.id, missing.rate = 0.1, maf=0.05,  autosome.only = F)

###################################
pca <- pca5

#Here's the percent variance explained for each eigenvector
pc.percent <- pca$varprop*100
round(pc.percent, 2)
# pca5: 11.29  4.73  3.71  3.22  2.90  2.62  2.48  2.00  1.84  1.74  1.67  1.60  1.53  1.52  1.45  1.40
# pca0: 10.32  4.56  3.95  3.52  2.97  2.56  2.27  2.14  1.87  1.84  1.67  1.62  1.59  1.49  1.42  1.37

#Make a dataframe of your PCA results
PCA_tab <- data.frame(sample = pca$sample.id,
                      PC1 = pca$eigenvect[,1],    # the first eigenvector
                      PC2 = pca$eigenvect[,2],    # the second eigenvector
                      PC3 = pca$eigenvect[,3],    # the first eigenvector
                      PC4 = pca$eigenvect[,4],    # the second eigenvector
                      PC5 = pca$eigenvect[,5],    # the first eigenvector
                      PC6 = pca$eigenvect[,6],    # the second eigenvector
                      stringsAsFactors = FALSE)

# Use strsplit to split by underscore '_'
split_filenames <- strsplit(as.character(PCA_tab$sample), "_")

# Extract the first part (e.g., 'Asian1' or 'C')
Treatment <- sapply(split_filenames, function(x) x[1])
Site  <- sapply(split_filenames, function(x) substr(x[2], 1, 3))
Unique.Random.ID <- sapply(split_filenames, function(x) x[4])

PCA_tab_data <- cbind(Site, Treatment, Unique.Random.ID,  PCA_tab)

colnames(PCA_tab_data) <- c("Site", "Treatment", "Unique.Random.ID", "sample", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6")

#####################################
# Join all the data

# Sample specific information

# samples_list
# Admix_tbl
# sample_information
# Het_data
# PCA_tab_data

# Join Samples information
samples_list <- PCA_tab_data[,c(1:4)]
sample = samples_list[,4]

Admix_tbl1 <- cbind(sample, Admix_tbl)
All_samples_data <- left_join(Admix_tbl1, samples_list, by="sample")
colnames(samples_imiss)[1] <- "sample"
All_samples_data <- left_join(All_samples_data, samples_imiss, by="sample")
All_samples_data <- left_join(All_samples_data, PCA_tab_data, by="sample")

#------------------------------
# join all Population information

# population_information
# Pi_wind_data
# Pi_site_data
# Tajima_data
#
# colnames(Pi_site_data) <- c("Pop" ,   "SumSitePi" , "Site_SD_Pi",  "Site_var_Pi")
# colnames(Pi_wind_data) <- c("Pop" ,   "SumWindPi",  "Wind_SD_Pi",  "Wind_var_Pi")
# colnames(Tajima_data ) <- c("Pop" ,   "TajimaAvg",  "TajimaSD",  "TajimaVar")
# population_information$Pop[which(population_information$Pop=="25")] <- "025"
#
# All_pop_data <- left_join(population_information, Tajima_data, by="Pop")
# All_pop_data <- left_join(All_pop_data, Pi_site_data, by="Pop")
# All_pop_data <- left_join(All_pop_data, Pi_wind_data, by="Pop")

#------------------------------
# average the sample information by population
Pop_avg <- All_samples_data  %>%
  group_by(Site.x) %>%
  dplyr::summarise(V1 = mean(V1, na.rm=TRUE),
                   V2 = mean(V2, na.rm=TRUE),
                   V3 = mean(V3, na.rm=TRUE),
                   V4 = mean(V4, na.rm=TRUE),
                   #V5 = mean(V5, na.rm=TRUE),
                   # V6 = mean(V6, na.rm=TRUE),
                   # V7 = mean(V7, na.rm=TRUE),
                   # V8 = mean(V8, na.rm=TRUE),
                   # V9 = mean(V9, na.rm=TRUE),
                   # V10 = mean(V10, na.rm=TRUE),
                   # Lat = mean(Lat, na.rm=TRUE),
                   # Long = mean(Long, na.rm=TRUE),
                   # Lat_shift = mean(Lat_shift, na.rm=TRUE),
                   # Long_shift = mean(Long_shift, na.rm=TRUE),
                   # Avg.Temp = mean(Avg.Temp, na.rm=TRUE),
                   # Leaf_weight = mean(Leaf_weight, na.rm=TRUE),
                   # Ice_time = mean(Ice_time, na.rm=TRUE),
                   # O.HOM. = mean(O.HOM., na.rm=TRUE),
                   # E.HOM. = mean(E.HOM., na.rm=TRUE),
                   # N_SITES = mean(N_SITES, na.rm=TRUE),
                   # FIS = mean(F, na.rm=TRUE),
                   N_MISS = mean(N_MISS, na.rm=TRUE),
                   F_MISS = mean(F_MISS, na.rm=TRUE),
                   PC1 = mean(PC1, na.rm=TRUE),
                   PC2 = mean(PC2, na.rm=TRUE),
                   PC3 = mean(PC3, na.rm=TRUE),
                   PC4 = mean(PC4, na.rm=TRUE),
                   PC5 = mean(PC5, na.rm=TRUE),
                   PC6 = mean(PC6, na.rm=TRUE),
                   Treatment = list(Treatment)
  )

Pop_avg$Pop <- Pop_avg$Site.x

All_pop_data <- Pop_avg # add lat long here for maps

############################################
# calculate the max Admix group for each location
Admix_groups <- gather(All_pop_data, Group, Amount, V1:V4, factor_key=TRUE)

maxAdmixgroup <- Admix_groups  %>%
  group_by(Site.x) %>%
  dplyr::summarise(Group = Group[which(Amount==max(Amount, na.rm=TRUE))])

Admix_groups$Mix <- paste0(Admix_groups$Group, Admix_groups$Site.x)
maxAdmixgroup$Mix <- paste0(maxAdmixgroup$Group, maxAdmixgroup$Site.x)

Final_AdmixGroups <- left_join(Admix_groups,maxAdmixgroup,  by="Mix")

Final_AdmixGroups$Site.x <- Final_AdmixGroups$Site.x.x

Max_Admix_Group <- c("V1", "V2", "V3", "V4")
Admix_Location <- c("North", "Center1", "Center2", "South")
Group_location <- as.data.frame(cbind(Max_Admix_Group, Admix_Location))

All_pop_data <- left_join(All_pop_data, Final_AdmixGroups, by="Site.x")

All_pop_data$Max_Admix_Group <- All_pop_data$Group.y
All_pop_data <- left_join(All_pop_data, Group_location, by="Max_Admix_Group")

#4 groups
Groups_summary_admix <- All_pop_data %>%
  group_by(Group.x) %>%
  dplyr::summarise(Locations = list(Site.x),
                   #Lat = mean(Lat.x, na.rm=TRUE),
                   #Long = mean(Long.x, na.rm=TRUE),
                   V1 = mean(V1, na.rm=TRUE),
                   V2 = mean(V2, na.rm=TRUE),
                   V3 = mean(V3, na.rm=TRUE),
                   V4 = mean(V4, na.rm=TRUE))


###################################
# Set up colour scheme

# map_colours_5g <- c("deepskyblue", "yellow",  "green", "red")
#
# Max_Admix_Group <- c("V1", "V2","V3", "V4")
#
# map_col <- cbind(Max_Admix_Group, map_colours_5g)
# map_col <- as.data.frame(map_col)
#
# All_pop_data <- left_join(All_pop_data, map_col, by="Max_Admix_Group")
#
# # join map_colours with individual sample data
# Map_col_pop <- as.data.frame(cbind(All_pop_data$Site.x, All_pop_data$map_colours_5g, All_pop_data$Max_Admix_Group))
#
# colnames(Map_col_pop) <- c("Site.x", "map_colours_5g", "Max_Admix_Group")
#
# All_samples_data$Pop <- All_samples_data$Site.x
#
# All_samples_data <- left_join(All_samples_data, Map_col_pop, by="Site.x")
#
# All_samples_data$map_colours_5g <- as.factor(All_samples_data$map_colours_5g)

#####################################

#Plot a PCA image

# PCA = 27.00  5.62  2.63  1.63  1.50  1.26  1.02  0.93  0.93  0.85  0.82  0.77  0.73  0.73  0.68  0.67
All_samples_data$Site.x <- as.factor(All_samples_data$Site.x)
pdf("./plots/PCA_PC1_PC2_maf0.pdf", width = 100, height = 90)
All_samples_data %>%
  ggplot(.,aes(x=PC1,y=PC2)) +
  geom_point(aes(color = Treatment.x), size=15)  +
  geom_text(aes(label = Site.x), color = "black", fontface = 2, size = 25.4/72.27*15)+
  theme_classic()+
  theme(legend.text = element_text(color = "black", size = 50),
        axis.text=element_text(size=60 ,face="bold"),
        axis.title=element_text(size=60,face="bold")) +
  #guides(colour = guide_legend(override.aes = list(size=60)))+
  scale_colour_manual(values =c("grey","deepskyblue",  "green",  "black", "orange", "red","purple"))+
  labs(y= "PC2", x = "PC1")
dev.off()

All_samples_data$Site.x <- as.factor(All_samples_data$Site.x)
pdf("./plots/PCA_PC1_PC3_maf0.pdf", width = 100, height = 90)
All_samples_data %>%
  ggplot(.,aes(x=PC1,y=PC3)) +
  geom_point(aes(color = Treatment.x), size=15)  +
  geom_text(aes(label = Site.x), color = "black", fontface = 2, size = 25.4/72.27*15)+
  theme_classic()+
  theme(legend.text = element_text(color = "black", size = 50),
        axis.text=element_text(size=60 ,face="bold"),
        axis.title=element_text(size=60,face="bold")) +
  #guides(colour = guide_legend(override.aes = list(size=60)))+
  scale_colour_manual(values =c( "grey","deepskyblue",  "green",  "black", "orange", "red","purple"))+
  labs(y= "PC3", x = "PC1")
dev.off()
#--------------------------------
W_C_samples_data <- All_samples_data[c(which(All_samples_data$Treatment.x=="W"),which(All_samples_data$Treatment.x=="C")),]

pdf("./plots/PCA_PC1_PC2_maf0_W_C.pdf", width = 100, height = 90)
W_C_samples_data %>%
  ggplot(.,aes(x=PC1,y=PC2)) +
  geom_point(aes(color = Treatment.x), size=20)  +
  geom_text(aes(label = Site.x), color = "black", fontface = 2, size = 25.4/72.27*20)+
  theme_classic()+
  theme(legend.text = element_text(color = "black", size = 50),
        axis.text=element_text(size=60 ,face="bold"),
        axis.title=element_text(size=60,face="bold")) +
  #guides(colour = guide_legend(override.aes = list(size=60)))+
  scale_colour_manual(values =c("deepskyblue",  "red"))+
  labs(y= "PC2", x = "PC1")
dev.off()


#-------------------------------

pdf("./plots/PCA_PC2_PC3_maf0_W_C.pdf", width = 100, height = 90)
W_C_samples_data %>%
  ggplot(.,aes(x=PC2,y=PC3)) +
  geom_point(aes(color = Treatment.x), size=20)  +
  geom_text(aes(label = Site.x), color = "black", fontface = 2, size = 25.4/72.27*20)+
  theme_classic()+
  theme(legend.text = element_text(color = "black", size = 50),
        axis.text=element_text(size=60 ,face="bold"),
        axis.title=element_text(size=60,face="bold")) +
  #guides(colour = guide_legend(override.aes = list(size=60)))+
  scale_colour_manual(values =c("deepskyblue",  "red"))+
  labs(y= "PC3", x = "PC2")
dev.off()


#-------------------------------

pdf("./plots/PCA_PC3_PC4_maf0_W_C.pdf", width = 100, height = 90)
W_C_samples_data %>%
  ggplot(.,aes(x=PC3,y=PC4)) +
  geom_point(aes(color = Treatment.x), size=20)  +
  geom_text(aes(label = Site.x), color = "black", fontface = 2, size = 25.4/72.27*20)+
  theme_classic()+
  theme(legend.text = element_text(color = "black", size = 50),
        axis.text=element_text(size=60 ,face="bold"),
        axis.title=element_text(size=60,face="bold")) +
  #guides(colour = guide_legend(override.aes = list(size=60)))+
  scale_colour_manual(values =c("deepskyblue",  "red"))+
  labs(y= "PC4", x = "PC3")
dev.off()


#-------------------------------

pdf("./plots/PCA_PC5_PC6_maf0_W_C.pdf", width = 100, height = 90)
W_C_samples_data %>%
  ggplot(.,aes(x=PC5,y=PC6)) +
  geom_point(aes(color = Treatment.x), size=20)  +
  geom_text(aes(label = Site.x), color = "black", fontface = 2, size = 25.4/72.27*20)+
  theme_classic()+
  theme(legend.text = element_text(color = "black", size = 50),
        axis.text=element_text(size=60 ,face="bold"),
        axis.title=element_text(size=60,face="bold")) +
  #guides(colour = guide_legend(override.aes = list(size=60)))+
  scale_colour_manual(values =c("deepskyblue",  "red"))+
  labs(y= "PC6", x = "PC5")
dev.off()



#----------------------
# Plot boxplot of PC2 and PC1 for mixed treatment and site factor


W_C_samples_data$Site_Alex <- gsub("CAS|WIL|FER|MEA|DRY", "ALEX", W_C_samples_data$Site.x)
W_C_samples_data$mix <- as.factor(paste0(W_C_samples_data$Site_Alex, "_", W_C_samples_data$Treatment.x))

pdf("./plots/PCA_PC1_W_C_diff.pdf",width = 30, height = 20)

# Create and save the plot
print(
  ggplot(data = W_C_samples_data, aes(x = mix, y = PC1, fill =Treatment.x)) +
    geom_boxplot() +
    geom_point(size = 3) +
    labs(x = "Site_Treatment", y = "PC1") +
    scale_fill_manual(values = c("skyblue", "red") ) +
    theme_classic() +
    scale_x_discrete(drop = TRUE) +  # Remove categories with NA or zero
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
)
dev.off()

pdf("./plots/PCA_PC3_W_C_diff.pdf", width = 30, height = 20)

# Create and save the plot
print(
  ggplot(data = W_C_samples_data, aes(x = mix, y = PC3, fill =Treatment.x)) +
    geom_boxplot() +
    geom_point(size = 3) +
    labs(x = "Site_Treatment", y = "PC3") +
    scale_fill_manual(values = c("skyblue", "red") ) +
    theme_classic() +
    scale_x_discrete(drop = TRUE) +  # Remove categories with NA or zero
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
)
dev.off()

pdf("./plots/PCA_PC1_W_C_treat.pdf", width = 30, height = 20)

# Create and save the plot
print(
  ggplot(data = W_C_samples_data, aes(x =  Treatment.x, y = PC1, fill =Treatment.x)) +
    geom_boxplot() +
    geom_point(size = 3) +
    labs(x = "Site_Treatment", y = "PC1") +
    scale_fill_manual(values = c("skyblue", "red") ) +
    theme_classic() +
    scale_x_discrete(drop = TRUE) +  # Remove categories with NA or zero
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
)
dev.off()

pdf("./plots/PCA_PC3_W_C_treat.pdf",  width = 30, height = 20)

# Create and save the plot
print(
  ggplot(data = W_C_samples_data, aes(x = Treatment.x, y = PC3, fill =Treatment.x)) +
    geom_boxplot() +
    geom_point(size = 3) +
    labs(x = "Site_Treatment", y = "PC3") +
    scale_fill_manual(values = c("skyblue", "red") ) +
    theme_classic() +
    scale_x_discrete(drop = TRUE) +  # Remove categories with NA or zero
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
)
dev.off()


#-----------------
# significant difference between warmed and control plants?

# Fit the MANOVA model
model <- manova(cbind(PC1, PC2, PC3, PC4, PC5, PC6) ~ Treatment.x+Site_Alex, data = W_C_samples_data)

# Print the summary of the MANOVA
summary(model)

# You can also get a detailed result using:
summary(model, test = "Wilks")  # Wilks' Lambda test is commonly used

# Df  Pillai approx F num Df den Df  Pr(>F)
# Treatment.x  1 0.12852    2.237      6     91 0.04652 *
#   Site_Alex    3 2.73980  163.205     18    279 < 2e-16 ***
#   Residuals   96
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#--------------------------
# Which PC is significant?

# Fit the ANOVA model
anova_model <- aov(PC6 ~ Treatment.x+Site_Alex, data = W_C_samples_data)

# Print the summary of the ANOVA
summary(anova_model)

# Treatment pvals
# PC1: 0.0249 *
# PC3: 0.00932 **
# PC5: 0.00880 **


pdf("./plots/PCA_PC1_PC3_maf0_W_C.pdf", width = 100, height = 90)
W_C_samples_data %>%
  ggplot(.,aes(x=PC1,y=PC3)) +
  geom_point(aes(color = Treatment.x), size=20)  +
  geom_text(aes(label = Site.x), color = "black", fontface = 2, size = 25.4/72.27*20)+
  theme_classic()+
  theme(legend.text = element_text(color = "black", size = 50),
        axis.text=element_text(size=60 ,face="bold"),
        axis.title=element_text(size=60,face="bold")) +
  #guides(colour = guide_legend(override.aes = list(size=60)))+
  scale_colour_manual(values =c("deepskyblue",  "red"))+
  labs(y= "PC3", x = "PC1")
dev.off()

pdf("./plots/PCA_PC5_PC3_maf0_W_C.pdf", width = 100, height = 90)
W_C_samples_data %>%
  ggplot(.,aes(x=PC5,y=PC3)) +
  geom_point(aes(color = Treatment.x), size=20)  +
  geom_text(aes(label = Site.x), color = "black", fontface = 2, size = 25.4/72.27*20)+
  theme_classic()+
  theme(legend.text = element_text(color = "black", size = 50),
        axis.text=element_text(size=60 ,face="bold"),
        axis.title=element_text(size=60,face="bold")) +
  #guides(colour = guide_legend(override.aes = list(size=60)))+
  scale_colour_manual(values =c("deepskyblue",  "red"))+
  labs(y= "PC3", x = "PC5")
dev.off()


############################################

#Make Admixture barplot

#Excel: =RIGHT(A2,LEN(A2) - SEARCH("_", A2, SEARCH("_", A2) + 1))

All_samples_data$Site <- All_samples_data$Site.x

All_samples_data_map <- All_samples_data

mergedAdmixTable <- All_samples_data_map[,c("Site", "Treatment.x", "sample","V1", "V2","V3", "V4")]

rownames(mergedAdmixTable) <- mergedAdmixTable$sample

barNaming <- function(vec) {
  retVec <- vec
  for(k in 2:length(vec)) {
    if(vec[k-1] == vec[k])
      retVec[k] <- ""
  }
  return(retVec)
}

#----------------------------------
#4 Groups

#ordered4 = mergedAdmixTable[order(mergedAdmixTable$Site),]

ordered0 = mergedAdmixTable[order(mergedAdmixTable$Site),]
ordered1 = ordered0[order(ordered0$V1),]
ordered2 = ordered1[order(ordered1$V2),]
ordered3 = ordered2[order(ordered2$V4),]
ordered4 = ordered3[order(ordered3$V3),]
ordered5 = ordered4#[order(ordered4$Site),]


#ordered4 = mergedAdmixTable[order(mergedAdmixTable$Latitude.x),]
ordered5$Site <-  as.factor(ordered4$Site)
ordered5$Site <-  as.factor(ordered5$Treatment.x)

#map_colours_5g <- c("deepskyblue", "yellow",  "green", "red")
map_colours_5g <- c("green",  "yellow", "deepskyblue", "red")

pdf("./plots/Site_Admix_4_bar.pdf", width = 90, height = 40)
# bottom, left, top, and right
par(mar=c(30,10,4,4))
barplot(t(as.matrix(ordered5[,c(4:7)])), col=map_colours_5g, border=NA,
        names.arg=barNaming(ordered5$sample), las=2, cex.names=2, cex.axis=6)
dev.off()

###################################
# Read in the FST data and plot along the genome

FST_Nunavut_W_C <- read.table("./Data_rmbiased/Nunavut_W_C.weir.fst", header = TRUE)
FST_Svalbard_W_C <- read.table("./Data_rmbiased/Svalbard_W_C.weir.fst", header = TRUE)
FST_Sweden_W_C <- read.table("./Data_rmbiased/Sweden_W_C.weir.fst", header = TRUE)
FST_Alaska_W_C <- read.table("./Data_rmbiased/Alaska_W_C.weir.fst", header = TRUE)


# rolling mean
bin_size=25
FST_Nunavut_W_C$mean_FST <- zoo::rollapply(FST_Nunavut_W_C$WEIR_AND_COCKERHAM_FST, width = bin_size, FUN=mean, fill=NA, na.rm=TRUE)
FST_Sweden_W_C$mean_FST <- zoo::rollapply(FST_Sweden_W_C$WEIR_AND_COCKERHAM_FST, width = bin_size, FUN=mean, fill=NA, na.rm=TRUE)
FST_Alaska_W_C$mean_FST <- zoo::rollapply(FST_Alaska_W_C$WEIR_AND_COCKERHAM_FST, width = bin_size, FUN=mean, fill=NA, na.rm=TRUE)
FST_Svalbard_W_C$mean_FST <- zoo::rollapply(FST_Svalbard_W_C$WEIR_AND_COCKERHAM_FST, width = bin_size, FUN=mean, fill=NA, na.rm=TRUE)

chromosome_order <- c("DoctH0-1", "DoctH0-2", "DoctH0-3", "DoctH0-4", "DoctH0-5", "DoctH0-6", "DoctH0-7", "DoctH0-8", "DoctH0-9")

# Ensure that 'CHROM' is a factor and set the order explicitly
FST_Svalbard_W_C$CHROM <- factor(FST_Svalbard_W_C$CHROM, levels = chromosome_order)
FST_Alaska_W_C$CHROM <- factor(FST_Alaska_W_C$CHROM,  levels = chromosome_order)
FST_Sweden_W_C$CHROM <- factor(FST_Sweden_W_C$CHROM, levels = chromosome_order)
FST_Nunavut_W_C$CHROM <- factor(FST_Nunavut_W_C$CHROM, levels = chromosome_order)

# plot
pdf("./plots/FST_Nunavut.pdf", width = 100, height = 90)
FST_Nunavut_W_C %>%
  ggplot(.,aes(x=POS,y=mean_FST)) +
  geom_point(size=15)  +
  theme_classic()+
  facet_wrap(~ CHROM, drop = TRUE, ncol = 1)
labs(y= "FST", x = "Genome_position")
dev.off()

pdf("./plots/FST_Alaska.pdf", width = 100, height = 90)
FST_Alaska_W_C %>%
  ggplot(.,aes(x=POS,y=mean_FST)) +
  geom_point(size=15)  +
  theme_classic()+
  facet_wrap(~ CHROM, drop = TRUE, ncol = 1)
labs(y= "FST", x = "Genome_position")
dev.off()

pdf("./plots/FST_Sweden.pdf", width = 100, height = 90)
FST_Sweden_W_C %>%
  ggplot(.,aes(x=POS,y=mean_FST)) +
  geom_point(size=15)  +
  theme_classic()+
  facet_wrap(~ CHROM, drop = TRUE, ncol = 1)
labs(y= "FST", x = "Genome_position")
dev.off()

pdf("./plots/FST_Svalbard.pdf", width = 100, height = 90)
FST_Svalbard_W_C %>%
  ggplot(.,aes(x=POS,y=mean_FST)) +
  geom_point(size=15)  +
  theme_classic()+
  facet_wrap(~ CHROM, drop = TRUE, ncol = 1)
labs(y= "FST", x = "Genome_position")
dev.off()

#-------------------
# plot
pdf("./plots/FST1_Nunavut.pdf", width = 100, height = 90)
FST_Nunavut_W_C %>%
  ggplot(.,aes(x=POS,y=WEIR_AND_COCKERHAM_FST)) +
  geom_point(size=15)  +
  theme_classic()+
  facet_wrap(~ CHROM, drop = TRUE, ncol = 1)
  labs(y= "FST", x = "Genome_position")
dev.off()

pdf("./plots/FST1_Alaska.pdf", width = 100, height = 90)
FST_Alaska_W_C %>%
  ggplot(.,aes(x=POS,y=WEIR_AND_COCKERHAM_FST)) +
  geom_point(size=15)  +
  theme_classic()+
  facet_wrap(~ CHROM, drop = TRUE, ncol = 1)
labs(y= "FST", x = "Genome_position")
dev.off()

pdf("./plots/FST1_Sweden.pdf", width = 100, height = 90)
FST_Sweden_W_C %>%
  ggplot(.,aes(x=POS,y=WEIR_AND_COCKERHAM_FST)) +
  geom_point(size=15)  +
  theme_classic()+
  facet_wrap(~ CHROM, drop = TRUE, ncol = 1)
labs(y= "FST", x = "Genome_position")
dev.off()

pdf("./plots/FST1_Svalbard.pdf", width = 100, height = 90)
FST_Svalbard_W_C %>%
  ggplot(.,aes(x=POS,y=WEIR_AND_COCKERHAM_FST)) +
  geom_point(size=15)  +
  theme_classic()+
  facet_wrap(~ CHROM, drop = TRUE, ncol = 1)
labs(y= "FST", x = "Genome_position")
dev.off()

#-------------------
# Join
FST_SVAL_SWED<- left_join(FST_Svalbard_W_C, FST_Sweden_W_C, by = c("CHROM", "POS"))
pdf("./plots/FST_SVal_SWED.pdf", width = 100, height = 90)
FST_SVAL_SWED %>%
  ggplot(.,aes(x=mean_FST.x,y=mean_FST.y)) +
  geom_point(size=15)  +
  theme_classic()+
  facet_wrap(~ CHROM, drop = TRUE)+
  labs(y= "Sweden FST", x = "Svalbard FST")
dev.off()

FST_ALAS_SWED<- left_join(FST_Alaska_W_C, FST_Sweden_W_C, by = c("CHROM", "POS"))
pdf("./plots/FST_ALAS_SWED.pdf", width = 100, height = 90)
FST_ALAS_SWED %>%
  ggplot(.,aes(x=mean_FST.x,y=mean_FST.y)) +
  geom_point(size=15)  +
  theme_classic()+
  facet_wrap(~ CHROM, drop = TRUE)+
  labs(y= "Sweden FST", x = "Alaska FST")
dev.off()

FST_ALAS_SVAL<- left_join(FST_Alaska_W_C, FST_Svalbard_W_C, by = c("CHROM", "POS"))
pdf("./plots/FST_ALAS_SVAL.pdf", width = 100, height = 90)
FST_ALAS_SVAL %>%
  ggplot(.,aes(x=mean_FST.x,y=mean_FST.y)) +
  geom_point(size=15)  +
  theme_classic()+
  facet_wrap(~ CHROM, drop = TRUE)+
  labs(y= "Svalbard FST", x = "Alaska FST")
dev.off()

#------------------------

FST_ALAS_SVAL<- left_join(FST_Alaska_W_C, FST_Svalbard_W_C, by = c("CHROM", "POS"))
FST_ALAS_SVAL_SWED<- left_join(FST_ALAS_SVAL, FST_Sweden_W_C, by = c("CHROM", "POS"))

FST_ALAS_SVAL_SWED[
  !is.na(FST_ALAS_SVAL_SWED$WEIR_AND_COCKERHAM_FST.x) &
    !is.na(FST_ALAS_SVAL_SWED$WEIR_AND_COCKERHAM_FST.y) &
    !is.na(FST_ALAS_SVAL_SWED$WEIR_AND_COCKERHAM_FST) &
    FST_ALAS_SVAL_SWED$WEIR_AND_COCKERHAM_FST.x > 0.01 &
    FST_ALAS_SVAL_SWED$WEIR_AND_COCKERHAM_FST.y > 0.01 &
    FST_ALAS_SVAL_SWED$WEIR_AND_COCKERHAM_FST > 0.01,
]

# CHROM      POS WEIR_AND_COCKERHAM_FST.x  mean_FST.x WEIR_AND_COCKERHAM_FST.y   mean_FST.y WEIR_AND_COCKERHAM_FST
# 6370  DoctH0-2  2280352                0.0555556 0.015279046                0.0181818  0.013643487              0.1666670
# 13790 DoctH0-3 12541693                0.1111110 0.008248637                0.1000000 -0.055504029              0.0370370
# 14559 DoctH0-3 21495626                0.1666670 0.017875012                0.4000000  0.039116880              0.1111110
# 14987 DoctH0-3 26421871                0.0277778 0.058142628                0.0333333  0.044383833              0.0370370
# 25146 DoctH0-6  9705642                0.1111110 0.086419667                0.0285714 -0.006457144              0.2222220
# 28851 DoctH0-7 10655266                0.0370370 0.011217039                0.0285714  0.017932880              0.0277778
# mean_FST
# 6370  0.12934804
# 13790 0.03779153
# 14559 0.05476821
# 14987 0.01712798
# 25146 0.20750941
# 28851 0.05656102
