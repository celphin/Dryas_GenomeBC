#import

setwd("~/GitHub/Dryas_GenomeBC/10_Phenotypes")

# load data
data_field_2019 <- read.csv("data_files/dryas_field_traits_2019.csv", header=TRUE)
Total_Methylation_contexts <- read.table("data_files/bismark_summary_report.txt", sep="\t", header=TRUE)
SE_methylation_contexts <-read.table("data_files/SE_bismark_summary_report.txt", sep="\t", header=TRUE)
pheno_2018 <- read.csv("data_files/2018pheno-data_219.csv", header=TRUE)[,c(1:16)]
seed_germ_2018 <- read.csv("data_files/GermData2018Dryas.csv", header=TRUE)
seed_germ_2019 <- read.csv("data_files/Survival_Jan2021.csv", header=TRUE)
DNA_extraction_traits <- read.csv("data_files/dryas_DNA_extraction_traits.csv", header=TRUE)
meta_data_2018_Alex <- read.csv("data_files/dryas_field_metadata_2018.csv", header=TRUE)

#----------------------
# # Check aphids  and phenotype data before join with sequencing data
# data_field_20191 <- data_field_2019
# data_field_2019$Treatment <- data_field_2019$C.W.G
# data_field_2019$mix <- paste0( data_field_2019$Community, "_",  data_field_2019$C.W.G)
#
# traits <- c(  "Weight", "Lf.hairs",
#             "Aphids", "Crystals.on.lv", "X..Cap.Twist.2019")
#
# # Loop over each trait and create a plot
# for (trait in traits) {
#
#   # Filter data to remove NAs from the trait
#   Y2019_data_filtered <-  data_field_2019 %>%
#     filter(!is.na(.data[[trait]]))  # Make sure to filter NAs only
#
#   # Convert the trait to numeric (if it's not already)
#   Y2019_data_filtered[[trait]] <- as.numeric(as.character(Y2019_data_filtered[[trait]]))
#
#   # Running the Kruskal-Wallis test without using .data[[trait]]
#   kruskal_test_result <- kruskal.test(as.formula(paste(trait, "~ Treatment")), data = Y2019_data_filtered)
#
#   # Viewing the results
#   print(trait)
#   print(kruskal_test_result)
#
#   # Generate the filename for the plot
#   plot_filename <- paste0("./plots/2019_data_", trait, ".pdf")
#
#   # Save the plot to a pdf file
#   pdf(plot_filename)
#
#   # Create and save the plot
#   print(
#     ggplot(data = Y2019_data_filtered, aes(x = mix, y = .data[[trait]], fill = Treatment)) +
#       geom_boxplot() +
#       geom_point(size = 3) +
#       labs(x = "Site_Treatment", y = trait) +
#       scale_fill_manual(values = c("skyblue", "red")) +
#       theme_classic() +
#       scale_x_discrete(drop = TRUE) +  # Remove categories with NA or zero
#       theme(axis.text.x = element_text(angle = 45, hjust = 1))
#   )
#
#   # Close the pdf device
#   dev.off()
#
#   # Generate the filename for the second plot
#   plot_filename <- paste0("./plots/2019_Treatment_", trait, ".pdf")
#
#   # Save the plot to a pdf file
#   pdf(plot_filename)
#
#   # Create and save the plot
#   print(
#     ggplot(data = Y2019_data_filtered, aes(x = Treatment, y = .data[[trait]], fill = Treatment)) +
#       geom_boxplot() +
#       geom_point(size = 3) +
#       labs(x = "Site_Treatment", y = trait) +
#       scale_fill_manual(values = c("skyblue", "red"))  +
#       theme_classic() +
#       scale_x_discrete(drop = TRUE) +  # Remove categories with NA or zero
#       theme(axis.text.x = element_text(angle = 45, hjust = 1))
#   )
#
#   # Close the pdf device
#   dev.off()
# }
#-----------------------
# Format names and add IDs

Total_Methylation_contexts <- rbind(Total_Methylation_contexts, SE_methylation_contexts[c(1:2),])


# string split file names to get unique random IDs
split_filenames <- strsplit(Total_Methylation_contexts$File, "_")
extracted_values <- t( sapply(split_filenames, function(x) c(x[2], x[3], x[4])))
colnames(extracted_values) <- c("Site", "Tag", "Unique.Random.ID")
Total_Methylation_contexts1 <- cbind(extracted_values, Total_Methylation_contexts)

split_filenames1 <- strsplit(SE_methylation_contexts$File, "_")
extracted_values1 <- t( sapply(split_filenames1, function(x) c(x[1],x[2], x[3], x[4], x[5], x[6])))
SE_Methylation_contexts1 <- cbind(extracted_values1, SE_methylation_contexts)

Seed_meth <- SE_Methylation_contexts1[c(12:28),]
colnames(Seed_meth)[c(1:6)] <- c("Seedling", "Site", "W.C", "GrowChamber", "Unique.Random.ID", "Pot")
Mat_Flw_Meth <- SE_Methylation_contexts1[c(3:11),]
colnames(Mat_Flw_Meth)[c(1:5)] <- c("PhenoStage", "Site", "Plot", "Unique.Random.ID", "Tag")

# get Site_treat_plot_ID value for metadata and pheno2018 to join for unique ID 2018pheno
pheno_2018$mix <- as.factor(paste0(pheno_2018$Site, "_",pheno_2018$Plot, "_", pheno_2018$C.T, "_", pheno_2018$Plant.ID))
meta_data_2018_Alex$mix <- as.factor(paste0(meta_data_2018_Alex$Community, "_",meta_data_2018_Alex$Plot, "_", meta_data_2018_Alex$C.W.G, "_", meta_data_2018_Alex$Plant.ID))

pheno_2018_IDs <- dplyr::left_join(meta_data_2018_Alex, pheno_2018, by="mix")

#----------------------------
# join by Unique.Random.ID
Total_Methylation_contexts1$Unique.Random.ID <- as.factor(Total_Methylation_contexts1$Unique.Random.ID)
pheno_2018_IDs$Unique.Random.ID <- as.factor(pheno_2018_IDs$Unique.Random.ID)
data_field_2019$Unique.Random.ID <- as.factor(data_field_2019$Unique.Random.ID)
Mat_Flw_Meth$Unique.Random.ID <- as.factor(Mat_Flw_Meth$Unique.Random.ID)
Seed_meth$Unique.Random.ID <- as.factor(Seed_meth$Unique.Random.ID)
seed_germ_2018 $Unique.Random.ID <- as.factor(seed_germ_2018 $Unique.Random.ID)
seed_germ_2019$Unique.Random.ID <- as.factor(seed_germ_2019$Unique.Random.ID)
DNA_extraction_traits$Unique.Random.ID <- as.factor(DNA_extraction_traits$Unique.Random.ID)

# join
Total_data0 <- dplyr::left_join(Total_Methylation_contexts1, pheno_2018_IDs, by="Unique.Random.ID")
Total_data1 <- dplyr::left_join(Total_data0, data_field_2019, by="Unique.Random.ID")
Total_data2 <- dplyr::left_join(Total_data1, seed_germ_2018, by="Unique.Random.ID")
Total_data3 <- dplyr::left_join(Total_data2, seed_germ_2019, by="Unique.Random.ID")
Total_data <- dplyr::left_join(Total_data3, DNA_extraction_traits, by="Unique.Random.ID")

# join seedling data with parent data methylation
Seedling_data <- dplyr::left_join(Seed_meth, Total_data, by="Unique.Random.ID")

# join mature flower data with sen methylation data and phenology of plants
Phenology_data <- dplyr::left_join(Mat_Flw_Meth, Total_data, by="Unique.Random.ID")

colnames(Phenology_data)
colnames(Seedling_data)
colnames(Total_data)

write.table(Total_data, "Wild_Plant_Traits_Merged.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(Seedling_data, "Seedling_Parent_Traits_Merged.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(Phenology_data, "MatFlw_Sen_Plant_Traits_Merged.tsv", sep = "\t", quote = FALSE, row.names = FALSE)


#------------------------------
# group and summarize
library(dplyr)

Total_data$Site <- substr(Total_data$Site.x, 1, 3)
Total_data$Treatment <- substr(Total_data$File, 1, 1)
Total_data$mix <- as.factor(paste0( Total_data$Site, "_", Total_data$Treatment))

Total_data$mix_chamber <- as.factor(paste0( Total_data$Site, "_", Total_data$Treatment, "_", Total_data$Chamber))

# calculate %meth
Total_data$percent_CpG <- (Total_data$Methylated.CpGs/(Total_data$Methylated.CpGs + Total_data$Unmethylated.CpGs))*100
Total_data$percent_chg <- (Total_data$Methylated.CpGs/(Total_data$Methylated.chgs + Total_data$Unmethylated.chgs))*100
Total_data$percent_CHH <- (Total_data$Methylated.CpGs/(Total_data$Methylated.CHHs + Total_data$Unmethylated.CHHs))*100

Total_data$TotalCs <- (Total_data$Methylated.CHHs + Total_data$Unmethylated.CHHs)+(Total_data$Methylated.chgs + Total_data$Unmethylated.chgs)+(Total_data$Methylated.CpGs + Total_data$Unmethylated.CpGs)
Total_data$TotalmCs <- Total_data$Methylated.CHHs+Total_data$Methylated.chgs+Total_data$Methylated.CpGs
Total_data$TotalCHH <- (Total_data$Methylated.CHHs + Total_data$Unmethylated.CHHs)+(Total_data$Methylated.chgs + Total_data$Unmethylated.chgs)+(Total_data$Methylated.CpGs + Total_data$Unmethylated.CpGs)
Total_data$Totalchg <- (Total_data$Methylated.CHHs + Total_data$Unmethylated.CHHs)+(Total_data$Methylated.chgs + Total_data$Unmethylated.chgs)+(Total_data$Methylated.CpGs + Total_data$Unmethylated.CpGs)
Total_data$TotalCpG <- (Total_data$Methylated.CHHs + Total_data$Unmethylated.CHHs)+(Total_data$Methylated.chgs + Total_data$Unmethylated.chgs)+(Total_data$Methylated.CpGs + Total_data$Unmethylated.CpGs)


databySite <- group_by(Total_data, by=Treatment)
#summarize


sum_treatment_summary <- summarize(databySite,
                                   Aligned.Reads = sum(Aligned.Reads, na.rm = TRUE),
                                   Total.Reads = sum(Total.Reads, na.rm = TRUE),
                                   Methylated.CpGs = sum(Methylated.CpGs, na.rm = TRUE),
                                   Methylated.chgs = sum(Methylated.chgs, na.rm = TRUE),
                                   Methylated.CHHs = sum(Methylated.CHHs, na.rm = TRUE),
                                   Unmethylated.CpGs = sum(Unmethylated.CpGs, na.rm = TRUE),
                                   Unmethylated.chgs = sum(Unmethylated.chgs, na.rm = TRUE),
                                   Unmethylated.CHHs = sum(Unmethylated.CHHs, na.rm = TRUE),
                                   percent_CpG = sum(percent_CpG, na.rm = TRUE),
                                   percent_chg = sum(percent_chg, na.rm = TRUE),
                                   percent_CHH = sum(percent_CHH, na.rm = TRUE),
                                   Height = sum(Height, na.rm = TRUE),
                                   Flwr_total = sum(Flwr_total, na.rm = TRUE),
                                   Lf.Weight.2018 = sum(Weight.x, na.rm = TRUE),
                                   Lf.hairs.field = sum(Lf.hairs, na.rm = TRUE),
                                   Aphids = sum(Aphids, na.rm = TRUE),
                                   Crystals.on.lv = sum(Crystals.on.lv, na.rm = TRUE),
                                   Seed_count2018 = sum(X..seeds..all.seeds., na.rm = TRUE),
                                   Seed_weight2018 = sum(weight.all.seeds..g.mg., na.rm = TRUE),
                                   Dryas.seed.Ripeness = sum(Dryas.seed.Ripeness, na.rm = TRUE),
                                   Seed_count2019 = sum(Number.of.Seeds, na.rm = TRUE),
                                   Seed.weight2019 = sum(Seed.weight.mg, na.rm = TRUE),
                                   Percent_germ2019 = sum(Percent_germ, na.rm = TRUE),
                                   Numb.Lv = sum(Numb.Lv, na.rm = TRUE),
                                   Hair.on.leaves = sum(Hair.on.leaves, na.rm = TRUE)
)

mean_treatment_summary <- summarize(databySite,
                                    Aligned.Reads = mean(Aligned.Reads, na.rm = TRUE),
                                    Total.Reads = mean(Total.Reads, na.rm = TRUE),
                                    Methylated.CpGs = mean(Methylated.CpGs, na.rm = TRUE),
                                    Methylated.chgs = mean(Methylated.chgs, na.rm = TRUE),
                                    Methylated.CHHs = mean(Methylated.CHHs, na.rm = TRUE),
                                    Unmethylated.CpGs = mean(Unmethylated.CpGs, na.rm = TRUE),
                                    Unmethylated.chgs = mean(Unmethylated.chgs, na.rm = TRUE),
                                    Unmethylated.CHHs = mean(Unmethylated.CHHs, na.rm = TRUE),
                                    Height = mean(Height, na.rm = TRUE),
                                    Flwr_total = mean(Flwr_total, na.rm = TRUE),
                                    Lf.Weight.2018 = mean(Weight.x, na.rm = TRUE),
                                    Lf.hairs.field = mean(Lf.hairs, na.rm = TRUE),
                                    Aphids = mean(Aphids, na.rm = TRUE),
                                    Crystals.on.lv = mean(Crystals.on.lv, na.rm = TRUE),
                                    Seed_count2018 = mean(X..seeds..all.seeds., na.rm = TRUE),
                                    Seed_weight2018 = mean(weight.all.seeds..g.mg., na.rm = TRUE),
                                    Dryas.seed.Ripeness = mean(Dryas.seed.Ripeness, na.rm = TRUE),
                                    Seed_count2019 = mean(Number.of.Seeds, na.rm = TRUE),
                                    Seed.weight2019 = mean(Seed.weight.mg, na.rm = TRUE),
                                    Percent_germ2019 = mean(Percent_germ, na.rm = TRUE),
                                    Numb.Lv = mean(Numb.Lv, na.rm = TRUE),
                                    Hair.on.leaves = mean(Hair.on.leaves, na.rm = TRUE)
)
write.table(mean_treatment_summary, "mean_treatment_summary.tsv", sep = "\t", quote = FALSE, row.names = FALSE)


#----------------
#group
databySite <- group_by(Total_data, by=mix)
#summarize


sum_Site_treatment_summary <- summarize(databySite,
                                        Aligned.Reads = sum(Aligned.Reads, na.rm = TRUE),
                                        Total.Reads = sum(Total.Reads, na.rm = TRUE),
                                        Methylated.CpGs = sum(Methylated.CpGs, na.rm = TRUE),
                                        Methylated.chgs = sum(Methylated.chgs, na.rm = TRUE),
                                        Methylated.CHHs = sum(Methylated.CHHs, na.rm = TRUE),
                                        Unmethylated.CpGs = sum(Unmethylated.CpGs, na.rm = TRUE),
                                        Unmethylated.chgs = sum(Unmethylated.chgs, na.rm = TRUE),
                                        Unmethylated.CHHs = sum(Unmethylated.CHHs, na.rm = TRUE),
                                        Height = sum(Height, na.rm = TRUE),
                                        Flwr_total = sum(Flwr_total, na.rm = TRUE),
                                        Lf.Weight.2018 = sum(Weight.x, na.rm = TRUE),
                                        Lf.hairs.field = sum(Lf.hairs, na.rm = TRUE),
                                        Aphids = sum(Aphids, na.rm = TRUE),
                                        Crystals.on.lv = sum(Crystals.on.lv, na.rm = TRUE),
                                        Seed_count2018 = sum(X..seeds..all.seeds., na.rm = TRUE),
                                        Seed_weight2018 = sum(weight.all.seeds..g.mg., na.rm = TRUE),
                                        Dryas.seed.Ripeness = sum(Dryas.seed.Ripeness, na.rm = TRUE),
                                        Seed_count2019 = sum(Number.of.Seeds, na.rm = TRUE),
                                        Seed.weight2019 = sum(Seed.weight.mg, na.rm = TRUE),
                                        Percent_germ2019 = sum(Percent_germ, na.rm = TRUE),
                                        Numb.Lv = sum(Numb.Lv, na.rm = TRUE),
                                        Hair.on.leaves = sum(Hair.on.leaves, na.rm = TRUE)
)

mean_Site_treatment_summary <- summarize(databySite,
                                        Aligned.Reads = mean(Aligned.Reads, na.rm = TRUE),
                                        Total.Reads = mean(Total.Reads, na.rm = TRUE),
                                        Methylated.CpGs = mean(Methylated.CpGs, na.rm = TRUE),
                                        Methylated.chgs = mean(Methylated.chgs, na.rm = TRUE),
                                        Methylated.CHHs = mean(Methylated.CHHs, na.rm = TRUE),
                                        Unmethylated.CpGs = mean(Unmethylated.CpGs, na.rm = TRUE),
                                        Unmethylated.chgs = mean(Unmethylated.chgs, na.rm = TRUE),
                                        Unmethylated.CHHs = mean(Unmethylated.CHHs, na.rm = TRUE),
                                        Height = mean(Height, na.rm = TRUE),
                                        Flwr_total = mean(Flwr_total, na.rm = TRUE),
                                        Lf.Weight.2018 = mean(Weight.x, na.rm = TRUE),
                                        Lf.hairs.field = mean(Lf.hairs, na.rm = TRUE),
                                        Aphids = mean(Aphids, na.rm = TRUE),
                                        Crystals.on.lv = mean(Crystals.on.lv, na.rm = TRUE),
                                        Seed_count2018 = mean(X..seeds..all.seeds., na.rm = TRUE),
                                        Seed_weight2018 = mean(weight.all.seeds..g.mg., na.rm = TRUE),
                                        Dryas.seed.Ripeness = mean(Dryas.seed.Ripeness, na.rm = TRUE),
                                        Seed_count2019 = mean(Number.of.Seeds, na.rm = TRUE),
                                        Seed.weight2019 = mean(Seed.weight.mg, na.rm = TRUE),
                                        Percent_germ2019 = mean(Percent_germ, na.rm = TRUE),
                                        Numb.Lv = mean(Numb.Lv, na.rm = TRUE),
                                        Hair.on.leaves = mean(Hair.on.leaves, na.rm = TRUE)
)

write.table(mean_Site_treatment_summary, "mean_Site_treatment_summary.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

#---------------------------------
#ggplot

# remove the extra samples
Total_data <- Total_data[-c(101, 202:203),]


# link Alex sites

Total_data$Site_Alex <- gsub("CAS|WIL|FER|MEA|DRY", "ALEX", Total_data$Site)

Total_data$mix <- as.factor(paste0(Total_data$Site_Alex, "_", Total_data$Treatment))

#----------------------
library(dplyr)
library(ggplot2)
library(car)
library(lme4)

# Define the trait
traits <- c(#"Aligned.Reads", "Total.Reads", "Methylated.CpGs", "Methylated.chgs", "Methylated.CHHs",
            #"Unmethylated.CpGs", "Unmethylated.chgs", "Unmethylated.CHHs",
            #"percent_CpG", "percent_chg", "percent_CHH",
            #"TotalCpG", "Totalchg", "TotalCHH", "TotalmCs", "TotalCs")
             "Height", "Flwr_total", "Weight.x", "Lf.hairs",
             "Aphids",  "X..seeds..all.seeds.", "weight.all.seeds..g.mg.",
             "Dryas.seed.Ripeness", "Number.of.Seeds", "Seed.weight.mg",
             "Percent_germ", "Numb.Lv", "Hair.on.leaves")

Total_data$Methylated.CHHs[which(Total_data$Methylated.CHHs>1e8)] <- NA

# Loop over each trait and create a plot
for (trait in traits) {
  # Filter data to remove NAs and zeros from the trait
  Total_data_filtered <- Total_data %>%
    filter(!is.na(.data[[trait]]) & .data[[trait]] != 0)

  # Convert the trait to numeric (if it's not already)
  Total_data_filtered[[trait]] <- as.numeric(as.character(Total_data_filtered[[trait]]))

  # Running the Kruskal-Wallis test without using .data[[trait]]
  kruskal_test_result <- kruskal.test(as.formula(paste(trait, "~ Treatment")), data = Total_data_filtered)

  # Viewing the results
  print(trait)
  print(kruskal_test_result)

  #----------------------------
  model <- lmer(as.formula(paste(trait, "~ Treatment+(1 | Site)")), data = Total_data)
  summary(model)
  #plot(resid(model))
  # Type II ANOVA
  print(Anova(model, type = 2))

  #----------------------------
  # Generate the filename for the plot
  plot_filename <- paste0("./plots/", trait, ".pdf")

  # Save the plot to a pdf file
  pdf(plot_filename)

  # Create and save the plot
  print(
    ggplot(data = Total_data_filtered, aes(x = mix, y = .data[[trait]], fill =Treatment)) +
      geom_boxplot() +
      geom_point(size = 3) +
      labs(x = "Site_Treatment", y = trait) +
      scale_fill_manual(values = c("skyblue", "red") ) +
      theme_classic() +
      scale_x_discrete(drop = TRUE) +  # Remove categories with NA or zero
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  )

  # Close the pdf device
  dev.off()

  # Generate the filename for the plot
  plot_filename <- paste0("./plots/Treatment_", trait, ".pdf")

  # Save the plot to a pdf file
  pdf(plot_filename)

  # Create and save the plot
  print(
    ggplot(data = Total_data_filtered, aes(x =Treatment, y = .data[[trait]], fill = Treatment)) +
      geom_boxplot() +
      geom_point(size = 3) +
      labs(x = "Site_Treatment", y = trait) +
      scale_fill_manual(values = c("skyblue", "red"))  +
      theme_classic() +
      scale_x_discrete(drop = TRUE) +  # Remove categories with NA or zero
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  )

  # Close the pdf device
  dev.off()
}

####################################
# check significance of total methylation in each context

# Fit the MANOVA model
model <- manova(cbind(Methylated.CpGs, Methylated.chgs, Methylated.CHHs) ~ Treatment+Site, data = Total_data)

# Print the summary of the MANOVA
summary(model)

# You can also get a detailed result using:
summary(model, test = "Wilks")  # Wilks' Lambda test is commonly used

# Df   Wilks approx F num Df den Df  Pr(>F)
# Treatment   1 0.95319   3.0937      3 189.00 0.02819 *
#   Site        7 0.38312  10.2627     21 543.26 < 2e-16 ***
#   Residuals 191
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#--------------------
# Fit the ANOVA model
anova_model <- aov(Methylated.CpGs ~ Treatment +Site, data = Total_data)

# Print the summary of the ANOVA
summary(anova_model)

# Df    Sum Sq   Mean Sq F value Pr(>F)
# Treatment     1 9.954e+13 9.954e+13   3.605 0.0591 .
# Site          7 3.184e+15 4.548e+14  16.473 <2e-16 ***
#   Residuals   191 5.274e+15 2.761e+13

#--------------------------
# Fit the ANOVA model
anova_model <- aov(Methylated.chgs ~ Treatment + Site, data = Total_data)

# Print the summary of the ANOVA
summary(anova_model)

# Df    Sum Sq   Mean Sq F value   Pr(>F)
# Treatment     1 9.104e+13 9.104e+13   7.438  0.00698 **
#   Site          7 1.295e+15 1.850e+14  15.113 1.17e-15 ***
#   Residuals   191 2.338e+15 1.224e+13

#------------------------
# Example data (replace with your own dataset)

# Fit the ANOVA model
anova_model <- aov(Methylated.CHHs ~ Treatment +Site, data = Total_data)

# Print the summary of the ANOVA
summary(anova_model)

# Df    Sum Sq   Mean Sq F value  Pr(>F)
# Treatment     1 6.243e+14 6.243e+14   5.855 0.01647 *
#   Site          7 2.526e+15 3.609e+14   3.385 0.00198 **
#   Residuals   191 2.036e+16 1.066e+14

##############################
# if not normally distributed

# Running the Kruskal-Wallis test
kruskal_test_result <- kruskal.test(Methylated.CHHs ~ Treatment, data = Total_data)

# Viewing the results
print(kruskal_test_result)

# data:  Methylated.CHHs by Treatment
# Kruskal-Wallis chi-squared = 6.0665, df = 1, p-value = 0.01378

#----------------------
# Running the Kruskal-Wallis test
kruskal_test_result <- kruskal.test(Unmethylated.CHHs ~ Treatment, data = Total_data)

# Viewing the results
print(kruskal_test_result)

#data:  Unmethylated.CHHs by Treatment
#Kruskal-Wallis chi-squared = 0.9745, df = 1, p-value = 0.3236
#---------------------
# Running the Kruskal-Wallis test
kruskal_test_result <- kruskal.test(percent_CHH ~ Treatment, data = Total_data)

# Viewing the results
print(kruskal_test_result)

#data:  percent_CHH by Treatment
#Kruskal-Wallis chi-squared = 1.3986, df = 1, p-value = 0.237

###########################################
library(lme4)
# Fit a mixed-effects model
model <- lmer(Methylated.CHHs ~ Treatment + (1 | Site), data = Total_data)

# View model summary
summary(model)

plot(resid(model))

# get significance
install.packages("car")
library(car)

# Type II ANOVA
Anova(model, type = 2)

# Treatment 5.8553  1    0.01553 *

####################################
# Test percentage specifically

# Fit the MANOVA model
model <- manova(cbind(percent_CpG,percent_chg, percent_CHH) ~ Treatment+Site, data = Total_data)

# Print the summary of the MANOVA
summary(model)

# You can also get a detailed result using:
summary(model, test = "Wilks")  # Wilks' Lambda test is commonly used

# Df   Wilks approx F num Df den Df  Pr(>F)
# Treatment   1 0.96888    3.051      2    190 0.04964 *
#   Site        7 0.37351   17.270     14    380 < 2e-16 ***
#   Residuals 191

#--------------------------
# Fit the ANOVA model
anova_model <- aov(percent_CpG ~ Treatment + Site, data = Total_data)

# Print the summary of the ANOVA
summary(anova_model)

# Df Sum Sq Mean Sq F value Pr(>F)
# Treatment     1   15.0   14.98   4.305 0.0393 *
#   Site          7 1723.2  246.18  70.741 <2e-16 ***
#   Residuals   191  664.7    3.48

#--------------------------
# Fit the ANOVA model
anova_model <- aov(percent_chg ~ Treatment + Site, data = Total_data)

# Print the summary of the ANOVA
summary(anova_model)

# Df Sum Sq Mean Sq F value Pr(>F)
# Treatment     1   6.64    6.64   4.763 0.0303 *
#   Site          7 314.97   45.00  32.269 <2e-16 ***
#   Residuals   191 266.33    1.39

#------------------------
# Example data (replace with your own dataset)

# Fit the ANOVA model
anova_model <- aov(percent_CHH ~ Treatment +Site, data = Total_data)

# Print the summary of the ANOVA
summary(anova_model)

# Df Sum Sq Mean Sq F value Pr(>F)
# Treatment     1  0.176  0.1762   3.735 0.0548 .
# Site          7  8.560  1.2228  25.922 <2e-16 ***
#   Residuals   191  9.010  0.0472

###########################################

# Fit the MANOVA model
model <- manova(cbind(Unmethylated.CpGs, Unmethylated.chgs, Unmethylated.CHHs) ~ Treatment+Site, data = Total_data)

# Print the summary of the MANOVA
summary(model)

# You can also get a detailed result using:
summary(model, test = "Wilks")  # Wilks' Lambda test is commonly used

# Df   Wilks approx F num Df den Df Pr(>F)
# Treatment   1 0.97083   1.8927      3 189.00 0.1322
# Site        7 0.19677  19.6995     21 543.26 <2e-16 ***
#   Residuals 191
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#---------------
# Fit the ANOVA model
anova_model <- aov(TotalmCs ~ Treatment +Site, data = Total_data)

# Print the summary of the ANOVA
summary(anova_model)
# Treatment     1 1.981e+15 1.981e+15   6.830  0.00968 **
#   Site          7 1.710e+16 2.443e+15   8.427 5.68e-09 ***
#   Residuals   191 5.538e+16 2.900e+14
# ---

#------------------
# Fit the ANOVA model
anova_model <- aov(TotalCs ~ Treatment +Site, data = Total_data)

# Print the summary of the ANOVA
summary(anova_model)

# Df    Sum Sq   Mean Sq F value   Pr(>F)
# Treatment     1 2.041e+16 2.041e+16   0.780    0.378
# Site          7 1.202e+18 1.717e+17   6.559 6.01e-07 ***
#   Residuals   191 5.001e+18 2.618e+16
#




##############################
# Traits measured across all sites

# Fit the MANOVA model
model <- manova(cbind(Lf.hairs,Number.of.Seeds, Percent_germ, Weight.x, Seed.weight.mg) ~ Treatment+Site, data = Total_data)

# Print the summary of the MANOVA
summary(model)

# You can also get a detailed result using:
summary(model, test = "Wilks")  # Wilks' Lambda test is commonly used

# Df   Wilks approx F num Df den Df    Pr(>F)
# Treatment  1 0.43504   7.7919      5     30 8.499e-05 ***
#   Site       2 0.67394   1.3087     10     60    0.2467
# Residuals 34
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#############################

# Fit the ANOVA model
anova_model <- aov(Lf.hairs ~ Treatment + Site, data = Total_data)

# Print the summary of the ANOVA
summary(anova_model)

# Df Sum Sq Mean Sq F value   Pr(>F)
# Treatment    1  161.8  161.78  28.704 6.95e-07 ***
#   Site         4  205.7   51.43   9.125 3.42e-06 ***
#   Residuals   86  484.7    5.64
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 108 observations deleted due to missingness
#---------------------------

# Fit the ANOVA model
anova_model <- aov(Number.of.Seeds ~ Treatment + Site, data = Total_data)

# Print the summary of the ANOVA
summary(anova_model)

# Df Sum Sq Mean Sq F value   Pr(>F)
# Treatment     1   1356  1355.6  22.440 6.14e-06 ***
#   Site          5   1337   267.5   4.428 0.000993 ***
#   Residuals   117   7068    60.4
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 76 observations deleted due to missingness
#---------------------------
# Fit the ANOVA model
anova_model <- aov(Percent_germ ~ Treatment + Site, data = Total_data)

# Print the summary of the ANOVA
summary(anova_model)
# Df Sum Sq Mean Sq F value   Pr(>F)
# Treatment     1 0.0215 0.02150   2.445 0.120931
# Site          5 0.2375 0.04751   5.401 0.000184 ***
#   Residuals   105 0.9236 0.00880
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 88 observations deleted due to missingness

###################################


# Take a look at seedlings data
# Seedling_data

#------------------------------
# group and summarize

library(dplyr)

Seedling_data$Treatment <-Seedling_data$W.C
Seedling_data$mix <- as.factor(paste0( Seedling_data$Site, "_", Seedling_data$W.C))

Seedling_data$mix_chamber <- as.factor(paste0( Seedling_data$Site, "_", Seedling_data$Treatment, "_", Seedling_data$Chamber))

#group
databySite <- group_by(Seedling_data, by=mix_chamber)
#summarize


sum_Site_treatment_summary <- summarize(databySite,
                                        Aligned.Reads.x = sum(Aligned.Reads.x, na.rm = TRUE),
                                        Total.Reads.x = sum(Total.Reads.x, na.rm = TRUE),
                                        Methylated.CpGs.x = sum(Methylated.CpGs.x, na.rm = TRUE),
                                        Methylated.chgs.x = sum(Methylated.chgs.x, na.rm = TRUE),
                                        Methylated.CHHs.x = sum(Methylated.CHHs.x, na.rm = TRUE),
                                        Unmethylated.CpGs.x = sum(Unmethylated.CpGs.x, na.rm = TRUE),
                                        Unmethylated.chgs.x = sum(Unmethylated.chgs.x, na.rm = TRUE),
                                        Unmethylated.CHHs.x = sum(Unmethylated.CHHs.x, na.rm = TRUE),

                                        Aligned.Reads.y = sum(Aligned.Reads.y, na.rm = TRUE),
                                        Total.Reads.y = sum(Total.Reads.y, na.rm = TRUE),
                                        Methylated.CpGs.y = sum(Methylated.CpGs.y, na.rm = TRUE),
                                        Methylated.chgs.y = sum(Methylated.chgs.y, na.rm = TRUE),
                                        Methylated.CHHs.y = sum(Methylated.CHHs.y, na.rm = TRUE),
                                        Unmethylated.CpGs.y = sum(Unmethylated.CpGs.y, na.rm = TRUE),
                                        Unmethylated.chgs.y = sum(Unmethylated.chgs.y, na.rm = TRUE),
                                        Unmethylated.CHHs.y = sum(Unmethylated.CHHs.y, na.rm = TRUE),


                                        Height = sum(Height, na.rm = TRUE),
                                        Flwr_total = sum(Flwr_total, na.rm = TRUE),
                                        Lf.Weight.2018 = sum(Weight.x, na.rm = TRUE),
                                        Lf.hairs.field = sum(Lf.hairs, na.rm = TRUE),
                                        Aphids = sum(Aphids, na.rm = TRUE),
                                        Crystals.on.lv = sum(Crystals.on.lv, na.rm = TRUE),
                                        Seed_count2018 = sum(X..seeds..all.seeds., na.rm = TRUE),
                                        Seed_weight2018 = sum(weight.all.seeds..g.mg., na.rm = TRUE),
                                        Dryas.seed.Ripeness = sum(Dryas.seed.Ripeness, na.rm = TRUE),
                                        Seed_count2019 = sum(Number.of.Seeds, na.rm = TRUE),
                                        Seed.weight2019 = sum(Seed.weight.mg, na.rm = TRUE),
                                        Percent_germ2019 = sum(Percent_germ, na.rm = TRUE),
                                        Numb.Lv = sum(Numb.Lv, na.rm = TRUE),
                                        Hair.on.leaves = sum(Hair.on.leaves, na.rm = TRUE)
)

mean_Site_treatment_summary <- summarize(databySite,
                                         Aligned.Reads.x = mean(Aligned.Reads.x, na.rm = TRUE),
                                         Total.Reads.x = mean(Total.Reads.x, na.rm = TRUE),
                                         Methylated.CpGs.x = mean(Methylated.CpGs.x, na.rm = TRUE),
                                         Methylated.chgs.x = mean(Methylated.chgs.x, na.rm = TRUE),
                                         Methylated.CHHs.x = mean(Methylated.CHHs.x, na.rm = TRUE),
                                         Unmethylated.CpGs.x = mean(Unmethylated.CpGs.x, na.rm = TRUE),
                                         Unmethylated.chgs.x = mean(Unmethylated.chgs.x, na.rm = TRUE),
                                         Unmethylated.CHHs.x = mean(Unmethylated.CHHs.x, na.rm = TRUE),

                                         Aligned.Reads.y = mean(Aligned.Reads.y, na.rm = TRUE),
                                         Total.Reads.y = mean(Total.Reads.y, na.rm = TRUE),
                                         Methylated.CpGs.y = mean(Methylated.CpGs.y, na.rm = TRUE),
                                         Methylated.chgs.y = mean(Methylated.chgs.y, na.rm = TRUE),
                                         Methylated.CHHs.y = mean(Methylated.CHHs.y, na.rm = TRUE),
                                         Unmethylated.CpGs.y = mean(Unmethylated.CpGs.y, na.rm = TRUE),
                                         Unmethylated.chgs.y = mean(Unmethylated.chgs.y, na.rm = TRUE),
                                         Unmethylated.CHHs.y = mean(Unmethylated.CHHs.y, na.rm = TRUE),


                                         Height = mean(Height, na.rm = TRUE),
                                         Flwr_total = mean(Flwr_total, na.rm = TRUE),
                                         Lf.Weight.2018 = mean(Weight.x, na.rm = TRUE),
                                         Lf.hairs.field = mean(Lf.hairs, na.rm = TRUE),
                                         Aphids = mean(Aphids, na.rm = TRUE),
                                         Crystals.on.lv = mean(Crystals.on.lv, na.rm = TRUE),
                                         Seed_count2018 = mean(X..seeds..all.seeds., na.rm = TRUE),
                                         Seed_weight2018 = mean(weight.all.seeds..g.mg., na.rm = TRUE),
                                         Dryas.seed.Ripeness = mean(Dryas.seed.Ripeness, na.rm = TRUE),
                                         Seed_count2019 = mean(Number.of.Seeds, na.rm = TRUE),
                                         Seed.weight2019 = mean(Seed.weight.mg, na.rm = TRUE),
                                         Percent_germ2019 = mean(Percent_germ, na.rm = TRUE),
                                         Numb.Lv = mean(Numb.Lv, na.rm = TRUE),
                                         Hair.on.leaves = mean(Hair.on.leaves, na.rm = TRUE)
)

#----------------------
library(dplyr)
library(ggplot2)

# Define the trait
traits <- c("Aligned.Reads.x", "Total.Reads.x", "Methylated.CpGs.x", "Methylated.chgs.x", "Methylated.CHHs.x",
            "Unmethylated.CpGs.x", "Unmethylated.chgs.x", "Unmethylated.CHHs.x",
            "Height", "Flwr_total", "Weight.x", "Lf.hairs",
              "X..seeds..all.seeds.", "weight.all.seeds..g.mg.",
            "Dryas.seed.Ripeness", "Number.of.Seeds", "Seed.weight.mg",
            "Percent_germ", "Numb.Lv", "Hair.on.leaves")

#"Aphids""Crystals.on.lv",

# Loop over each trait and create a plot
for (trait in traits) {
  # Filter data to remove NAs and zeros from the trait
  Seedling_data_filtered <- Seedling_data %>%
    filter(!is.na(.data[[trait]]) & .data[[trait]] != 0)

  # Convert the trait to numeric (if it's not already)
  Seedling_data_filtered[[trait]] <- as.numeric(as.character(Seedling_data_filtered[[trait]]))

  # Generate the filename for the plot
  plot_filename <- paste0("./plots_seedlings/", trait, ".pdf")

  # Save the plot to a pdf file
  pdf(plot_filename)

  # Create and save the plot
  print(
    ggplot(data = Seedling_data_filtered, aes(x = mix, y = .data[[trait]], fill =Treatment)) +
      geom_boxplot() +
      geom_point(size = 3) +
      labs(x = "Site_Treatment", y = trait) +
      scale_fill_manual(values = c("skyblue", "red") ) +
      theme_classic() +
      scale_x_discrete(drop = TRUE) +  # Remove categories with NA or zero
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  )

  # Close the pdf device
  dev.off()

  # Generate the filename for the plot
  plot_filename <- paste0("./plots_seedlings/Treatment_", trait, ".pdf")

  # Save the plot to a pdf file
  pdf(plot_filename)

  # Create and save the plot
  print(
    ggplot(data = Seedling_data_filtered, aes(x =Treatment, y = .data[[trait]], fill = Treatment)) +
      geom_boxplot() +
      geom_point(size = 3) +
      labs(x = "Site_Treatment", y = trait) +
      scale_fill_manual(values = c("skyblue", "red"))  +
      theme_classic() +
      scale_x_discrete(drop = TRUE) +  # Remove categories with NA or zero
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  )

  # Close the pdf device
  dev.off()
}

####################################
# check significance of total methylation in each context

# Fit the MANOVA model
model <- manova(cbind(Methylated.CpGs.x, Methylated.chgs.x, Methylated.CHHs.x) ~ Treatment+Site, data = Seedling_data)

# Print the summary of the MANOVA
summary(model)

# You can also get a detailed result using:
summary(model, test = "Wilks")  # Wilks' Lambda test is commonly used

# Df   Pillai approx F num Df den Df Pr(>F)
# Treatment  1 0.027282  0.26177      3     28 0.8523
# Site       2 0.279923  1.57314      6     58 0.1714
# Residuals 30






