# Dryas DNA methylation bismark report results

library(tidyverse)
library(stringr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)

mydata <- read.csv("C:/Users/gaiaa/Documents/Cassandra/PhD/GenomeBC_Dryas/Wild_Dryas_WGBS/Analysis/bismark_summary_report_May2021.csv", header = TRUE)

mydata$perCpG <- (mydata$Methylated.CpGs/(mydata$Methylated.CpGs+mydata$Unmethylated.CpGs))*100
mydata$perCHH <- (mydata$Methylated.CHHs/(mydata$Methylated.CHHs+mydata$Unmethylated.CHHs))*100
mydata$perCHG <- (mydata$Methylated.chgs/(mydata$Methylated.chgs+mydata$Unmethylated.chgs))*100


sitemethylation <- mydata %>%
  group_by(Site) %>%
  dplyr::summarise(avgCpG = mean(perCpG, na.rm=TRUE), 
                   avgCHH = mean(perCHH, na.rm=TRUE), 
                   avgCHG = mean(perCHG, na.rm=TRUE),
                   sdCpG = sd(perCpG, na.rm=TRUE), 
                   sdCHH = sd(perCHH, na.rm=TRUE), 
                   sdCHG = sd(perCHG, na.rm=TRUE))


treatmethylation <- mydata %>%
  group_by(Treatment) %>%
  dplyr::summarise(avgCpG = mean(perCpG, na.rm=TRUE), 
                   avgCHH = mean(perCHH, na.rm=TRUE), 
                   avgCHG = mean(perCHG, na.rm=TRUE),
                   sdCpG = sd(perCpG, na.rm=TRUE), 
                   sdCHH = sd(perCHH, na.rm=TRUE), 
                   sdCHG = sd(perCHG, na.rm=TRUE))

Sitetreatmethylation <- mydata %>%
  group_by(Treatment, Site) %>%
  dplyr::summarise(avgCpG = mean(perCpG, na.rm=TRUE), 
                   avgCHH = mean(perCHH, na.rm=TRUE), 
                   avgCHG = mean(perCHG, na.rm=TRUE),
                   sdCpG = sd(perCpG, na.rm=TRUE), 
                   sdCHH = sd(perCHH, na.rm=TRUE), 
                   sdCHG = sd(perCHG, na.rm=TRUE))


#ANOVA
# https://www.datanovia.com/en/lessons/anova-in-r/#prerequisites

library(tidyverse)
library(ggpubr)
library(rstatix)

# one way - just treatment

####################
# test assumption of equal variances
mydata %>% levene_test(perCHH ~ Treatment)
# if above is signifcantly different use welch_anova_test()[rstatix package]. but it is not

resCHH.aov <- mydata %>% anova_test(perCHH ~ Treatment)
# not  significant

#   Effect     DFn DFd    F     p     p<.05   ges
# 1 Treatment   1 100    2.72  0.102         0.026  

####################
# test assumption of equal variances
mydata %>% levene_test(perCHG ~ Treatment)
# if above is signifcantly different use welch_anova_test()[rstatix package]. but it is not

resCHG.aov <- mydata %>% anova_test(perCHG ~ Treatment)
# not  significant

# Effect DFn DFd     F     p p<.05   ges
# 1 Treatment   1 100 2.861 0.094       0.028


library(nlme) #linear models with random effect
library(lme4) #linear models with random effect

summary(lme(perCHG ~ Treatment, random = ~1|Site, data=mydata, na.action=na.exclude))

# Fixed effects: perCHG ~ Treatment 
# Value Std.Error DF  t-value p-value
# (Intercept) 16.068797 0.3689733 93 43.55003  0.0000
# TreatmentW   0.546255 0.2682709 93  2.03620  0.0446


######################
mydata %>% levene_test(perCpG ~ Site)
# if above is signifcantly different use welch_anova_test()[rstatix package]. but it is not

resCpG.aov <- mydata %>% anova_test(perCpG ~ Site)
# 56% of the variance in CpG is explained by the treatment
#    Effect  DFn DFd     F        p       p<.05   ges
# 1   Site   7  94    16.99    2.47e-14     *    0.559

# Pairwise comparisons # all significant with Alaska
pwc <- mydata %>% tukey_hsd(perCpG ~ Site)
pwc


# two way - treatment and site

CPGts.aov <- mydata %>% anova_test(perCpG ~ Treatment * Site)
CPGts.aov

CHGts.aov <- mydata %>% anova_test(perCHG ~ Treatment * Site)
CHGts.aov

pwc <- mydata %>% tukey_hsd(perCHG ~ Treatment * Site, p.adjust.method = "bonferroni")
pwc

CHHts.aov <- mydata %>% anova_test(perCHH ~ Treatment * Site)
CHHts.aov

###################################
# plots

mydata[(ncol(mydata)+1)] <- paste(mydata$Site,mydata$Treat, sep=".")
names(mydata)[(ncol(mydata))] <- "MixPlot"
mydata$MixPlot <- as.factor(mydata$MixPlot)

#colour
colour <- c("blue", "red","blue", "red","blue", "red","blue", "red","blue", "red","blue", "red","blue", "red","blue", "red")

# AF_WILL.W has one outlier - remove them for the plots

mydata$perCHH[which(mydata$perCHH == 11.99544)] <- NA
mydata$perCHG[which(mydata$perCHG == max(mydata$perCHG, na.rm=TRUE))] <- NA

jpeg("C:/Users/gaiaa/Documents/Cassandra/PhD/GenomeBC_Dryas/Wild_Dryas_WGBS/Methylation_contexts/tsite_methylation_CpG.jpg", width = 1200, height = 907)
X= mydata$perCpG
Y= ' perCpG'
plot(X~MixPlot, data = mydata, cex=1.2, las=2, ylim =c(min(X, na.rm=TRUE),max(X, na.rm=TRUE)), xlab="", ylab=paste0("DNA Methylation ", as.character(Y)), col=colour)
stripchart(X~MixPlot, vertical = TRUE, method = "jitter", pch = 16,cex.lab=1.5,cex.axis=1.2,
           col = "black", data = mydata, add=TRUE)
dev.off()

jpeg("C:/Users/gaiaa/Documents/Cassandra/PhD/GenomeBC_Dryas/Wild_Dryas_WGBS/Methylation_contexts/tsite_methylation_CHH.jpg", width = 1200, height = 907)
X= mydata$perCHH
Y= ' perCHH'
plot(X~MixPlot, data = mydata, cex=1.2, las=2, ylim =c(min(X, na.rm=TRUE),max(X, na.rm=TRUE)), xlab="", ylab=paste0("DNA Methylation ", as.character(Y)), col=colour)
stripchart(X~MixPlot, vertical = TRUE, method = "jitter", pch = 16,cex.lab=1.5,cex.axis=1.2,
           col = "black", data = mydata, add=TRUE)
dev.off()

jpeg("C:/Users/gaiaa/Documents/Cassandra/PhD/GenomeBC_Dryas/Wild_Dryas_WGBS/Methylation_contexts/tsite_methylation_CHG.jpg", width = 1200, height = 907)
X= mydata$perCHG
Y= ' perCHG'
plot(X~MixPlot, data = mydata, cex=1.2, las=2, ylim =c(min(X, na.rm=TRUE),max(X, na.rm=TRUE)), xlab="", ylab=paste0("DNA Methylation ", as.character(Y)), col=colour)
stripchart(X~MixPlot, vertical = TRUE, method = "jitter", pch = 16,cex.lab=1.5,cex.axis=1.2,
           col = "black", data = mydata, add=TRUE)
dev.off()


