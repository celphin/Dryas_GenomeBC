#Try: 
cp /home/msandler/projects/def-rieseber/Dryas_shared_data/MS_Seedling_metilene_output_bedgraphs/SE_W_C_70_5_4_0.9_qval.0.001.bedgraph SE_W_C.bedGraph
cp /home/msandler/projects/def-rieseber/Dryas_shared_data/MS_Wild_metilene_output_bedgraphs/True_Parent_Warming_bedgraphs/P_W_C_70_5_4_0.9_qval.0.001.bedgraph P_W_C.bedGraph

#Intersect: Sweden, Parent, Seedling
module load bedtools/2.30.0
bedtools intersect -u -a SE_W_C.bedGraph -b Sweden_W_C.bedGraph| bedtools intersect -u -a stdin -b P_W_C.bedGraph > Swed_SE_P_W_C_intersect.bedGraph


bedtools intersect -u -a Swed_SE_P_W_C_intersect.bedGraph -b Alaska_W_C.bedGraph > Alas_Swed_SE_P_W_C_intersect.bedGraph
bedtools intersect -u -a Swed_SE_P_W_C_intersect.bedGraph -b Nunavut_W_C.bedGraph > Alex_Swed_SE_P_W_C_intersect.bedGraph
bedtools intersect -u -a Swed_SE_P_W_C_intersect.bedGraph -b Svalbard_W_C.bedGraph > Sval_Swed_SE_P_W_C_intersect.bedGraph

wc -l *Swed*intersect*
13 Alas_Swed_SE_P_W_C_intersect.bedGraph
8 Alex_Swed_SE_P_W_C_intersect.bedGraph
9 Sval_Swed_SE_P_W_C_intersect.bedGraph
74 Swed_SE_P_W_C_intersect.bedGraph


###############################################################


#Rerun metilene prep for Seedling W/C
#Rerun metilene prep for True Parent W/C

#Combine: 
#Subset intersections:
#Alas_Swed_SE_P_W_C_intersect.bedGraph
#awk -F '\t' 'FNR == NR { a[$1]; next } $1 in a' P_metilene_W_C.input Alas_Swed_SE_P_W_C_intersect.bedGraph  > subset_Parents.tsv
#awk -F "\t" -v start="$start" -v end="$end" '{ if(($2 >= start) && ($2 <= end)) { print } }' "data/${chrom}.txt"  > "data/${chrom}_${start}_${end}.txt"
#awk -F '\t' 'FNR == NR { a[$1]; next } $1 in a' Alas_Swed_SE_P_W_C_intersect.bedGraph SE_metilene_W_C.input > subset_Seedlings.tsv


#grep "^${chrom}" data/Sweden_metilene_W_C.input > "data/${chrom}.txt"
#awk -F "\t" -v start="$start" -v end="$end" '{ if(($2 >= start) && ($2 <= end)) { print } }' "data/${chrom}.txt"  > "data/${chrom}_${start}_${end}.txt"

while read dmr
do

chrom=$(awk '{print $1}' <<< $dmr)
start=$(awk '{print $2}' <<< $dmr)
end=$(awk '{print $3}' <<< $dmr)


grep "^${chrom}" data/P_metilene_W_C.input > "data/P_${chrom}.txt"
grep "^${chrom}" data/SE_metilene_W_C.input > "data/SE_${chrom}.txt"

awk -F "\t" -v start="$start" -v end="$end" '{ if(($2 >= start) && ($2 <= end)) { print } }' "data/P_${chrom}.txt"  > "data/P_${chrom}_${start}_${end}.txt"
awk -F "\t" -v start="$start" -v end="$end" '{ if(($2 >= start) && ($2 <= end)) { print } }' "data/SE_${chrom}.txt"  > "data/SE_${chrom}_${start}_${end}.txt"


paste "data/P_${chrom}_${start}_${end}.txt"  <(cut -f3- "data/SE_${chrom}_${start}_${end}.txt") > "data/${chrom}_${start}_${end}.txt"
Rscript Box_Script_SE_P.R ${chrom} $start $end

done < Alas_Swed_SE_P_W_C_intersect.bedGraph

head -n 1 *.input
==> P_metilene_W_C.input <==
chrom   pos     W2.E12.W2e5_LATC1W_12_219       W2.5.3f6_ALAS_00W_228   W2.H10.W2h1_LATC5W_18_190       C1.C10.C1f12_FERT31C_15F_170     C1.F08.C1a10_WILL7C_445_125     C2.C11.C2h2_LATD2C_6_198        C2.A11.C2a2_LATD5C_5_191 C2.25.3d2_LATJ_02C_194  C2.27.3e1_LATJ_00C_187  C2.C01.C2h7_ALAS0C_5_238

#For Parent:
/home/msandler/scratch/Parent_Metilene/P_W_C_input_files

sed 's/^W..........._/"W_/g' ./W_list.txt | \
sed 's/^W.........._/"W_/g' | \
sed 's/^W........_/"W_/g' | \
sed 's/^W......._/"_/g' | \
sed 's/_R1_val_1_bismark_bt2_pe.deduplicated.bedGraph_sorted_tab.bedGraph/", /g' | \
sed ':a;N;$!ba;s/\n/ /g'

"W_LATC1W_12_219",  "W_ALAS_00W_228",  "W_LATC5W_18_190",

#C_List
sed 's/^C..........._/"C_/g' ./C_list.txt | \
sed 's/^C.........._/"C_/g' | \
sed 's/^C........_/"C_/g' | \
sed 's/^C......._/"C_/g' | \
sed 's/_R1_val_1_bismark_bt2_pe.deduplicated.bedGraph_sorted_tab.bedGraph/", /g' | \
sed ':a;N;$!ba;s/\n/ /g'




==> SE_metilene_W_C.input <==

chrom   pos     W_SE_L_W_L_198_34_F112563       W_SE_T_W_H_228_95_F112570       W_SE_L_W_L_190_43_F112556W_SE_L_W_H_219_104_F112567      W_SE_L_W_H_219_105J_F112571     W_SE_L_W_H_219_105_F112564      W_SE_L_W_L_219_42_F112557        C_SE_L_C_L_198_35_F112558       C_SE_L_C_H_194_92_F112566       C_SE_A_C_L_170_5_F112555 C_SE_L_C_H_187_98_F112565       C_SE_L_C_H_191_102_F112568      C_SE_T_C_L_238_47_F112559       C_SE_L_C_H_194_92J_F112572       C_SE_L_C_H_191_100_F112569      C_SE_A_C_L_125_11_F112561       C_SE_L_C_L_187_38_F112560 |
echo "W_SE_L_W_L_198_34_F112563       W_SE_T_W_H_228_95_F112570       W_SE_L_W_L_190_43_F112556W_SE_L_W_H_219_104_F112567      W_SE_L_W_H_219_105J_F112571     W_SE_L_W_H_219_105_F112564      W_SE_L_W_L_219_42_F112557        C_SE_L_C_L_198_35_F112558       C_SE_L_C_H_194_92_F112566       C_SE_A_C_L_170_5_F112555 C_SE_L_C_H_187_98_F112565       C_SE_L_C_H_191_102_F112568      C_SE_T_C_L_238_47_F112559       C_SE_L_C_H_194_92J_F112572       C_SE_L_C_H_191_100_F112569      C_SE_A_C_L_125_11_F112561       C_SE_L_C_L_187_38_F112560" | sed -e 's/ \+/","/g' -e 's/^/"/' -e 's/$/"/' -e 's/,/","/g'



#!/bin/bash

sed -e 's/ \+/","/g' -e 's/^/"/' -e 's/$/"/' -e 's/,/","/g'



#Combine and remove 
"W_LATC1W_12_219",  "W_ALAS_00W_228",  "W_LATC5W_18_190", "C_FERT31C_15F_170",  "C_WILL7C_445_125",  "C_LATD2C_6_198",  "C_LATD5C_5_191",  "C_LATJ_02C_194",  "C_LATJ_00C_187",  "C_ALAS0C_5_238", "W_SE_L_W_L_198_34_F112563", "W_SE_T_W_H_228_95_F112570", "W_SE_L_W_L_190_43_F112556", "W_SE_L_W_H_219_104_F112567", "W_SE_L_W_H_219_105J_F112571", "W_SE_L_W_H_219_105_F112564", "W_SE_L_W_L_219_42_F112557, "C_SE_L_C_L_198_35_F112558", "C_SE_L_C_H_194_92_F112566", "C_SE_A_C_L_170_5_F112555", "C_SE_L_C_H_187_98_F112565", "C_SE_L_C_H_191_102_F112568, "C_SE_T_C_L_238_47_F112559", "C_SE_L_C_H_194_92J_F112572", "C_SE_L_C_H_191_100_F112569", "C_SE_A_C_L_125_11_F112561", "C_SE_L_C_L_187_38_F112560"



##############################################3
module load StdEnv/2020
module load r/4.2.1
module load gdal
module load udunits
module load python
export R_LIBS_USER=/home/msandler/R/x86_64-pc-linux-gnu-library/4.2.1/
#####################
#Testing for:Do1_01_a00004_7719865_7720140.txt 
#convert to script later

chrom="Do1_01_a00004"
start=7719865
end=7720140
directory="/home/msandler/scratch/Barplots"

setwd ("/home/msandler/scratch/Barplots")


library(ggplot2)
library(ggpubr)

#import data
data <- read.table(paste0(directory,"/data/",chrom,"_",start,"_",end,".txt"), header = FALSE)

#parents

colnames(data) <- c("chrom", "pos", "W_LATC1W_12_219", "W_ALAS_00W_228", "W_LATC5W_18_190", "C_FERT31C_15F_170", "C_WILL7C_445_125", "C_LATD2C_6_198", "C_LATD5C_5_191", "C_LATJ_02C_194", "C_LATJ_00C_187", "C_ALAS0C_5_238", "W_SE_L_W_L_198_34_F112563", "W_SE_T_W_H_228_95_F112570", "W_SE_L_W_L_190_43_F112556", "W_SE_L_W_H_219_104_F112567", "W_SE_L_W_H_219_105J_F112571", "W_SE_L_W_H_219_105_F112564", "W_SE_L_W_L_219_42_F112557", "C_SE_L_C_L_198_35_F112558", "C_SE_L_C_H_194_92_F112566", "C_SE_A_C_L_170_5_F112555", "C_SE_L_C_H_187_98_F112565", "C_SE_L_C_H_191_102_F112568", "C_SE_T_C_L_238_47_F112559", "C_SE_L_C_H_194_92J_F112572", "C_SE_L_C_H_191_100_F112569", "C_SE_A_C_L_125_11_F112561", "C_SE_L_C_L_187_38_F112560")
#> ncol(data)
#[1] 29
#> nrow(data)
#[1] 54


datat <-  as.data.frame(t(as.matrix(data)))


#> ncol(datat)
#[1] 54
#> nrow(datat)
#[1] 29


datat$ID <- c("chrom", "pos", "W_LATC1W_12_219", "W_ALAS_00W_228", "W_LATC5W_18_190", "C_FERT31C_15F_170", "C_WILL7C_445_125", "C_LATD2C_6_198", "C_LATD5C_5_191", "C_LATJ_02C_194", "C_LATJ_00C_187", "C_ALAS0C_5_238", "W_SE_L_W_L_198_34_F112563", "W_SE_T_W_H_228_95_F112570", "W_SE_L_W_L_190_43_F112556", "W_SE_L_W_H_219_104_F112567", "W_SE_L_W_H_219_105J_F112571", "W_SE_L_W_H_219_105_F112564", "W_SE_L_W_L_219_42_F112557", "C_SE_L_C_L_198_35_F112558", "C_SE_L_C_H_194_92_F112566", "C_SE_A_C_L_170_5_F112555", "C_SE_L_C_H_187_98_F112565", "C_SE_L_C_H_191_102_F112568", "C_SE_T_C_L_238_47_F112559", "C_SE_L_C_H_194_92J_F112572", "C_SE_L_C_H_191_100_F112569", "C_SE_A_C_L_125_11_F112561", "C_SE_L_C_L_187_38_F112560")
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

#install.packages('stringr')
library(stringr)
p_se <- c()
treat <- c()

for (id in mydata$ID) {
  if (str_detect(id, "SE")) {
    p_se <- c(p_se, "SE")
  } else  {
    p_se <- c(p_se, "P")
  }
  if (str_starts(id, "W")) {
    treat <- c(treat, "W")
  } else  {
    treat <- c(treat, "C")
  }
}

print(p_se)
print(treat)


mydata <- cbind(mydata, treat)
mydata <- cbind(mydata, p_se)
colnames(mydata) <- c("ID", "TotalMeth","Treat", "Seedling/Parent")

#q <- as.data.frame(t(as.matrix(as.data.frame((strsplit(as.character(mydata$ID), "_"))))))
#split so that q

#
mydata <- cbind(mydata, q)

##Works until here
#If SE -> add column that says SE else add column that says P




# plots
mydata$TotalMeth <- as.numeric(as.character(mydata$TotalMeth))

#make Plot unique to Site
mydata[(ncol(mydata)+1)] <- paste(mydata$SE_P,mydata$Treat, sep=".")
names(mydata)[(ncol(mydata))] <- "SE_P_Plot"
unique(mydata$SE_P_Plot)
length(unique(mydata$SE_P_Plot))


mydata$SE_P_Plot <- as.factor(mydata$SE_P_Plot)
#mydata$RegPlot <- as.factor(mydata$RegPlot)
mydata$Treat <- as.factor(mydata$Treat)

#colour
colour <- c("blue", "red","blue", "red","blue", "red","blue", "red")
colour_reg <- c("blue", "red","blue", "red","blue", "red","blue", "red")


jpeg(paste0(directory,"/plots/SE_P_Plot_",chrom,"_",start,"_",end,".jpg"), width = 1700, height = 1500)
par(mgp=c(3,3,0))
par(mar = c(35, 20, 4.1, 2.1))
plot(TotalMeth~MixPlot, data = mydata, cex=1.5, ylim =c(min(mydata$TotalMeth, na.rm=TRUE),max(mydata$TotalMeth, na.rm=TRUE)), xlab = "", ylab="", col=colour, las=2, cex.lab=5, cex.axis=5)
stripchart(TotalMeth~SE_P_Plot, vertical = TRUE, method = "jitter", pch = 16, cex = 3,
           col = "black", data = mydata, add=TRUE)
mtext(side=1, text="", line=30, cex=5)
mtext(side=2, text=paste0(chrom,":",start,"-",end," Methylation"), line=15, cex=5)
dev.off()

jpeg(paste0(directory,"/plots/RegPlot_",chrom,"_",start,"_",end,".jpg"), width = 1700, height = 1500)
par(mgp=c(3,3,0))
par(mar = c(15, 20, 4.1, 2.1))
plot(TotalMeth~RegPlot, data = mydata, cex=1.5, ylim =c(min(mydata$TotalMeth, na.rm=TRUE),max(mydata$TotalMeth, na.rm=TRUE)), xlab = "", ylab="", col=colour_reg, las=2, cex.lab=5, cex.axis=5)
stripchart(TotalMeth~RegPlot, vertical = TRUE, method = "jitter", pch = 16, cex = 3,
           col = "black", data = mydata, add=TRUE)
mtext(side=1, text="Site", line=10, cex=5)
mtext(side=2, text=paste0(chrom,":",start,"-",end," Methylation"), line=15, cex=5)
dev.off()

jpeg(paste0(directory,"/plots/TreatPlot_",chrom,"_",start,"_",end,".jpg"), width = 1700, height = 1500)
par(mgp=c(3,3,0))
par(mar = c(15, 20, 4.1, 2.1))
plot(TotalMeth~Treat, data = mydata, cex=1.5, ylim =c(min(mydata$TotalMeth, na.rm=TRUE),max(mydata$TotalMeth, na.rm=TRUE)), xlab = "", ylab="", col=colour_reg, las=2, cex.lab=5, cex.axis=5)
stripchart(TotalMeth~Treat, vertical = TRUE, method = "jitter", pch = 16, cex = 3,
           col = "black", data = mydata, add=TRUE)
mtext(side=1, text="Warming Treatment", line=10, cex=5)
mtext(side=2, text=paste0(chrom,":",start,"-",end," Methylation"), line=15, cex=5)
dev.off()

jpeg(paste0(directory,"/plots/SE_vs_P_Plot_",chrom,"_",start,"_",end,".jpg"), width = 1700, height = 1500)
par(mgp=c(3,3,0))
par(mar = c(15, 20, 4.1, 2.1))
plot(TotalMeth~Treat, data = mydata, cex=1.5, ylim =c(min(mydata$TotalMeth, na.rm=TRUE),max(mydata$TotalMeth, na.rm=TRUE)), xlab = "", ylab="", col=colour_reg, las=2, cex.lab=5, cex.axis=5)
stripchart(TotalMeth~Treat, vertical = TRUE, method = "jitter", pch = 16, cex = 3,
           col = "black", data = mydata, add=TRUE)
mtext(side=1, text="Warming Treatment", line=10, cex=5)
mtext(side=2, text=paste0(chrom,":",start,"-",end," Methylation"), line=15, cex=5)
dev.off()

#################

scp -v msandler@cedar.computecanada.ca:/home/msandler/scratch/Barplots/plots/* .