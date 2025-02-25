########################################
#Plotting inheritence bar plots
#In progress (needs to be cleaned)
########################################
args = commandArgs(trailingOnly=TRUE)

# test if there is at least two arguments: if not, return an error
if (length(args)<3) {
  stop("At least three arguments must be supplied ('chrom', start, end )", call.=FALSE)
}

chrom=args[1]                # read first argument as string
start=as.integer( args[2] )  # read second argument as integer
end=as.integer(args[3])

# chrom="Do1_01_a00004"
# start=7719865
# end=7720140

directory="/home/msandler/scratch/Barplots"
setwd ("/home/msandler/scratch/Barplots")


library(ggplot2)
library(ggpubr)

#import data
P_data <- read.table(paste0(directory,"/data/P_",chrom,"_",start,"_",end,".txt"), header = FALSE)
SE_data <- read.table(paste0(directory,"/data/SE_",chrom,"_",start,"_",end,".txt"), header = FALSE)
#cleaned_Do1_01_a00004_9409994_9410213.txt
#parents

#data <- read.table(paste0(directory,"/data/cleaned_Do1_01_a00004_9409994_9410213.txt"), header = FALSE)


colnames(P_data) <- c("chrom", "pos", "W_LATC1W_12_219", "W_ALAS_00W_228", "W_LATC5W_18_190", "C_FERT31C_15F_170", "C_WILL7C_445_125", "C_LATD2C_6_198", "C_LATD5C_5_191", "C_LATJ_02C_194", "C_LATJ_00C_187", "C_ALAS0C_5_238")
colnames(SE_data) <- c("chrom", "pos", "W_SE_L_W_L_198_34_F112563", "W_SE_T_W_H_228_95_F112570", "W_SE_L_W_L_190_43_F112556", "W_SE_L_W_H_219_104_F112567", "W_SE_L_W_H_219_105J_F112571", "W_SE_L_W_H_219_105_F112564", "W_SE_L_W_L_219_42_F112557", "C_SE_L_C_L_198_35_F112558", "C_SE_L_C_H_194_92_F112566", "C_SE_A_C_L_170_5_F112555", "C_SE_L_C_H_187_98_F112565", "C_SE_L_C_H_191_102_F112568", "C_SE_T_C_L_238_47_F112559", "C_SE_L_C_H_194_92J_F112572", "C_SE_L_C_H_191_100_F112569", "C_SE_A_C_L_125_11_F112561", "C_SE_L_C_L_187_38_F112560")

P_datat <-  as.data.frame(t(as.matrix(P_data)))
SE_datat <-  as.data.frame(t(as.matrix(SE_data)))

P_datat$ID <- c("chrom", "pos", "W_LATC1W_12_219", "W_ALAS_00W_228", "W_LATC5W_18_190", "C_FERT31C_15F_170", "C_WILL7C_445_125", "C_LATD2C_6_198", "C_LATD5C_5_191", "C_LATJ_02C_194", "C_LATJ_00C_187", "C_ALAS0C_5_238")
SE_datat$ID <- c("chrom", "pos", "W_SE_L_W_L_198_34_F112563", "W_SE_T_W_H_228_95_F112570", "W_SE_L_W_L_190_43_F112556", "W_SE_L_W_H_219_104_F112567", "W_SE_L_W_H_219_105J_F112571", "W_SE_L_W_H_219_105_F112564", "W_SE_L_W_L_219_42_F112557", "C_SE_L_C_L_198_35_F112558", "C_SE_L_C_H_194_92_F112566", "C_SE_A_C_L_170_5_F112555", "C_SE_L_C_H_187_98_F112565", "C_SE_L_C_H_191_102_F112568", "C_SE_T_C_L_238_47_F112559", "C_SE_L_C_H_194_92J_F112572", "C_SE_L_C_H_191_100_F112569", "C_SE_A_C_L_125_11_F112561", "C_SE_L_C_L_187_38_F112560")

P_datat <- P_datat[-c(1,2),]
SE_datat <- SE_datat[-c(1,2),]

#Factor to numeric 
X <- length(P_datat)
P_factorToNumeric <- function(f) as.numeric(as.character(f))
P_datat[1:X] <- lapply(P_datat[1:X], P_factorToNumeric)


X <- length(SE_datat)
SE_factorToNumeric <- function(f) as.numeric(as.character(f))
SE_datat[1:X] <- lapply(SE_datat[1:X], SE_factorToNumeric)

#Sum rows
P_datat$total <- rowSums(P_datat[1:X], na.rm=TRUE)
P_datat$ID <- rownames(P_datat)

SE_datat$total <- rowSums(SE_datat[1:X], na.rm=TRUE)
SE_datat$ID <- rownames(SE_datat)

P_data <- as.data.frame(cbind(P_datat$ID ,P_datat$total))
SE_data <- as.data.frame(cbind(SE_datat$ID ,SE_datat$total))

colnames(P_data) <- c("P_ID", "P_TotalMeth")
colnames(SE_data) <- c("SE_ID", "SE_TotalMeth")

#install.packages('stringr')
library(stringr)
p_se <- c()
treat <- c()
randID <-c()


for (id in P_data$P_ID) {
  p_se <- c(p_se, "P")
  randID <- c(randID, str_split(id, "_", simplify =TRUE)[,4])
  if (str_starts(id, "W")) {
    treat <- c(treat, "W")
  } else  {
    treat <- c(treat, "C")
  }
  
}

P_data <- cbind(P_data, treat)
P_data <- cbind(P_data, p_se)
P_data <- cbind(P_data, randID)
colnames(P_data) <- c("Parent_ID", "P_TotalMeth","Treat", "SE_P", "randID")

p_se <- c()
treat <- c()
randID <-c()
highlow <-c()



for (id in SE_data$SE_ID) {
   p_se <- c(p_se, "SE")
   randID <- c(randID, str_split(id, "_", simplify =TRUE)[,6])
    if(str_split(id, "_", simplify =TRUE)[,5]=="H") {
        highlow <- c(highlow, "H")
    } else {
        highlow <- c(highlow, "L")
    }

    if (str_starts(id, "W")) {
        treat <- c(treat, "W")
    } else  {
        treat <- c(treat, "C")
    }
   
}

SE_data <- cbind(SE_data, treat)
SE_data <- cbind(SE_data, p_se)
SE_data <- cbind(SE_data, randID)
SE_data <- cbind(SE_data, highlow)
colnames(SE_data) <- c("Seedling_ID", "SE_TotalMeth","Treat", "SE_P", "randID", "HighLow")



#q <- as.data.frame(t(as.matrix(as.data.frame((strsplit(as.character(mydata$ID), "_"))))))
#split so that q




############################################
#subsetting parent and seedling info
#Left merge by id
se_p_data <- merge(SE_data, P_data, by = "randID", all.x = TRUE)
se_p_data <- subset(se_p_data, select = -Treat.x)
se_p_data <- subset(se_p_data, select = -SE_P.x)
se_p_data <- subset(se_p_data, select = -SE_P.y)
colnames(se_p_data) <- (c("randID", "Seedling_ID","SE_Meth", "HighLow", "Parent_ID", "P_Meth","Treat"))

############################################
se_p_data$SE_Meth <- as.numeric(as.character(se_p_data$SE_Meth))
se_p_data$P_Meth <- as.numeric(as.character(se_p_data$P_Meth))

#############################################
#Scatter plot, Linear
#col=as.factor(se_p_data$Treat)
jpeg(paste0(directory,"/plots/HighLowScatter_",chrom,"_",start,"_",end,".jpg"), width = 3000, height = 1700)
model <- lm(P_Meth ~ SE_Meth, se_p_data)
summary(model)
par(mar=c(20,20,4,4))
plot(se_p_data$P_Meth, se_p_data$SE_Meth, xlab="Parent Methylation", ylab="Seedling Methylation", pch= ifelse(se_p_data$HighLow == "H", 19, 17), cex = 3, mgp=c(10,5,0), cex.lab=3, cex.axis=3, col = ifelse(se_p_data$Treat == "W", "red", "blue"))
abline(model)
dev.off()
##################################################
#Remove 
SE_data1 <- subset(SE_data, select = -HighLow)
colnames(SE_data1) <- c("ID", "TotalMeth","Treat", "SE_P", "randID")
P_data1 <- P_data
colnames(P_data1) <- c("ID", "TotalMeth","Treat", "SE_P", "randID")
mydata <- rbind(P_data1, SE_data1)

#
mydata$TotalMeth <- as.numeric(as.character(mydata$TotalMeth))

#make Plot unique to SE_P
mydata[(ncol(mydata)+1)] <- paste(mydata$SE_P,mydata$Treat, sep=".")
names(mydata)[(ncol(mydata))] <- "SE_P_Plot"
unique(mydata$SE_P_Plot)
length(unique(mydata$SE_P_Plot))




####################

# plots
mydata$TotalMeth <- as.numeric(as.character(mydata$TotalMeth))

#make Plot unique to SE_P
mydata[(ncol(mydata)+1)] <- paste(mydata$SE_P,mydata$Treat, sep=".")
names(mydata)[(ncol(mydata))] <- "SE_P_BarPlot"
unique(mydata$SE_P_Plot)
length(unique(mydata$SE_P_Plot))

mydata$SE_P_Plot <- as.factor(mydata$SE_P_Plot)
mydata$Treat <- as.factor(mydata$Treat)

colour <- c("blue", "red","blue", "red","blue", "red","blue", "red")

jpeg(paste0(directory,"/plots/SE_P_BarPlot_",chrom,"_",start,"_",end,".jpg"), width = 1700, height = 1500)
par(mgp=c(3,3,0))
par(mar = c(35, 20, 4.1, 2.1))
plot(TotalMeth~SE_P_Plot, data = mydata, cex=1.5, ylim =c(min(mydata$TotalMeth, na.rm=TRUE),max(mydata$TotalMeth, na.rm=TRUE)), xlab = "", ylab="", col=colour, las=2, cex.lab=5, cex.axis=5)
stripchart(TotalMeth~SE_P_Plot, vertical = TRUE, method = "jitter", pch = 16, cex = 3,
           col = "black", data = mydata, add=TRUE)
mtext(side=1, text="", line=30, cex=5)
mtext(side=2, text=paste0(chrom,":",start,"-",end," Methylation"), line=15, cex=5)
dev.off()


#########################

jpeg(paste0(directory,"/plots/TreatPlot_",chrom,"_",start,"_",end,".jpg"), width = 1700, height = 1500)
par(mgp=c(3,3,0))
par(mar = c(15, 20, 4.1, 2.1))
plot(TotalMeth~Treat, data = mydata, cex=1.5, ylim =c(min(mydata$TotalMeth, na.rm=TRUE),max(mydata$TotalMeth, na.rm=TRUE)), xlab = "", ylab="", col=colour, las=2, cex.lab=5, cex.axis=5)
stripchart(TotalMeth~Treat, vertical = TRUE, method = "jitter", pch = 16, cex = 3,
           col = "black", data = mydata, add=TRUE)
mtext(side=1, text="Warming Treatment", line=10, cex=5)
mtext(side=2, text=paste0(chrom,":",start,"-",end," Methylation"), line=15, cex=5)
dev.off()