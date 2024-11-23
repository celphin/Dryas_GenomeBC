

# # DESeq and EdgeR are very similar and both assume that no genes are differentially expressed. 
# DESeq uses a "geometric" normalisation strategy, whereas EdgeR is a weighted mean of log ratios-based
 # method. Both normalise data initially via the calculation of size / normalisation factors.

# # Limma / Voom is different in that it normalises via the very successful (for microarrays) 
# quantile nomalisation, where an attempt is made to match gene count distributions across 
# samples in your dataset. It can somewhat loosely be viewed as scaling each sample's values 
# to be between the min and max values (across all samples). Thus, the final distributions will be similar.

# # Note added August 19, 2020: for limma, if we are referring to microarrays and not RNA-seq, 
# then normalisation will be performed by affy or oligo for Affymetrix arrays, while limma has 
# functionality to normalise Illumina and Agilent arrays.

# # Here is further information (important parts in bold):
# # DESeq2

    # # DESeq: This normalization method [14] is included in the DESeq Bioconductor package 
	# (version 1.6.0) [14] and is based on the hypothesis that most genes are not DE. A DESeq 
	# scaling factor for a given lane is computed as the median of the ratio, for each gene, 
	# of its read count over its geometric mean across all lanes. The underlying idea is that 
	# non-DE genes should have similar read counts across samples, leading to a ratio of 1. 
	# Assuming most genes are not DE, the median of this ratio for the lane provides an estimate
	# of the correction factor that should be applied to all read counts of this lane to fulfill 
	# the hypothesis. By calling the estimateSizeFactors() and sizeFactors() functions in the DESeq 
	# Bioconductor package, this factor is computed for each lane, and raw read counts are divided 
	# by the factor associated with their sequencing lane.

# # [source: https://www.ncbi.nlm.nih.gov/pubmed/22988256]
# # EdgeR

    # # Trimmed Mean of M-values (TMM): This normalization method [17] is implemented in the 
	# edgeR Bioconductor package (version 2.4.0). It is also based on the hypothesis that most
	# genes are not DE. The TMM factor is computed for each lane, with one lane being considered 
	# as a reference sample and the others as test samples. For each test sample, TMM is computed 
	# as the weighted mean of log ratios between this test and the reference, after exclusion of 
	# the most expressed genes and the genes with the largest log ratios. According to the hypothesis
	# of low DE, this TMM should be close to 1. If it is not, its value provides an estimate of the 
	# correction factor that must be applied to the library sizes (and not the raw counts) in order 
	# to fulfill the hypothesis. The calcNormFactors() function in the edgeR Bioconductor package
	# provides these scaling factors. To obtain normalized read counts, these normalization factors
	# are re-scaled by the mean of the normalized library sizes. Normalized read counts are obtained 
	# by dividing raw read counts by these re-scaled normalization factors.

# # [source: https://www.ncbi.nlm.nih.gov/pubmed/22988256]
# # Limma / Voom

    # Quantile (Q): First proposed in the context of microarray data, this normalization method 
	# consists in matching distributions of gene counts across lanes [22, 23]. It is implemented in 
	# the Bioconductor package limma [31] by calling the normalizeQuantiles() function.

# # [source: https://www.ncbi.nlm.nih.gov/pubmed/22988256]


#####################################################

rm(list=ls())

## Commented out because I've already installed the program

#source("http://bioconductor.org/biocLite.R")
#BiocManager::install("DESeq2")

library("DESeq2")

# Set working directory - This will be different on your machine
directory <- "~/Desktop/read_counts/"
setwd(directory)

# Set the prefix for each output file name
outputPrefix <- "RNA_tutorial"

# List the names for each of the count files
sampleFiles<- c("warm_sample_01.read_counts.c.txt",
                "warm_sample_03.read_counts.c.txt",
                "warm_sample_02.read_counts.c.txt",
                "warm_sample_04.read_counts.c.txt",
                "warm_sample_05.read_counts.c.txt",
                "warm_sample_06.read_counts.c.txt",
                "cold_sample_07.read_counts.c.txt",
                "cold_sample_08.read_counts.c.txt",
                "cold_sample_09.read_counts.c.txt",
                "cold_sample_10.read_counts.c.txt",
                "cold_sample_11.read_counts.c.txt",
                "cold_sample_12.read_counts.c.txt")

# List the sample names
sampleNames<- c("sample_01",
                "sample_02",
                "sample_03",
                "sample_04",
                "sample_05",
                "sample_06",
                "sample_07",
                "sample_08",
                "sample_09",
                "sample_10",
                "sample_11",
                "sample_12")

# A list of treatments
sampleCondition <- c(rep("warm", 6),rep("cold", 6))

# Combine the sample info into a table for DESeq
sampleTable <- data.frame(sampleName = sampleNames, fileName = sampleFiles, condition = sampleCondition)

# A list of the factors for a future dataframe
treatments = c("warm","cold")



# Built the DESeq dataset object using the counts data
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = directory,
                                       design = ~ condition)

# Estimate size factors - normalize read counts
dds <- estimateSizeFactors(ddsHTSeq)

# Take a look at the raw counts
View(counts(ddsHTSeq))
# Take a look at the normalised counts
View(counts(dds, normalized=TRUE))

##########################################

# Convert treatment into factors
colData(ddsHTSeq)$condition <- factor(colData(ddsHTSeq)$condition,
                                      levels = treatments)


# The guts of the analysis - run the DESeq2 GLM on the count data
dds <- DESeq(ddsHTSeq)

res <- results(dds)
# Make a simple volcano plot
plot(res$log2FoldChange, -log10(res$padj))

library(ggplot2)

res_DF <- data.frame(res)

# A quick Volcano plot
ggplot(data = res_DF, aes( y = -log10(padj), x = log2FoldChange, fill = abs(log2FoldChange)))+
  geom_point(shape = 16)+
  geom_point(shape = 21)+
  scale_fill_gradient(low="white",high = "red")+
  theme_bw()

# A quick Volcano plot - log10 the y-axis
ggplot(data = res_DF, aes( y = -log10(padj), x = log2FoldChange, fill = abs(log2FoldChange)))+
  geom_point(shape = 16)+
  geom_point(shape = 21)+
  scale_fill_gradient(low="white",high = "red")+
  scale_y_log10()+
  theme_bw()



# transform raw counts into normalized values
# DESeq2 has two options:  1) rlog transformed and 2) variance stabilization
# variance stabilization is very good for heatmaps, etc.
rld <- rlogTransformation(dds, blind=T)
vsd <- varianceStabilizingTransformation(dds, blind=T)

#Principal components plot shows additional but rough clustering of samples
library("genefilter")
library("ggplot2")
library("grDevices")

rv <- rowVars(assay(rld))
select <- order(rv, decreasing=T)[seq_len(min(500,length(rv)))]
pc <- prcomp(t(assay(vsd)[select,]))


# set condition
condition <- treatments
scores <- data.frame(pc$x)
scores$condition <- rep(c("warm","cold"),each =6)
(pcaplot <- ggplot(scores, aes(x = PC1, y = PC2, col = (factor(condition))))
  + geom_point(size = 5)
  + ggtitle("Principal Components")
  + scale_colour_brewer(name = " ", palette = "Set1")
  + theme(
    plot.title = element_text(face = 'bold'),
    legend.key = element_rect(fill = 'NA'),
    legend.text = element_text(size = 10, face = "bold"),
    axis.text.y = element_text(colour = "Black"),
    axis.text.x = element_text(colour = "Black"),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = 'bold'),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.background = element_rect(color = 'black',fill = NA)
  ))

ggsave(pcaplot,file=paste0(outputPrefix, "-ggplot2.pdf"))

