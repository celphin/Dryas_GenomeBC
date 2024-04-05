###########################
# Gene ontology analysis
# Mar 2024
##########################

#######################
# to run on species without annotation

# try cluster profiler
# https://rdrr.io/bioc/clusterProfiler/man/enricher.html

cd /home/celphin/scratch/Dryas/MS_Dryas_Merged_Data

less ./original_data/interproscan_dryas_full.tsv

Do1_04_a00001G02558V1.1 53f5d6e35629feb799400172e4290a35        162     Gene3D  G3DSA:3.30.40.10               Zinc/RING finger domain, C3HC4 (zinc finger)    96      156     1.1E-5  T       29-09-2023             IPR013083       Zinc finger, RING/FYVE/PHD-type -
Do1_04_a00004G00129V1.1 777732e682d73302a04e01445672e633        213     PANTHER PTHR37610       CCHC-TYPE DOMAIN-CONTAINING PROTEIN    4       127     2.8E-12 T       29-09-2023      -       -
Do1_07_a00002G00409V1.1 aaaa9a5a2bb10efad31307eaee179473        224     PANTHER PTHR47895       CYSTEINE AND HISTIDINE-RICH DOMAIN-CONTAINING PROTEIN RAR1     7       223     2.5E-105        T       29-09-2023     IPR043316       Cysteine and histidine-rich domain-containing protein RAR1      -
Do1_07_a00002G00409V1.1 aaaa9a5a2bb10efad31307eaee179473        224     MobiDBLite      mobidb-lite            consensus disorder prediction   82      100     -       T       29-09-2023      -       -
Do1_07_a00002G00409V1.1 aaaa9a5a2bb10efad31307eaee179473        224     ProSiteProfiles PS51401 CHORD domain profile.  11      69      28.896263       T       29-09-2023      IPR007051       CHORD domain           -
Do1_07_a00002G00409V1.1 aaaa9a5a2bb10efad31307eaee179473        224     Gene3D  G3DSA:4.10.1130.20             -       151     221     6.5E-27 T       29-09-2023      -       -
Do1_07_a00002G00409V1.1 aaaa9a5a2bb10efad31307eaee179473        224     Gene3D  G3DSA:4.10.1130.20             -       2       73      8.1E-26 T       29-09-2023      -       -
Do1_07_a00002G00409V1.1 aaaa9a5a2bb10efad31307eaee179473        224     Pfam    PF04968 CHORD   10             69      1.3E-23 T       29-09-2023      IPR007051       CHORD domain    -

#--------------------------
# edit input for enricher

# use excel to open and extract GOTERM, database, NAME and GeneID from ./original_data/interproscan_dryas_full.tsv
# make GOTERM_NAME_dryas_full.tsv.txt

scp GOTERM_NAME_dryas_full.tsv.txt celphin@cedar.computecanada.ca:/home/celphin/scratch/Dryas/MS_Dryas_Merged_Data/

# remove duplicated rows
sort GOTERM_NAME_dryas_full.tsv.txt | uniq > GOTERM_NAME_dryas_uniq.txt

#------------------------------
module load StdEnv/2023
module load r/4.3.1

R

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("clusterProfiler")

#-----------
# add more dependancies
BiocManager::install("BiocParallel")
BiocManager::install("clusterProfiler")

#-------------------------
# read in data
setwd("/home/celphin/scratch/Dryas/MS_Dryas_Merged_Data/")

total_gotermsfile <- read.table("GOTERM_NAME_dryas_uniq.txt", sep = "\t", header=FALSE) 

DMRs <- read.table("Gene_DMR_Total_GO_Merged_table.tsv", sep = "\t", header=TRUE) 
DEGs <- read.table("Gene_RNA_Total_GO_Merged_table.tsv", sep = "\t", header=TRUE) 
DMR_RNA <- read.table("Gene_DMR_RNA_GO_Left_join_table.tsv", sep = "\t", header=TRUE) 

#----------------------------
# make TERM2GENE dataframe
# https://stackoverflow.com/questions/72728057/create-new-column-with-unique-id-to-handle-duplicates-from-2-columns-in-r

total_gotermsfile$V1 <- stringr::str_split(total_gotermsfile$V1, "V", simplify =TRUE)[,1]
TERM2GENE <- total_gotermsfile[,c(3,1)]
TERM2GENE_narm <- TERM2GENE[-which(is.na(TERM2GENE[,1])),]

# make TERM2NAME dataframe
TERM2NAME <- total_gotermsfile[,c(3,2)]
TERM2NAME_narm <- TERM2NAME[-which(is.na(TERM2NAME[,1])),]

# make background_genes_list 
background_genes <- total_gotermsfile[,c(1)]
uniq_background_genes<- unique(background_genes)

# make gene_id_vector
DMR_genes <- DMRs[,c(1,2)]
DEG_genes <- DEGs[,c(1)]
DMR_DEG_genes <- DMRs[,c(1,2)]

unique(DMR_DEG_genes$Origin)

DMR_Wild_W_C <- DMR_genes[which(DMR_genes$Origin=="Wild_W_C"),][,1]
DMR_Wild_Lat_L_H <- DMR_genes[which(DMR_genes$Origin=="Wild_Lat_L_H"),][,1]
DMR_Mat_Sen <- DMR_genes[which(DMR_genes$Origin=="Mat_Sen"),][,1]

DMR_DEG_Wild_W_C <- DMR_DEG_genes[which(DMR_DEG_genes$Origin=="Wild_W_C"),][,1]

#-------------------------
library(clusterProfiler)

# https://rdrr.io/bioc/clusterProfiler/man/enricher.html

DEG_genes_enrich <- enricher(
  gene=DEG_genes,
  pvalueCutoff = 1,
  universe=uniq_background_genes,
  TERM2GENE=TERM2GENE_narm,
  TERM2NAME=TERM2NAME_narm
)
#...0 enriched terms found
[1] "Glycoside_hydrolase_family_17"
[2] "Sreceptorlike_serine/threonineprotein_kinase"
[3] "Ascorbate_oxidase_third_cupredoxin_domain"
[4] "UDPglucuronosyl/UDPglucosyltransferase"
[5] "START_domain"
[6] "Phospholipid/glycerol_acyltransferase"


DMR_DEG_Wild_W_C_enrich <- enricher(
  gene=DMR_DEG_Wild_W_C,
    pvalueCutoff = 1,
  universe=uniq_background_genes,
  TERM2GENE=TERM2GENE_narm,
  TERM2NAME=TERM2NAME_narm
)

#...@organism    UNKNOWN
#...@ontology    UNKNOWN
#...@gene        chr [1:119] "Do1_00153G00006" "Do1_00153G00016" "Do1_00153G00020" ...
#...pvalues adjusted by 'BH' with cutoff <0.05
#...2 enriched terms found
'data.frame':   2 obs. of  9 variables:
 $ ID         : chr  "GO:0006486|GO:0016757" "GO:0008194"
 $ Description: chr  "Exostosinlike" "UDPglucuronosyl/UDPglucosyltransferase"
 $ GeneRatio  : chr  "5/23" "4/23"
 $ BgRatio    : chr  "21/8057" "57/8057"
 $ pvalue     : num  2.35e-09 1.80e-05
 $ p.adjust   : num  3.76e-08 1.44e-04
 $ qvalue     : num  2.23e-08 8.54e-05
 $ geneID     : chr  "Do1_01_a00001G02423/Do1_01_a00001G02424/Do1_01_a00001G02427/Do1_01_a00001G02428/Do1_01_a00001G02430" "Do1_01_a00001G02221/Do1_01_a00001G02223/Do1_06_a00002G00975/Do1_06_a00002G00976"
 $ Count      : int  5 4

[1] "Exostosinlike"
 [2] "UDPglucuronosyl/UDPglucosyltransferase"
 [3] "Integrase_catalytic_core"
 [4] "Serine/threonineprotein_kinase_active_site"
 [5] "Laccase"
 [6] "Receptorlike_protein_kinase_ANXUR1like"
 [7] "Protein_kinase_domain"
 [8] "Mrp_conserved_site"
 [9] "Multicopper_oxidase_Cterminal"
[10] "Ascorbate_oxidase_third_cupredoxin_domain"
[11] "Zinc_finger_CCHCtype_superfamily"
[12] "Cytochrome_P450_conserved_site"
[13] "FHY3/FAR1_family"
[14] "Cytochrome_P450_superfamily"



DMR_Wild_W_C_enrich<- enricher(
  gene=DMR_Wild_W_C,
  pvalueCutoff = 1,
  universe=uniq_background_genes,
  TERM2GENE=TERM2GENE_narm,
  TERM2NAME=TERM2NAME_narm
)

'data.frame':   2 obs. of  9 variables:
 $ ID         : chr  "GO:0006486|GO:0016757" "GO:0008194"
 $ Description: chr  "Exostosinlike" "UDPglucuronosyl/UDPglucosyltransferase"
 $ GeneRatio  : chr  "5/23" "4/23"
 $ BgRatio    : chr  "21/8057" "57/8057"
 $ pvalue     : num  2.35e-09 1.80e-05
 $ p.adjust   : num  3.76e-08 1.44e-04
 $ qvalue     : num  2.23e-08 8.54e-05
 $ geneID     : chr  "Do1_01_a00001G02423/Do1_01_a00001G02424/Do1_01_a00001G02427/Do1_01_a00001G02428/Do1_01_a00001G02430" "Do1_01_a00001G02221/Do1_01_a00001G02223/Do1_06_a00002G00975/Do1_06_a00002G00976"
 $ Count      : int  5 4

enricher(
  gene=DMR_Wild_Lat_L_H,
  universe=uniq_background_genes,
  TERM2GENE=TERM2GENE_narm,
  TERM2NAME=TERM2NAME_narm
)

#...@organism    UNKNOWN
#...@ontology    UNKNOWN
#...@gene        chr [1:2746] "Do1_00107G00001" "Do1_00107G00009" "Do1_00107G00010" ...
#...pvalues adjusted by 'BH' with cutoff <0.05
#...10 enriched terms found
'data.frame':   10 obs. of  9 variables:
 $ ID         : chr  "GO:0016491" "GO:0008270|GO:0016491" "GO:0016616" "GO:0030246" ...
 $ Description: chr  "Redoxin" "Alcohol_dehydrogenase_zinctype_conserved_site" "Lactate/malate_dehydrogenase_Cterminal" "Legume_lectin_domain" ...
 $ GeneRatio  : chr  "28/467" "9/467" "10/467" "12/467" ...
 $ BgRatio    : chr  "137/8057" "16/8057" "24/8057" "47/8057" ...
 $ pvalue     : num  3.26e-09 5.44e-08 3.65e-07 1.01e-05 1.32e-05 ...
 $ p.adjust   : num  3.74e-07 3.13e-06 1.40e-05 2.91e-04 3.04e-04 ...
 $ qvalue     : num  3.32e-07 2.78e-06 1.24e-05 2.59e-04 2.70e-04 ...
 $ geneID     : chr  "Do1_01_a00001G00580/Do1_01_a00001G00585/Do1_01_a00001G02089/Do1_01_a00001G02136/Do1_01_a00001G02144/Do1_01_a000"| __truncated__ "Do1_01_a00001G02197/Do1_01_a00001G02199/Do1_01_a00001G02202/Do1_01_a00001G02203/Do1_01_a00001G02204/Do1_01_a000"| __truncated__ "Do1_01_a00001G02136/Do1_01_a00001G02197/Do1_01_a00001G02199/Do1_01_a00001G02202/Do1_01_a00001G02203/Do1_01_a000"| __truncated__ "Do1_01_a00001G02187/Do1_01_a00001G02196/Do1_02_a00001G00555/Do1_02_a00001G00556/Do1_02_a00001G00557/Do1_02_a000"| __truncated__ ...
 $ Count      : int  28 9 10 12 6 8 6 8 19 6

enricher(
  gene=DMR_Mat_Sen,
  universe=uniq_background_genes,
  TERM2GENE=TERM2GENE_narm,
  TERM2NAME=TERM2NAME_narm
)

#...@organism    UNKNOWN
#...@ontology    UNKNOWN
#...@gene        chr [1:834] "Do1_00107G00001" "Do1_00107G00006" "Do1_00107G00040" ...
#...pvalues adjusted by 'BH' with cutoff <0.05
#...6 enriched terms found
'data.frame':   6 obs. of  9 variables:
 $ ID         : chr  "GO:0003676|GO:0004523" "GO:0015074" "GO:0006486|GO:0016757" "GO:0005507|GO:0016491" ...
 $ Description: chr  "Ribonuclease_H_domain" "Integrase_catalytic_core" "Exostosinlike" "Multicopper_oxidase_Cterminal" ...
 $ GeneRatio  : chr  "16/103" "18/103" "4/103" "4/103" ...
 $ BgRatio    : chr  "165/8057" "320/8057" "21/8057" "23/8057" ...
 $ pvalue     : num  2.22e-10 9.55e-08 1.28e-04 1.85e-04 3.11e-04 ...
 $ p.adjust   : num  9.54e-09 2.05e-06 1.83e-03 1.99e-03 2.67e-03 ...
 $ qvalue     : num  8.41e-09 1.81e-06 1.61e-03 1.75e-03 2.36e-03 ...
 $ geneID     : chr  "Do1_01_a00001G02008/Do1_01_a00001G02024/Do1_01_a00001G02033/Do1_01_a00001G02209/Do1_01_a00001G02215/Do1_01_a000"| __truncated__ "Do1_01_a00001G02024/Do1_01_a00001G02141/Do1_01_a00001G02209/Do1_01_a00001G02269/Do1_01_a00001G02438/Do1_01_a000"| __truncated__ "Do1_01_a00001G02423/Do1_01_a00001G02424/Do1_01_a00001G02427/Do1_01_a00001G02428" "Do1_02_a00004G00437/Do1_02_a00004G00438/Do1_02_a00004G00440/Do1_02_a00004G00509" ...
 $ Count      : int  16 18 4 4 3 4


#----------------------
# https://yulab-smu.top/biomedical-knowledge-mining-book/index.html

enrichplot

    barplot
    cnetplot
    dotplot
    emapplot
    gseaplot
    goplot
    upsetplot




#########################
# other options

# https://bioconductor.org/packages/release/bioc/vignettes/AnnotationForge/inst/doc/MakingNewOrganismPackages.html
# https://jokergoo.github.io/GSEAtraining/articles/topic3_02_local_GREAT.html#other-organisms-who-do-not-have-txdb-objects


#-----------------------
# rGREAT 
# installation

module load StdEnv/2023
module load r/4.3.1

R

if(!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
BiocManager::install("rGREAT")

#----------------------------
# to run online 

set.seed(123)
gr = randomRegions(nr = 1000, genome = "hg19")

job = submitGreatJob(gr)
tbl = getEnrichmentTables(job)

#---------------
# to run on any species
# https://useast.ensembl.org/info/data/biomart/how_to_use_biomart.html

res = great(gr, "MSigDB:H", "TxDb.Hsapiens.UCSC.hg19.knownGene")
tb = getEnrichmentTable(res)

# giant panda
great(gr, "GO:BP", biomart_dataset = "amelanoleuca_gene_ensembl")
