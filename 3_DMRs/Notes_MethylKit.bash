###########################
# methylKit
# https://github.com/al2na/methylKit
# Mar 2024
##############################

module load StdEnv/2023
module load r/4.3.1

R

library(devtools)
install_github("al2na/methylKit", build_vignettes=FALSE, 
  repos=BiocManager::repositories(),
  dependencies=TRUE)
  

#-----------------------------
# to run and explore clustering and DMRs
# https://bioconductor.org/packages/release/bioc/vignettes/methylKit/inst/doc/methylKit.html#1_Introduction

# input data
cd /home/celphin/projects/rpp-rieseber/celphin/Dryas/Methylation_calling/Aug2022_Seedlings/bismark_methylation_calls/methylation_coverage
more Asian1_F112573_R1_val_1_bismark_bt2_pe.deduplicated.bismark.cov
ID             start  End      
Do1_00107       9       9       0       0       1
Do1_00107       10      10      0       0       2
Do1_00107       138     138     5.26315789473684        1       18
Do1_00107       139     139     5       1       19
Do1_00107       248     248     0       0       20

