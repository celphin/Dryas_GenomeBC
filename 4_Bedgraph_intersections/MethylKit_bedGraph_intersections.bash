###########################
# methylKit bedGraph intersections
# https://github.com/al2na/methylKit
# Jan 2025
##############################

# back in bash compare with metilene 

############################################
# intersection of sites

# collect data
cd /lustre04/scratch/celphin/Dryas_large_folders/
mkdir intersections
cd /lustre04/scratch/celphin/Dryas_large_folders/intersections

# CpG
cp /lustre04/scratch/celphin/Dryas_large_folders/CpG/Wild_stranded_CpG_report/*DMRs.txt .
cp /lustre04/scratch/celphin/Dryas_large_folders/CpG/stranded_CpG_report/Seedling/*DMRs.txt .
cp /lustre04/scratch/celphin/Dryas_large_folders/CpG/stranded_CpG_report/Phenology/*DMRs.txt .

# CHH
cp /lustre04/scratch/celphin/Dryas_large_folders/CHG_CHH/data/coverage/CHH_cytosine_report/Methylkit_*.txt .
cp /lustre04/scratch/celphin/Dryas_large_folders/CHG_CHH/seedling_data/coverage/seedling/CHH_cytosine_report/*_DMRs_CHH.txt .
cp /lustre04/scratch/celphin/Dryas_large_folders/CHG_CHH/seedling_data/coverage/phenology/CHH_cytosine_report/*_DMRs_CHH.txt .

cp /lustre04/scratch/celphin/Dryas_large_folders/CHG_CHH/seedling_data/coverage/phenology/CHH_cytosine_report/Methylkit_Pheno_10_DMRs_CHH.txt .
cp /lustre04/scratch/celphin/Dryas_large_folders/CHG_CHH/seedling_data/coverage/phenology/CHH_cytosine_report/Methylkit_Pheno_25_DMRs_CHH.txt .

# random
# CpG
cp /lustre04/scratch/celphin/Dryas_large_folders/CpG/Wild_stranded_CpG_report/*DMRs_rand.txt .
cp /lustre04/scratch/celphin/Dryas_large_folders/CpG/stranded_CpG_report/Seedling/*DMRs_rand.txt .
cp /lustre04/scratch/celphin/Dryas_large_folders/CpG/stranded_CpG_report/Phenology/*DMRs_rand.txt .

# CHH
cp /lustre04/scratch/celphin/Dryas_large_folders/CHG_CHH/data/coverage/CHH_cytosine_report/Methylkit_*.txt .
cp /lustre04/scratch/celphin/Dryas_large_folders/CHG_CHH/seedling_data/coverage/seedling/CHH_cytosine_report/*_DMRs_CHH_rand.txt .
cp /lustre04/scratch/celphin/Dryas_large_folders/CHG_CHH/seedling_data/coverage/phenology/CHH_cytosine_report/*_DMRs_CHH_rand.txt .

# Metilene
mkdir metilene; cd metilene
cp /lustre04/scratch/celphin/Dryas/DMR_bedGraph_files/*.bedGraph .

# random
cd /lustre04/scratch/celphin/Dryas_large_folders/CpG/Metilene/MS_Wild_metilene_output_bedgraphs/Random_Site_bedgraphs
cp *_random_W_C.bedGraph /lustre04/scratch/celphin/Dryas_large_folders/intersections

mkdir metilene_rand; cd metilene_rand
cp /lustre04/scratch/celphin/Dryas/snpEff/Dryas_DMRs/*bedGraph .

# Overdispersion
cd /lustre04/scratch/celphin/Dryas_large_folders/intersections
cp /lustre04/scratch/celphin/Dryas_large_folders/CpG/Wild_stranded_CpG_report/*_overdisp_DMRs.txt .

#----------------------
# make methylkit text output into bedGraph file

for file in *.txt; do
    base_name=$(basename "$file" .txt)
    awk '{print $1, $2, $3, $NF}' "$file" > "${base_name}.bedGraph"
    sed -i 's/ /\t/g' "${base_name}.bedGraph"
    sed -i 's/chr\tstart\tend\tmeth.diff/track type=bedGraph/g' "${base_name}.bedGraph"
done

for file in Methylkit_Pheno_*_DMRs_CHH.txt; do
    base_name=$(basename "$file" .txt)
    awk '{print $1, $2, $3, $NF}' "$file" > "${base_name}.bedGraph"
    sed -i 's/ /\t/g' "${base_name}.bedGraph"
    sed -i 's/chr\tstart\tend\tmeth.diff/track type=bedGraph/g' "${base_name}.bedGraph"
done

for file in *_overdisp_DMRs.txt; do
    base_name=$(basename "$file" .txt)
    awk '{print $1, $2, $3, $NF}' "$file" > "${base_name}.bedGraph"
    sed -i 's/ /\t/g' "${base_name}.bedGraph"
    sed -i 's/chr\tstart\tend\tmeth.diff/track type=bedGraph/g' "${base_name}.bedGraph"
done

#----------------------------
# intersect with metilene output

tmux new-session -s snpEff
tmux attach-session -t snpEff

module load StdEnv/2023 bedtools/2.31.0

# warming data LAT/ALAS CpG
bedtools intersect -u -a Methylkit_ALAS_W_C_10plowcovDMRs.bedGraph \
-b Methylkit_LAT_W_C_10plowcovDMRs.bedGraph > intersect_ALAS_LAT_W_C_10plowcovDMRs.bedGraph
#2 ~1% or less

bedtools intersect -u -a Methylkit_SVAL_W_C_10plowcovDMRs.bedGraph \
-b Methylkit_LAT_W_C_10plowcovDMRs.bedGraph > intersect_SVAL_LAT_W_C_10plowcovDMRs.bedGraph
#0

bedtools intersect -u -a Methylkit_ALAS_W_C_10plowcovDMRs.bedGraph \
-b Methylkit_CASS_W_C_10plowcovDMRs.bedGraph > intersect_ALAS_CASS_W_C_10plowcovDMRs.bedGraph
#0

# compare High/Low and warming - CpG
bedtools intersect -u -a Methylkit_ALAS_W_C_10plowcovDMRs.bedGraph \
-b Methylkit_HL_10plowcovDMRs.bedGraph > intersect_ALAS_HL_10plowcovDMRs.bedGraph
# 28/10 ~19%
bedtools intersect -u -a Methylkit_HL_10plowcovDMRs.bedGraph \
-b Methylkit_LAT_W_C_10plowcovDMRs.bedGraph > intersect_HL_LAT_10plowcovDMRs.bedGraph
#152/900 ~17%

# compare phenology and warming - CpG
bedtools intersect -u -a Methylkit_ALAS_W_C_10plowcovDMRs.bedGraph \
-b Methylkit_Pheno_10_DMRs.bedGraph > intersect_ALAS_Pheno_10plowcovDMRs.bedGraph
# 6/100 ~6%
bedtools intersect -u -a Methylkit_Pheno_10_DMRs.bedGraph \
-b Methylkit_LAT_W_C_10plowcovDMRs.bedGraph > intersect_Pheno_LAT_10plowcovDMRs.bedGraph
# 76/900 ~8% 

# compare phenology and HL Arctic
bedtools intersect -u -a Methylkit_Pheno_10_DMRs.bedGraph \
-b Methylkit_HL_10plowcovDMRs.bedGraph > intersect_Pheno_HL_10plowcovDMRs.bedGraph
# ~10%

# compare seedling and warming - CpG
bedtools intersect -u -a Methylkit_LAT_W_C_10plowcovDMRs.bedGraph \
-b Methylkit_SE_W_C_10_DMRs.bedGraph > intersect_LAT_SE_W_C_10plowcovDMRs.bedGraph
# 34/900 overlap,  ~3%

# compare seedling warming and growth chamber - CpG
bedtools intersect -u -a Methylkit_SE_W_C_10_DMRs.bedGraph \
-b Methylkit_SE_HL_10_DMRs.bedGraph > intersect_SE_WC_HL_10plowcovDMRs.bedGraph

#----
# warming data LAT/ALAS CHH
bedtools intersect -u -a Methylkit_ALAS_rand_W_C_10plowcovDMRs.bedGraph \
-b Methylkit_LAT_rand_W_C_10plowcovDMRs.bedGraph > intersect_rand_ALAS_LAT_W_C_10plowcovDMRs.bedGraph

# compare High/Low and warming - CHH

# compare phenology and warming - CHH

# compare seedling and warming - CHH

# compare seedling warming and growth chamber - CHH

#---
# intersect CpG and CHH

bedtools intersect -u -a Methylkit_ALAS_rand_W_C_10plowcovDMRs.bedGraph \
-b Methylkit_ALAS_W_C_10_CHH_plowcovDMRs.bedGraph > intersect_ALAS_CpG_CHH_W_C_10plowcovDMRs.bedGraph

bedtools intersect -u -a Methylkit_LAT_rand_W_C_10plowcovDMRs.bedGraph \
-b Methylkit_LAT_W_C_10_CHH_plowcovDMRs.bedGraph > intersect_LAT_CpG_CHH_W_C_10plowcovDMRs.bedGraph

bedtools intersect -u -a Methylkit_Pheno_10_DMRs.bedGraph \
-b Methylkit_Pheno_10_DMRs_CHH.bedGraph > intersect_Pheno_CpG_CHH_10plowcovDMRs.bedGraph


#---
# intersect random and non-random

bedtools intersect -u -a Methylkit_LAT_W_C_10plowcovDMRs.bedGraph \
-b Methylkit_LAT_rand_W_C_10plowcovDMRs.bedGraph > intersect_rand-non_LAT_W_C_10plowcovDMRs.bedGraph

bedtools intersect -u -a Methylkit_ALAS_rand_W_C_10plowcovDMRs.bedGraph \
-b Methylkit_ALAS_W_C_10plowcovDMRs.bedGraph > intersect_rand-non_ALAS_W_C_10plowcovDMRs.bedGraph


#------
# random

# ALAS
bedtools intersect -u -a Methylkit_ALAS_rand_W_C_10plowcovDMRs.bedGraph \
-b Methylkit_LAT_rand_W_C_10plowcovDMRs.bedGraph > intersect_ALAS_LAT_rand_W_C_10plowcovDMRs.bedGraph

# LAT
bedtools intersect -u -a Methylkit_ALAS_rand_W_C_10plowcovDMRs.bedGraph \
-b Methylkit_LAT_rand_W_C_10plowcovDMRs.bedGraph > intersect_ALAS_LAT_rand_W_C_10plowcovDMRs.bedGraph

##################################
# three sites and all

# Methylkit10 overdispersion corrected

#Bedtools intersect all 4
module load StdEnv/2023 bedtools/2.31.0
bedtools intersect -u -a Methylkit_CASS_W_C_10plowcovDMRs.bedGraph -b Methylkit_SVAL_W_C_10plowcovDMRs.bedGraph | bedtools intersect -u -a stdin -b Methylkit_LAT_W_C_10plowcovDMRs.bedGraph | bedtools intersect -u -a stdin -b Methylkit_ALAS_W_C_10plowcovDMRs.bedGraph > ALL_Sites_intersect_DMRs_Methylkit10.bedgraph

#Bedtools intersect all but Nunavut
bedtools intersect -u -a Methylkit_SVAL_W_C_10plowcovDMRs.bedGraph -b Methylkit_LAT_W_C_10plowcovDMRs.bedGraph | bedtools intersect -u -a stdin -b Methylkit_ALAS_W_C_10plowcovDMRs.bedGraph > SVAL_SWED_ALAS_Sites_intersect_DMRs_Methylkit10.bedgraph

#Bedtools intersect all but Svalbard
bedtools intersect -u -a Methylkit_CASS_W_C_10plowcovDMRs.bedGraph -b Methylkit_LAT_W_C_10plowcovDMRs.bedGraph | bedtools intersect -u -a stdin -b Methylkit_ALAS_W_C_10plowcovDMRs.bedGraph > ALEX_SWED_ALAS_Sites_intersect_DMRs_Methylkit10.bedgraph

#Bedtools intersect all but Sweden
bedtools intersect -u -a Methylkit_CASS_W_C_10plowcovDMRs.bedGraph -b Methylkit_SVAL_W_C_10plowcovDMRs.bedGraph| bedtools intersect -u -a stdin -b Methylkit_ALAS_W_C_10plowcovDMRs.bedGraph > ALEX_SVAL_ALAS_Sites_intersect_DMRs_Methylkit10.bedgraph

#Bedtools intersect all but Alaska
bedtools intersect -u -a Methylkit_CASS_W_C_10plowcovDMRs.bedGraph -b Methylkit_SVAL_W_C_10plowcovDMRs.bedGraph| bedtools intersect -u -a stdin -b Methylkit_LAT_W_C_10plowcovDMRs.bedGraph > ALEX_SVAL_SWED_Sites_intersect_DMRs_Methylkit10.bedgraph

#Bedtools intersect all but Nunavut - random
bedtools intersect -u -a Methylkit_SVAL_rand_W_C_10plowcovDMRs.bedGraph -b Methylkit_LAT_rand_W_C_10plowcovDMRs.bedGraph | bedtools intersect -u -a stdin -b Methylkit_ALAS_rand_W_C_10plowcovDMRs.bedGraph > SVAL_SWED_ALAS_Sites_intersect_rand_DMRs_Methylkit10.bedgraph

# Counts 
wc -l *DMRs_Methylkit10.bedgraph
# 0 ALEX_SVAL_ALAS_Sites_intersect_DMRs_Methylkit10.bedgraph
# 0 ALEX_SVAL_SWED_Sites_intersect_DMRs_Methylkit10.bedgraph
# 0 ALEX_SWED_ALAS_Sites_intersect_DMRs_Methylkit10.bedgraph
# 0 ALL_Sites_intersect_DMRs_Methylkit10.bedgraph
# 0 SVAL_SWED_ALAS_Sites_intersect_DMRs_Methylkit10.bedgraph
# 0 SVAL_SWED_ALAS_Sites_intersect_rand_DMRs_Methylkit10.bedgraph
# 0 total

# none 

wc -l *10plowcovDMRs.bedGraph
      # 7 Methylkit_ALAS_rand_W_C_10plowcovDMRs.bedGraph
    # 159 Methylkit_ALAS_W_C_10plowcovDMRs.bedGraph
     # 94 Methylkit_CASS_W_C_10plowcovDMRs.bedGraph
      # 2 Methylkit_LAT_rand_W_C_10plowcovDMRs.bedGraph
    # 965 Methylkit_LAT_W_C_10plowcovDMRs.bedGraph
    # 129 Methylkit_SVAL_rand_W_C_10plowcovDMRs.bedGraph
     # 93 Methylkit_SVAL_W_C_10plowcovDMRs.bedGraph


#-------------------------
# # Methylkit10 without overdispersion corrected

#Bedtools intersect all 4
module load StdEnv/2023 bedtools/2.31.0
bedtools intersect -u -a Methylkit_Alex_W_C_10_overdisp_DMRs.bedGraph -b Methylkit_SVAL_W_C_10_overdisp_DMRs.bedGraph | bedtools intersect -u -a stdin -b Methylkit_LAT_W_C_10_overdisp_DMRs.bedGraph | bedtools intersect -u -a stdin -b Methylkit_ALAS_W_C_10_overdisp_DMRs.bedGraph > ALL_Sites_intersect_overdisp_DMRs_Methylkit10.bedgraph

#Bedtools intersect all but Nunavut
bedtools intersect -u -a Methylkit_SVAL_W_C_10_overdisp_DMRs.bedGraph -b Methylkit_LAT_W_C_10_overdisp_DMRs.bedGraph | bedtools intersect -u -a stdin -b Methylkit_ALAS_W_C_10_overdisp_DMRs.bedGraph > SVAL_SWED_ALAS_Sites_intersect_overdisp_DMRs_Methylkit10.bedgraph

#Bedtools intersect all but Svalbard
bedtools intersect -u -a Methylkit_Alex_W_C_10_overdisp_DMRs.bedGraph -b Methylkit_LAT_W_C_10_overdisp_DMRs.bedGraph | bedtools intersect -u -a stdin -b Methylkit_ALAS_W_C_10_overdisp_DMRs.bedGraph > ALEX_SWED_ALAS_Sites_intersect_overdisp_DMRs_Methylkit10.bedgraph

#Bedtools intersect all but Sweden
bedtools intersect -u -a Methylkit_Alex_W_C_10_overdisp_DMRs.bedGraph -b Methylkit_SVAL_W_C_10_overdisp_DMRs.bedGraph| bedtools intersect -u -a stdin -b Methylkit_ALAS_W_C_10_overdisp_DMRs.bedGraph > ALEX_SVAL_ALAS_Sites_intersect_overdisp_DMRs_Methylkit10.bedgraph

#Bedtools intersect all but Alaska
bedtools intersect -u -a Methylkit_Alex_W_C_10_overdisp_DMRs.bedGraph -b Methylkit_SVAL_W_C_10_overdisp_DMRs.bedGraph| bedtools intersect -u -a stdin -b Methylkit_LAT_W_C_10_overdisp_DMRs.bedGraph > ALEX_SVAL_SWED_Sites_intersect_overdisp_DMRs_Methylkit10.bedgraph

#------------------------
bedtools intersect -u -a Methylkit_Alex_rand_W_C_10_overdisp_DMRs.bedGraph -b Methylkit_SVAL_rand_W_C_10_overdisp_DMRs.bedGraph | bedtools intersect -u -a stdin -b Methylkit_LAT_rand_W_C_10_overdisp_DMRs.bedGraph | bedtools intersect -u -a stdin -b Methylkit_ALAS_rand_W_C_10_overdisp_DMRs.bedGraph > ALL_Sites_intersect_rand_overdisp_DMRs_Methylkit10.bedgraph

#Bedtools intersect all but Nunavut - random
bedtools intersect -u -a Methylkit_SVAL_rand_W_C_10_overdisp_DMRs.bedGraph -b Methylkit_LAT_rand_W_C_10_overdisp_DMRs.bedGraph | bedtools intersect -u -a stdin -b Methylkit_ALAS_rand_W_C_10_overdisp_DMRs.bedGraph > SVAL_SWED_ALAS_Sites_intersect_rand_overdisp_DMRs_Methylkit10.bedgraph

#Bedtools intersect all but Svalbard
bedtools intersect -u -a Methylkit_Alex_rand_W_C_10_overdisp_DMRs.bedGraph -b Methylkit_LAT_rand_W_C_10_overdisp_DMRs.bedGraph | bedtools intersect -u -a stdin -b Methylkit_ALAS_rand_W_C_10_overdisp_DMRs.bedGraph > ALEX_SWED_ALAS_Sites_intersect_rand_overdisp_DMRs_Methylkit10.bedgraph

#Bedtools intersect all but Sweden
bedtools intersect -u -a Methylkit_Alex_rand_W_C_10_overdisp_DMRs.bedGraph -b Methylkit_SVAL_rand_W_C_10_overdisp_DMRs.bedGraph| bedtools intersect -u -a stdin -b Methylkit_ALAS_rand_W_C_10_overdisp_DMRs.bedGraph > ALEX_SVAL_ALAS_Sites_intersect_rand_overdisp_DMRs_Methylkit10.bedgraph

#Bedtools intersect all but Alaska
bedtools intersect -u -a Methylkit_Alex_rand_W_C_10_overdisp_DMRs.bedGraph -b Methylkit_SVAL_rand_W_C_10_overdisp_DMRs.bedGraph| bedtools intersect -u -a stdin -b Methylkit_LAT_rand_W_C_10_overdisp_DMRs.bedGraph > ALEX_SVAL_SWED_Sites_intersect_rand_overdisp_DMRs_Methylkit10.bedgraph

#------------------------
# Counts 
wc -l *_overdisp_DMRs_Methylkit10.bedgraph
   # 29 ALEX_SVAL_ALAS_Sites_intersect_overdisp_DMRs_Methylkit10.bedgraph
   # 32 ALEX_SVAL_ALAS_Sites_intersect_rand_overdisp_DMRs_Methylkit10.bedgraph
   
   # 30 ALEX_SVAL_SWED_Sites_intersect_overdisp_DMRs_Methylkit10.bedgraph
   # 27 ALEX_SVAL_SWED_Sites_intersect_rand_overdisp_DMRs_Methylkit10.bedgraph
   
   # 32 ALEX_SWED_ALAS_Sites_intersect_overdisp_DMRs_Methylkit10.bedgraph
   # 33 ALEX_SWED_ALAS_Sites_intersect_rand_overdisp_DMRs_Methylkit10.bedgraph
   
   # 10 ALL_Sites_intersect_overdisp_DMRs_Methylkit10.bedgraph
   # 7 ALL_Sites_intersect_rand_overdisp_DMRs_Methylkit10.bedgraph

  # 142 SVAL_SWED_ALAS_Sites_intersect_overdisp_DMRs_Methylkit10.bedgraph
  # 106 SVAL_SWED_ALAS_Sites_intersect_rand_overdisp_DMRs_Methylkit10.bedgraph

# no pattern

wc -l *10_overdisp_DMRs.bedGraph
   # 4582 Methylkit_ALAS_rand_W_C_10_overdisp_DMRs.bedGraph
   # 5086 Methylkit_ALAS_W_C_10_overdisp_DMRs.bedGraph
    # 783 Methylkit_Alex_rand_W_C_10_overdisp_DMRs.bedGraph
    # 678 Methylkit_Alex_W_C_10_overdisp_DMRs.bedGraph
   # 4328 Methylkit_LAT_rand_W_C_10_overdisp_DMRs.bedGraph
   # 6944 Methylkit_LAT_W_C_10_overdisp_DMRs.bedGraph
   # 6670 Methylkit_SVAL_rand_W_C_10_overdisp_DMRs.bedGraph
   # 6714 Methylkit_SVAL_W_C_10_overdisp_DMRs.bedGraph
  # 35785 total

###############################
# CHH

#Bedtools intersect all 4
module load StdEnv/2023 bedtools/2.31.0

#Bedtools intersect all but Nunavut
bedtools intersect -u -a Methylkit_SVAL_W_C_10_CHH_plowcovDMRs.bedGraph -b Methylkit_LAT_W_C_10_CHH_plowcovDMRs300.bedGraph | bedtools intersect -u -a stdin -b Methylkit_ALAS_W_C_10_CHH_plowcovDMRs.bedGraph > SVAL_SWED_ALAS_Sites_intersect_DMRs_Methylkit_CHH_10.bedgraph

#----------------------------
#Bedtools intersect all but Nunavut - random
bedtools intersect -u -a Methylkit_SVAL_W_C_10_CHH_rand_plowcovDMRs.bedGraph -b Methylkit_LAT_W_C_10_CHH_rand_plowcovDMRs_300.bedGraph | bedtools intersect -u -a stdin -b Methylkit_ALAS_W_C_10_CHH_rand_plowcovDMRs.bedGraph > SVAL_SWED_ALAS_Sites_intersect_rand_DMRs_Methylkit_CHH_10.bedgraph

wc -l *DMRs_Methylkit_CHH_10.bedgraph
0 SVAL_SWED_ALAS_Sites_intersect_DMRs_Methylkit_CHH_10.bedgraph
0 SVAL_SWED_ALAS_Sites_intersect_rand_DMRs_Methylkit_CHH_10.bedgraph
0 total

###################################
# Methylkit and metilene overlap

bedtools intersect -u -a Methylkit_ALAS_W_C_10plowcovDMRs.bedGraph -b Methylkit_ALAS_W_C_10_overdisp_DMRs.bedGraph | bedtools intersect -u -a stdin -b ./metilene/Alaska_W_C.bedGraph > intersect_ALAS_W_C_DMRs_metilene_methylkit_CG.bedgraph
bedtools intersect -u -a Methylkit_SVAL_W_C_10plowcovDMRs.bedGraph -b Methylkit_SVAL_W_C_10_overdisp_DMRs.bedGraph | bedtools intersect -u -a stdin -b ./metilene/Svalbard_W_C.bedGraph > intersect_SVAL_W_C_DMRs_metilene_methylkit_CG.bedgraph
bedtools intersect -u -a Methylkit_LAT_W_C_10plowcovDMRs.bedGraph -b Methylkit_LAT_W_C_10_overdisp_DMRs.bedGraph | bedtools intersect -u -a stdin -b ./metilene/Sweden_W_C.bedGraph > intersect_LAT_W_C_DMRs_metilene_methylkit_CG.bedgraph
bedtools intersect -u -a Methylkit_CASS_W_C_10plowcovDMRs.bedGraph -b Methylkit_Alex_W_C_10_overdisp_DMRs.bedGraph | bedtools intersect -u -a stdin -b ./metilene/Nunavut_W_C.bedGraph > intersect_Alex_W_C_DMRs_metilene_methylkit_CG.bedgraph

bedtools intersect -u -a Methylkit_ALAS_rand_W_C_10plowcovDMRs.bedGraph  -b Alaska_random_W_C.bedGraph > intersect_ALAS_W_C_DMRs_rand_metilene_methylkit_CG.bedgraph
bedtools intersect -u -a Methylkit_SVAL_rand_W_C_10plowcovDMRs.bedGraph -b Svalbard_random_W_C.bedGraph > intersect_SVAL_W_C_DMRs_rand_metilene_methylkit_CG.bedgraph
bedtools intersect -u -a Methylkit_LAT_rand_W_C_10plowcovDMRs.bedGraph  -b Sweden_random_W_C.bedGraph > intersect_LAT_W_C_DMRs_rand_metilene_methylkit_CG.bedgraph
bedtools intersect -u -a Methylkit_CASS_rand_W_C_10plowcovDMRs.bedGraph -b Nunavut_random_W_C.bedGraph > intersect_Alex_W_C_DMRs_rand_metilene_methylkit_CG.bedgraph

wc -l *DMRs_metilene_methylkit_CG.bedgraph
   # 26 intersect_ALAS_W_C_DMRs_metilene_methylkit_CG.bedgraph
    # 3 intersect_Alex_W_C_DMRs_metilene_methylkit_CG.bedgraph
  # 206 intersect_LAT_W_C_DMRs_metilene_methylkit_CG.bedgraph
    # 4 intersect_SVAL_W_C_DMRs_metilene_methylkit_CG.bedgraph

wc -l *DMRs_rand_metilene_methylkit_CG.bedgraph
  0 intersect_ALAS_W_C_DMRs_rand_metilene_methylkit_CG.bedgraph
  0 intersect_Alex_W_C_DMRs_rand_metilene_methylkit_CG.bedgraph
  0 intersect_LAT_W_C_DMRs_rand_metilene_methylkit_CG.bedgraph
  3 intersect_SVAL_W_C_DMRs_rand_metilene_methylkit_CG.bedgraph


wc -l *10plowcovDMRs.bedGraph
    # 159 Methylkit_ALAS_W_C_10plowcovDMRs.bedGraph - 16%
     # 94 Methylkit_CASS_W_C_10plowcovDMRs.bedGraph - 3%
    # 965 Methylkit_LAT_W_C_10plowcovDMRs.bedGraph - 21%
     # 93 Methylkit_SVAL_W_C_10plowcovDMRs.bedGraph - 3%


####################################
# Count intersections
wc -l intersect_*.bedGraph
    0 intersect_ALAS_CASS_W_C_10plowcovDMRs.bedGraph
    0 intersect_ALAS_CpG_CHH_W_C_10plowcovDMRs.bedGraph
   28 intersect_ALAS_HL_10plowcovDMRs.bedGraph
    2 intersect_ALAS_LAT_W_C_10plowcovDMRs.bedGraph
    6 intersect_ALAS_Pheno_10plowcovDMRs.bedGraph
  152 intersect_HL_LAT_10plowcovDMRs.bedGraph
    0 intersect_LAT_CpG_CHH_W_C_10plowcovDMRs.bedGraph
   34 intersect_LAT_SE_W_C_10plowcovDMRs.bedGraph
 4984 intersect_Pheno_CpG_CHH_10plowcovDMRs.bedGraph
 1229 intersect_Pheno_HL_10plowcovDMRs.bedGraph
   76 intersect_Pheno_LAT_10plowcovDMRs.bedGraph
    0 intersect_rand_ALAS_LAT_W_C_10plowcovDMRs.bedGraph
    0 intersect_rand-non_ALAS_W_C_10plowcovDMRs.bedGraph
    1 intersect_rand-non_LAT_W_C_10plowcovDMRs.bedGraph
    0 intersect_SE_WC_HL_10plowcovDMRs.bedGraph
    3 intersect_SVAL_LAT_W_C_10plowcovDMRs.bedGraph


wc -l *.bedGraph
      7 Methylkit_ALAS_rand_W_C_10plowcovDMRs.bedGraph
     57 Methylkit_ALAS_rand_W_C_25plowcovDMRs.bedGraph
      6 Methylkit_ALAS_W_C_10_CHH_plowcovDMRs.bedGraph
      3 Methylkit_ALAS_W_C_10_CHH_rand_plowcovDMRs.bedGraph
    159 Methylkit_ALAS_W_C_10plowcovDMRs.bedGraph
     21 Methylkit_ALAS_W_C_25plowcovDMRs.bedGraph
     14 Methylkit_CASS_rand_W_C_25plowcovDMRs.bedGraph
     94 Methylkit_CASS_W_C_10plowcovDMRs.bedGraph
     12 Methylkit_CASS_W_C_25plowcovDMRs.bedGraph
    222 Methylkit_DRY_rand_W_C_25plowcovDMRs.bedGraph
    147 Methylkit_DRY_W_C_10plowcovDMRs.bedGraph
     19 Methylkit_DRY_W_C_25plowcovDMRs.bedGraph
    448 Methylkit_FERT_rand_W_C_25plowcovDMRs.bedGraph
   9437 Methylkit_FERT_W_C_10plowcovDMRs.bedGraph
    550 Methylkit_FERT_W_C_25plowcovDMRs.bedGraph
  12823 Methylkit_HL_10plowcovDMRs.bedGraph
   1465 Methylkit_HL_25plowcovDMRs.bedGraph
      2 Methylkit_HL_rand_10plowcovDMRs.bedGraph
      2 Methylkit_HL_rand_25plowcovDMRs.bedGraph
      2 Methylkit_LAT_rand_W_C_10plowcovDMRs.bedGraph
     21 Methylkit_LAT_rand_W_C_25plowcovDMRs.bedGraph
     28 Methylkit_LAT_W_C_10_CHH_plowcovDMRs300.bedGraph
     10 Methylkit_LAT_W_C_10_CHH_plowcovDMRs.bedGraph
      2 Methylkit_LAT_W_C_10_CHH_rand_plowcovDMRs_300.bedGraph
    965 Methylkit_LAT_W_C_10plowcovDMRs.bedGraph
    115 Methylkit_LAT_W_C_25plowcovDMRs.bedGraph
    424 Methylkit_MEAD_rand_W_C_25plowcovDMRs.bedGraph
   8527 Methylkit_MEAD_W_C_10plowcovDMRs.bedGraph
    478 Methylkit_MEAD_W_C_25plowcovDMRs.bedGraph
  14384 Methylkit_Pheno_10_DMRs.bedGraph
      5 Methylkit_Pheno_10_DMRs_CHH_rand.bedGraph
   3531 Methylkit_Pheno_25_DMRs.bedGraph
      3 Methylkit_Pheno_25_DMRs_CHH_rand.bedGraph
    416 Methylkit_SE_HL_10_DMRs.bedGraph
     53 Methylkit_SE_HL_25_DMRs.bedGraph
     12 Methylkit_SE_HL_rand_10_DMRs.bedGraph
      2 Methylkit_SE_HL_rand_25_DMRs.bedGraph
   1517 Methylkit_SE_W_C_10_DMRs.bedGraph
     40 Methylkit_SE_W_C_10_DMRs_CHH.bedGraph
     10 Methylkit_SE_W_C_10_DMRs_CHH_rand.bedGraph
    370 Methylkit_SE_W_C_25_DMRs.bedGraph
     15 Methylkit_SE_W_C_25_DMRs_CHH.bedGraph
    129 Methylkit_SVAL_rand_W_C_10plowcovDMRs.bedGraph
     15 Methylkit_SVAL_rand_W_C_25plowcovDMRs.bedGraph
      4 Methylkit_SVAL_W_C_10_CHH_plowcovDMRs.bedGraph
      4 Methylkit_SVAL_W_C_10_CHH_rand_plowcovDMRs.bedGraph
     93 Methylkit_SVAL_W_C_10plowcovDMRs.bedGraph
     12 Methylkit_SVAL_W_C_25plowcovDMRs.bedGraph
    406 Methylkit_WILL_rand_W_C_25plowcovDMRs.bedGraph
   7853 Methylkit_WILL_W_C_10plowcovDMRs.bedGraph
    395 Methylkit_WILL_W_C_25plowcovDMRs.bedGraph
      9 Rand_Methylkit_Pheno_10_DMRs.bedGraph
      2 Rand_Methylkit_Pheno_25_DMRs.bedGraph
    145 Rand_Methylkit_SE_W_C_10_DMRs.bedGraph
     30 Rand_Methylkit_SE_W_C_25_DMRs.bedGraph



#-------------
# Run all combinations
# Directory containing the bedGraph files
bedgraph_dir="/lustre04/scratch/celphin/Dryas_large_folders/intersections"

# Create an output directory if it doesn't exist
mkdir -p "$bedgraph_dir/intersect_results"

# Loop through all the bedGraph files in the directory
for a in "$bedgraph_dir"/*.bedGraph; do
  for b in "$bedgraph_dir"/*.bedGraph; do
    # Skip if the files are the same
    if [ "$a" != "$b" ]; then
      # Get the output file name based on the input files
      output_file="$bedgraph_dir/intersect_results/intersect_$(basename "$a" .bedGraph)_$(basename "$b" .bedGraph).bedGraph"
      
      # Run bedtools intersect
      bedtools intersect -u -a "$a" -b "$b" > "$output_file"
      
      echo "Intersected $a and $b, output saved to $output_file"
    fi
  done
done

##########################################


