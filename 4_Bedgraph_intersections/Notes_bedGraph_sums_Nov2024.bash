###########################################################
# Sum Bedgegraph files from methylseq for viewing in IGV
# Nov 2024
#########################################

tmux new-session -s bedgraph
tmux attach-session -t bedgraph

cd ~/scratch/Dryas/bedgraph_old_ref

cp -r /home/celphin/projects/rpp-rieseber/celphin/Dryas/Methylation_calling/May2021_Parents/bismark_methylation_calls/bedGraph* ~/scratch/Dryas/bedgraph_old_ref

gunzip *deduplicated.bedGraph.gz

module load StdEnv/2023 bedtools/2.31.0

# test
bedtools unionbedg -i W*SVAL*.deduplicated.bedGraph | \
awk 'OFS="\t" {sum=0; for (col=4; col<=NF; col++) sum += $col; print $1, $2, $3, sum/(NF-4+1); }' > W_SVAL.bedGraph

bedtools unionbedg -i W*SVAL*.deduplicated.bedGraph > W_SVAL_total.bedGraph

head W_SVAL_total.bedGraph
# Do1_00107       8       9       0       0       0       0       0       0
# Do1_00107       9       10      0       0       0       0       0       0
# Do1_00107       137     138     12.5    0       16.6666666666667        0       0       0
# Do1_00107       138     139     28.5714285714286        0       7.69230769230769        12.5 00
# Do1_00107       247     248     0       0       46.1538461538462        0       0       14.2857142857143
# Do1_00107       248     249     0       0       8.69565217391304        0       0       0
# Do1_00107       572     573     30      50      40      0       25      38.0952380952381
# Do1_00107       573     574     41.6666666666667        50      27.2727272727273        57.1428571428571      62.5    54.5454545454545
# Do1_00107       580     581     77.7777777777778        0       75      33.3333333333333     027.2727272727273
# Do1_00107       581     582     82.6086956521739        25      36.3636363636364        16.6666666666667      50      81.8181818181818

head W_SVAL.bedGraph
# Do1_00107       8       9       0
# Do1_00107       9       10      0
# Do1_00107       137     138     4.86111
# Do1_00107       138     139     8.12729
# Do1_00107       247     248     10.0733
# Do1_00107       248     249     1.44928
# Do1_00107       572     573     30.5159
# Do1_00107       573     574     48.8546
# Do1_00107       580     581     35.564
# Do1_00107       581     582     48.7429

# Working!

#############################
# Now for all combinations

module load StdEnv/2023 bedtools/2.31.0

bedtools unionbedg -i W*SVAL*.deduplicated.bedGraph > W_SVAL_total.bedGraph
awk 'OFS="\t" {sum=0; for (col=4; col<=NF; col++) sum += $col; print $1, $2, $3, sum/(NF-4+1); }' W_SVAL_total.bedGraph > W_SVAL.bedGraph

bedtools unionbedg -i C*SVAL*.deduplicated.bedGraph > C_SVAL_total.bedGraph
awk 'OFS="\t" {sum=0; for (col=4; col<=NF; col++) sum += $col; print $1, $2, $3, sum/(NF-4+1); }' C_SVAL_total.bedGraph > C_SVAL.bedGraph

#---------------
bedtools unionbedg -i W*LAT*.deduplicated.bedGraph > W_Swed_total.bedGraph
awk 'OFS="\t" {sum=0; for (col=4; col<=NF; col++) sum += $col; print $1, $2, $3, sum/(NF-4+1); }' W_Swed_total.bedGraph > W_SwedC.bedGraph

bedtools unionbedg -i C*LAT*.deduplicated.bedGraph > C_Swed_total.bedGraph
awk 'OFS="\t" {sum=0; for (col=4; col<=NF; col++) sum += $col; print $1, $2, $3, sum/(NF-4+1); }' C_Swed_total.bedGraph > C_SwedC.bedGraph

#---------------
bedtools unionbedg -i W*CASS*.deduplicated.bedGraph > W_CASS_total.bedGraph
awk 'OFS="\t" {sum=0; for (col=4; col<=NF; col++) sum += $col; print $1, $2, $3, sum/(NF-4+1); }' W_CASS_total.bedGraph > W_CASS.bedGraph

bedtools unionbedg -i C*CASS*.deduplicated.bedGraph  > C_CASS_total.bedGraph
awk 'OFS="\t" {sum=0; for (col=4; col<=NF; col++) sum += $col; print $1, $2, $3, sum/(NF-4+1); }' C_CASS_total.bedGraph > C_CASS.bedGraph

#---------------
bedtools unionbedg -i W*WILL*.deduplicated.bedGraph > W_WILL_total.bedGraph
awk 'OFS="\t" {sum=0; for (col=4; col<=NF; col++) sum += $col; print $1, $2, $3, sum/(NF-4+1); }' W_WILL_total.bedGraph > W_WILL.bedGraph

bedtools unionbedg -i C*WILL*.deduplicated.bedGraph > C_WILL_total.bedGraph
awk 'OFS="\t" {sum=0; for (col=4; col<=NF; col++) sum += $col; print $1, $2, $3, sum/(NF-4+1); }' C_WILL_total.bedGraph > C_WILL.bedGraph

#---------------
bedtools unionbedg -i W*FERT*.deduplicated.bedGraph > W_FERT_total.bedGraph
awk 'OFS="\t" {sum=0; for (col=4; col<=NF; col++) sum += $col; print $1, $2, $3, sum/(NF-4+1); }' W_FERT_total.bedGraph > W_FERT.bedGraph

bedtools unionbedg -i C*FERT*.deduplicated.bedGraph > C_FERT_total.bedGraph
awk 'OFS="\t" {sum=0; for (col=4; col<=NF; col++) sum += $col; print $1, $2, $3, sum/(NF-4+1); }' C_FERT_total.bedGraph > C_FERT.bedGraph

#---------------
bedtools unionbedg -i W*DRY*.deduplicated.bedGraph > W_DRY_total.bedGraph
awk 'OFS="\t" {sum=0; for (col=4; col<=NF; col++) sum += $col; print $1, $2, $3, sum/(NF-4+1); }' W_DRY_total.bedGraph > W_DRY.bedGraph

bedtools unionbedg -i C*DRY*.deduplicated.bedGraph > C_DRY_total.bedGraph
awk 'OFS="\t" {sum=0; for (col=4; col<=NF; col++) sum += $col; print $1, $2, $3, sum/(NF-4+1); }' C_DRY_total.bedGraph > C_DRY.bedGraph

#---------------
bedtools unionbedg -i W*MEAD*.deduplicated.bedGraph > W_MEAD_total.bedGraph
awk 'OFS="\t" {sum=0; for (col=4; col<=NF; col++) sum += $col; print $1, $2, $3, sum/(NF-4+1); }' W_MEAD_total.bedGraph > W_MEAD.bedGraph

bedtools unionbedg -i C*MEAD*.deduplicated.bedGraph > C_MEAD_total.bedGraph
awk 'OFS="\t" {sum=0; for (col=4; col<=NF; col++) sum += $col; print $1, $2, $3, sum/(NF-4+1); }' C_MEAD_total.bedGraph > C_MEAD.bedGraph

#---------------
bedtools unionbedg -i W*ALAS*.deduplicated.bedGraph  > W_ALAS_total.bedGraph
awk 'OFS="\t" {sum=0; for (col=4; col<=NF; col++) sum += $col; print $1, $2, $3, sum/(NF-4+1); }' W_ALAS_total.bedGraph > W_ALAS.bedGraph

bedtools unionbedg -i C*ALAS*.deduplicated.bedGraph > C_ALAS_total.bedGraph
awk 'OFS="\t" {sum=0; for (col=4; col<=NF; col++) sum += $col; print $1, $2, $3, sum/(NF-4+1); }' C_ALAS_total.bedGraph > C_ALAS.bedGraph
