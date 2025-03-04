############################################
#MS_Seedling_VS_Parent_Plotting.bash
    #Inheritence Plots Notes
    #Creates: 
        #Scatter plot of Parents vs Seedling
        #Bar plots for DMRS comparing parents and seedlings
        #Bar plots for WC for the whole SE/Parent individuals
#Requires: Box_Script_SE_P.R
########################################
#Create Intersections: 
cp /home/msandler/projects/def-rieseber/Dryas_shared_data/MS_Seedling_metilene_output_bedgraphs/SE_W_C_70_5_4_0.9_qval.0.001.bedgraph SE_W_C.bedGraph
cp /home/msandler/projects/def-rieseber/Dryas_shared_data/MS_Wild_metilene_output_bedgraphs/True_Parent_Warming_bedgraphs/P_W_C_70_5_4_0.9_qval.0.001.bedgraph P_W_C.bedGraph
#Copy Site specific ones 

#Intersect: Sweden, Parent, Seedling
module load bedtools/2.30.0
bedtools intersect -u -a SE_W_C.bedGraph -b Sweden_W_C.bedGraph| bedtools intersect -u -a stdin -b P_W_C.bedGraph > Swed_SE_P_W_C_intersect.bedGraph


bedtools intersect -u -a Swed_SE_P_W_C_intersect.bedGraph -b Alaska_W_C.bedGraph > Alas_Swed_SE_P_W_C_intersect.bedGraph
bedtools intersect -u -a Swed_SE_P_W_C_intersect.bedGraph -b Nunavut_W_C.bedGraph > Alex_Swed_SE_P_W_C_intersect.bedGraph
bedtools intersect -u -a Swed_SE_P_W_C_intersect.bedGraph -b Svalbard_W_C.bedGraph > Sval_Swed_SE_P_W_C_intersect.bedGraph

wc -l *Swed*intersect*
#Intersect: Sweden, Parent, Seedling
module load bedtools/2.30.0
bedtools intersect -u -a SE_W_C.bedGraph -b Sweden_W_C.bedGraph| bedtools intersect -u -a stdin -b P_W_C.bedGraph > Swed_SE_P_W_C_intersect.bedGraph


bedtools intersect -u -a Swed_SE_P_W_C_intersect.bedGraph -b Alaska_W_C.bedGraph > Alas_Swed_SE_P_W_C_intersect.bedGraph
bedtools intersect -u -a Swed_SE_P_W_C_intersect.bedGraph -b Nunavut_W_C.bedGraph > Alex_Swed_SE_P_W_C_intersect.bedGraph
bedtools intersect -u -a Swed_SE_P_W_C_intersect.bedGraph -b Svalbard_W_C.bedGraph > Sval_Swed_SE_P_W_C_intersect.bedGraph

##########################################
#Plotting: 
module load StdEnv/2020
module load r/4.2.1
module load gdal
module load udunits
module load python
export R_LIBS_USER=/home/msandler/R/x86_64-pc-linux-gnu-library/4.2.1/
#####################

while read dmr
do

chrom=$(awk '{print $1}' <<< $dmr)
start=$(awk '{print $2}' <<< $dmr)
end=$(awk '{print $3}' <<< $dmr)


grep "^${chrom}" data/P_metilene_W_C.input > "data/P_${chrom}.txt"
grep "^${chrom}" data/SE_metilene_W_C.input > "data/SE_${chrom}.txt"

awk -F "\t" -v start="$start" -v end="$end" '{ if(($2 >= start) && ($2 <= end)) { print } }' "data/P_${chrom}.txt"  > "data/P_${chrom}_${start}_${end}.txt"
awk -F "\t" -v start="$start" -v end="$end" '{ if(($2 >= start) && ($2 <= end)) { print } }' "data/SE_${chrom}.txt"  > "data/SE_${chrom}_${start}_${end}.txt"

Rscript Box_Script_SE_P.R ${chrom} $start $end

done < Alas_Swed_SE_P_W_C_intersect.bedGraph
#############################
#Plotting inheritane for WC-MatSen
cp ~/projects/def-rieseber/Dryas_shared_data/MS_blast_input_bedgraphs/total_subtract_W_C_Mat_Sen.bedGraph .


while read dmr
do

chrom=$(awk '{print $1}' <<< $dmr)
start=$(awk '{print $2}' <<< $dmr)
end=$(awk '{print $3}' <<< $dmr)


grep "^${chrom}" data/P_metilene_W_C.input > "data/P_${chrom}.txt"
grep "^${chrom}" data/SE_metilene_W_C.input > "data/SE_${chrom}.txt"

awk -F "\t" -v start="$start" -v end="$end" '{ if(($2 >= start) && ($2 <= end)) { print } }' "data/P_${chrom}.txt"  > "data/P_${chrom}_${start}_${end}.txt"
awk -F "\t" -v start="$start" -v end="$end" '{ if(($2 >= start) && ($2 <= end)) { print } }' "data/SE_${chrom}.txt"  > "data/SE_${chrom}_${start}_${end}.txt"

Rscript Box_Script_SE_P.R ${chrom} $start $end

done < total_subtract_W_C_Mat_Sen.bedGraph
######################################