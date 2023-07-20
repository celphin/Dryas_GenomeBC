###########################################################
#MS_Bedgraph_Intersection_Notes.bash
    #Partly from MS_Blast_Notes.bash
    #Partly from MS_Wild_Site_Specific_Metilene_Notes.bash
#Required input:
    #All output bedgraphs from metilene for statistics: 70_5_4_0.9_qval.0.001
#Output:
    #Series of intersections bedgraphs
    #Location: 
#Todo:
    #Update file copy location, add slurm/directory/other such notes
    #Maybe put to metilene folder
    #Remove from original files
###########################################################

module load bedtools/2.30.0

#Update this
cp /home/msandler/projects/def-rieseber/Dryas_shared_data/MS_Phenology_metilene_output_bedgraphs/metilene_Mat_Sen_70_5_4_0.9_qval.0.001.bedgraph Mat_Sen.bedGraph
cp /home/msandler/projects/def-rieseber/Dryas_shared_data/MS_Wild_metilene_output_bedgraphs/All_Wild_Warming_bedgraphs/Wild_W_C_70_5_4_0.9_qval.0.001.bedgraph Wild_W_C.bedGraph
cp /home/msandler/projects/def-rieseber/Dryas_shared_data/MS_Seedling_metilene_output_bedgraphs/SE_W_C_70_5_4_0.9_qval.0.001.bedgraph SE_W_C.bedGraph
cp /home/msandler/projects/def-rieseber/Dryas_shared_data/MS_Seedling_metilene_output_bedgraphs/SE_L_H_70_5_4_0.9_qval.0.001.bedgraph SE_L_H.bedGraph
cp /home/msandler/projects/def-rieseber/Dryas_shared_data/MS_Wild_metilene_output_bedgraphs/All_Wild_Species_bedgraphs/Wild_Species_DO_DI_70_5_4_0.9_qval.0.001.bedgraph Wild_DO_DI.bedGraph
cp /home/msandler/projects/def-rieseber/Dryas_shared_data/MS_Wild_metilene_output_bedgraphs/All_Wild_Latitude_bedgraphs/Wild_Lat_L_H_70_5_4_0.9_qval.0.001.bedgraph Wild_Lat_W_C.bedGraph
cp /home/msandler/projects/def-rieseber/Dryas_shared_data/MS_Wild_metilene_output_bedgraphs/True_Parent_Warming_bedgraphs/P_W_C_70_5_4_0.9_qval.0.001.bedgraph Parent_W_C.bedGraph
cp Nunavut_Metilene/Nunavut_W_C_70_5_4_0.9_qval.0.001.bedgraph Nunavut_W_C.bedGraph
cp Svalbard_Metilene/Svalbard_W_C_70_5_4_0.9_qval.0.001.bedgraph Svalbard_W_C.bedGraph
cp Sweden_Metilene/Sweden_W_C_70_5_4_0.9_qval.0.001.bedgraph Sweden_W_C.bedGraph
cp Alaska_Metilene/Alaska_W_C_70_5_4_0.9_qval.0.001.bedgraph Alaska_W_C.bedGraph

#Intersect Warming + Phenology DMRS:
bedtools intersect -u -a Wild_W_C.bedGraph -b Mat_Sen.bedGraph > intersect_Wild_W_C_Mat_Sen.bedGraph
#Total subtract: Warming-Phenology DMRS:
bedtools subtract -A -a Wild_W_C.bedGraph -b Mat_Sen.bedGraph > total_subtract_W_C_Mat_Sen.bedGraph
#Intersection: Seedling Warming, Seedling Low/High
bedtools intersect -u -a SE_W_C.bedGraph -b SE_L_H.bedGraph > intersect_SE_W_C_SE_L_H.bedGraph
#Total Subtract: Seedling Warming - Seedling Low High
bedtools subtract -A -a SE_W_C.bedGraph -b SE_L_H.bedGraph > total_subtract_SE_W_C_SE_L_H.bedGraph
#Include intersect warming Seedling Warming + Parents:
bedtools intersect -u -a SE_W_C.bedGraph -b Parent_W_C.bedGraph > intersect_SE_W_C_P_W_C.bedGraph
#Intersect Seedling Warming + Wild Warming
bedtools intersect -u -a SE_W_C.bedGraph -b Wild_W_C.bedGraph > intersect_SE_W_C_Wild_W_C.bedGraph
#Total subtracts Seedling Warming - Parents
bedtools subtract -A -a SE_W_C.bedGraph -b Parent_W_C.bedGraph > total_subtract_SE_W_C_P_W_C.bedGraph
#Intersect: Seedling Low High + Wild Site DMRs (Sweden and Alaska)  vs (Alex + Svalbard)
bedtools intersect -u -a SE_L_H.bedGraph -b Wild_Lat_L_H.bedGraph > intersect_SE_L_H_Wild_L_H.bedGraph
#Bedtools intersect all 4
bedtools intersect -u -a Nunavut_W_C.bedGraph -b Svalbard_W_C.bedGraph| bedtools intersect -u -a stdin -b Sweden_W_C.bedGraph | bedtools intersect -u -a stdin -b Alaska_W_C.bedGraph > ALL_Sites_intersect_DMRs.bedgraph
#Bedtools intersect all but Nunavut
bedtools intersect -u -a Svalbard_W_C.bedGraph -b Sweden_W_C.bedGraph | bedtools intersect -u -a stdin -b Alaska_W_C.bedGraph > SVAL_SWED_ALAS_Sites_intersect_DMRs.bedgraph
#Bedtools intersect all but Svalbard
bedtools intersect -u -a Nunavut_W_C.bedGraph -b Sweden_W_C.bedGraph | bedtools intersect -u -a stdin -b Alaska_W_C.bedGraph > ALEX_SWED_ALAS_Sites_intersect_DMRs.bedgraph
#Bedtools intersect all but Sweden
bedtools intersect -u -a Nunavut_W_C.bedGraph -b Svalbard_W_C.bedGraph| bedtools intersect -u -a stdin -b Alaska_W_C.bedGraph > ALEX_SVAL_ALAS_Sites_intersect_DMRs.bedgraph
#Bedtools intersect all but Alaska
bedtools intersect -u -a Nunavut_W_C.bedGraph -b Svalbard_W_C.bedGraph| bedtools intersect -u -a stdin -b Sweden_W_C.bedGraph > ALEX_SVAL_SWED_Sites_intersect_DMRs.bedgraph

