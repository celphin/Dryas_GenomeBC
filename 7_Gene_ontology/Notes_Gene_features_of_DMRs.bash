###########################
# Determining the gene features of each DMR with SNPEff
# Dec 2024
##########################
# https://github.com/RobertsLab/project-gigas-oa-meth/blob/master/code/10-Genomic-Location-of-DML.ipynb 
# https://www.biostars.org/p/432180/ 
# https://pcingola.github.io/SnpEff/snpeff/inputoutput/#bed-files 

# Download SNPEff
wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip
unzip snpEff_latest_core.zip

# Example annotating variants in vcf file
java -Xmx8g -jar snpEff.jar GRCh37.75 examples/test.chr22.vcf > test.chr22.ann.vcf

# Example annotating variants in bed file
java -Xmx8g -jar snpEff.jar -i bed BDGP5.69 chipSeq_peaks.bed

#--------------------------
# setting up an annotation 
# https://pcingola.github.io/SnpEff/snpeff/build_db/

# https://pcingola.github.io/SnpEff/snpeff/build_db/#step-1-configure-a-new-genome
mkdir /home/celphin/scratch/Dryas/snpEff/data/OldDoct/
cd /home/celphin/scratch/Dryas/snpEff/data/OldDoct/

# copy over whole genome and gff files from Cedar

# rename files
rename Dryas_octopetala_H1 OldDoct *

mv OldDoct.CDS.fa cds.fa
mv  OldDoct.gff3 genes.gff
mv OldDoct.genome.fa sequences.fa
mv OldDoct.protein.fa protein.fa

# edit config file to add
nano snpEffect.config

# Dryas octopetala genome, version Old
OldDoct.genome : OldDoct

#-------------------------------
tmux new-session -s snpEff
tmux attach-session -t snpEff

salloc -c1 --time 3:00:00 --mem 120000m --account def-cronk

module load StdEnv/2023 java/21.0.1
cd /home/celphin/scratch/Dryas/snpEff
java -Xmx8g -jar snpEff.jar build -gff3 -v OldDoct

# depending if the annotation has protein or cds files you can include:
# -noCheckProtein
# -noCheckCds

#  Protein check:  OldDoct OK: 40930       Not found: 0    Errors: 251     Error percentage: 0.6095043830892887%

# 00:01:23 [Optional] Reading regulation elements: GFF
# WARNING_FILE_NOT_FOUND: Cannot read optional regulation file '/lustre04/scratch/celphin/Dryas/snpEff/./data/OldDoct/regulation.gff', nothing done.
# 00:01:23 [Optional] Reading regulation elements: BED
# 00:01:23 Cannot find optional regulation dir '/lustre04/scratch/celphin/Dryas/snpEff/./data/OldDoct/regulation.bed/', nothing done.
# 00:01:23 [Optional] Reading motifs: GFF
# WARNING_FILE_NOT_FOUND: Cannot open PWMs file /lustre04/scratch/celphin/Dryas/snpEff/./data/OldDoct/pwms.bin. Nothing done

#-------------------------------
cd /home/celphin/scratch/Dryas/snpEff
mkdir Dryas_DMRs
# copy over DMRs

# convert all bedGraph files to bed files
grep -v track file.bedGraph | awk '{print $1 "\t" $2 "\t" $3}' > file.bed

# Alaska_W_C.bedGraph                           SE_L_H.bedGraph
# ALEX_SVAL_ALAS_Sites_intersect_DMRs.bedgraph  SE_W_C.bedGraph
# ALEX_SVAL_SWED_Sites_intersect_DMRs.bedgraph  Summary_Data_Bedgraph_Intersection_DMR_Count
# ALEX_SWED_ALAS_Sites_intersect_DMRs.bedgraph  Svalbard_W_C.bedGraph
# ALL_Sites_intersect_DMRs.bedgraph             SVAL_NUN_Sites_intersect_DMRs.bedgraph
# C_ALAS.bedGraph                               SVAL_SWED_ALAS_Sites_intersect_DMRs.bedgraph
# CASS5C_529_159_Non_CpG.bedGraph               SVAL_SWED_Sites_intersect_DMRs.bedgraph
# C_CASS.bedGraph                               SWED_ALAS_Sites_intersect_DMRs.bedgraph
# C_DRY.bedGraph                                Sweden_W_C.bedGraph
# C_FERT.bedGraph                               total_subtract_SE_W_C_P_W_C.bedGraph
# C_MEAD.bedGraph                               total_subtract_SE_W_C_SE_L_H.bedGraph
# C_SVAL.bedGraph                               total_subtract_W_C_Mat_Sen.bedGraph
# C_SwedC.bedGraph                              W_ALAS.bedGraph
# C_WILL.bedGraph                               W_CASS.bedGraph
# intersect_SE_L_H_Wild_L_H.bedGraph            W_DRY.bedGraph
# intersect_SE_W_C_P_W_C.bedGraph               W_FERT.bedGraph
# intersect_SE_W_C_SE_L_H.bedGraph              Wild_Lat_L_H.bedGraph
# intersect_SE_W_C_Wild_W_C.bedGraph            Wild_W_C.bedGraph
# intersect_Wild_W_C_Mat_Sen.bedGraph           W_MEAD.bedGraph
# Mat_Sen.bedGraph                              W_SVAL.bedGraph
# NUN_ALAS_Sites_intersect_DMRs.bedgraph        W_SwedC.bedGraph
# Nunavut_W_C.bedGraph                          W_WILL.bedGraph
# Parent_W_C.bedGraph

FILES=$(ls /home/celphin/scratch/Dryas/snpEff/Dryas_DMRs/*.bedGraph)
echo ${FILES}

for file in ${FILES}; do
    # Run the command on the file and create the .bed file
    grep -v track "$file" | awk '{print $1 "\t" $2 "\t" $3}' > "${file%.bedGraph}.bed"
    echo "Processed $file"
done


#----------------------------------
# annotating the DMRs
# https://pcingola.github.io/SnpEff/snpeff/inputoutput/#bed-files 

FILES=$(ls cd /home/celphin/scratch/Dryas/snpEff/Dryas_DMRs/*.bed)
echo ${FILES}

cd /home/celphin/scratch/Dryas/snpEff

for file in ${FILES}; do
java -Xmx8g -jar snpEff.jar -i bed OldDoct "$file" > "$file".out
echo "Processed $file"
done

# SnpEff version 5.2e (build 2024-10-04 18:09), by Pablo Cingolani
# Command line: SnpEff  -i bed OldDoct /home/celphin/scratch/Dryas/snpEff/Dryas_DMRs/Alaska_W_C.bed
# Chromo        Start   End     Name;Effect|Gene|BioType        Score
# Do1_00107       195322  195749  line_1;Intergenic:Do1_00107G00046-Do1_00107G00047;Upstream:4688|Transcript:Do1_00107G00044V1.1:protein_coding|Gene:Do1_00107G00044:protein_coding;Exon:1:1:RETAINED|Transcript:Do1_00107G00047V1.1:protein_coding|Gene:Do1_00107G00047:protein_coding;Downstream:89|Transcript:Do1_00107G00047V1.1:protein_coding|Gene:Do1_00107G00047:protein_coding;Upstream:1007|Transcript:Do1_00107G00046V1.1:protein_coding|Gene:Do1_00107G00046:protein_coding;Upstream:2661|Transcript:Do1_00107G00045V1.1:protein_coding|Gene:Do1_00107G00045:protein_coding;Intergenic:Do1_00107G00047-Do1_00107G00048;Upstream:0|Transcript:Do1_00107G00047V1.1:protein_coding|Gene:Do1_00107G00047:protein_coding
# Do1_00151       12865   12998   line_2
# Do1_00153       43724   44084   line_3;Intergenic:Do1_00153G00006-Do1_00153G00007;Downstream:4795|Transcript:Do1_00153G00005V1.1:protein_coding|Gene:Do1_00153G00005:protein_coding;Downstream:1898|Transcript:Do1_00153G00008V1.1:protein_coding|Gene:Do1_00153G00008:protein_coding;Upstream:163|Transcript:Do1_00153G00006V1.1:protein_coding|Gene:Do1_00153G00006:protein_coding;Upstream:1457|Transcript:Do1_00153G00007V1.1:protein_coding|Gene:Do1_00153G00007:protein_coding
# Do1_00165       44980   45178   line_4;Upstream:2447|Transcript:Do1_00165G00004V1.1:protein_coding|Gene:Do1_00165G00004:protein_coding;Intergenic:Do1_00165G00005-Do1_00165G00006;Upstream:999|Transcript:Do1_00165G00005V1.1:protein_coding|Gene:Do1_00165G00005:protein_coding
# Do1_00167       81807   81969   line_5;Intergenic:Do1_00167G00015-CHR_END
# Do1_00206       48138   48507   line_6;Downstream:2456|Transcript:Do1_00206G00012V1.1:protein_coding|Gene:Do1_00206G00012:protein_coding;Intergenic:Do1_00206G00009-Do1_00206G00010;Downstream:4120|Transcript:Do1_00206G00008V1.1:protein_coding|Gene:Do1_00206G00008:protein_coding;Downstream:2044|Transcript:Do1_00206G00011V1.1:protein_coding|Gene:Do1_00206G00011:protein_coding;Downstream:578|Transcript:Do1_00206G00010V1.1:protein_coding|Gene:Do1_00206G00010:protein_coding;Downstream:2932|Transcript:Do1_00206G00009V1.1:protein_coding|Gene:Do1_00206G00009:protein_coding
# Do1_00211       20021   20090   line_7
# Do1_00233       29596   29949   line_8;Upstream:3731|Transcript:Do1_00233G00004V1.1:protein_coding|Gene:Do1_00233G00004:protein_coding;Intergenic:Do1_00233G00003-Do1_00233G00004;Downstream:939|Transcript:Do1_00233G00003V1.1:protein_coding|Gene:Do1_00233G00003:protein_coding
# Do1_00276       22230   22535   line_9

# Do1_01_a00001   3955238 3955782 line_21;Gene:Do1_01_a00001G00761:protein_coding;Upstream:1141|Transcript:Do1_01_a00001G00760V1.1:protein_coding|Gene:Do1_01_a00001G00760:protein_coding;Downstream:3187|Transcript:Do1_01_a00001G00762V1.1:protein_coding|Gene:Do1_01_a00001G00762:protein_coding
# Do1_01_a00001   4248114 4248158 line_22;Downstream:1490|Transcript:Do1_01_a00001G00815V1.1:protein_coding|Gene:Do1_01_a00001G00815:protein_coding;Exon:1:1:RETAINED|Transcript:Do1_01_a00001G00814V1.1:protein_coding|Gene:Do1_01_a00001G00814:protein_coding;Downstream:3965|Transcript:Do1_01_a00001G00813V1.1:protein_coding|Gene:Do1_01_a00001G00813:protein_coding
# Do1_01_a00001   4303561 4303898 line_23;Upstream:3543|Transcript:Do1_01_a00001G00818V1.1:protein_coding|Gene:Do1_01_a00001G00818:protein_coding;Intergenic:Do1_01_a00001G00818-Do1_01_a00001G00819
# Do1_01_a00001   4640406 4640646 line_24;Gene:Do1_01_a00001G00871:protein_coding
# Do1_01_a00001   4698876 4699117 line_25;Intron:14:14:RETAINED-RETAINED|Transcript:Do1_01_a00001G00877V1.1:protein_coding|Gene:Do1_01_a00001G00877:protein_coding
# Do1_01_a00001   4806306 4806403 line_26;Downstream:269|Transcript:Do1_01_a00001G00896V1.1:protein_coding|Gene:Do1_01_a00001G00896:protein_coding;Intergenic:Do1_01_a00001G00896-Do1_01_a00001G00897;Downstream:2742|Transcript:Do1_01_a00001G00897V1.1:protein_coding|Gene:Do1_01_a00001G00897:protein_coding;Upstream:1148|Transcript:Do1_01_a00001G00895V1.1:protein_coding|Gene:Do1_01_a00001G00895:protein_coding
# Do1_01_a00001   4805921 4806232 line_27;Downstream:3127|Transcript:Do1_01_a00001G00897V1.1:protein_coding|Gene:Do1_01_a00001G00897:protein_coding;Exon:1:1:RETAINED|Transcript:Do1_01_a00001G00896V1.1:protein_coding|Gene:Do1_01_a00001G00896:protein_coding;Intergenic:Do1_01_a00001G00896-Do1_01_a00001G00897;Upstream:4748|Transcript:Do1_01_a00001G00894V1.1:protein_coding|Gene:Do1_01_a00001G00894:protein_coding;Upstream:763|Transcript:Do1_01_a00001G00895V1.1:protein_coding|Gene:Do1_01_a00001G00895:protein_coding;Downstream:0|Transcript:Do1_01_a00001G00896V1.1:protein_coding|Gene:Do1_01_a00001G00896:protein_coding
# Do1_01_a00001   4835883 4836073 line_28;Downstream:2671|Transcript:Do1_01_a00001G00902V1.1:protein_coding|Gene:Do1_01_a00001G00902:protein_coding;Downstream:4441|Transcript:Do1_01_a00001G00900V1.1:protein_coding|Gene:Do1_01_a00001G00900:protein_coding;Exon:2:2:RETAINED|Transcript:Do1_01_a00001G00901V1.1:protein_coding|Gene:Do1_01_a00001G00901:protein_coding
# Do1_01_a00001   4836083 4836353 line_29;Downstream:4641|Transcript:Do1_01_a00001G00900V1.1:protein_coding|Gene:Do1_01_a00001G00900:protein_coding;Downstream:2471|Transcript:Do1_01_a00001G00902V1.1:protein_coding|Gene:Do1_01_a00001G00902:protein_coding;Exon:2:2:RETAINED|Transcript:Do1_01_a00001G00901V1.1:protein_coding|Gene:Do1_01_a00001G00901:protein_coding
# Do1_01_a00001   5125657 5125884 line_30;Downstream:3010|Transcript:Do1_01_a00001G00955V1.1:protein_coding|Gene:Do1_01_a00001G00955:protein_coding;Upstream:2707|Transcript:Do1_01_a00001G00956V1.1:protein_coding|Gene:Do1_01_a00001G00956:protein_coding;Intergenic:Do1_


# extract the immediate feature types and count them

cd /home/celphin/scratch/Dryas/snpEff/Dryas_DMRs/
FILES=$(ls *.out)
echo ${FILES}

for file in ${FILES}; do
awk -F'[;:]' '{print $1 $2}' "$file" > "$file"1
echo 
echo "$file Upstream"
grep "Upstream" "$file"1 | wc -l 
echo "$file Downstream"
grep "Downstream" "$file"1 | wc -l 
echo "$file Gene"
grep "Gene" "$file"1 | wc -l 
echo "$file Intergenic"
grep "Intergenic" "$file"1 | wc -l 
echo "$file Intron"
grep "Intron" "$file"1 | wc -l 
echo "$file Exon"
grep "Exon" "$file"1 | wc -l 
done

Alaska_W_C.bed.out Upstream
282
Alaska_W_C.bed.out Downstream
244
Alaska_W_C.bed.out Gene
35
Alaska_W_C.bed.out Intergenic
212
Alaska_W_C.bed.out Intron
11
Alaska_W_C.bed.out Exon
80

intersect_SE_L_H_Wild_L_H.bed.out Upstream
116
intersect_SE_L_H_Wild_L_H.bed.out Downstream
83
intersect_SE_L_H_Wild_L_H.bed.out Gene
16
intersect_SE_L_H_Wild_L_H.bed.out Intergenic
74
intersect_SE_L_H_Wild_L_H.bed.out Intron
2
intersect_SE_L_H_Wild_L_H.bed.out Exon
48

intersect_SE_W_C_P_W_C.bed.out Upstream
48
intersect_SE_W_C_P_W_C.bed.out Downstream
39
intersect_SE_W_C_P_W_C.bed.out Gene
4
intersect_SE_W_C_P_W_C.bed.out Intergenic
41
intersect_SE_W_C_P_W_C.bed.out Intron
1
intersect_SE_W_C_P_W_C.bed.out Exon
20

intersect_SE_W_C_SE_L_H.bed.out Upstream
34
intersect_SE_W_C_SE_L_H.bed.out Downstream
26
intersect_SE_W_C_SE_L_H.bed.out Gene
8
intersect_SE_W_C_SE_L_H.bed.out Intergenic
25
intersect_SE_W_C_SE_L_H.bed.out Intron
1
intersect_SE_W_C_SE_L_H.bed.out Exon
15

intersect_SE_W_C_Wild_W_C.bed.out Upstream
24
intersect_SE_W_C_Wild_W_C.bed.out Downstream
23
intersect_SE_W_C_Wild_W_C.bed.out Gene
4
intersect_SE_W_C_Wild_W_C.bed.out Intergenic
24
intersect_SE_W_C_Wild_W_C.bed.out Intron
0
intersect_SE_W_C_Wild_W_C.bed.out Exon
10

intersect_Wild_W_C_Mat_Sen.bed.out Upstream
19
intersect_Wild_W_C_Mat_Sen.bed.out Downstream
23
intersect_Wild_W_C_Mat_Sen.bed.out Gene
3
intersect_Wild_W_C_Mat_Sen.bed.out Intergenic
16
intersect_Wild_W_C_Mat_Sen.bed.out Intron
0
intersect_Wild_W_C_Mat_Sen.bed.out Exon
9

Mat_Sen.bed.out Upstream
1468
Mat_Sen.bed.out Downstream
1267
Mat_Sen.bed.out Gene
177
Mat_Sen.bed.out Intergenic
1496
Mat_Sen.bed.out Intron
133
Mat_Sen.bed.out Exon
436

Nunavut_W_C.bed.out Upstream
104
Nunavut_W_C.bed.out Downstream
71
Nunavut_W_C.bed.out Gene
16
Nunavut_W_C.bed.out Intergenic
89
Nunavut_W_C.bed.out Intron
3
Nunavut_W_C.bed.out Exon
40

Parent_W_C.bed.out Upstream
138
Parent_W_C.bed.out Downstream
124
Parent_W_C.bed.out Gene
14
Parent_W_C.bed.out Intergenic
110
Parent_W_C.bed.out Intron
12
Parent_W_C.bed.out Exon
51

SE_L_H.bed.out Upstream
142
SE_L_H.bed.out Downstream
105
SE_L_H.bed.out Gene
18
SE_L_H.bed.out Intergenic
94
SE_L_H.bed.out Intron
6
SE_L_H.bed.out Exon
62

SE_W_C.bed.out Upstream
277
SE_W_C.bed.out Downstream
224
SE_W_C.bed.out Gene
36
SE_W_C.bed.out Intergenic
226
SE_W_C.bed.out Intron
12
SE_W_C.bed.out Exon
114

Svalbard_W_C.bed.out Upstream
115
Svalbard_W_C.bed.out Downstream
92
Svalbard_W_C.bed.out Gene
18
Svalbard_W_C.bed.out Intergenic
94
Svalbard_W_C.bed.out Intron
7
Svalbard_W_C.bed.out Exon
33

Sweden_W_C.bed.out Upstream
521
Sweden_W_C.bed.out Downstream
475
Sweden_W_C.bed.out Gene
88
Sweden_W_C.bed.out Intergenic
442
Sweden_W_C.bed.out Intron
21
Sweden_W_C.bed.out Exon
190

total_subtract_SE_W_C_P_W_C.bed.out Upstream
229
total_subtract_SE_W_C_P_W_C.bed.out Downstream
185
total_subtract_SE_W_C_P_W_C.bed.out Gene
33
total_subtract_SE_W_C_P_W_C.bed.out Intergenic
185
total_subtract_SE_W_C_P_W_C.bed.out Intron
11
total_subtract_SE_W_C_P_W_C.bed.out Exon
94

total_subtract_SE_W_C_SE_L_H.bed.out Upstream
243
total_subtract_SE_W_C_SE_L_H.bed.out Downstream
198
total_subtract_SE_W_C_SE_L_H.bed.out Gene
29
total_subtract_SE_W_C_SE_L_H.bed.out Intergenic
201
total_subtract_SE_W_C_SE_L_H.bed.out Intron
11
total_subtract_SE_W_C_SE_L_H.bed.out Exon
99

total_subtract_W_C_Mat_Sen.bed.out Upstream
101
total_subtract_W_C_Mat_Sen.bed.out Downstream
68
total_subtract_W_C_Mat_Sen.bed.out Gene
12
total_subtract_W_C_Mat_Sen.bed.out Intergenic
96
total_subtract_W_C_Mat_Sen.bed.out Intron
4
total_subtract_W_C_Mat_Sen.bed.out Exon
35

Wild_Lat_L_H.bed.out Upstream
4491
Wild_Lat_L_H.bed.out Downstream
4024
Wild_Lat_L_H.bed.out Gene
591
Wild_Lat_L_H.bed.out Intergenic
4067
Wild_Lat_L_H.bed.out Intron
373
Wild_Lat_L_H.bed.out Exon
1427

Wild_W_C.bed.out Upstream
120
Wild_W_C.bed.out Downstream
91
Wild_W_C.bed.out Gene
14
Wild_W_C.bed.out Intergenic
112
Wild_W_C.bed.out Intron
4
Wild_W_C.bed.out Exon
44




#-----------------------------
# Example annotating variants in vcf file
java -Xmx8g -jar snpEff.jar OldDoct examples/test.chr22.vcf > test.chr22.ann.vcf

#--------------------------

# Look at transposons
# https://github.com/RobertsLab/project-gigas-oa-meth/blob/master/code/10-Genomic-Location-of-DML.ipynb 

module load StdEnv/2023 bedtools/2.31.0

cd /home/celphin/scratch/Dryas/snpEff/Dryas_DMRs/
FILES=$(ls *.bed)
echo ${FILES}

for file in ${FILES}; do
intersectBed \
-u \
-a ${file} \
-b /home/celphin/scratch/Dryas/snpEff/data/OldDoct/OldDoct.genome.fa.mod.EDTA.TEanno.gff3 \
> ${file}-All-TE.bed

#head ${file}-All-TE.bed
wc -l ${file}-All-TE.bed

done

# 330 Alaska_W_C.bed-All-TE.bed
# 111 intersect_SE_L_H_Wild_L_H.bed-All-TE.bed
# 53 intersect_SE_W_C_P_W_C.bed-All-TE.bed
# 29 intersect_SE_W_C_SE_L_H.bed-All-TE.bed
# 30 intersect_SE_W_C_Wild_W_C.bed-All-TE.bed
# 40 intersect_Wild_W_C_Mat_Sen.bed-All-TE.bed
# 3406 Mat_Sen.bed-All-TE.bed
# 139 Nunavut_W_C.bed-All-TE.bed
# 176 Parent_W_C.bed-All-TE.bed
# 148 SE_L_H.bed-All-TE.bed
# 225 SE_W_C.bed-All-TE.bed
# 136 Svalbard_W_C.bed-All-TE.bed
# 581 Sweden_W_C.bed-All-TE.bed
# 172 total_subtract_SE_W_C_P_W_C.bed-All-TE.bed
# 196 total_subtract_SE_W_C_SE_L_H.bed-All-TE.bed
# 133 total_subtract_W_C_Mat_Sen.bed-All-TE.bed
# 5429 Wild_Lat_L_H.bed-All-TE.bed
# 173 Wild_W_C.bed-All-TE.bed

#---------------------
# get overall totals

wc -l *.bedGraph

       # 866 Alaska_W_C.bedGraph
       # 338 intersect_SE_L_H_Wild_L_H.bedGraph
       # 152 intersect_SE_W_C_P_W_C.bedGraph
       # 108 intersect_SE_W_C_SE_L_H.bedGraph
        # 84 intersect_SE_W_C_Wild_W_C.bedGraph
        # 69 intersect_Wild_W_C_Mat_Sen.bedGraph
      # 4989 Mat_Sen.bedGraph
       # 322 Nunavut_W_C.bedGraph
       # 449 Parent_W_C.bedGraph
       # 426 SE_L_H.bedGraph
       # 888 SE_W_C.bedGraph
       # 358 Svalbard_W_C.bedGraph
      # 1737 Sweden_W_C.bedGraph
       # 736 total_subtract_SE_W_C_P_W_C.bedGraph
       # 780 total_subtract_SE_W_C_SE_L_H.bedGraph
       # 316 total_subtract_W_C_Mat_Sen.bedGraph
     # 14986 Wild_Lat_L_H.bedGraph
       # 385 Wild_W_C.bedGraph

#-------------------------
# could subet the TEs by type and rerun for each type

awk '{print $3}' /home/celphin/scratch/Dryas/snpEff/data/OldDoct/OldDoct.genome.fa.mod.EDTA.TEanno.gff3 | sort | uniq

# CACTA_TIR_transposon
# Copia_LTR_retrotransposon
# Gypsy_LTR_retrotransposon
# hAT_TIR_transposon
# helitron
# identity
# long_terminal_repeat
# LTR_retrotransposon
# Mutator_TIR_transposon
# Nov
# PIF_Harbinger_TIR_transposon
# repeat_region
# sequence_ontology
# site
# target_site_duplication
# Tc1_Mariner_TIR_transposon

cd /home/celphin/scratch/Dryas/snpEff/Dryas_DMRs/
# Extract unique values from the third column (or any other column)
unique_values=$(awk '{print $3}' /home/celphin/scratch/Dryas/snpEff/data/OldDoct/OldDoct.genome.fa.mod.EDTA.TEanno.gff3 | sort | uniq)

# Loop through each unique value and extract lines to a new file
for value in $unique_values; do
    grep "$value" /home/celphin/scratch/Dryas/snpEff/data/OldDoct/OldDoct.genome.fa.mod.EDTA.TEanno.gff3 > "${value}.gff3"
    echo "Extracted lines for $value into ${value}.gff3"
done

mkdir gff3
mv *.gff3 gff3/

#-----------------
module load StdEnv/2023 bedtools/2.31.0

cd /home/celphin/scratch/Dryas/snpEff/Dryas_DMRs/
FILES=$(ls *.bed)
echo ${FILES}

for value in $unique_values; do
for file in ${FILES}; do
intersectBed \
-u \
-a ${file} \
-b ./gff3/${value}.gff3 \
> ${file}-${value}.bed

wc -l ${file}-${value}.bed

done
done


35 Alaska_W_C.bed-CACTA_TIR_transposon.bed
35 Alaska_W_C.bed-All-TE.bed-CACTA_TIR_transposon.bed
16 intersect_SE_L_H_Wild_L_H.bed-CACTA_TIR_transposon.bed
16 intersect_SE_L_H_Wild_L_H.bed-All-TE.bed-CACTA_TIR_transposon.bed
5 intersect_SE_W_C_P_W_C.bed-CACTA_TIR_transposon.bed
5 intersect_SE_W_C_P_W_C.bed-All-TE.bed-CACTA_TIR_transposon.bed
6 intersect_SE_W_C_SE_L_H.bed-CACTA_TIR_transposon.bed
6 intersect_SE_W_C_SE_L_H.bed-All-TE.bed-CACTA_TIR_transposon.bed
1 intersect_SE_W_C_Wild_W_C.bed-CACTA_TIR_transposon.bed
1 intersect_SE_W_C_Wild_W_C.bed-All-TE.bed-CACTA_TIR_transposon.bed
2 intersect_Wild_W_C_Mat_Sen.bed-CACTA_TIR_transposon.bed
2 intersect_Wild_W_C_Mat_Sen.bed-All-TE.bed-CACTA_TIR_transposon.bed
451 Mat_Sen.bed-CACTA_TIR_transposon.bed
451 Mat_Sen.bed-All-TE.bed-CACTA_TIR_transposon.bed
16 Nunavut_W_C.bed-CACTA_TIR_transposon.bed
16 Nunavut_W_C.bed-All-TE.bed-CACTA_TIR_transposon.bed
12 Parent_W_C.bed-CACTA_TIR_transposon.bed
12 Parent_W_C.bed-All-TE.bed-CACTA_TIR_transposon.bed
16 SE_L_H.bed-CACTA_TIR_transposon.bed
16 SE_L_H.bed-All-TE.bed-CACTA_TIR_transposon.bed
26 SE_W_C.bed-CACTA_TIR_transposon.bed
26 SE_W_C.bed-All-TE.bed-CACTA_TIR_transposon.bed
13 Svalbard_W_C.bed-CACTA_TIR_transposon.bed
13 Svalbard_W_C.bed-All-TE.bed-CACTA_TIR_transposon.bed
39 Sweden_W_C.bed-CACTA_TIR_transposon.bed
39 Sweden_W_C.bed-All-TE.bed-CACTA_TIR_transposon.bed
21 total_subtract_SE_W_C_P_W_C.bed-CACTA_TIR_transposon.bed
21 total_subtract_SE_W_C_P_W_C.bed-All-TE.bed-CACTA_TIR_transposon.bed
20 total_subtract_SE_W_C_SE_L_H.bed-CACTA_TIR_transposon.bed
20 total_subtract_SE_W_C_SE_L_H.bed-All-TE.bed-CACTA_TIR_transposon.bed
15 total_subtract_W_C_Mat_Sen.bed-CACTA_TIR_transposon.bed
15 total_subtract_W_C_Mat_Sen.bed-All-TE.bed-CACTA_TIR_transposon.bed
824 Wild_Lat_L_H.bed-CACTA_TIR_transposon.bed
824 Wild_Lat_L_H.bed-All-TE.bed-CACTA_TIR_transposon.bed
17 Wild_W_C.bed-CACTA_TIR_transposon.bed
17 Wild_W_C.bed-All-TE.bed-CACTA_TIR_transposon.bed


27 Alaska_W_C.bed-Copia_LTR_retrotransposon.bed
27 Alaska_W_C.bed-All-TE.bed-Copia_LTR_retrotransposon.bed
6 intersect_SE_L_H_Wild_L_H.bed-Copia_LTR_retrotransposon.bed
6 intersect_SE_L_H_Wild_L_H.bed-All-TE.bed-Copia_LTR_retrotransposon.bed
3 intersect_SE_W_C_P_W_C.bed-Copia_LTR_retrotransposon.bed
3 intersect_SE_W_C_P_W_C.bed-All-TE.bed-Copia_LTR_retrotransposon.bed
2 intersect_SE_W_C_SE_L_H.bed-Copia_LTR_retrotransposon.bed
2 intersect_SE_W_C_SE_L_H.bed-All-TE.bed-Copia_LTR_retrotransposon.bed
2 intersect_SE_W_C_Wild_W_C.bed-Copia_LTR_retrotransposon.bed
2 intersect_SE_W_C_Wild_W_C.bed-All-TE.bed-Copia_LTR_retrotransposon.bed
5 intersect_Wild_W_C_Mat_Sen.bed-Copia_LTR_retrotransposon.bed
5 intersect_Wild_W_C_Mat_Sen.bed-All-TE.bed-Copia_LTR_retrotransposon.bed
315 Mat_Sen.bed-Copia_LTR_retrotransposon.bed
315 Mat_Sen.bed-All-TE.bed-Copia_LTR_retrotransposon.bed
8 Nunavut_W_C.bed-Copia_LTR_retrotransposon.bed
8 Nunavut_W_C.bed-All-TE.bed-Copia_LTR_retrotransposon.bed
12 Parent_W_C.bed-Copia_LTR_retrotransposon.bed
12 Parent_W_C.bed-All-TE.bed-Copia_LTR_retrotransposon.bed
10 SE_L_H.bed-Copia_LTR_retrotransposon.bed
10 SE_L_H.bed-All-TE.bed-Copia_LTR_retrotransposon.bed
20 SE_W_C.bed-Copia_LTR_retrotransposon.bed
20 SE_W_C.bed-All-TE.bed-Copia_LTR_retrotransposon.bed
13 Svalbard_W_C.bed-Copia_LTR_retrotransposon.bed
13 Svalbard_W_C.bed-All-TE.bed-Copia_LTR_retrotransposon.bed
45 Sweden_W_C.bed-Copia_LTR_retrotransposon.bed
45 Sweden_W_C.bed-All-TE.bed-Copia_LTR_retrotransposon.bed
17 total_subtract_SE_W_C_P_W_C.bed-Copia_LTR_retrotransposon.bed
17 total_subtract_SE_W_C_P_W_C.bed-All-TE.bed-Copia_LTR_retrotransposon.bed
18 total_subtract_SE_W_C_SE_L_H.bed-Copia_LTR_retrotransposon.bed
18 total_subtract_SE_W_C_SE_L_H.bed-All-TE.bed-Copia_LTR_retrotransposon.bed
12 total_subtract_W_C_Mat_Sen.bed-Copia_LTR_retrotransposon.bed
12 total_subtract_W_C_Mat_Sen.bed-All-TE.bed-Copia_LTR_retrotransposon.bed
479 Wild_Lat_L_H.bed-Copia_LTR_retrotransposon.bed
479 Wild_Lat_L_H.bed-All-TE.bed-Copia_LTR_retrotransposon.bed
17 Wild_W_C.bed-Copia_LTR_retrotransposon.bed
17 Wild_W_C.bed-All-TE.bed-Copia_LTR_retrotransposon.bed


35 Alaska_W_C.bed-Gypsy_LTR_retrotransposon.bed
35 Alaska_W_C.bed-All-TE.bed-Gypsy_LTR_retrotransposon.bed
15 intersect_SE_L_H_Wild_L_H.bed-Gypsy_LTR_retrotransposon.bed
15 intersect_SE_L_H_Wild_L_H.bed-All-TE.bed-Gypsy_LTR_retrotransposon.bed
2 intersect_SE_W_C_P_W_C.bed-Gypsy_LTR_retrotransposon.bed
2 intersect_SE_W_C_P_W_C.bed-All-TE.bed-Gypsy_LTR_retrotransposon.bed
3 intersect_SE_W_C_SE_L_H.bed-Gypsy_LTR_retrotransposon.bed
3 intersect_SE_W_C_SE_L_H.bed-All-TE.bed-Gypsy_LTR_retrotransposon.bed
2 intersect_SE_W_C_Wild_W_C.bed-Gypsy_LTR_retrotransposon.bed
2 intersect_SE_W_C_Wild_W_C.bed-All-TE.bed-Gypsy_LTR_retrotransposon.bed
10 intersect_Wild_W_C_Mat_Sen.bed-Gypsy_LTR_retrotransposon.bed
10 intersect_Wild_W_C_Mat_Sen.bed-All-TE.bed-Gypsy_LTR_retrotransposon.bed
899 Mat_Sen.bed-Gypsy_LTR_retrotransposon.bed
899 Mat_Sen.bed-All-TE.bed-Gypsy_LTR_retrotransposon.bed
16 Nunavut_W_C.bed-Gypsy_LTR_retrotransposon.bed
16 Nunavut_W_C.bed-All-TE.bed-Gypsy_LTR_retrotransposon.bed
17 Parent_W_C.bed-Gypsy_LTR_retrotransposon.bed
17 Parent_W_C.bed-All-TE.bed-Gypsy_LTR_retrotransposon.bed
24 SE_L_H.bed-Gypsy_LTR_retrotransposon.bed
24 SE_L_H.bed-All-TE.bed-Gypsy_LTR_retrotransposon.bed
20 SE_W_C.bed-Gypsy_LTR_retrotransposon.bed
20 SE_W_C.bed-All-TE.bed-Gypsy_LTR_retrotransposon.bed
17 Svalbard_W_C.bed-Gypsy_LTR_retrotransposon.bed
17 Svalbard_W_C.bed-All-TE.bed-Gypsy_LTR_retrotransposon.bed
89 Sweden_W_C.bed-Gypsy_LTR_retrotransposon.bed
89 Sweden_W_C.bed-All-TE.bed-Gypsy_LTR_retrotransposon.bed
18 total_subtract_SE_W_C_P_W_C.bed-Gypsy_LTR_retrotransposon.bed
18 total_subtract_SE_W_C_P_W_C.bed-All-TE.bed-Gypsy_LTR_retrotransposon.bed
17 total_subtract_SE_W_C_SE_L_H.bed-Gypsy_LTR_retrotransposon.bed
17 total_subtract_SE_W_C_SE_L_H.bed-All-TE.bed-Gypsy_LTR_retrotransposon.bed
16 total_subtract_W_C_Mat_Sen.bed-Gypsy_LTR_retrotransposon.bed
16 total_subtract_W_C_Mat_Sen.bed-All-TE.bed-Gypsy_LTR_retrotransposon.bed
816 Wild_Lat_L_H.bed-Gypsy_LTR_retrotransposon.bed
816 Wild_Lat_L_H.bed-All-TE.bed-Gypsy_LTR_retrotransposon.bed
26 Wild_W_C.bed-Gypsy_LTR_retrotransposon.bed
26 Wild_W_C.bed-All-TE.bed-Gypsy_LTR_retrotransposon.bed


111 Alaska_W_C.bed-hAT_TIR_transposon.bed
111 Alaska_W_C.bed-All-TE.bed-hAT_TIR_transposon.bed
28 intersect_SE_L_H_Wild_L_H.bed-hAT_TIR_transposon.bed
28 intersect_SE_L_H_Wild_L_H.bed-All-TE.bed-hAT_TIR_transposon.bed
17 intersect_SE_W_C_P_W_C.bed-hAT_TIR_transposon.bed
17 intersect_SE_W_C_P_W_C.bed-All-TE.bed-hAT_TIR_transposon.bed
6 intersect_SE_W_C_SE_L_H.bed-hAT_TIR_transposon.bed
6 intersect_SE_W_C_SE_L_H.bed-All-TE.bed-hAT_TIR_transposon.bed
11 intersect_SE_W_C_Wild_W_C.bed-hAT_TIR_transposon.bed
11 intersect_SE_W_C_Wild_W_C.bed-All-TE.bed-hAT_TIR_transposon.bed
10 intersect_Wild_W_C_Mat_Sen.bed-hAT_TIR_transposon.bed
10 intersect_Wild_W_C_Mat_Sen.bed-All-TE.bed-hAT_TIR_transposon.bed
592 Mat_Sen.bed-hAT_TIR_transposon.bed
592 Mat_Sen.bed-All-TE.bed-hAT_TIR_transposon.bed
44 Nunavut_W_C.bed-hAT_TIR_transposon.bed
44 Nunavut_W_C.bed-All-TE.bed-hAT_TIR_transposon.bed
56 Parent_W_C.bed-hAT_TIR_transposon.bed
56 Parent_W_C.bed-All-TE.bed-hAT_TIR_transposon.bed
39 SE_L_H.bed-hAT_TIR_transposon.bed
39 SE_L_H.bed-All-TE.bed-hAT_TIR_transposon.bed
70 SE_W_C.bed-hAT_TIR_transposon.bed
70 SE_W_C.bed-All-TE.bed-hAT_TIR_transposon.bed
47 Svalbard_W_C.bed-hAT_TIR_transposon.bed
47 Svalbard_W_C.bed-All-TE.bed-hAT_TIR_transposon.bed
181 Sweden_W_C.bed-hAT_TIR_transposon.bed
181 Sweden_W_C.bed-All-TE.bed-hAT_TIR_transposon.bed
53 total_subtract_SE_W_C_P_W_C.bed-hAT_TIR_transposon.bed
53 total_subtract_SE_W_C_P_W_C.bed-All-TE.bed-hAT_TIR_transposon.bed
64 total_subtract_SE_W_C_SE_L_H.bed-hAT_TIR_transposon.bed
64 total_subtract_SE_W_C_SE_L_H.bed-All-TE.bed-hAT_TIR_transposon.bed
38 total_subtract_W_C_Mat_Sen.bed-hAT_TIR_transposon.bed
38 total_subtract_W_C_Mat_Sen.bed-All-TE.bed-hAT_TIR_transposon.bed
1237 Wild_Lat_L_H.bed-hAT_TIR_transposon.bed
1237 Wild_Lat_L_H.bed-All-TE.bed-hAT_TIR_transposon.bed
48 Wild_W_C.bed-hAT_TIR_transposon.bed
48 Wild_W_C.bed-All-TE.bed-hAT_TIR_transposon.bed


15 Alaska_W_C.bed-helitron.bed
15 Alaska_W_C.bed-All-TE.bed-helitron.bed
6 intersect_SE_L_H_Wild_L_H.bed-helitron.bed
6 intersect_SE_L_H_Wild_L_H.bed-All-TE.bed-helitron.bed
4 intersect_SE_W_C_P_W_C.bed-helitron.bed
4 intersect_SE_W_C_P_W_C.bed-All-TE.bed-helitron.bed
1 intersect_SE_W_C_SE_L_H.bed-helitron.bed
1 intersect_SE_W_C_SE_L_H.bed-All-TE.bed-helitron.bed
4 intersect_SE_W_C_Wild_W_C.bed-helitron.bed
4 intersect_SE_W_C_Wild_W_C.bed-All-TE.bed-helitron.bed
3 intersect_Wild_W_C_Mat_Sen.bed-helitron.bed
3 intersect_Wild_W_C_Mat_Sen.bed-All-TE.bed-helitron.bed
106 Mat_Sen.bed-helitron.bed
106 Mat_Sen.bed-All-TE.bed-helitron.bed
13 Nunavut_W_C.bed-helitron.bed
13 Nunavut_W_C.bed-All-TE.bed-helitron.bed
13 Parent_W_C.bed-helitron.bed
13 Parent_W_C.bed-All-TE.bed-helitron.bed
9 SE_L_H.bed-helitron.bed
9 SE_L_H.bed-All-TE.bed-helitron.bed
13 SE_W_C.bed-helitron.bed
13 SE_W_C.bed-All-TE.bed-helitron.bed
7 Svalbard_W_C.bed-helitron.bed
7 Svalbard_W_C.bed-All-TE.bed-helitron.bed
39 Sweden_W_C.bed-helitron.bed
39 Sweden_W_C.bed-All-TE.bed-helitron.bed
9 total_subtract_SE_W_C_P_W_C.bed-helitron.bed
9 total_subtract_SE_W_C_P_W_C.bed-All-TE.bed-helitron.bed
12 total_subtract_SE_W_C_SE_L_H.bed-helitron.bed
12 total_subtract_SE_W_C_SE_L_H.bed-All-TE.bed-helitron.bed
8 total_subtract_W_C_Mat_Sen.bed-helitron.bed
8 total_subtract_W_C_Mat_Sen.bed-All-TE.bed-helitron.bed
349 Wild_Lat_L_H.bed-helitron.bed
349 Wild_Lat_L_H.bed-All-TE.bed-helitron.bed
11 Wild_W_C.bed-helitron.bed
11 Wild_W_C.bed-All-TE.bed-helitron.bed


14 Alaska_W_C.bed-identity.bed
14 Alaska_W_C.bed-All-TE.bed-identity.bed
3 intersect_SE_L_H_Wild_L_H.bed-identity.bed
3 intersect_SE_L_H_Wild_L_H.bed-All-TE.bed-identity.bed
0 intersect_SE_W_C_P_W_C.bed-identity.bed
0 intersect_SE_W_C_P_W_C.bed-All-TE.bed-identity.bed
2 intersect_SE_W_C_SE_L_H.bed-identity.bed
2 intersect_SE_W_C_SE_L_H.bed-All-TE.bed-identity.bed
0 intersect_SE_W_C_Wild_W_C.bed-identity.bed
0 intersect_SE_W_C_Wild_W_C.bed-All-TE.bed-identity.bed
4 intersect_Wild_W_C_Mat_Sen.bed-identity.bed
4 intersect_Wild_W_C_Mat_Sen.bed-All-TE.bed-identity.bed
409 Mat_Sen.bed-identity.bed
409 Mat_Sen.bed-All-TE.bed-identity.bed
7 Nunavut_W_C.bed-identity.bed
7 Nunavut_W_C.bed-All-TE.bed-identity.bed
3 Parent_W_C.bed-identity.bed
3 Parent_W_C.bed-All-TE.bed-identity.bed
8 SE_L_H.bed-identity.bed
8 SE_L_H.bed-All-TE.bed-identity.bed
7 SE_W_C.bed-identity.bed
7 SE_W_C.bed-All-TE.bed-identity.bed
7 Svalbard_W_C.bed-identity.bed
7 Svalbard_W_C.bed-All-TE.bed-identity.bed
26 Sweden_W_C.bed-identity.bed
26 Sweden_W_C.bed-All-TE.bed-identity.bed
7 total_subtract_SE_W_C_P_W_C.bed-identity.bed
7 total_subtract_SE_W_C_P_W_C.bed-All-TE.bed-identity.bed
5 total_subtract_SE_W_C_SE_L_H.bed-identity.bed
5 total_subtract_SE_W_C_SE_L_H.bed-All-TE.bed-identity.bed
7 total_subtract_W_C_Mat_Sen.bed-identity.bed
7 total_subtract_W_C_Mat_Sen.bed-All-TE.bed-identity.bed
315 Wild_Lat_L_H.bed-identity.bed
315 Wild_Lat_L_H.bed-All-TE.bed-identity.bed
11 Wild_W_C.bed-identity.bed
11 Wild_W_C.bed-All-TE.bed-identity.bed

2 Alaska_W_C.bed-long_terminal_repeat.bed
2 Alaska_W_C.bed-All-TE.bed-long_terminal_repeat.bed
0 intersect_SE_L_H_Wild_L_H.bed-long_terminal_repeat.bed
0 intersect_SE_L_H_Wild_L_H.bed-All-TE.bed-long_terminal_repeat.bed
0 intersect_SE_W_C_P_W_C.bed-long_terminal_repeat.bed
0 intersect_SE_W_C_P_W_C.bed-All-TE.bed-long_terminal_repeat.bed
0 intersect_SE_W_C_SE_L_H.bed-long_terminal_repeat.bed
0 intersect_SE_W_C_SE_L_H.bed-All-TE.bed-long_terminal_repeat.bed
0 intersect_SE_W_C_Wild_W_C.bed-long_terminal_repeat.bed
0 intersect_SE_W_C_Wild_W_C.bed-All-TE.bed-long_terminal_repeat.bed
0 intersect_Wild_W_C_Mat_Sen.bed-long_terminal_repeat.bed
0 intersect_Wild_W_C_Mat_Sen.bed-All-TE.bed-long_terminal_repeat.bed
53 Mat_Sen.bed-long_terminal_repeat.bed
53 Mat_Sen.bed-All-TE.bed-long_terminal_repeat.bed
0 Nunavut_W_C.bed-long_terminal_repeat.bed
0 Nunavut_W_C.bed-All-TE.bed-long_terminal_repeat.bed
0 Parent_W_C.bed-long_terminal_repeat.bed
0 Parent_W_C.bed-All-TE.bed-long_terminal_repeat.bed
1 SE_L_H.bed-long_terminal_repeat.bed
1 SE_L_H.bed-All-TE.bed-long_terminal_repeat.bed
0 SE_W_C.bed-long_terminal_repeat.bed
0 SE_W_C.bed-All-TE.bed-long_terminal_repeat.bed
2 Svalbard_W_C.bed-long_terminal_repeat.bed
2 Svalbard_W_C.bed-All-TE.bed-long_terminal_repeat.bed
5 Sweden_W_C.bed-long_terminal_repeat.bed
5 Sweden_W_C.bed-All-TE.bed-long_terminal_repeat.bed
0 total_subtract_SE_W_C_P_W_C.bed-long_terminal_repeat.bed
0 total_subtract_SE_W_C_P_W_C.bed-All-TE.bed-long_terminal_repeat.bed
0 total_subtract_SE_W_C_SE_L_H.bed-long_terminal_repeat.bed
0 total_subtract_SE_W_C_SE_L_H.bed-All-TE.bed-long_terminal_repeat.bed
1 total_subtract_W_C_Mat_Sen.bed-long_terminal_repeat.bed
1 total_subtract_W_C_Mat_Sen.bed-All-TE.bed-long_terminal_repeat.bed
62 Wild_Lat_L_H.bed-long_terminal_repeat.bed
62 Wild_Lat_L_H.bed-All-TE.bed-long_terminal_repeat.bed
1 Wild_W_C.bed-long_terminal_repeat.bed
1 Wild_W_C.bed-All-TE.bed-long_terminal_repeat.bed


122 Alaska_W_C.bed-LTR_retrotransposon.bed
122 Alaska_W_C.bed-All-TE.bed-LTR_retrotransposon.bed
47 intersect_SE_L_H_Wild_L_H.bed-LTR_retrotransposon.bed
47 intersect_SE_L_H_Wild_L_H.bed-All-TE.bed-LTR_retrotransposon.bed
10 intersect_SE_W_C_P_W_C.bed-LTR_retrotransposon.bed
10 intersect_SE_W_C_P_W_C.bed-All-TE.bed-LTR_retrotransposon.bed
16 intersect_SE_W_C_SE_L_H.bed-LTR_retrotransposon.bed
16 intersect_SE_W_C_SE_L_H.bed-All-TE.bed-LTR_retrotransposon.bed
11 intersect_SE_W_C_Wild_W_C.bed-LTR_retrotransposon.bed
11 intersect_SE_W_C_Wild_W_C.bed-All-TE.bed-LTR_retrotransposon.bed
22 intersect_Wild_W_C_Mat_Sen.bed-LTR_retrotransposon.bed
22 intersect_Wild_W_C_Mat_Sen.bed-All-TE.bed-LTR_retrotransposon.bed
1872 Mat_Sen.bed-LTR_retrotransposon.bed
1872 Mat_Sen.bed-All-TE.bed-LTR_retrotransposon.bed
47 Nunavut_W_C.bed-LTR_retrotransposon.bed
47 Nunavut_W_C.bed-All-TE.bed-LTR_retrotransposon.bed
63 Parent_W_C.bed-LTR_retrotransposon.bed
63 Parent_W_C.bed-All-TE.bed-LTR_retrotransposon.bed
68 SE_L_H.bed-LTR_retrotransposon.bed
68 SE_L_H.bed-All-TE.bed-LTR_retrotransposon.bed
71 SE_W_C.bed-LTR_retrotransposon.bed
71 SE_W_C.bed-All-TE.bed-LTR_retrotransposon.bed
56 Svalbard_W_C.bed-LTR_retrotransposon.bed
56 Svalbard_W_C.bed-All-TE.bed-LTR_retrotransposon.bed
230 Sweden_W_C.bed-LTR_retrotransposon.bed
230 Sweden_W_C.bed-All-TE.bed-LTR_retrotransposon.bed
61 total_subtract_SE_W_C_P_W_C.bed-LTR_retrotransposon.bed
61 total_subtract_SE_W_C_P_W_C.bed-All-TE.bed-LTR_retrotransposon.bed
55 total_subtract_SE_W_C_SE_L_H.bed-LTR_retrotransposon.bed
55 total_subtract_SE_W_C_SE_L_H.bed-All-TE.bed-LTR_retrotransposon.bed
54 total_subtract_W_C_Mat_Sen.bed-LTR_retrotransposon.bed
54 total_subtract_W_C_Mat_Sen.bed-All-TE.bed-LTR_retrotransposon.bed
2217 Wild_Lat_L_H.bed-LTR_retrotransposon.bed
2217 Wild_Lat_L_H.bed-All-TE.bed-LTR_retrotransposon.bed
76 Wild_W_C.bed-LTR_retrotransposon.bed
76 Wild_W_C.bed-All-TE.bed-LTR_retrotransposon.bed
73 Alaska_W_C.bed-Mutator_TIR_transposon.bed
73 Alaska_W_C.bed-All-TE.bed-Mutator_TIR_transposon.bed
24 intersect_SE_L_H_Wild_L_H.bed-Mutator_TIR_transposon.bed
24 intersect_SE_L_H_Wild_L_H.bed-All-TE.bed-Mutator_TIR_transposon.bed
19 intersect_SE_W_C_P_W_C.bed-Mutator_TIR_transposon.bed
19 intersect_SE_W_C_P_W_C.bed-All-TE.bed-Mutator_TIR_transposon.bed
2 intersect_SE_W_C_SE_L_H.bed-Mutator_TIR_transposon.bed
2 intersect_SE_W_C_SE_L_H.bed-All-TE.bed-Mutator_TIR_transposon.bed
5 intersect_SE_W_C_Wild_W_C.bed-Mutator_TIR_transposon.bed
5 intersect_SE_W_C_Wild_W_C.bed-All-TE.bed-Mutator_TIR_transposon.bed
9 intersect_Wild_W_C_Mat_Sen.bed-Mutator_TIR_transposon.bed
9 intersect_Wild_W_C_Mat_Sen.bed-All-TE.bed-Mutator_TIR_transposon.bed
608 Mat_Sen.bed-Mutator_TIR_transposon.bed
608 Mat_Sen.bed-All-TE.bed-Mutator_TIR_transposon.bed
24 Nunavut_W_C.bed-Mutator_TIR_transposon.bed
24 Nunavut_W_C.bed-All-TE.bed-Mutator_TIR_transposon.bed
43 Parent_W_C.bed-Mutator_TIR_transposon.bed
43 Parent_W_C.bed-All-TE.bed-Mutator_TIR_transposon.bed
29 SE_L_H.bed-Mutator_TIR_transposon.bed
29 SE_L_H.bed-All-TE.bed-Mutator_TIR_transposon.bed
55 SE_W_C.bed-Mutator_TIR_transposon.bed
55 SE_W_C.bed-All-TE.bed-Mutator_TIR_transposon.bed
21 Svalbard_W_C.bed-Mutator_TIR_transposon.bed
21 Svalbard_W_C.bed-All-TE.bed-Mutator_TIR_transposon.bed
123 Sweden_W_C.bed-Mutator_TIR_transposon.bed
123 Sweden_W_C.bed-All-TE.bed-Mutator_TIR_transposon.bed
36 total_subtract_SE_W_C_P_W_C.bed-Mutator_TIR_transposon.bed
36 total_subtract_SE_W_C_P_W_C.bed-All-TE.bed-Mutator_TIR_transposon.bed
53 total_subtract_SE_W_C_SE_L_H.bed-Mutator_TIR_transposon.bed
53 total_subtract_SE_W_C_SE_L_H.bed-All-TE.bed-Mutator_TIR_transposon.bed
23 total_subtract_W_C_Mat_Sen.bed-Mutator_TIR_transposon.bed
23 total_subtract_W_C_Mat_Sen.bed-All-TE.bed-Mutator_TIR_transposon.bed
1192 Wild_Lat_L_H.bed-Mutator_TIR_transposon.bed
1192 Wild_Lat_L_H.bed-All-TE.bed-Mutator_TIR_transposon.bed
32 Wild_W_C.bed-Mutator_TIR_transposon.bed
32 Wild_W_C.bed-All-TE.bed-Mutator_TIR_transposon.bed

17 Alaska_W_C.bed-PIF_Harbinger_TIR_transposon.bed
17 Alaska_W_C.bed-All-TE.bed-PIF_Harbinger_TIR_transposon.bed
2 intersect_SE_L_H_Wild_L_H.bed-PIF_Harbinger_TIR_transposon.bed
2 intersect_SE_L_H_Wild_L_H.bed-All-TE.bed-PIF_Harbinger_TIR_transposon.bed
2 intersect_SE_W_C_P_W_C.bed-PIF_Harbinger_TIR_transposon.bed
2 intersect_SE_W_C_P_W_C.bed-All-TE.bed-PIF_Harbinger_TIR_transposon.bed
2 intersect_SE_W_C_SE_L_H.bed-PIF_Harbinger_TIR_transposon.bed
2 intersect_SE_W_C_SE_L_H.bed-All-TE.bed-PIF_Harbinger_TIR_transposon.bed
2 intersect_SE_W_C_Wild_W_C.bed-PIF_Harbinger_TIR_transposon.bed
2 intersect_SE_W_C_Wild_W_C.bed-All-TE.bed-PIF_Harbinger_TIR_transposon.bed
1 intersect_Wild_W_C_Mat_Sen.bed-PIF_Harbinger_TIR_transposon.bed
1 intersect_Wild_W_C_Mat_Sen.bed-All-TE.bed-PIF_Harbinger_TIR_transposon.bed
152 Mat_Sen.bed-PIF_Harbinger_TIR_transposon.bed
152 Mat_Sen.bed-All-TE.bed-PIF_Harbinger_TIR_transposon.bed
7 Nunavut_W_C.bed-PIF_Harbinger_TIR_transposon.bed
7 Nunavut_W_C.bed-All-TE.bed-PIF_Harbinger_TIR_transposon.bed
8 Parent_W_C.bed-PIF_Harbinger_TIR_transposon.bed
8 Parent_W_C.bed-All-TE.bed-PIF_Harbinger_TIR_transposon.bed
4 SE_L_H.bed-PIF_Harbinger_TIR_transposon.bed
4 SE_L_H.bed-All-TE.bed-PIF_Harbinger_TIR_transposon.bed
12 SE_W_C.bed-PIF_Harbinger_TIR_transposon.bed
12 SE_W_C.bed-All-TE.bed-PIF_Harbinger_TIR_transposon.bed
2 Svalbard_W_C.bed-PIF_Harbinger_TIR_transposon.bed
2 Svalbard_W_C.bed-All-TE.bed-PIF_Harbinger_TIR_transposon.bed
27 Sweden_W_C.bed-PIF_Harbinger_TIR_transposon.bed
27 Sweden_W_C.bed-All-TE.bed-PIF_Harbinger_TIR_transposon.bed
10 total_subtract_SE_W_C_P_W_C.bed-PIF_Harbinger_TIR_transposon.bed
10 total_subtract_SE_W_C_P_W_C.bed-All-TE.bed-PIF_Harbinger_TIR_transposon.bed
10 total_subtract_SE_W_C_SE_L_H.bed-PIF_Harbinger_TIR_transposon.bed
10 total_subtract_SE_W_C_SE_L_H.bed-All-TE.bed-PIF_Harbinger_TIR_transposon.bed
7 total_subtract_W_C_Mat_Sen.bed-PIF_Harbinger_TIR_transposon.bed
7 total_subtract_W_C_Mat_Sen.bed-All-TE.bed-PIF_Harbinger_TIR_transposon.bed
215 Wild_Lat_L_H.bed-PIF_Harbinger_TIR_transposon.bed
215 Wild_Lat_L_H.bed-All-TE.bed-PIF_Harbinger_TIR_transposon.bed
8 Wild_W_C.bed-PIF_Harbinger_TIR_transposon.bed
8 Wild_W_C.bed-All-TE.bed-PIF_Harbinger_TIR_transposon.bed


14 Alaska_W_C.bed-repeat_region.bed
14 Alaska_W_C.bed-All-TE.bed-repeat_region.bed
3 intersect_SE_L_H_Wild_L_H.bed-repeat_region.bed
3 intersect_SE_L_H_Wild_L_H.bed-All-TE.bed-repeat_region.bed
0 intersect_SE_W_C_P_W_C.bed-repeat_region.bed
0 intersect_SE_W_C_P_W_C.bed-All-TE.bed-repeat_region.bed
2 intersect_SE_W_C_SE_L_H.bed-repeat_region.bed
2 intersect_SE_W_C_SE_L_H.bed-All-TE.bed-repeat_region.bed
0 intersect_SE_W_C_Wild_W_C.bed-repeat_region.bed
0 intersect_SE_W_C_Wild_W_C.bed-All-TE.bed-repeat_region.bed
4 intersect_Wild_W_C_Mat_Sen.bed-repeat_region.bed
4 intersect_Wild_W_C_Mat_Sen.bed-All-TE.bed-repeat_region.bed
409 Mat_Sen.bed-repeat_region.bed
409 Mat_Sen.bed-All-TE.bed-repeat_region.bed
7 Nunavut_W_C.bed-repeat_region.bed
7 Nunavut_W_C.bed-All-TE.bed-repeat_region.bed
3 Parent_W_C.bed-repeat_region.bed
3 Parent_W_C.bed-All-TE.bed-repeat_region.bed
8 SE_L_H.bed-repeat_region.bed
8 SE_L_H.bed-All-TE.bed-repeat_region.bed
7 SE_W_C.bed-repeat_region.bed
7 SE_W_C.bed-All-TE.bed-repeat_region.bed
7 Svalbard_W_C.bed-repeat_region.bed
7 Svalbard_W_C.bed-All-TE.bed-repeat_region.bed
26 Sweden_W_C.bed-repeat_region.bed
26 Sweden_W_C.bed-All-TE.bed-repeat_region.bed
7 total_subtract_SE_W_C_P_W_C.bed-repeat_region.bed
7 total_subtract_SE_W_C_P_W_C.bed-All-TE.bed-repeat_region.bed
5 total_subtract_SE_W_C_SE_L_H.bed-repeat_region.bed
5 total_subtract_SE_W_C_SE_L_H.bed-All-TE.bed-repeat_region.bed
7 total_subtract_W_C_Mat_Sen.bed-repeat_region.bed
7 total_subtract_W_C_Mat_Sen.bed-All-TE.bed-repeat_region.bed
315 Wild_Lat_L_H.bed-repeat_region.bed
315 Wild_Lat_L_H.bed-All-TE.bed-repeat_region.bed
11 Wild_W_C.bed-repeat_region.bed
11 Wild_W_C.bed-All-TE.bed-repeat_region.bed



1 Alaska_W_C.bed-site.bed
1 Alaska_W_C.bed-All-TE.bed-site.bed
0 intersect_SE_L_H_Wild_L_H.bed-site.bed
0 intersect_SE_L_H_Wild_L_H.bed-All-TE.bed-site.bed
0 intersect_SE_W_C_P_W_C.bed-site.bed
0 intersect_SE_W_C_P_W_C.bed-All-TE.bed-site.bed
0 intersect_SE_W_C_SE_L_H.bed-site.bed
0 intersect_SE_W_C_SE_L_H.bed-All-TE.bed-site.bed
0 intersect_SE_W_C_Wild_W_C.bed-site.bed
0 intersect_SE_W_C_Wild_W_C.bed-All-TE.bed-site.bed
0 intersect_Wild_W_C_Mat_Sen.bed-site.bed
0 intersect_Wild_W_C_Mat_Sen.bed-All-TE.bed-site.bed
7 Mat_Sen.bed-site.bed
7 Mat_Sen.bed-All-TE.bed-site.bed
0 Nunavut_W_C.bed-site.bed
0 Nunavut_W_C.bed-All-TE.bed-site.bed
0 Parent_W_C.bed-site.bed
0 Parent_W_C.bed-All-TE.bed-site.bed
0 SE_L_H.bed-site.bed
0 SE_L_H.bed-All-TE.bed-site.bed
0 SE_W_C.bed-site.bed
0 SE_W_C.bed-All-TE.bed-site.bed
2 Svalbard_W_C.bed-site.bed
2 Svalbard_W_C.bed-All-TE.bed-site.bed
2 Sweden_W_C.bed-site.bed
2 Sweden_W_C.bed-All-TE.bed-site.bed
0 total_subtract_SE_W_C_P_W_C.bed-site.bed
0 total_subtract_SE_W_C_P_W_C.bed-All-TE.bed-site.bed
0 total_subtract_SE_W_C_SE_L_H.bed-site.bed
0 total_subtract_SE_W_C_SE_L_H.bed-All-TE.bed-site.bed
1 total_subtract_W_C_Mat_Sen.bed-site.bed
1 total_subtract_W_C_Mat_Sen.bed-All-TE.bed-site.bed
12 Wild_Lat_L_H.bed-site.bed
12 Wild_Lat_L_H.bed-All-TE.bed-site.bed
1 Wild_W_C.bed-site.bed
1 Wild_W_C.bed-All-TE.bed-site.bed


1 Alaska_W_C.bed-target_site_duplication.bed
1 Alaska_W_C.bed-All-TE.bed-target_site_duplication.bed
0 intersect_SE_L_H_Wild_L_H.bed-target_site_duplication.bed
0 intersect_SE_L_H_Wild_L_H.bed-All-TE.bed-target_site_duplication.bed
0 intersect_SE_W_C_P_W_C.bed-target_site_duplication.bed
0 intersect_SE_W_C_P_W_C.bed-All-TE.bed-target_site_duplication.bed
0 intersect_SE_W_C_SE_L_H.bed-target_site_duplication.bed
0 intersect_SE_W_C_SE_L_H.bed-All-TE.bed-target_site_duplication.bed
0 intersect_SE_W_C_Wild_W_C.bed-target_site_duplication.bed
0 intersect_SE_W_C_Wild_W_C.bed-All-TE.bed-target_site_duplication.bed
0 intersect_Wild_W_C_Mat_Sen.bed-target_site_duplication.bed
0 intersect_Wild_W_C_Mat_Sen.bed-All-TE.bed-target_site_duplication.bed
7 Mat_Sen.bed-target_site_duplication.bed
7 Mat_Sen.bed-All-TE.bed-target_site_duplication.bed
0 Nunavut_W_C.bed-target_site_duplication.bed
0 Nunavut_W_C.bed-All-TE.bed-target_site_duplication.bed
0 Parent_W_C.bed-target_site_duplication.bed
0 Parent_W_C.bed-All-TE.bed-target_site_duplication.bed
0 SE_L_H.bed-target_site_duplication.bed
0 SE_L_H.bed-All-TE.bed-target_site_duplication.bed
0 SE_W_C.bed-target_site_duplication.bed
0 SE_W_C.bed-All-TE.bed-target_site_duplication.bed
2 Svalbard_W_C.bed-target_site_duplication.bed
2 Svalbard_W_C.bed-All-TE.bed-target_site_duplication.bed
2 Sweden_W_C.bed-target_site_duplication.bed
2 Sweden_W_C.bed-All-TE.bed-target_site_duplication.bed
0 total_subtract_SE_W_C_P_W_C.bed-target_site_duplication.bed
0 total_subtract_SE_W_C_P_W_C.bed-All-TE.bed-target_site_duplication.bed
0 total_subtract_SE_W_C_SE_L_H.bed-target_site_duplication.bed
0 total_subtract_SE_W_C_SE_L_H.bed-All-TE.bed-target_site_duplication.bed
1 total_subtract_W_C_Mat_Sen.bed-target_site_duplication.bed
1 total_subtract_W_C_Mat_Sen.bed-All-TE.bed-target_site_duplication.bed
12 Wild_Lat_L_H.bed-target_site_duplication.bed
12 Wild_Lat_L_H.bed-All-TE.bed-target_site_duplication.bed
1 Wild_W_C.bed-target_site_duplication.bed
1 Wild_W_C.bed-All-TE.bed-target_site_duplication.bed


1 Alaska_W_C.bed-Tc1_Mariner_TIR_transposon.bed
1 Alaska_W_C.bed-All-TE.bed-Tc1_Mariner_TIR_transposon.bed
0 intersect_SE_L_H_Wild_L_H.bed-Tc1_Mariner_TIR_transposon.bed
0 intersect_SE_L_H_Wild_L_H.bed-All-TE.bed-Tc1_Mariner_TIR_transposon.bed
0 intersect_SE_W_C_P_W_C.bed-Tc1_Mariner_TIR_transposon.bed
0 intersect_SE_W_C_P_W_C.bed-All-TE.bed-Tc1_Mariner_TIR_transposon.bed
0 intersect_SE_W_C_SE_L_H.bed-Tc1_Mariner_TIR_transposon.bed
0 intersect_SE_W_C_SE_L_H.bed-All-TE.bed-Tc1_Mariner_TIR_transposon.bed
0 intersect_SE_W_C_Wild_W_C.bed-Tc1_Mariner_TIR_transposon.bed
0 intersect_SE_W_C_Wild_W_C.bed-All-TE.bed-Tc1_Mariner_TIR_transposon.bed
0 intersect_Wild_W_C_Mat_Sen.bed-Tc1_Mariner_TIR_transposon.bed
0 intersect_Wild_W_C_Mat_Sen.bed-All-TE.bed-Tc1_Mariner_TIR_transposon.bed
1 Mat_Sen.bed-Tc1_Mariner_TIR_transposon.bed
1 Mat_Sen.bed-All-TE.bed-Tc1_Mariner_TIR_transposon.bed
0 Nunavut_W_C.bed-Tc1_Mariner_TIR_transposon.bed
0 Nunavut_W_C.bed-All-TE.bed-Tc1_Mariner_TIR_transposon.bed
1 Parent_W_C.bed-Tc1_Mariner_TIR_transposon.bed
1 Parent_W_C.bed-All-TE.bed-Tc1_Mariner_TIR_transposon.bed
0 SE_L_H.bed-Tc1_Mariner_TIR_transposon.bed
0 SE_L_H.bed-All-TE.bed-Tc1_Mariner_TIR_transposon.bed
0 SE_W_C.bed-Tc1_Mariner_TIR_transposon.bed
0 SE_W_C.bed-All-TE.bed-Tc1_Mariner_TIR_transposon.bed
0 Svalbard_W_C.bed-Tc1_Mariner_TIR_transposon.bed
0 Svalbard_W_C.bed-All-TE.bed-Tc1_Mariner_TIR_transposon.bed
1 Sweden_W_C.bed-Tc1_Mariner_TIR_transposon.bed
1 Sweden_W_C.bed-All-TE.bed-Tc1_Mariner_TIR_transposon.bed
0 total_subtract_SE_W_C_P_W_C.bed-Tc1_Mariner_TIR_transposon.bed
0 total_subtract_SE_W_C_P_W_C.bed-All-TE.bed-Tc1_Mariner_TIR_transposon.bed
0 total_subtract_SE_W_C_SE_L_H.bed-Tc1_Mariner_TIR_transposon.bed
0 total_subtract_SE_W_C_SE_L_H.bed-All-TE.bed-Tc1_Mariner_TIR_transposon.bed
0 total_subtract_W_C_Mat_Sen.bed-Tc1_Mariner_TIR_transposon.bed
0 total_subtract_W_C_Mat_Sen.bed-All-TE.bed-Tc1_Mariner_TIR_transposon.bed
9 Wild_Lat_L_H.bed-Tc1_Mariner_TIR_transposon.bed
9 Wild_Lat_L_H.bed-All-TE.bed-Tc1_Mariner_TIR_transposon.bed
0 Wild_W_C.bed-Tc1_Mariner_TIR_transposon.bed
0 Wild_W_C.bed-All-TE.bed-Tc1_Mariner_TIR_transposon.bed


#############################################
# Run for sites

tmux new-session -s snpEff
tmux attach-session -t snpEff

cd /home/celphin/scratch/Dryas/snpEff/Dryas_DMRs/sites

FILES=$(ls /home/celphin/scratch/Dryas/snpEff/Dryas_DMRs/sites/*.bedGraph)
echo ${FILES}

for file in ${FILES}; do
    # Run the command on the file and create the .bed file
    grep -v track "$file" | awk '{print $1 "\t" $2 "\t" $3}' > "${file%.bedGraph}.bed"
    echo "Processed $file"
done

# run SNPEff
module load StdEnv/2023 java/21.0.1
cd /home/celphin/scratch/Dryas/snpEff/Dryas_DMRs/sites/
FILES=$(ls /home/celphin/scratch/Dryas/snpEff/Dryas_DMRs/sites/*.bed)
echo ${FILES}

cd /home/celphin/scratch/Dryas/snpEff

for file in ${FILES}; do
java -Xmx8g -jar snpEff.jar -i bed OldDoct "$file" > "$file".out
echo "Processed $file"
done


# extract the immediate feature types and count them
cd /home/celphin/scratch/Dryas/snpEff/Dryas_DMRs/sites
FILES=$(ls *.out)
echo ${FILES}

for file in ${FILES}; do
awk -F'[;:]' '{print $1 $2}' "$file" > "$file"1
echo 
echo "$file Upstream"
grep "Upstream" "$file"1 | wc -l 
echo "$file Downstream"
grep "Downstream" "$file"1 | wc -l 
echo "$file Gene"
grep "Gene" "$file"1 | wc -l 
echo "$file Intergenic"
grep "Intergenic" "$file"1 | wc -l 
echo "$file Intron"
grep "Intron" "$file"1 | wc -l 
echo "$file Exon"
grep "Exon" "$file"1 | wc -l 
done

#--------------------------

# Look at transposons
# https://github.com/RobertsLab/project-gigas-oa-meth/blob/master/code/10-Genomic-Location-of-DML.ipynb 

module load StdEnv/2023 bedtools/2.31.0

cd /home/celphin/scratch/Dryas/snpEff/Dryas_DMRs/
FILES=$(ls *.bed)
echo ${FILES}

for file in ${FILES}; do
intersectBed \
-u \
-a ${file} \
-b /home/celphin/scratch/Dryas/snpEff/data/OldDoct/OldDoct.genome.fa.mod.EDTA.TEanno.gff3 \
> ${file}-All-TE.bed

#head ${file}-All-TE.bed
wc -l ${file}-All-TE.bed

done


#-------------------------
# could subet the TEs by type and rerun for each type

awk '{print $3}' /home/celphin/scratch/Dryas/snpEff/data/OldDoct/OldDoct.genome.fa.mod.EDTA.TEanno.gff3 | sort | uniq

# CACTA_TIR_transposon
# Copia_LTR_retrotransposon
# Gypsy_LTR_retrotransposon
# hAT_TIR_transposon
# helitron
# identity
# long_terminal_repeat
# LTR_retrotransposon
# Mutator_TIR_transposon
# Nov
# PIF_Harbinger_TIR_transposon
# repeat_region
# sequence_ontology
# site
# target_site_duplication
# Tc1_Mariner_TIR_transposon

cd /home/celphin/scratch/Dryas/snpEff/Dryas_DMRs/
# Extract unique values from the third column (or any other column)
unique_values=$(awk '{print $3}' /home/celphin/scratch/Dryas/snpEff/data/OldDoct/OldDoct.genome.fa.mod.EDTA.TEanno.gff3 | sort | uniq)

# Loop through each unique value and extract lines to a new file
for value in $unique_values; do
    grep "$value" /home/celphin/scratch/Dryas/snpEff/data/OldDoct/OldDoct.genome.fa.mod.EDTA.TEanno.gff3 > "${value}.gff3"
    echo "Extracted lines for $value into ${value}.gff3"
done

mkdir gff3
mv *.gff3 gff3/

#-----------------
module load StdEnv/2023 bedtools/2.31.0

cd /home/celphin/scratch/Dryas/snpEff/Dryas_DMRs/
FILES=$(ls *.bed)
echo ${FILES}

for value in $unique_values; do
for file in ${FILES}; do
intersectBed \
-u \
-a ${file} \
-b ./gff3/${value}.gff3 \
> ${file}-${value}.bed

wc -l ${file}-${value}.bed

done
done
