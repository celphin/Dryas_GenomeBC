#################################
# Exploring Dryas DMRs
# Nov 2023
###############################

cd /home/celphin/projects/def-rieseber/Dryas_shared_data/MS_Dryas_Merged_Data

# search and count for the various terms
grep "GO:" Gene_DMR_Total_GO_Merged_table.tsv | grep "Mat_Sen"  | awk '{print $2 " " $13 " " $14}' | sort | uniq -c | sort -bgr  > go_terms_Mat_Sen.txt
grep "GO:" Gene_DMR_Total_GO_Merged_table.tsv | grep "_W_C"  | awk '{print $2 " " $13 " " $14}'  | sort | uniq -c | sort -bgr > go_terms_Warming.txt
grep "GO:" Gene_DMR_Total_GO_Merged_table.tsv | grep "Wild_Lat_L_H"  | awk '{print $2 " " $13 " " $14}'  | sort | uniq -c | sort -bgr > go_terms_Latitude.txt
grep "GO:" Gene_DMR_Total_GO_Merged_table.tsv | grep "SE_W_C"  | awk '{print $2 " " $13 " " $14}'  | sort | uniq -c | sort -bgr > go_terms_Seedling.txt
grep "GO:" Gene_DMR_Total_GO_Merged_table.tsv | grep "Wild_W_C"  | awk '{print $2 " " $13 " " $14}'  | sort | uniq -c | sort -bgr > go_terms_Wild.txt
grep "GO:" Gene_DMR_Total_GO_Merged_table.tsv | grep "Parent_W_C"  | awk '{print $2 " " $13 " " $14}'  | sort | uniq -c | sort -bgr > go_terms_Parent.txt
grep "GO:" Gene_DMR_Total_GO_Merged_table.tsv | grep "Sweden_W_C"  | awk '{print $2 " " $13 " " $14}'  | sort | uniq -c | sort -bgr > go_terms_Sweden.txt
grep "GO:" Gene_DMR_Total_GO_Merged_table.tsv | grep "Alaska_W_C"  | awk '{print $2 " " $13 " " $14}'  | sort | uniq -c | sort -bgr > go_terms_Alaska.txt
grep "GO:" Gene_DMR_Total_GO_Merged_table.tsv | grep "Nunavut_W_C"  | awk '{print $2 " " $13 " " $14}'  | sort | uniq -c | sort -bgr > go_terms_Nunavut.txt
grep "GO:" Gene_DMR_Total_GO_Merged_table.tsv | grep "Svalbard_W_C"  | awk '{print $2 " " $13 " " $14}'  | sort | uniq -c | sort -bgr > go_terms_Svalbard.txt

##########################################
# compare to total counts	 
cd ./original_data
awk '{print $2}' dryas_goterm_file.tsv | sort | uniq -c | sort -bgr  > ../go_terms_Dryas.txt

##########################################
# Genes of interest: monooxygenase, defense, UDP-glycosyltransferase, carbohydrate, lipid, ATP, ADP, polysaccharide, oxidoreductase, membrane, structural,
# hydrolase, DNA, RNA, serine-type, polygalacturonase, catalytic, signal 
 
 awk '{s+=$1}END{print s}' go_terms_Dryas.txt
60939

awk '{s+=$1}END{print s}' go_terms_Sweden.txt
1230

awk '{s+=$1}END{print s}' go_terms_Latitude.txt
12374

awk '{s+=$1}END{print s}'  go_terms_Mat_Sen.txt
2079

awk '{s+=$1}END{print s}'  go_terms_Wild.txt
210

awk '{s+=$1}END{print s}'  go_terms_Seedling.txt
545

#-----------------------------
 grep GO:monooxygenase *.txt
 
go_terms_Alaska.txt:     84 Alaska_W_C GO:0004497,GO:0005506,GO:0016705,GO:0020037 GO:monooxygenase
go_terms_Latitude.txt:    723 Wild_Lat_L_H GO:0004497,GO:0005506,GO:0016705,GO:0020037 GO:monooxygenase
go_terms_Mat_Sen.txt:     45 Mat_Sen GO:0004497,GO:0005506,GO:0016705,GO:0020037 GO:monooxygenase
go_terms_Nunavut.txt:     16 Nunavut_W_C GO:0004497,GO:0005506,GO:0016705,GO:0020037 GO:monooxygenase
go_terms_Parent.txt:     28 Parent_W_C GO:0004497,GO:0005506,GO:0016705,GO:0020037 GO:monooxygenase
go_terms_Seedling.txt:     45 SE_W_C GO:0004497,GO:0005506,GO:0016705,GO:0020037 GO:monooxygenase
go_terms_Svalbard.txt:     15 Svalbard_W_C GO:0004497,GO:0005506,GO:0016705,GO:0020037 GO:monooxygenase
go_terms_Sweden.txt:    107 Sweden_W_C GO:0004497,GO:0005506,GO:0016705,GO:0020037 GO:monooxygenase
go_terms_Warming.txt:    107 Sweden_W_C GO:0004497,GO:0005506,GO:0016705,GO:0020037 GO:monooxygenase
go_terms_Warming.txt:     84 Alaska_W_C GO:0004497,GO:0005506,GO:0016705,GO:0020037 GO:monooxygenase
go_terms_Warming.txt:     45 SE_W_C GO:0004497,GO:0005506,GO:0016705,GO:0020037 GO:monooxygenase
go_terms_Warming.txt:     28 Parent_W_C GO:0004497,GO:0005506,GO:0016705,GO:0020037 GO:monooxygenase
go_terms_Warming.txt:     16 Nunavut_W_C GO:0004497,GO:0005506,GO:0016705,GO:0020037 GO:monooxygenase
go_terms_Warming.txt:     15 Svalbard_W_C GO:0004497,GO:0005506,GO:0016705,GO:0020037 GO:monooxygenase
go_terms_Warming.txt:      7 Wild_W_C GO:0004497,GO:0005506,GO:0016705,GO:0020037 GO:monooxygenase
go_terms_Wild.txt:      7 Wild_W_C GO:0004497,GO:0005506,GO:0016705,GO:0020037 GO:monooxygenase



grep GO:0004497 go_terms_Dryas.txt
2826 GO:0004497|GO:0005506|GO:0016705|GO:0020037

Null  2826/60939 = 0.04
Mat Sen 45/2079 = 0.02
Wild 7/210 = 0.03

Lat = 723/12374 = 0.06
Seed = 45/ 545 =  0.08

#-----------------------------
# try defense
# GO:0006952

grep defense *.txt
grep GO:0006952 *.txt
 
go_terms_Dryas.txt:    311 GO:0006952
go_terms_Latitude.txt:    500 Wild_Lat_L_H GO:0006952,GO:0043531 GO:defense
go_terms_Alaska.txt:     40 Alaska_W_C GO:0006952,GO:0043531 GO:defense
go_terms_Mat_Sen.txt:      20 Mat_Sen GO:0006952,GO:0043531 GO:defense
go_terms_Sweden.txt:     100 Sweden_W_C GO:0006952,GO:0043531 GO:defense
go_terms_Seedling.txt:     26 SE_W_C GO:0006952,GO:0043531 GO:defense


Null  311/60939 = 0.005
Mat Sen 8/2079 = 0.003
Wild 14/210 = 0.07
Lat = 500/12374 = 0.04
Seed = 22/ 545 =  0.04

#-------------------------------
# try UDP-glycosyltransferase

grep UDP-glycosyltransferase *.txt

go_terms_Alaska.txt:     16 Alaska_W_C GO:0008194 GO:UDP-glycosyltransferase
go_terms_Latitude.txt:    235 Wild_Lat_L_H GO:0008194 GO:UDP-glycosyltransferase
go_terms_Mat_Sen.txt:     28 Mat_Sen GO:0008194 GO:UDP-glycosyltransferase
go_terms_Nunavut.txt:     14 Nunavut_W_C GO:0008194 GO:UDP-glycosyltransferase
go_terms_Parent.txt:     13 Parent_W_C GO:0008194 GO:UDP-glycosyltransferase
go_terms_Seedling.txt:     16 SE_W_C GO:0008194 GO:UDP-glycosyltransferase
go_terms_Sweden.txt:     19 Sweden_W_C GO:0008194 GO:UDP-glycosyltransferase
go_terms_Warming.txt:     19 Sweden_W_C GO:0008194 GO:UDP-glycosyltransferase
go_terms_Warming.txt:     16 SE_W_C GO:0008194 GO:UDP-glycosyltransferase
go_terms_Warming.txt:     16 Alaska_W_C GO:0008194 GO:UDP-glycosyltransferase
go_terms_Warming.txt:     14 Nunavut_W_C GO:0008194 GO:UDP-glycosyltransferase
go_terms_Warming.txt:     13 Parent_W_C GO:0008194 GO:UDP-glycosyltransferase
go_terms_Warming.txt:      8 Wild_W_C GO:0008194 GO:UDP-glycosyltransferase
go_terms_Wild.txt:      8 Wild_W_C GO:0008194 GO:UDP-glycosyltransferase

grep GO:0008194 *.txt
go_terms_Dryas.txt:    217 GO:0008194


Null  217/60939 = 0.004
Mat Sen 28/2079 = 0.01
Wild 8/210 = 0.04
Lat = 235/12374 = 0.02
Seed = 16/ 545 =  0.03

#------------------------------
# try GO:0043457 regulation

grep GO:0043457 *.txt

go_terms_Dryas.txt:     27 GO:0043457

grep regulation *.txt

go_terms_Latitude.txt:     15 Wild_Lat_L_H GO:0043457 GO:regulation
go_terms_Mat_Sen.txt:     32 Mat_Sen GO:0043457 GO:regulation


#-----------------------------
# try GO:0030247 GO:polysaccharide

grep GO:0030247 *.txt
go_terms_Dryas.txt:     59 GO:0030247

grep GO:polysaccharide *.txt
go_terms_Alaska.txt:      3 Alaska_W_C GO:0030247 GO:polysaccharide
go_terms_Latitude.txt:    120 Wild_Lat_L_H GO:0030247 GO:polysaccharide
go_terms_Parent.txt:      7 Parent_W_C GO:0030247 GO:polysaccharide
go_terms_Seedling.txt:      9 SE_W_C GO:0030247 GO:polysaccharide
go_terms_Sweden.txt:     21 Sweden_W_C GO:0030247 GO:polysaccharide
go_terms_Warming.txt:     21 Sweden_W_C GO:0030247 GO:polysaccharide
go_terms_Warming.txt:      9 Wild_W_C GO:0030247 GO:polysaccharide
go_terms_Warming.txt:      9 SE_W_C GO:0030247 GO:polysaccharide
go_terms_Warming.txt:      7 Parent_W_C GO:0030247 GO:polysaccharide
go_terms_Warming.txt:      3 Alaska_W_C GO:0030247 GO:polysaccharide
go_terms_Wild.txt:      9 Wild_W_C GO:0030247 GO:polysaccharide


Null  59/60939 = 0.001
Mat Sen 0/2079 = 0

Wild 9/210 = 0.04
Lat = 120/12374 = 0.01
Seed = 9/ 545 =  0.02


####################################################
# examples
 
   304 Mat_Sen GO:0003676,GO:0004523 GO:nucleic
    265 Mat_Sen GO:0003676,GO:0015074 GO:nucleic
    157 Mat_Sen GO:0003676,GO:0004523,GO:0015074 GO:nucleic
    111 Mat_Sen GO:0000166,GO:0004812,GO:0004827,GO:0005524,GO:0005737,GO:0006418,GO:0006433 GO:nucleotide
     98 Mat_Sen GO:0005515 GO:protein
     66 Mat_Sen GO:0004672,GO:0005524,GO:0006468 GO:protein
     64 Mat_Sen GO:0003676 GO:nucleic
     45 Mat_Sen GO:0004497,GO:0005506,GO:0016705,GO:0020037 GO:monooxygenase
     40 Mat_Sen GO:0003676,GO:0008270 GO:nucleic
     32 Mat_Sen GO:0043457 GO:regulation
     28 Mat_Sen GO:0008194 GO:UDP-glycosyltransferase
     23 Mat_Sen GO:0046983 GO:protein
     17 Mat_Sen GO:0010073,GO:0048507 GO:meristem
     15 Mat_Sen GO:0004650,GO:0005975 GO:polygalacturonase
     14 Mat_Sen GO:0006355 GO:regulation
     14 Mat_Sen GO:0003677 GO:DNA
     13 Mat_Sen GO:0005315,GO:0006817,GO:0016020,GO:0022857,GO:0055085 GO:inorganic
     12 Mat_Sen GO:0016491 GO:oxidoreductase
     12 Mat_Sen GO:0000723,GO:0003678,GO:0006281 GO:telomere
     11 Mat_Sen GO:0005509 GO:calcium
     10 Mat_Sen GO:0015144,GO:0016020,GO:0022857,GO:0055085 GO:carbohydrate
     10 Mat_Sen GO:0008270 GO:zinc
     10 Mat_Sen GO:0004672,GO:0004674,GO:0005524,GO:0006468 GO:protein
     10 Mat_Sen GO:0003735,GO:0005840,GO:0006412 GO:structural
     10 Mat_Sen GO:0003690,GO:0006355 GO:double-stranded



     74 SE_W_C GO:0005515 GO:protein
     45 SE_W_C GO:0004497,GO:0005506,GO:0016705,GO:0020037 GO:monooxygenase
     31 SE_W_C GO:0004672,GO:0005524,GO:0006468 GO:protein
     16 SE_W_C GO:0008194 GO:UDP-glycosyltransferase
     11 SE_W_C GO:0006952,GO:0043531 GO:defense
     10 SE_W_C GO:0004672,GO:0005515,GO:0005524,GO:0006468 GO:protein
      9 SE_W_C GO:0030247 GO:polysaccharide
      8 SE_W_C GO:0005515,GO:0006952,GO:0043531 GO:protein
      8 SE_W_C GO:0003677,GO:0006355 GO:DNA
      7 SE_W_C GO:0043531 GO:ADP
      7 SE_W_C GO:0004672,GO:0004674,GO:0005524,GO:0006468 GO:protein
      7 SE_W_C GO:0004553,GO:0005975 GO:hydrolase
      7 SE_W_C GO:0003676,GO:0015074 GO:nucleic


     24 Wild_W_C GO:0004672,GO:0005524,GO:0006468 GO:protein
     17 Wild_W_C GO:0005515 GO:protein
     10 Wild_W_C GO:0006952,GO:0043531 GO:defense
      9 Wild_W_C GO:0030247 GO:polysaccharide
      8 Wild_W_C GO:0008194 GO:UDP-glycosyltransferase
      8 Wild_W_C GO:0004672,GO:0005524,GO:0006468,GO:0030247 GO:protein
      7 Wild_W_C GO:0004497,GO:0005506,GO:0016705,GO:0020037 GO:monooxygenase
      5 Wild_W_C GO:0015144,GO:0016020,GO:0022857,GO:0055085 GO:carbohydrate
      5 Wild_W_C GO:0006486,GO:0016757 GO:protein
      5 Wild_W_C GO:0004672,GO:0005524,GO:0006468,GO:0030246 GO:protein


 3262 Wild_Lat_L_HGO:protein
   1006 Wild_Lat_L_HGO:nucleic
    810 Wild_Lat_L_HGO:defense
    724 Wild_Lat_L_HGO:monooxygenase
    415 Wild_Lat_L_HGO:DNA
    276 Wild_Lat_L_HGO:ATP
    235 Wild_Lat_L_HGO:UDP-glycosyltransferase
    224 Wild_Lat_L_HGO:nucleotide
    217 Wild_Lat_L_HGO:RNA
    204 Wild_Lat_L_HGO:hydrolase
    178 Wild_Lat_L_HGO:recognition
    167 Wild_Lat_L_HGO:oxidoreductase
    156 Wild_Lat_L_HGO:catalytic
    154 Wild_Lat_L_HGO:methyltransferase
    153 Wild_Lat_L_HGO:membrane
    150 Wild_Lat_L_HGO:serine-type
    149 Wild_Lat_L_HGO:peroxidase
    142 Wild_Lat_L_HGO:carbohydrate
    131 Wild_Lat_L_HGO:lipid
    120 Wild_Lat_L_HGO:polysaccharide
    116 Wild_Lat_L_HGO:tRNA
    114 Wild_Lat_L_HGO:ADP
    107 Wild_Lat_L_HGO:DNA-binding
    106 Wild_Lat_L_HGO:signal
    106 Wild_Lat_L_HGO:monoatomic
    100 Wild_Lat_L_HGO:calcium
     87 Wild_Lat_L_HGO:polygalacturonase
     80 Wild_Lat_L_HGO:zinc
     78 Wild_Lat_L_HGO:regulation
     77 Wild_Lat_L_HGO:structural
     72 Wild_Lat_L_HGO:oligopeptide
     65 Wild_Lat_L_HGO:manganese
     63 Wild_Lat_L_HGO:enzyme
     61 Wild_Lat_L_HGO:mRNA
     59 Wild_Lat_L_HGO:ubiquitin-protein
     57 Wild_Lat_L_HGO:GTPase
     51 Wild_Lat_L_HGO:ATP:ADP
     48 Wild_Lat_L_HGO:translation
     47 Wild_Lat_L_HGO:proteolysis
     47 Wild_Lat_L_HGO:microtubule
     44 Wild_Lat_L_HGO:small-subunit
     42 Wild_Lat_L_HGO:metal
     40 Wild_Lat_L_HGO:ligand-gated
     40 Wild_Lat_L_HGO:intracellular
     38 Wild_Lat_L_HGO:biosynthetic
     37 Wild_Lat_L_HGO:S-adenosylmethionine-dependent
     35 Wild_Lat_L_HGO:transcription
     35 Wild_Lat_L_HGO:meristem
     34 Wild_Lat_L_HGO:phosphatase
     33 Wild_Lat_L_HGO:electron
     33 Wild_Lat_L_HGO:acyl-CoA
     32 Wild_Lat_L_HGO:mitochondrial'
	 
	 
	 
