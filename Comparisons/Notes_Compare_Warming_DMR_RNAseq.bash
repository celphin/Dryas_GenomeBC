####################################################
# Compare genes/regions that show warming signal in RNAseq and DMRs
# March 2023
###############################################

# Compare regions between RNAseq and DMRs

# Step 1: add 1000bp to each end of DMR
# add 500bp on each side of the DMR
awk '{$3+=1000}1' Total_DMRs_100-10-10_July2022_qval.1e-20.bed > Total_DMRs_100-10-10_July2022_qval.1e-20_extended.bed
awk '{$2-=1000}1' Total_DMRs_100-10-10_July2022_qval.1e-20_extended.bed > Total_DMRs_100-10-10_July2022_qval.1e-20_extended2.bed

# Step 2: get gene positions of each differntially expressed gene
Do1_00107G00001V1.1 CDS=1-1578


# Step 3: Union of DMRs and RNAseq regions
# any overlap? - maybe add more/less to DMRs

bedtools unionbedg -header -names ${Site1}_${h1}name ${Site1}_${h2}name -filler $NA -i ${Site1}_${h1}files ${Site1}_${h2}files | cut -f1,3- | sed 's/end/pos/' > "$Site1_$in_metilene"
bedtools unionbedg -header -names ${Site2}_${h1}name ${Site2}_${h2}name -filler $NA -i ${Site2}_${h1}files ${Site2}_${h2}files | cut -f1,3- | sed 's/end/pos/' > "$Site2_$in_metilene"
bedtools unionbedg -header -names ${Site3}_${h1}name ${Site3}_${h2}name -filler $NA -i ${Site3}_${h1}files ${Site3}_${h2}files | cut -f1,3- | sed 's/end/pos/' > "$Site3_$in_metilene"
bedtools unionbedg -header -names ${Site4}_${h1}name ${Site4}_${h2}name -filler $NA -i ${Site4}_${h1}files ${Site4}_${h2}files | cut -f1,3- | sed 's/end/pos/' > "$Site4_$in_metilene"



# Step4: BLAST genes that match 


# Step 5: determine how to plot




####################
# Alternatively

# Get fasta of DMRs
sed -e 's/ /\t/g' Total_DMRs_100-10-10_July2022_qval.1e-20_extended2.bed > Total_DMRs_100-10-10_July2022_qval.1e-20_extended2_tab.bed
bedtools getfasta -fi /home/celphin/scratch/Dryas/Dryas_octopetala_reference/genomes/Dryas_octopetala_H1.supercontigs.fa -bed Total_DMRs_100-10-10_July2022_qval.1e-20_extended2_tab.bed > Total_DMRs_100-10-10_July2022_qval.1e-20_extended.fasta
more Total_DMRs_100-10-10_July2022_qval.1e-20_extended.fasta
blastn -db nt -query Total_DMRs_100-10-10_July2022_qval.1e-20_extended.fasta -out blast_Total_DMRs_100-10-10_July2022_qval.1e-20_extended.out -remote -outfmt "6 qseqid stitle" 
-evalue 0.05
-word_size 11
-gapopen 5
-gapextend 2
-penalty -3
-reward 2
-max_target_seqs 6

blastn -db nt -query Total_DMRs_100-10-10_July2022_qval.1e-20_extended.fasta -out blast_Total_DMRs_100-10-10_July2022_qval.1e-20_extended_settings_test.out -remote -outfmt "6 qseqid stitle" -evalue 0.05 -word_size 11 -gapopen 5 -gapextend 2 -penalty -3 -reward 2 -max_target_seqs 6

Error: [blastn] internal_error: (Severe Error) Blast search error: Details: search failed. # Informational Message: [blastsrv4.REAL]: Error: CPU usage limit was exceeded, resulting in SIGXCPU (24).

more blast_Total_DMRs_100-10-10_July2022_qval.1e-20_extended_settings_test.out

# take blast results and compare to RNAseq results

#################################3
# convert bedGraph files to bed files to fasta files

cd /home/celphin/projects/rpp-rieseber/celphin/Dryas/DMRs/June2022_Metilene_DMRs/
cat Total_DMRs_100-10-10_July2022_qval.1e-20.bedgraph | awk '{print $1 "\t" $2 "\t" $3}' > Total_DMRs_100-10-10_July2022_qval.1e-20.bed

#----------------------------
# add 500bp on each side of the DMR
awk '{$3+=1000}1' Total_DMRs_100-10-10_July2022_qval.1e-20.bed > Total_DMRs_100-10-10_July2022_qval.1e-20_extended.bed
awk '{$2-=1000}1' Total_DMRs_100-10-10_July2022_qval.1e-20_extended.bed > Total_DMRs_100-10-10_July2022_qval.1e-20_extended2.bed

Do1_01_a00004 4252188 4253712
Do1_03_a00002 5352537 5354712
Do1_03_a00002 12037475 12038826
Do1_03_a00002 12806450 12807803
Do1_03_a00004 1828629 1829973
Do1_05_a00001 17163975 17165496
Do1_05_a00003 1152346 1153487
Do1_05_a00003 1479741 1481401
Do1_05_a00004 323760 325304
Do1_06_a00002 2717520 2718959
Do1_07_a00002 11004673 11006313
Do1_07_a00002 12442051 12443627
Do1_07_a00004 2484455 2485900
Do1_a00028 195151 196295

Do1_01_a00004   4251688 4254212
Do1_03_a00002   5352037 5355212
Do1_03_a00002   12036975        12039326
Do1_03_a00002   12805950        12808303
Do1_03_a00004   1828129 1830473
Do1_05_a00001   17163475        17165996
Do1_05_a00003   1151846 1153987
Do1_05_a00003   1479241 1481901
Do1_05_a00004   323260  325804
Do1_06_a00002   2717020 2719459
Do1_07_a00002   11004173        11006813
Do1_07_a00004   2483955 2486400
Do1_a00028      194651  196795



#----------------------
# https://bedtools.readthedocs.io/en/latest/content/tools/getfasta.html
# bedtools getfasta [OPTIONS] -fi <input FASTA> -bed <BED/GFF/VCF>

sed -e 's/ /\t/g' Total_DMRs_100-10-10_July2022_qval.1e-20_extended2.bed > Total_DMRs_100-10-10_July2022_qval.1e-20_extended2_tab.bed

bedtools getfasta -fi /home/celphin/scratch/Dryas/Dryas_octopetala_reference/genomes/Dryas_octopetala_H1.supercontigs.fa -bed Total_DMRs_100-10-10_July2022_qval.1e-20_extended2_tab.bed > Total_DMRs_100-10-10_July2022_qval.1e-20_extended.fasta

more Total_DMRs_100-10-10_July2022_qval.1e-20_extended.fasta

# try the only blast
>Do1_a00028:194651-196795
GTATATGATAACATCAAAATCGAAGAGGTGAAGTGTGTACCAGAGTTGAAGTTAAAGGCGAGAGGAATGCTGATCCGAGAGACCTCCGGAGCAAACCAATTCCCTGTAATAAGAATTGGACGAGTCTTCCACTCTTTACAAAATGTAGGCTTGTTC
CCAAAGAACACACAATGGTGTCGAGGCGTAATATGGTAATACCCAATCTCATTCCTGACAGTAGCACAACAATGATAAGCAAAGAAATCGATAAGTCCCAAGGAGACACCCCAACACAGGCCGAGCTCGAGAAAAGCAGAAATTGTAAGAAAACTA
TTGGGGTTTACTTGCATTGGAGACAAGCCCAAGAAGCTGAAAACCAGATGGAATACGAGGTGTAAAGGAAACCTAACATGATCAAGATGATTAAAAGACAGAACTATAGAATTCCCAGGTAATCCTGCTCCAGTTCGAATGTTAGCATCACGAGGT
AAGCCACGGAAAACTACATCATCTGGGGCTTTGCACTTTTGTCGGGCTATGGCAAGGTCTTGGTCAGTTTGGTTCCACAAAGTTGTGTTATTTTTATAATTCTTTGCCATAAAAATAAAATGCAGAAAGGCACTAACCTAAACTAAGAGAAGATAA
GAGCTAATTTGAGAAGTGCCCAAGTTCTGCTAAACATCGCCGAACCCAAGTCTACTGAACCAAGGAAACGCCGCAGACCTCGCCGAAACCCAAGTCTGCTGAACCAGGGAAATGCCGCAAACCTTGCCGAACCCAAGTCTACTGAACCAAGGAAAC
GTCGCAGACCCCATCGAAAGCCACCAGCTTCGCTGAACAGCTGAACGCACGCAGCAAGTCACACCGAACCGCCGAGCGCCCAGGCCTTACTGAACTGCCGAACCCAGCAACCTTGCTGAATCGTCGAACTGCCGAAACCAGAGTTGCTAGATCGCC
AAGACTCGCCGAACAAGAAGACTCGAATTCTCGCAGCCGAAACGGCAAGCAAGTTGCATGTTGTCGTCACCGGACCTCAACTCTGGTGAGACGCGAAGAGGGTCGCGTTACTGAAATTCAGGAAACCCAAGCCTTGACAGCGTGTTTGGAGACCAA
AGCCGACTTCGATGCCGAGAGCTACTAGCAAACTGCGAGACAAACAGATGCGATGATCGTGCGCAAGTCTGAGTGTGGTGCACAGTCGGTGGCGGCGAACGCGGCTGTGCAAGAATTCAATGGTGAGGCAGCCCAGATGAGTAATAGGGAGATATG
GCACCGAGAAGGATAAGTGGTGAAATTTTATTCAAAACATGAAAGGAGTTCGCGCAGGTGATGATGTCTTCTCCTTAAAATTGCTTTTTAATGGTTTAAAGATGTTTATTTAATTTGAAAATTGGCTAATATTTATCCGGCCCAACTAGGTATGTG
ATAAAGTCTATCATATTAATTTTATTTGAATTAGTGAAAAAATAGGAAATCCAGATGATGAAGGAAAATAATTGATTTTTCTTTAATAAAATCTATAGAGATTAAATTAAAGAAAAAGGGAGCATTATAGGGGAGATTATTTTCGTGAGCCCAAAT
TTGAAGTTTGAAGCCTAGAATATCAAACAACATTTGGGTTGAACCCATAAGAATAATTTAACTGAGTTTGGGTTTGTGAAGCTATTATGTTGAAAAGTCCAGGATGTAGTTAAGTTTATCCTCGTGACTAAATTCAATTGAAAACGTGGGCGTGAA
GTGCAAGTAACACCACGAAACATCGGAATGGAAGACGCATGGATGGAGAAACTTGCGCCTTACGTCGTAAACTTTATCCACAAGAAATTTTTGACTATATATATATATATATACATAATCATAGTCATTAGAGAAGGTACACAATCCAACCATTTG
AGCATTTATTTGCATCAGACATAAAAATAACTTGAGCGTCGGAGTATCTTTTGTGGGTACCTCACCCGCCGTACAAAAACAAACTCAAATCAAAGTCATATATATAGCGAGCGTGATCTAGAGTTCAATTCCAAAGTCAACATCAACAACCTGCAG
GCCAAGTTATCTTCTTTGAAGATTTTGGGTCCAATAGGCAGCCAGCAAATGGAGGTTCAGGGGCGCAGCAGAGCCAGTACTGAGCAGTTGGGAGAGGTTGTATCATGGTGATTTTC


>Do1_05_a00003:1151846-1153987
GTCTTTCCCACAAACGGAGCCAATTGTTGAACCCAAAATTTCTCGAAGATAAATTAAGTTAACAATAAAATTATTTGGAGCTTGAACGTTTGACGTGAAGCGGATTGAACGGCTTGGTTTCTCTGTCTGACTTGTCAAAGCTTCTTCTCCTCAAAG
ATGGGTGAGGTTACCTGCAAAAAACCCTCTGATGCCTAAGTTAGCACAAGAAAACTCTCAAGTGTCTGTAAAATGATTTTCTACGTACCTTTCATTAATGACCTCCATGAGACTTTATAATAGTTTGATTCAACTAATTTCATGGGTTTTCTCCAC
CTTTCTTGTCGTCTAAGATGGTGGACTTTGTTGGAGAAAATCATGGGAATTTCTCTAGGAGATTTTCAGTGGAGACTTGGTCTTCTATTGAAGAGATTTTTCCCACGTTGTCTCTGCATTTACTTCCTTATTTAATTCCTTCATTATCCTCTTAAT
TCTCAAAAATCAAAGGTTAAGAAATAAATGTTAAATAATTATCATTAATTAATAACAATGAGTCTCAAATTGAATGGGAAACAAAGAGATTTTGACTGGGCTTGTATTCCACGTGAAAATTCTTCCCAATTGCCATAACTCCTATTTATCTCATTG
TTAATAAATAAAATATTATTGAAATGTAATTATCATGAAATTGACTCCTTATCAAGATAGGTTATCAATAAATATTTATAAAATTTAAATGACACATCATCTTTGATTGAGTCATTGTTCTTCAAATGACAAATGGGAAATTGCCATGTGTCGCGG
TAATTTTGCTCCAACATAGAATATTCTAGCTTTAATTTGTATGAGTTTTGAGGAAGTCTCTAGAATTTTGTAGTATTTTCCAATAATGTGTATTAATTTTCAACAAATGTGAACCCACCAAAACAACACTGCCTAAGCTTCATGGTATTCCTGTTT
TACCCAAGGACAGTATAGCAGCCAAGGACAAGCTGAACACCACTGCTGGGTCTGGGAGAAAAGCCGCCAGTATGACTCTGGAAACCGCGTCTGTGTCGATGTAGACCATATGGACGGTCGACGTCGAGTGTCCTATATAGTCGACAAAGCCAATAT
ATTATAGTTGACGTCGACCATATGGACAGTCAAAAATATATTATGGTCGACGTTTATACATCAGTTTGGCTTTGTAAACCAAATTTAGATTTCAAATCTGTAAGGGACTTGTTTACTTTAAAGGATCAATAAGGACTTGTCTGGAAATTTCTCTAT
TTCCAGATATACCAGAAGGTCTGCAAATTTCTCAGTTTCTAGATATACCAGTAGACAAAGATGAACAATTATGTAGATTAGCAGCAAGCCAACTGTCCTAGTCTTAGATGCAGAAGGTTCCCAGGTAGACCAATAGCTTTCCACATTCTCCTGGCT
CTACCTTGGATTTTTTCCAAATTGCCTGCATAAATAGCTCTATAAACTCTGGAATGAACCACATGCATTCCCAAATACTTACCCTTAGCTTTATAAACTTTGTTTATAAACTTGCTACATATGTTTAGAGAAGCTTTTGTTGGGACGCATCAGGGT
AGAGTTAATTTCTTCCATGTTGCGTTTGGGTTTTTCTTTTGTTGGAGAGAGAGAGAGAGAAAGCAGAGATGAGTGAGAGGTTATGTGGGAAGATGAGAGAGAGGTTATGTGGGAAAATGAGAAAGCAGGTATATTCATTGATTTCAAAATAAGATG
TCATGACAAAAGGTAATAACGAAAAAATGTAAGGGGTAATTGAGAGTTTTCCTGTCTACAAAGTGTCTTTTTCTGCAAAGATTAGTGATGGTTCTAAATTTGAAAGGGGGGGAGACTATTCAGACCACCCAAAACCACTTAGTTGACCACCCCATT
AATTTTTTGCCCATTTTTTTTTTCCTATTTTGCCCCTGCATTAAAAAAAAACTCAGTTTTCCAATTTTGGCCGGCCACCCCTCCCCGCTACCGGCGCAACCCCACCCCCTCCCCACTGCCGGCGCGACCCCACCCCTCCTCTTCGCCGGTCGACCC
CACGCCTCCCCATCGCCGGCGCGACCCACCCACCCCTCCCCTTCTCTCTCTCTCTCTCTTCCCTCTCTCTGCCATCTTCTCTCTCTCCCCTCTCTCTACCTCCCTCCCTCACC
PREDICTED: Rosa chinensis probable amidase At4g34880 (LOC112173655), mRNA 

###############################################
# Blast bedGraph files to ncbi database  - see Ancient Cassiope BLAST

# BLAST

# run blast
cd /home/celphin/projects/rpp-rieseber/celphin/Dryas/DMRs/June2022_Metilene_DMRs/

module load StdEnv/2020
module load gcc/9.3.0
module load blast+/2.12.0

# https://www.metagenomics.wiki/tools/blast/blastn-output-format-6

blastn -db nt -query Total_DMRs_100-10-10_July2022_qval.1e-20_extended.fasta -out blast_Total_DMRs_100-10-10_July2022_qval.1e-20_extended.out -remote -outfmt "6 qseqid stitle"

more blast_Total_DMRs_100-10-10_July2022_qval.1e-20_extended.out

Do1_05_a00004:323760-325304 - syntaxin-71 
Do1_06_a00002:2717520-2718959 - POLYGALACTURONASE INVOLVED IN EXPANSION3 - AT1G48100.1
Do1_07_a00002:12442051-12443627 - NADH dehydrogenase subunit 4 gene, complete cds; mitochondrial


#------------------
# try blast online - works to find same genes as before


# try more strict BLAST or less strict BLAST online
blastn -help

USAGE
  blastn [-h] [-help] [-import_search_strategy filename]
    [-export_search_strategy filename] [-task task_name] [-db database_name]
    [-dbsize num_letters] [-gilist filename] [-seqidlist filename]
    [-negative_gilist filename] [-negative_seqidlist filename]
    [-taxids taxids] [-negative_taxids taxids] [-taxidlist filename]
    [-negative_taxidlist filename] [-entrez_query entrez_query]
    [-db_soft_mask filtering_algorithm] [-db_hard_mask filtering_algorithm]
    [-subject subject_input_file] [-subject_loc range] [-query input_file]
    [-out output_file] [-evalue evalue] [-word_size int_value]
    [-gapopen open_penalty] [-gapextend extend_penalty]
    [-perc_identity float_value] [-qcov_hsp_perc float_value]
    [-max_hsps int_value] [-xdrop_ungap float_value] [-xdrop_gap float_value]
    [-xdrop_gap_final float_value] [-searchsp int_value]
    [-sum_stats bool_value] [-penalty penalty] [-reward reward] [-no_greedy]
    [-min_raw_gapped_score int_value] [-template_type type]
    [-template_length int_value] [-dust DUST_options]
    [-filtering_db filtering_database]
    [-window_masker_taxid window_masker_taxid]
    [-window_masker_db window_masker_db] [-soft_masking soft_masking]
    [-ungapped] [-culling_limit int_value] [-best_hit_overhang float_value]
    [-best_hit_score_edge float_value] [-subject_besthit]
    [-window_size int_value] [-off_diagonal_range int_value]
    [-use_index boolean] [-index_name string] [-lcase_masking]
    [-query_loc range] [-strand strand] [-parse_deflines] [-outfmt format]
    [-show_gis] [-num_descriptions int_value] [-num_alignments int_value]
    [-line_length line_length] [-html] [-sorthits sort_hits]
    [-sorthsps sort_hsps] [-max_target_seqs num_sequences]
    [-num_threads int_value] [-mt_mode int_value] [-remote] [-version]

blastn -db nt -query Total_DMRs_100-10-10_July2022_qval.1e-20_extended.fasta -out blast_Total_DMRs_100-10-10_July2022_qval.1e-20_extended.out -remote -outfmt "6 qseqid stitle" 
-evalue 0.05
-word_size 11
-gapopen 5
-gapextend 2
-penalty -3
-reward 2
-max_target_seqs 6

blastn -db nt -query Total_DMRs_100-10-10_July2022_qval.1e-20_extended.fasta -out blast_Total_DMRs_100-10-10_July2022_qval.1e-20_extended_settings_test.out -remote -outfmt "6 qseqid stitle" -evalue 0.05 -word_size 11 -gapopen 5 -gapextend 2 -penalty -3 -reward 2 -max_target_seqs 6

Error: [blastn] internal_error: (Severe Error) Blast search error: Details: search failed. # Informational Message: [blastsrv4.REAL]: Error: CPU usage limit was exceeded, resulting in SIGXCPU (24).

more blast_Total_DMRs_100-10-10_July2022_qval.1e-20_extended_settings_test.out


# finish Alex samples
- running now

# check for same Amidase gene
yes is there

######################################
# run sample/site overlaps and BLAST


################################
# Compare blast results with RNAseq data - any overlap?

# find matches between studies - e.g. the overlap/union between the warming and phenology or seedling and parents DMRs

# https://stackoverflow.com/questions/13272717/inner-join-on-two-text-files
#https://itectec.com/ubuntu/ubuntu-awk-compare-2-files-and-print-columns-from-both-files/ 

# this returns all rows from second file that match the first column of first file
awk 'NR==FNR {a[$1]; next} $1 in a {print $0, a[$2]}' OFS='\t' LAT_DMRs_qval.0.05.out metilene_qval.0.05.out > LAT_Total.txt
# returns all that match for the 'chrom'

bedtools unionbedg -header -names ${Site1}_${h1}name ${Site1}_${h2}name -filler $NA -i ${Site1}_${h1}files ${Site1}_${h2}files | cut -f1,3- | sed 's/end/pos/' > "$Site1_$in_metilene"
bedtools unionbedg -header -names ${Site2}_${h1}name ${Site2}_${h2}name -filler $NA -i ${Site2}_${h1}files ${Site2}_${h2}files | cut -f1,3- | sed 's/end/pos/' > "$Site2_$in_metilene"
bedtools unionbedg -header -names ${Site3}_${h1}name ${Site3}_${h2}name -filler $NA -i ${Site3}_${h1}files ${Site3}_${h2}files | cut -f1,3- | sed 's/end/pos/' > "$Site3_$in_metilene"
bedtools unionbedg -header -names ${Site4}_${h1}name ${Site4}_${h2}name -filler $NA -i ${Site4}_${h1}files ${Site4}_${h2}files | cut -f1,3- | sed 's/end/pos/' > "$Site4_$in_metilene"

#------------------------------
#header get rid of enter in first line

#header should be:
#chrom pos ${Wname} ${Cname}
# but the Cnames are one new line
W_ALAS0W_8_265
	C_ALAS_00C_227

#edit to remove extra enter
nano "$Site1_$in_metilene"
nano "$Site2_$in_metilene"
nano "$Site3_$in_metilene"
nano "$Site4_$in_metilene"

#-----------------------------
# check header
head -5 "$Site1_$in_metilene" > subset1.txt
head -5 "$Site2_$in_metilene" > subset2.txt
head -5 "$Site3_$in_metilene" > subset3.txt
head -5 "$Site4_$in_metilene" > subset4.txt

grep '^c.*' subset1.txt | wc -w
grep '^c.*' subset2.txt | wc -w
grep '^c.*' subset3.txt | wc -w
grep '^c.*' subset4.txt | wc -w
Alex -52
SV -14
LAT -22
ALAS -22

# grep '^Do.*' "$Site1_$in_metilene" | wc -l # too long
Alex - 8506281
SV - 8088914
LAT - 8318294
ALAS - 8315924

nano "$Site1_$in_metilene"
nano "$Site2_$in_metilene"
nano "$Site3_$in_metilene"
nano "$Site4_$in_metilene"
# delete enter on second line if still there



########################################################################
########################################################################


