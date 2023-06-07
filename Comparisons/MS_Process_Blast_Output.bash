##########################################
#Processing Blast output files - to have clean gene names
#Goal - remove any duplicates, leave only fragments with gene names
#Note: This file is mainly scratch work for now, ignore most things
##########################################
#Blast name abbreviation meanings:
    #XR, XM - predicted gene name, contains uncharacterized locus
    #OU, OX, OW - genome assembly (both with and without chromosome #)
    #AP - pseudo molecule
    #NC, OM, ON, MZ, MW - mitochondrial complete genome
    #OL - chromosome # and mitochondria
##########################################

cd scratch
mkdir process_blast_output
cd process_blast_output

tmux new-session -s Process_Blast
tmux attach-session -t Process_Blast
##################################################
#Testing filtering process:
#Want:
    #1) remove all OU/PX/AP/OW matches
    #2) Remove any identical copies for gene name
    #3) Isolate top 3 for each blast


cp ~/projects/def-rieseber/Dryas_shared_data/MS_blast_output/blast_ncbi_rosaceae_intersect_Wild_W_C_Mat_Sen.out intersect_Wild_W_C_Mat_Sen.out

#try removing duplicates:
awk '!a[$1]++' intersect_Wild_W_C_Mat_Sen.out > no_duplicates.out
# Problem : not all of them have the first one be the protien name 

#try isolate XM (potentially not getting name for all)
grep "XM" intersect_Wild_W_C_Mat_Sen.out > grep_xm.out
awk '!a[$1]++' grep_xm.out > no_duplicate_grep.out
grep -v "uncharacterized" no_duplicate_grep.out > genes_only.out

### Trying to pipe: 
grep "XM" intersect_Wild_W_C_Mat_Sen.out | awk '!a[$1]++'|\
    grep -v "uncharacterized" > genes_only2.out
#works correctly - do once done




#Try getting ALL genes (same wc as no_duplicates, and then cleaning)

#Want:
    #1) remove all OU/PX/AP/OW matches
    #2) Remove any identical copies for gene name
    #3) Isolate top 3 for each blast


#Removes OU, OX (genome assembly), and AP (pseudomolecules)
grep -v "OU\|AP\|OX\|OW" intersect_Wild_W_C_Mat_Sen.out > grep.out

awk -F ' ' '{ key = $1 } key && count[key]++ < 3' RS=" \n" grep.out > top3.out
awk -F ' ' '{ key = $1 } !seen[key]++ && count[key]++ < 3' RS=" \n" input_file > output_file

#First three lines of each blast
awk -F ' ' '!a[$1]++ && c[$1]<3' intersect_Wild_W_C_Mat_Sen.out > top3.out

awk -F ' ' '!a[$1]++' intersect_Wild_W_C_Mat_Sen.out > top3.out

awk -F '  ' -v RS='\n' '{count[$1]++; lines[$1][count[$1]]=$0} END {for (key in count) {for (i=1; i<=3 && i<=count[key]; i++) print lines[key][i]}}' intersect_Wild_W_C_Mat_Sen.out > top3.out

awk -F '  ' -v RS='\n' '{count[$1]++; lines[$1][count[$1]]=$0} END {for (key in count) {for (i=1; i<=count[key]; i++) {if (i <= 3) print lines[key][i]}}}' input_file > output_file

#try removing duplicates:

$ awk '{blast_id[$1]++;} END {for (var in id) print var, "id: ", Ip[var]," times" }' >  numtimes.txt

#Remove duplicates of the same exact blast (field #4):
awk '!a[$1$4$5$6$7]++ {print}' grep.out > no_duplicates2.out

awk '{print $4$5$6$7}' grep.out > onlyfield4.out
#remove duplicates
awk '!a[$4$5$6$7]++ {print}' grep.out > no_duplicates4.out
#########################################################################
#Actual filtering:
#For loop to create output for every file in blast folder
#Input: folder of ncbi blast outputs
#Output: Top 3 genes per blast, no repititions 
#########################################################################
#input_folder = 
#output_folder=Blast_rosaceae_genes





#output_folder=Blast_arabidopsis_genes


cd ~/scratch
mkdir Process_Blast_Output

#For loop to create output for every file in blast folder
#Input: folder of ncbi blast outputs
#Output: Top 3 genes per blast, no repititions 

#input_folder = 
#output_folder=Blast_rosaceae_genes
#output_folder=Blast_arabidopsis_genes

mkdir output_folder


#Body of for loop
grep -v "OU\|AP\|OX\|OW" ${filename} | awk '!a[$1$4$5$6$7]++ {print}' | #isolate top3 > "$genes_{filename}"

#done


#

