##########################################
#Processing Blast output files - to have clean gene names
#Goal - remove any duplicates, leave only fragments with gene names
#TODO: check sort order for loop
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
#Removes OU, OX (genome assembly), and AP (pseudomolecules)
grep -v "OU\|AP\|OX\|OW" intersect_Wild_W_C_Mat_Sen.out > grep.out
#try removing duplicates:

#Remove duplicates of the same exact blast (field #4):
awk '!a[$1$4$5$6$7]++ {print}' grep.out > no_duplicates2.out

#Check sort order
#Isolate top 3 for each
touch top3.out
key=""
count=0
num_copies=3
while read line ; 
do 
    name=${line%%g*}
    if [$name==$key]
    then 
        if [$count<3]
        then 
            echo $line >> top3.out
            count=$((count+1))
        fi
    else
        key=$name
        count=1
        echo $line >> top3.out
    fi
done < intersect_Wild_W_C_Mat_Sen.out

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

