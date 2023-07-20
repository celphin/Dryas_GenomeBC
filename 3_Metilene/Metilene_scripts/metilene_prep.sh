#metilene_prep.sh: script to prepare files for metilene
################################################
#Input needed:
    # Unzipped methylseq output files (in methylseq_output_dir)
    # category names
    # input names
#Possible adjustments to script: 
    # renaming files process (marked) will vary with input name format
    # eg, looking for SE_*_W_*... and SE_*_C_* has this specific input
#Output:
    # Union bed file ready for finding DMRS for groups h1, h2
    # metilene_"$h1"_"$h2".input 
    # will be located in input_dir

###############################################
# Making input directory
h1="W"
h2="C"
methylseq_output_dir="CE_Seedling_metilene_input_bedGraphs"

input_dir=SE_${h1}_${h2}_input_files 
in_metilene="SE_metilene_"$h1"_"$h2".input"
mkdir ${input_dir}

#-------------------------
# Copying files to input directory

cp -r ${methylseq_output_dir}/*.deduplicated.bedGraph ./${input_dir}

# Renaming files to include category names, update depending on format
cd ./${input_dir}
for bg in SE_*_${h1}*; do mv "$bg" "${h1}_${bg}"; done
for bg in SE_*_${h2}*; do mv "$bg" "${h2}_${bg}"; done

#-------------------------
#Make sorted bedgraph files:

find . -name "*.bedGraph" -print > file_list.txt
while read name
do
bedtools sort -i ${name} >  ${name}_sorted.bedGraph
sort -c -k1,1 -k2,2n ${name}_sorted.bedGraph
cat ${name}_sorted.bedGraph | tr ' ' '\t' > ${name}_sorted_tab.bedGraph
done < file_list.txt

#--------------------------
#Make list of files:

find . -name "${h1}*_sorted_tab.bedGraph" -print > ${h1}_list.txt
find . -name "${h2}*_sorted_tab.bedGraph" -print > ${h2}_list.txt

#---------------------------
#Reformatting files:
# replace ./ in file names
sed -i 's/\.\///g' ${h1}_list.txt
sed -i 's/\.\///g' ${h2}_list.txt


# make csv and space delimited files
sed -z 's/\n/ /g;s/,$/\n/'  ${h1}_list.txt > ${h1}_list_sp.txt
sed -z 's/\n/ /g;s/,$/\n/'  ${h2}_list.txt > ${h2}_list_sp.txt

sed -z 's/\n/,/g;s/,$/\n/'  ${h1}_list.txt > ${h1}_list_csv.csv
sed -z 's/\n/,/g;s/,$/\n/'  ${h2}_list.txt > ${h2}_list_csv.csv
#-----------------------------
#make union bedfile for metilene input
NA="NA"

sed ':a;N;$!ba;s/\n/ /g' ${h1}_list.txt > h1files
sed ':a;N;$!ba;s/\n/ /g' ${h2}_list.txt > h2files


h1filesinput=$(more h1files)
h2filesinput=$(more h2files)


sed ':a;N;$!ba;s/\n/ /g' ${h1}_list.txt | \
sed 's/_R1_val_1_bismark_bt2_pe.deduplicated.bedGraph_sorted_tab.bedGraph//g' > H1name

sed ':a;N;$!ba;s/\n/ /g' ${h2}_list.txt | \
sed -z 's/_R1_val_1_bismark_bt2_pe.deduplicated.bedGraph_sorted_tab.bedGraph//g' > H2name

h1name=$(more H1name)
h2name=$(more H2name)

bedtools unionbedg -header -names $h1name $h2name -filler $NA -i $h1filesinput $h2filesinput \
| cut -f1,3- | sed 's/end/pos/' > "$in_metilene"

########################################################

