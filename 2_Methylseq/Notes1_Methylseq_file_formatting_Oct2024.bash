#########################################################
# Dryas DNA methylation calling
# nf-core Methylseq
# Mapping WGBS data to reference
# Calling methylated sites
#############################################################

tmux new-session -s Dryas
tmux attach-session -t Dryas

cd /home/celphin/scratch/Dryas/methylseq

#salloc -c1 --time 23:00:00 --mem 120000m --account rrg-rieseber-ac

# copy from Globus
# unzip them first
tar -zxvf Seedlings_MatureFlwr_WGBS_raw_data.tar.gz
tar -zxvf Total_Wild_WGBS_raw_data.tar.gz

# done
#------------------------
# rename raw files Wild parents

# test
for file in *.fastq.gz; do
    # Extract everything after the first underscore and prepend it with the first letter
    new_name="${file#?}_"
    new_name="${new_name#*_}"
    echo "mv "$file" "${file:0:1}_$new_name""
done

# run 
for file in *.fastq.gz; do
    # Extract everything after the first underscore and prepend it with the first letter
    new_name="${file#?}_"
    new_name="${new_name#*_}"
    mv "$file" "${file:0:1}_$new_name"
done
rename .gz_ .gz *

#-----------------------
# rename seedling and Pheno files
#sh rename_files.sh

cd /home/celphin/scratch/Dryas/methylseq/input/Seedlings_MatureFlwr_raw_data
# first copy all files to top level
mv ./*/150bp/*_passed.fastq.gz /home/celphin/scratch/Dryas/methylseq/input/Seedlings_MatureFlwr_raw_data/

#test
while IFS=, read -r oldname newname
do
    echo "mv $oldname $newname"
done < newnames.csv

# rename
while IFS=, read -r oldname newname
do
    mv $oldname $newname
done < newnames.csv

# move all files to one shared directory
# Dryas_total_data
mv *.fastq.gz ../Dryas_total_data

#------------------------
# formatting input file for methylseq run
cd /home/celphin/scratch/Dryas/methylseq/input/Dryas_total_data

#sh make_input_file.sh
# get list of fastq.gz files
ls *fastq.gz >  /home/celphin/scratch/Dryas/methylseq/input/samples_list.csv

cd /home/celphin/scratch/Dryas/methylseq/input

#file to be inputted to nextflowconfig
touch input_files.csv
#header
echo "sample,fastq_1,fastq_2" >> input_files.csv

#adding columns with samplename,fastq1,fastq2
while read -r line
do
    ext="_R1.fastq.gz"
    sample=${line%$ext}
    echo -n $sample",input/Dryas_total_data/"$line"," >> input_files.csv
    read -r line
    echo "input/Dryas_total_data/"$line >> input_files.csv
done <samples_list.csv

#-----------------------------
