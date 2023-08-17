#######################################
#Script to create methylseq input file from fastq.gz files
#Requires: be in directory which contains all fastq.gz files containing sequence for methylseq
#Output: input_files.csv
#Modifications:
    #cd into desired directory containing files
    #file_dir location
    #sample_list path
#######################################
#list of fastq.gz files 
file_dir="~/scratch/testfilecopy/dir_list.txt"
sample_list=" ~/scratch/testfilecopy/samples_list.csv"

ls > $file_dir
grep "fastq.gz" ${file_dir} > ${sample_list}

#file to be inputted to nextflowconfig
touch input_files.csv
#header
echo "sample,fastq1,fastq2" >> input_files.csv

#adding columns with samplename,fastq1,fastq2
while read -r line
do 
    ext="_R1.fastq.gz"
    sample=${line%$ext}
    echo -n $sample",input/"$line"," >> input_files.csv
    read -r line
    echo "input/"$line >> input_files.csv
done <samples_list.csv

# removing all intermediate  files
rm dir_list.txt
rm samples_list.csv
