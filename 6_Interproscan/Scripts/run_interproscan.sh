################################################
#Run_interproscan.sh
#With required loaded modules, and input directory
#runs interproscan on a series of blasts
#Input directory: output of prep_interpsocan.sh, or any directory with fasta files
#Output: list of Interproscan tsv files
################################################
input_dir=interproscan_input
output_dir=interproscan_output

mkdir $output_dir


for f in $input_dir/*; do
    no_path="${f##*/}"
    file_name="${no_path%.fasta}"
    interproscan.sh -i $f -f tsv -o $output_dir/interproscan_$file_name.tsv  --goterms;
done
