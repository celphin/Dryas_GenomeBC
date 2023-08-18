#!/bin/bash
#SBATCH --account=rpp-rieseber
#SBATCH --job-name=dmr_interproscan
#SBATCH --time=1-23:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=120000m

module load StdEnv/2020
module load interproscan/5.63-95.0

input_dir=interproscan_input
output_dir=interproscan_output

mkdir $output_dir

for f in $input_dir/*; do
    no_path="${f##*/}"
    file_name="${no_path%.fasta}"
    srun interproscan.sh -i $f -f tsv -cpu 30 -o $output_dir/interproscan_$file_name.tsv  --goterms;
done

