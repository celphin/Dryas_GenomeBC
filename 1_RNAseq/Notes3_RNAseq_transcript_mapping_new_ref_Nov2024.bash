########################################
# Dryas differential gene expression data
# Mapping to transcriptome 
# Nov 2024 map to new refernece genome 
#####################################

# Try interactive runs
# https://github.com/agshumate/Liftoff

tmux new-session -s RNA
tmux attach-session -t RNA

salloc -c40 --time 2:55:00 --mem 120000m --account def-rieseber

source ~/miniconda2/bin/activate liftoff
module load StdEnv/2020
module load python/3.10.2
module load scipy-stack/2022a
module load gcc/9.3.0
module load cuda/11.7
module load hdf5/1.12.1
module load parasail/2.5
module load minimap2/2.24

cd /home/celphin/scratch/Dryas/Liftoff/

liftoff \
-g Dryas_octopetala_H1.gff3 -p 40 -o new_Dryocto_liftoffpolish.gff3 \
-dir Dryocto_intermed -u Dryocto_unmapped_features.txt \
-infer_genes -copies -a 0.95 -s 0.95 -d 5.0 -flank 0.8 -polish \
DoctH0_Main.fasta Dryas_octopetala_H1.supercontigs.fa

