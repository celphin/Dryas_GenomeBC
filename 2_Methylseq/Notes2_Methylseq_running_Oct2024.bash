#########################################################
# Dryas DNA methylation calling
# nf-core Methylseq
# Mapping WGBS data to reference
# Calling methylated sites
#############################################################
# organize folders

mv input/ nf-core-methylseq_2.7.0/2_7_0/
mv output/ nf-core-methylseq_2.7.0/2_7_0/
mv nextflow.config nf-core-methylseq_2.7.0/2_7_0/

#--------------
# Try run to setup reference genome  in Bismarkindex/
# on Cedar

tmux new-session -s Dryas
tmux attach-session -t Dryas

cd /home/celphin/scratch/Dryas/methylseq

salloc -c1 --time 23:00:00 --mem 120000m --account rrg-rieseber-ac

cd /home/celphin/scratch/Dryas/methylseq
source nf-core-env/bin/activate
module load nextflow/24.04.4
module load apptainer/1.3.4
export NXF_SINGULARITY_CACHEDIR=/project/rrg-rieseber-ac/NXF_SINGULARITY_CACHEDIR

nextflow run nf-core-methylseq_2.7.0/2_7_0/  -profile singularity,cedar 

# executor >  slurm (240)
# [b0/3103ee] NFC…SMARK_GENOMEPREPARATION (BismarkIndex/DoctH0_Main.fasta) | 0 of 1
# [-        ] NFCORE_METHYLSEQ:METHYLSEQ:CAT_FASTQ                         -
# [52/17c391] NFCORE_METHYLSEQ:METHYLSEQ:FASTQC (W_WILL7W_448_107)         | 0 of 131
# [36/1561e3] NFCORE_METHYLSEQ:METHYLSEQ:TRIMGALORE (W_LATC5W_18_190)      | 0 of 131
# [-        ] NFCORE_METHYLSEQ:METHYLSEQ:BISMARK:BISMARK_ALIGN             -
# [-        ] NFCORE_METHYLSEQ:METHYLSEQ:BISMARK:SAMTOOLS_SORT_ALIGNED     -
# [-        ] NFCORE_METHYLSEQ:METHYLSEQ:BISMARK:BISMARK_DEDUPLICATE       -
# [-        ] NFC…METHYLSEQ:METHYLSEQ:BISMARK:BISMARK_METHYLATIONEXTRACTOR -
# [-        ] NFCORE_METHYLSEQ:METHYLSEQ:BISMARK:BISMARK_COVERAGE2CYTOSINE -
# [-        ] NFCORE_METHYLSEQ:METHYLSEQ:BISMARK:BISMARK_REPORT            -
# [-        ] NFCORE_METHYLSEQ:METHYLSEQ:BISMARK:BISMARK_SUMMARY           -
# [-        ] NFC…E_METHYLSEQ:METHYLSEQ:BISMARK:SAMTOOLS_SORT_DEDUPLICATED -
# [-        ] NFC…_METHYLSEQ:METHYLSEQ:BISMARK:SAMTOOLS_INDEX_DEDUPLICATED -
# [-        ] NFCORE_METHYLSEQ:METHYLSEQ:QUALIMAP_BAMQC                    -
# [-        ] NFCORE_METHYLSEQ:METHYLSEQ:PRESEQ_LCEXTRAP                   -
# [-        ] NFCORE_METHYLSEQ:METHYLSEQ:MULTIQC                           -

########################################
# After run edit config for BismarkIndex

params{

    // Input options
    input                      = '/home/celphin/scratch/Dryas/methylseq/input/input_files.csv'
    fasta                      = '/home/celphin/scratch/Dryas/methylseq/input/DoctH0_Main.fasta'
    //add after second run:
    //bismark_index            = '/home/celphin/scratch/Dryas/methylseq/output/BismarkIndex/'
}

###################################################
# restart pipeline with edited config to show new reference in Bismarkindex/

tmux new-session -s Dryas
tmux attach-session -t Dryas

cd /home/celphin/scratch/Dryas/methylseq

salloc -c1 --time 23:00:00 --mem 120000m --account rrg-rieseber-ac

cd /home/celphin/scratch/Dryas/methylseq
source nf-core-env/bin/activate
module load nextflow/24.04.4
module load apptainer/1.3.4
export NXF_SINGULARITY_CACHEDIR=/project/rrg-rieseber-ac/NXF_SINGULARITY_CACHEDIR

nextflow run nf-core-methylseq_2.7.0/2_7_0/  -profile singularity,cedar -resume

#######################################################
# resume pipeline with more time

tmux new-session -s Dryas
tmux attach-session -t Dryas

cd /home/celphin/scratch/Dryas/methylseq

salloc -c1 --time 23:00:00 --mem 120000m --account rrg-rieseber-ac

cd /home/celphin/scratch/Dryas/methylseq
source nf-core-env/bin/activate
module load nextflow/24.04.4
module load apptainer/1.3.4
export NXF_SINGULARITY_CACHEDIR=/project/rrg-rieseber-ac/NXF_SINGULARITY_CACHEDIR

nextflow run nf-core-methylseq_2.7.0/2_7_0/  -profile singularity,cedar -resume

#######################################################
# resume pipeline with more time

tmux new-session -s Dryas
tmux attach-session -t Dryas

cd /home/celphin/scratch/Dryas/methylseq

salloc -c1 --time 23:00:00 --mem 120000m --account rrg-rieseber-ac

cd /home/celphin/scratch/Dryas/methylseq
source nf-core-env/bin/activate
module load nextflow/24.04.4
module load apptainer/1.3.4
export NXF_SINGULARITY_CACHEDIR=/project/rrg-rieseber-ac/NXF_SINGULARITY_CACHEDIR

nextflow run nf-core-methylseq_2.7.0/2_7_0/  -profile singularity,cedar -resume
