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

# some errors but launched more jobs
#                 /scratch (user celphin)              11T/20T            37k/1000k

#------------------------------------
# After run for future runs edit config for BismarkIndex
nano /home/celphin/scratch/Dryas/methylseq/nf-core-methylseq_2.7.0/2_7_0/nextflow.config

params{

    // Input options
    input                      = '/home/celphin/scratch/Dryas/methylseq/input/input_files.csv'
    fasta                      = '/home/celphin/scratch/Dryas/methylseq/input/DoctH0_Main.fasta'
    //add after second run:
    bismark_index              = '/home/celphin/scratch/Dryas/methylseq/output/bismark/reference_genome/BismarkIndex/'
}

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

# ERROR ~ Error executing process > 'NFCORE_METHYLSEQ:METHYLSEQ:BISMARK:BISMARK_ALIGN (C_CAS
# S17C_576_175)'

# Caused by:
  # Process `NFCORE_METHYLSEQ:METHYLSEQ:BISMARK:BISMARK_ALIGN (C_CASS17C_576_175)` terminate
# d with an error exit status (255)

# Command executed:

  # if [ ! -f BismarkIndex/DoctH0_Main.fasta ]; then
      # ln -s $(readlink DoctH0_Main.fasta) BismarkIndex/DoctH0_Main.fasta;
  # fi

  # bismark \
      # -1 C_CASS17C_576_175_1_val_1.fq.gz -2 C_CASS17C_576_175_2_val_2.fq.gz \
      # --genome BismarkIndex \
      # --bam \
      # --bowtie2    --non_directional  --unmapped  --score_min L,0,-0.6 --multicore 2

  # cat <<-END_VERSIONS > versions.yml
  # "NFCORE_METHYLSEQ:METHYLSEQ:BISMARK:BISMARK_ALIGN":
      # bismark: $(echo $(bismark -v 2>&1) | sed 's/^.*Bismark Version: v//; s/Copyright.*$/
# /')
  # END_VERSIONS

# Command exit status:
  # 255

# Command error:

  # Created C -> T as well as G -> A converted versions of the FastQ file C_CASS17C_576_175_
# 2_val_2.fq.gz.temp.1 (31399235 sequences in total)

  # Input files are C_CASS17C_576_175_1_val_1.fq.gz.temp.1_C_to_T.fastq and C_CASS17C_576_17
# 5_1_val_1.fq.gz.temp.1_G_to_A.fastq and C_CASS17C_576_175_2_val_2.fq.gz.temp.1_C_to_T.fast
# q and C_CASS17C_576_175_2_val_2.fq.gz.temp.1_G_to_A.fastq (FastQ)
  # Now running 4 individual instances of Bowtie 2 against the bisulfite genome of /scratch/
# celphin/Dryas/methylseq/output/bismark/reference_genome/BismarkIndex/ with the specified o
# ptions: -q --score-min L,0,-0.6 --ignore-quals --no-mixed --no-discordant --dovetail --max
# ins 500

  # Now starting a Bowtie 2 paired-end alignment for CTread1GAread2CTgenome (reading in sequ
# ences from C_CASS17C_576_175_1_val_1.fq.gz.temp.1_C_to_T.fastq and C_CASS17C_576_175_2_val
# _2.fq.gz.temp.1_G_to_A.fastq, with the options: -q --score-min L,0,-0.6 --ignore-quals --n
# o-mixed --no-discordant --dovetail --maxins 500 --norc))
  # Chromosomal sequence could not be extracted for       A00977:110:HFTTCDSXY:2:1119:32461:
# 35978_1:N:0:CAACTCCA+GAATCCGT   DoctH0-10       1
  # Found first alignment:
  # A00977:110:HFTTCDSXY:2:1101:1868:1000_1:N:0:CAACTCCA+GAATCCGT/1       77      *       00
# *       *       0       0       TGAAGAATAATAAGAAGAGGGTTTTTATGGTTATAATAATTATTAGTGTTTTTTTTTT
# TTTTTTTATAAT    F:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FF:FFFFFFFFFFFFFFFFFFFFF:::FFF  YT
# :Z:UP
  # A00977:110:HFTTCDSXY:2:1101:1868:1000_2:N:0:CAACTCCA+GAATCCGT/2       141     *       00
# *       *       0       0       ATTATAAAAAAAAAAAAAAAAACACTAATAATTATTATAACCAAAAAAACCCTCTACT
# TATTATTCTACA    FF:,,F:FFFFFFFFFFFF:F:F,FFFFFF:F,FFF,FFFFF,:FFFFF:F:F,F,:,:::,FF,::,:F  YT
# :Z:UP
  # Now starting a Bowtie 2 paired-end alignment for GAread1CTread2GAgenome (reading in sequ
# ences from C_CASS17C_576_175_1_val_1.fq.gz.temp.1_G_to_A.fastq and C_CASS17C_576_175_2_val
# _2.fq.gz.temp.1_C_to_T.fastq, with the options: -q --score-min L,0,-0.6 --ignore-quals --n
# o-mixed --no-discordant --dovetail --maxins 500 --norc))

  # >>> Writing bisulfite mapping results to C_CASS17C_576_175_1_val_1.fq.gz.temp.1_bismark_
# bt2_pe.bam <<<

  # Unmapped sequences will be written to C_CASS17C_576_175_1_val_1.fq.gz.temp.1_unmapped_re
# ads_1.fq and C_CASS17C_576_175_2_val_2.fq.gz.temp.1_unmapped_reads_2.fq

  # Reading in the sequence files C_CASS17C_576_175_1_val_1.fq.gz.temp.1 and C_CASS17C_576_1
# 75_2_val_2.fq.gz.temp.1
  # Chromosomal sequence could not be extracted for       A00977:110:HFTTCDSXY:2:1105:21992:
# 18161_1:N:0:CAACTCCA+GAATCCGT   DoctH0-10       1
  # Chromosomal sequence could not be extracted for       A00977:110:HFTTCDSXY:2:1106:12771:
# 19351_1:N:0:CAACTCCA+GAATCCGT   DoctH0-12       2
  # Processed 1000000 sequence pairs so far

#-----------------------
# 22122_1:N:0:CAACTCCA+GAATCCGT   DoctH0-12       1
  # Chromosomal sequence could not be extracted for       A00977:110:HFTTCDSXY:2:1163:8793:1
# 4747_1:N:0:CAACTCCA+GAATCCGT    DoctH0-10       3
  # Chromosomal sequence could not be extracted for       A00977:110:HFTTCDSXY:2:1175:31675:
# 2895_1:N:0:CAACTCCA+GAATCCGT    DoctH0-12       101899
  # (ERR): bowtie2-align died with signal 9 (KILL)
  # (ERR): bowtie2-align died with signal 9 (KILL)
  # Use of uninitialized value $flag_2 in numeric eq (==) at /usr/local/bin/bismark line 328
# 9, <__ANONIO__> line 6300638.
  # Chromosome number extraction failed for *
  # (ERR): bowtie2-align died with signal 13 (PIPE)
  # (ERR): bowtie2-align died with signal 13 (PIPE)

# Work dir:
  # /scratch/celphin/Dryas/methylseq/work/cc/ce5d35131e327b810735d283760055

# Tip: when you have fixed the problem you can continue the execution adding the option `-re
# sume` to the run command line

 # -- Check '.nextflow.log' file for details
# ERROR ~ Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage
# /troubleshooting

 # -- Check '.nextflow.log' file for details




# Chromosome number extraction failed
# Signal 9 (KILL): This typically indicates that the process was killed by the operating system, 
# often due to exceeding resource limits (like memory). 


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

nextflow run nf-core-methylseq_2.7.0/2_7_0/  -profile singularity,cedar --max_memory 200GB -resume

executor >  slurm (109)
[-        ] NFCORE_METHYLSEQ:METHYLSEQ:CAT_FASTQ     -
[52/17c391] NFC…:METHYLSEQ:FASTQC (W_WILL7W_448_107) | 131 of 131, cached: 131 ✔
[97/69db15] NFC…THYLSEQ:TRIMGALORE (W_CASS7W_600_19) | 113 of 131, cached: 113
[91/f991c0] NFC…ARK:BISMARK_ALIGN (W_CASS10W_544_60) | 28 of 113, cached: 27, failed: 1
[e5/239c3b] NFC…OLS_SORT_ALIGNED (W_CASS17W_574_137) | 0 of 27
[e7/ef8c39] NFC…BISMARK_DEDUPLICATE (W_ALAS0W_7_263) | 0 of 27
[-        ] NFC…BISMARK:BISMARK_METHYLATIONEXTRACTOR -
[-        ] NFC…EQ:BISMARK:BISMARK_COVERAGE2CYTOSINE -
[-        ] NFC…SEQ:METHYLSEQ:BISMARK:BISMARK_REPORT -
[-        ] NFC…EQ:METHYLSEQ:BISMARK:BISMARK_SUMMARY -
[-        ] NFC…Q:BISMARK:SAMTOOLS_SORT_DEDUPLICATED -
[-        ] NFC…:BISMARK:SAMTOOLS_INDEX_DEDUPLICATED -
[-        ] NFC…E_METHYLSEQ:METHYLSEQ:QUALIMAP_BAMQC -
[-        ] NFC…_METHYLSEQ:METHYLSEQ:PRESEQ_LCEXTRAP -
[-        ] NFCORE_METHYLSEQ:METHYLSEQ:MULTIQC       | 0 of 1



#ERROR ~ Error executing process > 'NFCORE_METHYLSEQ:METHYLSEQ:BISMARK:BISMARK_ALIGN (C_LATJ_00C_187)'

#     47344041  celphin rrg-rieseber nf-NFCORE_METH  PD 8-00:00:00     1   12        N/A     72G  (ReqNodeNotAvail, Reserved for maintenance)


scontrol update jobid=47343921 TimeLimit=7-00:00:00


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

nextflow run nf-core-methylseq_2.7.0/2_7_0/  -profile singularity,cedar --max_memory 200GB -resume
