#########################################################
# Dryas DNA methylation calling
# nf-core Methylseq
# Mapping WGBS data to reference
# Calling methylated sites
#############################################################

# Try run to setup reference genome  in Bismarkindex/
# try on narval but slow due to allocation running low
# try on Cedar with high allocation but poor disk connection

tmux new-session -s Dryas
tmux attach-session -t Dryas

cd /home/celphin/scratch/Dryas/methylseq

#salloc -c1 --time 1:00:00 --mem 120000m --account rrg-rieseber-ac

cd /home/celphin/scratch/Dryas/methylseq
source nf-core-env/bin/activate
module load nextflow/24.04.4
module load apptainer/1.3.4

export NXF_SINGULARITY_CACHEDIR=/project/def-rieseber/NXF_SINGULARITY_CACHEDIR
nextflow run nf-core-methylseq_2.7.1/2_7_1/ -profile singularity,narval 

#export NXF_SINGULARITY_CACHEDIR=/project/rrg-rieseber-ac/NXF_SINGULARITY_CACHEDIR
#nextflow run nf-core-methylseq_2.7.1/2_7_1/ -profile singularity,cedar

#######################
# Check pipeline info to get run times

# task_id  hash            native_id       name                                                        status           exit     submit                 duration         realtime        %cpu     peak_rss     peak_vmem       rchar   wchar
# 1       6c/07611a       48656697        NFCORE_METHYLSEQ:PREPARE_GENOME:BISMARK_GENOMEPREPARATION (BismarkIndex/DoctH0_Main.fasta)     COMPLETED       0       2024-11-18 14:28:21.524 46m 48s        35m 11s 43.4%   1.1 GB  1.3 GB  1.9 GB  1.1 GB
# 9       98/4c70d0       48656718        NFCORE_METHYLSEQ:METHYLSEQ:TRIMGALORE (C_ALAS0C_13_254)        COMPLETED       0       2024-11-18 14:28:32.612 57m 56s         30m 23s        600.2%  591.7 MB      3.4 GB   202.3 GB        195.7 GB
# 7       c7/1fb64b       48656713        NFCORE_METHYLSEQ:METHYLSEQ:TRIMGALORE (C_ALAS0C_12_256)        COMPLETED       0       2024-11-18 14:28:26.221 1h 3m 25s       35m 31s        576.2%  587.3 MB       3.4 GB  230.7 GB        223.2 GB
# 6       11/c631b8       48656693        NFCORE_METHYLSEQ:METHYLSEQ:TRIMGALORE (C_ALAS0C_10_246)        COMPLETED       0       2024-11-18 14:28:20.968 1h 7m 30s       1h 1m 20s       466.5% 584.7 MB       3.4 GB  286.2 GB        276.8 GB
# 5       f5/e63d65       48656700        NFCORE_METHYLSEQ:METHYLSEQ:TRIMGALORE (Asian1_F112573)         COMPLETED        0       2024-11-18 14:28:22.193 1h 14m 9s       46m 22s        437.8%  590.9 MB      3.4 GB   282 GB         274.6 GB
# 16      88/18ee78       48656752        NFCORE_METHYLSEQ:METHYLSEQ:TRIMGALORE (C_ALAS0C_4_240)         COMPLETED        0       2024-11-18 14:28:48.923 1h 16m 46s      28m 54s        537.0%  589.5 MB      3.4 GB   204.7 GB        198.1 GB
# 15      14/3e0097       48656744        NFCORE_METHYLSEQ:METHYLSEQ:TRIMGALORE (C_ALAS0C_3_258)         COMPLETED        0       2024-11-18 14:28:45.951 1h 19m 44s      35m 59s         629.2%  589.2 MB      3.4 GB   227.8 GB        220.5 GB
# 4       e9/23a371       48656711        NFCORE_METHYLSEQ:METHYLSEQ:FASTQC     (C_ALAS0C_12_256)       -COMPLETED       0       2024-11-18 14:28:25.287 3h 24m 39s      9m 30s           197.8%  7.3 GB         40.3 GB  7.4 GB         4.5 MB
# 2       e0/0ba4bf       48656705        NFCORE_METHYLSEQ:METHYLSEQ:FASTQC     (Asian1_F112573)        -COMPLETED       0       2024-11-18 14:28:23.221 3h 21m 10s      10m 26s          197.2%  6.8 GB         40.3 GB  7.4 GB         4.5 MB



####################################################
# Make copy of reference and now run without running ref
 
cp -r /home/celphin/scratch/Dryas/methylseq/output/bismark/reference_genome/BismarkIndex* /home/celphin/scratch/Dryas/methylseq/input/reference/

#------------------------------------
# After run for future runs edit config for BismarkIndex
nano /home/celphin/scratch/Dryas/methylseq/nf-core-methylseq_2.7.1/2_7_1/nextflow.config

params{

    // Input options
    input                      = '/home/celphin/scratch/Dryas/methylseq/input/input_files.csv'
    fasta                      = '/home/celphin/scratch/Dryas/methylseq/input/reference/DoctH0_Main.fasta'
    //add after second run:
    bismark_index            = '/home/celphin/scratch/Dryas/methylseq/input/reference/BismarkIndex/'
}

# Change 
save_reference             = true to false
save_trimmed               = true to false

#----------------
nano /home/celphin/scratch/Dryas/methylseq/nf-core-methylseq_2.7.1/2_7_1/conf/base.config

# change 5 hours to 1 hours
# change 7 hours to 2 hours
# change 4 days to 1 days

# change to 4 attempts as Error stratedgy

# add memory to alignment process 72G to 120Gb

#######################################################
# start pipeline again with reference

cd /home/celphin/scratch/Dryas/methylseq

# narval3
tmux new-session -s Dryas1
tmux attach-session -t Dryas1

#salloc -c1 --time 18:00:00 --mem 128000m --account def-rieseber

cd /home/celphin/scratch/Dryas/methylseq

source nf-core-env/bin/activate
module load nextflow/24.04.4
module load apptainer/1.3.4
export NXF_SINGULARITY_CACHEDIR=/project/def-rieseber/NXF_SINGULARITY_CACHEDIR

nextflow run nf-core-methylseq_2.7.1/2_7_1/  -profile singularity,narval 

# executor >  slurm (333)
# [-        ] NFCORE_METHYLSEQ:METHYLSEQ:CAT_FASTQ           -
# [b4/23fe1c] NFC…HYLSEQ:METHYLSEQ:FASTQC (W_WILL5W_421_154) | 131 of 131 ✔
# [6c/5352bd] NFC…LSEQ:METHYLSEQ:TRIMGALORE (W_SVAL8W_8_272) | 131 of 131, cached: 60 ✔
# [f7/1537df] NFC…SEQ:BISMARK:BISMARK_ALIGN (W_SVAL8W_8_272) | 51 of 51, failed: 12
# [-        ] NFC…EQ:METHYLSEQ:BISMARK:SAMTOOLS_SORT_ALIGNED -
# [-        ] NFC…LSEQ:METHYLSEQ:BISMARK:BISMARK_DEDUPLICATE -
# [-        ] NFC…YLSEQ:BISMARK:BISMARK_METHYLATIONEXTRACTOR -
# [-        ] NFC…ETHYLSEQ:BISMARK:BISMARK_COVERAGE2CYTOSINE -
# [-        ] NFC…METHYLSEQ:METHYLSEQ:BISMARK:BISMARK_REPORT -
# [-        ] NFC…ETHYLSEQ:METHYLSEQ:BISMARK:BISMARK_SUMMARY -
# [-        ] NFC…THYLSEQ:BISMARK:SAMTOOLS_SORT_DEDUPLICATED -
# [-        ] NFC…HYLSEQ:BISMARK:SAMTOOLS_INDEX_DEDUPLICATED -
# [-        ] NFCORE_METHYLSEQ:METHYLSEQ:QUALIMAP_BAMQC      -
# [-        ] NFCORE_METHYLSEQ:METHYLSEQ:PRESEQ_LCEXTRAP     -
# [-        ] NFCORE_METHYLSEQ:METHYLSEQ:MULTIQC             -


#---------------------------
# [agitated_spence]
# ERROR ~ Error executing process > 'NFCORE_METHYLSEQ:METHYLSEQ:BISMARK:BISMARK_ALIGN (C_FERT5C_1F_97)'

# Caused by:
  # Process `NFCORE_METHYLSEQ:METHYLSEQ:BISMARK:BISMARK_ALIGN (C_FERT5C_1F_97)` terminated with an error exit status (255)

cd /lustre07/scratch/celphin/Dryas/methylseq/work/dc/66383cd334b9b52ce22e71be15abcd
more .command.log

# Chromosome number extraction failed for *
# (ERR): bowtie2-align died with signal 13 (PIPE)
# (ERR): bowtie2-align died with signal 13 (PIPE)
# (ERR): bowtie2-align died with signal 13 (PIPE)
# slurmstepd: error: Detected 1 oom_kill event in StepId=37167039.batch. Some of the step tasks ha
# ve been OOM Killed.

# C_SVAL8C_8_275 (37167079) 
# C_ALAS0C_10_246

# check various output file sizes - maybe not completing due to short on memory

# final
# -rw-r-----. 1 celphin celphin 11112756284 Nov 24 10:49 Asian1_F112573_1_val_1_bismark_bt2_pe.bam
# -rw-r-----. 1 celphin celphin  9762814771 Nov 24 23:29 C_ALAS0C_18_229_1_val_1_bismark_bt2_pe.bam

#--------------------
# -rw-r----- 1 celphin celphin  9661453352 Nov 14 16:37 Asian1_F112573_1_val_1_bismark_bt2_pe.bam
# -rw-r----- 1 celphin celphin 10596301569 Nov 14 16:36 C_ALAS0C_10_246_1_val_1_bismark_bt2_pe.bam
# -rw-r----- 1 celphin celphin  7331628014 Nov 14 16:38 C_ALAS0C_12_256_1_val_1_bismark_bt2_pe.bam
# -rw-r----- 1 celphin celphin  9760597564 Nov 14 16:38 C_ALAS0C_18_229_1_val_1_bismark_bt2_pe.bam
# -rw-r----- 1 celphin celphin  9529231898 Nov 14 16:40 C_ALAS0C_19_261_1_val_1_bismark_bt2_pe.bam
# -rw-r----- 1 celphin celphin  6495125640 Nov 14 16:39 C_ALAS0C_4_240_1_val_1_bismark_bt2_pe.bam
# -rw-r----- 1 celphin celphin 10758344752 Nov 14 16:41 C_ALAS0C_5_238_1_val_1_bismark_bt2_pe.bam
# -rw-r----- 1 celphin celphin 10475024815 Nov 14 16:43 C_ALAS_00C_231_1_val_1_bismark_bt2_pe.bam
# -rw-r----- 1 celphin celphin  8055652897 Nov 14 16:43 C_CASS10C_548_144_1_val_1_bismark_bt2_pe.bam
# -rw-r----- 1 celphin celphin  8281149750 Nov 14 16:44 C_CASS8C_535_54_1_val_1_bismark_bt2_pe.bam
# -rw-r----- 1 celphin celphin  6906953612 Nov 14 16:45 C_DRY5C_28_92_1_val_1_bismark_bt2_pe.bam
# -rw-r----- 1 celphin celphin  8220094753 Nov 14 16:46 C_FERT31C_15F_170_1_val_1_bismark_bt2_pe.bam
# -rw-r----- 1 celphin celphin  7926455786 Nov 14 16:46 C_FERT39C_20F_71_1_val_1_bismark_bt2_pe.bam
# -rw-r----- 1 celphin celphin 10177843181 Nov 14 16:47 C_FERT5C_1F_97_1_val_1_bismark_bt2_pe.bam
# -rw-r----- 1 celphin celphin 10810601227 Nov 14 16:49 C_LATD1C_4_223_1_val_1_bismark_bt2_pe.bam
# -rw-r----- 1 celphin celphin  7524951589 Nov 14 16:48 C_LATD2C_1_203_1_val_1_bismark_bt2_pe.bam
# -rw-r----- 1 celphin celphin  9650252163 Nov 14 16:49 C_LATD5C_20_199_1_val_1_bismark_bt2_pe.bam
# -rw-r----- 1 celphin celphin  9520776418 Nov 14 16:50 C_MEAD1C_446_33_1_val_1_bismark_bt2_pe.bam
# -rw-r----- 1 celphin celphin  8293634258 Nov 14 16:51 C_MEAD6C_468_22_1_val_1_bismark_bt2_pe.bam
# -rw-r----- 1 celphin celphin  9637296980 Nov 14 16:51 C_MEAD7C_473_95_1_val_1_bismark_bt2_pe.bam
# -rw-r----- 1 celphin celphin  6342303744 Nov 14 16:52 C_SVAL16C_16_276_1_val_1_bismark_bt2_pe.bam

# #----------
# -rw-r-----. 1 celphin celphin  9727700151 Nov 23 12:22 Asian1_F112573_1_val_1_bismark_bt2_pe.bam
# -rw-r-----. 1 celphin celphin 11082997345 Nov 23 11:50 C_ALAS_00C_227_1_val_1_bismark_bt2_pe.bam
# -rw-r-----. 1 celphin celphin 10482924112 Nov 23 11:04 C_ALAS_00C_231_1_val_1_bismark_bt2_pe.bam
# -rw-r-----. 1 celphin celphin  7441835006 Nov 23 09:26 C_ALAS0C_13_254_1_val_1_bismark_bt2_pe.bam
# -rw-r-----. 1 celphin celphin  9762131252 Nov 23 10:32 C_ALAS0C_18_229_1_val_1_bismark_bt2_pe.bam
# -rw-r-----. 1 celphin celphin  8950226082 Nov 23 10:43 C_ALAS0C_19_261_1_val_1_bismark_bt2_pe.bam
# -rw-r-----. 1 celphin celphin  8431881588 Nov 23 09:15 C_ALAS0C_3_258_1_val_1_bismark_bt2_pe.bam
# -rw-r-----. 1 celphin celphin  7756507998 Nov 23 09:06 C_ALAS0C_4_240_1_val_1_bismark_bt2_pe.bam
# -rw-r-----. 1 celphin celphin 12091358709 Nov 23 11:54 C_CASS_09C_541_1_val_1_bismark_bt2_pe.bam
# -rw-r-----. 1 celphin celphin  8095671583 Nov 23 12:34 C_CASS17C_576_175_1_val_1_bismark_bt2_pe.bam
# -rw-r-----. 1 celphin celphin  8008259941 Nov 23 11:27 C_CASS4C_524_4_1_val_1_bismark_bt2_pe.bam
# -rw-r-----. 1 celphin celphin  8755044955 Nov 23 11:05 C_CASS5C_529_159_1_val_1_bismark_bt2_pe.bam
# -rw-r-----. 1 celphin celphin  8283471388 Nov 23 10:54 C_CASS8C_535_54_1_val_1_bismark_bt2_pe.bam
# -rw-r-----. 1 celphin celphin  8149344809 Nov 23 12:09 C_DRY2C_10_52_1_val_1_bismark_bt2_pe.bam
# -rw-r-----. 1 celphin celphin  8201724757 Nov 23 12:06 C_DRY4C_23_82_1_val_1_bismark_bt2_pe.bam
# -rw-r-----. 1 celphin celphin  7556929481 Nov 23 10:18 C_DRY5C_28_92_1_val_1_bismark_bt2_pe.bam
# -rw-r-----. 1 celphin celphin  8085317491 Nov 23 10:52 C_DRY9C_53_149_1_val_1_bismark_bt2_pe.bam
# -rw-r-----. 1 celphin celphin 12656804860 Nov 23 13:00 C_FERT13C_7F_112_1_val_1_bismark_bt2_pe.bam
# -rw-r-----. 1 celphin celphin  8222461474 Nov 23 11:56 C_FERT31C_15F_170_1_val_1_bismark_bt2_pe.bam
# -rw-r-----. 1 celphin celphin  7928228326 Nov 23 10:44 C_FERT39C_20F_71_1_val_1_bismark_bt2_pe.bam
# -rw-r-----. 1 celphin celphin  9091686390 Nov 23 12:05 C_FERT5C_1F_97_1_val_1_bismark_bt2_pe.bam
# -rw-r-----. 1 celphin celphin 10812002267 Nov 23 12:06 C_LATD1C_4_223_1_val_1_bismark_bt2_pe.bam
# -rw-r-----. 1 celphin celphin  8536382223 Nov 23 11:00 C_LATD2C_1_203_1_val_1_bismark_bt2_pe.bam
# -rw-r-----. 1 celphin celphin 11068553783 Nov 23 12:10 C_LATD2C_6_198_1_val_1_bismark_bt2_pe.bam
# -rw-r-----. 1 celphin celphin  9111750794 Nov 23 13:23 C_LATD2C_7_209_1_val_1_bismark_bt2_pe.bam
# -rw-r-----. 1 celphin celphin  9710697680 Nov 23 11:35 C_LATD4C_3_196_1_val_1_bismark_bt2_pe.bam
# -rw-r-----. 1 celphin celphin  9304470441 Nov 23 13:39 C_LATD5C_20_199_1_val_1_bismark_bt2_pe.bam
# -rw-r-----. 1 celphin celphin  9316620796 Nov 23 13:18 C_LATD5C_2_201_1_val_1_bismark_bt2_pe.bam
# -rw-r-----. 1 celphin celphin  9351360838 Nov 23 12:35 C_LATD5C_5_191_1_val_1_bismark_bt2_pe.bam
# -rw-r-----. 1 celphin celphin 12177072907 Nov 23 14:20 C_LATJ_00C_187_1_val_1_bismark_bt2_pe.bam
# -rw-r-----. 1 celphin celphin  9848917818 Nov 23 14:24 C_MEAD2C_451_76_1_val_1_bismark_bt2_pe.bam
# -rw-r-----. 1 celphin celphin  7371305657 Nov 23 13:09 C_MEAD6C_468_22_1_val_1_bismark_bt2_pe.bam
# -rw-r-----. 1 celphin celphin  8275295038 Nov 23 12:09 C_SVAL12C_12_273_1_val_1_bismark_bt2_pe.bam
# -rw-r-----. 1 celphin celphin  9916882420 Nov 23 13:18 C_SVAL16C_16_276_1_val_1_bismark_bt2_pe.bam
# -rw-r-----. 1 celphin celphin  8713308994 Nov 23 12:59 C_SVAL49C_49_278_1_val_1_bismark_bt2_pe.bam
# -rw-r-----. 1 celphin celphin  8987550736 Nov 23 13:21 C_SVAL8C_8_275_1_val_1_bismark_bt2_pe.bam
# -rw-r-----. 1 celphin celphin  7995768054 Nov 23 14:17 C_WILL1C_406_152_1_val_1_bismark_bt2_pe.bam
# -rw-r-----. 1 celphin celphin  8798756353 Nov 23 14:28 C_WILL3C_414_100_1_val_1_bismark_bt2_pe.bam
# -rw-r-----. 1 celphin celphin 10145793519 Nov 23 15:02 C_WILL5C_422_31_1_val_1_bismark_bt2_pe.ba

# #-----------------
# -rw-r----- 1 celphin celphin 4541786807 Nov 17 01:58 C_CASS4C_524_4.sorted.bam
# -rw-r----- 1 celphin celphin 8006069754 Nov 16 21:46 C_CASS4C_524_4_1_val_1_bismark_bt2_pe.bam


#---------------------------------
# check log of past runs

nextflow log
# TIMESTAMP               DURATION        RUN NAME        STATUS  REVISION ID     SESSION ID                              COMMAND                                                                 
# 2024-11-21 01:55:05     -               awesome_agnesi  -       78f1f08d90      ceaa73c5-ac10-4639-97e5-b8221f81bd20    nextflow run nf-core-methylseq_2.7.1/2_7_1/ -profile singularity,narval 
# 2024-11-22 02:05:41     1d 1h 47m 4s    cheeky_mercator ERR     78f1f08d90      ceaa73c5-ac10-4639-97e5-b8221f81bd20    nextflow run nf-core-methylseq_2.7.1/2_7_1/ -profile singularity,narval -resume
# 2024-11-23 03:54:02     11h 20m 46s     prickly_mcnulty ERR     78f1f08d90      ceaa73c5-ac10-4639-97e5-b8221f81bd20    nextflow run nf-core-methylseq_2.7.1/2_7_1/ -profile singularity,narval -resume

#---------------------------------
# resume pipeline with more memory 

nano /home/celphin/scratch/Dryas/methylseq/nf-core-methylseq_2.7.1/2_7_1/conf/base.config

# add memory to alignment process 120G to 200Gb

#-----------------------------
#narval3
tmux new-session -s Dryas
tmux attach-session -t Dryas

cd /home/celphin/scratch/Dryas/methylseq
source nf-core-env/bin/activate
module load nextflow/24.04.4
module load apptainer/1.3.4
export NXF_SINGULARITY_CACHEDIR=/project/def-rieseber/NXF_SINGULARITY_CACHEDIR

nextflow run nf-core-methylseq_2.7.1/2_7_1/  -profile singularity,narval --max_memory 200GB -resume

# executor >  slurm (131)
# [-        ] NFCORE_METHYLSEQ:METHYLSEQ:CAT_FASTQ           -
# [b4/23fe1c] NFC…HYLSEQ:METHYLSEQ:FASTQC (W_WILL5W_421_154) | 131 of 131, cached: 131 ✔
# [75/0151c8] NFC…EQ:METHYLSEQ:TRIMGALORE (W_WILL5W_421_154) | 131 of 131, cached: 131 ✔
# [26/de19c6] NFC…Q:BISMARK:BISMARK_ALIGN (W_WILL5W_421_154) | 0 of 131
# [-        ] NFC…EQ:METHYLSEQ:BISMARK:SAMTOOLS_SORT_ALIGNED -
# [-        ] NFC…LSEQ:METHYLSEQ:BISMARK:BISMARK_DEDUPLICATE -
# [-        ] NFC…YLSEQ:BISMARK:BISMARK_METHYLATIONEXTRACTOR -
# [-        ] NFC…ETHYLSEQ:BISMARK:BISMARK_COVERAGE2CYTOSINE -
# [-        ] NFC…METHYLSEQ:METHYLSEQ:BISMARK:BISMARK_REPORT -
# [-        ] NFC…ETHYLSEQ:METHYLSEQ:BISMARK:BISMARK_SUMMARY -
# [-        ] NFC…THYLSEQ:BISMARK:SAMTOOLS_SORT_DEDUPLICATED -
# [-        ] NFC…HYLSEQ:BISMARK:SAMTOOLS_INDEX_DEDUPLICATED -
# [-        ] NFCORE_METHYLSEQ:METHYLSEQ:QUALIMAP_BAMQC      -
# [-        ] NFCORE_METHYLSEQ:METHYLSEQ:PRESEQ_LCEXTRAP     -
# [-        ] NFCORE_METHYLSEQ:METHYLSEQ:MULTIQC             -

# [-        ] NFCORE_METHYLSEQ:METHYLSEQ:CAT_FASTQ         -
# [b4/23fe1c] NFC…LSEQ:METHYLSEQ:FASTQC (W_WILL5W_421_154) | 131 of 131, cached: 131 ✔
# [75/0151c8] NFC…:METHYLSEQ:TRIMGALORE (W_WILL5W_421_154) | 131 of 131, cached: 131 ✔
# [af/b8213e] NFC…:BISMARK:BISMARK_ALIGN (W_WILL1W_403_67) | 131 of 131 ✔
# [35/ad22d6] NFC…:SAMTOOLS_SORT_ALIGNED (W_WILL1W_403_67) | 129 of 174, failed: 43, retries: 43[a1/552079] NFC…RK:BISMARK_DEDUPLICATE (W_WILL1W_403_67) | 5 of 131
# [35/d0e8f1] NFC…K_METHYLATIONEXTRACTOR (C_ALAS0C_12_256) | 0 of 5
# [-        ] NFC…HYLSEQ:BISMARK:BISMARK_COVERAGE2CYTOSINE -
# [-        ] NFC…THYLSEQ:METHYLSEQ:BISMARK:BISMARK_REPORT -
# [-        ] NFC…HYLSEQ:METHYLSEQ:BISMARK:BISMARK_SUMMARY -
# [ba/b97c56] NFC…OOLS_SORT_DEDUPLICATED (C_ALAS0C_12_256) | 0 of 5
# [-        ] NFC…LSEQ:BISMARK:SAMTOOLS_INDEX_DEDUPLICATED -
# [-        ] NFCORE_METHYLSEQ:METHYLSEQ:QUALIMAP_BAMQC    -
# [a3/5b7617] NFC…THYLSEQ:PRESEQ_LCEXTRAP (C_LATJ_00C_187) | 81 of 90, failed: 4, retries: 4
# [-        ] NFCORE_METHYLSEQ:METHYLSEQ:MULTIQC           -


#---------------------------------
# update 
 narval {
        process.clusterOptions = "--account def-henryg"
        max_memory='249G'
        max_cpu=64
        max_time='168h'
    }

#---------------------------
# resume pipeline 
 
#Cedar1
tmux new-session -s Dryas
tmux attach-session -t Dryas

cd /home/celphin/scratch/Dryas/methylseq
source nf-core-env/bin/activate
module load nextflow/24.04.4
module load apptainer/1.3.4
export NXF_SINGULARITY_CACHEDIR=/project/def-rieseber/NXF_SINGULARITY_CACHEDIR

nextflow run nf-core-methylseq_2.7.1/2_7_1/  -profile singularity,narval --max_memory 200GB -resume






#######################################
# updating job times and memory

squeue -u celphin -t PD -o "%.18i"

for job_id in $(squeue -u celphin -t PD -o "%.18i %.10l" --noheader| awk '$2 == "5:00:00" {print $1}'); do
  scontrol update JobId=$job_id TimeLimit=1:00:00
done

for job_id in $(squeue -u celphin -t PD -o "%.18i %.10l" --noheader| awk '$2 == "7:00:00" {print $1}'); do
  scontrol update JobId=$job_id TimeLimit=2:00:00
done

for job_id in $(squeue -u celphin -t PD -o "%.18i %.10l" --noheader| awk '$2 == "4-00:00:00" {print $1}'); do
  scontrol update JobId=$job_id TimeLimit=1-00:00:00
done


for job_id in $(squeue -u celphin -t PD -o "%.18i %.10l" --noheader| awk '$2 == "2:00:00" {print $1}'); do
  scontrol update JobId=$job_id MinMemoryNode=72000
done

for job_id in $(squeue -u celphin -t PD -o "%.18i %.10m" --noheader | awk '$2 == "72000M" {print $1}'); do
  scancel $job_id 
done

for job_id in $(squeue -u celphin -t PD -o "%.18i %.10m" --noheader | awk '$2 == "200G" {print $1}'); do
  scontrol update JobId=$job_id MinMemoryNode=150000
done

for job_id in $(squeue -u celphin -t PD -o "%.18i %.10m" --noheader | awk '$2 == "200G" {print $1}'); do
  scontrol update JobId=$job_id MinMemoryNode=150000
done

for job_id in $(squeue -u celphin -t PD -o "%.18i %.10l" --noheader| awk '$2 == "1-00:00:00" {print $1}'); do
  scontrol update JobId=$job_id TimeLimit=12:00:00
done

for job_id in $(squeue -u celphin -t PD -o "%.18i %.10m" --noheader | awk '$2 == "120G" {print $1}'); do
  scontrol update JobId=$job_id MinMemoryNode=200000
done

for job_id in $(squeue -u celphin -t PD -o "%.18i %.10l" --noheader| awk '$2 == "12:00:00" {print $1}'); do
  scontrol update JobId=$job_id TimeLimit=10:00:00
done


for job_id in $(squeue -u celphin -t PD --noheader| awk '{print $1}'); do
  scontrol update JobId=$job_id Account=def-henryg
done

for job_id in $(squeue -u celphin -t PD -o "%.18i %.10l" --noheader| awk '$2 == "1:00:00" {print $1}'); do
  scontrol update JobId=$job_id TimeLimit=2:50:00
done

scontrol update JobId=37188221 TimeLimit=1-00:00:00
scontrol update JobId=37272312 TimeLimit=1-00:00:00


