########################################
# Dryas differential gene expression data
# Mapping to transcriptome - only get mRNA mapping here
# Sept 2023 - Sweden, Feb 2023 other sites
#####################################

#################################
# check RNAseq data
# https://www.frontiersin.org/articles/10.3389/fgene.2021.649619/full 
# https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0881-8/figures/1
# https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/lessons/02_assessing_quality.html

# quality check reads
# FASTQC

salloc -c32 --time 12:50:00 --mem 120000m --account def-henryg

cd /home/celphin/projects/rpp-rieseber/celphin/Dryas/RNAseq_analysis/input/fastq_files/

module load StdEnv/2020
module load fastqc/0.11.9

fastqc -t 32 *.fastq

mv *fastqc* /home/celphin/projects/rpp-rieseber/celphin/Dryas/RNAseq_analysis/input/fastq_files/fastqc/

# on local machine
cd /home/Owner/Desktop

scp -v celphin@cedar.computecanada.ca:/home/celphin/projects/rpp-rieseber/celphin/Dryas/RNAseq_analysis/input/fastq_files/fastqc/* .


##############################################
# data analysis
# https://github.com/owensgl/biol525D/tree/master/Topic_6 

cd /home/celphin/projects/rpp-rieseber/celphin/Dryas/RNAseq_raw_data
# Data I1, I2, R1, R2

tmux new-session -s RNA
tmux attach-session -t RNA

# salloc -c32 --time 02:50:00 --mem 120000m --account def-henryg
# module spider bowtie

module load StdEnv/2020
module load bowtie2/2.4.4
module load r/4.2.1
module load samtools/1.15.1

cd /home/celphin/projects/rpp-rieseber/celphin/Other/
tar -xvzf biol525d.tar.gz

cd /home/celphin/projects/rpp-rieseber/celphin/Other/BIOL525d/cassandra.elphinstone/biol525d/Topic_6/scripts
cp RSEM-1.2.31.tar.gz /home/celphin/projects/rpp-rieseber/celphin/Dryas/RNAseq_analysis/

#--------------------------------
# install program

cd /home/celphin/projects/rpp-rieseber/celphin/Dryas/RNAseq_analysis/

tar -xzf RSEM-1.2.31.tar.gz

cd RSEM-1.2.31
make
make install

install -d /usr/local/bin /usr/local/bin/samtools-1.3
/cvmfs/soft.computecanada.ca/gentoo/2020/usr/bin/install: cannot change permissions of ‘/usr/local/bin/samtools-1.3’: No such file or directory
make: *** [Makefile:162: install] Error 1

#------------------------------------
# copy over reference
# DoctH0_Main.fasta
# DoctH0.AED_0.6_protein.fasta
# FINAL_DoctH0.AED_0.6.sorted.gff3
# Dryasoct_interproscan_edited.tsv
# Dryasoct_GO_mappings.ermineJ.txt


#------------------------------
# try https://github.com/gpertea/gffread
# to get to protein conversion

# Install
 cd /home/celphin/scratch/
  git clone https://github.com/gpertea/gffread
  cd gffread
  make release

cd /home/celphin/projects/rrg-rieseber-ac/rpp-rieseber/celphin/Dryas/RNAseq_analysis/reference

/home/celphin/scratch/Annotation/gffread/gffread -h

# to run proteins
/home/celphin/scratch/gffread/gffread -g DoctH0_Main.fasta FINAL_DoctH0.AED_0.6.sorted.gff3 -y Dry-octo-H0_proteins.fa

# to run proteins
/home/celphin/scratch/gffread/gffread -g DoctH0_Main.fasta FINAL_DoctH0.AED_0.6.sorted.gff3 -x Dry-octo-H0_cds.fa

#-----------------------------------
# Prep reference
# http://deweylab.github.io/RSEM/rsem-prepare-reference.html
# needs fasta and GFF file

module load StdEnv/2020
module load bowtie/1.3.0
module load r/4.2.1
module load samtools/1.15.1

cd /home/celphin/projects/rpp-rieseber/celphin/Dryas/RNAseq_analysis

cp /home/celphin/projects/rpp-rieseber/celphin/Dryas/Reference_genomes/gff/Dryas_octopetala_H1.gff3 .
cp /home/celphin/projects/rpp-rieseber/celphin/Dryas/Reference_genomes/genomes/Dryas_octopetala_H1.supercontigs.fa .
cp /home/celphin/projects/rpp-rieseber/celphin/Dryas/Reference_genomes/fasta/Dryas_octopetala_H1.transcript.fa .

/home/celphin/projects/rpp-rieseber/celphin/Dryas/RNAseq_analysis/RSEM-1.2.31/rsem-prepare-reference \
Dryas_octopetala_H1.transcript.fa \
Dryas_octopetala_H1_rnaseq

module load StdEnv/2020
module load bowtie2/2.4.4
module load r/4.2.1
module load samtools/1.15.1

bowtie2-build -f Dryas_octopetala_H1.transcript.fa Dryas_octopetala_H1_rnaseq

grep ">" Dryas_octopetala_H1.transcript.fa | wc -l

#----------------------------
# prep another reference using gff3 file

/home/celphin/projects/rpp-rieseber/celphin/Dryas/RNAseq_analysis/RSEM-1.2.31/rsem-prepare-reference \
--gff3 Dryas_octopetala_H1.gff3 \
--gff3-RNA-patterns mRNA,rRNA \
Dryas_octopetala_H1.supercontigs.fa \
./gff3ref/Dryas_octopetala_H1_gff3ref

bowtie2-build -f Dryas_octopetala_H1.supercontigs.fa ./gff3ref/Dryas_octopetala_H1_gff3ref

#---------------------------
# try mapping
# http://deweylab.github.io/RSEM/README.html 
# http://deweylab.github.io/RSEM/rsem-calculate-expression.html

#--------------------------------------------------------
# taking too much time with full fasta - try again with fasta only of transcripts

#--------------------------------------------------
tmux new-session -s RNA
tmux attach-session -t RNA

salloc -c32 --time 2:50:00 --mem 120000m --account rrg-rieseber-ac

module load StdEnv/2020
module load bowtie2/2.4.4
module load r/4.2.1
module load samtools/1.15.1

/home/celphin/projects/rpp-rieseber/celphin/Dryas/RNAseq_analysis/RSEM-1.2.31/rsem-calculate-expression \
--bowtie2 -p 32 --paired-end \
/home/celphin/projects/rpp-rieseber/celphin/Dryas/RNAseq_analysis/input/NS.1889.002.NEBNext_dual_i7_128---NEBNext_dual_i5_128.SE1W_R1.fastq.gz \
/home/celphin/projects/rpp-rieseber/celphin/Dryas/RNAseq_analysis/input/NS.1889.002.NEBNext_dual_i7_128---NEBNext_dual_i5_128.SE1W_R2.fastq.gz \
Dryas_octopetala_H1_rnaseq \
/home/celphin/projects/rpp-rieseber/celphin/Dryas/RNAseq_analysis/output/W_SE1

# worked!!

# SE1W
# gene_id                transcript_id(s)       length  effective_length      expected_count  TPM     FPKM
# Do1_00107G00001V1.1     Do1_00107G00001V1.1     1578.00 1385.75 41.94   0.58    0.92
# Do1_00107G00002V1.1     Do1_00107G00002V1.1     279.00  92.03   0.00    0.00    0.00
# Do1_00107G00003V1.1     Do1_00107G00003V1.1     298.00  109.59  29.49   5.13    8.14
# Do1_00107G00004V1.1     Do1_00107G00004V1.1     267.00  81.21   64.38   15.13   24.00
# Do1_00107G00005V1.1     Do1_00107G00005V1.1     210.00  34.17   0.00    0.00    0.00
# Do1_00107G00006V1.1     Do1_00107G00006V1.1     726.00  533.77  324.36  11.60   18.39
# Do1_00107G00007V1.1     Do1_00107G00007V1.1     1857.00 1664.75 263.23  3.02    4.79
# Do1_00107G00008V1.1     Do1_00107G00008V1.1     358.00  167.16  0.00    0.00    0.00
# Do1_00107G00009V1.1     Do1_00107G00009V1.1     326.00  136.15  0.00    0.00    0.00
# Do1_00107G00010V1.1     Do1_00107G00010V1.1     205.00  30.57   0.00    0.00    0.00
# Do1_00107G00011V1.1     Do1_00107G00011V1.1     278.00  91.12   3.00    0.63    1.00
# Do1_00107G00012V1.1     Do1_00107G00012V1.1     373.00  181.86  5.93    0.62    0.99
# Do1_00107G00013V1.1     Do1_00107G00013V1.1     246.00  62.89   0.00    0.00    0.00
# Do1_00107G00014V1.1     Do1_00107G00014V1.1     335.00  144.81  0.00    0.00    0.00
# Do1_00107G00015V1.1     Do1_00107G00015V1.1     358.00  167.16  0.00    0.00    0.00
# Do1_00107G00016V1.1     Do1_00107G00016V1.1     326.00  136.15  0.40    0.06    0.09
# Do1_00107G00017V1.1     Do1_00107G00017V1.1     214.00  37.13   2.00    1.03    1.63

# https://translational-medicine.biomedcentral.com/articles/10.1186/s12967-021-02936-w

TPM (transcripts per million)
FPKM (fragments per kilobase of transcript per million fragments mapped)
normalized counts using coefficient of variation, intraclass correlation coefficient, and cluster analysis

#------------------------------------------------------------
# try another individual

/home/celphin/projects/rpp-rieseber/celphin/Dryas/RNAseq_analysis/RSEM-1.2.31/rsem-calculate-expression \
--bowtie2 -p 32 --paired-end \
/home/celphin/projects/rpp-rieseber/celphin/Dryas/RNAseq_analysis/input/NS.1889.003.NEBNext_dual_i7_E8---NEBNext_dual_i5_E8.L_C11_R1.fastq.gz \
/home/celphin/projects/rpp-rieseber/celphin/Dryas/RNAseq_analysis/input/NS.1889.003.NEBNext_dual_i7_E8---NEBNext_dual_i5_E8.L_C11_R2.fastq.gz \
Dryas_octopetala_H1_rnaseq \
/home/celphin/projects/rpp-rieseber/celphin/Dryas/RNAseq_analysis/output/C_L11


# L_C11
# gene_id transcript_id(s)        length  effective_length        expected_count  TPM     FPKM
# Do1_00107G00001V1.1     Do1_00107G00001V1.1     1578.00 1382.17 32.81   0.31    0.49
# Do1_00107G00002V1.1     Do1_00107G00002V1.1     279.00  89.17   3.00    0.44    0.69
# Do1_00107G00003V1.1     Do1_00107G00003V1.1     298.00  106.55  153.09  18.70   29.61
# Do1_00107G00004V1.1     Do1_00107G00004V1.1     267.00  78.47   210.04  34.83   55.16
# Do1_00107G00005V1.1     Do1_00107G00005V1.1     210.00  32.60   0.00    0.00    0.00
# Do1_00107G00006V1.1     Do1_00107G00006V1.1     726.00  530.18  740.04  18.16   28.77
# Do1_00107G00007V1.1     Do1_00107G00007V1.1     1857.00 1661.17 413.62  3.24    5.13
# Do1_00107G00008V1.1     Do1_00107G00008V1.1     358.00  163.77  1.05    0.08    0.13
# Do1_00107G00009V1.1     Do1_00107G00009V1.1     326.00  132.91  0.00    0.00    0.00
# Do1_00107G00010V1.1     Do1_00107G00010V1.1     205.00  29.21   0.00    0.00    0.00
# Do1_00107G00011V1.1     Do1_00107G00011V1.1     278.00  88.26   0.00    0.00    0.00
# Do1_00107G00012V1.1     Do1_00107G00012V1.1     373.00  178.42  0.00    0.00    0.00
# Do1_00107G00013V1.1     Do1_00107G00013V1.1     246.00  60.42   0.00    0.00    0.00
# Do1_00107G00014V1.1     Do1_00107G00014V1.1     335.00  141.53  0.00    0.00    0.00
# Do1_00107G00015V1.1     Do1_00107G00015V1.1     358.00  163.77  0.00    0.00    0.00
# Do1_00107G00016V1.1     Do1_00107G00016V1.1     326.00  132.91  0.03    0.00    0.01
# Do1_00107G00017V1.1     Do1_00107G00017V1.1     214.00  35.40   2.22    0.82    1.29
# Do1_00107G00018V1.1     Do1_00107G00018V1.1     278.00  88.26   0.00    0.00    0.00
# Do1_00107G00019V1.1     Do1_00107G00019V1.1     374.00  179.40  35.00   2.54    4.02
# Do1_00107G00020V1.1     Do1_00107G00020V1.1     246.00  60.42   0.00    0.00    0.00
# Do1_00107G00021V1.1     Do1_00107G00021V1.1     431.00  235.66  349.31  19.29   30.55
# Do1_00107G00022V1.1     Do1_00107G00022V1.1     425.00  229.70  0.00    0.00    0.00


###############################
# should think about
# http://combine-australia.github.io/RNAseq-R/slides/RNASeq_filtering_qc.pdf
# https://www.biostars.org/p/83604/
# composition expression bias - normalize the data
# filtering out all very low count genes

##################################
# need generalized code to loop through all the RNAseq files

tmux new-session -s R
tmux attach-session -t R

salloc -c32 --time 23:50:00 --mem 120000m --account rpp-rieseber

module load StdEnv/2020
module load bowtie2/2.4.4
module load r/4.2.1
module load samtools/1.15.1

#------------------------------------
# make script for each individual

cd /home/celphin/projects/rpp-rieseber/celphin/Dryas/RNAseq_analysis/input/input_data
arr=($(ls *R1.fastq.gz))
${arr[2]}       

for x in ${arr[@]}; do   # better way to iterate over bash array elements; 5 loop iterations

outname=$(echo ${x} | sed 's/NS\.1889\.[0-9]*\.NEBNext_dual_i7..[0-9]*---NEBNext_dual_i[0-9]*..[0-9]*\.//' | sed 's/_R1.fastq.gz//') 
Rname=$(echo ${x} | sed 's/_R1.fastq.gz//')

echo ${outname}
echo ${Rname}

cat << EOF > ../RNAscripts/${outname}.sh
#!/bin/bash
#SBATCH --account=rpp-rieseber
#SBATCH --time=0-05:50:00 # 3:50 hours needed for all but files listed below
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=120000m
#SBATCH --output=slurm-%x.%j.out
#SBATCH --error=slurm-%x.%j.err

module load StdEnv/2020
module load bowtie2/2.4.4
module load r/4.2.1
module load samtools/1.15.1

cd /home/celphin/projects/rpp-rieseber/celphin/Dryas/RNAseq_analysis

srun /home/celphin/projects/rpp-rieseber/celphin/Dryas/RNAseq_analysis/RSEM-1.2.31/rsem-calculate-expression \
--bowtie2 -p 32 --paired-end \
/home/celphin/projects/rpp-rieseber/celphin/Dryas/RNAseq_analysis/input/input_data/${Rname}_R1.fastq.gz \
/home/celphin/projects/rpp-rieseber/celphin/Dryas/RNAseq_analysis/input/input_data/${Rname}_R2.fastq.gz \
/home/celphin/projects/rpp-rieseber/celphin/Dryas/RNAseq_analysis/reference/Dryas_octopetala_H1_rnaseq \
/home/celphin/projects/rpp-rieseber/celphin/Dryas/RNAseq_analysis/output/${outname}
EOF

done

#-----------------------
# batch submit scripts made above
cd /home/celphin/projects/rpp-rieseber/celphin/Dryas/RNAseq_analysis/input/RNAscripts/slurm_logs

# test
sbatch ../A_C11.sh
sbatch ../T_C12.sh
# takes ~5 hours

# try submitting all Toolik
cd /home/celphin/projects/rpp-rieseber/celphin/Dryas/RNAseq_analysis/input/RNAscripts/slurm_logs

for x in ../T_*.sh
do
sbatch ${x}
done

sq
          JOBID     USER      ACCOUNT           NAME  ST  TIME_LEFT NODES CPUS TRES_PER_N MIN_MEM NODELIST (REASON)
       60270168  celphin rpp-rieseber       A_C11.sh   R    3:37:34     1   32        N/A 120000M cdr2073 (None)
       60270909  celphin rpp-rieseber       T_C14.sh   R    3:45:49     1   32        N/A 120000M cdr1548 (None)
       60270910  celphin rpp-rieseber       T_C15.sh   R    3:45:49     1   32        N/A 120000M cdr1561 (None)
       60270911  celphin rpp-rieseber       T_C16.sh   R    3:45:49     1   32        N/A 120000M cdr1567 (None)
       60270912  celphin rpp-rieseber        T_C1.sh   R    3:45:49     1   32        N/A 120000M cdr1571 (None)
       60270913  celphin rpp-rieseber       T_C20.sh   R    3:45:49     1   32        N/A 120000M cdr1577 (None)
       60270914  celphin rpp-rieseber        T_C5.sh   R    3:45:49     1   32        N/A 120000M cdr2055 (None)
       60270915  celphin rpp-rieseber        T_C7.sh   R    3:45:49     1   32        N/A 120000M cdr2056 (None)
       60270917  celphin rpp-rieseber       T_W12.sh   R    3:45:49     1   32        N/A 120000M cdr2062 (None)
       60270918  celphin rpp-rieseber       T_W13.sh   R    3:45:49     1   32        N/A 120000M cdr2065 (None)
       60270919  celphin rpp-rieseber       T_W15.sh   R    3:45:49     1   32        N/A 120000M cdr2070 (None)
       60270920  celphin rpp-rieseber        T_W1.sh   R    3:45:49     1   32        N/A 120000M cdr2071 (None)
       60270921  celphin rpp-rieseber        T_W2.sh   R    3:45:49     1   32        N/A 120000M cdr2072 (None)
       60270923  celphin rpp-rieseber        T_W3.sh   R    3:45:49     1   32        N/A 120000M cdr2040 (None)
       60270924  celphin rpp-rieseber        T_W5.sh   R    3:45:49     1   32        N/A 120000M cdr2042 (None)
       60270925  celphin rpp-rieseber        T_W9.sh   R    3:45:49     1   32        N/A 120000M cdr2045 (None)
       60271227  celphin rpp-rieseber       T_C12.sh   R    3:49:51     1   32        N/A 120000M cdr1605 (None)
ls
slurm-A_C11.sh.60270168.err  slurm-T_C15.sh.60270910.out  slurm-T_C5.sh.60270914.err   slurm-T_W13.sh.60270918.out  slurm-T_W3.sh.60270923.err
slurm-A_C11.sh.60270168.out  slurm-T_C16.sh.60270911.err  slurm-T_C5.sh.60270914.out   slurm-T_W15.sh.60270919.err  slurm-T_W3.sh.60270923.out
slurm-T_C12.sh.60271227.err  slurm-T_C16.sh.60270911.out  slurm-T_C7.sh.60270915.err   slurm-T_W15.sh.60270919.out  slurm-T_W5.sh.60270924.err
slurm-T_C12.sh.60271227.out  slurm-T_C1.sh.60270912.err   slurm-T_C7.sh.60270915.out   slurm-T_W1.sh.60270920.err   slurm-T_W5.sh.60270924.out
slurm-T_C14.sh.60270909.err  slurm-T_C1.sh.60270912.out   slurm-T_W12.sh.60270917.err  slurm-T_W1.sh.60270920.out   slurm-T_W9.sh.60270925.err
slurm-T_C14.sh.60270909.out  slurm-T_C20.sh.60270913.err  slurm-T_W12.sh.60270917.out  slurm-T_W2.sh.60270921.err   slurm-T_W9.sh.60270925.out
slurm-T_C15.sh.60270910.err  slurm-T_C20.sh.60270913.out  slurm-T_W13.sh.60270918.err  slurm-T_W2.sh.60270921.out

#-------------------------------
# submit all Norway

cd /home/celphin/projects/rpp-rieseber/celphin/Dryas/RNAseq_analysis/input/RNAscripts/slurm_logs

for x in ../N_*.sh
do
sbatch ${x}
done

#-------------------------------
# submit all Alex

cd /home/celphin/projects/rpp-rieseber/celphin/Dryas/RNAseq_analysis/input/RNAscripts/slurm_logs

for x in ../A_*.sh
do
sbatch ${x}
done

#-------------------------------

sq
          JOBID     USER      ACCOUNT           NAME  ST  TIME_LEFT NODES CPUS TRES_PER_N MIN_MEM NODELIST (REASON)
       60270909  celphin rpp-rieseber       T_C14.sh   R      12:53     1   32        N/A 120000M cdr1548 (None)
       60270910  celphin rpp-rieseber       T_C15.sh   R      12:53     1   32        N/A 120000M cdr1561 (None)
       60270915  celphin rpp-rieseber        T_C7.sh   R      12:53     1   32        N/A 120000M cdr2056 (None)
       60271227  celphin rpp-rieseber       T_C12.sh   R      16:55     1   32        N/A 120000M cdr1605 (None)
       60271452  celphin rpp-rieseber      N_C2_B.sh   R      18:46     1   32        N/A 120000M cdr1539 (None)
       60271457  celphin rpp-rieseber      N_C9_B.sh   R      18:47     1   32        N/A 120000M cdr2050 (None)
       60271460  celphin rpp-rieseber       N_W11.sh   R      18:47     1   32        N/A 120000M cdr1999 (None)
       60271466  celphin rpp-rieseber        N_W8.sh   R      24:29     1   32        N/A 120000M cdr1591 (None)
       60271637  celphin rpp-rieseber       A_W23.sh   R      36:13     1   32        N/A 120000M cdr1673 (None)
       60271639  celphin rpp-rieseber       A_W24.sh   R      37:12     1   32        N/A 120000M cdr1662 (None)
       60271647  celphin rpp-rieseber        A_W7.sh   R      39:16     1   32        N/A 120000M cdr1668 (None)


# some above may need to be resubmitted - since taking longer
cd /home/celphin/projects/rpp-rieseber/celphin/Dryas/RNAseq_analysis/input/RNAscripts/
sbatch T_C14.sh
sbatch T_C15.sh
sbatch T_C7.sh
sbatch T_C12.sh
sbatch N_C2_B.sh
sbatch N_C9_B.sh
sbatch N_W11.sh
sbatch N_W8.sh
sbatch A_W23.sh
sbatch A_W24.sh 
sbatch A_W7.sh

##############################
# Rename seedlings and Lat samples done in Sept 2022
cd /home/celphin/projects/rpp-rieseber/celphin/Dryas/RNAseq_analysis/June2022_output

rename W_L W_Swed_ *
rename C_L C_Swed_ *

rename W_SE W_Seed_ *
rename C_SE C_Seed_ *


# rename new files
rename A_W W_Alex_ *
rename A_C C_Alex_ *

rename N_W W_Norw_ *
rename N_C C_Norw_ *

rename T_W W_Alas_ *
rename T_C C_Alas_ *


#-------------
# one directory down rename too
rename W_L W_Swed_ */*
rename C_L C_Swed_ */*

rename W_SE W_Seed_ */*
rename C_SE C_Seed_ */*

rename N_W W_Norw_ */*
rename N_C C_Norw_ */*

rename T_W W_Alas_ */*
rename T_C C_Alas_ */*

rename A_W W_Alex_ */*
rename A_C C_Alex_ */*

############################################################
# check data

tmux new-session -s RNA
tmux attach-session -t RNA

module load StdEnv/2020
module load bowtie2/2.4.4
module load r/4.2.1
module load samtools/1.15.1

/home/celphin/projects/rpp-rieseber/celphin/Dryas/RNAseq_analysis/RSEM-1.2.31/rsem-plot-model C_Norw_5 C_Norw_5.pdf

# looks good!

################################
# join data in table 

#------------------------------
# find perl script from before
#cd /home/celphin/projects/rpp-rieseber/celphin/Other
#grep add_RSEM_data_to_table.pl biol525d.txt
#tar -zxvf biol525d.tar.gz home/cassandra.elphinstone/biol525d/data/Topic_6/scripts/add_RSEM_data_to_table.pl
#cp add_RSEM_data_to_table.pl /home/celphin/projects/rpp-rieseber/celphin/Dryas/RNAseq_analysis/

# run script on RSEM output
cd /home/celphin/projects/rpp-rieseber/celphin/Dryas/RNAseq_analysis/RNA_site_tables
ls ../output/*genes.results > list_to_add.txt
grep ">" ../reference/Dryas_octopetala_H1.transcript.fa | sed 's/>//' > gene_names.txt
perl ./add_RSEM_data_to_table.pl list_to_add.txt gene_names.txt _expression_table.txt

# edit header line for R
head gene_names_expression_table.txt

# remove ./output/ from header
sed 's/\.\.\/output\///g' gene_names_expression_table.txt > gene_names_expression_table1.txt
head gene_names_expression_table1.txt

#------------------------------------------------------
# run script on RSEM output on SE=seedlings only
cd /home/celphin/projects/rpp-rieseber/celphin/Dryas/RNAseq_analysis/RNA_site_tables
ls ./output/*SE*genes.results > SE_list_to_add.txt
grep ">" ../reference/Dryas_octopetala_H1.transcript.fa | sed 's/>//' > gene_names.txt
perl ./add_RSEM_data_to_table.pl SE_list_to_add.txt gene_names.txt _SE_expression_table.txt

# edit header line for R
more gene_names_SE_expression_table.txt
sequence ./output/C_SE1 ./output/C_SE2 ./output/C_SE3 ./output/W_SE1 ./output/W_SE2 ./output/W_SE3

# remove ./output/ from header
sed 's/\.\.\/output\///g' gene_names_SE_expression_table.txt
less gene_names_SE_expression_table.txt
sequence C_SE1 C_SE2 C_SE3 W_SE1 W_SE2 W_SE3

# needs at least 3 samples per type to check

#------------------------------------------------------
# run script on RSEM output on Seedlings divided by growth chamber

cd /home/celphin/projects/rpp-rieseber/celphin/Dryas/RNAseq_analysis/RNA_site_tables

cp gene_names_SE_expression_table.txt gene_names_SE_L_H_expression_table.txt

# remove ./output/ from header
sed 's/\.\/output\///g' gene_names_SE_L_H_expression_table.txt
nano gene_names_SE_L_H_expression_table.txt

sequence L_C_SE1 L_C_SE2 H_C_SE3 L_W_SE1 H_W_SE2 H_W_SE3


###########################################################