##############################
# SNP calling from WGBS data
# Nov 2024
# https://github.com/hellbelly/BS-Snper 
##############################

# Install

# Hardware requirements
	# One computing node equipped with at least 10 GB Memory
# Software requirements
	# GCC 4.6.0 or higher
	# Perl 5.16.3 or higher
	# zlib 1.2.8 or higher

module load StdEnv/2023 gcc/12.3 perl/5.36.1 r-bundle-bioconductor/3.18 

# Beluga
cd /home/celphin/scratch/Dryas/BS-Snper

git clone https://github.com/hellbelly/BS-Snper.git

module load StdEnv/2020 samtools/0.1.20 gcc/9.3.0 r-bundle-bioconductor/3.17 perl/5.30.2 bcftools

# GCC 4.6.0 or higher
# Perl 5.16.3 or higher
# zlib 1.2.8 or higher
# samtools

cd BS-Snper
sh BS-Snper.sh

# Make sure the executable file rrbsSnp is generated.

/cvmfs/soft.computecanada.ca/gentoo/2020/usr/x86_64-pc-linux-gnu/binutils-bin/2.33.1/ld.gold: internal error in read_header_prolog, at /cvmfs/soft.computecanada.ca/gentoo/2020/var/tmp/portage/sys-devel/binutils-2.33.1-r1/work/binutils-2.33.1/gold/dwarf_reader.cc:1678
collect2: error: ld returned 1 exit status
make[1]: *** [Makefile:47: samtools] Error 1
make[1]: Leaving directory '/lustre04/scratch/celphin/Dryas/BS-Snper/BS-Snper/samtools-0.1.19'
make: *** [Makefile:27: all-recur] Error 1
g++ -O2 -o rrbsSnp main.o sam_funcs.o hash_funcs.o chrome_funcs.o  -m64 -I./samtools-0.1.19/ -L./samtools-0.1.19/ -lbam -lz -lpthread
/cvmfs/soft.computecanada.ca/gentoo/2020/usr/x86_64-pc-linux-gnu/binutils-bin/2.33.1/ld.gold: warning: main.o: unknown program property type 0xc0010002 in .note.gnu.property section

#-----------------------------
# copy over the raw filtered bam data and the reference genome to run this on

# Cedar to beluga
tmux new-session -s BS-Snper
tmux attach-session -t BS-Snper

tar -xvf May2021_Methylseq_output.tar.gz

#---------------------------------
# run on one indvidual

tmux new-session -s BS-Snper1
tmux attach-session -t BS-Snper1

module load StdEnv/2020 samtools/0.1.20 gcc/9.3.0 r-bundle-bioconductor/3.17 perl/5.30.2 bcftools

cd /home/celphin/scratch/Dryas/BS-Snper

# Too many characters in one row! Try to split the long row into several short rows (fewer than 1000000 character s per row).
# Error! at /home/celphin/scratch/Dryas/BS-Snper/BS-Snper/BS-Snper.pl line 110.
seqkit seq DoctH0_Main.fasta -w 60 > DoctH0_Main60.fasta 

NAME=Asian1_F112573

    perl /home/celphin/scratch/Dryas/BS-Snper/BS-Snper/BS-Snper.pl \
    --fa /home/celphin/scratch/Dryas/BS-Snper/Bam_data/DoctH0_Main60.fasta \
    --input /home/celphin/scratch/Dryas/BS-Snper/Bam_data/${NAME}.sorted.bam \
    --output /home/celphin/scratch/Dryas/BS-Snper/out/${NAME}.SNP-candidates.txt \
    --methcg /home/celphin/scratch/Dryas/BS-Snper/out/${NAME}.CpG-meth-info.tab \
    --methchg /home/celphin/scratch/Dryas/BS-Snper/out/${NAME}.CHG-meth-info.tab \
    --methchh /home/celphin/scratch/Dryas/BS-Snper/out/${NAME}.CHH-meth-info.tab \
    --mincover 5 \
    > /home/celphin/scratch/Dryas/BS-Snper/out/${NAME}.SNP-results.vcf 2> /home/celphin/scratch/Dryas/BS-Snper/out/${NAME}.ERR.log


# need to fix samtools Install
tail ${NAME}.ERR.log
# Can't exec "samtools": No such file or directory at /home/celphin/scratch/Dryas/BS-Snper/BS-Snper/BS-Snper.pl line 129.
# No such file or directory at /home/celphin/scratch/Dryas/BS-Snper/BS-Snper/BS-Snper.pl line 129.

$PATH
/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx512/Core/samtools/1.17/bin:/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx512/Compiler/gcc9/bcftools/1.16/bin:/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx512/Core/htslib/1.16/bin:/cvmfs/soft.computecanada.ca/easybuild/software/2020/Core/perl/5.30.2/bin:/cvmfs/soft.computecanada.ca/easybuild/software/2020/Core/expat/2.2.9/bin:/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx512/Core/protobuf/3.21.3/bin:/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx512/Compiler/gcc9/fftw/3.3.8/bin:/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx512/Compiler/gcc9/netcdf/4.7.4/bin:/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx512/Compiler/gcc9/hdf5/1.10.6/bin:/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx512/Core/r/4.3.1/bin:/cvmfs/soft.computecanada.ca/easybuild/software/2020/Core/flexiblas/3.0.4/bin:/cvmfs/soft.computecanada.ca/easybuild/software/2020/Core/flexiblas/3.0.4/bin:/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx512/Compiler/gcc9/openmpi/4.0.3/bin:/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx512/Core/libfabric/1.10.1/bin:/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx512/Core/ucx/1.8.0/bin:/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx512/Compiler/gcc9/gsl/2.6/bin:/cvmfs/soft.computecanada.ca/easybuild/software/2020/Core/gcccore/9.3.0/bin:/cvmfs/soft.computecanada.ca/easybuild/software/2020/Core/java/13.0.2:/cvmfs/soft.computecanada.ca/easybuild/software/2020/Core/java/13.0.2/bin:/cvmfs/soft.computecanada.ca/easybuild/software/2020/Core/rust/1.70.0:/cvmfs/soft.computecanada.ca/easybuild/software/2020/Core/rust/1.70.0/bin:/cvmfs/soft.computecanada.ca/easybuild/software/2020/Core/mii/1.1.2/bin:/cvmfs/soft.computecanada.ca/easybuild/bin:/cvmfs/soft.computecanada.ca/custom/bin:/cvmfs/soft.computecanada.ca/gentoo/2020/usr/sbin:/cvmfs/soft.computecanada.ca/gentoo/2020/usr/bin:/cvmfs/soft.computecanada.ca/gentoo/2020/sbin:/cvmfs/soft.computecanada.ca/gentoo/2020/bin:/cvmfs/soft.computecanada.ca/custom/bin/computecanada:/home/celphin/scratch/Oxyria/DupGen_finder:/opt/software/slurm/bin:/home/celphin/scratch/Oxyria/DupGen_finder:/home/celphin/miniconda2/condabin:/opt/software/slurm/bin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/opt/puppetlabs/bin:/home/celphin/.local/bin:/home/celphin/bin:/home/celphin/.cargo/bin:/home/celphin/.local/bin:/home/

# add samtools path to line 129
nano /home/celphin/scratch/Dryas/BS-Snper/BS-Snper/BS-Snper.pl

open INTV,"/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx512/Core/samtools/1.17/bin/samtools view -H $bam|" or die $!;

# try again
# works



#---------------------------
# fix old ref
module load StdEnv/2020 samtools/0.1.20 gcc/9.3.0 r-bundle-bioconductor/3.17 perl/5.30.2 bcftools samtools/1.17 seqkit/2.3.1

seqkit seq \
/home/celphin/scratch/Dryas/BS-Snper/Bam_data/output/reference_genome/BismarkIndex/Dryas_octopetala_H1.supercontigs.fa \
-w 60 > /home/celphin/scratch/Dryas/BS-Snper/Bam_data/output/reference_genome/BismarkIndex/Dryas_octopetala_H1.supercontigs60.fa

seqkit seq \
/home/celphin/scratch/Dryas/BS-Snper/New_ref_Bams/DoctH0_Main.fasta \
-w 60 > /home/celphin/scratch/Dryas/BS-Snper/New_ref_Bams/DoctH0_Main60.fasta

#---------------------
# test
# https://bioinformatics.stackexchange.com/questions/18538/samtools-sort-most-efficient-memory-and-thread-settings-for-many-samples-on-a-c

tmux new-session -s BS-Snper2
tmux attach-session -t BS-Snper2

salloc -c7 --time 5:00:00 --mem 60000m --account def-cronk

module load StdEnv/2023 minimap2/2.28 samtools/1.20 perl/5.36.1
cd /home/celphin/scratch/Dryas/BS-Snper/Bam_data/output/bismark_deduplicated
samtools cat C1.A10.C1d12_CASS5C_529_159_R1_val_1_bismark_bt2_pe.deduplicated.bam | samtools sort -T C_CASS5C_529_159.sorted --threads 6 -m10g -o C_CASS5C_529_159.sorted.bam

#-------------------------

# loop and submit script for all on old ref
cd /home/celphin/scratch/Dryas/BS-Snper/Bam_data/output/bismark_deduplicated
FILES=$(ls *bismark_bt2_pe.deduplicated.bam)
echo ${FILES}

for file in ${FILES}
do
    NAME=$(echo ${file} | awk -F'_' '{print substr($1,1,1) "_" $2 "_" $3 "_" $4}')
    echo ${NAME}
    echo ${file}

cat << EOF > bamsort_${NAME}.sh
#!/bin/bash
#SBATCH --account=def-cronk
#SBATCH --time=0-5:00:00
#SBATCH --cpus-per-task=7
#SBATCH --mem=60000M

module load StdEnv/2023 minimap2/2.28 samtools/1.20 perl/5.36.1
cd /home/celphin/scratch/Dryas/BS-Snper/Bam_data/output/bismark_deduplicated
samtools cat ${file} | samtools sort -T ${NAME}.sorted --threads 6 -m10g -o ${NAME}.sorted.bam

EOF

sbatch bamsort_${NAME}.sh

done

set -e # Exit if a tool returns a non-zero status/exit code
set -u # Treat unset variables and parameters as an error
set -o pipefail # Returns the status of the last command to exit with a non-zero status or zero if all successfully execute
set -C # No clobber - prevent output redirection from overwriting files.

#---------------------------

# no allocation needed

tmux new-session -s BS-Snper2
tmux attach-session -t BS-Snper2

module load StdEnv/2020 samtools/0.1.20 gcc/9.3.0 r-bundle-bioconductor/3.17 perl/5.30.2 bcftools samtools/1.17

module load StdEnv/2020 seqkit/2.3.1
seqkit seq \
/home/celphin/scratch/Dryas/BS-Snper/Bam_data/output/reference_genome/BismarkIndex/Dryas_octopetala_H1.supercontigs.fa \
-w 60 > /home/celphin/scratch/Dryas/BS-Snper/Bam_data/output/reference_genome/BismarkIndex/Dryas_octopetala_H1.supercontigs60.fa

cd /home/celphin/scratch/Dryas/BS-Snper/Bam_data/output/bismark_deduplicated
FILES=$(ls *bismark_bt2_pe.deduplicated.bam)
echo ${FILES}

for file in ${FILES}
do
    NAME=$(echo ${file} | awk -F'_' '{print substr($1,1,1) "_" $2 "_" $3 "_" $4}')
    echo ${NAME}
	echo ${file}

    perl /home/celphin/scratch/Dryas/BS-Snper/BS-Snper/BS-Snper.pl \
    --fa /home/celphin/scratch/Dryas/BS-Snper/Bam_data/DoctH0_Main60.fasta \
    --input /home/celphin/scratch/Dryas/BS-Snper/Bam_data/${NAME}.sorted.bam \
    --output /home/celphin/scratch/Dryas/BS-Snper/oldref_out/${NAME}.SNP-candidates.txt \
    --methcg /home/celphin/scratch/Dryas/BS-Snper/oldref_out/${NAME}.CpG-meth-info.tab \
    --methchg /home/celphin/scratch/Dryas/BS-Snper/oldref_out/${NAME}.CHG-meth-info.tab \
    --methchh /home/celphin/scratch/Dryas/BS-Snper/oldref_out/${NAME}.CHH-meth-info.tab \
    --mincover 5 \
    > /home/celphin/scratch/Dryas/BS-Snper/oldref_out/${NAME}.SNP-results.vcf 2> /home/celphin/scratch/Dryas/BS-Snper/oldref_out/${NAME}.ERR.log

done

#-------------------------

# index and merge

for file in ${FILES}
do
    NAME=$(echo ${file} | awk -F'_' '{print substr($1,1,1) "_" $2 "_" $3 "_" $4}')
    echo ${NAME}
    echo ${file}

module load StdEnv/2023 vcftools/0.1.16 
bgzip /home/celphin/scratch/Dryas/BS-Snper/oldref_out/${NAME}.SNP-results.vcf
tabix -p vcf /home/celphin/scratch/Dryas/BS-Snper/oldref_out/${NAME}.SNP-results.vcf.gz

done


# Merge vcf files
# https://vcftools.github.io/perl_module.html#vcf-merge

module load StdEnv/2023 vcftools/0.1.16 
vcf-merge *.SNP-results.vcf.gz  > Total_Dryas_Oldref.vcf




####################################
# Also run  merged bam file and individual bam files for new reference genome
# https://www.htslib.org/doc/samtools-merge.html

# control
tmux new-session -s BS-Snper
tmux attach-session -t BS-Snper

#salloc -c20 --time 7:00:00 --mem 120000m --account def-cronk

cd /home/celphin/scratch/Dryas/BS-Snper/New_ref_Bams
module load StdEnv/2023 gcc/12.3 samtools/1.20
samtools merge C_merged_deduplicated_sorted.bam C_*

#-----------
# Warming
tmux new-session -s BS-Snper1
tmux attach-session -t BS-Snper1

#salloc -c20 --time 7:00:00 --mem 120000m --account def-cronk

cd /home/celphin/scratch/Dryas/BS-Snper/New_ref_Bams
module load StdEnv/2023 gcc/12.3 samtools/1.20
samtools merge W_merged_deduplicated_sorted.bam W_*.deduplicated.sorted.bam

samtools merge Wild_merged_deduplicated_sorted.bam W_merged_deduplicated_sorted.bam C_merged_deduplicated_sorted.bam

samtools view -H Wild_merged_deduplicated_sorted.bam

cd /home/celphin/scratch/Dryas/BS-Snper/BS-Snper/

#Defaults: --minhetfreq 0.1 --minhomfreq 0.85 --minquali 15 --maxcover 1000 --minread2 2 --errorate 0.02 --mapvalue 20
perl /home/celphin/scratch/Dryas/BS-Snper/BS-Snper/BS-Snper.pl \
    --fa /home/celphin/scratch/Dryas/BS-Snper/New_ref_Bams/DoctH0_Main60.fasta \
    --input /home/celphin/scratch/Dryas/BS-Snper/New_ref_Bams/Wild_merged_deduplicated_sorted.bam  \
    --output /home/celphin/scratch/Dryas/BS-Snper/out/NewRef_Wild.SNP-candidates.txt \
    --methcg /home/celphin/scratch/Dryas/BS-Snper/out/NewRef_Wild.CpG-meth-info.tab \
    --methchg /home/celphin/scratch/Dryas/BS-Snper/out/NewRef_Wild.CHG-meth-info.tab \
    --methchh /home/celphin/scratch/Dryas/BS-Snper/out/NewRef_Wild.CHH-meth-info.tab \
    --mincover 5 \
    > /home/celphin/scratch/Dryas/BS-Snper/out/NewRef_Wild.SNP-results.vcf \
	2> /home/celphin/scratch/Dryas/BS-Snper/out/NewRef_Wild.ERR.log

# names
samtools view -H Wild_merged_deduplicated_sorted.bam

C_ALAS_00C_227.deduplicated.sorted.bam C_ALAS_00C_231.deduplicated.sorted.bam C_ALAS0C_10_246.deduplicated.sorted.bam C_ALAS0C_12_256.deduplicated.sorted.bam C_ALAS0C_13_254.deduplicated.sorted.bam C_ALAS0C_18_229.deduplicated.sorted.bam C_ALAS0C_19_261.deduplicated.sorted.bam C_ALAS0C_3_258.deduplicated.sorted.bam C_ALAS0C_4_240.deduplicated.sorted.bam C_ALAS0C_5_238.deduplicated.sorted.bam C_CASS_09C_541.deduplicated.sorted.bam C_CASS10C_548_144.deduplicated.sorted.bam C_CASS17C_576_175.deduplicated.sorted.bam C_CASS4C_524_4.deduplicated.sorted.bam C_CASS5C_529_159.deduplicated.sorted.bam C_CASS8C_535_54.deduplicated.sorted.bam C_DRY10C_60_41.deduplicated.sorted.bam C_DRY2C_10_52.deduplicated.sorted.bam C_DRY4C_23_82.deduplicated.sorted.bam C_DRY5C_28_92.deduplicated.sorted.bam C_DRY9C_53_149.deduplicated.sorted.bam C_FERT13C_7F_112.deduplicated.sorted.bam C_FERT31C_15F_170.deduplicated.sorted.bam C_FERT39C_20F_71.deduplicated.sorted.bam C_FERT5C_1F_97.deduplicated.sorted.bam C_LATD1C_4_223.deduplicated.sorted.bam C_LATD2C_1_203.deduplicated.sorted.bam C_LATD2C_6_198.deduplicated.sorted.bam C_LATD2C_7_209.deduplicated.sorted.bam C_LATD4C_3_196.deduplicated.sorted.bam C_LATD5C_20_199.deduplicated.sorted.bam C_LATD5C_2_201.deduplicated.sorted.bam C_LATD5C_5_191.deduplicated.sorted.bam C_LATJ_00C_187.deduplicated.sorted.bam C_LATJ_02C_194.deduplicated.sorted.bam C_MEAD_03C_1R0.deduplicated.sorted.bam C_MEAD1C_446_33.deduplicated.sorted.bam C_MEAD2C_451_76.deduplicated.sorted.bam C_MEAD6C_468_22.deduplicated.sorted.bam C_MEAD7C_473_95.deduplicated.sorted.bam C_SVAL_0C_268.deduplicated.sorted.bam C_SVAL_0C_269.deduplicated.sorted.bam C_SVAL12C_12_273.deduplicated.sorted.bam C_SVAL16C_16_276.deduplicated.sorted.bam C_SVAL49C_49_278.deduplicated.sorted.bam C_SVAL8C_8_275.deduplicated.sorted.bam C_WILL10C_440_16.deduplicated.sorted.bam C_WILL1C_406_152.deduplicated.sorted.bam C_WILL3C_414_100.deduplicated.sorted.bam C_WILL5C_422_31.deduplicated.sorted.bam C_WILL7C_445_125.deduplicated.sorted.bam
W_ALAS_00W_228.deduplicated.sorted.bam W_ALAS_00W_232.deduplicated.sorted.bam W_ALAS0W_14_249.deduplicated.sorted.bam W_ALAS0W_15_242.deduplicated.sorted.bam W_ALAS0W_16_239.deduplicated.sorted.bam W_ALAS0W_17_235.deduplicated.sorted.bam W_ALAS0W_18_248.deduplicated.sorted.bam W_ALAS0W_3_236.deduplicated.sorted.bam W_ALAS0W_7_263.deduplicated.sorted.bam W_ALAS0W_8_265.deduplicated.sorted.bam W_CASS_04W_519.deduplicated.sorted.bam W_CASS10W_544_60.deduplicated.sorted.bam W_CASS17W_574_137.deduplicated.sorted.bam W_CASS5W_525_130.deduplicated.sorted.bam W_CASS7W_600_19.deduplicated.sorted.bam W_CASS9W_539_128.deduplicated.sorted.bam W_DRY1W_3_39.deduplicated.sorted.bam W_DRY3W_15_69.deduplicated.sorted.bam W_DRY6W_31_147.deduplicated.sorted.bam W_DRY8W_45_155.deduplicated.sorted.bam W_DRY9W_50_185.deduplicated.sorted.bam W_FERT14W_6F_126.deduplicated.sorted.bam W_FERT22W_12F_111.deduplicated.sorted.bam W_FERT30W_14F_40.deduplicated.sorted.bam W_FERT6W_3F_110.deduplicated.sorted.bam W_LATC1W_12_219.deduplicated.sorted.bam W_LATC3W_16_220.deduplicated.sorted.bam W_LATC5W_18_190.deduplicated.sorted.bam W_LATC9W_11_216.deduplicated.sorted.bam W_LATD2W_4_212.deduplicated.sorted.bam W_LATD2W_5_206.deduplicated.sorted.bam W_LATD4W_8_211.deduplicated.sorted.bam W_LATD4W_9_207.deduplicated.sorted.bam W_LATJ_02W_193.deduplicated.sorted.bam W_LATJ_04W_188.deduplicated.sorted.bam W_MEAD_08W_4R0.deduplicated.sorted.bam W_MEAD1W_444_116.deduplicated.sorted.bam W_MEAD2W_450_70.deduplicated.sorted.bam W_MEAD6W_466_163.deduplicated.sorted.bam W_MEAD7W_470_173.deduplicated.sorted.bam W_SVAL_0W_267.deduplicated.sorted.bam W_SVAL_0W_270.deduplicated.sorted.bam W_SVAL16W_16_277.deduplicated.sorted.bam W_SVAL18W_18_271.deduplicated.sorted.bam W_SVAL6W_6_274.deduplicated.sorted.bam W_SVAL8W_8_272.deduplicated.sorted.bam W_WILL1W_403_67.deduplicated.sorted.bam W_WILL4W_417_13.deduplicated.sorted.bam W_WILL5W_421_154.deduplicated.sorted.bam W_WILL7W_448_107.deduplicated.sorted.bam
> Wild_list.txt

#---------------
# Seedling/Parent
tmux new-session -s BS-Snper2
tmux attach-session -t BS-Snper2

#salloc -c20 --time 7:00:00 --mem 120000m --account def-cronk

cd /home/celphin/scratch/Dryas/BS-Snper/New_ref_Bams
module load StdEnv/2023 gcc/12.3 samtools/1.20
samtools merge SE_merged_deduplicated_sorted.bam SE_*.deduplicated.sorted.bam \
C_LATD2C_6_198.deduplicated.sorted.bam \
W_LATC1W_12_219.deduplicated.sorted.bam \
C_LATJ_00C_187.deduplicated.sorted.bam \
C_LATJ_02C_194.deduplicated.sorted.bam \
W_LATC1W_12_219.deduplicated.sorted.bam \
C_LATD5C_5_191.deduplicated.sorted.bam \
W_ALAS_00W_228.deduplicated.sorted.bam

perl /home/celphin/scratch/Dryas/BS-Snper/BS-Snper/BS-Snper.pl \
    --fa /home/celphin/scratch/Dryas/BS-Snper/New_ref_Bams/DoctH0_Main60.fasta \
    --input /home/celphin/scratch/Dryas/BS-Snper/New_ref_Bams/SE_merged_deduplicated_sorted.bam  \
    --output /home/celphin/scratch/Dryas/BS-Snper/out/NewRef_SE.SNP-candidates.txt \
    --methcg /home/celphin/scratch/Dryas/BS-Snper/out/NewRef_SE.CpG-meth-info.tab \
    --methchg /home/celphin/scratch/Dryas/BS-Snper/out/NewRef_SE.CHG-meth-info.tab \
    --methchh /home/celphin/scratch/Dryas/BS-Snper/out/NewRef_SE.CHH-meth-info.tab \
    --mincover 5 \
    > /home/celphin/scratch/Dryas/BS-Snper/out/NewRef_SE.SNP-results.vcf \
	2> /home/celphin/scratch/Dryas/BS-Snper/out/NewRef_SE.ERR.log

ls SE_*.deduplicated.sorted.bam \
C_LATD2C_6_198.deduplicated.sorted.bam \
W_LATC1W_12_219.deduplicated.sorted.bam \
C_LATJ_00C_187.deduplicated.sorted.bam \
C_LATJ_02C_194.deduplicated.sorted.bam \
W_LATC1W_12_219.deduplicated.sorted.bam \
C_LATD5C_5_191.deduplicated.sorted.bam \
W_ALAS_00W_228.deduplicated.sorted.bam > SE_list.txt

# names
samtools view -H SE_merged_deduplicated_sorted.bam

# SE_A_C_L_125_11_F112561 SE_A_C_L_170_5_F112555 SE_L_C_H_187_98_F112565 SE_L_C_H_191_100_F112569 SE_L_C_H_191_102_F112568 SE_L_C_H_194_92_F112566 SE_L_C_H_194_92J_F112572 SE_L_C_L_187_38_F112560 SE_L_C_L_198_35_F112558 SE_L_W_H_219_104_F112567 SE_L_W_H_219_105_F112564 SE_L_W_H_219_105J_F112571 SE_L_W_L_190_43_F112556 SE_L_W_L_198_34_F112563 SE_L_W_L_219_42_F112557 SE_T_C_L_238_47_F112559 SE_T_W_H_228_95_F112570 C_LATD2C_6_198 W_LATC1W_12_219 C_LATJ_00C_187 C_LATJ_02C_194 W_LATC1W_12_219 C_LATD5C_5_191 W_ALAS_00W_228

#-----------------
# Phenology
tmux new-session -s BS-Snper3
tmux attach-session -t BS-Snper3

#salloc -c20 --time 7:00:00 --mem 120000m --account def-cronk

cd /home/celphin/scratch/Dryas/BS-Snper/New_ref_Bams
module load StdEnv/2023 gcc/12.3 samtools/1.20
samtools merge MatFl_merged_deduplicated_sorted.bam MatFl_*.deduplicated.sorted.bam \
C_CASS4C_524_4.deduplicated.sorted.bam \
C_FERT5C_1F_97.deduplicated.sorted.bam \
C_MEAD1C_446_33.deduplicated.sorted.bam \
C_WILL3C_414_100.deduplicated.sorted.bam \
W_CASS10W_544_60.deduplicated.sorted.bam \
W_CASS5W_525_130.deduplicated.sorted.bam \
W_FERT6W_3F_110.deduplicated.sorted.bam \
W_MEAD1W_444_116.deduplicated.sorted.bam \
W_WILL4W_417_13.deduplicated.sorted.bam

perl /home/celphin/scratch/Dryas/BS-Snper/BS-Snper/BS-Snper.pl \
    --fa /home/celphin/scratch/Dryas/BS-Snper/New_ref_Bams/DoctH0_Main60.fasta \
    --input /home/celphin/scratch/Dryas/BS-Snper/New_ref_Bams/MatFl_merged_deduplicated_sorted.bam  \
    --output /home/celphin/scratch/Dryas/BS-Snper/out/NewRef_MatFl.SNP-candidates.txt \
    --methcg /home/celphin/scratch/Dryas/BS-Snper/out/NewRef_MatFl.CpG-meth-info.tab \
    --methchg /home/celphin/scratch/Dryas/BS-Snper/out/NewRef_MatFl.CHG-meth-info.tab \
    --methchh /home/celphin/scratch/Dryas/BS-Snper/out/NewRef_MatFl.CHH-meth-info.tab \
    --mincover 5 \
    > /home/celphin/scratch/Dryas/BS-Snper/out/NewRef_MatFl.SNP-results.vcf \
	2> /home/celphin/scratch/Dryas/BS-Snper/out/NewRef_MatFl.ERR.log

# names
samtools view -H MatFl_merged_deduplicated_sorted.bam

# MatFl_Cass_10W_60_544_F112581 MatFl_Cass_4C_4_524_F112579 MatFl_Cass_5W_130_525_F112580 MatFl_Fert_5C_97_1F_F112583 MatFl_Fert_6W_110_3F_F112584 MatFl_Mead_1C_33_446_F112577 MatFl_Mead_1W_116_444_F112578 MatFl_Will_3C_100_414_F112575 MatFl_Will_4W_13_417_F112576 C_CASS4C_524_4 C_FERT5C_1F_97 C_MEAD1C_446_33 C_WILL3C_414_100 W_CASS10W_544_60 W_CASS5W_525_130 W_FERT6W_3F_110 W_MEAD1W_444_116 W_WILL4W_417_13

#-----------------
tmux new-session -s BS-Snper4
tmux attach-session -t BS-Snper4

#salloc -c20 --time 7:00:00 --mem 120000m --account def-cronk

cd /home/celphin/scratch/Dryas/BS-Snper/New_ref_Bams
module load StdEnv/2023 gcc/12.3 samtools/1.20
samtools merge Other_merged_deduplicated_sorted.bam \
Asian1_F112573.deduplicated.sorted.bam \
L_sample5.deduplicated.sorted.bam \
Chilliwack1_F112574.deduplicated.sorted.bam \
C_ALAS_00C_227.deduplicated.sorted.bam \
C_ALAS0C_10_246.deduplicated.sorted.bam \
C_LATD1C_4_223.deduplicated.sorted.bam \
C_LATD5C_20_199.deduplicated.sorted.bam \
C_SVAL_0C_268.deduplicated.sorted.bam \
C_SVAL16C_16_276.deduplicated.sorted.bam \
C_CASS4C_524_4.deduplicated.sorted.bam \
C_FERT5C_1F_97.deduplicated.sorted.bam \
C_MEAD1C_446_33.deduplicated.sorted.bam \
C_DRY10C_60_41.deduplicated.sorted.bam \
C_WILL3C_414_100.deduplicated.sorted.bam 

perl /home/celphin/scratch/Dryas/BS-Snper/BS-Snper/BS-Snper.pl \
    --fa /home/celphin/scratch/Dryas/BS-Snper/New_ref_Bams/DoctH0_Main60.fasta \
    --input /home/celphin/scratch/Dryas/BS-Snper/New_ref_Bams/Other_merged_deduplicated_sorted.bam  \
    --output /home/celphin/scratch/Dryas/BS-Snper/out/NewRef_Other.SNP-candidates.txt \
    --methcg /home/celphin/scratch/Dryas/BS-Snper/out/NewRef_Other.CpG-meth-info.tab \
    --methchg /home/celphin/scratch/Dryas/BS-Snper/out/NewRef_Other.CHG-meth-info.tab \
    --methchh /home/celphin/scratch/Dryas/BS-Snper/out/NewRef_Other.CHH-meth-info.tab \
    --mincover 5 \
    > /home/celphin/scratch/Dryas/BS-Snper/out/NewRef_Other.SNP-results.vcf \
	2> /home/celphin/scratch/Dryas/BS-Snper/out/NewRef_Other.ERR.log

# names
samtools view -H Other_merged_deduplicated_sorted.bam

# Asian1_F112573 L_sample5 Chilliwack1_F112574 C_ALAS_00C_227 C_ALAS0C_10_246 C_LATD1C_4_223 C_LATD5C_20_199 C_SVAL_0C_268 C_SVAL16C_16_276 C_CASS4C_524_4 C_FERT5C_1F_97 C_MEAD1C_446_33 C_DRY10C_60_41 C_WILL3C_414_100

##################################
# Run for all samples
tmux new-session -s BS-Snper
tmux attach-session -t BS-Snper

cd /home/celphin/scratch/Dryas/BS-Snper/New_ref_Bams/
FILES=$(ls *.deduplicated.sorted.bam)
echo ${FILES}

cd /home/celphin/scratch/Dryas/BS-Snper/scripts
for file in ${FILES}
do
    NAME=$(echo ${file} | awk -F'.' '{print $1}')
    echo ${NAME}
    echo ${file}

cat << EOF > BS-SNPer_${NAME}.sh
#!/bin/bash
#SBATCH --account=def-cronk
#SBATCH --time=0-1:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8000M

    perl /home/celphin/scratch/Dryas/BS-Snper/BS-Snper/BS-Snper.pl \
    --fa /home/celphin/scratch/Dryas/BS-Snper/New_ref_Bams/DoctH0_Main60.fasta \
    --input //home/celphin/scratch/Dryas/BS-Snper/New_ref_Bams/${file} \
    --output /home/celphin/scratch/Dryas/BS-Snper/out/${NAME}.SNP-candidates.txt \
    --methcg /home/celphin/scratch/Dryas/BS-Snper/out/${NAME}.CpG-meth-info.tab \
    --methchg /home/celphin/scratch/Dryas/BS-Snper/out/${NAME}.CHG-meth-info.tab \
    --methchh /home/celphin/scratch/Dryas/BS-Snper/out/${NAME}.CHH-meth-info.tab \
    --mincover 5 \
    > /home/celphin/scratch/Dryas/BS-Snper/out/${NAME}.Indiv_SNP-results.vcf 2> /home/celphin/scratch/Dryas/BS-Snper/out/${NAME}.ERR.log

EOF

sbatch BS-SNPer_${NAME}.sh

done


for file in ${FILES}
do
    NAME=$(echo ${file} | awk -F'.' '{print $1}')
    echo ${NAME}
    echo ${file}

module load StdEnv/2023 vcftools/0.1.16 
bgzip /home/celphin/scratch/Dryas/BS-Snper/out/${NAME}.Indiv_SNP-results.vcf
tabix -p vcf /home/celphin/scratch/Dryas/BS-Snper/out/${NAME}.Indiv_SNP-results.vcf.gz

done

#------------------------
# Merge vcf files
# https://vcftools.github.io/perl_module.html#vcf-merge

module load StdEnv/2023 vcftools/0.1.16 
vcf-merge *.Indiv_SNP-results.vcf.gz  > Total_Dryas_Newref.vcf


##############################
# Basic pop gen 

# SNP filtering (vcftools)
## quality, minor allele counts, % found accross individuals
## filter out individuals with lots of missing data
## A script to count the number of potential genotyping errors due to low read depth

# Plink
## make a plink file from vcf
## make .bed file from Plink

# Admixture
## run admixture for different k values with loop
## bootstrap admixture results

###########################
# General SNP filtering starts

mkdir pop_gen
cp ./out/Total_Dryas_Newref.vcf ./pop_gen/
cd /home/celphin/scratch/Dryas/BS-Snper/pop_gen

#load vcftools
module load StdEnv/2020 
module load vcftools/0.1.16
module load gnuplot/5.4.2
module load gcc/9.3.0
module load bcftools/1.16
module load plink/1.9b_6.21-x86_64

#---------------------------
# rename samples in vcf file
#https://www.biostars.org/p/279195/ 

salloc -c1 --time 2:00:00 --mem 70000m --account def-henryg

sed 's/\/\/home\/celphin\/scratch\/Dryas\/BS\-Snper\/New_ref_Bams\// /g' Total_Dryas_Newref.vcf > Total_Dryas_Newref_rename.vcf
sed 's/\.deduplicated\.sorted\.bam//g' Total_Dryas_Newref_rename.vcf > Total_Dryas_Newref_rename1.vcf

#-------------
# # not used 
# bcftools query -l Total_Dryas_Newref.vcf > sample_names.txt

# sed 's/\/\/home\/celphin\/scratch\/Dryas\/BS\-Snper\/New_ref_Bams\// /g' ample_names.txt > samples_renamed0.txt

# sed 's/\.deduplicated\.sorted\.bam//g' samples_renamed0.txt > samples_renamed.txt

# bcftools reheader -s samples_renamed.txt -o Total_Dryas_Newref_rename.vcf Total_Dryas_Newref.vcf

# bcftools query -l Total_Dryas_Newref_rename.vcf


#-----------------------------
# all individuals
# filter for quality, indels, biallelic, missing in less than 90%, monomorphic, LD, and minor allele count of 10
vcftools --vcf Total_Dryas_Newref_rename1.vcf \
--minGQ 30 \
--remove-indels \
--max-alleles 2 \
--max-missing 0.90 \
--min-alleles 2 \
--thin 1000 \
--mac 4 \
--recode --recode-INFO-all --out Dryas_filtered
# After filtering, kept 130 out of 130 Individuals
# After filtering, kept 134 905 out of a possible 14860573 Sites

vcftools --vcf Total_Dryas_Newref_rename1.vcf \
--minGQ 30 \
--remove-indels \
--max-alleles 2 \
--max-missing 0.90 \
--min-alleles 2 \
--thin 1000 \
--mac 10 \
--recode --recode-INFO-all --out Dryas_filtered1
# After filtering, kept 130 out of 130 Individuals
# After filtering, kept 69002 out of a possible 14860573 Sites


# Warning: Expected at least 2 parts in INFO entry: ID=DP4,Number=4,Type=Integer,Description="Number of high-quality ref-forward , ref-reverse, alt-forward and alt-reverse bases">
# Warning: Expected at least 2 parts in INFO entry: ID=DP4,Number=4,Type=Integer,Description="Number of high-quality ref-forward , ref-reverse, alt-forward and alt-reverse bases">
# Warning: Expected at least 2 parts in FORMAT entry: ID=BSD,Number=8,Type=Integer,Description="Depth, ATCG in watson strand and crick strand">
# Warning: Expected at least 2 parts in FORMAT entry: ID=BSQ,Number=8,Type=Integer,Description="Avarage Base Quality, ATCG in watson strand and crick strand">
# Warning: Expected at least 2 parts in INFO entry: ID=SF,Number=.,Type=String,Description="Source File (index to sourceFiles, f when filtered)">

################################

#list the amount of missing data per individual - find indivdiuals with no reads mapped
vcftools --vcf Dryas_filtered1.recode.vcf --missing-indv --out Dryas_filtered

#filter out the individuals with greater than 40% missing SNPs 
#mawk '$5 > 0.40' Dryas_filtered.imiss | cut -f1 > lowDP.indv
#vcftools --vcf Dryas_filtered.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out Dryas_filtered_rm20

# get list of all individuals left
#vcftools --vcf Dryas_filtered_rm20.recode.vcf --missing-indv --out Dryas_filtered_rm20

#############################
# Filter vcf to include only some individuals

# grep -E '2006|2007|2008|2009|2010' Dryas_filtered_rm20.imiss | cut -f1  > baseline_indv.txt

# vcftools --vcf Dryas_filtered_rm20.recode.vcf --keep baseline_indv.txt --recode --recode-INFO-all --out Dryas_filtered_baseline
# # After filtering, kept 285 out of 576 Individuals

# vcftools --vcf Dryas_filtered_baseline.recode.vcf --missing-indv --out Dryas_filtered_baseline

################################
# Dryas RUN ADMIXTURE

tmux new-session -s BS-Snper
tmux attach-session -t BS-Snper

# get allocation
salloc -c40 --time 02:50:00 --mem 187G --account def-cronk

# make Plink file to be able to run Admixture
mkdir Admixture_Dryas_Baseline
cp Dryas_filtered1.recode.vcf Admixture_Dryas_Baseline
cd Admixture_Dryas_Baseline

sed 's/DoctH0-//g' Dryas_filtered1.recode.vcf > Dryas_filtered_modified.vcf

# convert the vcf file to Plink
vcftools --vcf Dryas_filtered_modified.vcf --plink --out Dryas_filtered

#MAKE A BED FILE 
plink --file Dryas_filtered --allow-no-sex --chr-set 27 --make-bed --out Dryas_filtered

# Unsupervised analysis K from 1 to 9 and 5 replicates

module load StdEnv/2020
module load admixture/1.3.0

# run loop through various K values
mkdir Take1

cd Take1
cp ../Dryas_filtered* .

for K in 1 2 3 4 5 6 7 8 9 10; \
do admixture --cv=10 -s time -j48 -C 0.0000000001  Dryas_filtered.bed $K | tee log${K}.out; done

#to get the CV errors and see which K value is the best model
cd /home/celphin/scratch/Dryas/BS-Snper/pop_gen/Admixture_Dryas_Baseline/
cd /home/celphin/scratch/Dryas/BS-Snper/pop_gen/Admixture_Dryas_Baseline/Take1
cd /home/celphin/scratch/Dryas/BS-Snper/pop_gen/Admixture_Dryas_Baseline/Take1
grep -h CV ./log*out


CV error (K=1): 0.24980
CV error (K=2): 0.14127
CV error (K=3): 0.13806
CV error (K=4): 0.13754
CV error (K=5): 0.13592
CV error (K=6): 0.13670
CV error (K=7): 0.13612
CV error (K=8): 0.13719



#Dryas low filter
CV error (K=1): 0.25255
CV error (K=2): 0.18321
CV error (K=3): 0.17840
CV error (K=4): 0.18205
CV error (K=5): 0.18009
CV error (K=6): 0.18797


# Mimulus
CV error (K=1): 0.17135
CV error (K=2): 0.13176
CV error (K=3): 0.12380
CV error (K=4): 0.12098
CV error (K=5): 0.12077
CV error (K=6): 0.12008
CV error (K=7): 0.11967
CV error (K=8): 0.11962
CV error (K=9): 0.12030
CV error (K=10): 0.12241


######################################
# Make PCA GDS file

tmux new-session -s BS-Snper1
tmux attach-session -t BS-Snper1

cd /home/celphin/scratch/Dryas/BS-Snper/pop_gen
mkdir PCA


module load r/4.0.2
module load gdal/3.0.1
module load udunits/2.2.26
module load python/3.8.2
export R_LIBS_USER=/home/celphin/R/x86_64-pc-linux-gnu-library/4.0/

R

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("SNPRelate")

library(SNPRelate)

snpgdsVCF2GDS("./Dryas_filtered1.recode.vcf",
              "./Dryas_filtered1_total.recode.gds",
              method="biallelic.only")

# Start file conversion from VCF to SNP GDS ...
# Method: exacting biallelic SNPs
# Number of samples: 130
# Parsing "./Dryas_filtered.recode.vcf" ...
        # import 134905 variants.
# + genotype   { Bit2 130x134905, 4.2M } *
# Optimize the access efficiency ...
# Clean up the fragments of GDS file:
    # open the file './Dryas_filtered_total.recode.gds' (4.9M)
    # # of fragments: 63
    # save to './Dryas_filtered_total.recode.gds.tmp'
    # rename './Dryas_filtered_total.recode.gds.tmp' (4.9M, reduced: 516B)
    # # of fragments: 20


############################
# Plotting

# PCA analysis

# Join all the data into population file and sample file

# Admixture
# - Determine Admixture groups order
# - Set groups colours
# - Plot PCA images with colours
# - Plot map of sites with dominant colours
# - Make Admixture barplot (total, Ancient, Ellemere)
# - Make Admixture maps (samples and population scales)
# - Kluane correlation and barplot

# PopStats
# - Icetime, Diversity and Temperature correlations
# - Plots of Pi and Tajima's D per population
# - Maps of all population avgeraged stats

# Table of PopStats for paper

# Site specific Admixture plots

#####################################
# download files from Cedar to computer

# Admixture
# .4.Q

# # PopStats
# Fst_mean.txt
# Fst_weighted.txt
# Fst_count.txt
# Fst_sites.txt
# population_TajimaD.txt
# summed_site_Pi.txt
# summed_wind_Pi.txt
# Het_data_by_pop.txt

# # PCA GDS file
# Cassiope_noMER_r10i.recode.gds


######################################

# # install packages
# install.packages("devtools") 
# devtools::install_github("MikkoVihtakari/PlotSvalbard", upgrade = "never")
# devtools::install_github("celphin/Equitable")
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("SNPRelate")
# 
# install.packages(c("tidyverse","dplyr","vctrs", "tidyr","reshape","geosphere","ade4","ape","scatterpie","data.table","maps","readr","stringr","ggplot2"))

#----------------------------------
# libraries

library(Equitable)
library(PlotSvalbard)
library(tidyverse)
library(dplyr)
library(tidyr)
library(reshape)
library(geosphere)
library(ade4)
library(ape)
library(scatterpie)
library(data.table)
library(maps)
library(readr) 
library(stringr)
library(SNPRelate)
library(ggplot2)

##############################
# Load in the data

setwd("~/GitHub/Dryas_GenomeBC/9_SNP_calling")

#Load the gds file - for PCA
genofile <- snpgdsOpen("./Figures_data/Mimulus_filtered_baseline.recode.gds")

# load list of samples in filtered vcf
samples_list <- read.table("./Figures_data/Mimulus_filtered_baseline.imiss", header = TRUE)

# load in Admixture data - code setup for 5 populations
Admix_tbl=read.table("./Figures_data/Admixture/Mimulus_filtered_baseline.4.Q")

# load detailed information about all 371 individuals
sample_information  <- read.csv("./Figures_data/M_caridnalis_genomics_meta_data.csv", header = TRUE)

# Population lat long
pop_lat_long <- read.csv("./Figures_data/56_pop_lat.csv", header = TRUE)

# PopStats
#Het_data <- read.table("./Figures_data/PopStats/Het_data_by_pop.txt", header = TRUE)
#Pi_wind_data <- read.table("./Figures_data/PopStats/summed_wind_Pi.txt", header=TRUE, sep=" ")
#Pi_site_data <- read.table("./Figures_data/PopStats/summed_site_Pi.txt", header=TRUE, sep=" ")
#Tajima_data <- read.table("./Figures_data/PopStats/population_TajimaD.txt", header=TRUE, sep=" ")


#################################
# PCA
# https://owensgl.github.io/biol525D/Topic_8-9/pca.html

# create GDS file on server
# snpgdsVCF2GDS("/scratch/celphin/GBS_Cassiope/Mar2020_SNPFiltering_PopStats/R_plots/FiltXXg9mac5minq30r60i_chrom_rLD.vcf.gz",
# "/scratch/celphin/GBS_Cassiope/Mar2020_SNPFiltering_PopStats/R_plots/FiltXXg9mac5minq30r60i_chrom_rLD.gds",
# method="biallelic.only")

#Prune for linkage
snpset_pruned <- snpgdsLDpruning(genofile, autosome.only=F)

#Make a list of sites we're keeping.
snpset.id <- unlist(snpset_pruned)

#Run the PCA
pca0 <- snpgdsPCA(genofile, num.thread = 1, eigen.cnt = 16, snp.id = snpset.id, missing.rate = 0.1, autosome.only = F)
pca5 <- snpgdsPCA(genofile, num.thread = 1, eigen.cnt = 16, snp.id = snpset.id, missing.rate = 0.1, maf=0.05,  autosome.only = F)

###################################
pca <- pca5

#Here's the percent variance explained for each eigenvector
pc.percent <- pca$varprop*100
round(pc.percent, 2)
# 19.99 18.55 15.99 13.10 11.25  6.45  4.42  2.81  2.58  1.11  1.05  0.84  0.67  0.42  0.39  0.26 

#Make a dataframe of your PCA results
PCA_tab <- data.frame(sample = pca$sample.id,
                      PC1 = pca$eigenvect[,1],    # the first eigenvector
                      PC2 = pca$eigenvect[,2],    # the second eigenvector
                      PC3 = pca$eigenvect[,3],    # the first eigenvector
                      PC4 = pca$eigenvect[,4],    # the second eigenvector
                      PC5 = pca$eigenvect[,5],    # the first eigenvector
                      PC6 = pca$eigenvect[,6],    # the second eigenvector
                      stringsAsFactors = FALSE)

str_split(as.character(PCA_tab$sample), "_")

q <- as.data.frame(t(as.matrix(as.data.frame((strsplit(as.character(PCA_tab$sample), "_"))))))
PCA_tab_data <- cbind(PCA_tab, q[,1])

colnames(PCA_tab_data) <- c("ID_code", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6","Pop")

#####################################
# Join all the data

# Sample specific information

# samples_list 
# Admix_tbl
# sample_information 
# Het_data
# PCA_tab_data

# Join Samples information
samples_list$ID_code <- as.factor(samples_list$INDV)
ID_code = samples_list$ID_code

Admix_tbl1 <- cbind(ID_code, Admix_tbl)
All_samples_data <- left_join(Admix_tbl1, sample_information, by="ID_code")
#All_samples_data <- left_join(All_samples_data, Het_data, by="ID_code")
All_samples_data <- left_join(All_samples_data, samples_list, by="ID_code")
All_samples_data <- left_join(All_samples_data, PCA_tab_data, by="ID_code")

#------------------------------
# join all Population information

# population_information  
# Pi_wind_data 
# Pi_site_data 
# Tajima_data 
# 
# colnames(Pi_site_data) <- c("Pop" ,   "SumSitePi" , "Site_SD_Pi",  "Site_var_Pi")
# colnames(Pi_wind_data) <- c("Pop" ,   "SumWindPi",  "Wind_SD_Pi",  "Wind_var_Pi")
# colnames(Tajima_data ) <- c("Pop" ,   "TajimaAvg",  "TajimaSD",  "TajimaVar")
# population_information$Pop[which(population_information$Pop=="25")] <- "025"
# 
# All_pop_data <- left_join(population_information, Tajima_data, by="Pop")
# All_pop_data <- left_join(All_pop_data, Pi_site_data, by="Pop")
# All_pop_data <- left_join(All_pop_data, Pi_wind_data, by="Pop")

#------------------------------
# average the sample information by population
Pop_avg <- All_samples_data  %>%
  group_by(Site) %>%
  dplyr::summarise(V1 = mean(V1, na.rm=TRUE), 
                   V2 = mean(V2, na.rm=TRUE), 
                   V3 = mean(V3, na.rm=TRUE), 
                   V4 = mean(V4, na.rm=TRUE), 
                   #V5 = mean(V5, na.rm=TRUE), 
                   # V6 = mean(V6, na.rm=TRUE), 
                   # V7 = mean(V7, na.rm=TRUE), 
                   # V8 = mean(V8, na.rm=TRUE), 
                   # V9 = mean(V9, na.rm=TRUE), 
                   # V10 = mean(V10, na.rm=TRUE),
                   # Lat = mean(Lat, na.rm=TRUE),
                   # Long = mean(Long, na.rm=TRUE),
                   # Lat_shift = mean(Lat_shift, na.rm=TRUE),
                   # Long_shift = mean(Long_shift, na.rm=TRUE),
                   # Avg.Temp = mean(Avg.Temp, na.rm=TRUE),
                   # Leaf_weight = mean(Leaf_weight, na.rm=TRUE),
                   # Ice_time = mean(Ice_time, na.rm=TRUE),
                   # O.HOM. = mean(O.HOM., na.rm=TRUE),
                   # E.HOM. = mean(E.HOM., na.rm=TRUE),
                   # N_SITES = mean(N_SITES, na.rm=TRUE),
                   # FIS = mean(F, na.rm=TRUE),
                   N_MISS = mean(N_MISS, na.rm=TRUE),
                   F_MISS = mean(F_MISS, na.rm=TRUE),
                   PC1 = mean(PC1, na.rm=TRUE),
                   PC2 = mean(PC2, na.rm=TRUE),
                   PC3 = mean(PC3, na.rm=TRUE),
                   PC4 = mean(PC4, na.rm=TRUE),
                   PC5 = mean(PC5, na.rm=TRUE),
                   PC6 = mean(PC6, na.rm=TRUE), 
                   Group = list(Geo_Region),
                   Paper_ID = Paper_ID
  )

Pop_avg$Pop <- Pop_avg$Site
All_pop_data <- left_join(pop_lat_long, Pop_avg, by="Paper_ID")

All_pop_data <- All_pop_data[-which(duplicated(All_pop_data)),] 

############################################
# calculate the max Admix group for each location
Admix_groups <- gather(All_pop_data, Group, Amount, V1:V4, factor_key=TRUE)

maxAdmixgroup <- Admix_groups  %>%
  group_by(Pop) %>%
  dplyr::summarise(Group = Group[which(Amount==max(Amount, na.rm=TRUE))])

Admix_groups$Mix <- paste0(Admix_groups$Group, Admix_groups$Pop)
maxAdmixgroup$Mix <- paste0(maxAdmixgroup$Group, maxAdmixgroup$Pop)

Final_AdmixGroups <- left_join(maxAdmixgroup, Admix_groups, by="Mix")


Max_Admix_Group <- c("V1", "V2", "V3", "V4")
Admix_Location <- c("North", "Center1", "Center2", "South")
Group_location <- as.data.frame(cbind(Max_Admix_Group, Admix_Location))

All_pop_data <- left_join(All_pop_data, Final_AdmixGroups, by="Site")

All_pop_data$Max_Admix_Group <- All_pop_data$Group.y
All_pop_data <- left_join(All_pop_data, Group_location, by="Max_Admix_Group")

#4 groups
Groups_summary_admix <- All_pop_data %>%
  group_by(Admix_Location) %>%
  dplyr::summarise(Locations = list(Pop.x), 
                   #Lat = mean(Lat.x, na.rm=TRUE),
                   #Long = mean(Long.x, na.rm=TRUE),
                   V1 = mean(V1, na.rm=TRUE), 
                   V2 = mean(V2, na.rm=TRUE), 
                   V3 = mean(V3, na.rm=TRUE), 
                   V4 = mean(V4, na.rm=TRUE))


###################################
# Set up colour scheme 

map_colours_5g <- c("deepskyblue", "yellow",  "green", "red")

Max_Admix_Group <- c("V1", "V2","V3", "V4")

map_col <- cbind(Max_Admix_Group, map_colours_5g)
map_col <- as.data.frame(map_col)

All_pop_data <- left_join(All_pop_data, map_col, by="Max_Admix_Group")

# join map_colours with individual sample data
Map_col_pop <- as.data.frame(cbind(All_pop_data$Pop, All_pop_data$map_colours_5g, All_pop_data$Max_Admix_Group))

colnames(Map_col_pop) <- c("Pop", "map_colours_5g", "Max_Admix_Group")

All_samples_data$Pop <- All_samples_data$Site

All_samples_data <- left_join(All_samples_data, Map_col_pop, by="Pop")

All_samples_data$map_colours_5g <- as.factor(All_samples_data$map_colours_5g)

#####################################

#Plot a PCA image

# PCA = 27.00  5.62  2.63  1.63  1.50  1.26  1.02  0.93  0.93  0.85  0.82  0.77  0.73  0.73  0.68  0.67
All_samples_data$Site <- as.factor(All_samples_data$Site)
jpeg("./Figures_data/Plots/PCA_PC1_PC2_maf1.jpg", width = 3000, height = 2700)
All_samples_data %>%
  ggplot(.,aes(x=PC1,y=PC2)) +
  geom_point(aes(color = Max_Admix_Group), size=15)  +
  geom_text(aes(label = ID_code), color = "black", fontface = 2, size = 25.4/72.27*15)+
  theme_classic()+
  theme(legend.text = element_text(color = "black", size = 50),
        axis.text=element_text(size=60 ,face="bold"),
        axis.title=element_text(size=60,face="bold")) +
  #guides(colour = guide_legend(override.aes = list(size=60)))+
  scale_colour_manual(values = map_colours_5g)+
  labs(y= "PC2", x = "PC1")
dev.off()

jpeg("./Figures_data/Plots/PCA_PC2_PC3_maf1.jpg", width = 3000, height = 2700)

All_samples_data %>%
  ggplot(.,aes(x=PC2,y=PC3)) + 
  geom_point(aes(color = Max_Admix_Group), size=15)  + 
  theme_classic()+
  theme(legend.text = element_text(color = "black", size = 70),
        axis.text=element_text(size=60 ,face="bold"),
        axis.title=element_text(size=60,face="bold")) + 
  guides(colour = guide_legend(override.aes = list(size=60)))+
  scale_colour_manual(values = map_colours_5g)+
  labs(y= "PC3 (2.63%)", x = "PC2 (5.62%)")
dev.off()

jpeg("./Figures_data/Plots/PCA_PC3_PC4_maf1.jpg", width = 3000, height = 2700)

All_samples_data %>%
  ggplot(.,aes(x=PC3,y=PC4)) + 
  geom_point(aes(color = Max_Admix_Group), size=15)  + 
  theme_classic()+
  theme(legend.text = element_text(color = "black", size = 70),
        axis.text=element_text(size=60 ,face="bold"),
        axis.title=element_text(size=60,face="bold")) + 
  guides(colour = guide_legend(override.aes = list(size=60)))+
  scale_colour_manual(values = map_colours_5g)+
  labs(y= "PC4 (1.63%)", x = "PC3 (2.63%)")
dev.off()

jpeg("./Figures_data/Plots/PCA_PC5_PC6_maf1.jpg", width = 3000, height = 2700)

All_samples_data %>%
  ggplot(.,aes(x=PC5,y=PC6)) + 
  geom_point(aes(color = Max_Admix_Group), size=15)  + 
  theme_classic()+
  theme(legend.text = element_text(color = "black", size = 70),
        axis.text=element_text(size=60 ,face="bold"),
        axis.title=element_text(size=60,face="bold")) + 
  guides(colour = guide_legend(override.aes = list(size=60)))+
  scale_colour_manual(values = map_colours_5g)+
  labs(y= "PC6 (1.26%)", x = "PC5 (1.50%)")
dev.off()

############################################

#Make Admixture barplot

#Excel: =RIGHT(A2,LEN(A2) - SEARCH("_", A2, SEARCH("_", A2) + 1))

All_lat_long_data <-  All_pop_data[,c(1:5)]
All_samples_data$Site <- All_samples_data$Pop

All_samples_data_map <- left_join(All_samples_data, All_lat_long_data, by="Site")

mergedAdmixTable <- All_samples_data_map[,c("Site", "ID_code","Latitude.x", "Longitude.x","V1", "V2","V3", "V4")]

rownames(mergedAdmixTable) <- mergedAdmixTable$ID_code

barNaming <- function(vec) {
  retVec <- vec
  for(k in 2:length(vec)) {
    if(vec[k-1] == vec[k])
      retVec[k] <- ""
  }
  return(retVec)
}

#----------------------------------
#4 Groups

ordered4 = mergedAdmixTable[order(mergedAdmixTable$Geo_Region),]

ordered0 = mergedAdmixTable[order(mergedAdmixTable$Site),]
ordered1 = ordered0[order(ordered0$V3),]
ordered2 = ordered1[order(ordered1$V2),]
ordered3 = ordered2[order(ordered2$V1),]
ordered4 = ordered3[order(ordered3$V4),]

ordered4 = mergedAdmixTable[order(mergedAdmixTable$Latitude.x),]


ordered4$Site <-  as.factor(ordered4$Site)

#map_colours_5g <- c("deepskyblue", "yellow",  "green", "red")
map_colours_5g <- c("green",  "yellow", "deepskyblue", "red")

jpeg("./Figures_data/Plots/Site_Admix_4_bar.jpg", width = 6000, height = 2000)
# bottom, left, top, and right
par(mar=c(30,10,4,4))
barplot(t(as.matrix(ordered4[,c(5:8)])), col=map_colours_5g, border=NA,
        names.arg=barNaming(ordered4$ID_code), las=2, cex.names=1, cex.axis=6)
dev.off()



#################################
# Explore CT SNPs and methylation
# https://yaaminiv.github.io/WGBS-Analysis-Part28/


#Intersect VCF with SNP locations and CG motif track
#BEDtools output has C/T and non-C/T SNPs
#Use grep to isolate C/T SNPs (tab between C and T required)

intersectBed \
-u \
-a /Volumes/web/spartina/project-gigas-oa-meth/output/BS-Snper/out/SNP-results.vcf \
-b CGmotif.bedGraph \
| grep "C T" \
> /Volumes/web/spartina/project-gigas-oa-meth/output/BS-Snper/out/CT-SNPs.tab


# for each sample
for f in zr3616*vcf
do
    /Users/Shared/bioinformatics/bedtools2/bin/intersectBed \
    -u \
    -a ${f} \
    -b /Users/yaamini/Documents/project-gigas-oa-meth/genome-feature-files/cgigas_uk_roslin_v1_fuzznuc_CGmotif.gff \
    | grep "C	T" \
    > ${f}_CT-SNPs.tab
    head ${f}_CT-SNPs.tab
    wc -l ${f}_CT-SNPs.tab
done


#######################################
# EWAS with GEM
# https://www.bioconductor.org/packages/devel/bioc/vignettes/GEM/inst/doc/user_guide.html

## install dependent packages
if(!require("methods")){
    install.packages("methods")
}
if(!require("ggplot2")){
    install.packages("ggplot2")
}
if(!require("Rcpp")){
    install.packages("Rcpp")
}
if(!require("digest")){
    install.packages("digest")
}
if(!require("devtools")){
    install.packages("devtools")
}

## install gem package
## suggest devtools version: devtools >= 1.11.1.9000
devtools::install_github("fastGEM/GEM")

#--------------------
# Setup input files

# covariates
# Sites

##        S1 S2 S3 S4 S5 S6 S7 S8 S9 S10
## Gender  2  2  2  1  2  1  1  2  1   1




#----------------------
# Run
library(GEM)
#GEM_GUI()

library(GEM)
DATADIR = system.file('extdata',package='GEM')
RESULTDIR = getwd()

env_file_name = paste(DATADIR, "env.txt", sep = .Platform$file.sep)
covariates_file_name = paste(DATADIR, "cov.txt", sep = .Platform$file.sep)
covariates_file_name_gxe = paste(DATADIR, "gxe.txt", sep = .Platform$file.sep)
methylation_file_name = paste(DATADIR, "methylation.txt", sep = .Platform$file.sep)
snp_file_name = paste(DATADIR, "snp.txt", sep = .Platform$file.sep)

Emodel_pv = 1
Gmodel_pv = 1e-04
GxEmodel_pv = 1

Emodel_result_file_name = paste(RESULTDIR, "Result_Emodel.txt", sep = .Platform$file.sep)
Gmodel_result_file_name = paste(RESULTDIR, "Result_Gmodel.txt", sep = .Platform$file.sep)
GxEmodel_result_file_name = paste(RESULTDIR, "Result_GxEmodel.txt", sep = .Platform$file.sep)

Emodel_qqplot_file_name = paste(RESULTDIR, "QQplot_Emodel.jpg", sep = .Platform$file.sep)

GEM_Emodel(env_file_name, covariates_file_name, methylation_file_name, Emodel_pv, Emodel_result_file_name,Emodel_qqplot_file_name)

GEM_Gmodel(snp_file_name, covariates_file_name, methylation_file_name, Gmodel_pv, Gmodel_result_file_name)

GEM_GxEmodel(snp_file_name, covariates_file_name_gxe, methylation_file_name, GxEmodel_pv, GxEmodel_result_file_name)





















#########################################
# Map of samples
# https://mikkovihtakari.github.io/PlotSvalbard/articles/PlotSvalbard_user_manual.html

# https://ggplot2-book.org/maps.html
# https://jakob.schwalb-willmann.de/basemaps/
All_pop_data  <- transform_coord(All_pop_data, lon = "Longitude.x", lat = "Latitude.x", bind = TRUE, proj.og = "+proj=longlat +datum=WGS84", proj.out = "+init=epsg:3857")
All_pop_data$Pop <- as.factor(All_pop_data$Pop)
#All_pop_data <- All_pop_data[,-c(21:22)]

library(basemaps)
library(mapedit)
library(ggmap)
library(ggnewscale)
# view all available maps
#get_maptypes()

# use draw_ext() to interactively draw an extent yourself
# ext <- draw_ext()

# load and return basemap map as class of choice, e.g. as image using magick:
#basemap_magick(ext, map_service = "osm", map_type = "topographic")

# set defaults for the basemap
#set_defaults(map_service = "osm_stamen", map_type = "toner")
#set_defaults(map_service = "osm_stamen", map_type = "terrain_bg")
#set_defaults(map_service = "osm", map_type = "topographic")
set_defaults(map_service = "osm_stamen", map_type = "terrain")

#All_pop_data <- left_join(All_pop_data, Map_col_pop, by="Pop")

jpeg("./Figures_data/Plots/Mimulus_samples_map_colours.jpg",width = 2700, height = 3300)
ggplot() + 
  basemap_gglayer(ext) +
  scale_fill_identity() + 
  coord_sf() +
  geom_point(data = All_pop_data, aes(x = lon.utm,  y = lat.utm, col = "Population"), color="yellow", size = 20) +
  geom_text(data = All_pop_data, aes(x = lon.utm, y = lat.utm, label = Site), color = "black", fontface = 2, size = 25.4/72.27*30)
dev.off()

#https://eliocamp.github.io/codigo-r/2018/09/multiple-color-and-fill-scales-with-ggplot2/

jpeg("./Figures_data/Plots/Admix_map4_take1_lat_long.jpg", width = 900, height = 1500)
ggplot(data = All_pop_data, aes(x = lon.utm, y = lat.utm)) +
  #basemap_gglayer(ext) +
  #coord_sf() +
  #scale_fill_identity() +
  #new_scale_colour() +
  geom_scatterpie(aes(x = lon.utm, y = lat.utm, group=Site, r=20000), data = All_pop_data, cols=colnames(All_pop_data[,c(6:9)]), size = 0.1)+
  scale_fill_manual(values=map_colours_5g) +
  geom_text(data = All_pop_data, aes(x = lon.utm, y = lat.utm, label = Site), color = "black", fontface = 2, size = 25.4/72.27*10) 
dev.off()

#cannot get the two scales for the map and the pie charts to work

*************************
  
##############################
# Write out samples file

write.table(All_samples_data, file = "./Figures_data/All_samples_data.txt", quote = FALSE, row.names=FALSE, col.names=TRUE, sep="\t")

#####################################
# Admixture map - shifted to see all groups

#library(devtools)
#devtools::install_github("MikkoVihtakari/PlotSvalbard", upgrade = "never")
#library(PlotSvalbard)
# https://mikkovihtakari.github.io/PlotSvalbard/articles/PlotSvalbard_user_manual.html
#library(scatterpie)

# shifted coordinates
All_samples_data <- transform_coord(All_samples_data, lon = "Long_shift", lat = "Lat_shift", bind = TRUE, proj.og = "+proj=longlat +datum=WGS84", proj.out = "+init=epsg:3995")
All_samples_data$lon_shift.utm <- All_samples_data$lon.utm 
All_samples_data$lat_shift.utm <- All_samples_data$lat.utm 
All_samples_data <- subset(All_samples_data, select = -c(lon.utm, lat.utm))

# regular coordinates
All_samples_data <- transform_coord(All_samples_data, lon = "Long", lat = "Lat", bind = TRUE, proj.og = "+proj=longlat +datum=WGS84", proj.out = "+init=epsg:3995")

# make dataset of just what needs to be plotted
shift_Admix_lat_long <- subset(All_samples_data, select = c(ID_code, lon.utm, lat.utm, lon_shift.utm, lat_shift.utm, Pop, V1, V2, V3, V4, V5))

Pop <- c("V1", "V2", "V3", "V4", "V5")

jpeg("./Figures_data/Plots/Admix_map5_take1_shiftedlat_long.jpg", width = 2000, height = 1700)
basemap("panarctic", limits=50) + 
  geom_scatterpie(aes(x = lon_shift.utm, y = lat_shift.utm, group = ID_code, r=200000), data = shift_Admix_lat_long, cols = Pop, size = 0.9) +
  geom_point(aes(x = lon_shift.utm, y = lat_shift.utm), data = shift_Admix_lat_long, col="white", size=17) +
  geom_point(aes(x = lon.utm, y = lat.utm), data = shift_Admix_lat_long, col="green", size=6) +
  scale_fill_manual(values=map_colours_5g)+
  geom_text(data = shift_Admix_lat_long, aes(x = lon_shift.utm, y = lat_shift.utm, label = Pop), color = "black", fontface = 2, size = 25.4/72.27*20)
dev.off()

#---------------------------------
# Averaged admixture map - not plotting each individual

Pop <- c("V1", "V2", "V3", "V4", "V5")

# map showing population structure makeup
jpeg("./Figures_data/Plots/Admix_map5g_take1_avg.jpg", width = 2000, height = 1700)
basemap("panarctic", limits=50) + 
  geom_scatterpie(aes(x = lon.utm, y = lat.utm, group = Pop, r=200000), data = All_pop_data, cols = Pop, size = 0.5) +
  scale_fill_manual(values=map_colours_5g)
dev.off()

#----------------------------
All_pop_data <- subset(All_pop_data, select = -c(lon.utm, lat.utm))

All_pop_data <- transform_coord(All_pop_data, lon = "Long_shift", lat = "Lat_shift", bind = TRUE, proj.og = "+proj=longlat +datum=WGS84", proj.out = "+init=epsg:3995")
All_pop_data$lon_shift.utm <- All_pop_data$lon.utm 
All_pop_data$lat_shift.utm <- All_pop_data$lat.utm 
All_pop_data <- subset(All_pop_data, select = -c(lon.utm, lat.utm))

All_pop_data <- transform_coord(All_pop_data, lon = "Long.x", lat = "Lat.x", bind = TRUE, proj.og = "+proj=longlat +datum=WGS84", proj.out = "+init=epsg:3995")

jpeg("./Figures_data/Plots/Admix_map5_take1_shiftedlat_long_avg.jpg", width = 2000, height = 1700)
basemap("panarctic", limits=50) + 
  geom_scatterpie(aes(x = lon_shift.utm, y = lat_shift.utm, group = Pop, r=200000), data = All_pop_data, cols = Pop, size = 0.9) +
  geom_point(aes(x = lon_shift.utm, y = lat_shift.utm), data = All_pop_data, col="white", size=17) +
  geom_point(aes(x = lon.utm, y = lat.utm), data = All_pop_data, col="green", size=6) +
  scale_fill_manual(values=map_colours_5g)+
  geom_text(data = All_pop_data[-which(All_pop_data$Pop=="AlexOld"), ], aes(x = lon_shift.utm, y = lat_shift.utm, label = Pop), color = "black", fontface = 2, size = 25.4/72.27*20)
dev.off()

######################################

# Kluane
#check correlation with %BC and elevation

All_samples_data
Kluane <- All_samples_data[which(All_samples_data$Pop.x=="KL"|All_samples_data$Pop.x=="PC"),]

# Saximontana is V5 for Take 1

mod1 = lm(V3~Elevation, data = Kluane)
modsum = summary(mod1)

# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  1.292e+00  1.440e-01   8.971 7.42e-08 ***
#   Elevation   -6.859e-04  9.243e-05  -7.421 1.00e-06 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.07808 on 17 degrees of freedom
# Multiple R-squared:  0.7641,	Adjusted R-squared:  0.7502 
# F-statistic: 55.07 on 1 and 17 DF,  p-value: 9.995e-07

jpeg("./Figures_data/Plots/Elevation_BC_5.jpg", width = 3000, height = 1400)
par(mar=c(20,20,4,4))
plot(Kluane$Elevation, Kluane$V3, pch=20, cex = 10, mgp=c(10,5,0), cex.lab=5, cex.axis=5, xlab="Elevation(m)", ylab="BC Admixture Proportion")
abline(mod1)
dev.off()


orderedKluane = Kluane[order(Kluane$Elevation),]

barNaming <- function(vec) {
  retVec <- vec
  for(k in 2:length(vec)) {
    if(vec[k-1] == vec[k])
      retVec[k] <- ""
  }
  return(retVec)
}

jpeg("./Figures_data/Plots/Admix_Kluane_bar5.jpg", width = 1000, height = 707)
barplot(t(as.matrix(orderedKluane[,c(3:7)])), col=map_colours_5g, border=NA,
        names.arg=barNaming(orderedKluane$ID_code), las=2, cex.names=1.7)
dev.off()

#################################################
# Climate, ice and nucleotide diversity
# Icetime.x
# tavglong
# SumWindPi

# check site and window pi correlate
jpeg("./Figures_data/Plots/SitevsWIndPi.jpg", width = 1000, height = 707)
plot(All_pop_data$SumWindPi ~ All_pop_data$SumSitePi, pch=20, xlab="Sites Pi", ylab="Window Pi")
dev.off()

# Models
mod3 = lm(SumWindPi ~ tavglong, data = All_pop_data)
summary(mod3)
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -14.2604  -1.5667   0.1728   1.9179  11.5762 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  16.8080     1.7931   9.374 5.95e-11 ***
#   tavglong      0.5029     0.2390   2.104   0.0428 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 4.205 on 34 degrees of freedom
# (6 observations deleted due to missingness)
# Multiple R-squared:  0.1152,	Adjusted R-squared:  0.08918 
# F-statistic: 4.427 on 1 and 34 DF,  p-value: 0.04285


mod4 = lm(Ice_time.x ~ tavglong, data = All_pop_data)
summary(mod4)

# Residuals:
#   Min     1Q Median     3Q    Max 
# -9.730 -2.978 -1.069  1.762  8.380 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  10.9930     1.5866   6.928 3.07e-08 ***
#   tavglong      0.1624     0.1863   0.872    0.389    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 4.565 on 38 degrees of freedom
# (2 observations deleted due to missingness)
# Multiple R-squared:  0.01962,	Adjusted R-squared:  -0.006184 
# F-statistic: 0.7603 on 1 and 38 DF,  p-value: 0.3887


mod5 = lm(SumWindPi ~ Ice_time.x, data = All_pop_data)
summary(mod5)

# Residuals:
#   Min       1Q   Median       3Q      Max 
# -11.3903  -1.5145   0.0607   1.0588  14.5881 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  15.7977     1.8528   8.526 3.66e-10 ***
#   Ice_time.x    0.3520     0.1417   2.484   0.0178 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 4.082 on 36 degrees of freedom
# (4 observations deleted due to missingness)
# Multiple R-squared:  0.1463,	Adjusted R-squared:  0.1226 
# F-statistic: 6.171 on 1 and 36 DF,  p-value: 0.01777


# Plots
jpeg("./Figures_data/Plots/WindPivsTemp.jpg", width = 1000, height = 707)
plot(All_pop_data$SumWindPi ~ All_pop_data$tavglong, col="lightblue", pch=19, cex=2, xlab="Average Temperature (degrees Celsius)", ylab="Sites Pi")
abline(mod3, col="red", lwd=3)
text(SumWindPi ~ tavglong, labels=Pop, data=All_pop_data, cex=0.9, font=2)
dev.off()

jpeg("./Figures_data/Plots/Ice_temp.jpg", width = 1000, height = 707)
plot(All_pop_data$Ice_time.x ~ All_pop_data$tavglong, col="lightblue", pch=19, cex=2, xlab="Average Temperature (degrees Celsius)", ylab="Icetime")
abline(mod4, col="red", lwd=3)
text(Ice_time.x ~ tavglong, labels=Pop,data=All_pop_data, cex=0.9, font=2)
dev.off()

jpeg("./Figures_data/Plots/IcetimesWIndPi.jpg", width = 1000, height = 707)
plot(All_pop_data$SumWindPi ~ All_pop_data$Ice_time.x, col="lightblue", pch=19, cex=2, xlab="Ice retreat time (thousands of years)", ylab="Sites Pi")
abline(mod5, col="red", lwd=3)
text(SumWindPi ~ Ice_time.x, labels=Pop, data=All_pop_data, cex=0.9, font=2)
dev.off()

##################
# make plots of Pi, Tajimas D

# order population averages by longitude
All_pop_data_ord <- All_pop_data[order(All_pop_data$Long.x),]

# remove Sverdrup New
All_pop_data_ord_sub <- All_pop_data_ord[-which(All_pop_data_ord$Pop=="SVN"),]

# remove NA values
All_pop_data_ord_sub <- All_pop_data_ord_sub[-which(is.na(All_pop_data_ord_sub$SumWindPi)),]

# make Pop a factor
All_pop_data_ord_sub$Pop <- as.factor(All_pop_data_ord_sub$Pop)

# plot Nucleotide diversity
jpeg("./Figures_data/Plots/AvgPi_noSverdrup.jpg", width = 2200, height = 707)
plot(All_pop_data_ord_sub$Pop, All_pop_data_ord_sub$SumWindPi)
dev.off()

jpeg("./Figures_data/Plots/AvgTajima_noSverdrup.jpg", width = 1700, height = 707)
plot(All_pop_data_ord_sub$Pop, All_pop_data_ord_sub$TajimaAvg, cex.axis=2, cex.lab=1.5, las=2)
dev.off()

##################
# make barplots of heterozygotsity with per individual data

# All_samples_data$F, All_samples_data$O.HOM. , All_samples_data$E.HOM.

All_samples_data$Pop <- as.factor(All_samples_data$Pop)

# jpeg("./Figures_data/Plots/FIS.jpg", width = 2200, height = 707)
# plot(All_samples_data$Pop, All_samples_data$F, col="blue")
# dev.off()

jpeg("./Figures_data/Plots/Obser_Homozygosity.jpg", width = 2200, height = 707)
plot(All_samples_data$Pop, All_samples_data$O.HOM., col="blue")
dev.off()

jpeg("./Figures_data/Plots/Exp_Homozygosity.jpg", width = 2200, height = 707)
plot(All_samples_data$Pop, All_samples_data$E.HOM., col="blue")
dev.off()


#-------
# order by longitude
Het_data <-  as.data.frame(cbind(All_samples_data$Pop.x, All_samples_data$F, All_samples_data$Long, All_samples_data$map_colours_5g))
colnames(Het_data) <- c("Pop", "F", "Long", "map_colours_5g")
Het_data$F <- as.numeric(Het_data$F)
Het_data$Long <- as.numeric(Het_data$Long)

Het_data_ord <- Het_data[order(Het_data$Long),]

Het_data_ord$Pop <- as.factor(Het_data_ord$Pop)
levels(Het_data_ord$Pop)
Het_data_ord$Pop <- factor(Het_data_ord$Pop, levels = c("GEN", "ATQ", "BARD", "MAT", "MNT", "IMN", "SAG" , "DEN", "MIL",  "QHI",  "KL", "PC",  "PEA", 
                                                        "KUQ",   "YAM", "AXE", "EUR",  "FOS", "GF", "SVO" , "BY" , "AlexNew", "AlexOld" ,  "HAZ",  "Iq",
                                                        "CR", "Kik", "DLG", "DQG",  "025", "IG", "ZAC", "LON", "PET", "LAJ", "SW", "SAM", "YED"))

Het_cols <-  c("purple", "deepskyblue", "deepskyblue", "deepskyblue", "deepskyblue", "deepskyblue", "deepskyblue" , "deepskyblue", "deepskyblue",  "deepskyblue",  "deepskyblue", "deepskyblue",  "deepskyblue", 
               "deepskyblue",   "deepskyblue", "deepskyblue", "yellow",  "yellow", "yellow", "yellow" , "yellow" , "yellow", "yellow" ,  "yellow",  "orange",
               "orange", "orange", "orange", "orange",  "orange", "orange", "yellow", "yellow", "yellow", "yellow", "yellow", "red", "red")

# jpeg("./Figures_data/Plots/FIS_ordered.jpg", width = 2200, height = 707)
# plot(Het_data_ord$Pop, Het_data_ord$F, col=Het_cols)
# dev.off()


# par(mar = c(bottom, left, top, right)) 
jpeg("./Figures_data/Plots/Fis_col_box.jpg", width = 4000, height = 1700)
par(mar= c(20,15,4,4))
plot(F ~ Pop,  data = Het_data_ord, xlab="", ylab="", cex=3, cex.lab=6, cex.axis=5, las=2, col = Het_cols)
stripchart(F ~ Pop, vertical = TRUE, method = "jitter", pch = 16,las = 2, col = "blue", cex = 1.5, data = Het_data_ord, add=TRUE)
lines(Het_data_ord$Pop, y=rep(0, length(Het_data_ord$Pop)), lwd=3)
dev.off()


##########################################
# Estimate haplotypes for seedlings?

# https://odelaneau.github.io/shapeit5/
# https://odelaneau.github.io/shapeit5/docs/tutorials/simulated/

#####################################
# Genes associated with methylation differences based on natural variation
# https://epidiverse.gitbook.io/project/-MfxkdBDZggX_vc_sG5l/lectures/natural-variation-of-methylation-detlef-weigel
# https://github.com/EpiDiverse/ewas
# https://github.com/EpiDiverse/snp



