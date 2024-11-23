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

cd /home/celphin/scratch/Dryas/BS-Snper

git clone https://github.com/hellbelly/BS-Snper.git

cd BS-Snper
sh BS-Snper.sh

# Make sure the executable file rrbsSnp is generated.

# /cvmfs/soft.computecanada.ca/gentoo/2023/x86-64-v3/usr/x86_64-pc-linux-gnu/binutils-bin/2.41/ld: cannot find -lcurses: No such file or directory
# collect2: error: ld returned 1 exit status
# make[1]: *** [Makefile:47: samtools] Error 1
# make[1]: Leaving directory '/scratch/celphin/Dryas/BS-Snper/BS-Snper/samtools-0.1.19'
# make: *** [Makefile:27: all-recur] Error 1
# g++ -O2 -o rrbsSnp main.o sam_funcs.o hash_funcs.o chrome_funcs.o  -m64 -I./samtools-0.1.19/ -L./samtools-0.1.19/ -lbam -lz -lpthread

#-----------------------------
# copy over the raw filtered bam data to run this on



# copy over the new reference genome



#---------------------------------
# run on one indvidual

salloc 

module load 

cd 

perl BS-Snper.pl <sorted_bam_file> --fa <reference_file> --output <snp_result_file> --methcg <meth_cg_result_file> --methchg <meth_chg_result_file> --methchh <meth_chh_result_file> --minhetfreq 0.1 --minhomfreq 0.85 --minquali 15 --mincover 10 --maxcover 1000 --minread2 2 --errorate 0.02 --mapvalue 20 >SNP.out 2>ERR.log

