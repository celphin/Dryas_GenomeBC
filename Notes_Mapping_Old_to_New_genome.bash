######################################
# RagTag Dryas genomes
# June 2024
###################################
# download genomes
cd ~/scratch/Dryas/Dryas_genomes/

tmux new-session -s CrossMap
tmux attach-session -t CrossMap

salloc -c40 --time 2:55:00 --mem 190000m --account def-cronk

cd /home/celphin/scratch/RagTag/
source ~/RagTag/bin/activate

module load StdEnv/2020
module load python/3.10.2
module load minimap2/2.24
module load scipy-stack/2022a

pip install RagTag

cd ~/scratch/Dryas/Dryas_genomes/

# Try running
/home/celphin/scratch/RagTag/RagTag-master/ragtag.py \
scaffold DoctH0_Main.fasta Dryas_octopetala_H1.supercontigs.fa \
-t 39 -u -o ./Old_Dryas_octopetala_H1_ragtag_output/


######################################
# check stats
# https://github.com/sanger-pathogens/assembly-stats

cd ~
git clone https://github.com/sanger-pathogens/assembly-stats.git
cd assembly-stats
mkdir build
cd build
cmake -DINSTALL_DIR:PATH=/home/celphin/assembly-stats/ ..
make
make test
make install

#-----------------------
module load StdEnv/2020 
module load gcc/9.3.0
module load r-bundle-bioconductor/3.14
cd ~/assembly-stats/build

./assembly-stats ~/scratch/Dryas/Dryas_genomes/Dryas_octopetala_H1.supercontigs.fa
# stats for /home/celphin/scratch/Dryas/Dryas_genomes/Dryas_octopetala_H1.supercontigs.fa
# sum = 233557465, n = 179, ave = 1304790.31, largest = 19685959
# N50 = 12105649, n = 8
# N60 = 11145057, n = 10
# N70 = 8228389, n = 12
# N80 = 4503589, n = 16
# N90 = 2104560, n = 24
# N100 = 9458, n = 179
# N_count = 134679
# Gaps = 190

./assembly-stats ~/scratch/Dryas/Dryas_genomes/Old_Dryas_octopetala_H1_ragtag_output/ragtag.scaffold.fasta
# stats for /home/celphin/scratch/Dryas/Dryas_genomes/Old_Dryas_octopetala_H1_ragtag_output/ragtag.scaffold.fasta
# sum = 233562465, n = 129, ave = 1810561.74, largest = 29810504
# N50 = 26945421, n = 5
# N60 = 26945421, n = 5
# N70 = 22319059, n = 6
# N80 = 20495270, n = 8
# N90 = 19685959, n = 9
# N100 = 9458, n = 129
# N_count = 139679
# Gaps = 227


#########################
# RepeatOBserver on new genome
cd /home/celphin/scratch/repeats/auto_script/
cp ~/scratch/Dryas/Dryas_genomes/Old_Dryas_octopetala_H1_ragtag_output/ragtag.scaffold.fasta OldRagtagDryOcto.fasta

cat << EOF > Auto_OldRagtagDryOcto.sh
#!/bin/bash
#SBATCH --account=def-cronk
#SBATCH --time=5:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=15
#SBATCH --mem=128000M

module load StdEnv/2020
module load seqkit/2.3.1
module load emboss/6.6.0
module load r/4.1.2

srun Setup_Run_Repeats.sh -i OldRagtagDryOcto -f OldRagtagDryOcto.fasta -h H0 -c 15 -m 128000M -g FALSE

EOF

sbatch Auto_OldRagtagDryOcto.sh

###############################
# Mapping all other files between genomes
# make chain file
# https://genome.ucsc.edu/goldenPath/help/chain.html
# Generate alignments with the cigar string attached to the cg:Z: tag. 
# These can be made by several aligners, including minimap2 -c

module load minimap2
minimap2 -c  DoctH0_Main.fasta Dryas_old.fa > results.paf
# Real time: 313.507 sec; CPU: 897.067 sec; Peak RSS: 6.508 GB
minimap2 -c Dryas_old.fa  DoctH0_Main.fasta > rev_results.paf
minimap2 -c Dryas_old.fa OldRagtagDryOcto.fasta > Ragtag_results.paf

module load rust
git clone https://github.com/AndreaGuarracino/paf2chain
cd paf2chain
cargo install --force --path .
export PATH=$PATH:/home/celphin/.cargo/bin

paf2chain -i results.paf > Dryas_Old_New.chain
paf2chain -i rev_results.paf > Dryas_New_Old.chain
paf2chain -i Ragtag_results.paf > Dryas_RagTag_Old.chain


########################################
# CrossMap: https://crossmap.readthedocs.io/en/latest/#installation

tmux new-session -s CrossMap
tmux attach-session -t CrossMap

cd /home/celphin/scratch/Dryas/CrossMap/

module load StdEnv/2023 python/3.12.4
avail_wheels CrossMap
pip install CrossMap

#---------------------
# Get files
# Cedar

mkdir file_for_converting

cd /home/celphin/scratch/Dryas/CrossMap/file_for_converting

cp /home/celphin/scratch/Dryas/bedgraph_old_ref/W_*.bedGraph .
cp /home/celphin/scratch/Dryas/bedgraph_old_ref/C_*.bedGraph .

cp /home/celphin/scratch/Dryas/CHG_CHH/data/bedgraph/W_*.bedGraph .
cp /home/celphin/scratch/Dryas/CHG_CHH/data/bedgraph/C_*.bedGraph .

cp /home/celphin/scratch/Dryas/MS_Dryas_Merged_Data/original_data/Dryas_octopetala_H1.gff3 . 
cp /home/celphin/scratch/Oxyria/EDTA/Dryas_octopetala_H1.supercontigs.fa.mod.EDTA.TEanno.gff3 .

cp /home/celphin/scratch/Dryas/snpEff/Dryas_DMRs/CHH/*_CHH.bedGraph ./CHH/

cp /home/celphin/scratch/Dryas/BS-Snper/pop_gen/Dryas_filtered_biasSNPs.vcf .

#---------------------
# Map to RagTag reference
# to run
cd /home/celphin/scratch/Dryas/CrossMap/file_for_converting

# genes
CrossMap gff ../Dryas_RagTag_Old.chain Dryas_octopetala_H1.gff3 RagTag_Dryas_genes.gff3

# TEs
CrossMap gff ../Dryas_RagTag_Old.chain Dryas_octopetala_H1.supercontigs.fa.mod.EDTA.TEanno.gff3 RagTag_Dryas_TEs.gff3

# Bedgraphs CpG per site_treatment
CrossMap bed ../Dryas_RagTag_Old.chain C_ALAS.bedGraph RagTag_C_ALAS.bedGraph
CrossMap bed ../Dryas_RagTag_Old.chain C_CASS.bedGraph RagTag_C_CASS.bedGraph
CrossMap bed ../Dryas_RagTag_Old.chain C_DRY.bedGraph RagTag_C_DRY.bedGraph
CrossMap bed ../Dryas_RagTag_Old.chain C_FERT.bedGraph RagTag_C_FERT.bedGraph
CrossMap bed ../Dryas_RagTag_Old.chain C_MEAD.bedGraph RagTag_C_MEAD.bedGraph
CrossMap bed ../Dryas_RagTag_Old.chain C_SVAL.bedGraph RagTag_C_SVAL.bedGraph
CrossMap bed ../Dryas_RagTag_Old.chain C_SwedC.bedGraph RagTag_C_SwedC.bedGraph
CrossMap bed ../Dryas_RagTag_Old.chain C_WILL.bedGraph RagTag_C_WILL.bedGraph

CrossMap bed ../Dryas_RagTag_Old.chain W_ALAS.bedGraph RagTag_W_ALAS.bedGraph
CrossMap bed ../Dryas_RagTag_Old.chain W_CASS.bedGraph RagTag_W_CASS.bedGraph
CrossMap bed ../Dryas_RagTag_Old.chain W_DRY.bedGraph RagTag_W_DRY.bedGraph
CrossMap bed ../Dryas_RagTag_Old.chain W_FERT.bedGraph RagTag_W_FERT.bedGraph
CrossMap bed ../Dryas_RagTag_Old.chain W_MEAD.bedGraph RagTag_W_MEAD.bedGraph
CrossMap bed ../Dryas_RagTag_Old.chain W_SVAL.bedGraph RagTag_W_SVAL.bedGraph
CrossMap bed ../Dryas_RagTag_Old.chain W_SwedC.bedGraph RagTag_W_SwedC.bedGraph
CrossMap bed ../Dryas_RagTag_Old.chain W_WILL.bedGraph RagTag_W_WILL.bedGraph

# Bedgraphs CHH per site_treatment
CrossMap bed ../Dryas_RagTag_Old.chain C_ALAS_CHH.bedGraph RagTag_C_ALAS_CHH.bedGraph
CrossMap bed ../Dryas_RagTag_Old.chain C_CASS_CHH.bedGraph RagTag_C_CASS_CHH.bedGraph
CrossMap bed ../Dryas_RagTag_Old.chain C_DRY_CHH.bedGraph RagTag_C_DRY_CHH.bedGraph
CrossMap bed ../Dryas_RagTag_Old.chain C_FERT_CHH.bedGraph RagTag_C_FERT_CHH.bedGraph
CrossMap bed ../Dryas_RagTag_Old.chain C_MEAD_CHH.bedGraph RagTag_C_MEAD_CHH.bedGraph
CrossMap bed ../Dryas_RagTag_Old.chain C_SVAL_CHH.bedGraph RagTag_C_SVAL_CHH.bedGraph
CrossMap bed ../Dryas_RagTag_Old.chain C_SwedC_CHH.bedGraph RagTag_C_SwedC_CHH.bedGraph
CrossMap bed ../Dryas_RagTag_Old.chain C_WILL_CHH.bedGraph RagTag_C_WILL_CHH.bedGraph

CrossMap bed ../Dryas_RagTag_Old.chain W_ALAS_CHH.bedGraph RagTag_W_ALAS_CHH.bedGraph
CrossMap bed ../Dryas_RagTag_Old.chain W_CASS_CHH.bedGraph RagTag_W_CASS_CHH.bedGraph
CrossMap bed ../Dryas_RagTag_Old.chain W_DRY_CHH.bedGraph RagTag_W_DRY_CHH.bedGraph
CrossMap bed ../Dryas_RagTag_Old.chain W_FERT_CHH.bedGraph RagTag_W_FERT_CHH.bedGraph
CrossMap bed ../Dryas_RagTag_Old.chain W_MEAD_CHH.bedGraph RagTag_W_MEAD_CHH.bedGraph
CrossMap bed ../Dryas_RagTag_Old.chain W_SVAL_CHH.bedGraph RagTag_W_SVAL_CHH.bedGraph
CrossMap bed ../Dryas_RagTag_Old.chain W_SwedC_CHH.bedGraph RagTag_W_SwedC_CHH.bedGraph
CrossMap bed ../Dryas_RagTag_Old.chain W_WILL_CHH.bedGraph RagTag_W_WILL_CHH.bedGraph

# CpG DMRs (W/C, Low/High ...)
CrossMap bed ../Dryas_RagTag_Old.chain ./BedGraphs_From_Metilene/Alaska_W_C.bedGraph ./BedGraphs_From_Metilene/RagTag_Alaska_W_C.bedGraph
CrossMap bed ../Dryas_RagTag_Old.chain ./BedGraphs_From_Metilene/Mat_Sen.bedGraph ./BedGraphs_From_Metilene/RagTag_Mat_Sen.bedGraph
CrossMap bed ../Dryas_RagTag_Old.chain ./BedGraphs_From_Metilene/Svalbard_W_C.bedGraph ./BedGraphs_From_Metilene/RagTag_Svalbard_W_C.bedGraph
CrossMap bed ../Dryas_RagTag_Old.chain ./BedGraphs_From_Metilene/SVAL_NUN_Sites_intersect_DMRs.bedgraph ./BedGraphs_From_Metilene/RagTag_SVAL_NUN_Sites_intersect_DMRs.bedgraph
CrossMap bed ../Dryas_RagTag_Old.chain ./BedGraphs_From_Metilene/NUN_ALAS_Sites_intersect_DMRs.bedgraph ./BedGraphs_From_Metilene/RagTag_NUN_ALAS_Sites_intersect_DMRs.bedgraph
CrossMap bed ../Dryas_RagTag_Old.chain ./BedGraphs_From_Metilene/SVAL_SWED_Sites_intersect_DMRs.bedgraph ./BedGraphs_From_Metilene/RagTag_SVAL_SWED_Sites_intersect_DMRs.bedgraph
CrossMap bed ../Dryas_RagTag_Old.chain ./BedGraphs_From_Metilene/Nunavut_W_C.bedGraph ./BedGraphs_From_Metilene/RagTag_Nunavut_W_C.bedGraph
CrossMap bed ../Dryas_RagTag_Old.chain ./BedGraphs_From_Metilene/SWED_ALAS_Sites_intersect_DMRs.bedgraph ./BedGraphs_From_Metilene/RagTag_SWED_ALAS_Sites_intersect_DMRs.bedgraph
CrossMap bed ../Dryas_RagTag_Old.chain ./BedGraphs_From_Metilene/Parent_W_C.bedGraph ./BedGraphs_From_Metilene/RagTag_Parent_W_C.bedGraph
CrossMap bed ../Dryas_RagTag_Old.chain ./BedGraphs_From_Metilene/Sweden_W_C.bedGraph ./BedGraphs_From_Metilene/RagTag_Sweden_W_C.bedGraph
CrossMap bed ../Dryas_RagTag_Old.chain ./BedGraphs_From_Metilene/SE_L_H.bedGraph ./BedGraphs_From_Metilene/RagTag_SE_L_H.bedGraph
CrossMap bed ../Dryas_RagTag_Old.chain ./BedGraphs_From_Metilene/Wild_Lat_L_H.bedGraph ./BedGraphs_From_Metilene/RagTag_Wild_Lat_L_H.bedGraph
CrossMap bed ../Dryas_RagTag_Old.chain ./BedGraphs_From_Metilene/SE_W_C.bedGraph ./BedGraphs_From_Metilene/RagTag_SE_W_C.bedGraph
CrossMap bed ../Dryas_RagTag_Old.chain ./BedGraphs_From_Metilene/Wild_W_C.bedGraph ./BedGraphs_From_Metilene/RagTag_Wild_W_C.bedGraph


# Bedgraphs_Intersected_Subtracted
CrossMap bed ../Dryas_RagTag_Old.chain ./Bedgraphs_Intersected_Subtracted/ALEX_SVAL_ALAS_Sites_intersect_DMRs.bedgraph ./Bedgraphs_Intersected_Subtracted/RagTag_ALEX_SVAL_ALAS_Sites_intersect_DMRs.bedgraph
CrossMap bed ../Dryas_RagTag_Old.chain ./Bedgraphs_Intersected_Subtracted/ALEX_SVAL_ALAS_Sites_intersect_DMRs.bedgraph ./Bedgraphs_Intersected_Subtracted/RagTag_ALEX_SVAL_ALAS_Sites_intersect_DMRs.bedgraph
CrossMap bed ../Dryas_RagTag_Old.chain ./Bedgraphs_Intersected_Subtracted/intersect_SE_W_C_Wild_W_C.bedGraph ./Bedgraphs_Intersected_Subtracted/RagTag_intersect_SE_W_C_Wild_W_C.bedGraph
CrossMap bed ../Dryas_RagTag_Old.chain ./Bedgraphs_Intersected_Subtracted/ALEX_SVAL_SWED_Sites_intersect_DMRs.bedgraph ./Bedgraphs_Intersected_Subtracted/RagTag_ALEX_SVAL_SWED_Sites_intersect_DMRs.bedgraph
CrossMap bed ../Dryas_RagTag_Old.chain ./Bedgraphs_Intersected_Subtracted/intersect_Wild_W_C_Mat_Sen.bedGraph ./Bedgraphs_Intersected_Subtracted/RagTag_intersect_Wild_W_C_Mat_Sen.bedGraph
CrossMap bed ../Dryas_RagTag_Old.chain ./Bedgraphs_Intersected_Subtracted/ALEX_SWED_ALAS_Sites_intersect_DMRs.bedgraph ./Bedgraphs_Intersected_Subtracted/RagTag_ALEX_SWED_ALAS_Sites_intersect_DMRs.bedgraph
CrossMap bed ../Dryas_RagTag_Old.chain ./Bedgraphs_Intersected_Subtracted/SVAL_SWED_ALAS_Sites_intersect_DMRs.bedgraph ./Bedgraphs_Intersected_Subtracted/RagTag_SVAL_SWED_ALAS_Sites_intersect_DMRs.bedgraph
CrossMap bed ../Dryas_RagTag_Old.chain ./Bedgraphs_Intersected_Subtracted/ALL_Sites_intersect_DMRs.bedgraph ./Bedgraphs_Intersected_Subtracted/RagTag_ALL_Sites_intersect_DMRs.bedgraph
CrossMap bed ../Dryas_RagTag_Old.chain ./Bedgraphs_Intersected_Subtracted/total_subtract_SE_W_C_P_W_C.bedGraph ./Bedgraphs_Intersected_Subtracted/RagTag_total_subtract_SE_W_C_P_W_C.bedGraph
CrossMap bed ../Dryas_RagTag_Old.chain ./Bedgraphs_Intersected_Subtracted/intersect_SE_L_H_Wild_L_H.bedGraph ./Bedgraphs_Intersected_Subtracted/RagTag_intersect_SE_L_H_Wild_L_H.bedGraph
CrossMap bed ../Dryas_RagTag_Old.chain ./Bedgraphs_Intersected_Subtracted/total_subtract_SE_W_C_SE_L_H.bedGraph ./Bedgraphs_Intersected_Subtracted/RagTag_total_subtract_SE_W_C_SE_L_H.bedGraph
CrossMap bed ../Dryas_RagTag_Old.chain ./Bedgraphs_Intersected_Subtracted/intersect_SE_W_C_P_W_C.bedGraph ./Bedgraphs_Intersected_Subtracted/RagTag_intersect_SE_W_C_P_W_C.bedGraph
CrossMap bed ../Dryas_RagTag_Old.chain ./Bedgraphs_Intersected_Subtracted/total_subtract_W_C_Mat_Sen.bedGraph ./Bedgraphs_Intersected_Subtracted/RagTag_total_subtract_W_C_Mat_Sen.bedGraph
CrossMap bed ../Dryas_RagTag_Old.chain ./Bedgraphs_Intersected_Subtracted/intersect_SE_W_C_SE_L_H.bedGraph ./Bedgraphs_Intersected_Subtracted/RagTag_intersect_SE_W_C_SE_L_H.bedGraph


# CHH DMRs
CrossMap bed ../Dryas_RagTag_Old.chain ./CHH/Alaska_CHH_W_C_CHH.bedGraph ./CHH/RagTag_Alaska_CHH_W_C_CHH.bedGraph
CrossMap bed ../Dryas_RagTag_Old.chain ./CHH/Svalbard_CHH_W_C_CHH.bedGraph ./CHH/RagTag_Svalbard_CHH_W_C_CHH.bedGraph
CrossMap bed ../Dryas_RagTag_Old.chain ./CHH/Nunavut_CHH_W_C_CHH.bedGraph ./CHH/RagTag_Nunavut_CHH_W_C_CHH.bedGraph
CrossMap bed ../Dryas_RagTag_Old.chain ./CHH/Sweden_CHH_W_C_CHH.bedGraph ./CHH/RagTag_Sweden_CHH_W_C_CHH.bedGraph

CrossMap bed ../Dryas_RagTag_Old.chain ./CHH/ALEX_SVAL_ALAS_Sites_intersect_DMRs_CHH.bedGraph ./CHH/RagTag_ALEX_SVAL_ALAS_Sites_intersect_DMRs_CHH.bedGraph
CrossMap bed ../Dryas_RagTag_Old.chain ./CHH/ALEX_SVAL_SWED_Sites_intersect_DMRs_CHH.bedGraph ./CHH/RagTag_ALEX_SVAL_SWED_Sites_intersect_DMRs_CHH.bedGraph
CrossMap bed ../Dryas_RagTag_Old.chain ./CHH/ALEX_SWED_ALAS_Sites_intersect_DMRs_CHH.bedGraph ./CHH/RagTag_ALEX_SWED_ALAS_Sites_intersect_DMRs_CHH.bedGraph
CrossMap bed ../Dryas_RagTag_Old.chain ./CHH/ALL_Sites_intersect_DMRs_CHH.bedGraph ./CHH/RagTag_ALL_Sites_intersect_DMRs_CHH.bedGraph
CrossMap bed ../Dryas_RagTag_Old.chain ./CHH/SVAL_SWED_ALAS_Sites_intersect_DMRs_CHH.bedGraph ./CHH/RagTag_SVAL_SWED_ALAS_Sites_intersect_DMRs_CHH.bedGraph

# SNPs
CrossMap vcf ../Dryas_RagTag_Old.chain Dryas_filtered_biasSNPs.vcf ../OldRagtagDryOcto.fasta RagTag_Dryas_filtered_biasSNPs.vcf
# 2024-12-07 09:36:22 [INFO]  Total entries: 77866
# 2024-12-07 09:36:22 [INFO]  Failed to map: 77866
# SNPs called on the old Newest Ref!


#---------------------------------------------
# Map to new reference
# genes
CrossMap gff ../Dryas_New_Old.chain Dryas_octopetala_H1.gff3 CrossMap_Dryas_genes.gff3

# TEs
CrossMap gff ../Dryas_New_Old.chain Dryas_octopetala_H1.supercontigs.fa.mod.EDTA.TEanno.gff3 CrossMap_Dryas_TEs.gff3

# Bedgraphs per site_treatment
CrossMap bed ../Dryas_New_Old.chain C_ALAS.bedGraph CrossMap_C_ALAS.bedGraph
CrossMap bed ../Dryas_New_Old.chain C_CASS.bedGraph CrossMap_C_CASS.bedGraph
CrossMap bed ../Dryas_New_Old.chain C_DRY.bedGraph CrossMap_C_DRY.bedGraph
CrossMap bed ../Dryas_New_Old.chain C_FERT.bedGraph CrossMap_C_FERT.bedGraph
CrossMap bed ../Dryas_New_Old.chain C_MEAD.bedGraph CrossMap_C_MEAD.bedGraph

CrossMap bed ../Dryas_New_Old.chain C_SVAL.bedGraph CrossMap_C_SVAL.bedGraph
CrossMap bed ../Dryas_New_Old.chain C_SwedC.bedGraph CrossMap_C_SwedC.bedGraph
CrossMap bed ../Dryas_New_Old.chain C_WILL.bedGraph CrossMap_C_WILL.bedGraph

CrossMap bed ../Dryas_New_Old.chain W_ALAS.bedGraph CrossMap_W_ALAS.bedGraph
CrossMap bed ../Dryas_New_Old.chain W_CASS.bedGraph CrossMap_W_CASS.bedGraph
CrossMap bed ../Dryas_New_Old.chain W_DRY.bedGraph CrossMap_W_DRY.bedGraph
CrossMap bed ../Dryas_New_Old.chain W_FERT.bedGraph CrossMap_W_FERT.bedGraph
CrossMap bed ../Dryas_New_Old.chain W_MEAD.bedGraph CrossMap_W_MEAD.bedGraph
CrossMap bed ../Dryas_New_Old.chain W_SVAL.bedGraph CrossMap_W_SVAL.bedGraph
CrossMap bed ../Dryas_New_Old.chain W_SwedC.bedGraph CrossMap_W_SwedC.bedGraph
CrossMap bed ../Dryas_New_Old.chain W_WILL.bedGraph CrossMap_W_WILL.bedGraph

# DMRs (W/C, Low/High ...)
CrossMap bed ../Dryas_New_Old.chain C_ALAS.bedGraph CrossMap_C_ALAS.bedGraph

# SNPs
CrossMap vcf ../Dryas_New_Old.chain <vcffile> <outfile>































###############################
# Other attempts at chain files

# needs chain file for Old genome to new Ragtag 
# https://www.biostars.org/p/391080/

# create .chain file
module load StdEnv/2020 blat/3.5 kentutils/401

ID1=Dryas_octopetala_H1.supercontigs
ID2=OldRagtagDryOcto

cd $ID1
faToTwoBit $ID1.fa $ID1.2bit
twoBitInfo $ID1.2bit chrom.sizes
cd ..
cd $ID2
faToTwoBit $ID2.fasta $ID2.2bit
twoBitInfo $ID2.2bit chrom.sizes
cd ..

blat $ID1/$ID1.2bit $ID2/$ID2.fasta $ID1\to$ID2.psl -tileSize=12 -minScore=100 -minIdentity=98 
axtChain -linearGap=medium -psl $ID1\to$ID2.psl $ID1/$ID1.2bit $ID2/$ID2.2bit $ID1\to$ID2.chain

#--------------------------
# taking too long

# try
git clone https://github.com/icebert/pblat.git
make


cat << EOF > CrossMap_make_chain.sh
#!/bin/bash
#SBATCH --account=def-cronk
#SBATCH --time=0-5:00:00
#SBATCH --cpus-per-task=40
#SBATCH --mem=191000M

module load StdEnv/2020 blat/3.5 kentutils/401

cd /home/celphin/scratch/Dryas/CrossMap/

ID1=Dryas_octopetala_H1.supercontigs
ID2=OldRagtagDryOcto

/home/celphin/scratch/Dryas/CrossMap/pblat/pblat -threads=40 $ID1/$ID1.fa $ID2/$ID2.fasta NewtoOld.psl 
axtChain -linearGap=medium -psl NewtoOld.psl $ID1/$ID1.2bit $ID2/$ID2.2bit NewtoOld.chain
EOF

sbatch CrossMap_make_chain.sh

#-------------------------
# Too slow too
# try splitting chain file??
# https://genomewiki.ucsc.edu/index.php/Minimal_Steps_For_LiftOver
# try https://genomewiki.ucsc.edu/index.php/DoBlastzChainNet.pl












