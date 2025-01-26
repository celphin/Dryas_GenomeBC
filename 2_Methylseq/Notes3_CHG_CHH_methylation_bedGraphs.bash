#########################################################
# Dryas DNA CHG CHH methylation calling
# nf-core Methylseq
# Convert Bismark text to bedgraph bismark2bedGraph
# https://github.com/FelixKrueger/Bismark
# Nov 2024
#############################################################

# Beluga

cd /home/celphin/projects/rpp-rieseber/celphin/Dryas/Methylation_calling/May2021_Parents/bismark_methylation_calls

cp ./methylation_calls/Non_CpG_context* ~/scratch/Dryas/CHG_CHH/

cd ~/scratch/Dryas/CHG_CHH/data

gunzip *

#-------------------------
# Install Bismark

git clone https://github.com/FelixKrueger/Bismark.git

#############################
# Run for one sample
# https://felixkrueger.github.io/Bismark/bismark/methylation_extraction/

tmux new-session -s CHH
tmux attach-session -t CHH

cd ~/scratch/Dryas/CHG_CHH/
module load StdEnv/2023 minimap2/2.28 samtools/1.20 perl/5.36.1

/home/celphin/scratch/Dryas/CHG_CHH/Bismark/bismark2bedGraph --CX Non_CpG_context_C1.A10.C1d12_CASS5C_529_159_R1_val_1_bismark_bt2_pe.deduplicated.txt -o CASS5C_529_159_Non_CpG.bedGraph

# takes about 3 hour per sample

# Note try to run with new ref genome version in future

###################################
# Run for all Wild samples

# Change to the directory where your files are located
cd /home/celphin/scratch/Dryas/CHG_CHH/data

FILES=$(ls *bismark_bt2_pe.deduplicated.txt)
echo ${FILES}

cd /home/celphin/scratch/Dryas/CHG_CHH/
for file in ${FILES}
do
output_name=$(echo ${file} | awk -F 'Non_CpG_context_|_R1_val_1_bismark_bt2_pe.deduplicated.txt' '{print $2}' | sed 's/_/./g')

cat << EOF > CHH_${output_name}.sh
#!/bin/bash
#SBATCH --account=def-cronk
#SBATCH --time=0-5:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=120000M

module load StdEnv/2023 minimap2/2.28 samtools/1.20 perl/5.36.1

cd /home/celphin/scratch/Dryas/CHG_CHH/data
/home/celphin/scratch/Dryas/CHG_CHH/Bismark/bismark2bedGraph --CX "$file" -o "$output_name"
EOF

sbatch CHH_${output_name}.sh

done

sq
          JOBID     USER      ACCOUNT           NAME  ST  TIME_LEFT NODES CPUS TRES_PER_N MIN_MEM NODELIST (REASON)
       52914166  celphin def-cronk_cp CHH_C1.A10.C1d   R    4:59:37     1    1        N/A 120000M bl12402 (None)
       52914167  celphin def-cronk_cp CHH_C1.B05.C1d   R    4:59:37     1    1        N/A 120000M bl12402 (None)
       52914168  celphin def-cronk_cp CHH_C1.B07.C1b   R    4:59:37     1    1        N/A 120000M bl12437 (None)
       52914169  celphin def-cronk_cp CHH_C1.C06.C1d   R    4:59:37     1    1        N/A 120000M bl12437 (None)
       52914170  celphin def-cronk_cp CHH_C1.C07.C1c   R    4:59:37     1    1        N/A 120000M bl12438 (None)
       52914171  celphin def-cronk_cp CHH_C1.C09.C1d   R    4:59:37     1    1        N/A 120000M bl12438 (None)
       52914172  celphin def-cronk_cp CHH_C1.C10.C1f  PD    5:00:00     1    1        N/A 120000M  (Priority)
       52914173  celphin def-cronk_cp CHH_C1.D05.C1a  PD    5:00:00     1    1        N/A 120000M  (Priority)
       52914174  celphin def-cronk_cp CHH_C1.D06.C1h  PD    5:00:00     1    1        N/A 120000M  (Priority)
       52914175  celphin def-cronk_cp CHH_C1.D07.C1f  PD    5:00:00     1    1        N/A 120000M  (Priority)
       52914176  celphin def-cronk_cp CHH_C1.D08.C1h  PD    5:00:00     1    1        N/A 120000M  (Priority)
       52914177  celphin def-cronk_cp CHH_C1.E06.C1b  PD    5:00:00     1    1        N/A 120000M  (Priority)
       52914178  celphin def-cronk_cp CHH_C1.E09.C1g  PD    5:00:00     1    1        N/A 120000M  (Priority)
       52914179  celphin def-cronk_cp CHH_C1.F05.C1e  PD    5:00:00     1    1        N/A 120000M  (Priority)
       52914180  celphin def-cronk_cp CHH_C1.F07.C1f  PD    5:00:00     1    1        N/A 120000M  (Priority)
       52914181  celphin def-cronk_cp CHH_C1.F08.C1a  PD    5:00:00     1    1        N/A 120000M  (Priority)
       52914182  celphin def-cronk_cp CHH_C1.F09.C1a  PD    5:00:00     1    1        N/A 120000M  (Priority)
       52914183  celphin def-cronk_cp CHH_C1.G05.C1d  PD    5:00:00     1    1        N/A 120000M  (Priority)
       52914184  celphin def-cronk_cp CHH_C1.G07.C1h  PD    5:00:00     1    1        N/A 120000M  (Priority)
       52914185  celphin def-cronk_cp CHH_C1.H05.C1f  PD    5:00:00     1    1        N/A 120000M  (Priority)
       52914186  celphin def-cronk_cp CHH_C1.H07.C1b  PD    5:00:00     1    1        N/A 120000M  (Priority)
       52914187  celphin def-cronk_cp CHH_C2.19.1a1.  PD    5:00:00     1    1        N/A 120000M  (Priority)
       52914188  celphin def-cronk_cp CHH_C2.21.1f4.  PD    5:00:00     1    1        N/A 120000M  (Priority)
       52914189  celphin def-cronk_cp CHH_C2.2.3a7.A  PD    5:00:00     1    1        N/A 120000M  (Priority)
       52914190  celphin def-cronk_cp CHH_C2.25.3d2.  PD    5:00:00     1    1        N/A 120000M  (Priority)
       52914191  celphin def-cronk_cp CHH_C2.27.3e1.  PD    5:00:00     1    1        N/A 120000M  (Priority)
       52914192  celphin def-cronk_cp CHH_C2.4.3e6.A  PD    5:00:00     1    1        N/A 120000M  (Priority)
       52914193  celphin def-cronk_cp CHH_C2.6.3g11.  PD    5:00:00     1    1        N/A 120000M  (Priority)
       52914194  celphin def-cronk_cp CHH_C2.9.3f11.  PD    5:00:00     1    1        N/A 120000M  (Priority)
       52914195  celphin def-cronk_cp CHH_C2.A03.C2e  PD    5:00:00     1    1        N/A 120000M  (Priority)
       52914196  celphin def-cronk_cp CHH_C2.A11.C2a  PD    5:00:00     1    1        N/A 120000M  (Priority)
       52914198  celphin def-cronk_cp CHH_C2.A12.C2c  PD    5:00:00     1    1        N/A 120000M  (Priority)
       52914199  celphin def-cronk_cp CHH_C2.B02.C2h  PD    5:00:00     1    1        N/A 120000M  (Priority)
       52914200  celphin def-cronk_cp CHH_C2.B03.C2f  PD    5:00:00     1    1        N/A 120000M  (Priority)
       52914201  celphin def-cronk_cp CHH_C2.B11.C2f  PD    5:00:00     1    1        N/A 120000M  (Priority)
       52914202  celphin def-cronk_cp CHH_C2.C01.C2h  PD    5:00:00     1    1        N/A 120000M  (Priority)
       52914203  celphin def-cronk_cp CHH_C2.C02.C2b  PD    5:00:00     1    1        N/A 120000M  (Priority)
       52914204  celphin def-cronk_cp CHH_C2.C11.C2h  PD    5:00:00     1    1        N/A 120000M  (Priority)
       52914205  celphin def-cronk_cp CHH_C2.D02.C2d  PD    5:00:00     1    1        N/A 120000M  (Priority)
       52914206  celphin def-cronk_cp CHH_C2.D03.C2h  PD    5:00:00     1    1        N/A 120000M  (Priority)
       52914207  celphin def-cronk_cp CHH_C2.D11.C2a  PD    5:00:00     1    1        N/A 120000M  (Priority)
       52914208  celphin def-cronk_cp CHH_C2.E01.C2b  PD    5:00:00     1    1        N/A 120000M  (Priority)
       52914209  celphin def-cronk_cp CHH_C2.E02.C2g  PD    5:00:00     1    1        N/A 120000M  (Priority)
       52914210  celphin def-cronk_cp CHH_C2.E10.C1h  PD    5:00:00     1    1        N/A 120000M  (Priority)
       52914211  celphin def-cronk_cp CHH_C2.E11.C2c  PD    5:00:00     1    1        N/A 120000M  (Priority)
       52914212  celphin def-cronk_cp CHH_C2.F10.C1d  PD    5:00:00     1    1        N/A 120000M  (Priority)
       52914213  celphin def-cronk_cp CHH_C2.F11.C2e  PD    5:00:00     1    1        N/A 120000M  (Priority)
       52914214  celphin def-cronk_cp CHH_C2.G01.C2h  PD    5:00:00     1    1        N/A 120000M  (Priority)
       52914215  celphin def-cronk_cp CHH_C2.G04.C2c  PD    5:00:00     1    1        N/A 120000M  (Priority)
       52914216  celphin def-cronk_cp CHH_C2.G12.C2a  PD    5:00:00     1    1        N/A 120000M  (Priority)
       52914217  celphin def-cronk_cp CHH_C2.H12.C2g  PD    5:00:00     1    1        N/A 120000M  (Priority)
       52914218  celphin def-cronk_cp CHH_Lab.sample  PD    5:00:00     1    1        N/A 120000M  (Priority)
       52914219  celphin def-cronk_cp CHH_W1.A06.W1b  PD    5:00:00     1    1        N/A 120000M  (Priority)
       52914220  celphin def-cronk_cp CHH_W1.A07.W1a  PD    5:00:00     1    1        N/A 120000M  (Priority)
       52914221  celphin def-cronk_cp CHH_W1.A08.W1e  PD    5:00:00     1    1        N/A 120000M  (Priority)
       52914222  celphin def-cronk_cp CHH_W1.A09.W1f  PD    5:00:00     1    1        N/A 120000M  (Priority)
       52914223  celphin def-cronk_cp CHH_W1.B06.W1c  PD    5:00:00     1    1        N/A 120000M  (Priority)
       52914224  celphin def-cronk_cp CHH_W1.B08.W1f  PD    5:00:00     1    1        N/A 120000M  (Priority)
       52914225  celphin def-cronk_cp CHH_W1.B09.W1b  PD    5:00:00     1    1        N/A 120000M  (Priority)
       52914226  celphin def-cronk_cp CHH_W1.B10.W1e  PD    5:00:00     1    1        N/A 120000M  (Priority)
       52914227  celphin def-cronk_cp CHH_W1.C05.W1g  PD    5:00:00     1    1        N/A 120000M  (Priority)
       52914228  celphin def-cronk_cp CHH_W1.C08.W1g  PD    5:00:00     1    1        N/A 120000M  (Priority)
       52914229  celphin def-cronk_cp CHH_W1.D09.W1f  PD    5:00:00     1    1        N/A 120000M  (Priority)
       52914230  celphin def-cronk_cp CHH_W1.D10.W1g  PD    5:00:00     1    1        N/A 120000M  (Priority)
       52914231  celphin def-cronk_cp CHH_W1.E05.W1c  PD    5:00:00     1    1        N/A 120000M  (Priority)
       52914232  celphin def-cronk_cp CHH_W1.E07.W1g  PD    5:00:00     1    1        N/A 120000M  (Priority)
       52914233  celphin def-cronk_cp CHH_W1.E08.W1c  PD    5:00:00     1    1        N/A 120000M  (Priority)
       52914234  celphin def-cronk_cp CHH_W1.F06.W1e  PD    5:00:00     1    1        N/A 120000M  (Priority)
       52914235  celphin def-cronk_cp CHH_W1.G06.W1g  PD    5:00:00     1    1        N/A 120000M  (Priority)
       52914236  celphin def-cronk_cp CHH_W1.G08.W1b  PD    5:00:00     1    1        N/A 120000M  (Priority)
       52914237  celphin def-cronk_cp CHH_W1.G09.W1b  PD    5:00:00     1    1        N/A 120000M  (Priority)
       52914238  celphin def-cronk_cp CHH_W1.H06.W1h  PD    5:00:00     1    1        N/A 120000M  (Priority)
       52914239  celphin def-cronk_cp CHH_W1.H08.W1d  PD    5:00:00     1    1        N/A 120000M  (Priority)
       52914240  celphin def-cronk_cp CHH_W1.H09.W1c  PD    5:00:00     1    1        N/A 120000M  (Priority)
       52914241  celphin def-cronk_cp CHH_W2.1.3f1.L  PD    5:00:00     1    1        N/A 120000M  (Priority)
       52914242  celphin def-cronk_cp CHH_W2.20.1b3.  PD    5:00:00     1    1        N/A 120000M  (Priority)
       52914243  celphin def-cronk_cp CHH_W2.22.1e1.  PD    5:00:00     1    1        N/A 120000M  (Priority)
       52914244  celphin def-cronk_cp CHH_W2.23.3c2.  PD    5:00:00     1    1        N/A 120000M  (Priority)
       52914245  celphin def-cronk_cp CHH_W2.3.3b7.A  PD    5:00:00     1    1        N/A 120000M  (Priority)
       52914246  celphin def-cronk_cp CHH_W2.5.3f6.A  PD    5:00:00     1    1        N/A 120000M  (Priority)
       52914247  celphin def-cronk_cp CHH_W2.7.3h11.  PD    5:00:00     1    1        N/A 120000M  (Priority)
       52914248  celphin def-cronk_cp CHH_W2.8.3e11.  PD    5:00:00     1    1        N/A 120000M  (Priority)
       52914249  celphin def-cronk_cp CHH_W2.A01.W2e  PD    5:00:00     1    1        N/A 120000M  (Priority)
       52914250  celphin def-cronk_cp CHH_W2.A02.W2c  PD    5:00:00     1    1        N/A 120000M  (Priority)
       52914251  celphin def-cronk_cp CHH_W2.B01.W2f  PD    5:00:00     1    1        N/A 120000M  (Priority)
       52914252  celphin def-cronk_cp CHH_W2.B12.W2e  PD    5:00:00     1    1        N/A 120000M  (Priority)
       52914253  celphin def-cronk_cp CHH_W2.C03.W2g  PD    5:00:00     1    1        N/A 120000M  (Priority)
       52914254  celphin def-cronk_cp CHH_W2.C12.W2f  PD    5:00:00     1    1        N/A 120000M  (Priority)
       52914255  celphin def-cronk_cp CHH_W2.D01.W2a  PD    5:00:00     1    1        N/A 120000M  (Priority)
       52914256  celphin def-cronk_cp CHH_W2.D12.W2b  PD    5:00:00     1    1        N/A 120000M  (Priority)
       52914257  celphin def-cronk_cp CHH_W2.E12.W2e  PD    5:00:00     1    1        N/A 120000M  (Priority)
       52914258  celphin def-cronk_cp CHH_W2.F01.W2d  PD    5:00:00     1    1        N/A 120000M  (Priority)
       52914259  celphin def-cronk_cp CHH_W2.F02.W2a  PD    5:00:00     1    1        N/A 120000M  (Priority)
       52914261  celphin def-cronk_cp CHH_W2.F04.W2b  PD    5:00:00     1    1        N/A 120000M  (Priority)
       52914262  celphin def-cronk_cp CHH_W2.F12.W2f  PD    5:00:00     1    1        N/A 120000M  (Priority)
       52914263  celphin def-cronk_cp CHH_W2.G02.W2c  PD    5:00:00     1    1        N/A 120000M  (Priority)
       52914264  celphin def-cronk_cp CHH_W2.G11.W2h  PD    5:00:00     1    1        N/A 120000M  (Priority)
       52914265  celphin def-cronk_cp CHH_W2.H01.W2b  PD    5:00:00     1    1        N/A 120000M  (Priority)
       52914266  celphin def-cronk_cp CHH_W2.H02.W2a  PD    5:00:00     1    1        N/A 120000M  (Priority)
       52914267  celphin def-cronk_cp CHH_W2.H04.W2d  PD    5:00:00     1    1        N/A 120000M  (Priority)
       52914268  celphin def-cronk_cp CHH_W2.H10.W2h  PD    5:00:00     1    1        N/A 120000M  (Priority)
       52914269  celphin def-cronk_cp CHH_W2.H11.W2a  PD    5:00:00     1    1        N/A 120000M  (Priority)
       52914270  celphin def-cronk_cp CHH_W4.G10.W2d  PD    5:00:00     1    1        N/A 120000M  (Priority)


mkdir coverage
mv *.bismark.cov.gz coverage/
mkdir bedgraph
mv *.gz bedgraph/

cd bedgraph
rename .gz .bedGraph.gz *
gunzip *bedGraph.gz


###################################
# Run for all Seedling and Phenology samples

tmux new-session -s CHH
tmux attach-session -t CHH

module load StdEnv/2023 minimap2/2.28 samtools/1.20 perl/5.36.1

# Change to the directory where your files are located
cd /home/celphin/scratch/Dryas/CHG_CHH/seedling_data

#gunzip *.gz

FILES=$(ls *bismark_bt2_pe.deduplicated.txt)
echo ${FILES}

for file in ${FILES}
do
output_name=$(echo ${file} | awk -F 'Non_CpG_context_|_R1_val_1_bismark_bt2_pe.deduplicated.txt' '{print $2}' | sed 's/_/./g')

echo $file
echo $output_name

cat << EOF > CHH_${output_name}.sh
#!/bin/bash
#SBATCH --account=def-cronk
#SBATCH --time=0-7:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=200000M

module load StdEnv/2023 minimap2/2.28 samtools/1.20 perl/5.36.1

cd /home/celphin/scratch/Dryas/CHG_CHH/seedling_data
/home/celphin/scratch/Dryas/CHG_CHH/Bismark/bismark2bedGraph --CX "$file" -o "$output_name"
EOF

sbatch CHH_${output_name}.sh

done

mv *bismark.cov.gz coverage/
rename .gz .bedGraph.gz *.gz
mv *.gz ../bedgraph

cd ../bedgraph

gunzip *.gz

mkdir seedling
mkdir Phenology

mv SE* seedling/
mv MatFl.* Phenology/
cd Phenology/

# check files worked
ls -la

# MatFl.Cass.10W.60.544.F112581  MatFl.Fert.5C.97.1F.F112583   MatFl.Mead.1W.116.444.F112578
# MatFl.Cass.4C.4.524.F112579    MatFl.Fert.6W.110.3F.F112584  MatFl.Will.3C.100.414.F112575
# MatFl.Cass.5W.130.525.F112580  MatFl.Mead.1C.33.446.F112577  MatFl.Will.4W.13.417.F112576

cp W1.F06.W1e5.CASS10W.544.60_CHH.bedGraph ./Phenology/Sen.CASS10W.544.60_CHH.bedGraph
cp C1.B05.C1d1.CASS4C.524.4_CHH.bedGraph ./Phenology/Sen.CASS4C.524.4_CHH.bedGraph
cp W1.A09.W1f10.CASS5W.525.130_CHH.bedGraph ./Phenology/Sen.CASS5W.525.130_CHH.bedGraph
cp C1.G07.C1h7.FERT5C.1F.97_CHH.bedGraph ./Phenology/Sen.FERT5C.1F.97_CHH.bedGraph
cp W1.B08.W1f8.FERT6W.3F.110_CHH.bedGraph ./Phenology/Sen.FERT6W.3F.110_CHH.bedGraph
cp C1.H05.C1f3.MEAD1C.446.33_CHH.bedGraph ./Phenology/Sen.MEAD1C.446.33_CHH.bedGraph
cp W1.E08.W1c9.MEAD1W.444.116_CHH.bedGraph ./Phenology/Sen.MEAD1W.444.116_CHH.bedGraph
cp C1.H07.C1b8.WILL3C.414.100_CHH.bedGraph ./Phenology/Sen.WILL3C.414.100_CHH.bedGraph
cp W1.C05.W1g1.WILL4W.417.13_CHH.bedGraph ./Phenology/Sen.WILL4W.417.13_CHH.bedGraph

#############################
# Make cytosine reports from each CHG/CHH coverage file

# on narval1

#-----------------
# to make stranded cytosine reports for CHH
cd /lustre04/scratch/celphin/Dryas_large_folders/CHG_CHH/data/coverage

#############################
# Run for one sample
# https://felixkrueger.github.io/Bismark/bismark/methylation_extraction/

tmux new-session -s CHH
tmux attach-session -t CHH

cd /lustre04/scratch/celphin/Dryas_large_folders/CHG_CHH/data/coverage
module load StdEnv/2023 minimap2/2.28 samtools/1.20 perl/5.36.1

/lustre04/scratch/celphin/Dryas_large_folders/CHG_CHH/Bismark/coverage2cytosine \
--CX C1.A10.C1d12.CASS5C.529.159.gz.bismark.cov \
--genome_folder /lustre04/scratch/celphin/Dryas_large_folders/BismarkIndex \
-o CASS5C_529_159_Non_CpG.report.txt

# works

# Note try to run with new ref genome version in future

#--------------------------------
# Run for all Wild samples

# Change to the directory where your files are located
cd /lustre04/scratch/celphin/Dryas_large_folders/CHG_CHH/data/coverage

FILES=$(ls *.gz.bismark.cov)
echo ${FILES}

cd /lustre04/scratch/celphin/Dryas_large_folders/CHG_CHH/data/coverage
for file in ${FILES}
do
output_name=$(echo ${file} | sed 's/\.gz\.bismark\.cov//g')
echo "$output_name"


cat << EOF > CHH_report_${output_name}.sh
#!/bin/bash
#SBATCH --account=def-cronk
#SBATCH --time=0-5:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=120000M

module load StdEnv/2023 minimap2/2.28 samtools/1.20 perl/5.36.1

cd /lustre04/scratch/celphin/Dryas_large_folders/CHG_CHH/data/coverage

/lustre04/scratch/celphin/Dryas_large_folders/CHG_CHH/Bismark/coverage2cytosine \
--CX "$file" \
--genome_folder /lustre04/scratch/celphin/Dryas_large_folders/BismarkIndex \
-o ${output_name}_Non_CpG.report.txt


EOF

sbatch CHH_report_${output_name}.sh

done

mkdir CHH_cytosine_report
mv *CX_report.txt CHH_cytosine_report


##################################

# Run for all Seedling and Phenology samples

# Change to the directory where your files are located
cd /lustre04/scratch/celphin/Dryas_large_folders/CHG_CHH/seedling_data/coverage

FILES=$(ls *.gz.bismark.cov)
echo ${FILES}

for file in ${FILES}
do
output_name=$(echo ${file} | sed 's/\.gz\.bismark\.cov//g')

cat << EOF > CHH_report_${output_name}.sh
#!/bin/bash
#SBATCH --account=def-cronk
#SBATCH --time=0-5:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=120000M

module load StdEnv/2023 minimap2/2.28 samtools/1.20 perl/5.36.1

cd /lustre04/scratch/celphin/Dryas_large_folders/CHG_CHH/seedling_data/coverage

/lustre04/scratch/celphin/Dryas_large_folders/CHG_CHH/Bismark/coverage2cytosine \
--CX "$file" \
--genome_folder /lustre04/scratch/celphin/Dryas_large_folders/BismarkIndex \
-o ${output_name}_Non_CpG.report.txt


EOF

sbatch CHH_report_${output_name}.sh

done

mkdir CHH_cytosine_report
mv *CX_report.txt CHH_cytosine_report

#-------------------------
# For seedlings

# Change to the directory where your files are located
cd /lustre04/scratch/celphin/Dryas_large_folders/CHG_CHH/seedling_data/coverage/seedling

FILES=$(ls *.gz.bismark.cov)
echo ${FILES}

for file in ${FILES}
do
output_name=$(echo ${file} | sed 's/\.gz\.bismark\.cov//g')

cat << EOF > CHH_report_${output_name}.sh
#!/bin/bash
#SBATCH --account=def-cronk
#SBATCH --time=0-5:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=120000M

module load StdEnv/2023 minimap2/2.28 samtools/1.20 perl/5.36.1

cd /lustre04/scratch/celphin/Dryas_large_folders/CHG_CHH/seedling_data/coverage/seedling

/lustre04/scratch/celphin/Dryas_large_folders/CHG_CHH/Bismark/coverage2cytosine \
--CX "$file" \
--genome_folder /lustre04/scratch/celphin/Dryas_large_folders/BismarkIndex \
-o ${output_name}_Non_CpG.report.txt


EOF

sbatch CHH_report_${output_name}.sh

done

mkdir CHH_cytosine_report
mv *CX_report.txt CHH_cytosine_report



##################################
# Make average bedGraph file for each site_treatment

tmux attach-session -t CHH

cd /home/celphin/scratch/Dryas/CHG_CHH/data/bedgraph
module load StdEnv/2023 bedtools/2.31.0 

bedtools unionbedg -i W*SVAL*.bedGraph > W_SVAL_total.bedGraph
awk 'OFS="\t" {sum=0; for (col=4; col<=NF; col++) sum += $col; print $1, $2, $3, sum/(NF-4+1); }' W_SVAL_total.bedGraph > W_SVAL.bedGraph

bedtools unionbedg -i C*SVAL*.bedGraph > C_SVAL_total.bedGraph
awk 'OFS="\t" {sum=0; for (col=4; col<=NF; col++) sum += $col; print $1, $2, $3, sum/(NF-4+1); }' C_SVAL_total.bedGraph > C_SVAL.bedGraph

#---------------
bedtools unionbedg -i W*LAT*.bedGraph > W_Swed_total.bedGraph
awk 'OFS="\t" {sum=0; for (col=4; col<=NF; col++) sum += $col; print $1, $2, $3, sum/(NF-4+1); }' W_Swed_total.bedGraph > W_SwedC.bedGraph

bedtools unionbedg -i C*LAT*.bedGraph > C_Swed_total.bedGraph
awk 'OFS="\t" {sum=0; for (col=4; col<=NF; col++) sum += $col; print $1, $2, $3, sum/(NF-4+1); }' C_Swed_total.bedGraph > C_SwedC.bedGraph

#---------------
bedtools unionbedg -i W*CASS*.bedGraph > W_CASS_total.bedGraph
awk 'OFS="\t" {sum=0; for (col=4; col<=NF; col++) sum += $col; print $1, $2, $3, sum/(NF-4+1); }' W_CASS_total.bedGraph > W_CASS.bedGraph

bedtools unionbedg -i C*CASS*.bedGraph  > C_CASS_total.bedGraph
awk 'OFS="\t" {sum=0; for (col=4; col<=NF; col++) sum += $col; print $1, $2, $3, sum/(NF-4+1); }' C_CASS_total.bedGraph > C_CASS.bedGraph

#---------------
bedtools unionbedg -i W*WILL*.bedGraph > W_WILL_total.bedGraph
awk 'OFS="\t" {sum=0; for (col=4; col<=NF; col++) sum += $col; print $1, $2, $3, sum/(NF-4+1); }' W_WILL_total.bedGraph > W_WILL.bedGraph

bedtools unionbedg -i C*WILL*.bedGraph > C_WILL_total.bedGraph
awk 'OFS="\t" {sum=0; for (col=4; col<=NF; col++) sum += $col; print $1, $2, $3, sum/(NF-4+1); }' C_WILL_total.bedGraph > C_WILL.bedGraph

#---------------
bedtools unionbedg -i W*FERT*.bedGraph > W_FERT_total.bedGraph
awk 'OFS="\t" {sum=0; for (col=4; col<=NF; col++) sum += $col; print $1, $2, $3, sum/(NF-4+1); }' W_FERT_total.bedGraph > W_FERT.bedGraph

bedtools unionbedg -i C*FERT*.bedGraph > C_FERT_total.bedGraph
awk 'OFS="\t" {sum=0; for (col=4; col<=NF; col++) sum += $col; print $1, $2, $3, sum/(NF-4+1); }' C_FERT_total.bedGraph > C_FERT.bedGraph

#---------------
bedtools unionbedg -i W*DRY*.bedGraph > W_DRY_total.bedGraph
awk 'OFS="\t" {sum=0; for (col=4; col<=NF; col++) sum += $col; print $1, $2, $3, sum/(NF-4+1); }' W_DRY_total.bedGraph > W_DRY.bedGraph

bedtools unionbedg -i C*DRY*.bedGraph > C_DRY_total.bedGraph
awk 'OFS="\t" {sum=0; for (col=4; col<=NF; col++) sum += $col; print $1, $2, $3, sum/(NF-4+1); }' C_DRY_total.bedGraph > C_DRY.bedGraph

#---------------
bedtools unionbedg -i W*MEAD*.bedGraph > W_MEAD_total.bedGraph
awk 'OFS="\t" {sum=0; for (col=4; col<=NF; col++) sum += $col; print $1, $2, $3, sum/(NF-4+1); }' W_MEAD_total.bedGraph > W_MEAD.bedGraph

bedtools unionbedg -i C*MEAD*.bedGraph > C_MEAD_total.bedGraph
awk 'OFS="\t" {sum=0; for (col=4; col<=NF; col++) sum += $col; print $1, $2, $3, sum/(NF-4+1); }' C_MEAD_total.bedGraph > C_MEAD.bedGraph

#---------------
bedtools unionbedg -i W*ALAS*.bedGraph  > W_ALAS_total.bedGraph
awk 'OFS="\t" {sum=0; for (col=4; col<=NF; col++) sum += $col; print $1, $2, $3, sum/(NF-4+1); }' W_ALAS_total.bedGraph > W_ALAS.bedGraph

bedtools unionbedg -i C*ALAS*.bedGraph > C_ALAS_total.bedGraph
awk 'OFS="\t" {sum=0; for (col=4; col<=NF; col++) sum += $col; print $1, $2, $3, sum/(NF-4+1); }' C_ALAS_total.bedGraph > C_ALAS.bedGraph

# rename to add CHH to all

rename .bedGraph _CHH.bedGraph *

###############################
