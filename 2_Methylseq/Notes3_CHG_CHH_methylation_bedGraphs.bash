#########################################################
# Dryas DNA methylation calling
# nf-core Methylseq
# Convert Bismark text to bedgraph bismark2bedGraph
# https://github.com/FelixKrueger/Bismark
# Nov 2024
#############################################################


cd /home/celphin/projects/rpp-rieseber/celphin/Dryas/Methylation_calling/May2021_Parents/bismark_methylation_calls

cp ./methylation_calls/* ~/scratch/Dryas/CHG_CHH/

cd ~/scratch/Dryas/CHG_CHH/

gunzip *

#-------------------------
# Install Bismark

git clone https://github.com/FelixKrueger/Bismark.git

# Needs 
 # minimap2
 # Samtools
 # perl

module load 




#############################
# Run 
# https://felixkrueger.github.io/Bismark/bismark/methylation_extraction/

cd ~/scratch/Dryas/CHG_CHH/

bismark2bedGraph --CX Non_CpG_context_W4.G10.W2d1_DRY9W_50_185_R1_val_1_bismark_bt2_pe.deduplicated.txt

# Note try to run with new ref genome version
