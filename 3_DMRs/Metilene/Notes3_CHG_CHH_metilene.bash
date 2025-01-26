#########################################################
# Dryas DNA methylation calling
# nf-core Methylseq
# Convert Bismark text to bedgraph bismark2bedGraph
# https://github.com/FelixKrueger/Bismark
# Nov 2024
#############################################################

# Run metilene on CHG and CHH contexts

# install metilene

#Installing ggplot: https://docs.alliancecan.ca/wiki/
module load StdEnv/2023
module load r/4.4.0
module load gdal
module load udunits
module load python

#Install R ggplot2 : https://docs.alliancecan.ca/wiki/R
mkdir /home/celphin/R/x86_64-pc-linux-gnu-library/4.4/
export R_LIBS_USER=/home/celphin/R/x86_64-pc-linux-gnu-library/4.4/
R -e 'install.packages("ggplot2", repos="https://cloud.r-project.org/")'

#####################################################
#Installing metilene
#In desired directory (specified in script later)
cd /home/celphin/scratch/Dryas/
wget http://www.bioinf.uni-leipzig.de/Software/metilene/metilene_v02-8.tar.gz
tar -xvzf metilene_v02-8.tar.gz
cd metilene_v0.2-8
make

############################
cd /home/celphin/scratch/Dryas/CHG_CHH/

#----------------------
#Nunavut:
mkdir Nunavut_Metilene
cd Nunavut_Metilene
cp /home/celphin/scratch/Dryas/MS_scripts/metilene*.sh .

tmux new-session -s Nunavut_DMRS
tmux attach-session -t Nunavut_DMRS
salloc -c32 --time 7:00:00 --mem 120000m --account def-cronk

#----------------------------------------
#Metilene input prep:

#Nunavut adjustments
nano metilene_prep.sh
#h1,h2: 
h1="W", h2="C"
#Input directories:
input_dir=Nunavut_${h1}_${h2}_input_files 
in_metilene="Nunavut_metilene_"$h1"_"$h2".input"
methylseq_output_dir="/home/celphin/scratch/Dryas/CHG_CHH/bedgraph"
#Copy Loop:
cp ${methylseq_output_dir}/*FERT*_CHH.bedGraph ./${input_dir}
cp ${methylseq_output_dir}/*CASS*_CHH.bedGraph ./${input_dir}
cp ${methylseq_output_dir}/*WILL*_CHH.bedGraph ./${input_dir}
cp ${methylseq_output_dir}/*DRY*_CHH.bedGraph ./${input_dir}
cp ${methylseq_output_dir}/*MEAD*_CHH.bedGraph ./${input_dir}

#Comment out rename for loops 

module load StdEnv/2020
module load bedtools/2.30.0
cd /home/celphin/scratch/Dryas/CHG_CHH/Nunavut_Metilene
sh metilene_prep.sh 

sed -i 's/Nunavut_W_C_input_files\///g' /home/celphin/scratch/Dryas/CHG_CHH/Nunavut_Metilene/Nunavut_metilene_W_C.input
sed -i 's/_CHH.bedGraph_sorted_tab.bedGraph//g' /home/celphin/scratch/Dryas/CHG_CHH/Nunavut_Metilene/Nunavut_metilene_W_C.input

#------------------------------------------------------------------

#Nunavut specific adjustments,
nano metilene_run.sh

h1='W'
h2='C'
metilene_dir=/home/celphin/scratch/Dryas
input_dir="/home/celphin/scratch/Dryas/CHG_CHH/Nunavut_Metilene"

maxdist=$1
mincpgs=$2
mindiff=$3

output_name=Nunavut_CHH_"$h1"_"$h2"_"${maxdist}"_"${mincpgs}"_"${mindiff}"
in_metilene="Nunavut_metilene_"$h1"_"$h2".input"
threads=32


#metilene:
# params: maxdist, mincpgs, mindiff
salloc -c32 --time 7:00:00 --mem 200000m --account def-cronk

module load StdEnv/2020
module load bedtools/2.30.0 
sh metilene_run.sh 70 5 4
sh metilene_run.sh 150 5 4 
sh metilene_run.sh 70 5 0.7

# slurmstepd: error: Detected 2 oom_kill events in StepId=53023285.interactive. Some of the step tasks have been OOM Killed.
# srun: error: bc11945: task 0: Out Of Memory
# salloc: Relinquishing job allocation 53023285
# salloc: Job allocation 53023285 has been revoked.


#-------------------------------------------------
#Filter based on qval
#params: maxdist, mincpgs, mindiff, minmeandif, q-value

nano  metilene_filter_qval.sh
h1='W'
h2='C'
metilene_dir=/home/celphin/scratch/Dryas
input_dir="/home/celphin/scratch/Dryas/CHG_CHH/Nunavut_Metilene"
outputname=Nunavut_CHH_"$h1"_"$h2"_"${maxdist}"_"${mincpgs}"_"${mindiff}"

# to run
module load StdEnv/2023
module load r/4.4.0
module load gdal
module load udunits
module load python
export R_LIBS_USER=/home/celphin/R/x86_64-pc-linux-gnu-library/4.4/

sh metilene_filter_qval.sh 70 5 4 0.9 0.001
#   Wrote 1841 DMRs with adj. p-value<0.001, a minimum absolute difference>=0.9, a minimum length [CpG]>=10 and a minimum length [nt]>=0
sh metilene_filter_qval.sh 150 5 4 0.9 1e-5
#  Wrote 1186 DMRs with adj. p-value<1e-5, a minimum absolute difference>=0.9, a minimum length [CpG]>=10 and a minimum length [nt]>=0
sh metilene_filter_qval.sh 70 5 0.7 0.7 0.001
#   Wrote 2588 DMRs with adj. p-value<0.001, a minimum absolute difference>=0.7, a minimum length [CpG]>=10 and a minimum length [nt]>=0

####################################################################################
#Svalbard:
mkdir Svalbard_Metilene
cd Svalbard_Metilene
cp /home/celphin/scratch/Dryas/MS_scripts/metilene*.sh .

tmux new-session -s Svalbard_DMRS
tmux attach-session -t Svalbard_DMRS
salloc -c32 --time 4:50:00 --mem 120000m --account def-rieseber

#----------------------------------------
#Metilene input prep:
module load StdEnv/2020
module load bedtools/2.30.0

#Svalbard adjustments
nano metilene_prep.sh 

#h1,h2: 
h1="W", h2="C"
#Input directories:
input_dir=Svalbard_${h1}_${h2}_input_files 
in_metilene="Svalbard_metilene_"$h1"_"$h2".input"
methylseq_output_dir="/home/celphin/scratch/Dryas/CHG_CHH/bedgraph"
#Copy Loop:
cp ${methylseq_output_dir}/*SVAL*_CHH.bedGraph ./${input_dir}
#Comment out rename for loops 

module load StdEnv/2020
module load bedtools/2.30.0
cd /home/celphin/scratch/Dryas/CHG_CHH/Svalbard_Metilene
sh metilene_prep.sh 

sed -i 's/Svalbard_W_C_input_files\///g' /home/celphin/scratch/Dryas/CHG_CHH/Svalbard_Metilene/Svalbard_metilene_W_C.input
sed -i 's/_CHH.bedGraph_sorted_tab.bedGraph//g' /home/celphin/scratch/Dryas/CHG_CHH/Svalbard_Metilene/Svalbard_metilene_W_C.input

#---------------------------------------------

#Svalbard specific adjustments,
nano metilene_run.sh

h1='W'
h2='C'
metilene_dir=/home/celphin/scratch/Dryas
input_dir="/home/celphin/scratch/Dryas/CHG_CHH/Svalbard_Metilene"

maxdist=$1
mincpgs=$2
mindiff=$3

output_name=Svalbard_CHH_"$h1"_"$h2"_"${maxdist}"_"${mincpgs}"_"${mindiff}"
in_metilene="Svalbard_metilene_"$h1"_"$h2".input"
threads=32


#metilene:
# params: maxdist, mincpgs, mindiff
salloc -c15 --time 7:00:00 --mem 120000m --account def-cronk

module load StdEnv/2020
module load bedtools/2.30.0 
sh metilene_run.sh 70 5 4
sh metilene_run.sh 150 5 4 
sh metilene_run.sh 70 5 0.7

# slurmstepd: error: Detected 3 oom_kill events in StepId=53035864.interactive. Some of the step tasks have been OOM Killed.
# srun: error: bl12440: task 0: Out Of Memory
# salloc: Relinquishing job allocation 53035864
# salloc: Job allocation 53035864 has been revoked.
#-------------------------------------------------
#Filter based on qval
#params: maxdist, mincpgs, mindiff, minmeandif, q-value

nano  metilene_filter_qval.sh
h1='W'
h2='C'
metilene_dir=/home/celphin/scratch/Dryas
input_dir="/home/celphin/scratch/Dryas/CHG_CHH/Svalbard_Metilene"
outputname=Svalbard_CHH_"$h1"_"$h2"_"${maxdist}"_"${mincpgs}"_"${mindiff}"

# to run
module load StdEnv/2023
module load r/4.4.0
module load gdal
module load udunits
module load python
export R_LIBS_USER=/home/celphin/R/x86_64-pc-linux-gnu-library/4.4/


sh metilene_filter_qval.sh 70 5 4 0.9 0.001
# Wrote 2177 DMRs with adj. p-value<0.001, a minimum absolute difference>=0.9, a minimum length [CpG]>=10 and a minimum length [nt]>=0
sh metilene_filter_qval.sh 150 5 4 0.9 1e-5
# Wrote 1397 DMRs with adj. p-value<1e-5, a minimum absolute difference>=0.9, a minimum length [CpG]>=10 and a minimum length [nt]>=0
sh metilene_filter_qval.sh 70 5 0.7 0.7 0.001
# Wrote 2270 DMRs with adj. p-value<0.001, a minimum absolute difference>=0.7, a minimum length [CpG]>=10 and a minimum length [nt]>=0

###################################################################################
#Sweden:
mkdir Sweden_Metilene
cd Sweden_Metilene
cp /home/celphin/scratch/Dryas/MS_scripts/metilene*.sh .

tmux new-session -s Sweden_DMRS
tmux attach-session -t Sweden_DMRS
salloc -c1 --time 7:00:00 --mem 120000m --account def-henryg

#----------------------------------------
#Metilene input prep:
module load StdEnv/2020
module load bedtools/2.30.0

#Sweden adjustments
nano metilene_prep.sh 

#h1,h2: 
h1="W", h2="C"
#Input directories:
input_dir=Sweden_${h1}_${h2}_input_files 
in_metilene="Sweden_metilene_"$h1"_"$h2".input"
methylseq_output_dir="/home/celphin/scratch/Dryas/CHG_CHH/bedgraph"
#Copy Loop:
cp ${methylseq_output_dir}/*LAT*_CHH.bedGraph ./${input_dir}
#Comment out rename for loops 

module load StdEnv/2020
module load bedtools/2.30.0
cd /home/celphin/scratch/Dryas/CHG_CHH/Sweden_Metilene
sh metilene_prep.sh 

sed -i 's/Sweden_W_C_input_files\///g' /home/celphin/scratch/Dryas/CHG_CHH/Sweden_Metilene/Sweden_metilene_W_C.input
sed -i 's/_CHH.bedGraph_sorted_tab.bedGraph//g' /home/celphin/scratch/Dryas/CHG_CHH/Sweden_Metilene/Sweden_metilene_W_C.input

#---------------------------------------------

#Sweden  specific adjustments,
nano metilene_run.sh

h1='W'
h2='C'
metilene_dir=/home/celphin/scratch/Dryas
input_dir="/home/celphin/scratch/Dryas/CHG_CHH/Sweden_Metilene"

maxdist=$1
mincpgs=$2
mindiff=$3

output_name=Sweden_CHH_"$h1"_"$h2"_"${maxdist}"_"${mincpgs}"_"${mindiff}"
in_metilene="Sweden_metilene_"$h1"_"$h2".input"
threads=15


#metilene:
# params: maxdist, mincpgs, mindiff
salloc -c15 --time 7:00:00 --mem 200000m --account def-cronk
# needs more Memory

module load StdEnv/2020
module load bedtools/2.30.0 
sh metilene_run.sh 70 5 4
sh metilene_run.sh 150 5 4 
sh metilene_run.sh 70 5 0.7

#-------------------------------------------------
#Filter based on qval
#params: maxdist, mincpgs, mindiff, minmeandif, q-value

nano  metilene_filter_qval.sh
h1='W'
h2='C'
metilene_dir=/home/celphin/scratch/Dryas
input_dir="/home/celphin/scratch/Dryas/CHG_CHH/Sweden_Metilene"
outputname=Sweden_CHH_"$h1"_"$h2"_"${maxdist}"_"${mincpgs}"_"${mindiff}"

# to run
module load StdEnv/2023
module load r/4.4.0
module load gdal
module load udunits
module load python
export R_LIBS_USER=/home/celphin/R/x86_64-pc-linux-gnu-library/4.4/

sh metilene_filter_qval.sh 70 5 4 0.9 0.001
#  Wrote 3336 DMRs with adj. p-value<0.001, a minimum absolute difference>=0.9, a minimum length [CpG]>=10 and a minimum length [nt]>=0
sh metilene_filter_qval.sh 150 5 4 0.9 1e-5
# Wrote 2206 DMRs with adj. p-value<1e-5, a minimum absolute difference>=0.9, a minimum length [CpG]>=10 and a minimum length [nt]>=0
sh metilene_filter_qval.sh 70 5 0.7 0.7 0.001
# Wrote 3737 DMRs with adj. p-value<0.001, a minimum absolute difference>=0.7, a minimum length [CpG]>=10 and a minimum length [nt]>=0


#################################################################################
#Alaska 
mkdir Alaska_Metilene
cd Alaska_Metilene
cp /home/celphin/scratch/Dryas/MS_scripts/metilene*.sh .

tmux new-session -s Alaska_DMRS
tmux attach-session -t Alaska_DMRS
salloc -c32 --time 4:50:00 --mem 120000m --account def-rieseber

#----------------------------------------
#Metilene input prep:
module load StdEnv/2020
module load bedtools/2.30.0

#Alaska adjustments
nano metilene_prep.sh 

#h1,h2: 
h1="W", h2="C"
#Input directories:
input_dir=Alaska_${h1}_${h2}_input_files 
in_metilene="Alaska_metilene_"$h1"_"$h2".input"
methylseq_output_dir="/home/celphin/scratch/Dryas/CHG_CHH/bedgraph"
#Copy Loop:
cp ${methylseq_output_dir}/*ALAS*_CHH.bedGraph ./${input_dir}
#Comment out rename for loops 

module load StdEnv/2020
module load bedtools/2.30.0
cd /home/celphin/scratch/Dryas/CHG_CHH/Alaska_Metilene
sh metilene_prep.sh 

sed -i 's/Alaska_W_C_input_files\///g' /home/celphin/scratch/Dryas/CHG_CHH/Alaska_Metilene/Alaska_metilene_W_C.input
sed -i 's/_CHH.bedGraph_sorted_tab.bedGraph//g' /home/celphin/scratch/Dryas/CHG_CHH/Alaska_Metilene/Alaska_metilene_W_C.input

#----------------------------------------
#Alaska specific adjustments,
nano metilene_run.sh

h1='W'
h2='C'
metilene_dir=/home/celphin/scratch/Dryas
input_dir="/home/celphin/scratch/Dryas/CHG_CHH/Alaska_Metilene"

maxdist=$1
mincpgs=$2
mindiff=$3

output_name=Alaska_CHH_"$h1"_"$h2"_"${maxdist}"_"${mincpgs}"_"${mindiff}"
in_metilene="Alaska_metilene_"$h1"_"$h2".input"
threads=32


#metilene:
# params: maxdist, mincpgs, mindiff
salloc -c32 --time 7:00:00 --mem 200000m --account def-cronk

module load StdEnv/2020
module load bedtools/2.30.0 
sh metilene_run.sh 70 5 4
sh metilene_run.sh 150 5 4 
sh metilene_run.sh 70 5 0.7

#-------------------------------------------------
#Filter based on qval
#params: maxdist, mincpgs, mindiff, minmeandif, q-value

nano  metilene_filter_qval.sh
h1='W'
h2='C'
metilene_dir=/home/celphin/scratch/Dryas
input_dir="/home/celphin/scratch/Dryas/CHG_CHH/Alaska_Metilene"
outputname=Alaska_CHH_"$h1"_"$h2"_"${maxdist}"_"${mincpgs}"_"${mindiff}"

# to run
module load StdEnv/2023
module load r/4.4.0
module load gdal
module load udunits
module load python
export R_LIBS_USER=/home/celphin/R/x86_64-pc-linux-gnu-library/4.4/

sh metilene_filter_qval.sh 70 5 4 0.9 0.001
# Wrote 3097 DMRs with adj. p-value<0.001, a minimum absolute difference>=0.9, a minimum length [CpG]>=10 and a minimum length [nt]>=0
sh metilene_filter_qval.sh 150 5 4 0.9 1e-5
# Wrote 1937 DMRs with adj. p-value<1e-5, a minimum absolute difference>=0.9, a minimum length [CpG]>=10 and a minimum length [nt]>=0
sh metilene_filter_qval.sh 70 5 0.7 0.7 0.001
# Wrote 3381 DMRs with adj. p-value<0.001, a minimum absolute difference>=0.7, a minimum length [CpG]>=10 and a minimum length [nt]>=0

##############################################################
cd ..
#Copy bedGraph files to directory above
cp Nunavut_Metilene/Nunavut_CHH_W_C_150_5_4_0.9_qval.1e-5.bedgraph Nunavut_CHH_W_C.bedGraph
cp Alaska_Metilene/Alaska_CHH_W_C_150_5_4_0.9_qval.1e-5.bedgraph Alaska_CHH_W_C.bedGraph
cp Svalbard_Metilene/Svalbard_CHH_W_C_150_5_4_0.9_qval.1e-5.bedgraph Svalbard_CHH_W_C.bedGraph
cp Sweden_Metilene/Sweden_CHH_W_C_150_5_4_0.9_qval.1e-5.bedgraph Sweden_CHH_W_C.bedGraph

#Bedtools intersect all 4
module load StdEnv/2023 bedtools/2.31.0
bedtools intersect -u -a Nunavut_CHH_W_C.bedGraph -b Svalbard_CHH_W_C.bedGraph| bedtools intersect -u -a stdin -b Sweden_CHH_W_C.bedGraph | bedtools intersect -u -a stdin -b Alaska_CHH_W_C.bedGraph > ALL_Sites_intersect_DMRs.bedgraph

#Bedtools intersect all but Nunavut
bedtools intersect -u -a Svalbard_CHH_W_C.bedGraph -b Sweden_CHH_W_C.bedGraph | bedtools intersect -u -a stdin -b Alaska_CHH_W_C.bedGraph > SVAL_SWED_ALAS_Sites_intersect_DMRs.bedgraph
# Do1_01_a00001   7699414 7699600 -8.213185
# Do1_02_a00001   8510664 8510827 12.079776
# Do1_02_a00004   6400235 6400478 8.718939
# Do1_04_a00005   405319  405494  -5.244285
# Do1_06_a00001   1983414 1983677 4.578070
# Do1_06_a00001   6538326 6538612 4.997452

#Bedtools intersect all but Svalbard
bedtools intersect -u -a Nunavut_CHH_W_C.bedGraph -b Sweden_CHH_W_C.bedGraph | bedtools intersect -u -a stdin -b Alaska_CHH_W_C.bedGraph > ALEX_SWED_ALAS_Sites_intersect_DMRs.bedgraph

#Bedtools intersect all but Sweden
bedtools intersect -u -a Nunavut_CHH_W_C.bedGraph -b Svalbard_CHH_W_C.bedGraph| bedtools intersect -u -a stdin -b Alaska_CHH_W_C.bedGraph > ALEX_SVAL_ALAS_Sites_intersect_DMRs.bedgraph

#Bedtools intersect all but Alaska
bedtools intersect -u -a Nunavut_CHH_W_C.bedGraph -b Svalbard_CHH_W_C.bedGraph| bedtools intersect -u -a stdin -b Sweden_CHH_W_C.bedGraph > ALEX_SVAL_SWED_Sites_intersect_DMRs.bedgraph

wc -l ALL_Sites_intersect_DMRs.bedgraph
# 6
wc -l SVAL_SWED_ALAS_Sites_intersect_DMRs.bedgraph
# 37
wc -l ALEX_SWED_ALAS_Sites_intersect_DMRs.bedgraph
# 19
wc -l ALEX_SVAL_ALAS_Sites_intersect_DMRs.bedgraph
#20
wc -l ALEX_SVAL_SWED_Sites_intersect_DMRs.bedgraph
# 26

###############################
# Phenology
cd /home/celphin/scratch/Dryas/CHG_CHH
mkdir Phenology_Metilene
cd Phenology_Metilene
cp /home/celphin/scratch/Dryas/MS_scripts/metilene*.sh .

tmux new-session -s Phenology_DMRS
tmux attach-session -t Phenology_DMRS

# copy and rename the Sen 

# MatFl.Cass.4C.4.524.F112579    
# MatFl.Will.3C.100.414.F112575    
# MatFl.Fert.5C.97.1F.F112583 
# MatFl.Cass.5W.130.525.F112580  
# MatFl.Will.4W.13.417.F112576     
# MatFl.Fert.6W.110.3F.F112584 
# MatFl.Mead.1C.33.446.F112577   
# Sen.FERT5C.1F.97_CHH.bedGraph
# Sen.FERT6W.3F.110_CHH.bedGraph
# Sen.CASS10W.544.60_CHH.bedGraph  
# Sen.MEAD1W.444.116_CHH.bedGraph
# Sen.CASS4C.524.4_CHH.bedGraph    
# Sen.WILL3C.414.100_CHH.bedGraph
# Sen.CASS5W.525.130_CHH.bedGraph  
# Sen.WILL4W.417.13_CHH.bedGraph
# Sen.MEAD1C.446.33_CHH.bedGraph

cd /home/celphin/scratch/Dryas/CHG_CHH/bedgraph/Phenology
# add bedGraph
for file in MatFl*; do
  if [ -f "$file" ]; then
    mv -- "$file" "$file.bedGraph"
  fi
done

#----------------------------------------
#Metilene input prep:

#Phenology specific adjustments:
nano  metilene_prep.sh

#h1,h2: 
h1="Mat"
h2="Sen"
#Input directories:
input_dir=${h1}_${h2}_input_files 
in_metilene="metilene_"$h1"_"$h2".input"
methylseq_output_dir="/home/celphin/scratch/Dryas/CHG_CHH/bedgraph/Phenology"


#change copy process because Senesence bedgraphs don't have deduplicated:
cp -r ${methylseq_output_dir}/*.bedGraph ./${input_dir}

#Comment out rename for loops shown below (already with proper prefix)
    ##for bg in _*_${h1}*; do mv "$bg" "${h1}_${bg}"; done
    ##for bg in _*_${h2}*; do mv "$bg" "${h2}_${bg}"; done

tmux new-session -s Phenology_DMRS
tmux attach-session -t Phenology_DMRS

salloc -c1 --time 7:00:00 --mem 120000m --account def-cronk
module load StdEnv/2020
module load bedtools/2.30.0

sh metilene_prep.sh 


#----------------------------------------
#Run metilene:

#Phenology specific adjustments,
nano metilene_run.sh

h1='Mat'
h2='Sen'
metilene_dir=/home/celphin/scratch/Dryas
input_dir="/home/celphin/scratch/Dryas/CHG_CHH/Phenology_Metilene"

maxdist=$1
mincpgs=$2
mindiff=$3

output_name=Phenology_CHH_"$h1"_"$h2"_"${maxdist}"_"${mincpgs}"_"${mindiff}"
in_metilene=metilene_Mat_Sen.input
threads=15

#-------
# remove the folder name from the metilene input file

sed -i 's/Mat_Sen_input_files\///g'  metilene_Mat_Sen.input

#metilene:
# params: maxdist, mincpgs, mindiff
salloc -c15 --time 7:00:00 --mem 200000m --account def-cronk

module load StdEnv/2020
module load bedtools/2.30.0 
sh metilene_run.sh 70 5 4
sh metilene_run.sh 150 5 4 
sh metilene_run.sh 70 5 0.7
# done

#-------------------------------------
#Filter based on qval
#params: maxdist, mincpgs, mindiff, minmeandif, q-value

nano  metilene_filter_qval.sh
h1='Mat'
h2='Sen'
metilene_dir=/home/celphin/scratch/Dryas
input_dir="/home/celphin/scratch/Dryas/CHG_CHH/Phenology_Metilene"
outputname=Phenology_CHH_"$h1"_"$h2"_"${maxdist}"_"${mincpgs}"_"${mindiff}"

# to run
module load StdEnv/2023
module load r/4.4.0
module load gdal
module load udunits
module load python
export R_LIBS_USER=/home/celphin/R/x86_64-pc-linux-gnu-library/4.4/

sh metilene_filter_qval.sh 70 5 4 0.9 0.001
sh metilene_filter_qval.sh 150 5 4 0.9 1e-5
sh metilene_filter_qval.sh 70 5 0.7 0.7 0.001


#######################
# Seedlings

cd /home/celphin/scratch/Dryas/CHG_CHH
mkdir Seedling_Metilene
cd Seedling_Metilene
cp /home/celphin/scratch/Dryas/MS_scripts/metilene*.sh .

tmux new-session -s Seedling_Warming_DMRS
tmux attach-session -t Seedling_Warming_DMRS

cd /home/celphin/scratch/Dryas/CHG_CHH/bedgraph/seedling

rename SE C_SE SE.*.C.*
rename SE W_SE SE.*.W.*


cd /home/celphin/scratch/Dryas/CHG_CHH/Seedling_Metilene

#----------------------------------------
#Metilene input prep:

#Seedling, Warming control specific adjustments:
nano metilene_prep.sh

#h1,h2: 
h1="W"
h2="C"
#Input directories:
input_dir=SE_${h1}_${h2}_input_files 
in_metilene="SE_metilene_"$h1"_"$h2".input"
methylseq_output_dir="/home/celphin/scratch/Dryas/CHG_CHH/bedgraph/seedling"

#Copy Loop:
cp ${methylseq_output_dir}/*SE* ./${input_dir}

# make sbatch?
salloc -c1 --time 7:00:00 --mem 120000m --account def-cronk
module load StdEnv/2020
module load bedtools/2.30.0
cd /home/celphin/scratch/Dryas/CHG_CHH/Seedling_Metilene
sh metilene_prep.sh 


#----------------------------------------
#Run metilene:
#Seedling,Warming control specific adjustments,
nano metilene_run.sh

h1='W'
h2='C'
metilene_dir=/home/celphin/scratch/Dryas
input_dir="/home/celphin/scratch/Dryas/CHG_CHH/Seedling_Metilene"

maxdist=$1
mincpgs=$2
mindiff=$3

output_name=SE_CHH_"$h1"_"$h2"_"${maxdist}"_"${mincpgs}"_"${mindiff}"
in_metilene="SE_metilene_"$h1"_"$h2".input"
threads=15

# remove the folder name from the metilene input file

sed -i 's/SE_W_C_input_files\///g' SE_metilene_W_C.input

#metilene:
# params: maxdist, mincpgs, mindiff
salloc -c15 --time 7:00:00 --mem 200000m --account def-cronk

module load StdEnv/2020
module load bedtools/2.30.0 
sh metilene_run.sh 70 5 4
sh metilene_run.sh 150 5 4 
sh metilene_run.sh 70 5 0.7
# done

#----------------------------------------------------
#Filter based on qval
#params: maxdist, mincpgs, mindiff, minmeandif, q-value

nano  metilene_filter_qval.sh
h1="W" 
h2="C"
metilene_dir=/home/celphin/scratch/Dryas
input_dir="/home/celphin/scratch/Dryas/CHG_CHH/Seedling_Metilene"
outputname=Phenology_CHH_"$h1"_"$h2"_"${maxdist}"_"${mincpgs}"_"${mindiff}"

# to run
module load StdEnv/2023
module load r/4.4.0
module load gdal
module load udunits
module load python
export R_LIBS_USER=/home/celphin/R/x86_64-pc-linux-gnu-library/4.4/

sh metilene_filter_qval.sh 70 5 4 0.9 0.001
# Wrote 3951 DMRs with adj. p-value<0.001, a minimum absolute difference>=0.9, a minimum length [CpG]>=10 and a minimum length [nt]>=0
sh metilene_filter_qval.sh 150 5 4 0.9 1e-5
# Wrote 2356 DMRs with adj. p-value<1e-5, a minimum absolute difference>=0.9, a minimum length [CpG]>=10 and a minimum length [nt]>=0
sh metilene_filter_qval.sh 70 5 0.7 0.7 0.001
# Wrote 4008 DMRs with adj. p-value<0.001, a minimum absolute difference>=0.7, a minimum length [CpG]>=10 and a minimum length [nt]>=0

############################
# Intersect Seedlings

cp Seedling_Metilene/SE_CHH_W_C_150_5_4_0.9_qval.1e-5.bedgraph SE_CHH_W_C.bedGraph
cp Phenology_Metilene/Phenology_CHH_Mat_Sen_150_5_4_0.9_qval.1e-5.bedgraph Phenology_CHH_Mat_Sen.bedGraph

#Bedtools intersect all 4
module load StdEnv/2023 bedtools/2.31.0
bedtools intersect -u -a SE_CHH_W_C.bedGraph -b Sweden_CHH_W_C.bedGraph > intersect_SE_Sweden_W_C_CHH.bedGraph
# 191 intersect_SE_Sweden_W_C_CHH.bedGraph
bedtools intersect -u -a Phenology_CHH_Mat_Sen.bedGraph -b Nunavut_CHH_W_C.bedGraph > intersect_Pheno_Nunavut_Mat_Sen_CHH.bedGraph
# 631 intersect_Pheno_Nunavut_Mat_Sen_CHH.bedGraph

#############################
