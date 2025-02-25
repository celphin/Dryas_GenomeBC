#########################################################
# Dryas DMR calling
# Dec 2024
#############################################################
# Run metilene on CHG and CHH contexts randomly

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


#########################
# Nunavut
cd /home/celphin/scratch/Dryas/CHG_CHH/

mkdir Nunavut_Metilene_random
cd Nunavut_Metilene_random
cp /home/celphin/scratch/Dryas/CHG_CHH/Nunavut_Metilene/* /home/celphin/scratch/Dryas/CHG_CHH/Nunavut_Metilene_random/

tmux new-session -s Nunavut_DMRS
tmux attach-session -t Nunavut_DMRS
salloc -c15 --time 7:00:00 --mem 200000m --account def-cronk

# Define the path to the input file 
cd /home/celphin/scratch/Dryas/CHG_CHH/Nunavut_Metilene_random/

cp ../Nunavut_Metilene/Nunavut_metilene_W_C.input .

head Nunavut_metilene_W_C.input > Nunavut_metilene_W_C1.input

# check
head -n 1 Nunavut_metilene_W_C1.input

awk '
BEGIN {
    OFS="\t"
}
NR == 1 {
    # This is the header line, process it separately
    header = $1 "\t" $2   # Preserve the first two columns as they are

    # Create an array for all columns after the first two
    for (i = 3; i <= NF; i++) {
        cols[i-2] = $i
    }

    # Calculate the number of columns
    n = length(cols)

    # Shuffle the columns except the first two
    srand()
    for (i = n; i > 1; i--) {
        j = int(rand() * i) + 1
        tmp = cols[i]
        cols[i] = cols[j]
        cols[j] = tmp
    }

    # Add "W_" to the first half of the shuffled columns
    # and "C_" to the second half
    half = int(n / 2)
    for (i = 1; i <= n; i++) {
        if (i <= half) {
            header = header "\tW_" cols[i]  # First half gets W_
        } else {
            header = header "\tC_" cols[i]  # Second half gets C_
        }
    }
    
    print header  # Print the modified header
    next  # Skip to the next line (data rows)
}

# This is for processing the data rows
{
    # Print the first two columns unchanged
    printf "%s\t%s", $1, $2

    # Save the remaining columns
    for (i = 3; i <= NF; i++) {
        cols[i-2] = $i
    }

    # Shuffle the columns except the first two
    n = length(cols)
    srand()
    for (i = n; i > 1; i--) {
        j = int(rand() * i) + 1
        tmp = cols[i]
        cols[i] = cols[j]
        cols[j] = tmp
    }

    # Print the shuffled columns, but without adding W_ or C_
    for (i = 1; i <= n; i++) {
        printf "\t%s", cols[i]
    }
    print ""  # End the row
}
' Nunavut_metilene_W_C.input > Nunavut_metilene_W_C2.input

#-------------------------------------
# replace NA with . or - 

sed 's/NA/-/g' Nunavut_metilene_W_C2.input > Nunavut_metilene_W_C3.input

#------------------------------------------------------------------

#Nunavut specific adjustments,
nano metilene_run.sh

input_dir="/home/celphin/scratch/Dryas/CHG_CHH/Nunavut_Metilene_random"
...
output_name=Nunavut_CHH_random_"$h1"_"$h2"_"${maxdist}"_"${mincpgs}"_"${mindiff}"

in_metilene="Nunavut_metilene_"$h1"_"$h2"2.input"

threads=15

#-----------------------------------
#metilene:
# params: maxdist, mincpgs, mindiff
salloc -c15 --time 7:00:00 --mem 200000m --account def-cronk

module load StdEnv/2020
module load bedtools/2.30.0 

# with in_metilene="Nunavut_metilene_"$h1"_"$h2"2.input"
sh metilene_run.sh 150 5 4 

# with in_metilene="Nunavut_metilene_"$h1"_"$h2"3.input"
sh metilene_run3.sh 150 5 4 

# run with Minimal # of non-missing values for estimating missing values in g1 and g2
sh metilene_run_XY.sh 150 5 4 

#  $metilene_dir/metilene_v0.2-8/metilene "${input_dir}/${in_metilene}" \--maxdist ${maxdist} --mincpgs ${mincpgs} --minMethDiff ${mindiff} --mode 1 --threads ${threads} -a ${h1} -b ${h2} -v 0.2 --minNoA 20 --minNoB 20 > "${output_name}"

#-------------------------------------------------
#Filter based on qval
#params: maxdist, mincpgs, mindiff, minmeandif, q-value

nano  metilene_filter_qval.sh

input_dir="/home/celphin/scratch/Dryas/CHG_CHH/Nunavut_Metilene_random"
outputname=Nunavut_CHH_random_"$h1"_"$h2"_"${maxdist}"_"${mincpgs}"_"${mindiff}"

# to run
module load StdEnv/2023
module load r/4.4.0
module load gdal
module load udunits
module load python
export R_LIBS_USER=/home/celphin/R/x86_64-pc-linux-gnu-library/4.4/

sh metilene_filter_qval.sh 150 5 4 0.9 1e-5
#  Wrote 1278 DMRs with adj. p-value<1e-5, a minimum absolute difference>=0.9, a minimum length [CpG]>=10 and a minimum length [nt]>=0
cp metilene_filter_qval.sh metilene_filter_qval3.sh

nano metilene_filter_qval3.sh
sh metilene_filter_qval3.sh 150 5 4 0.9 1e-5
# Wrote 1278 DMRs with adj. p-value<1e-5, a minimum absolute difference>=0.9, a minimum length [CpG]>=10 and a minimum length [nt]>=0

# original data but filling in 80%
sh metilene_filter_qval_XY.sh 150 5 4 0.9 1e-5
# Wrote 1186 DMRs with adj. p-value<1e-5, a minimum absolute difference>=0.9, a minimum length [CpG]>=10 and a minimum length [nt]>=0


####################################################################################
#Svalbard:

cd /home/celphin/scratch/Dryas/CHG_CHH/

mkdir Svalbard_Metilene_random
cd Svalbard_Metilene_random
cp /home/celphin/scratch/Dryas/CHG_CHH/Svalbard_Metilene/* /home/celphin/scratch/Dryas/CHG_CHH/Svalbard_Metilene_random/

tmux new-session -s Svalbard_DMRS
tmux attach-session -t Svalbard_DMRS

# Define the path to the input file 
cd /home/celphin/scratch/Dryas/CHG_CHH/Svalbard_Metilene_random/

cp ../Svalbard_Metilene/Svalbard_metilene_W_C.input .

head Svalbard_metilene_W_C.input > Svalbard_metilene_W_C1.input

# check
head -n 1 Svalbard_metilene_W_C1.input

awk '
BEGIN {
    OFS="\t"
}
NR == 1 {
    # This is the header line, process it separately
    header = $1 "\t" $2   # Preserve the first two columns as they are

    # Create an array for all columns after the first two
    for (i = 3; i <= NF; i++) {
        cols[i-2] = $i
    }

    # Calculate the number of columns
    n = length(cols)

    # Shuffle the columns except the first two
    srand()
    for (i = n; i > 1; i--) {
        j = int(rand() * i) + 1
        tmp = cols[i]
        cols[i] = cols[j]
        cols[j] = tmp
    }

    # Add "W_" to the first half of the shuffled columns
    # and "C_" to the second half
    half = int(n / 2)
    for (i = 1; i <= n; i++) {
        if (i <= half) {
            header = header "\tW_" cols[i]  # First half gets W_
        } else {
            header = header "\tC_" cols[i]  # Second half gets C_
        }
    }
    
    print header  # Print the modified header
    next  # Skip to the next line (data rows)
}

# This is for processing the data rows
{
    # Print the first two columns unchanged
    printf "%s\t%s", $1, $2

    # Save the remaining columns
    for (i = 3; i <= NF; i++) {
        cols[i-2] = $i
    }

    # Shuffle the columns except the first two
    n = length(cols)
    srand()
    for (i = n; i > 1; i--) {
        j = int(rand() * i) + 1
        tmp = cols[i]
        cols[i] = cols[j]
        cols[j] = tmp
    }

    # Print the shuffled columns, but without adding W_ or C_
    for (i = 1; i <= n; i++) {
        printf "\t%s", cols[i]
    }
    print ""  # End the row
}
' Svalbard_metilene_W_C.input > Svalbard_metilene_W_C2.input

#-------------------------------------
# replace NA with . or - 

sed 's/NA/-/g' Svalbard_metilene_W_C2.input > Svalbard_metilene_W_C3.input

#------------------------------------------------------------------

#Svalbard specific adjustments,
nano metilene_run.sh

input_dir="/home/celphin/scratch/Dryas/CHG_CHH/Svalbard_Metilene_random"
...
output_name=Svalbard_CHH_random_"$h1"_"$h2"_"${maxdist}"_"${mincpgs}"_"${mindiff}"

threads=15


#metilene:
# params: maxdist, mincpgs, mindiff
salloc -c7 --time 7:00:00 --mem 200000m --account def-cronk

module load StdEnv/2020
module load bedtools/2.30.0 
cd /home/celphin/scratch/Dryas/CHG_CHH/Svalbard_Metilene_random
# with in_metilene="Svalbard_metilene_"$h1"_"$h2"2.input"
sh metilene_run.sh 150 5 4 


#-------------------------------------------------
#Filter based on qval
#params: maxdist, mincpgs, mindiff, minmeandif, q-value

nano  metilene_filter_qval.sh

input_dir="/home/celphin/scratch/Dryas/CHG_CHH/Svalbard_Metilene_random"
outputname=Svalbard_CHH_random_"$h1"_"$h2"_"${maxdist}"_"${mincpgs}"_"${mindiff}"

# to run
module load StdEnv/2023
module load r/4.4.0
module load gdal
module load udunits
module load python
export R_LIBS_USER=/home/celphin/R/x86_64-pc-linux-gnu-library/4.4/

sh metilene_filter_qval.sh 150 5 4 0.9 1e-5
# Wrote 2002 DMRs with adj. p-value<1e-5, a minimum absolute difference>=0.9, a minimum length [CpG]>=10 and a minimum length [nt]>=0

###############################################################
# Alaska
cd /home/celphin/scratch/Dryas/CHG_CHH/

mkdir Alaska_Metilene_random
cd Alaska_Metilene_random
cp /home/celphin/scratch/Dryas/CHG_CHH/Alaska_Metilene/* /home/celphin/scratch/Dryas/CHG_CHH/Alaska_Metilene_random/

tmux new-session -s Alaska_DMRS
tmux attach-session -t Alaska_DMRS
salloc -c15 --time 7:00:00 --mem 200000m --account def-cronk

# Define the path to the input file 
cd /home/celphin/scratch/Dryas/CHG_CHH/Alaska_Metilene_random/

cp ../Alaska_Metilene/Alaska_metilene_W_C.input .

head Alaska_metilene_W_C.input > Alaska_metilene_W_C1.input

# check
head -n 1 Alaska_metilene_W_C1.input

awk '
BEGIN {
    OFS="\t"
}
NR == 1 {
    # This is the header line, process it separately
    header = $1 "\t" $2   # Preserve the first two columns as they are

    # Create an array for all columns after the first two
    for (i = 3; i <= NF; i++) {
        cols[i-2] = $i
    }

    # Calculate the number of columns
    n = length(cols)

    # Shuffle the columns except the first two
    srand()
    for (i = n; i > 1; i--) {
        j = int(rand() * i) + 1
        tmp = cols[i]
        cols[i] = cols[j]
        cols[j] = tmp
    }

    # Add "W_" to the first half of the shuffled columns
    # and "C_" to the second half
    half = int(n / 2)
    for (i = 1; i <= n; i++) {
        if (i <= half) {
            header = header "\tW_" cols[i]  # First half gets W_
        } else {
            header = header "\tC_" cols[i]  # Second half gets C_
        }
    }
    
    print header  # Print the modified header
    next  # Skip to the next line (data rows)
}

# This is for processing the data rows
{
    # Print the first two columns unchanged
    printf "%s\t%s", $1, $2

    # Save the remaining columns
    for (i = 3; i <= NF; i++) {
        cols[i-2] = $i
    }

    # Shuffle the columns except the first two
    n = length(cols)
    srand()
    for (i = n; i > 1; i--) {
        j = int(rand() * i) + 1
        tmp = cols[i]
        cols[i] = cols[j]
        cols[j] = tmp
    }

    # Print the shuffled columns, but without adding W_ or C_
    for (i = 1; i <= n; i++) {
        printf "\t%s", cols[i]
    }
    print ""  # End the row
}
' Alaska_metilene_W_C.input > Alaska_metilene_W_C2.input

#-------------------------------------
# replace NA with . or - 

sed 's/NA/-/g' Alaska_metilene_W_C2.input > Alaska_metilene_W_C3.input

#------------------------------------------------------------------

#Alaska specific adjustments,
nano metilene_run.sh

input_dir="/home/celphin/scratch/Dryas/CHG_CHH/Alaska_Metilene_random"
...
output_name=Alaska_CHH_random_"$h1"_"$h2"_"${maxdist}"_"${mincpgs}"_"${mindiff}"

threads=15


#metilene:
# params: maxdist, mincpgs, mindiff
salloc -c15 --time 7:00:00 --mem 200000m --account def-cronk

module load StdEnv/2020
module load bedtools/2.30.0 
cd /home/celphin/scratch/Dryas/CHG_CHH/Alaska_Metilene_random
# with in_metilene="Alaska_metilene_"$h1"_"$h2"2.input"
sh metilene_run.sh 150 5 4 

#-------------------------------------------------
#Filter based on qval
#params: maxdist, mincpgs, mindiff, minmeandif, q-value

nano  metilene_filter_qval.sh

input_dir="/home/celphin/scratch/Dryas/CHG_CHH/Alaska_Metilene_random"
outputname=Alaska_CHH_random_"$h1"_"$h2"_"${maxdist}"_"${mincpgs}"_"${mindiff}"

# to run
module load StdEnv/2023
module load r/4.4.0
module load gdal
module load udunits
module load python
export R_LIBS_USER=/home/celphin/R/x86_64-pc-linux-gnu-library/4.4/

sh metilene_filter_qval.sh 150 5 4 0.9 1e-5
#  Wrote 1780 DMRs with adj. p-value<1e-5, a minimum absolute difference>=0.9, a minimum length [CpG]>=10 and a minimum length [nt]>=0

###############################################################
# Sweden
cd /home/celphin/scratch/Dryas/CHG_CHH/

mkdir Sweden_Metilene_random
cd Sweden_Metilene_random
cp /home/celphin/scratch/Dryas/CHG_CHH/Sweden_Metilene/* /home/celphin/scratch/Dryas/CHG_CHH/Sweden_Metilene_random/

tmux new-session -s Sweden_DMRS
tmux attach-session -t Sweden_DMRS
salloc -c15 --time 7:00:00 --mem 200000m --account def-cronk

# Define the path to the input file 
cd /home/celphin/scratch/Dryas/CHG_CHH/Sweden_Metilene_random/

cp ../Sweden_Metilene/Sweden_metilene_W_C.input .

head Sweden_metilene_W_C.input > Sweden_metilene_W_C1.input

# check
head -n 1 Sweden_metilene_W_C1.input

awk '
BEGIN {
    OFS="\t"
}
NR == 1 {
    # This is the header line, process it separately
    header = $1 "\t" $2   # Preserve the first two columns as they are

    # Create an array for all columns after the first two
    for (i = 3; i <= NF; i++) {
        cols[i-2] = $i
    }

    # Calculate the number of columns
    n = length(cols)

    # Shuffle the columns except the first two
    srand()
    for (i = n; i > 1; i--) {
        j = int(rand() * i) + 1
        tmp = cols[i]
        cols[i] = cols[j]
        cols[j] = tmp
    }

    # Add "W_" to the first half of the shuffled columns
    # and "C_" to the second half
    half = int(n / 2)
    for (i = 1; i <= n; i++) {
        if (i <= half) {
            header = header "\tW_" cols[i]  # First half gets W_
        } else {
            header = header "\tC_" cols[i]  # Second half gets C_
        }
    }
    
    print header  # Print the modified header
    next  # Skip to the next line (data rows)
}

# This is for processing the data rows
{
    # Print the first two columns unchanged
    printf "%s\t%s", $1, $2

    # Save the remaining columns
    for (i = 3; i <= NF; i++) {
        cols[i-2] = $i
    }

    # Shuffle the columns except the first two
    n = length(cols)
    srand()
    for (i = n; i > 1; i--) {
        j = int(rand() * i) + 1
        tmp = cols[i]
        cols[i] = cols[j]
        cols[j] = tmp
    }

    # Print the shuffled columns, but without adding W_ or C_
    for (i = 1; i <= n; i++) {
        printf "\t%s", cols[i]
    }
    print ""  # End the row
}
' Sweden_metilene_W_C.input > Sweden_metilene_W_C2.input

#-------------------------------------
# replace NA with . or - 

sed 's/NA/-/g' Sweden_metilene_W_C2.input > Sweden_metilene_W_C3.input

#------------------------------------------------------------------

#Sweden specific adjustments,
nano metilene_run.sh

input_dir="/home/celphin/scratch/Dryas/CHG_CHH/Sweden_Metilene_random"
...
output_name=Sweden_CHH_random_"$h1"_"$h2"_"${maxdist}"_"${mincpgs}"_"${mindiff}"

threads=7


#metilene:
# params: maxdist, mincpgs, mindiff
salloc -c7 --time 7:00:00 --mem 200000m --account def-cronk

module load StdEnv/2020
module load bedtools/2.30.0 
cd /home/celphin/scratch/Dryas/CHG_CHH/Sweden_Metilene_random

# with in_metilene="Sweden_metilene_"$h1"_"$h2"2.input"
sh metilene_run.sh 150 5 4 

#-------------------------------------------------
#Filter based on qval
#params: maxdist, mincpgs, mindiff, minmeandif, q-value

nano  metilene_filter_qval.sh

input_dir="/home/celphin/scratch/Dryas/CHG_CHH/Sweden_Metilene_random"
outputname=Sweden_CHH_random_"$h1"_"$h2"_"${maxdist}"_"${mincpgs}"_"${mindiff}"

# to run
module load StdEnv/2023
module load r/4.4.0
module load gdal
module load udunits
module load python
export R_LIBS_USER=/home/celphin/R/x86_64-pc-linux-gnu-library/4.4/

sh metilene_filter_qval.sh 150 5 4 0.9 1e-5
# 1653

##############################################################
cd ..
#Copy bedGraph files to directory above
cp Nunavut_Metilene_random/Nunavut_CHH_random_W_C_150_5_4_0.9_qval.1e-5.bedgraph Nunavut_CHH_random_W_C.bedGraph
cp Alaska_Metilene_random/Alaska_CHH_random_W_C_150_5_4_0.9_qval.1e-5.bedgraph Alaska_CHH_random_W_C.bedGraph
cp Svalbard_Metilene_random/Svalbard_CHH_random_W_C_150_5_4_0.9_qval.1e-5.bedgraph Svalbard_CHH_random_W_C.bedGraph
cp Sweden_Metilene_random/Sweden_CHH_random_W_C_150_5_4_0.9_qval.1e-5.bedgraph Sweden_CHH_random_W_C.bedGraph

#Bedtools intersect all 4
module load StdEnv/2023 bedtools/2.31.0
bedtools intersect -u -a Nunavut_CHH_random_W_C.bedGraph -b Svalbard_CHH_random_W_C.bedGraph| bedtools intersect -u -a stdin -b Sweden_CHH_random_W_C.bedGraph | bedtools intersect -u -a stdin -b Alaska_CHH_random_W_C.bedGraph > ALL_Sites_intersect_random_DMRs.bedgraph

#Bedtools intersect all but Nunavut
bedtools intersect -u -a Svalbard_CHH_random_W_C.bedGraph -b Sweden_CHH_random_W_C.bedGraph | bedtools intersect -u -a stdin -b Alaska_CHH_random_W_C.bedGraph > SVAL_SWED_ALAS_Sites_intersect_random_DMRs.bedgraph
# Do1_01_a00001   7699414 7699600 -8.213185
# Do1_02_a00001   8510664 8510827 12.079776
# Do1_02_a00004   6400235 6400478 8.718939
# Do1_04_a00005   405319  405494  -5.244285
# Do1_06_a00001   1983414 1983677 4.578070
# Do1_06_a00001   6538326 6538612 4.997452

#Bedtools intersect all but Svalbard
bedtools intersect -u -a Nunavut_CHH_random_W_C.bedGraph -b Sweden_CHH_random_W_C.bedGraph | bedtools intersect -u -a stdin -b Alaska_CHH_random_W_C.bedGraph > ALEX_SWED_ALAS_Sites_intersect_random_DMRs.bedgraph

#Bedtools intersect all but Sweden
bedtools intersect -u -a Nunavut_CHH_random_W_C.bedGraph -b Svalbard_CHH_random_W_C.bedGraph| bedtools intersect -u -a stdin -b Alaska_CHH_random_W_C.bedGraph > ALEX_SVAL_ALAS_Sites_intersect_random_DMRs.bedgraph

#Bedtools intersect all but Alaska
bedtools intersect -u -a Nunavut_CHH_random_W_C.bedGraph -b Svalbard_CHH_random_W_C.bedGraph| bedtools intersect -u -a stdin -b Sweden_CHH_random_W_C.bedGraph > ALEX_SVAL_SWED_Sites_intersect_random_DMRs.bedgraph

wc -l ALL_Sites_intersect_DMRs.bedgraph
# 6
wc -l SVAL_SWED_ALAS_Sites_intersect_DMRs.bedgraph
# 37
wc -l ALEX_SWED_ALAS_Sites_intersect_DMRs.bedgraph
# 19
wc -l ALEX_SVAL_ALAS_Sites_intersect_DMRs.bedgraph
#20/ 60 - no pattern related to warming

wc -l ALEX_SVAL_SWED_Sites_intersect_DMRs.bedgraph
# 26
#----------------------------
# Are the same regions showing up?

bedtools intersect -u -a ALEX_SVAL_ALAS_Sites_intersect_random_DMRs.bedgraph -b ALEX_SVAL_ALAS_Sites_intersect_DMRs_CHH.bedGraph | wc -l 
# 0
bedtools intersect -u -a Alaska_CHH_random_W_C.bedGraph -b Alaska_CHH_W_C.bedGraph | wc -l 
#236
bedtools intersect -u -a Svalbard_CHH_random_W_C.bedGraph -b Svalbard_CHH_W_C.bedGraph | wc -l 
#269
bedtools intersect -u -a Nunavut_CHH_random_W_C.bedGraph -b Nunavut_CHH_W_C.bedGraph | wc -l 
# 142
bedtools intersect -u -a Sweden_CHH_random_W_C.bedGraph -b Sweden_CHH_W_C.bedGraph | wc -l 


###############################
# Phenology
cd /home/celphin/scratch/Dryas/CHG_CHH_random
mkdir Phenology_Metilene
cd Phenology_Metilene
cp /home/celphin/scratch/Dryas/MS_scripts/metilene*.sh .

tmux new-session -s Phenology_DMRS
tmux attach-session -t Phenology_DMRS



#######################
# Seedlings

cd /home/celphin/scratch/Dryas/CHG_CHH
mkdir Seedling_Metilene
cd Seedling_Metilene
cp /home/celphin/scratch/Dryas/MS_scripts/metilene*.sh .

tmux new-session -s Seedling_Warming_DMRS
tmux attach-session -t Seedling_Warming_DMRS

cd /home/celphin/scratch/Dryas/CHG_CHH/Seedling_Metilene

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
# Exploring non-random with more filtering

sh metilene_filter_qval.sh 150 5 4 10 1e-5
# Phenology 11665
# Svalbard 1020

sh metilene_filter_qval.sh 150 5 4 25 1e-5
# Svalbard 169 DMRs

sh metilene_filter_qval.sh 150 5 4 40 1e-5
# 22 Svalbard