#########################################################
# Dryas DMR calling
# Dec 2024
###############################
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

##########################
# Use data on Cedar
cd /home/celphin/projects/def-rieseber/Dryas_shared_data/Dryas/MS_Metilene_Output

#########################
# All wild warming and control

cd /home/celphin/projects/def-rieseber/Dryas_shared_data/Dryas/MS_Metilene_Output/

mkdir Wild_W_C_Metilene_random
cd Wild_W_C_Metilene_random
cp /home/celphin/projects/def-rieseber/Dryas_shared_data/Dryas/MS_Metilene_Output/Wild_W_C_Metilene/* /home/celphin/projects/def-rieseber/Dryas_shared_data/Dryas/MS_Metilene_Output/Wild_W_C_Metilene_random/

tmux new-session -s Wild_DMRS
tmux attach-session -t Wild_DMRS
salloc -c15 --time 7:00:00 --mem 200000m --account def-cronk

# Define the path to the input file 
cd /home/celphin/projects/def-rieseber/Dryas_shared_data/Dryas/MS_Metilene_Output/Wild_W_C_Metilene_random/

cp ../Wild_W_C_Metilene/Wild_W_C_input_files/*.input .

head Wild_metilene_W_C.input > Wild_metilene_W_C1.input

# check
head -n 1 Wild_metilene_W_C1.input

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
' Wild_metilene_W_C.input > Wild_metilene_W_C2.input

#-------------------------------------
# replace NA with . or - 

sed 's/NA/-/g' Wild_metilene_W_C2.input > Wild_metilene_W_C3.input

#------------------------------------------------------------------

#Wild specific adjustments,
nano metilene_run.sh

h1='W'
h2='C'
metilene_dir="/home/celphin/projects/rpp-rieseber/celphin/Dryas/DMRs/June2022_Metilene_DMRs"
input_dir="/home/celphin/projects/def-rieseber/Dryas_shared_data/Dryas/MS_Metilene_Output/Wild_W_C_Metilene_random"
maxdist=$1
mincpgs=$2
mindiff=$3

output_name=Wild_random_"$h1"_"$h2"_"${maxdist}"_"${mincpgs}"_"${mindiff}"
in_metilene="Wild_metilene_"$h1"_"$h2"2.input"
threads=15


#---------------------------------------
#metilene:
# params: maxdist, mincpgs, mindiff
salloc -c15 --time 3:00:00 --mem 192000m --account rrg-rieseber-ac

module load StdEnv/2020
module load bedtools/2.30.0 

# with in_metilene="Wild_metilene_"$h1"_"$h2"2.input"
sh metilene_run.sh 150 5 4 


#-------------------------------------------------
#Filter based on qval
#params: maxdist, mincpgs, mindiff, minmeandif, q-value

nano  metilene_filter_qval.sh

metilene_dir="/home/celphin/scratch/Dryas"
input_dir="/home/celphin/projects/def-rieseber/Dryas_shared_data/Dryas/MS_Metilene_Output/Wild_W_C_Metilene_random"
outputname=Wild_random_"$h1"_"$h2"_"${maxdist}"_"${mincpgs}"_"${mindiff}"

# to run
module load StdEnv/2023
module load r/4.4.0
module load gdal
module load udunits
module load python
export R_LIBS_USER=/home/celphin/R/x86_64-pc-linux-gnu-library/4.4/

sh metilene_filter_qval.sh 150 5 4 25 1e-5
# Wild 490 DMRs with 0.9
# Wild 1 DMR with 25
# Wild 237 with 10

# real 109 with 10
# real 0 with 25



####################################################################################
#Svalbard:

cd /lustre04/scratch/celphin/Dryas_large_folders/CpG/Metilene/

mkdir Svalbard_Metilene_random
cd Svalbard_Metilene_random
cp /lustre04/scratch/celphin/Dryas_large_folders/CpG/Metilene/Svalbard_Metilene/*sh /lustre04/scratch/celphin/Dryas_large_folders/CpG/Metilene/Svalbard_Metilene_random/

tmux new-session -s Svalbard_DMRS
tmux attach-session -t Svalbard_DMRS
salloc -c15 --time 7:00:00 --mem 200000m --account def-henryg

# Define the path to the input file 
cd /lustre04/scratch/celphin/Dryas_large_folders/CpG/Metilene/Svalbard_Metilene_random/

cp ../Svalbard_Metilene/Svalbard_W_C_input_files/Svalbard_metilene_W_C.input .

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
metilene_dir="/lustre04/scratch/celphin/Dryas"
input_dir="/lustre04/scratch/celphin/Dryas_large_folders/CpG/Metilene/Svalbard_Metilene_random"
...
output_name=Svalbard_random_"$h1"_"$h2"_"${maxdist}"_"${mincpgs}"_"${mindiff}"

threads=15


#metilene:
# params: maxdist, mincpgs, mindiff
salloc -c7 --time 7:00:00 --mem 200000m --account def-cronk

module load StdEnv/2020
module load bedtools/2.30.0 
cd /lustre04/scratch/celphin/Dryas_large_folders/CpG/Metilene/Svalbard_Metilene_random
# with in_metilene="Svalbard_metilene_"$h1"_"$h2"2.input"
sh metilene_run.sh 150 5 4 


#-------------------------------------------------
#Filter based on qval
#params: maxdist, mincpgs, mindiff, minmeandif, q-value

nano  metilene_filter_qval.sh

maxdist=$1
mincpgs=$2
mindiff=$3
minmeandiff=$4
qval=$5

h1='W'
h2='C'
metilene_dir="/lustre04/scratch/celphin/Dryas"
input_dir="/lustre04/scratch/celphin/Dryas_large_folders/CpG/Metilene/Svalbard_Metilene_random"
outputname=Svalbard_random_"$h1"_"$h2"_"${maxdist}"_"${mincpgs}"_"${mindiff}"

mincpgs=10

perl $metilene_dir/metilene_v0.2-8/metilene_output.pl \
-q "${outputname}" -p ${qval} -c ${mincpgs} -d ${minmeandiff} -a ${h1} -b ${h2} \
-o "${outputname}_${minmeandiff}"

#-------------
# to run
module load StdEnv/2023
module load r/4.4.0
module load gdal
module load udunits
module load python
export R_LIBS_USER=/home/celphin/R/x86_64-pc-linux-gnu-library/4.4/

sh metilene_filter_qval.sh 150 5 4 0.9 1e-5
# Wrote 386 DMRs 

###############################################################
# Sweden
cd /lustre04/scratch/celphin/Dryas_large_folders/CpG/Metilene/

mkdir Sweden_Metilene_random
cd Sweden_Metilene_random
cp /lustre04/scratch/celphin/Dryas_large_folders/CpG/Metilene/Sweden_Metilene/*sh /lustre04/scratch/celphin/Dryas_large_folders/CpG/Metilene/Sweden_Metilene_random/

tmux new-session -s Sweden_DMRS
tmux attach-session -t Sweden_DMRS
salloc -c15 --time 7:00:00 --mem 200000m --account def-henryg

# Define the path to the input file 
cd /lustre04/scratch/celphin/Dryas_large_folders/CpG/Metilene/Sweden_Metilene_random/

cp ../Sweden_Metilene/Sweden_W_C_input_files/Sweden_metilene_W_C.input .

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
' Sweden_metilene_W_C.input > Sweden_metilene_W_C4.input
# ' Sweden_metilene_W_C.input > Sweden_metilene_W_C2.input

#-------------------------------------
# replace NA with . or - 

sed 's/NA/-/g' Sweden_metilene_W_C2.input > Sweden_metilene_W_C3.input

#------------------------------------------------------------------

#Sweden specific adjustments,
nano metilene_run.sh

h1='W'
h2='C'
metilene_dir="/lustre04/scratch/celphin/Dryas"
input_dir="/lustre04/scratch/celphin/Dryas_large_folders/CpG/Metilene/Sweden_Metilene_random"

maxdist=$1
mincpgs=$2
mindiff=$3

output_name=Sweden_random2_"$h1"_"$h2"_"${maxdist}"_"${mincpgs}"_"${mindiff}"
in_metilene="Sweden_metilene_"$h1"_"$h2"4.input"
threads=15



#metilene:
# params: maxdist, mincpgs, mindiff
salloc -c15 --time 7:00:00 --mem 200000m --account def-cronk

module load StdEnv/2020
module load bedtools/2.30.0 
cd /lustre04/scratch/celphin/Dryas_large_folders/CpG/Metilene/Sweden_Metilene_random

# with in_metilene="Sweden_metilene_"$h1"_"$h2"2.input"
sh metilene_run.sh 150 5 4 


#-------------------------------------------------
#Filter based on qval
#params: maxdist, mincpgs, mindiff, minmeandif, q-value

nano  metilene_filter_qval.sh

maxdist=$1
mincpgs=$2
mindiff=$3
minmeandiff=$4
qval=$5

h1='W'
h2='C'
metilene_dir="/lustre04/scratch/celphin/Dryas"
input_dir="/lustre04/scratch/celphin/Dryas_large_folders/CpG/Metilene/Sweden_Metilene_random"
outputname=Sweden_random2_"$h1"_"$h2"_"${maxdist}"_"${mincpgs}"_"${mindiff}"

mincpgs=10

perl $metilene_dir/metilene_v0.2-8/metilene_output.pl \
-q "${outputname}" -p ${qval} -c ${mincpgs} -d ${minmeandiff} -a ${h1} -b ${h2} \
-o "${outputname}_${minmeandiff}"

#--------
# to run
module load StdEnv/2023
module load r/4.4.0
module load gdal
module load udunits
module load python
export R_LIBS_USER=/home/celphin/R/x86_64-pc-linux-gnu-library/4.4/

sh metilene_filter_qval.sh 150 5 4 0.9 1e-5
# Wrote 429 DMRs
# Wrote 388 DMRs

###############################################################
# Alaska
cd /lustre04/scratch/celphin/Dryas_large_folders/CpG/Metilene/

mkdir Alaska_Metilene_random
cd Alaska_Metilene_random
cp /lustre04/scratch/celphin/Dryas_large_folders/CpG/Metilene/Alaska_Metilene/*.sh /lustre04/scratch/celphin/Dryas_large_folders/CpG/Metilene/Alaska_Metilene_random/

tmux new-session -s Alaska_DMRS
tmux attach-session -t Alaska_DMRS
salloc -c15 --time 7:00:00 --mem 200000m --account def-henryg

# Define the path to the input file 
cd /lustre04/scratch/celphin/Dryas_large_folders/CpG/Metilene/Alaska_Metilene_random/

cp ../Alaska_Metilene/Alaska_W_C_input_files/Alaska_metilene_W_C.input .

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
' Alaska_metilene_W_C.input > Alaska_metilene_W_C4.input
#' Alaska_metilene_W_C.input > Alaska_metilene_W_C2.input

#-------------------------------------
# replace NA with . or - 

sed 's/NA/-/g' Alaska_metilene_W_C2.input > Alaska_metilene_W_C3.input

#------------------------------------------------------------------

#Alaska specific adjustments,
nano metilene_run.sh

h1='W'
h2='C'
metilene_dir="/lustre04/scratch/celphin/Dryas"
input_dir="/lustre04/scratch/celphin/Dryas_large_folders/CpG/Metilene/Alaska_Metilene_random"
maxdist=$1
mincpgs=$2
mindiff=$3

output_name=Alaska_random2_"$h1"_"$h2"_"${maxdist}"_"${mincpgs}"_"${mindiff}"
in_metilene="Alaska_metilene_"$h1"_"$h2"4.input"
threads=15

$metilene_dir/metilene_v0.2-8/metilene \
"${input_dir}/${in_metilene}" \--maxdist ${maxdist} --mincpgs ${mincpgs} --minMethDiff ${mindiff} ->
"${output_name}"



#metilene:
# params: maxdist, mincpgs, mindiff
salloc -c15 --time 7:00:00 --mem 200000m --account def-cronk

module load StdEnv/2020
module load bedtools/2.30.0 
cd /lustre04/scratch/celphin/Dryas_large_folders/CpG/Metilene/Alaska_Metilene_random
# with in_metilene="Alaska_metilene_"$h1"_"$h2"2.input"
sh metilene_run.sh 150 5 4 




#-------------------------------------------------
#Filter based on qval
#params: maxdist, mincpgs, mindiff, minmeandif, q-value

nano  metilene_filter_qval.sh

maxdist=$1
mincpgs=$2
mindiff=$3
minmeandiff=$4
qval=$5

h1='W'
h2='C'
metilene_dir="/lustre04/scratch/celphin/Dryas"
input_dir="/lustre04/scratch/celphin/Dryas_large_folders/CpG/Metilene/Alaska_Metilene_random"
outputname=Alaska_random2_"$h1"_"$h2"_"${maxdist}"_"${mincpgs}"_"${mindiff}"

mincpgs=10

perl $metilene_dir/metilene_v0.2-8/metilene_output.pl \
-q "${outputname}" -p ${qval} -c ${mincpgs} -d ${minmeandiff} -a ${h1} -b ${h2} \
-o "${outputname}_${minmeandiff}"

#-------
# to run
module load StdEnv/2023
module load r/4.4.0
module load gdal
module load udunits
module load python
export R_LIBS_USER=/home/celphin/R/x86_64-pc-linux-gnu-library/4.4/

sh metilene_filter_qval.sh 150 5 4 0.9 1e-5
# Wrote 392 DMRs
# Wrote 408 DMRs

#########################
# Nunavut
cd /lustre04/scratch/celphin/Dryas_large_folders/CpG/Metilene/

mkdir Nunavut_Metilene_random
cd Nunavut_Metilene_random
cp /lustre04/scratch/celphin/Dryas_large_folders/CpG/Metilene/Nunavut_Metilene/* /lustre04/scratch/celphin/Dryas_large_folders/CpG/Metilene/Nunavut_Metilene_random/

tmux new-session -s Nunavut_DMRS
tmux attach-session -t Nunavut_DMRS
salloc -c15 --time 7:00:00 --mem 200000m --account def-henryg

# Define the path to the input file 
cd /lustre04/scratch/celphin/Dryas_large_folders/CpG/Metilene/Nunavut_Metilene_random/

cp ../Nunavut_Metilene/Nunavut_W_C_input_files/Nunavut_metilene_W_C.input .

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
' Nunavut_metilene_W_C.input > Nunavut_metilene_W_C4.input

#' Nunavut_metilene_W_C.input > Nunavut_metilene_W_C2.input

#-------------------------------------
# replace NA with . or - 

sed 's/NA/-/g' Nunavut_metilene_W_C2.input > Nunavut_metilene_W_C3.input

#------------------------------------------------------------------

#Nunavut specific adjustments,
nano metilene_run.sh

h1='W'
h2='C'
metilene_dir="/lustre04/scratch/celphin/Dryas"
input_dir="/lustre04/scratch/celphin/Dryas_large_folders/CpG/Metilene/Nunavut_Metilene_random"
maxdist=$1
mincpgs=$2
mindiff=$3

output_name=Nunavut_random2_"$h1"_"$h2"_"${maxdist}"_"${mincpgs}"_"${mindiff}"
in_metilene="Nunavut_metilene_"$h1"_"$h2"4.input"
threads=15


#---------------------------------------
#metilene:
# params: maxdist, mincpgs, mindiff
salloc -c15 --time 7:00:00 --mem 200000m --account def-henryg

module load StdEnv/2020
module load bedtools/2.30.0 

# with in_metilene="Nunavut_metilene_"$h1"_"$h2"2.input"
sh metilene_run.sh 150 5 4 


#-------------------------------------------------
#Filter based on qval
#params: maxdist, mincpgs, mindiff, minmeandif, q-value

nano  metilene_filter_qval.sh

maxdist=$1
mincpgs=$2
mindiff=$3
minmeandiff=$4
qval=$5

h1='W'
h2='C'
metilene_dir="/lustre04/scratch/celphin/Dryas"
input_dir="/lustre04/scratch/celphin/Dryas_large_folders/CpG/Metilene/Nunavut_Metilene_random"
outputname=Nunavut_random2_"$h1"_"$h2"_"${maxdist}"_"${mincpgs}"_"${mindiff}"

mincpgs=10

perl $metilene_dir/metilene_v0.2-8/metilene_output.pl \
-q "${outputname}" -p ${qval} -c ${mincpgs} -d ${minmeandiff} -a ${h1} -b ${h2} \
-o "${outputname}_${minmeandiff}"

#-------------------------
# to run
module load StdEnv/2023
module load r/4.4.0
module load gdal
module load udunits
module load python
export R_LIBS_USER=/home/celphin/R/x86_64-pc-linux-gnu-library/4.4/

sh metilene_filter_qval.sh 150 5 4 0.9 1e-5
# Wrote 291 DMRs
# Wrote 309 DMRs

##############################################################
mkdir /lustre04/scratch/celphin/Dryas_large_folders/CpG/Metilene/MS_Wild_metilene_output_bedgraphs/Random_Site_bedgraphs
cd /lustre04/scratch/celphin/Dryas_large_folders/CpG/Metilene/MS_Wild_metilene_output_bedgraphs/Random_Site_bedgraphs

#Copy bedGraph files to directory above
cp /lustre04/scratch/celphin/Dryas_large_folders/CpG/Metilene/Nunavut_Metilene_random/Nunavut_random2_W_C_150_5_4_0.9_qval.1e-5.bedgraph ./Nunavut_random_W_C.bedGraph
cp /lustre04/scratch/celphin/Dryas_large_folders/CpG/Metilene/Alaska_Metilene_random/Alaska_random2_W_C_150_5_4_0.9_qval.1e-5.bedgraph Alaska_random_W_C.bedGraph
cp /lustre04/scratch/celphin/Dryas_large_folders/CpG/Metilene/Svalbard_Metilene_random/Svalbard_random_W_C_150_5_4_0.9_qval.1e-5.bedgraph Svalbard_random_W_C.bedGraph
cp /lustre04/scratch/celphin/Dryas_large_folders/CpG/Metilene/Sweden_Metilene_random/Sweden_random_W_C_150_5_4_0.9_qval.1e-5.bedgraph Sweden_random_W_C.bedGraph

#Bedtools intersect all 4
module load StdEnv/2023 bedtools/2.31.0
bedtools intersect -u -a Nunavut_random_W_C.bedGraph -b Svalbard_random_W_C.bedGraph| bedtools intersect -u -a stdin -b Sweden_random_W_C.bedGraph | bedtools intersect -u -a stdin -b Alaska_random_W_C.bedGraph > ALL_Sites_intersect_DMRs.bedgraph

#Bedtools intersect all but Nunavut
bedtools intersect -u -a Svalbard_random_W_C.bedGraph -b Sweden_random_W_C.bedGraph | bedtools intersect -u -a stdin -b Alaska_random_W_C.bedGraph > SVAL_SWED_ALAS_Sites_intersect_DMRs.bedgraph
# Do1_01_a00001   7699414 7699600 -8.213185
# Do1_02_a00001   8510664 8510827 12.079776
# Do1_02_a00004   6400235 6400478 8.718939
# Do1_04_a00005   405319  405494  -5.244285
# Do1_06_a00001   1983414 1983677 4.578070
# Do1_06_a00001   6538326 6538612 4.997452

#Bedtools intersect all but Svalbard
bedtools intersect -u -a Nunavut_random_W_C.bedGraph -b Sweden_random_W_C.bedGraph | bedtools intersect -u -a stdin -b Alaska_random_W_C.bedGraph > ALEX_SWED_ALAS_Sites_intersect_DMRs.bedgraph

#Bedtools intersect all but Sweden
bedtools intersect -u -a Nunavut_random_W_C.bedGraph -b Svalbard_random_W_C.bedGraph| bedtools intersect -u -a stdin -b Alaska_random_W_C.bedGraph > ALEX_SVAL_ALAS_Sites_intersect_DMRs.bedgraph

#Bedtools intersect all but Alaska
bedtools intersect -u -a Nunavut_random_W_C.bedGraph -b Svalbard_random_W_C.bedGraph| bedtools intersect -u -a stdin -b Sweden_random_W_C.bedGraph > ALEX_SVAL_SWED_Sites_intersect_DMRs.bedgraph

# Random counts
wc -l *intersect_DMRs.bedgraph
  # 4 ALEX_SVAL_ALAS_Sites_intersect_DMRs.bedgraph
  # 2 ALEX_SVAL_SWED_Sites_intersect_DMRs.bedgraph
  # 4 ALEX_SWED_ALAS_Sites_intersect_DMRs.bedgraph
  # 1 ALL_Sites_intersect_DMRs.bedgraph
  # 4 SVAL_SWED_ALAS_Sites_intersect_DMRs.bedgraph
 # 15 total

wc -l *_W_C.bedGraph
  # 408 Alaska_random_W_C.bedGraph
  # 309 Nunavut_random_W_C.bedGraph
  # 386 Svalbard_random_W_C.bedGraph
  # 429 Sweden_random_W_C.bedGraph
 # 1532 total

grep "-" Sweden_random_W_C.bedGraph |wc -l
grep -v "-" Sweden_random_W_C.bedGraph |wc -l
# 241 188

grep "-" Alaska_random_W_C.bedGraph |wc -l
grep -v "-" Alaska_random_W_C.bedGraph |wc -l
# 201 207

grep "-" Svalbard_random_W_C.bedGraph |wc -l
grep -v "-" Svalbard_random_W_C.bedGraph |wc -l
# 205 181

grep "-" Nunavut_random_W_C.bedGraph |wc -l
grep -v "-" Nunavut_random_W_C.bedGraph |wc -l
# 158 151

#############################
# Non-Random

cd /lustre04/scratch/celphin/Dryas_large_folders/CpG/Metilene/MS_Wild_metilene_output_bedgraphs/Site_bedgraphs

#Bedtools intersect all 4
module load StdEnv/2023 bedtools/2.31.0
bedtools intersect -u -a Nunavut_W_C.bedGraph -b Svalbard_W_C.bedGraph| bedtools intersect -u -a stdin -b Sweden_W_C.bedGraph | bedtools intersect -u -a stdin -b Alaska_W_C.bedGraph > ALL_Sites_intersect_DMRs.bedgraph

#Bedtools intersect all but Nunavut
bedtools intersect -u -a Svalbard_W_C.bedGraph -b Sweden_W_C.bedGraph | bedtools intersect -u -a stdin -b Alaska_W_C.bedGraph > SVAL_SWED_ALAS_Sites_intersect_DMRs.bedgraph
# Do1_01_a00001   7699414 7699600 -8.213185
# Do1_02_a00001   8510664 8510827 12.079776
# Do1_02_a00004   6400235 6400478 8.718939
# Do1_04_a00005   405319  405494  -5.244285
# Do1_06_a00001   1983414 1983677 4.578070
# Do1_06_a00001   6538326 6538612 4.997452

#Bedtools intersect all but Svalbard
bedtools intersect -u -a Nunavut_W_C.bedGraph -b Sweden_W_C.bedGraph | bedtools intersect -u -a stdin -b Alaska_W_C.bedGraph > ALEX_SWED_ALAS_Sites_intersect_DMRs.bedgraph

#Bedtools intersect all but Sweden
bedtools intersect -u -a Nunavut_W_C.bedGraph -b Svalbard_W_C.bedGraph| bedtools intersect -u -a stdin -b Alaska_W_C.bedGraph > ALEX_SVAL_ALAS_Sites_intersect_DMRs.bedgraph

#Bedtools intersect all but Alaska
bedtools intersect -u -a Nunavut_W_C.bedGraph -b Svalbard_W_C.bedGraph| bedtools intersect -u -a stdin -b Sweden_W_C.bedGraph > ALEX_SVAL_SWED_Sites_intersect_DMRs.bedgraph


wc -l *intersect_DMRs.bedgraph
   # 7 ALEX_SVAL_ALAS_Sites_intersect_DMRs.bedgraph
   # 9 ALEX_SVAL_SWED_Sites_intersect_DMRs.bedgraph
  # 19 ALEX_SWED_ALAS_Sites_intersect_DMRs.bedgraph
   # 4 ALL_Sites_intersect_DMRs.bedgraph
  # 15 SVAL_SWED_ALAS_Sites_intersect_DMRs.bedgraph
  # 54 total

wc -l *_W_C.bedGraph
   # 866 Alaska_W_C.bedGraph
   # 322 Nunavut_W_C.bedGraph
   # 358 Svalbard_W_C.bedGraph
  # 1737 Sweden_W_C.bedGraph

#------------
# Count hypo/hyper

grep "-" Sweden_W_C.bedGraph |wc -l
grep -v "-" Sweden_W_C.bedGraph |wc -l
# 958 779

grep "-" Alaska_W_C.bedGraph |wc -l
grep -v "-" Alaska_W_C.bedGraph |wc -l
# 471 395

grep "-" Svalbard_W_C.bedGraph |wc -l
grep -v "-" Svalbard_W_C.bedGraph |wc -l
#  158 200

grep "-" Nunavut_W_C.bedGraph |wc -l
grep -v "-" Nunavut_W_C.bedGraph |wc -l
# 169 153


###########################
# CHG/CHH 
cd /lustre04/scratch/celphin/Dryas_large_folders/CHG_CHH/Metilene


#Bedtools intersect all 4
module load StdEnv/2023 bedtools/2.31.0
bedtools intersect -u -a Nunavut_CHH_W_C.bedGraph -b Svalbard_CHH_W_C.bedGraph| bedtools intersect -u -a stdin -b Sweden_CHH_W_C.bedGraph | bedtools intersect -u -a stdin -b Alaska_CHH_W_C.bedGraph > ALL_Sites_intersect_DMRs_CHH.bedgraph

#Bedtools intersect all but Nunavut
bedtools intersect -u -a Svalbard_CHH_W_C.bedGraph -b Sweden_CHH_W_C.bedGraph | bedtools intersect -u -a stdin -b Alaska_CHH_W_C.bedGraph > SVAL_SWED_ALAS_Sites_intersect_DMRs_CHH.bedgraph
# Do1_01_a00001   7699414 7699600 -8.213185
# Do1_02_a00001   8510664 8510827 12.079776
# Do1_02_a00004   6400235 6400478 8.718939
# Do1_04_a00005   405319  405494  -5.244285
# Do1_06_a00001   1983414 1983677 4.578070
# Do1_06_a00001   6538326 6538612 4.997452

#Bedtools intersect all but Svalbard
bedtools intersect -u -a Nunavut_CHH_W_C.bedGraph -b Sweden_CHH_W_C.bedGraph | bedtools intersect -u -a stdin -b Alaska_CHH_W_C.bedGraph > ALEX_SWED_ALAS_Sites_intersect_DMRs_CHH.bedgraph

#Bedtools intersect all but Sweden
bedtools intersect -u -a Nunavut_CHH_W_C.bedGraph -b Svalbard_CHH_W_C.bedGraph| bedtools intersect -u -a stdin -b Alaska_CHH_W_C.bedGraph > ALEX_SVAL_ALAS_Sites_intersect_DMRs_CHH.bedgraph

#Bedtools intersect all but Alaska
bedtools intersect -u -a Nunavut_CHH_W_C.bedGraph -b Svalbard_CHH_W_C.bedGraph| bedtools intersect -u -a stdin -b Sweden_CHH_W_C.bedGraph > ALEX_SVAL_SWED_Sites_intersect_DMRs_CHH.bedgraph

# Random counts
wc -l *intersect_DMRs_CHH.bedgraph
  # 20 ALEX_SVAL_ALAS_Sites_intersect_DMRs_CHH.bedgraph
  # 26 ALEX_SVAL_SWED_Sites_intersect_DMRs_CHH.bedgraph
  # 19 ALEX_SWED_ALAS_Sites_intersect_DMRs_CHH.bedgraph
   # 6 ALL_Sites_intersect_DMRs_CHH.bedgraph
  # 37 SVAL_SWED_ALAS_Sites_intersect_DMRs_CHH.bedgraph
 # 108 total


wc -l *CHH_W_C.bedGraph
 # 1937 Alaska_CHH_W_C.bedGraph
  # 1186 Nunavut_CHH_W_C.bedGraph
  # 2356 SE_CHH_W_C.bedGraph
  # 1397 Svalbard_CHH_W_C.bedGraph
  # 2206 Sweden_CHH_W_C.bedGraph
  # 9082 total


grep "-" Sweden_CHH_W_C.bedGraph |wc -l
grep -v "-" Sweden_CHH_W_C.bedGraph |wc -l
# 1138 1068

grep "-" Alaska_CHH_W_C.bedGraph |wc -l
grep -v "-" Alaska_CHH_W_C.bedGraph |wc -l
# 436 1501

grep "-" Svalbard_CHH_W_C.bedGraph |wc -l
grep -v "-" Svalbard_CHH_W_C.bedGraph |wc -l
# 725 672

grep "-" Nunavut_CHH_W_C.bedGraph |wc -l
grep -v "-" Nunavut_CHH_W_C.bedGraph |wc -l
# 320 866 

#########################
# Random CHH

#Bedtools intersect all 4
module load StdEnv/2023 bedtools/2.31.0
bedtools intersect -u -a Nunavut_CHH_random_W_C.bedGraph -b Svalbard_CHH_random_W_C.bedGraph| bedtools intersect -u -a stdin -b Sweden_CHH_random_W_C.bedGraph | bedtools intersect -u -a stdin -b Alaska_CHH_random_W_C.bedGraph > ALL_Sites_intersect_DMRs_CHH_random.bedgraph

#Bedtools intersect all but Nunavut
bedtools intersect -u -a Svalbard_CHH_random_W_C.bedGraph -b Sweden_CHH_random_W_C.bedGraph | bedtools intersect -u -a stdin -b Alaska_CHH_random_W_C.bedGraph > SVAL_SWED_ALAS_Sites_intersect_DMRs_CHH_random.bedgraph
# Do1_01_a00001   7699414 7699600 -8.213185
# Do1_02_a00001   8510664 8510827 12.079776
# Do1_02_a00004   6400235 6400478 8.718939
# Do1_04_a00005   405319  405494  -5.244285
# Do1_06_a00001   1983414 1983677 4.578070
# Do1_06_a00001   6538326 6538612 4.997452

#Bedtools intersect all but Svalbard
bedtools intersect -u -a Nunavut_CHH_random_W_C.bedGraph -b Sweden_CHH_random_W_C.bedGraph | bedtools intersect -u -a stdin -b Alaska_CHH_random_W_C.bedGraph > ALEX_SWED_ALAS_Sites_intersect_DMRs_CHH_random.bedgraph

#Bedtools intersect all but Sweden
bedtools intersect -u -a Nunavut_CHH_random_W_C.bedGraph -b Svalbard_CHH_random_W_C.bedGraph| bedtools intersect -u -a stdin -b Alaska_CHH_random_W_C.bedGraph > ALEX_SVAL_ALAS_Sites_intersect_DMRs_CHH_random.bedgraph

#Bedtools intersect all but Alaska
bedtools intersect -u -a Nunavut_CHH_random_W_C.bedGraph -b Svalbard_CHH_random_W_C.bedGraph| bedtools intersect -u -a stdin -b Sweden_CHH_random_W_C.bedGraph > ALEX_SVAL_SWED_Sites_intersect_DMRs_CHH_random.bedgraph

# Random counts
wc -l *intersect_DMRs_CHH_random.bedgraph
  # 60 ALEX_SVAL_ALAS_Sites_intersect_DMRs_CHH_random.bedgraph
  # 66 ALEX_SVAL_SWED_Sites_intersect_DMRs_CHH_random.bedgraph
  # 43 ALEX_SWED_ALAS_Sites_intersect_DMRs_CHH_random.bedgraph
  # 17 ALL_Sites_intersect_DMRs_CHH_random.bedgraph
  # 67 SVAL_SWED_ALAS_Sites_intersect_DMRs_CHH_random.bedgraph
 # 253 total

wc -l *_CHH_random_W_C.bedGraph
  # 1780 Alaska_CHH_random_W_C.bedGraph
  # 1278 Nunavut_CHH_random_W_C.bedGraph
  # 2002 Svalbard_CHH_random_W_C.bedGraph
  # 1653 Sweden_CHH_random_W_C.bedGraph


grep "-" Sweden_CHH_random_W_C.bedGraph |wc -l
grep -v "-" Sweden_CHH_random_W_C.bedGraph |wc -l
# 799 854

grep "-" Alaska_CHH_random_W_C.bedGraph |wc -l
grep -v "-" Alaska_CHH_random_W_C.bedGraph |wc -l
# 861 919

grep "-" Svalbard_CHH_random_W_C.bedGraph |wc -l
grep -v "-" Svalbard_CHH_random_W_C.bedGraph |wc -l
# 994 1008

grep "-" Nunavut_CHH_random_W_C.bedGraph |wc -l
grep -v "-" Nunavut_CHH_random_W_C.bedGraph |wc -l
# 646 632
