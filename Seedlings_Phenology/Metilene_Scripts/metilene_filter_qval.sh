#metilene_filter_qval.sh: script to prepare files for metilene
################################################
#Input needed:
    #parameters:
        #1: maxdist
        #2: mincpgs
        #3: mindiff
        #4: minmeandiff
        #5: q-value to use
    #h1, h2, metilenedir, inputdir
    #finished metilene output file, matching outputname format
#Possible adjustments to script: 
    # inmetilene format
    # outputmetilene format
#Output:
    # DMRS, filtered by q-value

###############################################
maxdist=$1
mincpgs=$2
mindiff=$3
minmeandiff=$4
qval=$5

h1='W'
h2='C'
metilene_dir=/home/msandler/projects/rpp-rieseber/msandler/Dryas/Dryas_metilene_run
input_dir="/home/msandler/scratch/Seedling_Metilene/SE_W_C_input_files"
outputname=SE_"$h1"_"$h2"_"${maxdist}"_"${mincpgs}"_"${mindiff}"

mincpgs=10

perl $metilene_dir/metilene_v0.2-8/metilene_output.pl \
-q "${outputname}" -p ${qval} -c ${mincpgs} -d ${minmeandiff} -a ${h1} -b ${h2} \
-o "${outputname}_${minmeandiff}"
