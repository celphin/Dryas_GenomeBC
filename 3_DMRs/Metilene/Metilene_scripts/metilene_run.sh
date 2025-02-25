#metilene_run.sh: script to prepare files for metilene
################################################
#Input needed:
    #parameters:
        #1: maxdist
        #2: mincpgs
        #3: mindiff
    #h1, h2, metilenedir, inputdir
    # .input file from unionbedg before
#Possible adjustments to script: 
    # inmetilene format
    # outputmetilene format
#Output:
    # in_metilene.input (name specifed in parameters)

###############################################
h1='W'
h2='C'
metilene_dir=/home/msandler/projects/rpp-rieseber/msandler/Dryas/Dryas_metilene_run
input_dir="/home/msandler/scratch/Seedling_Metilene/SE_${h1}_${h2}_input_files"

maxdist=$1
mincpgs=$2
mindiff=$3

output_name=SE_"$h1"_"$h2"_"${maxdist}"_"${mincpgs}"_"${mindiff}"
in_metilene="SE_metilene_"$h1"_"$h2".input"
threads=32


$metilene_dir/metilene_v0.2-8/metilene \
"${input_dir}/${in_metilene}" \--maxdist ${maxdist} --mincpgs ${mincpgs} --minMethDiff ${mindiff} --mode 1 --threads ${threads} -a ${h1} -b ${h2} -v 0.2 > \
"${output_name}"
