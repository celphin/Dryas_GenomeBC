#########################################################
# Dryas DNA methylation calling
# nf-core Methylseq
# Mapping WGBS data to reference
# Calling methylated sites
#############################################################

#Installing nextflow/methylseq
# following instructions here: 
# https://nf-co.re/methylseq/2.7.0

# First Setup nextflow profile
# https://docs.alliancecan.ca/wiki/Nextflow

# install pip for setup
module purge # Make sure that previously loaded modules are not polluting the installation 
module load python/3.11
module load rust # New nf-core installations will err out if rust hasn't been loaded
module load postgresql # Will not use PostgresSQL here, but some Python modules which list psycopg2 as a dependency in the installation would crash without it.
python -m venv nf-core-env
source nf-core-env/bin/activate
python -m pip install nf_core==2.13
deactivate


module load nextflow/24.04.4
module load apptainer/1.3.4

# download images for nextflow generally
mkdir /project/rrg-rieseber-ac/NXF_SINGULARITY_CACHEDIR
export NXF_SINGULARITY_CACHEDIR=/project/rrg-rieseber-ac/NXF_SINGULARITY_CACHEDIR

#-----------------
# download the methylseq pipeline to scratch
# need allocation
tmux new-session -s Dryas
tmux attach-session -t Dryas

cd /home/celphin/scratch/Dryas/methylseq

salloc -c1 --time 23:00:00 --mem 120000m --account rrg-rieseber-ac

cd /home/celphin/scratch/Dryas/methylseq
source nf-core-env/bin/activate
module load nextflow/24.04.4
module load apptainer/1.3.4
export NXF_SINGULARITY_CACHEDIR=/project/rrg-rieseber-ac/NXF_SINGULARITY_CACHEDIR
# set name and version of pipeline
export NFCORE_PL=methylseq
export PL_VERSION=2.7.0

nf-core download --container-cache-utilisation amend --container-system  singularity   --compress none -r ${PL_VERSION}  -p 6  ${NFCORE_PL}

# check everything is there
ls nf-core-${NFCORE_PL}_${PL_VERSION}


#-------------------------------
# edit the nextflow config for CC
# https://docs.alliancecan.ca/wiki/Cedar
# https://docs.alliancecan.ca/wiki/Job_scheduling_policies

nano ~/.nextflow/config

params {
    config_profile_description = 'Alliance HPC config'
    config_profile_contact = 'support@alliancecan.ca'
    config_profile_url = 'docs.alliancecan.ca/mediawiki/index.php?title=Nextflow'
}

singularity {
  enabled = true
  autoMounts = true
}

apptainer {
  autoMounts = true
}

process {
  executor = 'slurm' 
  clusterOptions = '--account=rrg-rieseber-ac'
  maxRetries = 1
  errorStrategy = { task.exitStatus in [125,139] ? 'retry' : 'finish' }
  memory = { check_max( 4.GB * task.attempt, 'memory' ) }
  cpu = 1  
  time = '3h' 
}

executor {
  pollInterval = '60 sec'
  submitRateLimit = '60/1min'
  queueSize = 100 
}

profiles {
  cedar {
    max_memory='186G'
    max_cpu=48
    max_time='168h'
  }
  beluga {
    max_memory='186G'
    max_cpu=40
    max_time='168h'
  }
  narval {
    max_memory='249G'
    max_cpu=64
    max_time='168h'
  }
}


##################################
# Editing config files for setup

# installed here 
cd /home/celphin/scratch/Dryas/methylseq/nf-core-methylseq_2.7.0/2_7_0

# Note: need to copy the nf-core/methylseq/conf/base.config to the directory running the program
cp nextflow.config /home/celphin/scratch/Dryas/methylseq/
cp -r ./conf /home/celphin/scratch/Dryas/methylseq/

cd /home/celphin/scratch/Dryas/methylseq/

nano nextflow.config
# https://nf-co.re/methylseq/2.7.0/parameters/

/*
 * -------------------------------------------------
 *  nf-core/methylseq Nextflow config file
 * -------------------------------------------------
 * Default config options for all environments.
 */

singularity {
    enabled = true
    autoMounts = true
    runOptions = '-B /project:/project -B /scratch:/scratch'
}

process {
    executor = 'slurm'
}

executor {
    queueSize = 500
}

// Global default params, used in configs
params {

    // Input options
    input                      = '/home/celphin/scratch/Dryas/methylseq/input/input_files.csv'
    fasta                      = '/home/celphin/scratch/Dryas/methylseq/input/DoctH0_Main.fasta'
    //add after second run:
    //bismark_index            = '/home/celphin/scratch/Dryas/methylseq/output/BismarkIndex/'


 // Intermediate files
save_reference             = true
save_trimmed               = true

 // Alignment options
comprehensive              = true

  // Library presets
 zymo                      = true

 // Bismark options
 non_directional           = true
 cytosine_report           = true
 relax_mismatches          = true
 no_overlap                = true

 // Boilerplate options
outdir                      = 'output/' //creates and writes to output

}

profiles {
cedar {
        process.clusterOptions = "--account rrg-rieseber-ac"
        max_memory='186G'
        max_cpu=48
        max_time='168h'
    }

}
'

##############################################