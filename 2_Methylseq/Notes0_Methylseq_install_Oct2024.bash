#########################################################
# Dryas DNA methylation calling
# nf-core Methylseq
# Mapping WGBS data to reference
# Calling methylated sites
#############################################################
# https://nbis-workshop-epigenomics.readthedocs.io/en/latest/content/tutorials/nextflow.html#setup-nf-core

#Installing nextflow/methylseq
# following instructions here: 
# https://nf-co.re/methylseq/2.7.0

# First Setup nextflow profile
# https://docs.alliancecan.ca/wiki/Nextflow

cd /home/celphin/scratch/Dryas/methylseq

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
mkdir /project/def-rieseber/NXF_SINGULARITY_CACHEDIR
export NXF_SINGULARITY_CACHEDIR=/project/def-rieseber/NXF_SINGULARITY_CACHEDIR

# check 
# https://nbis-workshop-epigenomics.readthedocs.io/en/latest/content/tutorials/nextflow.html#setup-nf-core
# mkdir nextflow-hello-test
# cd nextflow-hello-test
# nextflow run hello

#-----------------
# download the methylseq pipeline to scratch
# need allocation
tmux new-session -s Dryas
tmux attach-session -t Dryas

cd /home/celphin/scratch/Dryas/methylseq

#salloc -c6 --time 3:00:00 --mem 120000m --account rrg-rieseber-ac

cd /home/celphin/scratch/Dryas/methylseq
source nf-core-env/bin/activate
module load nextflow/24.04.4
module load apptainer/1.3.4
export NXF_SINGULARITY_CACHEDIR=/project/def-rieseber/NXF_SINGULARITY_CACHEDIR
# set name and version of pipeline
export NFCORE_PL=methylseq
export PL_VERSION=2.7.1

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
  clusterOptions = '--account=def-rieseber'
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
cd /home/celphin/scratch/Dryas/methylseq/nf-core-methylseq_2.7.1/2_7_1

# Note: need to copy the nf-core/methylseq/conf/base.config to the directory running the program
cp nextflow.config /home/celphin/scratch/Dryas/methylseq/nextflow_raw.config
cp -r ./conf /home/celphin/scratch/Dryas/methylseq/

cd /home/celphin/scratch/Dryas/methylseq/

nano /home/celphin/scratch/Dryas/methylseq/nf-core-methylseq_2.7.1/2_7_1/nextflow.config

# https://nf-co.re/methylseq/2.7.1/parameters/

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
    fasta                      = '/home/celphin/scratch/Dryas/methylseq/input/reference/DoctH0_Main.fasta'
    //add after second run:
    //bismark_index            = '/home/celphin/scratch/Dryas/methylseq/input/reference/BismarkIndex/'

    // Intermediate files
    save_reference             = true
    save_align_intermeds       = true
    unmapped                   = true
    save_trimmed               = true

    // Alignment options
    aligner                    = 'bismark'
    comprehensive              = true

    // Library presets
    zymo                       = true

    // Bismark options
    non_directional            = true
    cytosine_report            = true
    relax_mismatches           = true
    num_mismatches             = 0.6
    // 0.6 will allow a penalty of bp * -0.6
    // For 100bp reads, this is -60. Mismatches cost -6, gap opening -5 and gap extension -2
    // So -60 would allow 10 mismatches or ~ 8 x 1-2bp indels
    // Bismark default is 0.2 (L,0,-0.2), Bowtie2 default is 0.6 (L,0,-0.6)
    meth_cutoff                = null
    no_overlap                 = true
    ignore_r1                  = 0
    ignore_r2                  = 2
    ignore_3prime_r1           = 0
    ignore_3prime_r2           = 2
    known_splices              = null
    local_alignment            = false
    minins                     = null
    maxins                     = null
    nomeseq                    = false

    // Boilerplate options
    outdir                       = 'output/'

}

profiles {
cedar {
        process.clusterOptions = "--account rrg-rieseber-ac"
        max_memory='186G'
        max_cpu=48
        max_time='168h'
    }
narval {
        process.clusterOptions = "--account def-rieseber"
        max_memory='249G'
        max_cpu=64
        max_time='168h'
    }

}

#------------------------
cd /home/celphin/scratch/Dryas/methylseq/nf-core-methylseq_2.7.1/2_7_1/conf
 nano base.config
 
# change times to shorter to help run faster

    withLabel:process_single {
        cpus   = { 1                   }
        memory = { 6.GB * task.attempt }
        time   = { 3.h  * task.attempt }
    }
    withLabel:process_low {
        cpus   = { 2     * task.attempt }
        memory = { 12.GB * task.attempt }
        time   = { 3.h   * task.attempt }
    }
    withLabel:process_medium {
        cpus   = { 6     * task.attempt }
        memory = { 36.GB * task.attempt }
        time   = { 5.h   * task.attempt }
    }
    withLabel:process_high {
        cpus   = { 12    * task.attempt }
        memory = { 72.GB * task.attempt }
        time   = { 7.h  * task.attempt }
    }
    withLabel:process_long {
        time   = { 20.h  * task.attempt }
    }
    withLabel:process_high_memory {
        memory = { 200.GB * task.attempt }
    }
    withLabel:error_ignore {
        errorStrategy = 'ignore'
    }
    withLabel:error_retry {
        errorStrategy = 'retry'
        maxRetries    = 2
    }
    withName: BISMARK_ALIGN {
        time = { 4.d * task.attempt }
    }
    withName: BISMARK_DEDUPLICATE {
        time = { 2.d * task.attempt }
    }
    withName: BISMARK_METHYLATIONEXTRACTOR {
        time = { 1.d * task.attempt }
    }
    withName: BWAMETH_ALIGN {
        time = { 3.d * task.attempt }
    }

#---------------------------
# add this to bashrc to prevent java from needing too much memory

nano ~/.bashrc
NXF_OPTS='-Xms1g -Xmx4g'


##############################################