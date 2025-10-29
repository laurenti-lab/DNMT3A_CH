#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH -p icelake
#SBATCH -A LAURENTI-SL2-CPU
#SBATCH --time=0-11:59:59
#SBATCH --output=logs/cellranger_%A_%a_%j.out
#SBATCH --error=logs/cellranger_%A_%a_%j.err


# -------------------------------------------------------
#title           :cellranger_count_array.sh
#author          :Hugo P Bastos
#date            :20220922
#version         :0.9

#description     : sbatch --array=1-6 cellranger_count_array.sbatch
# --------------------------------------------------------


refdir=/rds/project/rds-.../references/10x/


# reference genome to use:
# ------------------------------------------

# cellranger3x/
ref=refdata-cellranger-hg19-3.0.0

# path to fastq dir
# -------------------------------------------------------------------------------
fastqdir=/rds/project/rds-../fastq/

outdir=/rds/project/rds-.../processed/


echo "Running on hosts: $SLURM_NODELIST"
echo "Running on $SLURM_NNODES nodes."
echo "Running on $SLURM_NPROCS processors." 

source /etc/profile.d/modules.sh

module load cellranger/6.0.1

cellwrangler () {

  echo "Started processing library ${library[$1]} on `date`"
  
  cd ${outdir}

  cellranger count --id=${library[$1]} \
                   --sample=${library[$1]} \
                   --transcriptome=${refdir}${ref} \
                   --fastqs=${fastqdir} \
		   --include-introns \
		   --force-cells=${expected[$1]}
                   # \                 
                   # --chemistry=auto
                   # SC3Pv3
                   # --expect-cells=$4\

  echo "Finishing processing library ${library[$1]} on `date`"

}


declare -A library

library[1]=SITTA2
library[2]=SITTB2
library[3]=SITTC2
library[4]=SITTD2
library[5]=SITTE2
library[6]=SITTF2

declare -A expected

expected[1]=12000
expected[2]=12000
expected[3]=5500
expected[4]=12000
expected[5]=9000
expected[6]=10000

cellwrangler $SLURM_ARRAY_TASK_ID

