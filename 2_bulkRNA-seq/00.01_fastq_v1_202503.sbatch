#!/bin/bash
#SBATCH -N 1
#SBATCH -n 6 
#SBATCH -p icelake
#SBATCH -A LAURENTI-SL2-CPU
#SBATCH --time=0-03:00:00
#SBATCH --output=logs/fastqc_%A_%a_%j.out
#SBATCH --error=logs/fastqc_%A_%a_%j.err


# -------------------------------------------------------
#title           :fastq_v1_202503.sbatch
#author          :Giovanna Mantica
#date            :20250325
#version         :v1

#description     : sbatch fastq_v1_202503.sbatch #for fastqc of the file
# --------------------------------------------------------


projdir=/rds/project/rds-0p1wRpSNFIw/users/gm686/CHIP30/



# path to dir where fastq files are
# -------------------------------------------------------------------------------
fastqdir=/rds/project/rds-0p1wRpSNFIw/users/gm686/CHIP30/20250110_LH00442_0078_B22VCLFLT3/fastq



outdir=/rds/project/rds-0p1wRpSNFIw/users/gm686/CHIP30/processed/fastqc_before_processing #remember to create directory before running it



echo "Running on hosts: $SLURM_NODELIST"
echo "Running on $SLURM_NNODES nodes."
echo "Running on $SLURM_NPROCS processors." 
echo "Started processing library ${library[$1]} on `date`"

source ~/.bashrc
micromamba activate bulkRNA_trimming_b202503


echo "Running FastQC on $projdir"

fastqc "$fastqdir"/* -o "$outdir" #remember always forward slash when applicable

echo "FastQC complete. Output saved to $outdir."


