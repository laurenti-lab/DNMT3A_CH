#!/bin/bash
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -p icelake
#SBATCH -A LAURENTI-SL2-CPU
#SBATCH --time=0-5:00:00
#SBATCH --output=logs/fastqc_%A_%a_%j.out #edit
#SBATCH --error=logs/fastqc_%A_%a_%j.err #edit


# -------------------------------------------------------
#title           :star_createref_v1_202503.sbatch
#author          :Giovanna Mantica
#date            :20250325
#version         :v1

#description     : sbatch star_createref_v1_202503.sbatch
# --------------------------------------------------------


#specify directories
projdir=/rds/project/rds-0p1wRpSNFIw/users/gm686/CHIP30/
refdir=/rds/project/rds-0p1wRpSNFIw/users/gm686/CHIP30/ref

#micromamba activate env
source ~/.bashrc
micromamba activate bulkRNA_aligning_b202503

#job info
echo "Running on hosts: $SLURM_NODELIST"
echo "Running on $SLURM_NNODES nodes."
echo "Running on $SLURM_NPROCS processors." 


# ---
# CREATE INDEXED GENOME FROM GENCODE FILES
# ---
echo "Started processing library on `date`"

STAR --runMode genomeGenerate --genomeDir "$refdir" --genomeFastaFiles "$refdir/GRCh37.p13.genome.fa" --sjdbGTFfile "$refdir/gencode.v19.annotation.gtf" --sjdbOverhang 99 --runThreadN 16

echo "Finished creating ref genome on $(date)"




