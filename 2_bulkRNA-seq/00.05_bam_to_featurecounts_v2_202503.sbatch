#!/bin/bash
#SBATCH -N 1
#SBATCH -n 24
#SBATCH -p icelake-himem
#SBATCH -A LAURENTI-SL2-CPU
#SBATCH --time=0-1:00:00
#SBATCH --output=logs/bam_to_featurecounts_%A_%a_%j.out
#SBATCH --error=logs/bam_to_featurecounts_%A_%a_%j.err


# -------------------------------------------------------
#title           :bam_to_featurecounts_v2_202503.sbatch
#author          :Giovanna Mantica
#date            :20250325
#version         :v2 exon only

#description	 : copy script in the folder where the bams are stored and run it from there.
# --------------------------------------------------------


#specify directories
projdir=/rds/project/rds-0p1wRpSNFIw/users/gm686/CHIP30/
bamdir=/rds/project/rds-0p1wRpSNFIw/users/gm686/CHIP30/processed/star_bam
outdir=/rds/project/rds-0p1wRpSNFIw/users/gm686/CHIP30/output/tables
refdir=/rds/project/rds-0p1wRpSNFIw/users/gm686/CHIP30/ref

#micromamba activate env
source ~/.bashrc
micromamba activate bulkRNA_featurecount_b202503

#job info
echo "Running on hosts: $SLURM_NODELIST"
echo "Running on $SLURM_NNODES nodes."
echo "Running on $SLURM_NPROCS processors." 

echo "Started processing library on `date`"


#get sample name from CHIP30_sample.txt obtained from ls /rds/project/rds-0p1wRpSNFIw/users/gm686/CHIP30/20250110_LH00442_0078_B22VCLFLT3/fastq/*R1_001.fastq.gz | sed 's/_R1_001.fastq.gz//' | awk -F'/' '{print $NF}' > CHIP30_samples.txt

#sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" CHIP30_samples.txt)

# ---
# ALIGN TO BAMS
# ---

featureCounts -a "${refdir}/gencode.v19.annotation.gtf" -o "${outdir}/202503_CHIP30_countmatrix.txt" -t exon -T 24 -Q 30 -p --countReadPairs -B -C -g gene_name *.bam

echo "Finished processing on $(date)"


