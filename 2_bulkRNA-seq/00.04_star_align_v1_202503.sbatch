#!/bin/bash
#SBATCH -N 1
#SBATCH -n 6
#SBATCH -p icelake-himem
#SBATCH -A LAURENTI-SL2-CPU
#SBATCH --time=0-6:00:00
#SBATCH --output=logs/star_align_%A_%a_%j.out
#SBATCH --error=logs/star_align_%A_%a_%j.err


# -------------------------------------------------------
#title           :star_align_v1_202503.sbatch
#author          :Giovanna Mantica
#date            :20250325
#version         :v1

#description     : sbatch --array=1-28 star_align_v1_202503.sbatch #parallelising for number of files --array=1-$(wc -l < CHIP30_samples.txt
# --------------------------------------------------------


#specify directories
projdir=/rds/project/rds-0p1wRpSNFIw/users/gm686/CHIP30/
fastqdir=/rds/project/rds-0p1wRpSNFIw/users/gm686/CHIP30/processed/post_trim/fastq #trimmed fastq
bam_outdir=/rds/project/rds-0p1wRpSNFIw/users/gm686/CHIP30/processed/star_bam
refdir=/rds/project/rds-0p1wRpSNFIw/users/gm686/CHIP30/ref

#micromamba activate env
source ~/.bashrc
micromamba activate bulkRNA_aligning_b202503

#job info
echo "Running on hosts: $SLURM_NODELIST"
echo "Running on $SLURM_NNODES nodes."
echo "Running on $SLURM_NPROCS processors." 

echo "Started processing library on `date`"


#get sample name from CHIP30_sample.txt obtained from ls /rds/project/rds-0p1wRpSNFIw/users/gm686/CHIP30/20250110_LH00442_0078_B22VCLFLT3/fastq/*R1_001.fastq.gz | sed 's/_R1_001.fastq.gz//' | awk -F'/' '{print $NF}' > CHIP30_samples.txt

sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" CHIP30_samples.txt)

# ---
# ALIGN TO BAMS
# ---

STAR  --genomeDir "$refdir" --outSAMattributes NH HI NM MD AS --outSAMtype BAM SortedByCoordinate --readFilesIn "$fastqdir/${sample}_R1_001_val_1.fq.gz" "$fastqdir/${sample}_R2_001_val_2.fq.gz" --readFilesCommand gunzip -c --outFileNamePrefix "$bam_outdir/${sample}" --runThreadN 6 --limitBAMsortRAM 31532137230

echo "Finished processing sample ${sample} on $(date)"


