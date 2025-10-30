#!/bin/bash
#SBATCH -N 1
#SBATCH -n 14
#SBATCH -p icelake
#SBATCH -A LAURENTI-SL2-CPU
#SBATCH --time=0-12:00:00
#SBATCH --output=logs/fastqc_%A_%a_%j.out #edit
#SBATCH --error=logs/fastqc_%A_%a_%j.err #edit


# -------------------------------------------------------
#title           :trimgalore_v1_202503.sbatch
#author          :Giovanna Mantica
#date            :20250325
#version         :v1

#description     : sbatch --array=1-28 trimgalore_v1_202503.sbatch #parallelising for number of files --array=1-$(wc -l < CHIP30_samples.txt
# --------------------------------------------------------


#specify directories
projdir=/rds/project/rds-0p1wRpSNFIw/users/gm686/CHIP30/
fastqdir=/rds/project/rds-0p1wRpSNFIw/users/gm686/CHIP30/20250110_LH00442_0078_B22VCLFLT3/fastq
fq_outdir=/rds/project/rds-0p1wRpSNFIw/users/gm686/CHIP30/processed/fastq_posttrim/
#qc_outdir=/rds/project/rds-0p1wRpSNFIw/users/gm686/CHIP30/processed/fastqc_posttrim/ #remember to create directory before running it

#job info
echo "Running on hosts: $SLURM_NODELIST"
echo "Running on $SLURM_NNODES nodes."
echo "Running on $SLURM_NPROCS processors." 
echo "Started processing library on `date`"

#micromamba activate env
source ~/.bashrc
micromamba activate bulkRNA_trimming_b202503

#get sample name from CHIP30_sample.txt obtained from ls /rds/project/rds-0p1wRpSNFIw/users/gm686/CHIP30/20250110_LH00442_0078_B22VCLFLT3/fastq/*R1_001.fastq.gz | sed 's/_R1_001.fastq.gz//' | awk -F'/' '{print $NF}' > CHIP30_samples.txt

sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" CHIP30_samples.txt)

# ---
# STEP1: RUN TRIM GALORE Adapter trimming
# ---

trim_galore -q 20 -fastqc --illumina -stringency 1 -paired -length 20 -e 0.1 -o "$fq_outdir" -trim-n "$fastqdir/${sample}_R1_001.fastq.gz" "$fastqdir/${sample}_R2_001.fastq.gz"

echo "Finished processing sample ${sample} on $(date)"


# ----
# STEP2: RUN FASTQC POST TRIMMING not needed because -fastqc option added above
# ----
#echo "Running FASTQC post adapter trimming "

#fastqc "${fq_outdir}${sample}*.fq*" -o "$qc_outdir"

#echo "FastQC complete. Output saved to $qc_outdir."


