#!/bin/bash
#SBATCH --job-name="VvImp"
#SBATCH --account="dbiopast2"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=100G
#SBATCH --time=0-0:30:00
### Load the conda environment
module load Miniconda3/22.11.1-1
conda init bash
source ~/.bashrc
conda activate rustenv
DIR='/group/pasture/Jeff/imputef'
DIR_DATA='/group/pasture/Jeff/imputef/misc' ### contains grape.vcf (2n), cocksfoot.vcf (2n=4x), and soybean.vcf (2n x 42)
NREPS=3
cd ${DIR}/res
echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
echo ${SLURM_ARRAY_TASK_ID}
time \
Rscript perf.R \
    ${SLURM_ARRAY_TASK_ID} \
    ${DIR} \
    ${DIR_DATA} \
    ${NREPS} \
    ${SLURM_CPUS_PER_TASK}
### Execute:
# sbatch --array=1-60 perf.slurm
