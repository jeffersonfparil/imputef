#!/bin/bash
#SBATCH --job-name="imputef"
#SBATCH --account="dbiopast2"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=250G
#SBATCH --time=10-0:0:00
### Load the conda environment
module load Miniconda3/22.11.1-1
conda init bash
source ~/.bashrc
conda activate rustenv
DIR='/group/pasture/Jeff/imputef/res'
DIR_DATA='/group/pasture/Jeff/imputef/misc' ### contains grape.vcf (2n), lucerne.vcf (2n=4x), and soybean.vcf (2n x 42)
NREPS=3
cd $DIR
echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
echo ${SLURM_ARRAY_TASK_ID}
time \
for i in $(seq 1 $(cat ploidy_vcf-${SLURM_ARRAY_TASK_ID}.tmp | wc -l))
do
    # i=4; NREPS=1; SLURM_ARRAY_TASK_ID=40; SLURM_CPUS_PER_TASK=32
    time \
    Rscript perf.R \
        ${DIR_DATA} \
        ${SLURM_ARRAY_TASK_ID} \
        ${NREPS} \
        ${SLURM_CPUS_PER_TASK}
done
rm ploidy_vcf-${SLURM_ARRAY_TASK_ID}.tmp

### Execute:
# sbatch --array=1-44 perf.slurm