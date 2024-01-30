#!/bin/bash
#SBATCH --job-name="sensimp"
#SBATCH --account="dbiopast2"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=100G
#SBATCH --time=14-0:0:00
### Load the conda environment
module load Miniconda3/22.11.1-1
conda init bash
source ~/.bashrc
conda activate rustenv
DIR='/group/pasture/Jeff/imputef/res'
DIR_DATA='/group/pasture/Jeff/imputef/misc'
NREPS=3
PLOIDY=4
STRICTBOUNDARIES=FALSE
cd $DIR
time \
Rscript perf_misc.R \
    ${DIR} \
    ${DIR_DATA}/lucerne_chromosome1_population2.vcf \
    ${NREPS} \
    ${PLOIDY} \
    ${STRICTBOUNDARIES} \
    ${SLURM_CPUS_PER_TASK} \
    ${DIR}/sensitivity_analysis_output.csv
### Execute:
# sbatch perf_misc.slurm
