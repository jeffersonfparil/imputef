#!/bin/bash
#SBATCH --job-name="imputef"
#SBATCH --account="dbiopast2"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=500G
#SBATCH --time=14-0:0:00
### Load the conda environment
module load Miniconda3/22.11.1-1
conda init bash
source ~/.bashrc
conda activate rustenv
DIR='/group/pasture/Jeff/imputef/res'
NREPS=3
cd $DIR
echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
echo ${SLURM_ARRAY_TASK_ID}
DIR_DATA='/group/pasture/Jeff/imputef/misc'
echo "2,${DIR_DATA}/grape.vcf" > ploidy_vcf-${SLURM_ARRAY_TASK_ID}.tmp      ###   2.9 Mb
echo "2,${DIR_DATA}/apple.vcf" >> ploidy_vcf-${SLURM_ARRAY_TASK_ID}.tmp     ###  24.0 Mb
echo "84,${DIR_DATA}/soybean.vcf" >> ploidy_vcf-${SLURM_ARRAY_TASK_ID}.tmp  ###  36.0 Mb
echo "2,${DIR_DATA}/zucchini.vcf" >> ploidy_vcf-${SLURM_ARRAY_TASK_ID}.tmp  ### 135.0 Mb
echo "4,${DIR_DATA}/lucerne.vcf" >> ploidy_vcf-${SLURM_ARRAY_TASK_ID}.tmp   ### 400.0 Mb
time \
for i in $(seq 1 $(cat ploidy_vcf-${SLURM_ARRAY_TASK_ID}.tmp | wc -l))
do
    # i=4; NREPS=1; SLURM_ARRAY_TASK_ID=40; SLURM_CPUS_PER_TASK=32
    PLOIDY=$(head -n${i} ploidy_vcf-${SLURM_ARRAY_TASK_ID}.tmp | tail -n1 | cut -d',' -f1)
    VCF=$(head -n${i} ploidy_vcf-${SLURM_ARRAY_TASK_ID}.tmp | tail -n1 | cut -d',' -f2)
    echo "*********************************************************************************************"
    echo ${VCF}-${PLOIDY}X
    time \
    Rscript perf.R \
        ${VCF} \
        ${PLOIDY} \
        ${SLURM_ARRAY_TASK_ID} \
        ${NREPS} \
        ${SLURM_CPUS_PER_TASK}
done
rm ploidy_vcf-${SLURM_ARRAY_TASK_ID}.tmp

### Execute:
# sbatch --array=1-24%12 perf.slurm
