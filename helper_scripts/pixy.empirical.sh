#!/bin/bash
#SBATCH --job-name=pixy
#SBATCH --partition=sixhour
#SBATCH --time=4:00:00
#SBATCH --mem=10G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=output_pixy.%j.txt

#module load anaconda
module load conda
#conda activate /panfs/pfs.local/work/bi/l944d476/conda/envs/pixy
conda activate /kuhpc/work/colella/lab_software/pixy

bgzip snps_3_100_recode.vcf
tabix snps_3_100_recode.vcf.gz

pixy --stats pi fst dxy \
--vcf snps_3_100_recode.vcf.gz \
--populations passerella_samples_3.txt \
--bed_file average.bed \
--n_cores 12 \
--bypass_invariant_check 'yes' \
--output_prefix snps_3_100_recode
