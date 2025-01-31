#!/bin/bash
#SBATCH --job-name=pixy
#SBATCH --partition=colella
#SBATCH --time=4:00:00
#SBATCH --mem=120G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=output_pixy.%j.txt
#SBATCH --reservation=b686w673_43

module load anaconda
conda activate /panfs/pfs.local/work/bi/l944d476/conda/envs/pixy

bgzip D1000.end.of.div.vcf
tabix D1000.end.of.div.vcf.gz

pixy --stats pi fst dxy \
--vcf D1000.end.of.div.vcf.gz \
--populations D1000.end.of.div.pm.txt \
--bed_file D1000.bed \
--n_cores 12 \
--bypass_invariant_check 'yes' \
--output_prefix D1000.end.of.div
