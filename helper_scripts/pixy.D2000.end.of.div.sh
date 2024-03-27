#!/bin/bash
#SBATCH --job-name=pixy
#SBATCH --partition=sixhour
#SBATCH --time=6:00:00
#SBATCH --mem=120G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=output_pixy.%j.txt

module load anaconda
conda activate /panfs/pfs.local/work/bi/l944d476/conda/envs/pixy

bgzip D2000.end.of.div.invariants.vcf
tabix D2000.end.of.div.invariants.vcf.gz

pixy --stats pi fst dxy \
--vcf D2000.end.of.div.invariants.vcf.gz \
--populations D2000.end.of.div.pm.txt \
--bed_file D2000.bed \
--n_cores 12 \
--output_prefix D2000.end.of.div.invariants
