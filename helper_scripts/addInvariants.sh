#!/bin/bash
#SBATCH --job-name=addInvariants
#SBATCH --partition=sixhour
#SBATCH --time=1:00:00
#SBATCH --mem=120G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=output_addInvariants.%j.txt

module load python
module load anaconda

for d in D750 D1000 D2000
do
  gunzip ../$d/$d.end.of.div.vcf.gz
  ./addInvariants.py ../$d/$d.end.of.div.vcf $d/$d.end.of.div.invariants.vcf
  gzip ../$d/$d.end.of.div.vcf
done


