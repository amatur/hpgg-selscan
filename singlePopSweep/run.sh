#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=8GB
#SBATCH --time=5:00:00
#SBATCH --account=zps5164_sc_default
#SBATCH --mail-user=tqs5778@psu.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END

set -uex

# Generate our neutral reps.
for i in {1..100}
do 
    mkdir -p "neutralRep${i}"
    slim -d s=0 -d out=\"./neutralRep${i}/\" sim.slim
    tree="neutralRep${i}/gen1000.trees"
    python3 ../recap.py 1 "${tree}" "./neutralRep${i}/singlePopSweep"
    vcftools --vcf "./neutralRep${i}/singlePopSweep.vcf" --min-alleles 2 --max-alleles 2 --recode --out "neutralRep${i}/tmp"
    mv "neutralRep${i}/tmp.recode.vcf" "neutralRep${i}/p1.vcf"
    selscan --vcf "neutralRep${i}/p1.vcf" --out "neutralRep${i}/singlePopSweep" --ihs --nsl --pmap  --trunc-ok
    selscan --vcf "neutralRep${i}/p1.vcf" --out "neutralRep${i}/singlePopSweep.unphased" --unphased --ihs --nsl --pmap  --trunc-ok
done

# Our non-neutral replicate.
mkdir sweep
slim -d s=0.01 -d out=\"./sweep/\" sim.slim
tree=$(ls sweep/*_final.trees | grep -E '[0-9]+' | sort -t_ -k2,2n | tail -n 1)
python3 ../recap.py 1 "${tree}" "./sweep/singlePopSweep"
vcftools --vcf "./sweep/singlePopSweep.vcf" --min-alleles 2 --max-alleles 2 --recode --out "sweep/tmp"
mv "sweep/tmp.recode.vcf" "sweep/p1.vcf"
selscan --vcf "sweep/p1.vcf" --out "sweep/singlePopSweep" --ihs --nsl --pmap --trunc-ok
selscan --vcf "sweep/p1.vcf" --out "sweep/singlePopSweep.unphased" --unphased --ihs --nsl --pmap --trunc-ok

# Normalize w.r.t. neutral sims.
norm --ihs --files neutralRep*/singlePopSweep.ihs.out sweep/singlePopSweep.ihs.out --bins 100
norm --nsl --files neutralRep*/singlePopSweep.nsl.out sweep/singlePopSweep.nsl.out --bins 100
norm --ihs --files neutralRep*/singlePopSweep.unphased.ihs.out sweep/singlePopSweep.unphased.ihs.out --bins 100
norm --nsl --files neutralRep*/singlePopSweep.unphased.nsl.out sweep/singlePopSweep.unphased.nsl.out --bins 100

# Clean up
find "." -type f | egrep "log|trees|txt" | rm