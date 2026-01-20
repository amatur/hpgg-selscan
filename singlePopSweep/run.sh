#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=100GB
#SBATCH --time=72:00:00
#SBATCH --account=zps5164_sc_default
#SBATCH --mail-user=tqs5778@psu.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END

set -uex

# Generate our neutral models.
mkdir -p neutralRep{1..100}
parallel -j 20 slim -d s=0 -d 'out=\"./neutralRep{}/singlePopSweep\"' sim.slim ::: {1..100}

# Run selscan.
for i in {1..100}
do 
    vcftools --vcf "./neutralRep${i}/singlePopSweep.vcf" --min-alleles 2 --max-alleles 2 --recode --out "neutralRep${i}/tmp"
    mv "neutralRep${i}/tmp.recode.vcf" "neutralRep${i}/p1.vcf"
    selscan --vcf "neutralRep${i}/p1.vcf" --out "neutralRep${i}/singlePopSweep" --ihs --nsl --pmap  --trunc-ok
    selscan --vcf "neutralRep${i}/p1.vcf" --out "neutralRep${i}/singlePopSweep.unphased" --unphased --ihs --nsl --pmap  --trunc-ok
    selscan --vcf "neutralRep${i}/p1.vcf" --out "neutralRep${i}/singlePopSweep" --ihh12 --pmap
done

# Our non-neutral replicate.
mkdir sweep
slim -d s=0.1 -d 'out="./sweep/"' sim.slim
vcftools --vcf "./sweep/singlePopSweep.vcf" --min-alleles 2 --max-alleles 2 --recode --out "sweep/tmp"
mv "sweep/tmp.recode.vcf" "sweep/p1.vcf"
selscan --vcf "sweep/p1.vcf" --out "sweep/singlePopSweep" --ihs --nsl --pmap --trunc-ok
selscan --vcf "sweep/p1.vcf" --out "sweep/singlePopSweep.unphased" --unphased --ihs --nsl --pmap --trunc-ok
selscan --vcf "sweep/p1.vcf" --out "sweep/singlePopSweep" --ihh12 --pmap --trunc-ok

# Normalize w.r.t. neutral sims.
selscan norm --ihs --files neutralRep*/singlePopSweep.ihs.out sweep/singlePopSweep.ihs.out --bins 100
selscan norm --nsl --files neutralRep*/singlePopSweep.nsl.out sweep/singlePopSweep.nsl.out --bins 100
selscan norm --ihs --files neutralRep*/singlePopSweep.unphased.ihs.out sweep/singlePopSweep.unphased.ihs.out --bins 100
selscan norm --nsl --files neutralRep*/singlePopSweep.unphased.nsl.out sweep/singlePopSweep.unphased.nsl.out --bins 100
selscan norm --ihh12 --files neutralRep*/singlePopSweep.ihh12.out sweep/singlePopSweep.ihh12.out

# Clean up
#find "." -type f | egrep "log|trees|txt" | xargs rm