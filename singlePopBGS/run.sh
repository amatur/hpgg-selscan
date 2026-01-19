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

# Create enutral reps.
mkdir -p neutralReps{1..100}
for i in {1..100}
do
    echo "slim -d s=0 -d out=\"./neutralRep${i}/singlePopBGS\" sim.slim"
done | parallel -j 20

# Run selscan.
for i in {1..100}
do 
    vcftools --vcf "./neutralRep${i}/singlePopBGS.vcf" --min-alleles 2 --max-alleles 2 --recode --out "neutralRep${i}/tmp"
    mv "neutralRep${i}/tmp.recode.vcf" "neutralRep${i}/p1.vcf"
    selscan --vcf "neutralRep${i}/p1.vcf" --out "neutralRep${i}/singlePopBGS" --ihs --nsl --pmap  --trunc-ok
    selscan --vcf "neutralRep${i}/p1.vcf" --out "neutralRep${i}/singlePopBGS.unphased" --unphased --ihs --nsl --pmap  --trunc-ok
    selscan --vcf "neutralRep${i}/p1.vcf" --out "neutralRep${i}/singlePopBGS" --ihh12 --pmap  --trunc-ok
done

# Our non-neutral replicate.
mkdir sweep
slim -d s=0.1 -d out=\"./sweep/singlePopBGS\" sim.slim
vcftools --vcf "./sweep/singlePopBGS.vcf" --min-alleles 2 --max-alleles 2 --recode --out "sweep/tmp"
mv "sweep/tmp.recode.vcf" "sweep/p1.vcf"
selscan --vcf "sweep/p1.vcf" --out "sweep/singlePopBGS" --ihs --nsl --pmap --trunc-ok
selscan --vcf "sweep/p1.vcf" --out "sweep/singlePopBGS.unphased" --unphased --ihs --nsl --pmap --trunc-ok
selscan --vcf "sweep/p1.vcf" --out "sweep/singlePopBGS" --ihh12 --pmap --trunc-ok

# Normalize w.r.t. neutral sims.
selscan norm --ihs --files neutralRep*/singlePopBGS.ihs.out sweep/singlePopBGS.ihs.out --bins 100
selscan norm --nsl --files neutralRep*/singlePopBGS.nsl.out sweep/singlePopBGS.nsl.out --bins 100
selscan norm --ihs --files neutralRep*/singlePopBGS.unphased.ihs.out sweep/singlePopBGS.unphased.ihs.out --bins 100
selscan norm --nsl --files neutralRep*/singlePopBGS.unphased.nsl.out sweep/singlePopBGS.unphased.nsl.out --bins 100
selscan norm --ihh12 --files neutralRep*/singlePopBGS.ihh12.out sweep/singlePopBGS.ihh12.out

# Clean up
#find "." -type f | egrep "log|trees|txt" | xargs rm

selscan --ehh 503623 --vcf neutralRep1/p1.vcf --pmap --out neutral
selscan --ehh 500000 --vcf sweep/p1.vcf --pmap --out sweep