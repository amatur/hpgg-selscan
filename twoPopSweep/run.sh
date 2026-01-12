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
    tree="neutralRep${i}/gen4000.trees"
    python3 ../recap.py 2 "${tree}" "./neutralRep${i}/twoPopSweep"
    vcftools --vcf "./neutralRep${i}/twoPopSweep_p1.vcf" --min-alleles 2 --max-alleles 2 --recode --out "neutralRep${i}/tmp"
    mv "neutralRep${i}/tmp.recode.vcf" "neutralRep${i}/p1.vcf"
    bgzip "neutralRep${i}/p1.vcf"
    bcftools index -f "neutralRep${i}/p1.vcf.gz"
    vcftools --vcf "./neutralRep${i}/twoPopSweep_p2.vcf" --min-alleles 2 --max-alleles 2 --recode --out "neutralRep${i}/tmp"
    mv "neutralRep${i}/tmp.recode.vcf" "neutralRep${i}/p2.vcf"
    bgzip "neutralRep${i}/p2.vcf"
    bcftools index -f "neutralRep${i}/p2.vcf.gz"
    bcftools isec -p "neutralRep${i}/pp" "neutralRep${i}/p1.vcf.gz" "neutralRep${i}/p2.vcf.gz"
    selscan --vcf "neutralRep${i}/pp/0002.vcf" --vcf-ref "neutralRep${i}/pp/0003.vcf" --out "neutralRep${i}/twoPopSweep" --xpehh --pmap  --trunc-ok
    selscan --vcf "neutralRep${i}/pp/0002.vcf" --vcf-ref "neutralRep${i}/pp/0003.vcf" --out "neutralRep${i}/twoPopSweep.unphased" --unphased --xpehh --pmap  --trunc-ok
    selscan --vcf "neutralRep${i}/pp/0002.vcf" --vcf-ref "neutralRep${i}/pp/0003.vcf" --out "neutralRep${i}/twoPopSweep" --xpnsl --pmap --trunc-ok
    selscan --vcf "neutralRep${i}/pp/0002.vcf" --vcf-ref "neutralRep${i}/pp/0003.vcf" --out "neutralRep${i}/twoPopSweep.unphased" --unphased --xpnsl --pmap --trunc-ok
    rm -rf "neutralRep${i}/pp"
done

# Our non-neutral replicate.
mkdir sweep
slim -d s=0.1 -d out=\"./sweep/\" sim.slim
tree=$(ls sweep/*_final.trees | grep -E '[0-9]+' | sort -t_ -k2,2n | tail -n 1)
python3 ../recap.py 2 "${tree}" "./sweep/twoPopSweep"
vcftools --vcf "./sweep/twoPopSweep_p1.vcf" --min-alleles 2 --max-alleles 2 --recode --out "sweep/tmp"
mv "sweep/tmp.recode.vcf" "sweep/p1.vcf"
bgzip "sweep/p1.vcf"
bcftools index -f "sweep/p1.vcf.gz"
vcftools --vcf "./sweep/twoPopSweep_p2.vcf" --min-alleles 2 --max-alleles 2 --recode --out "sweep/tmp"
mv "sweep/tmp.recode.vcf" "sweep/p2.vcf"
bgzip "sweep/p2.vcf"
bcftools index -f "sweep/p2.vcf.gz"
bcftools isec -p "sweep/pp" "sweep/p1.vcf.gz" "sweep/p2.vcf.gz"
selscan --vcf "sweep/pp/0002.vcf" --vcf-ref "sweep/pp/0003.vcf" --out "sweep/twoPopSweep" --xpehh --pmap --trunc-ok
selscan --vcf "sweep/pp/0002.vcf" --vcf-ref "sweep/pp/0003.vcf" --out "sweep/twoPopSweep.unphased" --unphased --xpehh --pmap --trunc-ok
selscan --vcf "sweep/pp/0002.vcf" --vcf-ref "sweep/pp/0003.vcf" --out "sweep/twoPopSweep" --xpnsl --pmap --trunc-ok
selscan --vcf "sweep/pp/0002.vcf" --vcf-ref "sweep/pp/0003.vcf" --out "sweep/twoPopSweep.unphased" --unphased --xpnsl --pmap --trunc-ok
rm -rf "sweep/pp"

# Normalize w.r.t. neutral sims.
selscan norm --xpehh --files neutralRep*/twoPopSweep.xpehh.out sweep/twoPopSweep.xpehh.out
selscan norm --xpnsl --files neutralRep*/twoPopSweep.xpnsl.out sweep/twoPopSweep.xpnsl.out
selscan norm --xpehh --files neutralRep*/twoPopSweep.unphased.xpehh.out sweep/twoPopSweep.unphased.xpehh.out
selscan norm --xpnsl --files neutralRep*/twoPopSweep.unphased.xpnsl.out sweep/twoPopSweep.unphased.xpnsl.out

# Clean up
#find "." -type f | egrep "log|trees|txt" | xargs rm
