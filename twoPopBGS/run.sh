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

# Generate neutral reps.
mkdir -p neutralRep{1..100}
parallel -j 20 slim -d s=0 -d 'out="./neutralRep{}/twoPopBGS"' sim.slim ::: {1..100}

# Run selscan.
for i in {1..100}
do 
    vcftools --vcf "./neutralRep${i}/twoPopBGS_p1.vcf" --min-alleles 2 --max-alleles 2 --recode --out "neutralRep${i}/tmp"
    mv "neutralRep${i}/tmp.recode.vcf" "neutralRep${i}/p1.vcf"
    bgzip "neutralRep${i}/p1.vcf"
    bcftools index -f "neutralRep${i}/p1.vcf.gz"
    vcftools --vcf "./neutralRep${i}/twoPopBGS_p2.vcf" --min-alleles 2 --max-alleles 2 --recode --out "neutralRep${i}/tmp"
    mv "neutralRep${i}/tmp.recode.vcf" "neutralRep${i}/p2.vcf"
    bgzip "neutralRep${i}/p2.vcf"
    bcftools index -f "neutralRep${i}/p2.vcf.gz"
    bcftools isec -p "neutralRep${i}/pp" "neutralRep${i}/p1.vcf.gz" "neutralRep${i}/p2.vcf.gz"
    selscan --vcf "neutralRep${i}/pp/0002.vcf" --vcf-ref "neutralRep${i}/pp/0003.vcf" --out "neutralRep${i}/twoPopBGS" --xpehh --pmap  --trunc-ok
    selscan --vcf "neutralRep${i}/pp/0002.vcf" --vcf-ref "neutralRep${i}/pp/0003.vcf" --out "neutralRep${i}/twoPopBGS.unphased" --unphased --xpehh --pmap  --trunc-ok
    selscan --vcf "neutralRep${i}/pp/0002.vcf" --vcf-ref "neutralRep${i}/pp/0003.vcf" --out "neutralRep${i}/twoPopBGS" --xpnsl --pmap  --trunc-ok
    selscan --vcf "neutralRep${i}/pp/0002.vcf" --vcf-ref "neutralRep${i}/pp/0003.vcf" --out "neutralRep${i}/twoPopBGS.unphased" --unphased --xpnsl --pmap  --trunc-ok
    #rm -rf "neutralRep${i}/pp"
done

# Our non-neutral replicate.
mkdir sweep
slim -d s=0.1 -d 'out="./sweep/twoPopBGS"' sim.slim
vcftools --vcf "./sweep/twoPopBGS_p1.vcf" --min-alleles 2 --max-alleles 2 --recode --out "sweep/tmp"
mv "sweep/tmp.recode.vcf" "sweep/p1.vcf"
bgzip "sweep/p1.vcf"
bcftools index -f "sweep/p1.vcf.gz"
vcftools --vcf "./sweep/twoPopBGS_p2.vcf" --min-alleles 2 --max-alleles 2 --recode --out "sweep/tmp"
mv "sweep/tmp.recode.vcf" "sweep/p2.vcf"
bgzip "sweep/p2.vcf"
bcftools index -f "sweep/p2.vcf.gz"
bcftools isec -p "sweep/pp" "sweep/p1.vcf.gz" "sweep/p2.vcf.gz"
selscan --vcf "sweep/pp/0002.vcf" --vcf-ref "sweep/pp/0003.vcf" --out "sweep/twoPopBGS" --xpehh --pmap --trunc-ok
selscan --vcf "sweep/pp/0002.vcf" --vcf-ref "sweep/pp/0003.vcf" --out "sweep/twoPopBGS" --xpnsl --pmap --trunc-ok
selscan --vcf "sweep/pp/0002.vcf" --vcf-ref "sweep/pp/0003.vcf" --out "sweep/twoPopBGS.unphased" --unphased --xpehh --pmap --trunc-ok
selscan --vcf "sweep/pp/0002.vcf" --vcf-ref "sweep/pp/0003.vcf" --out "sweep/twoPopBGS.unphased" --unphased --xpnsl --trunc-ok
#rm -rf "sweep/pp"

# Normalize w.r.t. neutral sims.
selscan norm --xpehh --files neutralRep*/twoPopBGS.xpehh.out sweep/twoPopBGS.xpehh.out 
selscan norm --xpnsl --files neutralRep*/twoPopBGS.xpnsl.out sweep/twoPopBGS.xpnsl.out 
selscan norm --xpehh --files neutralRep*/twoPopBGS.unphased.xpehh.out sweep/twoPopBGS.unphased.xpehh.out 
selscan norm --xpnsl --files neutralRep*/twoPopBGS.unphased.xpnsl.out sweep/twoPopBGS.unphased.xpnsl.out 


# Clean up
#find "." -type f | egrep "log|trees|txt" | xargs rm

