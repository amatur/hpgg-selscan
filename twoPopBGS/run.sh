#!/bin/bash

set -uex

# Generate our neutral reps.
for i in {1..100}
do 
    mkdir -p "neutralRep${i}"
    slim -d s=0 -d out=\"./neutralRep${i}/\" sim.slim
    tree="neutralRep${i}/gen1000.trees"
    python3 ../recap.py 2 "${tree}" "./neutralRep${i}/twoPopBGS"
    vcftools --vcf "./neutralRep${i}/twoPopBGS_p1.vcf" --min-alleles 2 --max-alleles 2 --recode --out "neutralRep${i}/tmp"
    mv "neutralRep${i}/tmp.recode.vcf" "neutralRep${i}/p1.vcf"
    vcftools --vcf "./neutralRep${i}/twoPopBGS_p2.vcf" --min-alleles 2 --max-alleles 2 --recode --out "neutralRep${i}/tmp"
    mv "neutralRep${i}/tmp.recode.vcf" "neutralRep${i}/p2.vcf"
    selscan --vcf "neutralRep${i}/p1.vcf" --vcf-ref "neutralRep${i}/p2.vcf" --out "neutralRep${i}/twoPopBGS" --xpehh --xpnsl --pmap  --trunc-ok
done

# Our non-neutral replicate.
mkdir sweep
slim -d s=0.01 -d out=\"./sweep/\" sim.slim
tree=$(ls sweep/*_final.trees | grep -E '[0-9]+' | sort -t_ -k2,2n | tail -n 1)
python3 ../recap.py 1 "${tree}" "./sweep/twoPopBGS"
vcftools --vcf "./sweep/twoPopBGS_p1.vcf" --min-alleles 2 --max-alleles 2 --recode --out "sweep/tmp"
mv "sweep/tmp.recode.vcf" "sweep/p1.vcf"
vcftools --vcf "./sweep/twoPopBGS_p2.vcf" --min-alleles 2 --max-alleles 2 --recode --out "sweep/tmp"
mv "sweep/tmp.recode.vcf" "sweep/p2.vcf"
selscan --vcf "sweep/p1.vcf" --vcf-ref "sweep/p2.vcf" --out "sweep/twoPopBGS" --ihs --nsl --pmap --trunc-ok

# Normalize w.r.t. neutral sims.
selscan norm --xpehh --files neutralRep*/twoPopBGS.xpehh.out sweep/twoPopBGS.xpehh.out --bins 100
selscan norm --xpnsl --files neutralRep*/twoPopBGS.xpnsl.out sweep/twoPopBGS.xpnsl.out --bins 100