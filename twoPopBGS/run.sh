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
    python3 ../recap.py 2 "${tree}" "./neutralRep${i}/twoPopBGS"
    vcftools --vcf "./neutralRep${i}/twoPopBGS_p1.vcf" --min-alleles 2 --max-alleles 2 --recode --out "neutralRep${i}/tmp"
    mv "neutralRep${i}/tmp.recode.vcf" "neutralRep${i}/p1.vcf"
    bgzip "neutralRep${i}/p1.vcf"
    bcftools index "neutralRep${i}/p1.vcf.gz"
    vcftools --vcf "./neutralRep${i}/twoPopBGS_p2.vcf" --min-alleles 2 --max-alleles 2 --recode --out "neutralRep${i}/tmp"
    mv "neutralRep${i}/tmp.recode.vcf" "neutralRep${i}/p2.vcf"
    bgzip "neutralRep${i}/p2.vcf"
    bcftools index "neutralRep${i}/p2.vcf.gz"
    bcftools isec -p "neutralRep${i}/pp" "neutralRep${i}/p1.vcf.gz" "neutralRep${i}/p2.vcf.gz"
    ~/bin/selscan --vcf "neutralRep${i}/pp/0002.vcf" --vcf-ref "neutralRep${i}/pp/0003.vcf" --out "neutralRep${i}/twoPopBGS" --xpehh --xpnsl --pmap  --trunc-ok
    ~/bin/selscan --vcf "neutralRep${i}/pp/0002.vcf" --vcf-ref "neutralRep${i}/pp/0003.vcf" --out "neutralRep${i}/twoPopBGS.unphased" --unphased --xpehh --xpnsl --pmap  --trunc-ok
    rm -rf "neutralRep${i}/pp"
done

# Our non-neutral replicate.
mkdir sweep
slim -d s=0.01 -d out=\"./sweep/\" sim.slim
tree=$(ls sweep/*_final.trees | grep -E '[0-9]+' | sort -t_ -k2,2n | tail -n 1)
python3 ../recap.py 2 "${tree}" "./sweep/twoPopBGS"
vcftools --vcf "./sweep/twoPopBGS_p1.vcf" --min-alleles 2 --max-alleles 2 --recode --out "sweep/tmp"
mv "sweep/tmp.recode.vcf" "sweep/p1.vcf"
bgzip "sweep/p1.vcf"
bcftools index "sweep/p1.vcf.gz"
vcftools --vcf "./sweep/twoPopBGS_p2.vcf" --min-alleles 2 --max-alleles 2 --recode --out "sweep/tmp"
mv "sweep/tmp.recode.vcf" "sweep/p2.vcf"
bgzip "sweep/p2.vcf"
bcftools index "sweep/p2.vcf.gz"
bcftools isec -p "sweep/pp" "sweep/p1.vcf.gz" "sweep/p2.vcf.gz"
~/bin/selscan --vcf "sweep/pp/0002.vcf" --vcf-ref "sweep/pp/0003.vcf" --out "sweep/twoPopBGS" --xpehh --xpnsl --pmap --trunc-ok
~/bin/selscan --vcf "sweep/pp/0002.vcf" --vcf-ref "sweep/pp/0003.vcf" --out "sweep/twoPopBGS.unphased" --unphased --xpehh --xpnsl --pmap --trunc-ok
rm -rf "sweep/pp"

# Normalize w.r.t. neutral sims.
~/bin/norm --xpehh --files neutralRep*/twoPopBGS.xpehh.out sweep/twoPopBGS.xpehh.out --bins 100
~/bin/norm --xpnsl --files neutralRep*/twoPopBGS.xpnsl.out sweep/twoPopBGS.xpnsl.out --bins 100
~/bin/norm --xpehh --files neutralRep*/twoPopBGS.unphased.xpehh.out sweep/twoPopBGS.unphased.xpehh.out --bins 100
~/bin/norm --xpnsl --files neutralRep*/twoPopBGS.unphased.xpnsl.out sweep/twoPopBGS.unphased.xpnsl.out --bins 100


# Clean up
find "." -type f | egrep "log|trees|txt" | xargs rm