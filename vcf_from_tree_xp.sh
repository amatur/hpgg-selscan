#!/bin/bash
# === Script Overview ===

# Input: SLiM-generated .trees file

# Steps:
# Recapitate tree sequence using msprime
# Add neutral mutations
# Sample haplotypes (e.g., 100)
# Write VCF from mutated tree
# Filter VCF to keep only biallelic SNPs
# Run iHS scan using selscan 

# Outputs (in same directory as input):
# - *.recap.trees         # recapitated tree
# - *.recap.vcf           # mutated VCF
# - *.biallelic.vcf       # filtered VCF
# - *.pconfig*.ihs.out    # iHS scan results


# ======= Check for input =======
if [ $# -lt 1 ]; then
    echo "Usage: $0 <SLiM_tree_file.trees>"
    exit 1
fi

# ======= User-defined constants =======
TREE_FILE=$1   # SLiM-generated .trees file
SCRIPT_RECAP_MUT_SAMPLE="recap_mut_sample_xp.py"
SELSCAN="/storage/home/aur1111/s/transfer/selscan-bugfix/selscan/src/selscan"
SELSCAN="/Users/amatur/code/selscan_bug/src/selscan"
#SELSCAN="/storage/home/aur1111/s/transfer/selscan/src/selscan"
NE=10000
RECOMB_RATE=1e-8
MUT_RATE=1.29e-8
D_SAMP=100

# ======= Derived values =======
#DIR=$(dirname "$TREE_FILE")
#DIR=$(cd "$(dirname "$TREE_FILE")" && pwd)
DIR=$(realpath "$(dirname "$TREE_FILE")")
BASENAME=$(basename "$TREE_FILE" .trees)
DEST_PREFIX="${DIR}/${BASENAME}.recap"
RECAP_TREE="${DEST_PREFIX}.trees"


VCF="${DEST_PREFIX}.vcf"

VCF1="${DEST_PREFIX}_p1.vcf"
VCF2="${DEST_PREFIX}_p2.vcf"

BIALLELIC_VCF1="${DIR}/${BASENAME}_p1.biallelic.vcf"
BIALLELIC_VCF2="${DIR}/${BASENAME}_p2.biallelic.vcf"

LOG="${DEST_PREFIX}.log"

# ======= Step 1: Recapitate + Mutate + Sample =======
mkdir -p "$DIR"
python "$SCRIPT_RECAP_MUT_SAMPLE" \
    --source "$TREE_FILE" \
    --dest_prefix "$DEST_PREFIX" \
    --mu "$MUT_RATE" \
    --ne "$NE" \
    --recomb "$RECOMB_RATE" \
    --sample_size "$D_SAMP" \
    --sample_size_p2 "$D_SAMP" \
    --random \
    --vcf \
    --tree 

# ======= Step 2: Filter for biallelic SNPs =======
TMP_OUT="${DIR}/${BASENAME}_tmp"
vcftools --vcf "$VCF1" \
         --min-alleles 2 --max-alleles 2 \
         --recode --out "$TMP_OUT"
mv "${TMP_OUT}.recode.vcf" "$BIALLELIC_VCF1"

vcftools --vcf "$VCF2" \
         --min-alleles 2 --max-alleles 2 \
         --recode --out "$TMP_OUT"
mv "${TMP_OUT}.recode.vcf" "$BIALLELIC_VCF2"

# ======= Step 3: Run selscan =======
$SELSCAN \
    --vcf "$BIALLELIC_VCF1" \
    --out "${DIR}/${BASENAME}_p1" \
    --ihs --nsl --pmap  \
    --trunc-ok


$SELSCAN \
    --vcf "$BIALLELIC_VCF2" \
    --out "${DIR}/${BASENAME}_p2" \
    --ihs --nsl --pmap  \
    --trunc-ok
