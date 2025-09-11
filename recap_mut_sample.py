#!/usr/bin/env python3

import argparse
import collections
import numpy as np
import tskit
import msprime
import pyslim
import random

parser = argparse.ArgumentParser()
parser.add_argument("--source", type=str, help="Path to the input .trees file (including .trees extension)")
parser.add_argument("--dest_prefix", type=str, required=True, help="Prefix for output files")
parser.add_argument("--mu", type=float, required=True, help="Mutation rate per base per generation")
parser.add_argument("--random", action="store_true", help="Sample a random subset of diploid individuals")
parser.add_argument("--sample_size", type=int, help="Number of samples to draw (required if using --random)")
parser.add_argument("--seed", type=int, help="Random seed for reproducibility (required if using --random)")
parser.add_argument("--vcf", action="store_true", help="Write VCF output")
parser.add_argument("--tree", action="store_true", help="Write .trees output")
parser.add_argument("--recomb", type=float, default=1e-8, help="Recombination rate per base per generation")
parser.add_argument("--ne", type=float, default=1e4, help="Ancestral effective population size")

parser.add_argument("--norecap", action="store_true",
                    help="If set, do not recapitate the tree sequence before mutation overlay")

#popname = p1


args = parser.parse_args()

source_file = args.source
ts = tskit.load(source_file)


r_seed = random.randint(0, 2**31 - 1)
print(f"Generated random seed for recap: {r_seed}")


# demography = msprime.Demography()
# # ancestral population
# demography.add_population(name="pop_0", initial_size=10000)
# # split populations
# demography.add_population(name="p1", initial_size=10000)
# demography.add_population(name="p2", initial_size=10000)
# demography.add_population_split(time=500, derived=["p1", "p2"], ancestral="pop_0")


# Recapitate the tree sequence
if not args.norecap:
    ts = pyslim.recapitate(ts, recombination_rate=args.recomb, ancestral_Ne=args.ne, random_seed=args.seed if args.seed is not None else r_seed)
    print("Recapitation completed.")
    # ts = msprime.sim_ancestry(
    # recombination_rate=args.recomb,
    # sequence_length=ts.sequence_length,
    # initial_state=ts, demography=demography)
else:
    print("Skipping recapitation as --norecap is set.")


if args.random:
    # if args.sample_size is None or args.seed is None:
    #     raise ValueError("--sample_size and --seed are required when using --random")

    if args.sample_size is None:
        raise ValueError("--sample_size required when using --random")


    a_seed = random.randint(0, 2**31 - 1)
    print(f"Generated random seed for sample: {a_seed}")
    np.random.seed(a_seed)

    num_inds = ts.num_individuals
    inds = np.random.choice(num_inds, args.sample_size // 2, replace=False)

    samples = []
    for i in inds:
        samples.extend(ts.individual(i).nodes)

    subsample_nodes = np.sort(np.array(samples))
    ts = ts.simplify(subsample_nodes)

# Mutation overlay section
# Prepare to overlay new mutations on the tree sequence using SLiM-compatible format.
# - Gets next available SLiM mutation ID.
# - Applies mutations under the SLiMMutationModel with the provided mutation rate.
# - Keeps existing mutations.
# - Uses a fixed random seed for reproducibility.

m_seed = random.randint(0, 2**31 - 1)
print(f"Generated random seed for mut: {m_seed}")


next_id = pyslim.next_slim_mutation_id(ts)
ts = msprime.sim_mutations(
    ts,
    rate=args.mu,
    model=msprime.SLiMMutationModel(type=0, next_id=next_id),
    keep=True,
    random_seed=m_seed,
)

output_prefix = args.dest_prefix

# Write VCF output if requested
if args.vcf:
    print("Writing VCF file:", output_prefix + ".vcf")
    vcf_ts = pyslim.generate_nucleotides(ts)
    vcf_ts = pyslim.convert_alleles(vcf_ts)

    # Use unique individual IDs for naming
    inds = np.unique([ts.node(i).individual for i in ts.samples()])
    indv_names = [f"tsk_{i}indv" for i in range(len(inds))]

    # Write VCF with all samples included
    # 'isolated_as_missing=False' ensures all sample nodes are included even if they are isolated in the graph
    with open(output_prefix + ".vcf", "w") as vcf_file:
        vcf_ts.write_vcf(vcf_file, individual_names=indv_names, isolated_as_missing=False, position_transform=lambda x: np.fmax(1, x))
    print("VCF file written successfully.")

# Write .trees output if requested
if args.tree:
    print("Writing .trees file:", output_prefix + ".trees")
    ts.dump(output_prefix + ".trees")
    print(".trees file written successfully.")




if args.vcf:
    print("Writing VCF files for each subpopulation")

    vcf_ts = pyslim.generate_nucleotides(ts)
    vcf_ts = pyslim.convert_alleles(vcf_ts)

    # Collect sample nodes by population
    nodes_by_pop = {0: [], 1: []}  # Assuming p1=0, p2=1
    for node in vcf_ts.samples():
        ind_id = vcf_ts.node(node).individual
        if ind_id != tskit.NULL:
            pop_id = vcf_ts.individual(ind_id).population
            if pop_id in nodes_by_pop:
                nodes_by_pop[pop_id].append(node)

    for pop_id, nodes in nodes_by_pop.items():
        if not nodes:
            print(f"No samples found for population {pop_id}, skipping.")
            continue

        pop_name = f"p{pop_id+1}"  # SLiM: p1 = pop 0, p2 = pop 1
        vcf_filename = f"{output_prefix}_{pop_name}.vcf"
        print(f"Writing VCF for {pop_name} to {vcf_filename}")

        with open(vcf_filename, "w") as vcf_file:
            sub_ts = vcf_ts.simplify(samples=nodes)
            ind_names = [f"{pop_name}_ind{i}" for i in range(sub_ts.num_individuals)]
            sub_ts.write_vcf(
                vcf_file,
                individual_names=ind_names,
                isolated_as_missing=False,
                position_transform=lambda x: np.fmax(1, x),
            )

    print("VCF files written for all subpopulations.")
