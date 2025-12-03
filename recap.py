import numpy as np
import tskit
import msprime
import pyslim
import random
import sys
import collections

def simplify_and_subsample(ts, sample_size):
    a_seed = random.randint(0, 2**31 - 1)
    print(f"Generated random seed for sample: {a_seed}")
    np.random.seed(a_seed)

    num_inds = ts.num_individuals
    print(f"Number of individuals in the population: {num_inds}")
    print(f"Requested sample size: {sample_size}")
    inds = np.random.choice(num_inds, sample_size, replace=False)

    samples = []
    for i in inds:
        samples.extend(ts.individual(i).nodes)
    subsample_nodes = np.sort(np.array(samples))
    ts = ts.simplify(subsample_nodes, keep_unary=True, filter_sites=False)
    return ts

numPops = int(sys.argv[1])
inFile = sys.argv[2]
output_prefix = sys.argv[3]

ts = tskit.load(inFile)

r_seed = random.randint(0, 2**31 - 1)

if numPops == 1:
    ts = pyslim.recapitate(ts, recombination_rate=1e-8, ancestral_Ne=10000, random_seed=r_seed)
    a_seed = random.randint(0, 2**31 - 1)
    np.random.seed(a_seed)
    num_inds = ts.num_individuals
    inds = np.random.choice(num_inds, 50, replace=False)
    samples = []
    for i in inds:
        samples.extend(ts.individual(i).nodes)
    subsample_nodes = np.sort(np.array(samples))
    ts = ts.simplify(subsample_nodes)
    next_id = pyslim.next_slim_mutation_id(ts)
    ts = msprime.sim_mutations(
        ts,
        rate=1.29e-8,
        model=msprime.SLiMMutationModel(type=0, next_id=next_id),
        keep=True
    )
    vcf_ts = pyslim.generate_nucleotides(ts)
    vcf_ts = pyslim.convert_alleles(vcf_ts)
    inds = np.unique([ts.node(i).individual for i in ts.samples()])
    indv_names = [f"tsk_{i}indv" for i in range(len(inds))]
    with open(output_prefix + ".vcf", "w") as vcf_file:
        vcf_ts.write_vcf(vcf_file, individual_names=indv_names, isolated_as_missing=False, position_transform=lambda x: np.fmax(1, x))
elif numPops == 2:
    demography = msprime.Demography()
    # ancestral population
    demography.add_population(name="pop_0", initial_size=10000)
    # split populations
    demography.add_population(name="p1", initial_size=10000)
    demography.add_population(name="p2", initial_size=10000)
    demography.add_population_split(time=1500, derived=["p1", "p2"], ancestral="pop_0")

    ts = msprime.sim_ancestry(
        recombination_rate=1e-8,
        sequence_length=ts.sequence_length,
        initial_state=ts, 
        demography=demography
    )

    # Counter for individuals by population ID
    pop_counts = collections.Counter()
    for ind in ts.individuals():
        pop_counts[ind.population] += 1

    # # Print names and sizes
    # print("Populations and their sizes:")
    for pop_id, count in pop_counts.items():
        pop_metadata = ts.population(pop_id).metadata
        pop_name = pop_metadata.get("name", f"pop_{pop_id}")
        print(f"{pop_name} (ID: {pop_id}): {count} individuals") #pop_name=p1/p2, pop_id=1/2

    # Get individuals for each population
    pop1_id = [pop_id for pop_id in pop_counts if ts.population(pop_id).metadata.get("name") == "p1"]
    pop2_id = [pop_id for pop_id in pop_counts if ts.population(pop_id).metadata.get("name") == "p2"]

    # Get sample nodes for individuals in p1 and p2
    samples_p1 = []
    samples_p2 = []

    for ind in ts.individuals():
        if ind.population in pop1_id:
            samples_p1.extend(ind.nodes)
        elif ind.population in pop2_id:
            samples_p2.extend(ind.nodes)

    # Simplify and write VCF for p1
    ts_p1 = ts.simplify(samples=samples_p1, keep_unary=True,   filter_sites=False)
    ts_p2 = ts.simplify(samples=samples_p2, keep_unary=True,  filter_sites=False)

    next_id = pyslim.next_slim_mutation_id(ts)
    ts = msprime.sim_mutations(
        ts,
        rate=1.29e-8,
        model=msprime.SLiMMutationModel(type=0, next_id=next_id),
        keep=True
    )

    vcf_ts = pyslim.generate_nucleotides(ts)
    vcf_ts = pyslim.convert_alleles(vcf_ts)

    # Use unique individual IDs for naming
    inds = np.unique([ts.node(i).individual for i in ts.samples()])
    indv_names = [f"tsk_{i}indv" for i in range(len(inds))]

    # Collect sample nodes by population
    nodes_by_pop = {1: [], 2: []}  # Assuming p1=0, p2=1
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

        pop_name = f"p{pop_id}"  # SLiM: p1 = pop 1, p2 = pop 2
        vcf_filename = f"{output_prefix}_{pop_name}.vcf"
        print(f"Writing VCF for {pop_name} to {vcf_filename}")

        with open(vcf_filename, "w") as vcf_file:
            sub_ts = vcf_ts.simplify(samples=nodes,  filter_sites=False)
            sub_ts = simplify_and_subsample(sub_ts, 50)

            #ind_names = [f"{pop_name}_ind{i}" for i in range(sub_ts.num_individuals)]


            sub_ts.write_vcf(
                vcf_file,
                #individual_names=ind_names,
                isolated_as_missing=False,
                position_transform=lambda x: np.fmax(1, x)
            )