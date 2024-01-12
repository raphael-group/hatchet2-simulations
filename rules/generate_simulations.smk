
def _get_genome_file(wildcards):
    row = manifest[manifest.id == wildcards.id]
    assert len(row) == 1, (wildcards.id, len(row))
    return row.iloc[0].genome_file

def _get_snps_file(wildcards):
    row = manifest[manifest.id == wildcards.id]
    assert len(row) == 1, (wildcards.id, len(row))
    return row.iloc[0].snps_file

def _get_mixture_file(wildcards):
    row = manifest[manifest.id == wildcards.id]
    assert len(row) == 1, (wildcards.id, len(row))
    return row.iloc[0].mixture_file

def _get_events_string(wildcards):
    row = manifest[manifest.id == wildcards.id]
    assert len(row) == 1, (wildcards.id, len(row))
    return row.iloc[0].events_string

def _get_seed(wildcards):
    row = manifest[manifest.id == wildcards.id]
    assert len(row) == 1, (wildcards.id, len(row))
    return row.iloc[0].seed

def _get_prop_mirrored(wildcards):
    row = manifest[manifest.id == wildcards.id]
    assert len(row) == 1, (wildcards.id, len(row))
    return row.iloc[0].prop_mirrored

def _get_prop_focal_mirrored(wildcards):
    row = manifest[manifest.id == wildcards.id]
    assert len(row) == 1, (wildcards.id, len(row))
    return row.iloc[0].prop_focal_mirrored

def _get_prop_simplify(wildcards):
    row = manifest[manifest.id == wildcards.id]
    assert len(row) == 1, (wildcards.id, len(row))
    row = row.iloc[0]
    if 'prop_simplify' in row:
        return row.prop_simplify
    else:
        return 0

rule simulate_genome:
    input:
        genome_file=_get_genome_file,
    params:
        events_string=_get_events_string,
        seed=_get_seed,
        prop_mirrored=_get_prop_mirrored,
        prop_focal_mirrored=_get_prop_focal_mirrored,
        prop_simplify=_get_prop_simplify,
    output:
        os.path.join(sim_data_dir, 'simulated_genome{id}.tsv'),
    shell: 
        """
        python3 scripts/simulate_genome.py {input.genome_file} --n_clones {n_clones} --event_sizes {params.events_string} --prop_mirrored {params.prop_mirrored} --prop_focal_mirrored {params.prop_focal_mirrored} --prop_simplify {params.prop_simplify} --genome_version {genome_version} --seed {params.seed} --output {output}
        """

rule sample_chromosome:
    input:
        snps_file=_get_snps_file,
        mixture_file=_get_mixture_file,
        genome=os.path.join(sim_data_dir, 'simulated_genome{id}.tsv'),        
    params:
        outdir=os.path.join(intermediate_dir, '{id}'),
        seed=_get_seed,
    output:
        h1_rdr_tumor=os.path.join(intermediate_dir, '{id}', 'hatchet1_rd', '{chromosome}.tumor.1bed'),
        h1_rdr_normal=os.path.join(intermediate_dir, '{id}', 'hatchet1_rd', '{chromosome}.normal.1bed'),
        h1_baf_tumor=os.path.join(intermediate_dir, '{id}', 'baf', '{chromosome}.tumor_baf.1bed'),
        h1_baf_normal=os.path.join(intermediate_dir, '{id}', 'baf', '{chromosome}.normal_baf.1bed'),
        h2_arr=os.path.join(intermediate_dir, '{id}', 'hatchet2_rd', '{chromosome}.total.gz'),
        h2_thresholds=os.path.join(intermediate_dir, '{id}', 'hatchet2_rd', '{chromosome}.thresholds.gz'),
        totalcounts=os.path.join(intermediate_dir, '{id}', 'totalcounts', '{chromosome}.total.tsv'),
        raw_counts=os.path.join(intermediate_dir, '{id}', 'raw_counts', '{chromosome}.raw_counts.tsv'),
    shell:
        """
        mkdir -p {params.outdir}
        python scripts/sample_chromosome.py --genome {input.genome} --snps {input.snps_file} --chromosome {wildcards.chromosome} --sample_mixtures {input.mixture_file} --tumor_coverage {tumor_coverage} --normal_coverage {normal_coverage} --seed {params.seed} --genome_version {genome_version} --bin_width {bin_width} --normal_name {normal_name} --outdir {params.outdir}
        """

rule combine_intermediate_files:
    input:
        expand(os.path.join(intermediate_dir, '{id}', 'hatchet1_rd', '{chromosome}.tumor.1bed'), chromosome=chromosomes, allow_missing=True),
        expand(os.path.join(intermediate_dir, '{id}', 'hatchet1_rd', '{chromosome}.normal.1bed'), chromosome=chromosomes, allow_missing=True),
        expand(os.path.join(intermediate_dir, '{id}', 'baf', '{chromosome}.tumor_baf.1bed'), chromosome=chromosomes, allow_missing=True),
        expand(os.path.join(intermediate_dir, '{id}', 'baf', '{chromosome}.normal_baf.1bed'), chromosome=chromosomes, allow_missing=True),
        expand(os.path.join(intermediate_dir, '{id}', 'hatchet2_rd', '{chromosome}.total.gz'), chromosome=chromosomes, allow_missing=True),
        expand(os.path.join(intermediate_dir, '{id}', 'hatchet2_rd', '{chromosome}.thresholds.gz'), chromosome=chromosomes, allow_missing=True),
        expand(os.path.join(intermediate_dir, '{id}', 'totalcounts', '{chromosome}.total.tsv'), chromosome=chromosomes, allow_missing=True),
        expand(os.path.join(intermediate_dir, '{id}', 'raw_counts', '{chromosome}.raw_counts.tsv'), chromosome=chromosomes, allow_missing=True), 
    params:
        intermediate_dir=os.path.join(intermediate_dir, '{id}'),
        hatchet1_dir=os.path.join(sim_results_dir, 'hatchet1', '{id}'),
        hatchet2_dir=os.path.join(sim_results_dir, 'hatchet2', '{id}'),
    output:
        h1_baf_tumor=os.path.join(sim_results_dir, 'hatchet1', '{id}', 'baf', 'tumor.1bed'),
        h1_baf_normal=os.path.join(sim_results_dir, 'hatchet1', '{id}', 'baf', 'normal.1bed'),
        h1_totalcounts=os.path.join(sim_results_dir, 'hatchet1', '{id}', 'rdr', 'total.tsv'),
        h2_baf_tumor=os.path.join(sim_results_dir, 'hatchet2', '{id}', 'baf', 'tumor.1bed'),
        h2_baf_normal=os.path.join(sim_results_dir, 'hatchet2', '{id}', 'baf', 'normal.1bed'),
        h2_totalcounts=os.path.join(sim_results_dir, 'hatchet2', '{id}', 'rdr', 'total.tsv'),
        
        h1_rdr_tumor=os.path.join(sim_results_dir, 'hatchet1', '{id}', 'rdr', 'tumor.1bed'),
        h1_rdr_normal=os.path.join(sim_results_dir, 'hatchet1', '{id}', 'rdr', 'normal.1bed'),
        
        h2_samples=os.path.join(sim_results_dir, 'hatchet2', '{id}', 'rdr', 'samples.txt'),
        h2_totals=expand(os.path.join(sim_results_dir, 'hatchet2', '{id}', 'rdr', '{chromosome}.total.gz'), chromosome=chromosomes, allow_missing=True),
        h2_threshodls=expand(os.path.join(sim_results_dir, 'hatchet2', '{id}', 'rdr', '{chromosome}.thresholds.gz'), chromosome=chromosomes, allow_missing=True),
    shell:
        """
        python scripts/combine_chromosomes.py --intermediate_dir {params.intermediate_dir} --hatchet1_dir {params.hatchet1_dir} --hatchet2_dir {params.hatchet2_dir}
        """