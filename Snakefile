import pandas as pd
import os

chromosomes= config['chromosomes']
manifest = pd.read_table(config['manifest'])
manifest.id = manifest.id.astype(str)

sim_data_dir = config['sim_data_dir']
sim_results_dir = config['sim_results_dir']
intermediate_dir = os.path.join(sim_data_dir, 'intermediate')
analysis_dir = config['analysis_dir']
n_clones = config['n_clones']

# HATCHet2 binning parameters (used when running HATCHet2)
min_snpcov_reads = 500
min_total_reads = 500

# solving parameters
min_clones = 2 # actual number of clones, so nclones argument gets this value +1 for normal
max_clones = 2
diploid_cmax = 5
tetraploid_cmax = 10


# global parameters
genome_version = 'hg19'
normal_name = 'normal'
bin_width = 100000 # this is used when generating HATCHet1 input
normal_coverage = 30
tumor_coverage = 120

simulation_ids = manifest['id'].unique()
hatchet_versions = ['hatchet1', 'hatchet2']
#hatchet_versions = ['hatchet2']

rule all:
    input:
        # input files
        # expand(os.path.join(sim_results_dir, 'hatchet1', '{id}', 'baf', 'tumor.1bed'), id=simulation_ids),
        # expand(os.path.join(sim_results_dir, 'hatchet1', '{id}', 'baf', 'normal.1bed'), id=simulation_ids),
        # expand(os.path.join(sim_results_dir, 'hatchet1', '{id}', 'rdr', 'total.tsv'), id=simulation_ids),
        # expand(os.path.join(sim_results_dir, 'hatchet1', '{id}', 'rdr', 'tumor.1bed'), id=simulation_ids),
        # expand(os.path.join(sim_results_dir, 'hatchet1', '{id}', 'rdr', 'normal.1bed'), id=simulation_ids),
        expand(os.path.join(sim_results_dir, 'hatchet2', '{id}', 'baf', 'tumor.1bed'), id=simulation_ids),
        expand(os.path.join(sim_results_dir, 'hatchet2', '{id}', 'baf', 'normal.1bed'), id=simulation_ids),
        expand(os.path.join(sim_results_dir, 'hatchet2', '{id}', 'rdr', 'total.tsv'), id=simulation_ids),
        
        expand(os.path.join(sim_results_dir, 'hatchet2', '{id}', 'rdr', 'samples.txt'), id=simulation_ids),
        expand(os.path.join(sim_results_dir, 'hatchet2', '{id}', 'rdr', '{chromosome}.total.gz'), chromosome=chromosomes, id=simulation_ids),
        expand(os.path.join(sim_results_dir, 'hatchet2', '{id}', 'rdr', '{chromosome}.thresholds.gz'), chromosome=chromosomes, id=simulation_ids),

        # final results
        expand(os.path.join(sim_results_dir, '{hatchet_version}', '{simulation_id}', 'results', 'best.bbc.ucn'),
            hatchet_version=hatchet_versions, simulation_id=simulation_ids),

        # final plots
        #expand(os.path.join(sim_results_dir, '{hatchet_version}', '{id}', 'summary', 'intratumor-mixtures.pdf'), hatchet_version = hatchet_versions, id=simulation_ids),
        #expand(os.path.join(sim_results_dir, '{hatchet_version}', '{id}', 'summary', 'intratumor-subclonality.pdf'), hatchet_version = hatchet_versions, id=simulation_ids),
        expand(os.path.join(sim_results_dir, '{hatchet_version}', '{id}', 'summary', 'intratumor-clones-allelecn.pdf'), hatchet_version = hatchet_versions, id=simulation_ids),
        expand(os.path.join(sim_results_dir, '{hatchet_version}', '{id}', 'summary', 'intratumor-clones-totalcn.pdf'), hatchet_version = hatchet_versions, id=simulation_ids),
        expand(os.path.join(sim_results_dir, '{hatchet_version}', '{id}', 'summary', 'intratumor-profiles.pdf'), hatchet_version = hatchet_versions, id=simulation_ids),
        expand(os.path.join(sim_results_dir, '{hatchet_version}', '{id}', 'summary', 'intratumor-profilesreduced.pdf'), hatchet_version = hatchet_versions, id=simulation_ids),
        expand(os.path.join(sim_results_dir, '{hatchet_version}', '{id}', 'summary', 'intratumor-copynumber-allelecn.pdf'), hatchet_version = hatchet_versions, id=simulation_ids),
        expand(os.path.join(sim_results_dir, '{hatchet_version}', '{id}', 'summary', 'intratumor-copynumber-totalcn.pdf'), hatchet_version = hatchet_versions, id=simulation_ids),

        # analysis files
        expand(os.path.join(analysis_dir, 'joint_tables', '{hatchet_version}__{id}.tsv'), hatchet_version = hatchet_versions, id=simulation_ids),
        expand(os.path.join(analysis_dir, 'metrics', '{hatchet_version}__{id}.json'), hatchet_version = hatchet_versions, id=simulation_ids),


include: "rules/generate_simulations.smk"
include: "rules/run_methods.smk"