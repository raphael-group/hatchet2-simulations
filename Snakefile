import pandas as pd
import os

chromosomes= ['chr20', 'chr21']
manifest = pd.read_table('manifest.tsv')
manifest.id = manifest.id.astype(str)

sim_data_dir = '/n/fs/ragr-research/projects/hatchet2-results/newsims/data'
sim_results_dir = '/n/fs/ragr-research/projects/hatchet2-results/newsims/results/'
intermediate_dir = os.path.join(sim_data_dir, 'intermediate')

# global parameters
genome_version = 'hg19'
normal_name = 'normal'
bin_width = 100000
normal_coverage = 30
tumor_coverage = 80

simulation_ids = manifest['id'].unique()

rule all:
    input:
        expand(os.path.join(sim_results_dir, 'hatchet1', '{id}', 'baf', 'tumor.1bed'), id=simulation_ids),
        expand(os.path.join(sim_results_dir, 'hatchet1', '{id}', 'baf', 'normal.1bed'), id=simulation_ids),
        expand(os.path.join(sim_results_dir, 'hatchet1', '{id}', 'rdr', 'total.tsv'), id=simulation_ids),
        expand(os.path.join(sim_results_dir, 'hatchet2', '{id}', 'baf', 'tumor.1bed'), id=simulation_ids),
        expand(os.path.join(sim_results_dir, 'hatchet2', '{id}', 'baf', 'normal.1bed'), id=simulation_ids),
        expand(os.path.join(sim_results_dir, 'hatchet2', '{id}', 'rdr', 'total.tsv'), id=simulation_ids),
        
        expand(os.path.join(sim_results_dir, 'hatchet1', '{id}', 'rdr', 'tumor.1bed'), id=simulation_ids),
        expand(os.path.join(sim_results_dir, 'hatchet1', '{id}', 'rdr', 'normal.1bed'), id=simulation_ids),
        
        expand(os.path.join(sim_results_dir, 'hatchet2', '{id}', 'rdr', 'samples.txt'), id=simulation_ids),
        expand(os.path.join(sim_results_dir, 'hatchet2', '{id}', 'rdr', '{chromosome}.total.gz'), chromosome=chromosomes, id=simulation_ids),
        expand(os.path.join(sim_results_dir, 'hatchet2', '{id}', 'rdr', '{chromosome}.thresholds.gz'), chromosome=chromosomes, id=simulation_ids),

include: "rules/generate_simulations.smk"
include: "rules/run_methods.smk"