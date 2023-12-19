def _get_phase_file(wildcards):
    row = manifest[manifest.id == wildcards.id]
    assert len(row) == 1, (wildcards.id, len(row))
    return row.iloc[0].phase_file

rule write_ini_hatchet1:
    input:
        tumor_1bed = os.path.join(sim_results_dir, 'hatchet1', '{id}', 'baf', 'tumor.1bed'),
    params:
        work_dir = os.path.join(sim_results_dir, 'hatchet1', '{id}'),
    output:
        ini_file = os.path.join(sim_results_dir, 'hatchet1', '{id}', 'hatchet.ini'),
    shell:
        """ 
        python scripts/write_ini_hatchet1.py --tumor_1bed {input.tumor_1bed} --work_dir {params.work_dir} --ini_filename {output.ini_file} --min_clones {min_clones} --max_clones {max_clones} --maxcn_diploid {diploid_cmax} --maxcn_tetraploid {tetraploid_cmax}
        """

rule write_ini_hatchet2:
    input:
        tumor_1bed = os.path.join(sim_results_dir, 'hatchet2', '{id}', 'baf', 'tumor.1bed'),
        phase_file = _get_phase_file
    params:
        work_dir = os.path.join(sim_results_dir, 'hatchet2', '{id}'),
    output:
        ini_file = os.path.join(sim_results_dir, 'hatchet2', '{id}', 'hatchet.ini'),
    shell:
        """ 
        python scripts/write_ini_hatchet2.py --tumor_1bed {input.tumor_1bed} --work_dir {params.work_dir} --ini_filename {output.ini_file} --phase_file {input.phase_file} --min_clones {min_clones} --max_clones {max_clones} --maxcn_diploid {diploid_cmax} --maxcn_tetraploid {tetraploid_cmax}
        """

rule run_hatchet:
    input:
        os.path.join(sim_results_dir, '{hatchet_version}', '{id}', 'hatchet.ini'),
    params:
        work_dir = os.path.join(sim_results_dir, '{hatchet_version}', '{id}'),
        stdout = os.path.join(sim_results_dir, '{hatchet_version}', '{id}', 'log.out'),
        stderr = os.path.join(sim_results_dir, '{hatchet_version}', '{id}', 'log.err'),
    output:
        os.path.join(sim_results_dir, '{hatchet_version}', '{id}', 'results', 'best.bbc.ucn'),
        os.path.join(sim_results_dir, '{hatchet_version}', '{id}', 'summary', 'intratumor-clones-allelecn.pdf'),
        os.path.join(sim_results_dir, '{hatchet_version}', '{id}', 'summary', 'intratumor-profiles.pdf'),
        os.path.join(sim_results_dir, '{hatchet_version}', '{id}', 'summary', 'intratumor-clones-totalcn.pdf'),
        os.path.join(sim_results_dir, '{hatchet_version}', '{id}', 'summary', 'intratumor-profilesreduced.pdf'),
        os.path.join(sim_results_dir, '{hatchet_version}', '{id}', 'summary', 'intratumor-copynumber-allelecn.pdf'),
        os.path.join(sim_results_dir, '{hatchet_version}', '{id}', 'summary', 'intratumor-copynumber-totalcn.pdf'),
    resources:
        mem_mb=64000,
        threads=24,
        time_min=1320
    shell:
        """
        cd {params.work_dir}
        module load gurobi
        hatchet run {input} > {params.stdout} 2> {params.stderr}
        """
