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
        python scripts/write_ini_hatchet1.py --tumor_1bed {input.tumor_1bed} --work_dir {params.work_dir} --ini_filename {output.ini_file}
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
        python scripts/write_ini_hatchet2.py --tumor_1bed {input.tumor_1bed} --work_dir {params.work_dir} --ini_filename {output.ini_file} --phase_file {input.phase_file}
        """

rule run_hatchet:
    input:
        os.path.join(sim_results_dir, '{hatchet_version}', '{id}', 'hatchet.ini'),
    params:
        work_dir = os.path.join(sim_results_dir, '{hatchet_version}', '{id}'),
    output:
        result=os.path.join(sim_results_dir, '{hatchet_version}', '{id}', 'results', 'best.bbc.ucn'),
    resources:
        mem_mb=32000,
        threads=24,
        time_min=1320
    shell:
        """
        cd {params.work_dir}
        hatchet run {input} > {output.stdout} 2> {output.stderr}
        """
