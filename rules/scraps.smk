
rule create_solve_ini:
    input:
        os.path.join(sim_results_dir, '{hatchet_version}', '{id}', 'hatchet.ini'),
    params:
        work_dir=os.path.join(sim_results_dir, '{hatchet_version}', '{id}', 'results'),
    output:
        os.path.join(sim_results_dir, '{hatchet_version}', '{id}', 'results', 'hatchet.ini'),
    run:    
        work_dir = params[0]
        allowed_true = ['fixed_width', 'loc_clust', 'compute_cn', 'plot_cn']
        if not os.path.exists(workdir):
            os.makedirs(work_dir)
        with open(input[0], 'r') as infile:
            with open(output[0], 'w') as outfile:
                for line in infile:
                    if line.endswith('True'):
                        first_token = line.split('=')[0].strip()
                        if not first_token in allowed_true:
                            outfile.write(line.replace('True', 'False'))
                        else:
                            outfile.write(line)
                    else:
                        outfile.write(line)

rule run_solver:
    input:
        bbc=os.path.join(sim_results_dir, '{hatchet_version}', '{id}', 'bbc', 'bulk.bbc'),
        seg=os.path.join(sim_results_dir, '{hatchet_version}', '{id}', 'bbc', 'bulk.seg'),
        ini=os.path.join(sim_results_dir, '{hatchet_version}', '{id}', 'results', 'hatchet.ini'),
    params:
        work_dir=os.path.join(sim_results_dir, '{hatchet_version}', '{id}', 'results')
    resources:
        mem_mb=32000,
        threads=24,
        time_min=1320
    shell:
        """
        mkdir -p {params.work_dir}
        cd {params.work_dir}
        
        module load gurobi
        hatchet run {input} > {params.stdout} 2> {params.stderr}
        """