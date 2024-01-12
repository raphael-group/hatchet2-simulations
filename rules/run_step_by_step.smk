

rule run_binning_hatchet1:
    input:
        tumor_1bed=os.path.join(sim_results_dir, 'hatchet1', '{id}', 'rdr', 'tumor.1bed'),
        normal_1bed=os.path.join(sim_results_dir, 'hatchet1', '{id}', 'rdr', 'normal.1bed'),
        total_counts=os.path.join(sim_results_dir, 'hatchet1', '{id}', 'rdr', 'total.tsv'),
        baf_file=os.path.join(sim_results_dir, 'hatchet1', '{id}', 'baf', 'tumor.1bed'),
    params:
        work_dir = os.path.join(sim_results_dir, 'hatchet1', '{id}', 'bb'),
    output:
        os.path.join(sim_results_dir, 'hatchet1', '{id}', 'bb', 'bulk.bb'),
    threads: 24
    resources:
        mem_mb=32000
    shell:
        """
        mkdir -p {params.work_dir}
        cd {params.work_dir}
        hatchet combine-counts-fw -c {input.normal_1bed} -C {input.tumor_1bed} -B {input.baf_file} \
            -t {input.total_counts} > {output} 2> log.err
        """

rule run_gmm_clustering:
    input:
        os.path.join(sim_results_dir, 'hatchet1', '{id}', 'bb', 'bulk.bb'),
    params:
        work_dir = os.path.join(sim_results_dir, 'hatchet1', '{id}', 'bbc'),
    output:
        bbc_file=os.path.join(sim_results_dir, 'hatchet1', '{id}', 'bb', 'bulk.bbc'),
        seg_file=os.path.join(sim_results_dir, 'hatchet1', '{id}', 'bb', 'bulk.seg'),
    threads: 16
    resources:
        mem_mb=32000
    shell:
        """
        mkdir -p {params.work_dir}
        cd {params.work_dir}
        hatchet cluster-bins-gmm {input} -o {output.seg_file} -O {output.bbc_file} > log.out 2> log.err
        """

rule run_binning_hatchet2:
    input:
        expand(os.path.join(sim_results_dir, 'hatchet2', '{id}', 'rdr', '{chr}.thresholds.gz'), chr=chromosomes, allow_missing=True),
        expand(os.path.join(sim_results_dir, 'hatchet2', '{id}', 'rdr', '{chr}.total.gz'), chr=chromosomes, allow_missing=True),
        os.path.join(sim_results_dir, 'hatchet2', '{id}', 'rdr', 'samples.txt'),
        total_counts=os.path.join(sim_results_dir, 'hatchet2', '{id}', 'rdr', 'total.tsv'),
        baf_file=os.path.join(sim_results_dir, 'hatchet2', '{id}', 'baf', 'tumor.1bed'),
    params:
        work_dir = os.path.join(sim_results_dir, 'hatchet2', '{id}', 'bb'),
        rdr_dir = os.path.join(sim_results_dir, 'hatchet2', '{id}', 'rdr'),
    output:
        os.path.join(sim_results_dir, 'hatchet2', '{id}', 'bb', 'bulk.bb'),
    threads: 24
    resources:
        mem_mb=32000
    shell:
        """
        mkdir -p {params.work_dir}
        cd {params.work_dir}
        hatchet combine-counts -A {params.rdr_dir} -t {input.total_counts} -b {input.baf_file} -o {output} \
            --msr {min_snpcov_reads} --mtr {min_total_reads} -V hg19 > log.out 2> log.err
        """

rule run_loc_clustering:
    input:
        os.path.join(sim_results_dir, 'hatchet2', '{id}', 'bb', 'bulk.bb'),
    params:
        work_dir = os.path.join(sim_results_dir, 'hatchet2', '{id}', 'bbc'),
        plots_dir = os.path.join(sim_results_dir, 'hatchet2', '{id}', 'plots'),
        plots_1d2d_dir = os.path.join(sim_results_dir, 'hatchet2', '{id}', 'plots', '1d2d'),
    output:
        bbc_file=os.path.join(sim_results_dir, 'hatchet2', '{id}', 'bb', 'bulk.bbc'),
        seg_file=os.path.join(sim_results_dir, 'hatchet2', '{id}', 'bb', 'bulk.seg'),
    threads: 24
    resources:
        mem_mb=32000
    shell:
        """
        mkdir -p {params.work_dir}
        cd {params.work_dir}
        hatchet cluster-bins {input} -o {output.seg_file} -O {output.bbc_file} > log.out 2> log.err
        """

rule plot_clusters:
    input:
        bbc_file=os.path.join(sim_results_dir, '{hatchet_version}', '{id}', 'bb', 'bulk.bbc'),
        seg_file=os.path.join(sim_results_dir, 'hatchet2', '{id}', 'bb', 'bulk.seg'),
    params:
        plots_dir = os.path.join(sim_results_dir, '{hatchet_version}', '{id}', 'plots'),
        plots_1d2d_dir = os.path.join(sim_results_dir, '{hatchet_version}', '{id}', 'plots', '1d2d'),
    output:
        os.path.join(sim_results_dir, '{hatchet_version}', '{id}', 'plots', 'bb.pdf'),
        os.path.join(sim_results_dir, '{hatchet_version}', '{id}', 'plots', 'bb_clustered.png'),
        os.path.join(sim_results_dir, '{hatchet_version}', '{id}', 'plots', 'ballelefrequency.pdf'),
        os.path.join(sim_results_dir, '{hatchet_version}', '{id}', 'plots', 'ballelefrequency_clustered.pdf'),
        os.path.join(sim_results_dir, '{hatchet_version}', '{id}', 'plots', 'readdepthratio.pdf'),
        os.path.join(sim_results_dir, '{hatchet_version}', '{id}', 'plots', 'readdepthratio_clustered.pdf'),
    threads: 1
    resources:
        mem_mb=32000
    shell:
        """
        mkdir -p {params.plots_dir}
        cd {params.plots_dir}

        hatchet plot-bins {input.bbc_file} --rundir {params.plots_dir} --ymin 0 --ymax 3 \
             > log_plot.out 2> log_plot.err

        mkdir -p {params.plots_1d2d_dir}
        hatchet plot-bins-1d2d -s {input.seg_file} -b {input.bbc_file} --outdir {params.plots_1d2d_dir} \
         --centers --centromeres > log_1d2d.out 2> log_1d2d.err
        """

rule run_solver:
    input:
        bbc_file=os.path.join(sim_results_dir, '{hatchet_version}', '{id}', 'bb', 'bulk.bbc'),
        seg_file=os.path.join(sim_results_dir, '{hatchet_version}', '{id}', 'bb', 'bulk.seg'),
    params:
        bbc_dir=os.path.join(sim_results_dir, '{hatchet_version}', '{id}', 'bbc'),
        work_dir = os.path.join(sim_results_dir, '{hatchet_version}', '{id}', 'results'),
    output:
        cn_file=os.path.join(sim_results_dir, '{hatchet_version}', '{id}', 'results', 'best.bbc.ucn'),
    threads: 24
    resources:
        mem_mb=64000
    shell:
        """
        mkdir -p {params.work_dir}
        cd {params.work_dir}
        hatchet compute-cn -x {params.work_dir} -i {params.bbc_dir} --clones {min_clones},{max_clones} \
            --diploidcmax {diploid_cmax} --tetraploidcmax {tetraploid_cmax} > log.out 2> log.err
        """

rule plot_solution:
    input:
        os.path.join(sim_results_dir, '{hatchet_version}', '{id}', 'results', 'best.bbc.ucn'),
    params:
        summary_dir = os.path.join(sim_results_dir, '{hatchet_version}', '{id}', 'summary'),
        summary_1d2d_dir = os.path.join(sim_results_dir, '{hatchet_version}', '{id}', 'summary', '1d2d'),
    output:
        os.path.join(sim_results_dir, '{hatchet_version}', '{id}', 'summary', 'intratumor-mixtures.pdf'),
        os.path.join(sim_results_dir, '{hatchet_version}', '{id}', 'summary', 'intratumor-clones-allelecn.pdf'),
        os.path.join(sim_results_dir, '{hatchet_version}', '{id}', 'summary', 'intratumor-profiles.pdf'),
        os.path.join(sim_results_dir, '{hatchet_version}', '{id}', 'summary', 'intratumor-clones-totalcn.pdf'),
        os.path.join(sim_results_dir, '{hatchet_version}', '{id}', 'summary', 'intratumor-profilesreduced.pdf'),
        os.path.join(sim_results_dir, '{hatchet_version}', '{id}', 'summary', 'intratumor-copynumber-allelecn.pdf'),
        os.path.join(sim_results_dir, '{hatchet_version}', '{id}', 'summary', 'intratumor-copynumber-totalcn.pdf'),
        os.path.join(sim_results_dir, '{hatchet_version}', '{id}', 'summary', 'intratumor-subclonality.pdf'),
    threads: 1
    resources:
        mem_mb=32000
    shell:
        """
        mkdir -p {params.summary_dir}
        cd {params.summary_dir}
        hatchet plot-cn {input} --rundir {params.summary_dir} > log_plot.out 2> log_plot.err

        mkdir -p {params.summary_1d2d_dir}
        hatchet plot-cn-1d2d {input} --outdir {params.summary_1d2d_dir} --bysample --centromeres \
            > log_1d2d.out 2> log_1d2d.err
        """
