
import click
import pandas as pd

@click.command()
@click.option('--work_dir', help='Working directory for HATCHet1 run')
@click.option('--tumor_1bed', help='TSV table that includes "sample" column with all tumor samples')
@click.option('--ini_filename', help='Filename for output hatchet.ini files')
@click.option('--min_clones', help="Minimum number of clones for solver", default=2)
@click.option('--max_clones', help="Maximum number of clones for solver", default=4)
@click.option('--maxcn_diploid', help="Minimum copy number for diploid solutions", default=6)
@click.option('--maxcn_tetraploid', help="Maximum number of clones for tetraploid solutions", default=14)
def main(work_dir, tumor_1bed, ini_filename, min_clones, max_clones, maxcn_diploid, maxcn_tetraploid):
    assert min_clones <= max_clones, (min_clones, max_clones)
    assert min_clones > 0, min_clones
    assert maxcn_diploid >= 2, maxcn_diploid
    assert maxcn_tetraploid >= 4, maxcn_tetraploid
    
    if 'tetraploid' in work_dir:
        is_tetraploid = 1
    elif 'diploid' in work_dir:
        is_tetraploid = -1
    else:
        is_tetraploid = 0
    
    df = pd.read_table(tumor_1bed, names = ['chr', 'pos', 'sample', 'ref', 'alt'])    
    
    with open(ini_filename, 'w') as f:
        f.write('[run]\n')
        f.write('download_panel=False\n')
        f.write('count_reads=False\n')
        f.write('genotype_snps=False\n')
        f.write('phase_snps = False\n')
        f.write('fixed_width = True\n')
        f.write('count_alleles=False\n')
        f.write('combine_counts=True\n')
        f.write('cluster_bins=True\n')
        f.write('loc_clust = False\n')
        f.write('plot_bins=True\n')
        f.write('compute_cn=True\n')
        f.write('plot_cn=True\n\n')

        f.write('reference=/n/fs/ragr-data/datasets/ref-genomes/GRCh37_NCBI/GRCh37.p13.chr.fa\n')
        f.write('processes=24\n')
        f.write(f'samples={" ".join(sorted(df["sample"].unique()))}\n')
        f.write(f'output={work_dir}\n\n')
        
        f.write('[cluster_bins_gmm]\n')
        f.write('initclusters=100\n\n')
        
        
        f.write('[compute_cn]\n')
        f.write('solver=cpp\n')
        f.write(f'diploidcmax={maxcn_diploid}\n')
        f.write(f'tetraploidcmax={maxcn_tetraploid}\n')
        f.write(f'clones={min_clones+1},{max_clones+1}\n')
        if is_tetraploid == -1:
            f.write('diploid=True\n')
        elif is_tetraploid == 1:
             f.write('tetraploid=True\n')
             
if __name__ == '__main__':
    main()
