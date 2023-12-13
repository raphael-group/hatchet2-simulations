
import click
import pandas as pd

@click.command()
@click.option('--work_dir', help='Working directory for HATCHet1 run')
@click.option('--tumor_1bed', help='TSV table that includes "sample" column with all tumor samples')
@click.option('--ini_filename', help='Filename for output hatchet.ini files')
def main(work_dir, tumor_1bed, ini_filename):
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
        f.write('solver=cbc\n')
        f.write('diploidcmax=6\n')
        f.write('tetraploidcmax=14\n')
        f.write('clones=2,4\n')
        
if __name__ == '__main__':
    main()
