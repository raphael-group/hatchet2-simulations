
import click
import pandas as pd
import os
from shutil import copy

@click.command()
@click.option('--work_dir', help='Working directory for HATCHet1 run')
@click.option('--tumor_1bed', help='TSV table that includes "sample" column with all tumor samples')
@click.option('--ini_filename', help='Filename for output hatchet.ini files')
@click.option('--phase_file', help='File containing phasing results for these SNPs')
def main(work_dir, tumor_1bed, ini_filename, phase_file):
    df = pd.read_table(tumor_1bed, names = ['chr', 'pos', 'sample', 'ref', 'alt'])    
    
    with open(ini_filename, 'w') as f:
        f.write('[run]\n')
        f.write('download_panel=False\n')
        f.write('count_reads=False\n')
        f.write('genotype_snps=False\n')
        f.write('phase_snps = False\n')
        f.write('fixed_width = False\n')
        f.write('count_alleles=False\n')
        f.write('combine_counts=True\n')
        f.write('cluster_bins=True\n')
        f.write('loc_clust = True\n')
        f.write('plot_bins=True\n')
        f.write('compute_cn=True\n')
        f.write('plot_cn=True\n\n')

        f.write('reference=/n/fs/ragr-data/datasets/ref-genomes/GRCh37_NCBI/GRCh37.p13.chr.fa\n')
        f.write('reference_version=hg19\n')
        f.write('processes=24\n')
        f.write(f'samples={" ".join(sorted(df["sample"].unique()))}\n')
        f.write(f'output={work_dir}\n\n')
        
        f.write('[combine_counts]\n')
        f.write('msr=1500\n')
        f.write('mtr=1500\n\n')

        f.write('[genotype_snps]\n')
        f.write('reference_version=hg19\n')
        f.write('chr_notation=True\n')
        f.write('snps=/n/fs/ragr-data/datasets/dbSNP/hg37.vcf.gz\n\n')

        f.write('[cluster_bins]\n')
        f.write('minK=35\n')
        f.write('maxK=100\n')
        f.write('diploidbaf=0.07\n')
        f.write('tau=1e-4\n\n')

        f.write('[download_panel]\n')
        f.write('refpaneldir=/n/fs/ragr-data/datasets/phasing-panels/\n\n')

        f.write('[compute_cn]\n')
        f.write('solver=cbc\n')
        f.write('diploidcmax=6\n')
        f.write('tetraploidcmax=14\n')
        f.write('clones=2,4\n')
        
    if not os.path.exists(os.path.join(work_dir, 'phase')):
        os.makedirs(os.path.join(work_dir, 'phase'))
    if not os.path.exists(os.path.join(work_dir, 'phase', 'phased.vcf.gz')):
        copy(phase_file, os.path.join(work_dir, 'phase', 'phased.vcf.gz'))
        
if __name__ == '__main__':
    main()
