
import click
import pandas as pd
import os
from shutil import copy

@click.command()
@click.option('--work_dir', help='Working directory for HATCHet1 run')
@click.option('--tumor_1bed', help='TSV table that includes "sample" column with all tumor samples')
@click.option('--ini_filename', help='Filename for output hatchet.ini files')
@click.option('--phase_file', help='File containing phasing results for these SNPs')
@click.option('--min_clones', help="Minimum number of clones for solver", default=2)
@click.option('--msr', help="Minimum number of SNP-covering reads per bin", default=1000)
@click.option('--mtr', help="Minimum number of total reads per bin", default=1000)
@click.option('--max_clones', help="Maximum number of clones for solver", default=4)
@click.option('--maxcn_diploid', help="Minimum copy number for diploid solutions", default=6)
@click.option('--maxcn_tetraploid', help="Maximum number of clones for tetraploid solutions", default=14)
def main(work_dir, tumor_1bed, ini_filename, phase_file, min_clones, max_clones, maxcn_diploid, maxcn_tetraploid,
         msr, mtr):
    assert min_clones <= max_clones, (min_clones, max_clones)
    assert min_clones > 0, min_clones
    assert maxcn_diploid >= 2, maxcn_diploid
    assert maxcn_tetraploid >= 4, maxcn_tetraploid
    df = pd.read_table(tumor_1bed, names = ['chr', 'pos', 'sample', 'ref', 'alt'])    
    
    if 'tetraploid' in work_dir:
        is_tetraploid = 1
    elif 'diploid' in work_dir:
        is_tetraploid = -1
    else:
        is_tetraploid = 0
    
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
        f.write(f'msr={msr}\n')
        f.write(f'mtr={mtr}\n\n')

        f.write('[genotype_snps]\n')
        f.write('reference_version=hg19\n')
        f.write('chr_notation=True\n')
        f.write('snps=/n/fs/ragr-data/datasets/dbSNP/hg37.vcf.gz\n\n')

        f.write('[cluster_bins]\n')
        f.write('minK=30\n')
        f.write('maxK=45\n')
        f.write('diploidbaf=0.025\n')
        f.write('tau=1e-12\n')

        f.write('[download_panel]\n')
        f.write('refpaneldir=/n/fs/ragr-data/datasets/phasing-panels/\n\n')

        f.write('[compute_cn]\n')
        f.write('solver=cpp\n')
        f.write(f'diploidcmax={maxcn_diploid}\n')
        f.write(f'tetraploidcmax={maxcn_tetraploid}\n')
        f.write(f'clones={min_clones+1},{max_clones+1}\n')
        if is_tetraploid == -1:
            f.write('diploid=True\n')
        elif is_tetraploid == 1:
             f.write('tetraploid=True\n')
           
        
    if not os.path.exists(os.path.join(work_dir, 'phase')):
        os.makedirs(os.path.join(work_dir, 'phase'))
    if not os.path.exists(os.path.join(work_dir, 'phase', 'phased.vcf.gz')):
        copy(phase_file, os.path.join(work_dir, 'phase', 'phased.vcf.gz'))
        
if __name__ == '__main__':
    main()
