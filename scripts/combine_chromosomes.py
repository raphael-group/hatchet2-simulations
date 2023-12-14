import click
import os
from shutil import copy
import pandas as pd
import numpy as np

@click.command()
@click.option('--intermediate_dir', help='Path to intermediate simulated files')
@click.option('--hatchet1_dir', help='Path to intermediate simulated files')
@click.option('--hatchet2_dir', help='Path to intermediate simulated files')
def main(intermediate_dir, hatchet1_dir, hatchet2_dir):
    totalcount_files = os.listdir(os.path.join(intermediate_dir, 'totalcounts'))
    baf_files = os.listdir(os.path.join(intermediate_dir, 'baf'))
    rdr1_files = os.listdir(os.path.join(intermediate_dir, 'hatchet1_rd'))
    rdr2_files = os.listdir(os.path.join(intermediate_dir, 'hatchet2_rd'))

    chromosomes1 = set([a.split('.')[0] for a in totalcount_files])
    chromosomes2 = set([a.split('.')[0] for a in baf_files])
    chromosomes3 = set([a.split('.')[0] for a in rdr1_files])
    chromosomes4 = set([a.split('.')[0] for a in rdr2_files if not a == 'samples.txt'])
    assert chromosomes1 == chromosomes2 == chromosomes3 == chromosomes4

    ## make output directories
    if not os.path.exists(os.path.join(hatchet1_dir, 'baf')):
        os.makedirs(os.path.join(hatchet1_dir, 'baf'))
    if not os.path.exists(os.path.join(hatchet2_dir, 'baf')):
        os.makedirs(os.path.join(hatchet2_dir, 'baf'))

    if not os.path.exists(os.path.join(hatchet1_dir, 'rdr')):
        os.makedirs(os.path.join(hatchet1_dir, 'rdr'))
    if not os.path.exists(os.path.join(hatchet2_dir, 'rdr')):
        os.makedirs(os.path.join(hatchet2_dir, 'rdr'))


    # copy HATCHet2 RDR files and check array sizes
    sample_names = [a.strip() for a in open(os.path.join(os.path.join(intermediate_dir, 'hatchet2_rd', 'samples.txt')))]

    for ch in chromosomes4:
        my_thresholds = np.loadtxt(os.path.join(intermediate_dir, 'hatchet2_rd', f'{ch}.thresholds.gz'))
        my_array = np.loadtxt(os.path.join(intermediate_dir, 'hatchet2_rd', f'{ch}.total.gz'))
        assert my_array.shape == (len(my_thresholds), 2 * len(sample_names))
        
        copy(os.path.join(intermediate_dir, 'hatchet2_rd', f'{ch}.thresholds.gz'),
                os.path.join(hatchet2_dir, 'rdr', f'{ch}.thresholds.gz'))
        copy(os.path.join(intermediate_dir, 'hatchet2_rd', f'{ch}.total.gz'),
                os.path.join(hatchet2_dir, 'rdr', f'{ch}.total.gz'))
    copy(os.path.join(intermediate_dir, 'hatchet2_rd', 'samples.txt'),
            os.path.join(hatchet2_dir, 'rdr', 'samples.txt'))

   # merge files that are the same for both methods
    ## total read counts per sample
    tc = [pd.read_table(os.path.join(intermediate_dir, 'totalcounts', f), header = None, index_col = 0) for f in totalcount_files]
    totalcounts = tc[0].copy()
    if len(tc) > 0:
        for x in tc[1:]:
            assert totalcounts.index.equals(x.index), (totalcounts.index, x.index)
            totalcounts.iloc[:, 0] += x.iloc[:, 0]
    totalcounts.to_csv(os.path.join(hatchet1_dir, 'rdr', 'total.tsv'), header = False, sep = '\t')
    totalcounts.to_csv(os.path.join(hatchet2_dir, 'rdr', 'total.tsv'), header = False, sep = '\t')

    joint_tumor_baf = pd.concat([pd.read_table(os.path.join(intermediate_dir, 'baf', f)) for f in baf_files if 'normal' not in f])[['chr', 'pos', 'sample', 'ref', 'alt']].sort_values(by=['sample', 'chr', 'pos'])
    joint_normal_baf = pd.concat([pd.read_table(os.path.join(intermediate_dir, 'baf', f)) for f in baf_files if 'normal' in f])[['chr', 'pos', 'sample', 'ref', 'alt']].sort_values(by=['sample', 'chr', 'pos'])

    joint_tumor_baf.to_csv(os.path.join(hatchet1_dir, 'baf', 'tumor.1bed'), index = False, header = False, sep = '\t')
    joint_normal_baf.to_csv(os.path.join(hatchet1_dir, 'baf', 'normal.1bed'), index = False, header = False, sep = '\t')
    joint_tumor_baf.to_csv(os.path.join(hatchet2_dir, 'baf', 'tumor.1bed'), index = False, header = False, sep = '\t')
    joint_normal_baf.to_csv(os.path.join(hatchet2_dir, 'baf', 'normal.1bed'), index = False, header = False, sep = '\t')

    # merge HATCHet1 RDR files 
    joint_tumor_rdr = pd.concat([pd.read_table(os.path.join(intermediate_dir, 'hatchet1_rd', f), header = None) for f in rdr1_files if 'normal' not in f])
    joint_normal_rdr = pd.concat([pd.read_table(os.path.join(intermediate_dir, 'hatchet1_rd', f), header = None) for f in rdr1_files if 'normal' in f])

    joint_tumor_rdr.to_csv(os.path.join(hatchet1_dir, 'rdr', 'tumor.1bed'), index = False, header = False, sep = '\t')
    joint_normal_rdr.to_csv(os.path.join(hatchet1_dir, 'rdr', 'normal.1bed'), index = False, header = False, sep = '\t')

if __name__ == '__main__':
    main()
