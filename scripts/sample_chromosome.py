import click
import numpy as np
import pandas as pd
from scipy.stats import poisson, binom
import os 

def load_centromeres(centromeres_file):
    use_chr = True
    centromeres = pd.read_table(centromeres_file, header = None, 
                                    names = ['CHR', 'START', 'END', 'NAME', 'gieStain'])
    chr2centro = {}
    for ch in centromeres.CHR.unique():
        my_df = centromeres[centromeres.CHR == ch]
        assert (my_df.gieStain == 'acen').all()
        # Each centromere should consist of 2 adjacent segments
        assert len(my_df == 2)
        assert my_df.START.max() == my_df.END.min()
        if use_chr:
            if ch.startswith('chr'):
                chr2centro[ch] = my_df.START.min(), my_df.END.max()
            else:
                chr2centro['chr' + ch] = my_df.START.min(), my_df.END.max()
        else:
            if ch.startswith('chr'):
                chr2centro[ch[3:]] = my_df.START.min(), my_df.END.max()
            else:
                chr2centro[ch] 
    return chr2centro

def sample_genomic_interval(chr, start, end, sample, haploid_coverage, read_length=150,
                           return_fcn = False, fcn_only = False):    
    # Identify overlapping bins
    my_sample = sample[(sample.chr == chr) & (
        ((sample.start >= start) & (sample.start < end)) | 
        ((sample.end > start) & (sample.end <= end)) | 
        ((sample.start <= start) & (sample.end >= end))
    )].copy()
    assert len(my_sample) >= 1, (chr, start, end)

    # Compute overlap lengths between desired interval and genomic segments
    my_sample['overlap_length'] = my_sample.end - my_sample.start + 1
    first_overlap = min(end, my_sample.iloc[0].end) - max(start, my_sample.iloc[0].start)
    last_overlap = min(end, my_sample.iloc[-1].end) - max(start, my_sample.iloc[-1].start)
    assert first_overlap > 0, (chr, start, end, my_sample)
    assert last_overlap > 0, (chr, start, end, my_sample)
    my_sample.iloc[0, -1] = first_overlap
    my_sample.iloc[-1, -1] = last_overlap

    # Sample bin-covering reads 
    average_fcn = (my_sample.fcn_total * my_sample.overlap_length).sum() / my_sample.overlap_length.sum()
    if fcn_only:
        return average_fcn
    else:
        bin_reads =  poisson.rvs(average_fcn * haploid_coverage * (end - start) / read_length)

        # Sample threshold-covering reads
        head_reads = poisson.rvs(my_sample.iloc[0].fcn_total * haploid_coverage)
        tail_reads = poisson.rvs(my_sample.iloc[-1].fcn_total * haploid_coverage)

        if return_fcn:
            return head_reads, bin_reads, tail_reads, average_fcn
        else:
            return head_reads, bin_reads, tail_reads

def sample_raw_counts_chrom(samples, coverages, sample_names, chromosome, snps, cent_start, cent_end, bin_width):
    chr_end = samples[0][samples[0].chr == chromosome].end.max()
    
     ## Compute SNP-based thresholds
    # Identify thresholds for sub-bins (sampling intervals)
    p_snps = snps[snps.pos <= cent_start]
    p_snp_positions = np.array(sorted(p_snps.pos))
    
    # add 0.01 to ensure traditional rounding (0.5 goes up) rather than default rounding (closest even number)
    p_snp_thresholds = np.concatenate([np.trunc(np.vstack([p_snp_positions[:-1], p_snp_positions[1:]]).mean(axis=0)),
                                       [cent_start]]).astype(np.uint32)                                
    
    q_snps = snps[snps.pos >= cent_end]
    q_snp_positions = np.array(sorted(q_snps.pos))
    q_snp_thresholds = np.concatenate([[cent_end], np.trunc(np.vstack([q_snp_positions[:-1], q_snp_positions[1:]]).mean(axis=0)),
                                       [chr_end]]).astype(np.uint32)
    
    ## Combine with fixed-width thresholds
    p_thresholds = [int(a) for a in sorted(set(p_snp_thresholds).union(set([cent_start])).union(set(np.arange(0, cent_start, bin_width))))]
    q_thresholds = [int(a) for a in sorted(set(q_snp_thresholds).union(set([chr_end])).union(set(np.arange(cent_end, chr_end, bin_width))))]
    
    raw_rows = []
    for sample, coverage, sample_name in zip(samples, coverages, sample_names):
        for i in range(len(p_thresholds) - 1):
            start = p_thresholds[i]
            end = p_thresholds[i + 1]
            
            head_counts, bin_counts, tail_counts, average_fcn = sample_genomic_interval(chromosome, start, end, sample, coverage, return_fcn = True)
            
            raw_rows.append([chromosome, start, end, sample_name, head_counts, bin_counts, tail_counts, average_fcn])
    
    
    for sample, coverage, sample_name in zip(samples, coverages, sample_names):
        for i in range(len(q_thresholds) - 1):
            start = q_thresholds[i]
            end = q_thresholds[i + 1]
            
            head_counts, bin_counts, tail_counts, average_fcn = sample_genomic_interval(chromosome, start, end, sample, coverage, return_fcn = True)
            
            raw_rows.append([chromosome, start, end, sample_name, head_counts, bin_counts, tail_counts, average_fcn])
            
    raw_counts = pd.DataFrame(raw_rows, columns = ['chr', 'start', 'end', 'sample', 'head_counts', 'bin_counts', 'tail_counts', 'average_fcn'])
    return raw_counts, p_snp_thresholds, q_snp_thresholds

def mixture2name(mixture):
    '''
    Assumes that the normal sample mixture proportion is first, and the remainder are tumor clones
    numbered from 0
    '''
    normal_token = f'0{str(int(mixture[0] * 100)).rstrip("0")}normal'
    return 'bulk_' + '_'.join([f'0{str(int(a * 100)).rstrip("0")}clone{i}' for i, a in enumerate(mixture[1:])] +
                            [normal_token])

def aggregate_counts_hatchet1(raw_counts, sample_names, chromosome, cent_start, cent_end, chr_end, bin_width):
    """
    Aggregate raw counts into a table with fixed-width bins for HATCHet
    """
    h1_p_thresholds = np.arange(0, cent_start, bin_width)
    h1_q_thresholds = np.arange(cent_end, chr_end, bin_width)
    h1_rows = []
    for sample_name in sample_names:
        my_raw_counts = raw_counts[raw_counts['sample'] == sample_name]
        for i in range(len(h1_p_thresholds) - 1):
            start = h1_p_thresholds[i]
            end = h1_p_thresholds[i + 1]

            my_bins = my_raw_counts[(my_raw_counts.start >= start) & (my_raw_counts.end <= end - 1)]

            if start == 0:
                start = 1
            
            h1_rows.append([chromosome, int(start), int(end), sample_name, my_bins.bin_counts.sum()])
            
    h1_rdr = pd.DataFrame(h1_rows, columns = ['chr', 'start', 'end', 'sample', 'reads'])    
    return h1_rdr
        
def split_cn(row):
    # Convert from columns with CN strings to integerCN arrays 
    a = []
    b = []
    n_clones = max([i for i in range(20) if f'clone{i}' in row]) + 1
    for i in range(n_clones):
        tkns = row[f'clone{i}'].split('|')
        assert len(tkns) == 2
        a.append(int(tkns[0]))
        b.append(int(tkns[1]))
    return np.array([1] + a), np.array([1] + b)
        
def add_fcn(genome, mixture):
    # Add columns containing fractional copy numbers for each segment of the genome
    df = genome.copy()
    acn = []
    bcn = []
    for _, r in df.iterrows():
        a, b = split_cn(r)
        acn.append(np.sum(mixture * a))
        bcn.append(np.sum(mixture * b))
    df['fcn_a'] = acn
    df['fcn_b'] = bcn
    df['fcn_total'] = df.fcn_a + df.fcn_b
    return df
        
def aggregate_reads_hatchet1(raw_counts, sample_names, cent_start, cent_end, chromosome, bin_width):
    chr_end = raw_counts.end.max()
    h1_p_thresholds = np.arange(0, cent_start, bin_width)
    h1_q_thresholds = np.arange(cent_end, chr_end, bin_width)
    h1_rows = []
    for sample_name in sample_names:
        my_raw_counts = raw_counts[raw_counts['sample'] == sample_name]
        for i in range(len(h1_p_thresholds) - 1):
            start = h1_p_thresholds[i]
            end = h1_p_thresholds[i + 1]

            my_bins = my_raw_counts[(my_raw_counts.start >= start) & (my_raw_counts.end <= end - 1)]

            if start == 0:
                start = 1
            
            h1_rows.append([chromosome, int(start), int(end), sample_name, my_bins.bin_counts.sum()])
            
        for i in range(len(h1_q_thresholds) - 1):
            start = h1_q_thresholds[i]
            end = h1_q_thresholds[i + 1]

            my_bins = my_raw_counts[(my_raw_counts.start >= start) & (my_raw_counts.end <= end - 1)]

            if start == 0:
                start = 1
            
            h1_rows.append([chromosome, int(start), int(end), sample_name, my_bins.bin_counts.sum()])
    h1df = pd.DataFrame(h1_rows, columns = ['chromosome', 'start', 'end', 'sample', 'reads'])
    return h1df

def aggregate_reads_hatchet2(raw_counts, sample_names, p_snp_thresholds, q_snp_thresholds):    
    h2_p_thresholds = np.concatenate([[1], p_snp_thresholds])    
    h2_q_thresholds = q_snp_thresholds

    h2_p_arr = np.zeros((len(h2_p_thresholds), len(sample_names) * 2), dtype = int)
    for j, sample_name in enumerate(sample_names):
        my_raw_counts = raw_counts[raw_counts['sample'] == sample_name]
        for i in range(len(h2_p_thresholds) - 1):
            start = h2_p_thresholds[i]
            end = h2_p_thresholds[i + 1]

            my_bins = my_raw_counts[(my_raw_counts.start >= start - 1) & (my_raw_counts.end <= end)]
            assert len(my_bins) >= 1, (sample_name, start, end)

            total_counts = my_bins.bin_counts.sum()
            head_counts = my_bins.iloc[0].head_counts
            tail_counts = my_bins.iloc[-1].tail_counts
            h2_p_arr[i, j * 2] = total_counts
            h2_p_arr[i, j * 2 + 1] = head_counts # most of these will be overwritten which is fine
            h2_p_arr[i + 1, j * 2 + 1] = tail_counts

    h2_q_arr = np.zeros((len(h2_q_thresholds), len(sample_names) * 2), dtype = int)
    for j, sample_name in enumerate(sample_names):
        my_raw_counts = raw_counts[raw_counts['sample'] == sample_name]
        for i in range(len(h2_q_thresholds) - 1):
            start = h2_q_thresholds[i]
            end = h2_q_thresholds[i + 1]

            my_bins = my_raw_counts[(my_raw_counts.start >= start - 1) & (my_raw_counts.end <= end)]
            assert len(my_bins) >= 1, (sample_name, start, end)

            total_counts = my_bins.bin_counts.sum()
            head_counts = my_bins.iloc[0].head_counts
            tail_counts = my_bins.iloc[-1].tail_counts
            h2_q_arr[i, j * 2] = total_counts
            h2_q_arr[i, j * 2 + 1] = head_counts # most of these will be overwritten which is fine
            h2_q_arr[i + 1, j * 2 + 1] = tail_counts
    return h2_p_arr, h2_q_arr

def snp2fcn(r, mixture, genome):
    bin = genome[(genome.chr == r.chr) & (genome.start <= r.pos) & (genome.end >= r.pos)]
    assert len(bin) == 1, (r.chr, r.pos, len(bin))
    a, b = split_cn(bin.iloc[0])

    fcn_a = np.sum(mixture * a)
    fcn_b = np.sum(mixture * b)
    return f'{fcn_a}_{fcn_b}'

def sample_snps(chrpos, samples, mixtures, sample_names, coverages):
    # Sample SNP counts for each sample and bin
    l = [len(samples), len(mixtures), len(sample_names), len(coverages)]
    assert l[0] == l[1] == l[2] == l[3], l

    bafdfs = {}
    for (sample, mixture, sample_name, coverage) in zip(samples, mixtures, sample_names, coverages):
        # Compute FCN for each SNP conditioned on sample mixture proportions
        fcn_df = chrpos.apply(snp2fcn, mixture=mixture, genome=sample, axis=1).str.split('_', expand = True)
        bafdf = chrpos.copy()
        bafdf['sample'] = sample_name
        bafdf['fcn_a'] = fcn_df.iloc[:, 0].astype(float)
        bafdf['fcn_b'] = fcn_df.iloc[:, 1].astype(float)
        bafdf['total_fcn'] = bafdf.fcn_a + bafdf.fcn_b
        bafdf['true_baf'] = bafdf.fcn_b / (bafdf.fcn_a + bafdf.fcn_b)

        # Sample SNP counts
        total_reads = poisson.rvs(bafdf.total_fcn * coverage)
        b_reads = binom.rvs(total_reads, bafdf.true_baf)
        bafdf['ref'] = total_reads - b_reads
        bafdf['alt'] = b_reads
        bafdfs[sample_name] = bafdf.set_index(['chr', 'pos'])
        
    return bafdfs

@click.command()
@click.option('--genome', help='Simulated genome file path')
@click.option('--snps', help='SNP positions file path')
@click.option('--chromosome', help='Chromosome to analyze in this script')
@click.option('--sample_mixtures', help='Sample mixture proportions file path')
@click.option('--tumor_coverage', help='Expected coverage for tumor samples (default 80)', type=int, default=20)
@click.option('--normal_coverage', help='Expected coverage for normal sample (default 30)', type=int, default=15)
@click.option('--seed',  help='Random seed', default = 0)
@click.option('--genome_version', help = '"hg19" or "hg38"', default="hg19")
@click.option('--bin_width', help='Bin width', default = 100000)
@click.option('--normal_name', help='Name of normal sample', default = 'normal')
@click.option('--outdir', help='Output directory for intermediate files')
def main(genome, chromosome, snps, sample_mixtures, tumor_coverage, normal_coverage, genome_version, outdir, seed, bin_width, normal_name):
    np.random.seed(seed)
    
    genome = pd.read_table(genome)
    genome = genome[genome.chr == chromosome]
    assert len(genome) > 0, f'Found no bins in genome for chromosome: {chromosome}'
    
    centromeres_file = f'/n/fs/ragr-data/datasets/ref-genomes/centromeres/{genome_version}.centromeres.txt'
    chr2centro = load_centromeres(centromeres_file)
    
    snps = pd.read_table(snps, names=['chr', 'pos', 'sample', 'alt', 'ref'],
                         usecols = ['chr', 'pos', 'sample'],
                         dtype={'chr': object, 'pos': np.uint32, 'sample': object, 
                                      'ref': np.uint32, 'alt': np.uint32})
    snps = snps[snps.chr == chromosome]
    assert len(snps) > 0, f'Found no SNPs in SNPs table for chromosome: {chromosome}'
    chrpos = snps[['chr', 'pos']].drop_duplicates()
    
    mixtures = []
    with open(sample_mixtures, 'r') as f:
        for line in f:
            mix = np.array([float(a) for a in line.strip().split(',')])
            assert np.sum(mix) == 1, f"Mixture must sum to 1: {mix}"
            assert np.all(mix >= 0), f"Mixture must be non-negative: {mix}"
            mixtures.append(mix)
    sample_names = [mixture2name(m) for m in mixtures]
    
    sorted_order = np.argsort(sample_names)
    # hacky workaround to keep these as lists
    sn = []
    ms = []
    for idx in sorted_order:
        sn.append(sample_names[idx])    
        ms.append(mixtures[idx])
    sample_names = sn
    mixtures = ms

    # Generate sample dataframes
    samples = [add_fcn(genome, mixture) for mixture in mixtures]
    normal_sample = samples[0].copy()
    normal_sample.fcn_a = 1
    normal_sample.fcn_b = 1
    normal_sample.fcn_total = 2
    samples = [normal_sample] + samples
    sample_names = [normal_name] + sample_names

    # Compute expected sample haploid coverages
    sample_ploidies = []
    for sample in samples:
        sample['width'] = sample.end - sample.start
        sample_ploidies.append((sample.fcn_total * sample.width).sum() / sample.width.sum())
    
    expected_coverage = np.concatenate([[normal_coverage / sample_ploidies[0]],
                                 tumor_coverage / np.array(sample_ploidies[1:])])
    
    # Generate raw read counts
    raw_counts, p_snp_thresholds, q_snp_thresholds = sample_raw_counts_chrom(
        samples, expected_coverage, sample_names, chromosome, chrpos, *chr2centro[chromosome], bin_width)

    if not os.path.exists(os.path.join(outdir, 'raw_counts')):
        os.mkdir(os.path.join(outdir, 'raw_counts'))
    raw_counts.to_csv(os.path.join(outdir, 'raw_counts', f'{chromosome}.raw_counts.tsv'), sep = '\t', index = False)
    for sample, sample_name in zip(samples, sample_names):
        sample.to_csv(os.path.join(outdir, 'raw_counts', f'{chromosome}.{sample_name}.tsv'), sep = '\t', index = False)

    # Aggregate into fixed-width bins for HATCHet1
    h1df = aggregate_reads_hatchet1(raw_counts, sample_names, *chr2centro[chromosome], chromosome, bin_width)
    
    if not os.path.exists(os.path.join(outdir, 'hatchet1_rd')):
        os.mkdir(os.path.join(outdir, 'hatchet1_rd'))
    h1_normal = h1df[h1df['sample'] == normal_name].reset_index(drop = True)
    h1_normal.to_csv(os.path.join(outdir, 'hatchet1_rd', f'{chromosome}.normal.1bed'), sep = '\t', index = False, header = False)
    h1_tumor = h1df[h1df['sample'] != normal_name].reset_index(drop = True)
    h1_tumor.to_csv(os.path.join(outdir, 'hatchet1_rd', f'{chromosome}.tumor.1bed'), sep = '\t', index = False, header = False)

    if not os.path.exists(os.path.join(outdir, 'totalcounts')):
        os.mkdir(os.path.join(outdir, 'totalcounts'))
    raw_counts.groupby('sample').sum()['bin_counts'].to_csv(os.path.join(outdir, 'totalcounts', f'{chromosome}.total.tsv'), sep='\t', header = False)
    
    # Aggregate into SNP-based intervals for HATCHet2
    h2p, h2q = aggregate_reads_hatchet2(raw_counts, sample_names, p_snp_thresholds, q_snp_thresholds)
    my_thresholds = np.concatenate([[1], p_snp_thresholds, q_snp_thresholds]).astype(int)
    h2_arr = np.concatenate([h2p, h2q], axis = 0)
    
    if not os.path.exists(os.path.join(outdir, 'hatchet2_rd')):
        os.mkdir(os.path.join(outdir, 'hatchet2_rd'))
    np.savetxt(os.path.join(outdir, 'hatchet2_rd', f'{chromosome}.total.gz'), h2_arr, fmt='%s')
    np.savetxt(os.path.join(outdir, 'hatchet2_rd', f'{chromosome}.thresholds.gz'), my_thresholds, fmt='%s')
    with open(os.path.join(outdir, 'hatchet2_rd', 'samples.txt'), 'w') as f:
        for sample in sample_names:
            f.write(sample + '\n')

    # Sample SNP counts (same for both methods)
    bafdfs = sample_snps(chrpos, samples, [[1] + [0] * (len(mixtures[0]) - 1)] + mixtures, sample_names, expected_coverage)
    normal_df = bafdfs[normal_name]\
        
    if not os.path.exists(os.path.join(outdir, 'baf')):
        os.mkdir(os.path.join(outdir, 'baf'))
    tumor_df = pd.concat([df for sample_name, df in bafdfs.items() if sample_name != normal_name])
    normal_df.to_csv(os.path.join(outdir, 'baf', f'{chromosome}.normal_baf.1bed'), sep = '\t')
    tumor_df.to_csv(os.path.join(outdir, 'baf', f'{chromosome}.tumor_baf.1bed'), sep = '\t')     


if __name__ == '__main__':
    main()
