import click
import pandas as pd
import numpy as np


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

def inject_events(genome, event_sizes, chr2centro, seed = 0,
                 prop_mirrored = 0.2, prop_focal_mirrored = 0.5):
    
    # determine where events should be placed
    np.random.seed(seed)

    chroms = sorted(genome.chr.unique())
    n_chroms = len(chroms)
    chrom_placement = np.random.randint(0, n_chroms, size = len(event_sizes))
    event_sizes = np.array(event_sizes)
    n_clones = max([i for i in range(50) if f'clone{i}' in genome.columns]) + 1
    clone_cols = [f'clone{i}' for i in range(n_clones)]

    unique_ascns = genome[clone_cols].value_counts().keys().values
    genome_by_chrom = {ch:df.reset_index(drop = True) for ch, df in genome.groupby('chr')}
    new_genome = {}

    # spike in new events
    for ch in chroms:
        cent_start, cent_end = chr2centro[ch]
        my_event_indices = np.where(chrom_placement == chroms.index(ch))[0]
        df = genome_by_chrom[ch].copy()

        for size in event_sizes[my_event_indices]:
            # insert an event of this size
            thresholds = sorted(np.concatenate([df.start.values, [df.end.max()]]))

            # insert new event as the prefix of an existing segment
            while True:
                start_idx = np.random.randint(len(thresholds) - 1)
                start = thresholds[start_idx]
                end = thresholds[start_idx + 1]
                new_end = int(start + size)

                # reroll if it overlaps the centromere, or the new event doesn't fit within the segment
                if not ((start >= cent_start and start <= cent_end) or 
                   (new_end >= cent_start and new_end <= cent_end)) and new_end < end:
                    break

            # randomly sample cn state for the event
            prohib_states = df[(df.start == start) | (df.end == start)][['clone0', 'clone1']].values
            assert len(prohib_states) == 2 or start == 0

            while True:
                state_idx = np.random.randint(len(unique_ascns))
                state = list(unique_ascns[state_idx])

                # reroll if it's the same as its housing segment or neighboring segment
                if state not in prohib_states:
                    break

            # randomly mirror event w.p. prop_mirrored        
            if np.random.random() < prop_focal_mirrored:
                clone_to_flip = np.random.randint(n_clones)
                state[clone_to_flip] = state[clone_to_flip][::-1]

            new_rows = []

            # add a row for the new event
            new_rows.append([ch, start, new_end, *state])

            # add a row for the housing segment that had its prefix removed
            new_rows.append([ch, new_end, end, *prohib_states[0]])

            # delete old row
            to_delete = df[df.start == start].index
            assert len(to_delete) == 1
            df = df.drop(index = to_delete)

            df = pd.concat([df, pd.DataFrame(new_rows, columns = df.columns)]).reset_index(drop = True)
            thresholds = sorted(thresholds + [int(start + size)])

        # randomly mirror all events on the chromosome
        for i in df.index:
            if np.random.random() < prop_mirrored:                
                clone_to_flip = np.random.randint(n_clones)
                df.loc[i, f'clone{clone_to_flip}'] =  df.loc[i, f'clone{clone_to_flip}'][::-1]

        # store modified genome
        new_genome[ch] = df.sort_values(by = 'start').reset_index(drop = True)
        
    return pd.concat(new_genome.values()).reset_index(drop = True)

@click.command()
@click.argument('orig_genome')
@click.option('--event_sizes', help = 'Comma-separated list of focal event sizes to inject')
@click.option('--prop_mirrored', help = 'Proportion of ALL segments that are mirrored', type=float, default=0.2)
@click.option('--prop_focal_mirrored', help = 'Proportion of injected focal segments that are mirrored', type=float, default=0.5)
@click.option('--genome_version', help = '"hg19" or "hg38"', default="hg19")
@click.option('--seed', help = 'Random seed', type=int, default=0)
@click.option('--output', '-o', help='Output filename')
def main(orig_genome, event_sizes, prop_mirrored, prop_focal_mirrored, genome_version, seed, output):
    genome = pd.read_table(orig_genome).rename(columns={'#CHR':'chr', 'START':'start', 'END':'end'})
    event_sizes = [int(a) for a in event_sizes.split(',')]
    
    centromeres_file = f'/n/fs/ragr-data/datasets/ref-genomes/centromeres/{genome_version}.centromeres.txt'
    chr2centro = load_centromeres(centromeres_file)
    
    result = inject_events(genome, event_sizes, chr2centro,
                           seed=seed, prop_mirrored=prop_mirrored,
                           prop_focal_mirrored=prop_focal_mirrored)
    result.to_csv(output, index = False, sep = '\t')

if __name__ == '__main__':
    main()
