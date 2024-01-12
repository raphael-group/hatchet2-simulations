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
    print(clone_cols)
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
            prohib_states = df[(df.start == start) | (df.end == start)][clone_cols].values
            assert len(prohib_states) == 2 or start == 0

            print(prohib_states)
            while True:
                state_idx = np.random.randint(len(unique_ascns))
                state = list(unique_ascns[state_idx])
                print(state)

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

def simplify_genome(genome, prop_remove = 0.5):
    # Simplify the genome by removing a proportion of segments and extending neighboring segments
    cn_cols = genome.columns[3:]

    new_genome = []
    for ch, df in genome.groupby('chr'):
        if len(df) <= 3:
            new_genome.append(df)
        else:
            to_remove = sorted(np.random.choice(np.arange(len(df)), replace=False, size=int(len(df) * prop_remove)))

            new_segments = []
            to_delete = set()
            for i in to_remove:
                if i in to_delete:
                    continue

                df = df.reset_index(drop = True)
                to_delete.add(i)
                current = df.iloc[i]
                if i < len(df) - 1 and not all(j in to_delete for j in range(i, len(df))):
                    j = i + 1
                    while j in to_delete:
                        # keep looking later
                        j += 1

                    # if it's not the last bin, extend the following segment
                    successor = df.iloc[j]
                    new_start = current.start
                    new_end = successor.end
                    new_states = successor[cn_cols]
                    to_delete.add(j)
                else:
                    # extend the preceding segment
                    j = i - 1
                    while j in to_delete:
                        # keep looking earlier
                        j -= 1

                    predecessor = df.iloc[j]
                    new_start = predecessor.start
                    new_end = current.end
                    new_states = predecessor[cn_cols]
                    to_delete.add(j)
                new_segments.append([ch, new_start, new_end, *new_states])
            new_chromosome = pd.concat([df.drop(index=list(to_delete)), 
                   pd.DataFrame(new_segments, columns = df.columns)]).sort_values(by = ['chr', 'start']).reset_index(drop = True)
            
            # check for adjacent segments with the same states and merge segments
            merge_sets = []
            for i in range(len(new_chromosome) - 1):
                my_set = set()
                j = i + 1
                while j < len(new_chromosome) and np.all(new_chromosome.iloc[i][cn_cols] == new_chromosome.iloc[j][cn_cols]):
                    my_set.add(i)
                    my_set.add(j)
                    j += 1
                if len(my_set) > 0:
                    merge_sets.append(my_set)

            collapsed_segments = []
            to_delete = set()
            for s in merge_sets:
                to_delete = to_delete.union(s)
                s = sorted(s)
                collapsed_segments.append([ch, new_chromosome.iloc[s[0]].start, new_chromosome.iloc[s[-1]].end, 
                               *new_chromosome.iloc[s[0]][cn_cols]])
            new_chromosome = pd.concat([new_chromosome.drop(index=list(to_delete)), 
                   pd.DataFrame(collapsed_segments, columns = df.columns)]).sort_values(by = ['chr', 'start'])
            new_genome.append(new_chromosome)
            
    return pd.concat(new_genome).reset_index(drop = True)

@click.command()
@click.argument('orig_genome')
@click.option('--event_sizes', help = 'Comma-separated list of focal event sizes to inject')
@click.option('--n_clones', help = 'Number of tumor clones (default: use number in tumor genome)', default=-1, type=int)
@click.option('--prop_mirrored', help = 'Proportion of ALL segments that are mirrored', type=float, default=0.2)
@click.option('--prop_focal_mirrored', help = 'Proportion of injected focal segments that are mirrored', type=float, default=0.5)
@click.option('--genome_version', help = '"hg19" or "hg38"', default="hg19")
@click.option('--prop_simplify', default=0, type = float)
@click.option('--seed', help = 'Random seed', type=int, default=0)
@click.option('--output', '-o', help='Output filename')
def main(orig_genome, n_clones, event_sizes, prop_mirrored, prop_focal_mirrored, prop_simplify, genome_version, seed, output):
    genome = pd.read_table(orig_genome).rename(columns={'#CHR':'chr', 'START':'start', 'END':'end'})
    if n_clones is not None and n_clones > 0:
        assert n_clones <= len(genome.columns) - 3
        genome = genome.iloc[:, :n_clones + 3]
    event_sizes = [int(a) for a in event_sizes.split(',')]
    
    centromeres_file = f'/n/fs/ragr-data/datasets/ref-genomes/centromeres/{genome_version}.centromeres.txt'
    chr2centro = load_centromeres(centromeres_file)
    
    if prop_simplify > 0:
        genome = simplify_genome(genome, prop_remove = prop_simplify)
    result = inject_events(genome, event_sizes, chr2centro,
                           seed=seed, prop_mirrored=prop_mirrored,
                           prop_focal_mirrored=prop_focal_mirrored)
    result.to_csv(output, index = False, sep = '\t')

if __name__ == '__main__':
    main()
