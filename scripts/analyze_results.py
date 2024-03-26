import click
import pandas as pd
import numpy as np
import os 
from collections import defaultdict
import json

def sample2props(samplename):
    tkns = samplename.split('_')
    tumor_props = []

    for tkn in tkns[1:]:
        num = []
        for c in tkn:
            if c.isdigit():
                num.append(c)
            else:
                break
        if len(num) == 0:
            assert tkn.startswith('None'), tkn
            num = '0000'
            val = 0
        elif len(num) == 1 and num[0] == '0':
            val = 0
        else:
            try:
                val = np.round(int(''.join(num[1:])) * 10 ** (-1 * (len(num) - 1)), 3)
            except ValueError as e:
                print(samplename, tkn, num)
                raise e
        if tkn[len(num):].startswith('clone'):
            tumor_props.append(val)
        else:
            normal_prop = val
    return [normal_prop] + tumor_props

def ascn_error_per_base(gt, inf, max_gt_segment_size=np.inf, mirrored_only = False, debug = False):
    if mirrored_only:
        gt = gt[gt.is_mirrored]
        if len(gt) == 0:
            return -1
            raise ValueError("Flag 'mirrored_only' was set but there are no mirrored segments in this genome")
    chromosomes = gt['#CHR'].unique()
    n_clones_inf = max([i for i in range(20) if f'cn_clone{i}' in inf])

    errsum = 0
    errdenom = 0

    for ch in chromosomes:
        gt_ = gt[gt['#CHR'] == ch]
        inf_ = inf[inf['#CHR'] == ch]

        for _, r in gt_.iterrows():
            if r.END - r.START + 1 > max_gt_segment_size:
                continue
            
            my_inf = inf_[((inf_.START >= r.START + 1) & (inf_.START < r.END + 1)) | 
                          ((inf_.END >= r.START + 1) & (inf_.END < r.END + 1))]

            for _, r2 in my_inf.iterrows():
                true_u = sample2props(r2.SAMPLE)
                if r2.START >= r.START + 1 and r2.END < r.END + 1:
                    # inferred segment is entirely contained within gt segment
                    overlap_length = r2.END - r2.START

                elif r2.END >= r.END + 1:
                    # prefix of inferred segment is in gt segment
                    overlap_length = (r.END + 1) - r2.START

                elif r2.START < r.START + 1:
                    # suffix of inferred segment is in gt segment
                    overlap_length = r2.END - (r.START + 1)
                else:
                    print("Correspondence unclear between bins:", (r.START, r.END), (r2.START, r2.END))
                    break

                true_cn_props = defaultdict(lambda:0)
                for i,p in enumerate(true_u[1:]):
                    true_cn_props[r[f'clone{i}']] += p 
                true_cn_props['1|1'] += true_u[0]
                
                inf_cn_props = defaultdict(lambda:0)
                for i in range(n_clones_inf):
                    inf_cn_props[r2[f'cn_clone{i+1}']] += r2[f'u_clone{i+1}']
                inf_cn_props['1|1'] += r2.u_normal

                if debug:
                    print(true_cn_props.keys(), inf_cn_props.keys())
                
                biggest_diff = 0
                for state, prop1 in true_cn_props.items():
                    prop2 = inf_cn_props[state] if state in inf_cn_props else 0
                    my_diff = abs(prop1 - prop2)
                    if my_diff > biggest_diff:
                        biggest_diff = my_diff

                for state, prop2 in inf_cn_props.items():
                    prop1 = true_cn_props[state] if state in true_cn_props else 0
                    my_diff = abs(prop1 - prop2)
                    if my_diff > biggest_diff:
                        biggest_diff = my_diff      

                errsum += biggest_diff * overlap_length
                errdenom += overlap_length
    return errsum / errdenom

def resolve_segments(gt, inf):
    chromosomes = gt['#CHR'].unique()
    n_clones_inf = max([i for i in range(20) if f'cn_clone{i}' in inf])
    n_clones_gt = max([i for i in range(20) if f'clone{i}' in gt]) + 1

    errsum = 0
    errdenom = 0

    my_rows = []
    missing_calls = []
    for ch in chromosomes:
        gt_ = gt[gt['#CHR'] == ch]
        inf_ = inf[inf['#CHR'] == ch]

        for _, r in gt_.iterrows():
            my_inf = inf_[((inf_.START >= r.START + 1) & (inf_.START < r.END + 1)) | 
                          ((inf_.END >= r.START + 1) & (inf_.END < r.END + 1))]

            total_overlap = 0
            for _, r2 in my_inf.iterrows():
                true_u = sample2props(r2.SAMPLE)
                
                if r2.START >= r.START + 1 and r2.END < r.END + 1:
                    # inferred segment is entirely contained within gt segment
                    overlap_length = r2.END - r2.START

                elif r2.END >= r.END + 1:
                    # prefix of inferred segment is in gt segment
                    overlap_length = (r.END + 1) - r2.START

                elif r2.START < r.START + 1:
                    # suffix of inferred segment is in gt segment
                    overlap_length = r2.END - (r.START + 1)
                else:
                    print("Correspondence unclear between bins:", (r.START, r.END), (r2.START, r2.END))
                    break

                    
                true_cn_props = defaultdict(lambda:0)
                for i,p in enumerate(true_u[1:]):
                    true_cn_props[r[f'clone{i}']] += p 
                true_cn_props['1|1'] += true_u[0]

                inf_cn_props = defaultdict(lambda:0)
                for i in range(n_clones_inf):
                    inf_cn_props[r2[f'cn_clone{i+1}']] += r2[f'u_clone{i+1}']
                inf_cn_props['1|1'] += r2.u_normal

                my_row = [r2['#CHR'], r2.START, r2.END, overlap_length, r2.SAMPLE]
                my_row += list(r2.iloc[4:])
                my_row += ['1|1', true_u[0]]
                for i in range(n_clones_gt):
                    my_row += [r[f'clone{i}'], true_u[i + 1]]
                my_row.append(r.END - r.START)
                my_row.append(r.is_mirrored)
                
                my_rows.append(my_row)
                total_overlap += overlap_length
                
            nocall_bases = (r.END - r.START) - total_overlap
            if nocall_bases > 0:
                missing_calls.append([r['#CHR'], r.START, r.END, nocall_bases])

    header = ['#CHR', 'START', 'END', 'OVERLAP', 'SAMPLE']
    header += list(inf.columns[4:])
    header += ['gt_cn_normal', 'gt_u_normal']
    for i in range(n_clones_gt):
        header += [f'gt_cn_clone{i + 1}', f'gt_u_clone{i + 1}']  
    header.append('gt_segment_size')
    header.append('gt_is_mirrored')
    my_df = pd.DataFrame(my_rows, columns = header)
    return my_df, pd.DataFrame(missing_calls, columns = header[:3] + ['MISSING'])

def precision_recall(joint_seg, mirrored_only = False, max_gt_segment_size = np.inf):
    if mirrored_only:
        joint_seg = joint_seg[joint_seg.gt_is_mirrored]
        if len(joint_seg) == 0:
            return -1, -1, -1
            raise ValueError("Flag 'mirrored_only' was set but there are no mirrored segments in this genome")
            
    n_clones_inf = max([i for i in range(20) if f'cn_clone{i}' in joint_seg])
    n_clones_gt = max([i for i in range(20) if f'gt_cn_clone{i}' in joint_seg])

    p_sum = 0
    
    r_sum = 0
    total_length = 0
    
    acc_sum = 0
    
    for _, r in joint_seg.iterrows():
        if r.gt_segment_size <= max_gt_segment_size:
            inf_cns = set([r[f'cn_clone{i + 1}'] for i in range(n_clones_inf)])
            gt_cns = set([r[f'gt_cn_clone{i + 1}'] for i in range(n_clones_gt)])

            recovered = len(inf_cns.intersection(gt_cns))
            all_states = len(inf_cns.union(gt_cns))
            p_denom = len(inf_cns)
            r_denom = len(gt_cns)

            p_sum += r.OVERLAP * (recovered / p_denom)
            r_sum += r.OVERLAP * (recovered / r_denom)
            acc_sum += r.OVERLAP * (recovered / all_states)

            total_length += r.OVERLAP

    
    return p_sum / total_length, r_sum / total_length, acc_sum / total_length

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
        
def is_mirrored(row):
    acn, bcn = split_cn(row)
    agreater = acn > bcn
    bgreater = bcn > acn
    return agreater.any() and bgreater.any()

@click.command()
@click.argument('results_file')
@click.option('--simdata_dir', help = 'Path to simulated data', default='/n/fs/ragr-research/projects/hatchet2-results/newsims/data/')
@click.option('--joint_df_out', help='Filename of output joint dataframe')
@click.option('--stats_out', help='Filename of output metrics file')
def main(results_file, simdata_dir, joint_df_out, stats_out):
    # read results
    bbc = pd.read_table(results_file)
    
    # read in ground truth simulated data
    simkey = results_file.split(os.sep)[-3]
    
    #genome = pd.read_table(os.path.join(simdata_dir, f'simulated_genome{simkey}.tsv')).rename(
    #    columns = {'chr':'#CHR', 'start':'START', 'end':'END'})
    genome = pd.read_table(f'/n/fs/ragr-research/projects/hatchet2-results/newsims/mirrored6/genomes/{simkey[:-4]}.tsv').rename(
                                            columns = {'chr':'#CHR', 'start':'START', 'end':'END'})
    
    genome['is_mirrored'] = genome.apply(is_mirrored, axis = 1)

    # construct joint df
    joint = resolve_segments(genome, bbc)[0]
    joint.to_csv(joint_df_out, index = False, sep = '\t') 
    
    # compute metrics
    e = ascn_error_per_base(genome, bbc)
    pra = precision_recall(joint)
    e_mirrored = ascn_error_per_base(genome, bbc, mirrored_only = True)
    pra_mirrored = precision_recall(joint, mirrored_only = True)
    e_small = ascn_error_per_base(genome, bbc, max_gt_segment_size = 1e6)
    pra_small = precision_recall(joint, max_gt_segment_size = 1e6)

    #s25177s025_liquid_3clones_purity5_eventsA_0
    '''
    tkns = [a for a in simkey.split('_') if len(a) > 0] # account for double-underscore and single-underscore
    if tkns[0][0] == 'n':
        tkns = tkns[1:]
    dataset_id = tkns[0]
    liquidsolid = tkns[1]    
    n_clones = int(tkns[2][0])
    if tkns[3].startswith('purity'):
        mixture_id = int(tkns[3][-1])
        n_samples = 2
        events_id = tkns[4][-1]
        seed = int(tkns[5])
    else:
        n_samples = int(tkns[3][0])
        mixture_id = int(tkns[4])
        events_id = tkns[5][-1]
        seed = int(tkns[6])
    '''
    
    '''
    #dataset_n3_s4467_2_2
    tkns = simkey.split('_')
    dataset_id = simkey[2]
    liquidsolid = ''
    n_clones = 2
    mixture_id = int(tkns[3])
    seed = int(results_file.split(os.sep)[-3][-1])
    n_samples = 3 if mixture_id == 1 else 2
    events_id = 'A'
    '''
    #dataset_n2_tetraploid_2_0_1
    tkns = simkey.split('_')
    dataset_id = '_'.join(tkns[2:4])
    liquidsolid = ''
    n_clones = 2
    mixture_id = int(tkns[4])
    seed = int(tkns[5])
    n_samples = 3 if mixture_id == 1 else 2
    events_id = 'A'
    
    method = results_file.split(os.sep)[-4]
    
    result = {'simkey':simkey,
              'dataset_id':dataset_id,
                'liquidsolid':liquidsolid,
                'n_clones':n_clones,
                'n_samples':n_samples,
                'mixture_id':mixture_id,
                'events_id':events_id,
                'seed':seed,
                'error':e,
                'precision':pra[0],
                'recall':pra[1],
                'accuracy':pra[2],
                'error_mirrored':e_mirrored,
                'precision_mirrored':pra_mirrored[0],
                'recall_mirrored':pra_mirrored[1],
                'accuracy_mirrored':pra_mirrored[2],
                'error_small':e_small,
                'precision_small':pra_small[0],
                'recall_small':pra_small[1],
                'accuracy_small':pra_small[2],
                'method':method
                }
    json.dump(result, open(stats_out, 'w'))

if __name__ == '__main__':
    main()
