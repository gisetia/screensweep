# %%
import re
import os
from math import log10
import pandas as pd
from utils import timer

# Define parameters of screen to read
params = {'screen_name': 'PDL1_IFNg',
        'assembly': 'hg38',
        'trim_length': '50',
        'mode': 'collapse',
        'start': 'tx',
        'end': 'tx',
        'overlap': 'both',
        'direction': 'sense',
        'step': 500}

# indata_dir = 'sample-data/screen-analyzer-data'
indata_dir = 'data/screen-analyzer-data'

@timer
def get_sweep_data(data_dir: str, params: dict) -> pd.DataFrame:
    '''Get dataframe with all data resulting from a double parameter sweep of
    'screen-analyzer analyze'.
    ( dataframe columns are gene_name, low_counts, high_counts, p, fdr_p,
    log2_mi, 'start', position, offset)
    '''

    data_path = (f'''{data_dir}/{params['screen_name']}/'''
                f'''{params['assembly']}/{params['trim_length']}/'''
                f'''mode={params['mode']}_direction={params['direction']}'''
                f'''_overlap={params['overlap']}/double-sweep'''
                f'''_step={params['step']}/''')

    files = os.listdir(data_path)
    df_list = []  # Initialize list to append data
    for filename in files:

        if filename.startswith('out'):
            # Get param values from filename
            match = re.search(rf'_start=(.*?)_', filename)
            start = match.group(1)
            match = re.search(rf'_end=(.*?)_', filename)
            end = match.group(1)
        
            # File to df, add columns for each param value and add to list
            df = pd.read_csv(data_path + filename, sep='\t',
                                index_col=None, header=None)
            df[6] = start
            df[7] = end
            df_list.append(df)

    # Create single dataframe from list
    sweep_data = pd.concat(df_list, axis=0, ignore_index=True)
    sweep_data = sweep_data.rename(columns={0: 'gene_name', 1: 'low_counts',
                                            2: 'high_counts', 3: 'p',
                                            4: 'p_fdr', 5: 'log2_mi',
                                            6: 'start', 7: 'end'})

    # Split parameter values into txt and numbers, creating one column for each
    # eg. 'tx-100' is split in column 'position' with value 'tx':str and column
    # 'offset' with value -100:int
    sweep_data[['srt_pos', 'end_pos']] = sweep_data[['start', 'end']].applymap(
        lambda x: re.split('(-|\\+)', x)[0])
    sweep_data[['srt_off', 'end_off']] = sweep_data[['start', 'end']].applymap(
        lambda x: ''.join(re.split('(-|\\+)', x)[1:3])).astype(int)

    # Sort by name and offsets
    sweep_data = sweep_data.sort_values(by=['gene_name', 'srt_off', 'end_off'])

    # sweep_data.to_csv(f'{data_path}/all.csv')

    return sweep_data

sweep = get_sweep_data(indata_dir, params)

# %%

# sweep = sweep.head(10000)
# outdata_dir = 'sample-data/analyzed-data'
outdata_dir = 'data/analyzed-data'

@timer
def get_gene_info(data_dir: str, sweep_data: pd.DataFrame, 
                  params: dict) -> dict:

    data_path = (f'''{data_dir}/{params['screen_name']}/'''
                f'''{params['assembly']}/{params['trim_length']}/'''
                f'''mode={params['mode']}_direction={params['direction']}'''
                f'''_overlap={params['overlap']}/double-sweep'''
                f'''_step={params['step']}/''')
    
    if not os.path.exists(data_path):
        os.makedirs(data_path)
        print('Creating analyzed directory.')
    
    gene_info = dict()

    grouped = sweep_data.groupby('gene_name')
    for name, group in grouped:

        print(name)
        
        gene_data = group.pivot(index='srt_off', columns='end_off', 
                    values=['gene_name', 'low_counts', 'high_counts', 'p',
                            'p_fdr','log2_mi', 'srt_pos', 'end_pos'])

        # Get slopes of log2 MI when changing parameters in both directions 
        # (delta log2 MI per 1,000 bp)
        slope_sdir = (gene_data['log2_mi'] 
                      - gene_data['log2_mi'].shift(1))/(params['step']*0.001)

        slope_edir = (gene_data['log2_mi'] 
                      - gene_data['log2_mi'].shift(1, 
                                             axis=1))/(params['step']*0.001)

        # Get lof of ratio of p values when changing parameters in 
        # both directions 
        gene_data['p_fdr'] = gene_data['p_fdr'].replace(0, 1e-300)
        p_ratio_sdir = (gene_data['p_fdr']/
                        gene_data['p_fdr'].shift(1)).applymap(log10)
        p_ratio_edir = (gene_data['p_fdr']/
                        gene_data['p_fdr'].shift(1, axis=1)).applymap(log10)

        gene_data = gene_data.stack()
        gene_data['sl_sdir'] = slope_sdir.stack()
        gene_data['sl_edir'] = slope_edir.stack()
        gene_data['p_ratio_sdir'] = p_ratio_sdir.stack()
        gene_data['p_ratio_edir'] = p_ratio_edir.stack()
    
        gene_info[name] = gene_data

        # gene_data.to_csv(f'{data_path}{name}.csv')
    all_info = pd.concat(gene_info.values(), ignore_index=False)
    all_info.to_csv(f'{data_path}all_gene_info.csv')
    all_info.to_csv(f'{data_path}all_gene_info.gz', 
                    compression='gzip')

    return gene_info

gene_info = get_gene_info(outdata_dir, sweep, params)

# %%
