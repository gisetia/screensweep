# %%

# Define parameters of screen to read
params = {'screen_name': 'PDL1_IFNg',
          'assembly': 'hg38',
          'trim_length': '50',
          'mode': 'collapse',
          'start': 'tx',
          'end': 'tx',
          'overlap': 'both',
          'direction': 'sense',
          'step': 1000}

# data_dir = 'sample_analyzed'
data_dir = 'analyzed'

import re
import pandas as pd
import time

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

    return sweep_data

cols = ['gene_name', 'p_fdr', 'log2_mi', 'srt_off', 'end_off']
sweep = get_sweep_data(data_dir, params)[cols]
# %%

slope_thr = 2
p_thr = 0.001

genes = sweep.gene_name.unique()

start_time = time.time()

genes = genes[0:1000]
l_genes = len(genes)
flagged_genes = []
for idx, g in enumerate(genes):
# for g in genes[0:100]:
# for g in ['CD274']:

    if not idx % 50:
        print(f'Analyzing gene {idx+1} of {l_genes}.')

    gene = sweep.query('gene_name == @g')
    gene = gene.pivot(index='srt_off', columns='end_off', 
                      values=['p_fdr','log2_mi'])

    # Get slopes when changing parameters in both directions 
    # (delta log2 MI per 1,000 bp)
    slope_sdir = (gene['log2_mi'] 
                  - gene['log2_mi'].shift(1))/(params['step']*0.001)
    slope_edir = (gene['log2_mi'] 
                  - gene['log2_mi'].shift(1, axis=1))/(params['step']*0.001)

    gene = gene.stack()
    gene['sl_sdir'] = slope_sdir.stack()
    gene['sl_edir'] = slope_edir.stack()

    flags = gene.query('(sl_sdir > @slope_thr | sl_edir > @slope_thr) '
                       '& p_fdr < @p_thr')

    if not flags.empty:
        flagged_genes.append(g)

print('--- %s seconds ---\n' % (time.time() - start_time))

print(flagged_genes)

# %%
