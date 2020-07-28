# %%

from analyze import get_sweep_data, write_sweep_data

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

sweep = get_sweep_data(indata_dir, params)

# %%

# sweep = sweep.head(1000)
# outdata_dir = 'sample-data/analyzed-data'
outdata_dir = 'data/analyzed-data'

gene_info = write_sweep_data(outdata_dir, sweep, params)

# %%
