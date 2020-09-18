# %%
from importlib import reload
import sweeptools
from sweeptools.analyzesweep import get_sweep_data, write_sweep_data

# Define parameters of screen to read
params = {'screen_name': 'H3K27-but',
          'assembly': 'hg38',
          'trim_length': '50',
          'mode': 'collapse',
          'start': 'tx',
          'end': 'tx',
          'overlap': 'both',
          'direction': 'sense',
          'step': 500}

indata_dir = '../data/sweeps-screen-analyzer'

sweep = get_sweep_data(indata_dir, params)

# %%

outdata_dir = '../data/sweeps-analyzed'
gene_info = write_sweep_data(outdata_dir, sweep, params)

# %%
# from importlib import reload
# import sweeptools
# reload(sweeptools)
# from sweeptools.analyzesweep import read_analyzed_sweep

# data_dir = 'sample-data/analyzed-data'
# a = read_analyzed_sweep(data_dir, params)


# %%
