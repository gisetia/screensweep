# %%
from sweeptools.analyzesweep import get_sweep_data, write_sweep_data
from importlib import reload
import sweeptools
reload(sweeptools)

# Define parameters of screen to read
params = {'screen_name': 'LAD5-GFP_1',
          'assembly': 'hg38',
          'trim_length': '50',
          'mode': 'collapse',
          'start': 'tx',
          'end': 'tx',
          'overlap': 'both',
          'direction': 'sense',
          'step': 500}

indata_dir = '../data/sweeps-screen-analyzer_2020-09-21'
# indata_dir = '../data/sweeps-screen-analyzer_2020-06-25'

print(f'Reading sweep for screen {params["screen_name"]}')
sweep = sweeptools.analyzesweep.get_sweep_data(indata_dir, params)


# %%
reload(sweeptools)

outdata_dir = '../data/sweeps-analyzed_2020-09-21'
gene_info = sweeptools.analyzesweep.write_sweep_data(outdata_dir, sweep,
                                                     params)

# %%
# from importlib import reload
# import sweeptools
# reload(sweeptools)
# from sweeptools.analyzesweep import read_analyzed_sweep

# data_dir = 'sample-data/analyzed-data'
# a = read_analyzed_sweep(data_dir, params)


# %%
