# %%
from importlib import reload
import tools
from tools.analyzesweep import get_sweep_data, write_sweep_data

# Define parameters of screen to read
params = {'screen_name': 'Ac-gamma-Actin',
          'assembly': 'hg38',
          'trim_length': '50',
          'mode': 'collapse',
          'start': 'tx',
          'end': 'tx',
          'overlap': 'both',
          'direction': 'sense',
          'step': 500}

# indata_dir = 'sample-data/screen-analyzer-data'
# indata_dir = 'data/screen-analyzer-data'
indata_dir = '../data/sweeps-screen-analyzer'

sweep = get_sweep_data(indata_dir, params)

# %%

# reload(tools)
# from tools.analyzesweep import write_sweep_data

# sample_sweep = sweep.query('gene_name == "CD274"')
# sample_sweep = sweep.head(100000)
# outdata_dir = 'sample-data/analyzed-data'
# gene_info = write_sweep_data(outdata_dir, sample_sweep, params)

# outdata_dir = 'data/analyzed-data'
outdata_dir = '../data/sweeps-analyzed'
gene_info = write_sweep_data(outdata_dir, sweep, params)

# %%
# from importlib import reload
# import tools
# reload(tools)
# from tools.analyzesweep import read_analyzed_sweep

# data_dir = 'sample-data/analyzed-data'
# a = read_analyzed_sweep(data_dir, params)


# %%
