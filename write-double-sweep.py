# %%
import os
from importlib import reload
import sweeptools
# reload(sweeptools)

# Define parameters of screen to read
params = {'screen_name': '',
          'assembly': 'hg38',
          'trim_length': '50',
          'mode': 'collapse',
          'start': 'tx',
          'end': 'tx',
          'overlap': 'both',
          'direction': 'sense',
          'step': 500}

indata_dir = '../data/sweeps/sweeps-screen-analyzer_2020-09-21'
# indata_dir = '../data/sweeps/sweeps-screen-analyzer_2020-06-25'
outdata_dir = '../data/sweeps/sweeps-analyzed_2020-09-21'

in_files = set(os.listdir(indata_dir))
out_files = set(os.listdir(outdata_dir))
missing_screens = list(in_files.difference(out_files))
missing_screens.sort()

screen_list = missing_screens

# %%

for idx, screen in enumerate(screen_list):

    params['screen_name'] = screen
    print(f'\nReading sweep for screen {params["screen_name"]} - {idx+1} '
          f'of {len(screen_list)}')

    sweep = sweeptools.analyzesweep.get_sweep_data(indata_dir, params)

    gene_info = sweeptools.analyzesweep.write_sweep_data(outdata_dir, sweep,
                                                         params)

# %%
