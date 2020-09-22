# %%
import os
from sweeptools.analyzesweep import read_analyzed_sweep
import sys
import sweeptools as tls
from importlib import reload
reload(tls)

# Define parameters of screen to read
params = {'screen_name': 'c-MYC_HAP1',
          'assembly': 'hg38',
          'trim_length': '50',
          'mode': 'collapse',
          'start': 'tx',
          'end': 'tx',
          'overlap': 'both',
          'direction': 'sense',
          'step': 500}
# Flagging parameters
slope_thr = 1
mi_thr = 2
p_thr = 1e-5

sweep_data_dir = '../data/sweeps-analyzed_2020-09-21'

grouped_sweep = read_analyzed_sweep(sweep_data_dir, params)
# %%

out_path = (f'''flags/{params['screen_name']}/'''
            f'''{params['assembly']}/{params['trim_length']}/'''
            f'''mode={params['mode']}_direction={params['direction']}'''
            f'''_overlap={params['overlap']}/double-sweep'''
            f'''_step={params['step']}/''')

if not os.path.exists(out_path):
    os.makedirs(out_path)
    print('Creating analyzed directory.')

filename = (f'{out_path}flags_sl-thr={slope_thr}_p-thr={p_thr}'
            f'_mi-thr={mi_thr}.txt')

# Prevent overwriting
if os.path.exists(filename):
    raise SystemExit('File already exists, delete or rename to continue.')
    # print('File already exists, creating new file name.')
    # filename = ''.join([''.join([filename.split('.')[0], '_.']),
    #                     filename.split('.')[1]])

print(f'Finding flagged genes for screen {params["screen_name"]} - '
      f'{params["assembly"]} - sl-thr={slope_thr} - p_thr={p_thr} - '
      f'mi_thr={mi_thr}')

original_stdout = sys.stdout
with open(filename, 'w') as f:
    sys.stdout = f  # Change the standard output to the file we created.

    print(f'Flagged genes for screen {params["screen_name"]} - '
          f'{params["assembly"]} - sl-thr={slope_thr} - p_thr={p_thr} - '
          f'mi_thr={mi_thr}')
    flagged_genes = tls.analyzesweep.flag_by_slope(grouped_sweep, p_thr,
                                                   slope_thr, mi_thr)
    sys.stdout = original_stdout

print('Done.')


# %%
