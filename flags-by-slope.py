# %%
import os
from sweeptools.analyzesweep import read_analyzed_sweep
import sweeptools as tls
from importlib import reload
reload(tls)

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

# Flagging parameters
slope_thr = 1
p_thr = 1e-5

sweep_data_dir = '../data/sweeps/sweeps-analyzed_2020-09-21'
flag_data_dir = '../data/sweeps/sweep-flags_2020-10-02'

in_files = set(os.listdir(sweep_data_dir))
out_files = set(os.listdir(flag_data_dir))
missing_screens = list(in_files.difference(out_files))
missing_screens.sort()

screen_list = missing_screens[77+77:77+77+77]

# %%
# reload(tls)
for idx, screen in enumerate(screen_list):

    params['screen_name'] = screen
    print(f'\nFinding flagged genes for screen {idx+1} of {len(screen_list)}'
          f'\n{params["screen_name"]} - '
          f'{params["assembly"]} - sl_thr={slope_thr} - p_thr={p_thr}')
    grouped_sweep = read_analyzed_sweep(sweep_data_dir, params)
    flagged = tls.analyzesweep.flag_by_slope(grouped_sweep, p_thr, slope_thr,
                                             #  mi_seen_thr,p_seen_thr
                                             )

    out_path = (f'''{flag_data_dir}/{params['screen_name']}/'''
                f'''{params['assembly']}/{params['trim_length']}/'''
                f'''mode={params['mode']}_direction={params['direction']}'''
                f'''_overlap={params['overlap']}/double-sweep'''
                f'''_step={params['step']}/''')

    if not os.path.exists(out_path):
        os.makedirs(out_path)
        print('Creating analyzed directory.')

    filename = (f'{out_path}flags_sl-thr={slope_thr}_p-thr={p_thr}.csv')

    # Prevent overwriting
    if os.path.exists(filename):
        raise SystemExit('File already exists, delete or rename to continue.')

    flagged.to_csv(filename, index=False)

# %%
