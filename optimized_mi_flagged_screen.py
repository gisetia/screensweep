# %%
import os
import pandas as pd
from bokeh.plotting import output_file, show
from bokeh.layouts import gridplot

from sweeptools.analyzesweep import read_analyzed_sweep, optimize_flagged_genes

import sweeptools as tls
from importlib import reload
reload(tls)

# Define parameters of screen to read
params = {
    'screen_name': 'p-RPA',
    # 'screen_name': 'Ytub_in_TTL-VASH1-2-TKO',
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

screens = [x for x in os.listdir(flag_data_dir) if 'readme' not in x
           and 'DS_Store' not in x]

# %%
for screen in screens:

    params['screen_name'] = screen
    print(f'Getting optimized MI for screen {screen}.')
    plot_filename = f'plots/2020-10-13/{params["screen_name"]}.html'

    # Get sweep and flag data
    grouped_sweep = read_analyzed_sweep(sweep_data_dir, params)
    flag_path = (f'''{flag_data_dir}/{params['screen_name']}/'''
                 f'''{params['assembly']}/{params['trim_length']}/'''
                 f'''mode={params['mode']}_direction={params['direction']}'''
                 f'''_overlap={params['overlap']}/double-sweep'''
                 f'''_step={params['step']}/''')
    flagged = pd.read_csv(f'{flag_path}flags_sl-thr={slope_thr}_'
                          f'p-thr={p_thr}.csv')

    def background_genes(grouped_sweep):
        # Get background points (all values at tx start and end)
        all_list = []
        for gene, group in grouped_sweep:
            all_at_tx = pd.Series()
            try:
                all_at_tx['high_counts'] = group.loc[0, 0]['high_counts']
                all_at_tx['low_counts'] = group.loc[0, 0]['low_counts']
                all_at_tx['log2_mi'] = group.loc[0, 0]['log2_mi']
                all_at_tx['p'] = group.loc[0, 0]['p']
            except KeyError:
                pass
            all_list.append(all_at_tx)
        background = pd.concat(all_list, axis=1).T

        return background

    background = background_genes(grouped_sweep)
    try:
        optimized_mi = optimize_flagged_genes(flagged, grouped_sweep,
                                              delta_mi_thr=0.5)
    except ValueError:
        print(f'No flags for {params["screen_name"]}')
        continue

    opt_mi = tls.plotting.optimized_mi.OptimizedPlot(params['screen_name'],
                                                     params['assembly'],
                                                     optimized_mi, mode='mi')

    select = tls.plotting.optimized_mi.OptimizedPlot(params['screen_name'],
                                                     params['assembly'],
                                                     optimized_mi, mode='mi',
                                                     src=opt_mi.src)
    select.for_range_selection(opt_mi.plt)
    select.add_title()

    opt_p = tls.plotting.optimized_mi.OptimizedPlot(params['screen_name'],
                                                    params['assembly'],
                                                    optimized_mi, mode='p',
                                                    src=opt_mi.src,
                                                    x_range=opt_mi.plt.x_range,
                                                    tools=opt_mi.tools)
    opt_p.set_x_axis(opt_mi.src)

    volc = tls.plotting.optimized_mi.OptimizedPlot(params['screen_name'],
                                                   params['assembly'],
                                                   optimized_mi,
                                                   mode='volcano',
                                                   src=opt_mi.src)
    volc.add_background(background)
    fish = tls.plotting.optimized_mi.OptimizedPlot(params['screen_name'],
                                                   params['assembly'],
                                                   optimized_mi,
                                                   mode='fishtail',
                                                   src=opt_mi.src)
    fish.add_background(background)

    output_file(plot_filename)

    layout = gridplot([[select.plt, None], [opt_mi.plt, volc.plt],
                       [opt_p.plt, fish.plt]], merge_tools=False)

    show(layout)

# %%
