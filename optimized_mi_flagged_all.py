# %%
from bokeh.plotting import figure
from bokeh.io import export_svgs
from bokeh.models import (BasicTickFormatter, HoverTool, RangeTool,
                          TapTool, ColumnDataSource, Range1d, LinearAxis,
                          Tool, SingleIntervalTicker, PanTool, ResetTool)
from bokeh.transform import linear_cmap
from matplotlib import cm, colors
import os
from sweeptools.plotting.optimized_mi import OptimizedPlot
import pandas as pd
from bokeh.plotting import output_file, show
from bokeh.layouts import gridplot, column, row
import numpy as np

from sweeptools.analyzesweep import read_analyzed_sweep, optimize_flagged_genes

import sweeptools as tls
from importlib import reload
reload(tls)

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

# Flagging parameters
slope_thr = 1
p_thr = 1e-5

sweep_data_dir = '../data/sweeps/sweeps-analyzed_2020-09-21'
flag_data_dir = '../data/sweeps/sweep-flags_2020-10-02'

screens = [x for x in os.listdir(flag_data_dir) if 'readme' not in x
           and 'DS_Store' not in x]

# screens = screens[0:5]
opt_list = []
for idx, screen in enumerate(screens):
    print(f'Optimizing {screen}: {idx + 1} of {len(screens)}')

    params['screen_name'] = screen

    # Get sweep and flag data
    grouped_sweep = read_analyzed_sweep(sweep_data_dir, params)
    flag_path = (f'''{flag_data_dir}/{params['screen_name']}/'''
                 f'''{params['assembly']}/{params['trim_length']}/'''
                 f'''mode={params['mode']}_direction={params['direction']}'''
                 f'''_overlap={params['overlap']}/double-sweep'''
                 f'''_step={params['step']}/''')
    flagged = pd.read_csv(f'{flag_path}flags_sl-thr={slope_thr}_'
                          f'p-thr={p_thr}.csv')

    try:
        optimized_mi = optimize_flagged_genes(flagged, grouped_sweep,
                                              delta_mi_thr=0)
        optimized_mi['screen'] = screen
        opt_list.append(optimized_mi)
    except ValueError:
        print('Empty flags for', screen)

optimized_all = pd.concat(opt_list)

# %% Sort and save excel file
opt_all_sorted = optimized_all.sort_values(by=['p_fdr_at_tx'],
                                           ascending=False)
opt_all_sorted = opt_all_sorted.rename(columns={'log2_mi': 'mi_opt',
                                                'p': 'p_opt',
                                                'low_counts': 'low_opt',
                                                'high_counts': 'high_opt',
                                                'srt_off': 'srt_off_opt',
                                                'end_off': 'end_off_opt'})

cols = ['gene_name', 'screen', 'mi_at_tx', 'p_at_tx', 'p_fdr_at_tx', 'mi_opt',
        'p_opt', 'srt_off_opt', 'end_off_opt']
opt_all_sorted = opt_all_sorted[cols]
opt_all_sorted.to_excel('optimized_all.xlsx', index=False)

# %% Filter by original pfdr
p_fdr_thr = 0.05
optimized_all = optimized_all.query('p_fdr_at_tx > @p_fdr_thr')

# %%
plot_filename = 'plots/test_optimized_all.html'


class DeltaOptimizedPlot():

    def __init__(self, assembly: str,
                 optimized_all: pd.DataFrame) -> None:

        self.optimized = optimized_all
        self.screen = screen
        self.assembly = assembly

        self.set_colors()
        self.set_source()

        size = 800
        max_x = max(abs(self.optimized.delta_log2_mi)) + 0.4
        max_y = max(abs(self.optimized.delta_log10_p)) + 4

        self.plt = figure(plot_width=size, plot_height=int(np.ceil(size/1.5)),
                          x_range=Range1d(-max_x, max_x, bounds='auto'),
                          y_range=Range1d(-max_y, max_y, bounds='auto'),
                          margin=(10, 0, 0, 0),
                          tools=('wheel_zoom, pan, reset, box_zoom, tap, '
                                 'lasso_select'),
                          active_scroll='wheel_zoom'
                          )
        self.plt.toolbar.logo = None
        ticker = SingleIntervalTicker(interval=1, num_minor_ticks=5)
        self.plt.xaxis.ticker = ticker
        self.plt.xaxis.axis_label = '\u0394 log2 MI'
        self.plt.yaxis.axis_label = '\u0394 log10 p'

        self.plt.circle(x='delta_log2_mi', y='delta_log10_p', source=self.src,
                        color=self.opt_color, size='opt_p_size',
                        name='to_inspect',
                        alpha=0.8, line_alpha=0.8,
                        line_color='white', line_width=0.5,
                        selection_line_color='#4B98E5',
                        selection_line_width=6,
                        selection_line_alpha=1,
                        nonselection_alpha=0.8
                        )

        self.plt.add_tools(HoverTool(
            tooltips=[('Screen', '@screen'),
                      ('Gene', '@gene_name'),
                      #   ('Delta MI', '@delta_log2_mi'),
                      #   ('Delta p', '@delta_log10_p'),
                      ('Opt log2 MI', '@log2_mi'),
                      #   ('Tx log2 MI', '@mi_at_tx'),
                      ('Opt p-value', '@p'),
                      #   ('Tx p-value', '@p_at_tx'),
                      ('Opt start offset', '@srt_off'),
                      ('Opt end offset', '@end_off'),
                      #   ('Opt high counts', '@high_counts'),
                      #   ('Opt low counts', '@low_counts'),
                      ],
            names=['to_inspect'],
            mode='mouse', point_policy='snap_to_data'))

    def set_source(self) -> None:
        # Set sizes according to log10 of p-value
        max_val = 30  # Max -log10(p) value that will still increase size
        min_size = 5
        max_size = 15

        self.optimized['log_p_tx'] = [-np.log10(x) if x > 0 else 50
                                      for x in self.optimized.p_at_tx]
        self.optimized['log_opt_p'] = [-np.log10(x) if x > 0 else 50
                                       for x in self.optimized.p]
        self.optimized['opt_p_size'] = np.interp(self.optimized.log_opt_p,
                                                 (0, max_val),
                                                 (min_size, max_size))
        self.optimized['delta_log10_p'] = (self.optimized['log_opt_p']
                                           - self.optimized['log_p_tx'])
        self.optimized['delta_log2_mi'] = (self.optimized.log2_mi
                                           - self.optimized.mi_at_tx)

        self.src = ColumnDataSource(self.optimized)

    def set_colors(self) -> None:
        # Set colors according to log2MI value
        cmap = cm.get_cmap('PiYG', 256)
        PiYG256 = tuple(colors.rgb2hex(cmap(i)[:3]) for i in range(cmap.N))
        palette = PiYG256
        self.opt_color = linear_cmap(field_name='log2_mi', palette=palette,
                                     low=-8, high=8)

    def for_range_selection(self, other_plot: figure) -> None:
        # (a, b) = other_plot.x_range.bounds
        self.plt.x_range = Range1d(
            deltas.plt.x_range.start, deltas.plt.x_range.end, bounds='auto')
        range_tool = RangeTool(x_range=other_plot.x_range,
                               y_range=other_plot.y_range)
        range_tool.overlay.fill_color = 'navy'
        range_tool.overlay.fill_alpha = 0.1
        self.plt.add_tools(range_tool)
        self.plt.toolbar.active_inspect = None
        self.plt.toolbar_location = None
        self.plt.plot_height = 150
        self.plt.plot_width = 250
        self.plt.xgrid.grid_line_color = None
        self.plt.ygrid.grid_line_color = None
        self.plt.xaxis.visible = False
        self.plt.yaxis.visible = False

        self.plt.select(name='to_inspect').glyph.size = 3


class FishtailPlot():

    def __init__(self, src) -> None:

        size = 550

        max_x = max(max(src.data['counts']),
                    max(src.data['counts_at_tx'])) + 1000
        max_y = max(max(abs(src.data['log2_mi'])),
                    max(abs(src.data['mi_at_tx']))) + 0.4

        self.plt = figure(x_axis_label='Insertions', y_axis_label='log2 MI',
                          plot_width=size, plot_height=int(np.ceil(size/1.7)),
                          x_axis_type='log',
                          margin=(30, 0, 0, 50),
                          x_range=Range1d(0.9, max_x, bounds='auto'),
                          y_range=Range1d(-max_y, max_y, bounds='auto'),
                          tools=('wheel_zoom, pan, reset, box_zoom, tap'),
                          active_scroll='wheel_zoom')
        self.plt.toolbar.logo = None

        self.plt.xaxis[0].ticker.base = 10
        self.plt.xaxis.formatter = BasicTickFormatter(use_scientific=False)

        self.set_colors()

        self.plt.circle(x='counts', y='log2_mi', source=src,
                        color=self.opt_color, size=8, line_color='white',
                        alpha=0.1, selection_alpha=1, nonselection_alpha=0.03,
                        selection_line_color='#4B98E5',
                        selection_line_width=2, line_alpha=0.1,
                        line_width=0.5)
        self.plt.square(x='counts_at_tx', y='mi_at_tx', source=src,
                        color=self.tx_color, size=8, line_color='white',
                        alpha=0, selection_alpha=1, nonselection_alpha=0,
                        selection_line_color='#4B98E5',
                        selection_line_width=2, line_alpha=0.1,
                        line_width=0.5)

        self.plt.add_tools(HoverTool(
            tooltips=[('Screen', '@screen'),
                      ('Gene', '@gene_name'),
                      #   ('Delta MI', '@delta_log2_mi'),
                      #   ('Delta p', '@delta_log10_p'),
                      ('Opt log2 MI', '@log2_mi'),
                      #   ('Tx log2 MI', '@mi_at_tx'),
                      ('Opt p-value', '@p'),
                      #   ('Tx p-value', '@p_at_tx'),
                      ('Opt start offset', '@srt_off'),
                      ('Opt end offset', '@end_off'),
                      #   ('Opt high counts', '@high_counts'),
                      #   ('Opt low counts', '@low_counts'),
                      ],
            # names=['to_inspect'],
            mode='mouse', point_policy='snap_to_data'))

    def set_colors(self) -> None:
        # Set colors according to log2MI value
        cmap = cm.get_cmap('PiYG', 256)
        PiYG256 = tuple(colors.rgb2hex(cmap(i)[:3]) for i in range(cmap.N))
        palette = PiYG256
        self.opt_color = linear_cmap(field_name='log2_mi', palette=palette,
                                     low=-8, high=8)
        self.tx_color = linear_cmap(field_name='mi_at_tx', palette=palette,
                                    low=-8, high=8)


class VolcanoPlot():

    def __init__(self, src) -> None:

        size = 550
        max_y = max(max(src.data['log_opt_p']),
                    max(src.data['log_p_tx'])) + 1.5
        max_x = max(max(abs(src.data['log2_mi'])),
                    max(abs(src.data['mi_at_tx']))) + 0.4

        self.plt = figure(x_axis_label='log2 MI', y_axis_label='-log10 p',
                          plot_width=size, plot_height=int(np.ceil(size/1.7)),
                          margin=(30, 0, 0, 50),
                          x_range=Range1d(-max_x, max_x, bounds='auto'),
                          y_range=Range1d(0, max_y, bounds='auto'),
                          tools=('wheel_zoom, pan, reset, box_zoom, tap'),
                          active_scroll='wheel_zoom')

        self.set_colors()
        self.plt.toolbar.logo = None

        self.plt.circle(x='log2_mi', y='log_opt_p', source=src,
                        color=self.opt_color, size=8, line_color='white',
                        alpha=.1, selection_alpha=1, nonselection_alpha=0.03,
                        selection_line_color='#4B98E5',
                        selection_line_width=2, line_alpha=0.1,
                        line_width=0.5)
        self.plt.square(x='mi_at_tx', y='log_p_tx', source=src,
                        color=self.tx_color, size=8, line_color='white',
                        alpha=0, selection_alpha=1, nonselection_alpha=0,
                        selection_line_color='#4B98E5',
                        selection_line_width=2, line_alpha=0.1,
                        line_width=0.5)

    def set_colors(self) -> None:
        # Set colors according to log2MI value
        cmap = cm.get_cmap('PiYG', 256)
        PiYG256 = tuple(colors.rgb2hex(cmap(i)[:3]) for i in range(cmap.N))
        palette = PiYG256
        self.opt_color = linear_cmap(field_name='log2_mi', palette=palette,
                                     low=-8, high=8)
        self.tx_color = linear_cmap(field_name='mi_at_tx', palette=palette,
                                    low=-8, high=8)


deltas = DeltaOptimizedPlot(params['assembly'], optimized_all)
range_select = DeltaOptimizedPlot(params['assembly'], optimized_all)
range_select.for_range_selection(deltas.plt)

fish = FishtailPlot(deltas.src)
volc = VolcanoPlot(deltas.src)

output_file(plot_filename)
layout = row(column(range_select.plt, deltas.plt), column(volc.plt, fish.plt))
show(layout)

deltas.plt.output_backend = 'svg'
export_svgs(deltas.plt, filename='plots/delta-optimized.svg')
# %%
