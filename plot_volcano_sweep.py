# %%import pandas as pd
from bokeh.io import export_svgs
from math import pi, ceil
from numpy import interp
from scipy import stats
from matplotlib import cm, colors
from math import log10

from bokeh.plotting import figure, output_file, show, curdoc

from bokeh.plotting import figure
from bokeh.models import (ColumnDataSource, Arrow, OpenHead, HoverTool, Text,
                          CustomJS, TapTool, SquarePin)
from bokeh.palettes import PiYG8
from bokeh.transform import linear_cmap
from bokeh.layouts import row, column
from bokeh.models.annotations import Title

from sweeptools.analyzesweep import read_analyzed_sweep, flags_query
import sweeptools as tls
from importlib import reload
reload(tls)

# Define parameters of screen to read
params = {'screen_name': 'Ac-beta-actin_WT',
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
flag_data_dir = 'test_flags'

grouped_sweep = read_analyzed_sweep(sweep_data_dir, params)

gene = 'CRAMP1'
group = grouped_sweep.get_group(gene)
group['log10_p'] = group.p.apply(lambda x: abs(-log10(x)) if x > 0 else 50)
group = group.reset_index()

# %

cmap = cm.get_cmap('PiYG', 256)
PiYG256 = tuple(colors.rgb2hex(cmap(i)[:3]) for i in range(cmap.N))
palette = PiYG256
color = linear_cmap(field_name='log2_mi', palette=palette,
                    low=-8, high=8)

src = ColumnDataSource(group)
plt = figure(title=f'Gene: {gene} - Screen: {params["screen_name"]}'
             f' - Assembly: {params["assembly"]}',
             x_axis_label='Log2 MI',
             y_axis_label='-log10 p',
             plot_width=600,
             plot_height=300,
             x_range=(-10, 10),
             y_range=(0, 55),
             toolbar_location=None,
             #  margin=(0, 0, 0, 0)
             )

plt.line(x=[plt.x_range.start, plt.x_range.end],
         y=[-log10(0.05), -log10(0.05)],
         dash=[3, 2], color='#B8B8BB')

plt.circle(x='log2_mi', y='log10_p', color=color, size=10,
           line_color='white', alpha=0.8, line_alpha=0.3, source=src)
plt.add_tools(HoverTool(tooltips=[('Log2 MI', '@log2_mi'),
                                  ('Start offset', '@srt_off'),
                                  ('End offset', '@end_off'),
                                  ('High counts', '@high_counts'),
                                  ('Low counts', '@low_counts'),
                                  ('p-value', '@p'),
                                  ('fdr p', '@p_fdr'), ]))

output_file('plots/volcano_test4.html')
show(plt)

plt.output_backend = 'svg'
export_svgs(plt, filename='plots/volcano.svg')
# %%
