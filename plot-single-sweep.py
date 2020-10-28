# %%
# import pandas as pd
from bokeh.io import export_svgs
from math import pi, ceil
from numpy import interp
from scipy import stats
from matplotlib import cm, colors
from math import log10

from bokeh.plotting import figure, output_file, show, curdoc

from bokeh.plotting import figure
from bokeh.models import (ColumnDataSource, Range1d, LinearAxis, HoverTool, Text,
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
params = {'screen_name': 'PDL1_IFNg',
          'assembly': 'hg38',
          'trim_length': '50',
          'mode': 'collapse',
          'start': 'tx',
          'end': 'tx',
          'overlap': 'both',
          'direction': 'sense',
          'step': 500}

sweep_data_dir = '../data/sweeps/sweeps-analyzed_2020-09-21'

grouped_sweep = read_analyzed_sweep(sweep_data_dir, params)

gene = 'CD274'
group = grouped_sweep.get_group(gene)

single = group.query('srt_off == 0').reset_index()
single['slope'] = (single['log2_mi']
                   - single['log2_mi'].shift(1))*1000/params['step']

# %%
cmap = cm.get_cmap('PiYG', 256)
PiYG256 = tuple(colors.rgb2hex(cmap(i)[:3]) for i in range(cmap.N))
palette = PiYG256
color = linear_cmap(field_name='log2_mi', palette=palette,
                    low=-8, high=8)

src = ColumnDataSource(single)
plt = figure(title=f'Gene: {gene} - Screen: {params["screen_name"]}'
             f' - Assembly: {params["assembly"]}',
             x_axis_label='End offset (bp)',
             y_axis_label='log2 MI',
             plot_width=600,
             plot_height=300,
             x_range=(-10000, 0),
             y_range=(-8, 8),
             toolbar_location=None,
             #  margin=(0, 0, 0, 0)
             )
plt.xgrid.visible = False
plt.ygrid.visible = False
# plt.yaxis[0].axis_label_text_color = '#1F77B4'
# plt.yaxis[0].major_label_text_color = '#1F77B4'


plt.line(x=[plt.x_range.start, plt.x_range.end],
         y=[0, 0], color='#B8B8BB')

plt.line(x='end_off', y='log2_mi', source=src, color='#E2E2E4')
plt.circle(x='end_off', y='log2_mi', source=src, size=10, color=color)

plt.extra_y_ranges = {'slope': Range1d(start=-5, end=5)}
plt.add_layout(LinearAxis(y_range_name='slope',
                          axis_label='\u0394 log2 MI / 1 kb',
                          # axis_label_text_font_size=label_font_size,
                          # major_label_text_font_size=label_font_size),
                          ), 'right')
plt.yaxis[1].axis_label_text_color = '#11A0CE'
plt.yaxis[1].major_label_text_color = '#11A0CE'
plt.line(x='end_off', y='slope', source=src, color='#11A0CE',
         y_range_name='slope')
# plt.circle(x='end_off', y='slope', source=src, color='#11A0CE', size=10,
#            y_range_name='slope')

output_file('plots/single-sweep.html')
show(plt)

plt.output_backend = 'svg'
export_svgs(plt, filename='plots/single-sweep.svg')

# %%
