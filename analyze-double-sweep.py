# %%

import pandas as pd

from analyze import read_analyzed_sweep, get_flagged_genes

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

data_dir = 'data/analyzed-data'
# data_dir = 'sample-data/analyzed-data'

grouped_sweep = read_analyzed_sweep(data_dir, params)

# %%

slope_thr = 3
p_ratio_thr= 4
flagged_genes = get_flagged_genes(grouped_sweep, slope_thr, p_ratio_thr)

# %%

# Plot

gene = 'ZBTB14'

gene_info = pd.DataFrame()
for name, group in grouped_sweep:
    if name == gene.upper():
        
        gene_info = group.unstack()
        
from math import ceil
from numpy import interp, where
from bokeh.plotting import figure, output_file, show
from bokeh.models import ColumnDataSource, ColorBar
from bokeh.transform import linear_cmap
from bokeh.palettes import Spectral6, Viridis8, GnBu8, PiYG8, RdYlBu8


palette = PiYG8[::-1]
plot_width = 600

# Plot

# Axes ranges
padd = ceil(params['step']/1.5)
xlim = (-10000 - padd, 2000 + padd)
ylim = (10000 + padd, -2000 - padd)

plt = figure(title=f'Gene: {gene.upper()} - Screen: {params["screen_name"]}',
               x_axis_label='End offset (bp)',
               y_axis_label='Start offset (bp)',
               plot_width=plot_width,
               plot_height=plot_width-15,
               match_aspect=True,
               aspect_scale=1,
               x_range=xlim,
               y_range=ylim
               )

# Rename 0 to position (tx, cds end,start)
plt.xaxis.major_label_overrides = {0: (f'End '
                                     f'{gene_info["end_pos"].iloc[0,0]}')}
plt.yaxis.major_label_overrides = {0: (f'Start '
                                     f'{gene_info["srt_pos"].iloc[0,0]}')}
# Lines to mark 0,0
plt.line(x=[0, 0], y=[ylim[0], ylim[1]], color='#5f5f61')
plt.line(x=[xlim[0], xlim[1]], y=[0, 0], color='#5f5f61')

s = gene_info['log2_mi'].stack().reset_index()
s = s.rename(columns={0: 'log2_mi'})

p_thr = 0.00001
p = gene_info['p_fdr'].stack().reset_index()
p = p.rename(columns={0: 'p_fdr'})
p['p_masked'] = where(p.p_fdr < p_thr, 1, 0)
s['log2_mi_masked'] = s.log2_mi.where(p.p_masked == 1)

# Rescale log2 MI to get nice sizes in plot
min_size = plot_width/100
max_size = plot_width/38
s['log2_mi_rescaled'] = s['log2_mi'].map(lambda x: interp(abs(x), (0, 8), 
                                                         (min_size, 
                                                          max_size)))

source = ColumnDataSource(s)

line_color = linear_cmap(field_name='log2_mi', palette=palette,
                         low=-8 ,high=8)
fill_color = linear_cmap(field_name='log2_mi_masked', palette=palette,
                         low=-8 ,high=8, nan_color='white')

plt.circle(x='end_off', y='srt_off', source=source, size='log2_mi_rescaled',
line_color=line_color, color=fill_color, line_width=3)

color_bar = ColorBar(color_mapper=line_color['transform'], width=15, 
                     location=(0,0))
plt.add_layout(color_bar, 'right')

output_file('test_plot.html')
show(plt)

# %%
