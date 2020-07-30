# %%

from tools.analyzeinsertions import (read_all_insertions, get_gene_insertions, 
                                    read_gene_insertions, get_gene_positions)

# from importlib import reload
# reload(tools)

params = {'screen_name': 'PDL1_IFNg',
        'assembly': 'hg38',
        'trim_length': '50',
        'mode': 'collapse',
        'start': 'tx',
        'end': 'tx',
        'overlap': 'both',
        'direction': 'sense',
        'step': 500}

data_dir = 'data/screen-analyzer-data'

gene = 'SOCS1'
padd = 2000
gene_pos = get_gene_positions(gene, params['assembly'])
insertions = read_gene_insertions(gene, data_dir, params, 
                                  gene_pos=gene_pos, padding=padd)

# %%
from math import pi
from bokeh.plotting import figure, output_file, show
from bokeh.models import ColumnDataSource, HoverTool, Text
from bokeh.models.annotations import Title
from bokeh.transform import linear_cmap
from bokeh.palettes import Spectral11, Viridis8, GnBu8, PiYG8, RdYlBu8
from bokeh.layouts import row, column

# Plot setup --------------------------------------------------------------
#region -------------------------------------------------------------------

zero = min(gene_pos['txStart'])
xlim = (min(gene_pos['txStart']) - padd - zero,
        max(gene_pos['txEnd']) + padd - zero)
ylim = (0, 8 + len(gene_pos))

plt = figure(title=f'Insertions in gene: {gene.upper()} - '
                   f'Screen: {params["screen_name"]}',
               x_axis_label='Position (bp)',
               plot_width=1000,
               plot_height=300,

               x_range=(xlim[0]-padd*0.1, xlim[1] + padd*0.1),
               y_range=ylim
               )

plt.ygrid.grid_line_color = None
plt.xgrid.grid_line_color = None
plt.yaxis.major_tick_line_color = None
plt.yaxis.minor_tick_line_color = None
plt.yaxis.axis_line_color = None
plt.outline_line_color = None
plt.xaxis.major_label_overrides = {0: 'tx Start'}
plt.yaxis.ticker = [1, 2, 4, 5, 7]
plt.yaxis.major_label_overrides = {1: 'Low anti-sense',
                                   2: 'High anti-sense',
                                   4: 'Low sense',
                                   5: 'High sense',
                                   7: 'Transcript(s)'}

plt.line(x=xlim, y=[1, 1], color='#BCBCBF', line_width=2)
plt.line(x=xlim, y=[2, 2], color='#BCBCBF', line_width=2)
plt.line(x=xlim, y=[4, 4], color='#BCBCBF', line_width=2)
plt.line(x=xlim, y=[5, 5], color='#BCBCBF', line_width=2)

#endregion

# Plot transcripts --------------------------------------------------------
#region -------------------------------------------------------------------

# Plot dashed lines where transcript(s) end(s)
plt.line(x=(xlim[0]+padd, xlim[0]+padd), y=ylim, color='#E2E2E4', 
         line_width=1, line_dash=[5,2])
plt.line(x=(xlim[1]-padd, xlim[1]-padd), y=ylim, color='#E2E2E4',
         line_width=1, line_dash=[5,2])

from bokeh.models import Arrow, OpenHead, Rect
import pandas as pd
from numpy import NaN

# Plot transcripts
gene_pos = gene_pos.reset_index(drop=True)
for idx, tx in gene_pos.iterrows():

    # Plot line and arrow through transcript
    tx_length = tx['txEnd'] - tx['txStart']

    # Position of arrow
    if tx['strand'] == '+':
        x_end = tx['txEnd'] - zero + tx_length/30
        x_srt = x_end - tx_length/30
    else:
        x_end = tx['txStart'] - zero - tx_length/30
        x_srt = x_end + tx_length/30
    ypos = idx + 7

    plt.add_layout(Arrow(end=OpenHead(size=6, line_width=2, 
                         line_color='#5f5f61'), 
                         line_color='#5f5f61', line_width=2,
                         x_start=x_srt, y_start=ypos, 
                         x_end=x_end, y_end=ypos))
    plt.line(x=(tx['txStart']-zero, tx['txEnd']-zero), y=(ypos, ypos), 
            color='#5f5f61', line_width=2)






    cds_range = range(tx['cdsStart'], tx['cdsEnd'] + 1)
    



    

    exons = pd.DataFrame(columns=['start', 'end'])
    exons['start'] = tx['exonStarts'].split(',')
    exons['end'] = tx['exonEnds'].split(',')
    exons = exons.replace('', NaN).dropna().astype(int)

    exons['width'] = exons.eval('end - start')
    exons['height'] = 0.8
    exons['ypos'] = ypos
    exons['start_plot'] = exons.eval('start - @zero + width/2')

    ex_source = ColumnDataSource(exons)
    glyph = Rect(x='start_plot', y='ypos', width='width', height='height', 
                fill_color='#fdae61', line_color=None)
    plt.add_glyph(ex_source, glyph)
    
#endregion

# Plot insertions ---------------------------------------------------------
#region -------------------------------------------------------------------

# Transform data source for insertions
def get_ins_color(row):
    ins_colors = {  'hs': PiYG8[0],
                    'ls': PiYG8[-1],
                    'ha': PiYG8[3],
                    'la': PiYG8[-3]}
    return ins_colors[row['chan'][0]+row['dir'][0]]

insertions['color'] = insertions.apply(lambda row: get_ins_color(row),
                                       axis=1)
def get_ins_ypos(row):
    ins_pos = { 'hs': 5,
                'ls': 4,
                'ha': 2,
                'la': 1}
    return ins_pos[row['chan'][0]+row['dir'][0]]

insertions['ypos'] = insertions.apply(lambda row: get_ins_ypos(row),
                                       axis=1)
insertions['xpos'] = insertions.apply(lambda row: row['pos']-zero,
                                       axis=1)
source = ColumnDataSource(insertions)

# Plot insertions
plt.dash(x='xpos', y='ypos', color='color', source=source,
         angle=pi/2, line_width = 2, size=15,)

#endregion

output_file('plots/test_insertions.html')
# legend = column(col_leg, p_leg, flag_leg)
# layout = row(plt, legend)
show(plt)

# %%
