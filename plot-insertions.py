# %%
# from bokeh.layouts import row, column
from bokeh.models import NumeralTickFormatter, Range1d, LinearAxis
from bokeh.transform import linear_cmap
from bokeh.palettes import PiYG8
from bokeh.models import (ColumnDataSource, Arrow, OpenHead, HoverTool, Text,
                          CustomJS, TapTool, SquarePin, CrosshairTool,
                          Rect, RangeTool, Text)
from bokeh.plotting import figure, output_file, show, curdoc
from bokeh.layouts import column, row
from bokeh.transform import linear_cmap

from bokeh.palettes import PiYG8
from matplotlib import cm, colors

from typing import List, Union, Optional

from math import pi
import pandas as pd
from math import log2

# from tools.plots import plot_insertions, ins_select_range

import tools as tls
from importlib import reload
reload(tls)

data_dir = 'data/screen-analyzer-data'
chrom = '9'

# super enhancer
# start = 5495000
# end = 5505000


# start = 5580709
# end = 5591016

# https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0166626
# start = 5455433
# end = 5455446

# # # pdl1
# start = 5450503
# end = 5470567

#jak2
start = 4985033
end = 5128183

# start = 5492001
# end = 5495000

padd = 50

screen_name = 'PDL1_IFNg'
assembly = 'hg19'
trim_length = 50


ins = tls.plots.plot_insertions(data_dir, screen_name, assembly,
                                trim_length, chrom, start, end, 500000)
select = tls.plots.ins_select_range(ins)


def get_ratios(insertions):
    counts_per_strand = insertions.groupby(['chan', 'strand']).size()

    # Get counts for each channel and strand, set to 1 if no counts
    cnts_per_strand = {'+h': counts_per_strand.get(('high', '+'), 1),
                       '-h': counts_per_strand.get(('high', '-'), 1),
                       '+l': counts_per_strand.get(('low', '+'), 1),
                       '-l': counts_per_strand.get(('low', '-'), 1)}

    # print(cnts_per_strand)

    counts_both_strands = insertions.groupby(['chan']).size()
    cnts_both_strands = {'h': counts_both_strands.get(('high'), 1),
                         'l': counts_both_strands.get(('low'), 1)}
    # print(cnts_both_strands)

    ratios = {'+': log2(cnts_per_strand['+h']/cnts_per_strand['+l']),
              '-': log2(cnts_per_strand['-h']/cnts_per_strand['-l']),
              'both': log2(cnts_both_strands['h']/cnts_both_strands['l'])}

    # print(ratios)
    return ratios


def update_ratios(attr, old, new):

    # Get selected end and start
    start_sel = int(range_tool.x_range.start)
    end_sel = int(range_tool.x_range.end)

    # Filter insertions between new start and end
    ins_sel = insertions[insertions.pos.between(start_sel, end_sel)]

    ratios = get_ratios(ins_sel)

    rat_source.data['log2rat'] = [ratios['+'], ratios['-'], ratios['both']]
    rat_source.data['log2rat_str'] = list(map(lambda x: f'{x:.1f}',
                                              rat_source.data['log2rat']))

    ins.title.text = (f'Insertions of screen {screen_name} in '
                        f'chr{chrom}:{start_sel:,} - {end_sel:,}')

    return ratios


range_tool = select.tools[0]
ins_source = [x.data_source for x in ins.renderers
              if x.name == 'insertions'][0]
insertions = pd.DataFrame(ins_source.data)[['chan', 'chr', 'strand', 'pos']]


range_tool.x_range.on_change('start', update_ratios)
range_tool.x_range.on_change('end', update_ratios)


# Initialize ratios figure
rat = figure(title='Log2(high/low)', plot_width=200, plot_height=500,
             toolbar_location=None, y_range=(0, 4), x_range=(0, 10), 
             tools='')
rat.xaxis.visible = False
rat.yaxis.visible = False
rat.xgrid.visible = False
rat.ygrid.visible = False
rat.outline_line_color = None

ins_sel = insertions[insertions.pos.between(ins.x_range.start,
                                            ins.x_range.end)]

ratios = get_ratios(ins_sel)

rat_df = pd.DataFrame.from_dict(ratios, orient='index', columns=['log2rat'])
rat_df['log2rat_str'] = rat_df.log2rat.apply(lambda x: f'{x:.1f}')
rat_df['ypos'] = [3, 2, 1]
rat_df['xpos'] = [4, 4, 4]
rat_df['p_val'] = [0, 0, 0]

rat_source = ColumnDataSource(rat_df)

# Set colors
cmap = cm.get_cmap('PiYG', 256)
PiYG256 = tuple(colors.rgb2hex(cmap(i)[:3]) for i in range(cmap.N))
palette = PiYG256
line_color = linear_cmap(field_name='log2rat', palette=palette,
                         low=-8, high=8)
fill_color = 'white'

# Plot
rat.circle(x='xpos', y='ypos', size=50, line_color=line_color, line_width=10,
           fill_color=fill_color, source=rat_source)


# Add log2(ratio) text
rat_text = Text(x='xpos', y='ypos', text='log2rat_str',
                x_offset=40, y_offset=10, text_font_size='14pt')
rat.add_glyph(rat_source, rat_text)


# Add legend text
leg_source = ColumnDataSource(dict(x=[5, 5, 5], y=[3.3, 2.3, 1.3],
                                   text=['+ strand', '- strand',
                                         'Both strands']))
leg_glyph = Text(x='x', y='y', text='text', text_font_size='10pt',
                 text_align='center')
rat.add_glyph(leg_source, leg_glyph)


layout = row(column(ins, select), rat)

curdoc().title = 'Insertions per region'
curdoc().add_root(layout)

output_file('plots/test_insertions.html')
show(layout)
# %%
