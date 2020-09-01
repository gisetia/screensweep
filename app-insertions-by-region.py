# %%
import pandas as pd
from bokeh.transform import linear_cmap
from bokeh.models import (ColumnDataSource, Text)
from bokeh.plotting import figure, output_file, show, curdoc
from bokeh.layouts import column, row
from matplotlib import cm, colors
from math import log2
from scipy.stats import fisher_exact

from tools.plotting.insertionsrange import plot_insertions, ins_select_range

# import tools as tls
# from importlib import reload
# reload(tls)

data_dir = 'data/screen-analyzer-data'
screen_name = 'PDL1_IFNg'
assembly = 'hg38'
trim_length = 50

filename = 'plots/test_insertions1.html'

# position = 'chr9:4,984,390-5,129,948'

position = 'chr6:127,027,001-127,027,500'

chrom = position.split(':')[0][3:]
start = int(position.split(':')[1].split('-')[0].replace(',',''))
end = int(position.split(':')[1].split('-')[1].replace(',',''))



# pdl1 super enhancer
# chrom = '9'
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

# jak2
# start = 4985033
# end = 5128183

# start = 5492001
# end = 5495000

ins = plot_insertions(data_dir, screen_name, assembly,
                      trim_length, chrom, start, end, 200000)
select = ins_select_range(ins)


def get_ratios(insertions, start, end):

    ins_sel = insertions[insertions.pos.between(start, end)]
    counts_per_strand = ins_sel.groupby(['chan', 'strand']).size()

    # Get counts for each channel and strand, set to 1 if no counts
    cnts_per_strand = {'+h': counts_per_strand.get(('high', '+'), 1),
                       '-h': counts_per_strand.get(('high', '-'), 1),
                       '+l': counts_per_strand.get(('low', '+'), 1),
                       '-l': counts_per_strand.get(('low', '-'), 1)}

    # print(cnts_per_strand)

    counts_both_strands = ins_sel.groupby(['chan']).size()
    cnts_both_strands = {'h': counts_both_strands.get(('high'), 1),
                         'l': counts_both_strands.get(('low'), 1)}
    # print(cnts_both_strands)

    ratios = {'+': log2(cnts_per_strand['+h']/cnts_per_strand['+l']),
              '-': log2(cnts_per_strand['-h']/cnts_per_strand['-l']),
              'both': log2(cnts_both_strands['h']/cnts_both_strands['l'])}

    rat_df = pd.DataFrame.from_dict(ratios, orient='index',
                                    columns=['log2rat'])

    tot_space = end - start
    _, p_plus = fisher_exact([[cnts_per_strand['+h'], tot_space],
                              [cnts_per_strand['+l'], tot_space]],
                             alternative='two-sided')
    _, p_minus = fisher_exact([[cnts_per_strand['-h'], tot_space],
                               [cnts_per_strand['-l'], tot_space]],
                              alternative='two-sided')
    _, p_both = fisher_exact([[cnts_both_strands['h'], tot_space],
                              [cnts_both_strands['l'], tot_space]],
                             alternative='two-sided')

    rat_df['p'] = [p_plus, p_minus, p_both]
    rat_df['log2rat_masked'] = rat_df.log2rat.where(rat_df.p < 0.00001)

    return rat_df


def update_ratios(attr, old, new):

    # Get selected end and start
    start_sel = int(range_tool.x_range.start)
    end_sel = int(range_tool.x_range.end)

    ratios = get_ratios(insertions, start_sel, end_sel)

    rat_source.data['log2rat'] = ratios['log2rat']
    rat_source.data['p'] = ratios['p']
    rat_source.data['log2rat_masked'] = ratios['log2rat_masked']

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

ratios = get_ratios(insertions, ins.x_range.start, ins.x_range.end)

ratios['log2rat_str'] = ratios.log2rat.apply(lambda x: f'{x:.1f}')
ratios['ypos'] = [3.5, 2.5, 1.5]
ratios['xpos'] = [4, 4, 4]

rat_source = ColumnDataSource(ratios)

# Set colors
cmap = cm.get_cmap('PiYG', 256)
PiYG256 = tuple(colors.rgb2hex(cmap(i)[:3]) for i in range(cmap.N))
palette = PiYG256
line_color = linear_cmap(field_name='log2rat', palette=palette,
                         low=-8, high=8)
fill_color = linear_cmap(field_name='log2rat_masked', palette=palette,
                         low=-8, high=8, nan_color='white')

# Plot
rat.circle(x='xpos', y='ypos', size=50, line_color=line_color, line_width=10,
           fill_color=fill_color, source=rat_source)


# Add log2(ratio) text
rat_text = Text(x='xpos', y='ypos', text='log2rat_str',
                x_offset=40, y_offset=10, text_font_size='14pt')
rat.add_glyph(rat_source, rat_text)


# Add legend text
leg_source = ColumnDataSource(dict(x=[5, 5, 5], y=[3.8, 2.8, 1.8],
                                   text=['+ strand', '- strand',
                                         'Both strands']))
leg_glyph = Text(x='x', y='y', text='text', text_font_size='10pt',
                 text_align='center')
rat.add_glyph(leg_source, leg_glyph)


layout = row(column(ins, select), rat)

curdoc().title = 'Insertions per region'
curdoc().add_root(layout)

output_file(filename)
show(layout)
# %%
