# %%

from bokeh.plotting import figure, output_file, show
from bokeh.layouts import row, column
from bokeh.models import CustomJS, TapTool, ColumnDataSource, Circle, Square
from bokeh.models.annotations import Title

# from tools.plots import plot_sweep, plot_insertions
# from tools.analyzesweep import read_analyzed_sweep
# from tools.analyzeinsertions import get_gene_positions

import tools as tls
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

data_dir = 'data/analyzed-data'
ins_data_dir = 'data/screen-analyzer-data'
gene = 'socs1'

gene = gene.upper()

if not 'grouped_sweep' in locals():
    print('Loading sweep data.')
    grouped_sweep = tls.analyzesweep.read_analyzed_sweep(data_dir, params)

gene_pos = tls.analyzeinsertions.get_gene_positions(gene, params['assembly'])
sweep_layout, sweep_src, sweep_plt = tls.plots.plot_sweep(gene, params, 
                                                data_dir, grouped_sweep)
ins = tls.plots.plot_insertions(gene, params, ins_data_dir, gene_pos)

t = Title()
t.text = ''
ins.title = t

# Get actual end position (rather than just end offset)
tx_pos = (min(gene_pos['txStart']), max(gene_pos['txEnd']))

if gene_pos['strand'].iloc[0] == '+':
    sweep_src['start_pos'] = sweep_src['srt_off']
    sweep_src['end_pos'] = sweep_src.eval('@tx_pos[1] + end_off - @tx_pos[0]')
else:
    sweep_src['end_pos'] = -sweep_src['end_off']
    sweep_src['start_pos'] = sweep_src.eval('@tx_pos[1] - srt_off '
                                            '- @tx_pos[0]')

rect_src = ColumnDataSource({'x': [0, tx_pos[1]-tx_pos[0],
                                  tx_pos[1]-tx_pos[0], 0], 
                             'y': [0,0,6,6]})
rect = ins.patch(x='x', y='y', fill_color='gray', source=rect_src, 
                 line_color=None, alpha=0.15)
rect.visible = False


dummy = figure(title='', plot_width=230, plot_height=100,
             toolbar_location=None)
dummy.outline_line_color = None


code = '''
rect.visible = true;
var idx = cb_data.source.selected.indices;

var s = sweep.data['start_pos'][idx]
var e = sweep.data['end_pos'][idx]

rect.data_source.data['x'] = [e, s, s, e]
rect.data_source.change.emit(); 

console.log('Tap: ' + s + ' ' + e)
'''

sweep_source = ColumnDataSource(sweep_src[['start_pos', 'end_pos']])

sweep_plt.add_tools(TapTool(names=['datapoints']))
sweep_plt.select(TapTool).callback = CustomJS(args={'sweep': sweep_source,
                                                    'rect': rect}, 
                                              code=code)

fig = row(dummy, sweep_layout)
col = column(fig, ins)

output_file('plots/test_plot.html')
show(col)


# %%
