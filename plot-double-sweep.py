# %%

from tools.plots import plot_sweep, plot_insertions
from tools.analyzesweep import read_analyzed_sweep

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
gene = 'CD274'

from bokeh.plotting import figure, output_file, show
from bokeh.layouts import row, column

if not 'grouped_sweep' in locals():
    print('Loading sweep data.')
    grouped_sweep = read_analyzed_sweep(data_dir, params)

sweep = plot_sweep(gene, params, data_dir, grouped_sweep)
ins = plot_insertions(gene, params, ins_data_dir)

output_file('plots/test_plot.html')
fig = row(sweep, ins)
show(fig)

# %%
