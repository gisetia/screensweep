# %%
from bokeh.plotting import output_file, show
from bokeh.layouts import row, column
from bokeh.models import AutocompleteInput, Div

from sweeptools.analyzesweep import read_analyzed_sweep
from sweeptools.analyzeinsertions import read_insertions, read_refseq

from sweeptools.plotting.sweepplots import link_sweep_and_ins

filename = 'plots/PDL1-sweep.html'

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
ins_data_dir = data_dir
gene = 'CD274'

# Menus
menu_margins = (20, 50, 0, 10)
screen_opts = [params['screen_name']]
screen_menu = AutocompleteInput(title='Screen', value=params['screen_name'],
                                completions=screen_opts, width=150,
                                min_characters=1, case_sensitive=False,
                                margin=menu_margins)
assembly_opts = [params['assembly']]
assembly_menu = AutocompleteInput(title='Assembly', value=params['assembly'],
                                  completions=assembly_opts, width=150,
                                  min_characters=1, case_sensitive=False,
                                  margin=menu_margins)
gene_opts = [gene]
gene_menu = AutocompleteInput(title='Gene', value=gene,
                              completions=gene_opts, width=150,
                              min_characters=1, case_sensitive=False,
                              margin=menu_margins)

txt_out = Div(text='', margin=menu_margins, width=150)

refseq = read_refseq('hg38')

grouped_sweep = read_analyzed_sweep(data_dir, params)
insertions = read_insertions(data_dir, params['screen_name'],
                             params['assembly'], params['trim_length'])

sweep, ins = link_sweep_and_ins(gene, grouped_sweep, params,
                                data_dir, insertions, refseq)

menus = column(screen_menu, assembly_menu, gene_menu, txt_out)
layout = column(row(menus, sweep), ins)

# curdoc().add_root(layout)
output_file(filename)
show(layout)

# %%
