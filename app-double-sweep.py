# %%

from bokeh.plotting import figure, output_file, show, curdoc
from bokeh.layouts import row, column
from bokeh.models import AutocompleteInput, Div

# import tools as tls

from tools.analyzesweep import read_analyzed_sweep
from tools.analyzeinsertions import read_insertions, read_refseq
# from tools.plots import link_sweep_and_ins

from tools.plotting.sweepplots import link_sweep_and_ins


# from tools.utils import timer
# from importlib import reload
# reload(tls)

# Define parameters of screen to read
params = {'screen_name': 'Ac-gamma-Actin',
          'assembly': 'hg38',
          'trim_length': '50',
          'mode': 'collapse',
          'start': 'tx',
          'end': 'tx',
          'overlap': 'both',
          'direction': 'sense',
          'step': 500}

data_dir = '../data/sweeps-analyzed'
ins_data_dir = '../data/screen-insertions'

# data_dir = 'data/analyzed-data'
# ins_data_dir = 'data/screen-analyzer-data'
# ins_data_dir = data_dir
# gene = 'SOAT1'

gene_opts = []

# Menus
menu_margins = (20, 50, 0, 10)
screen_opts = ['PDL1_IFNg', 'p-AKT', 'Ac-gamma-Actin', 'IkBa']
screen_menu = AutocompleteInput(title='Screen', value='',
                                completions=screen_opts, width=150,
                                min_characters=1, case_sensitive=False,
                                margin=menu_margins)
assembly_opts = ['hg38', 'hg19']
assembly_menu = AutocompleteInput(title='Assembly', value='hg38',
                                  completions=assembly_opts, width=150,
                                  min_characters=1, case_sensitive=False,
                                  margin=menu_margins)
gene_menu = AutocompleteInput(title='Gene', value='',
                              completions=gene_opts, width=150,
                              min_characters=1, case_sensitive=False,
                              margin=menu_margins)

# Callbacks


def load_screen(attr, old, new):
    txt_out.text = 'Loading screen...'
    curdoc().add_next_tick_callback(update_screen)


def update_screen():
    screen = screen_menu.value
    assembly = assembly_menu.value
    params['screen_name'] = screen
    params['assembly'] = assembly
    global grouped_sweep
    global insertions

    try:
        grouped_sweep = read_analyzed_sweep(data_dir, params)

        # need to update reading insertions
        if params['screen_name'] == 'PDL1_IFNg':
            insertions = read_insertions(data_dir, params['screen_name'],
                                     params['assembly'], params['trim_length'])
        curdoc().add_next_tick_callback(update_gene_menu)
    except OSError:
        txt_out.text = 'No data found for these parameters.'


def update_gene_menu():
    txt_out.text = 'Finished loading.'
    gene_opts = list(grouped_sweep.groups.keys())
    gene_menu.completions = gene_opts


def load_assembly(attr, old, new):
    txt_out.text = 'Loading assembly...'
    assembly = assembly_menu.value
    global refseq
    refseq = read_refseq('../data/refseq', assembly)
    curdoc().add_next_tick_callback(update_screen)


def load_gene(attr, old, new):
    txt_out.text = 'Loading gene...'
    curdoc().add_next_tick_callback(update_gene)


def update_gene():
    txt_out.text = 'Finished loading gene.'
    gene = gene_menu.value
    sweep, ins = link_sweep_and_ins(gene, grouped_sweep, params,
                                    data_dir, insertions, refseq)
    layout.children[0].children[1] = sweep
    layout.children[1] = ins


screen_menu.on_change('value', load_screen)
assembly_menu.on_change('value', load_assembly)
gene_menu.on_change('value', load_gene)

# Load hg38 refseq by default
refseq = read_refseq('../data/refseq', 'hg38')

# Initialize empty figures
sweep = figure(plot_width=600, plot_height=600, toolbar_location=None)
sweep.outline_line_color = None
ins = figure(plot_width=1000, plot_height=400, toolbar_location=None)
ins.outline_line_color = None

txt_out = Div(text='', margin=menu_margins, width=150)

menus = column(screen_menu, assembly_menu, gene_menu, txt_out)
layout = column(row(menus, sweep), ins)

curdoc().add_root(layout)
output_file('plots/test_plot.html')
show(layout)

# %%
