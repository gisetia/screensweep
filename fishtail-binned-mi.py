# %%
import pandas as pd
import numpy as np
from bokeh.models import (BasicTickFormatter, Legend, HoverTool,
                          TapTool, ColumnDataSource, Div)
from bokeh.plotting import figure, output_file, show
from bokeh.models.callbacks import CustomJS
from bokeh.layouts import row, column

p_thr = 5e-2

screen_name = 'PDL1_IFNg'
assembly = 'hg38'
trim_length = 50
step = 500
data_dir = 'data/analyzed-data'

filename = (f'data/analyzed-data/{screen_name}/{assembly}/{trim_length}'
            f'/out-gene-mi_step={step}')

data = pd.read_parquet(f'{filename}.parquet')


# Manipulate data source
data_plot = data.stack()
data_plot = data_plot.query('p < 0.01')
data_plot = data_plot.reset_index()
data_plot = data_plot.query('strand == "+" | strand == "-"')
data_plot['tot_ins'] = data_plot.apply(lambda x: int(x.low + x.high), axis=1)
data_plot['color'] = data_plot.apply(
    lambda x: 'purple' if x.p_fdr <= p_thr else 'gray', axis=1)
data_plot['alpha'] = data_plot.apply(
    lambda x: 0.5 if x.p_fdr <= p_thr else 0.1, axis=1)

source = ColumnDataSource(data_plot)

# Plot
size = 800
title_font_size = f'{int(size*.018)}pt'
label_font_size = f'{int(size*.015)}pt'

max_ylim = np.ceil(max(abs(data_plot['log2mi']))+0.1)
max_xlim = 10**np.ceil(np.log10(max(data_plot['tot_ins'])))

# create a new plot with a title and axis labels
plt = figure(title=(f'Screen name: {screen_name} - Assembly: {assembly} - '
                    f'Window size: {step:,} bp'),
             x_axis_label='Insertions', y_axis_label='log2 MI',
             plot_width=size, plot_height=int(np.ceil(size/1.5)),
             x_axis_type='log',
             x_range=(1, max_xlim), y_range=(-max_ylim, max_ylim),
             tools='')
plt.xaxis[0].ticker.base = 10
plt.xaxis.formatter = BasicTickFormatter(use_scientific=False)
plt.title.text_font_size = title_font_size
plt.xaxis.axis_label_text_font_size = label_font_size
plt.yaxis.axis_label_text_font_size = label_font_size
plt.xaxis.major_label_text_font_size = label_font_size
plt.yaxis.major_label_text_font_size = label_font_size

points = plt.circle(x='tot_ins',
                    y='log2mi',
                    source=source,
                    color='color',
                    nonselection_alpha='alpha',
                    selection_color='skyblue',
                    selection_alpha=1,
                    alpha='alpha',
                    size=size/100,
                    line_color=None)

plt.add_tools(HoverTool(tooltips=[('Position', '@pos'),
                                  ('Strand', '@strand'),
                                  ('Log2(MI)', '@log2mi'),
                                  ('fcpv', '@p_fdr'),
                                  ('High counts', '@high'),
                                  ('Low counts', '@low')]))

txt_out = Div(text='', margin=(40, 0, 0, 50), width=200)

code = '''
var idx = cb_data.source.selected.indices;
div.text = source.data['pos'][idx]
console.log('Tap: ')
'''

plt.add_tools(TapTool())
plt.select(TapTool).callback = CustomJS(args={'div': txt_out,
                                              'source': source}, code=code)

output_file('plots/test_fishtail.html')
show(row(plt, txt_out))

# %%
