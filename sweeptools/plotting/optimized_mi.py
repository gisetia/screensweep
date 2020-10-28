import pandas as pd
import numpy as np
from math import pi
from typing import Optional, List
from matplotlib import cm, colors
from bokeh.transform import linear_cmap
from bokeh.models import (CrosshairTool, HoverTool, BoxZoomTool, RangeTool,
                          TapTool, ColumnDataSource, Range1d, LinearAxis,
                          Tool, SingleIntervalTicker, PanTool, ResetTool)
from bokeh.plotting import figure


class OptimizedPlot():

    def __init__(self, screen: str, assembly: str, optimized_mi: pd.DataFrame,
                 mode, src: Optional[ColumnDataSource] = None,
                 size: Optional[int] = 750,
                 x_range: Optional[Range1d] = None,
                 tools: Optional[List[Tool]] = None) -> None:

        self.optimized_mi = optimized_mi
        self.screen = screen
        self.assembly = assembly
        self.mode = mode

        if not x_range:
            self.genes = optimized_mi.gene_name
            max_x = min(len(self.genes) - 1, 30)
            x_range = Range1d(-0.5, max_x + 0.5,
                              bounds=(-0.5, len(self.genes) - 0.5),
                              min_interval=max_x)

        self.plt = figure(plot_width=size, plot_height=int(np.ceil(size/3)),
                          x_range=x_range,
                          margin=(20, 80, 0, 0),
                          tools='xwheel_pan, reset',
                          active_scroll='xwheel_pan')

        self.plt.toolbar.logo = None
        self.plt.xaxis.visible = False
        self.plt.xgrid.visible = False

        self.set_colors()

        if not src:
            self.set_source()
            src = self.src

        if not tools:
            self.create_tools()
        else:
            self.plt.add_tools(*tools)

        self.plot_opt(src, self.mode)

    def create_tools(self) -> None:
        crosshair = CrosshairTool(dimensions='height',
                                  line_color='#E9E9EC',
                                  line_width=1)
        taptool = TapTool()
        self.tools = [crosshair, taptool]
        self.plt.add_tools(crosshair, taptool)

    def set_x_axis(self, src) -> None:

        ticker = SingleIntervalTicker(interval=1, num_minor_ticks=0)
        xaxis = LinearAxis(ticker=ticker)
        self.plt.add_layout(xaxis, 'below')
        self.plt.xaxis.major_label_orientation = pi/2

        gene_dict = {int(gene_id): gene_name for gene_id, gene_name in
                     zip(src.data['gene_id'], src.data['gene_name'])}
        self.plt.xaxis.major_label_overrides = gene_dict

        self.plt.height = self.plt.plot_height + 100
        self.plt.min_border_bottom = 150

    def set_colors(self) -> None:
        # Set colors according to log2MI value
        cmap = cm.get_cmap('PiYG', 256)
        PiYG256 = tuple(colors.rgb2hex(cmap(i)[:3]) for i in range(cmap.N))
        palette = PiYG256
        self.tx_color = linear_cmap(field_name='mi_at_tx', palette=palette,
                                    low=-8, high=8)
        self.opt_color = linear_cmap(field_name='log2_mi', palette=palette,
                                     low=-8, high=8)

    def set_source(self) -> None:
        # Set sizes according to log10 of p-value
        max_val = 30  # Max -log10(p) value that will still increase size
        min_size = 5
        max_size = 20

        self.optimized_mi['log_p_tx'] = [-np.log10(x) if x > 0 else 50
                                         for x in self.optimized_mi.p_at_tx]
        self.optimized_mi['p_tx_size'] = np.interp(self.optimized_mi.log_p_tx,
                                                   (0, max_val),
                                                   (min_size, max_size))

        self.optimized_mi['log_opt_p'] = [-np.log10(x) if x > 0 else 50
                                          for x in self.optimized_mi.p]
        self.optimized_mi['opt_p_size'] = np.interp(self.optimized_mi.
                                                    log_opt_p,
                                                    (0, max_val),
                                                    (min_size, max_size))

        self.optimized_mi['gene_id'] = self.optimized_mi.index

        self.optimized_mi['log_counts_tx'] = [np.log10(x) for x in self
                                              .optimized_mi.counts_at_tx]
        self.optimized_mi['log_counts'] = [np.log10(x) for x in self
                                           .optimized_mi.counts]

        self.src = ColumnDataSource(self.optimized_mi)

    def plot_opt(self, src, mode) -> None:

        if mode == 'mi':
            x_tx = 'gene_id'
            x_opt = 'gene_id'
            y_opt = 'log2_mi'
            y_tx = 'mi_at_tx'
            line_y = ('log2_mi', 'mi_at_tx')
            line_x = ('gene_id', 'gene_id')
            self.plt.y_range = Range1d(-10, 10)
            self.plt.yaxis.axis_label = 'Log2 MI'

        elif mode == 'p':
            x_tx = 'gene_id'
            x_opt = 'gene_id'
            y_opt = 'log_opt_p'
            y_tx = 'log_p_tx'
            line_y = ('log_opt_p', 'log_p_tx')
            line_x = ('gene_id', 'gene_id')
            self.plt.y_range = Range1d(0, 53)
            self.plt.yaxis.axis_label = '-log10 p'
            self.plt.line(x=[x for x in self.plt.x_range.bounds],
                          y=[-np.log10(0.05), -np.log10(0.05)],
                          dash=[3, 2], color='#B8B8BB')

        elif mode == 'volcano':
            x_tx = 'mi_at_tx'
            x_opt = 'log2_mi'
            y_opt = 'log_opt_p'
            y_tx = 'log_p_tx'
            line_y = ('log_opt_p', 'log_p_tx')
            line_x = ('log2_mi', 'mi_at_tx')
            self.plt.y_range = Range1d(0, 53, bounds='auto')
            self.plt.x_range = Range1d(-10, 10, bounds='auto')
            self.plt.width = int(self.plt.plot_width/1.5)
            self.plt.xaxis.visible = True

            self.plt.yaxis.axis_label = '-log10 p'
            self.plt.xaxis.axis_label = 'Log2 MI'

        elif mode == 'fishtail':
            x_tx = 'log_counts_tx'
            x_opt = 'log_counts'
            y_opt = 'log2_mi'
            y_tx = 'mi_at_tx'
            line_y = ('log2_mi', 'mi_at_tx')
            line_x = ('log_counts', 'log_counts_tx')

            self.plt.y_range = Range1d(-10, 10, bounds='auto')

            max_x = max(max(self.optimized_mi.log_counts),
                        max(self.optimized_mi.log_counts_tx))

            self.plt.x_range = Range1d(0, round(max_x, 1) + 0.1,
                                       bounds='auto')
            ticker = SingleIntervalTicker(interval=1, num_minor_ticks=0)
            self.plt.xaxis.ticker = ticker
            self.plt.xaxis.major_label_overrides = {0: '1', 1: '10',
                                                    2: '100', 3: '1,000',
                                                    4: '10,000'}

            self.plt.width = int(self.plt.plot_width/1.5)
            self.plt.xaxis.visible = True
            self.plt.yaxis.axis_label = 'Log2 MI'
            self.plt.xaxis.axis_label = 'Insertions'

        def plot_lines(gene_flag):
            self.plt.line(y=[gene_flag[line_y[0]], gene_flag[line_y[1]]],
                          x=[gene_flag[line_x[0]], gene_flag[line_x[1]]],
                          color='#828283', name='opt_line')
        self.optimized_mi.apply(plot_lines, axis=1)

        self.plt.square(x=x_tx, y=y_tx, source=src,
                        color=self.tx_color, size='p_tx_size',
                        name='to_inspect', alpha=1,
                        line_color='white', line_width=1,
                        selection_line_color='gray', selection_line_width=2,
                        nonselection_alpha=1)
        self.plt.circle(x=x_opt, y=y_opt, source=src,
                        color=self.opt_color, size='opt_p_size',
                        name='to_inspect', alpha=1,
                        line_color='white', line_width=1,
                        selection_line_color='gray', selection_line_width=2,
                        nonselection_alpha=1)

        if mode in ['mi', 'p']:
            self.plt.add_tools(HoverTool(
                tooltips=[('Gene', '@gene_name'),
                          ('Opt log2 MI', '@log2_mi'),
                          ('Tx log2 MI', '@mi_at_tx'),
                          ('Opt p-value', '@p'),
                          ('Tx p-value', '@p_at_tx'),
                          ('Opt start offset', '@srt_off'),
                          ('Opt end offset', '@end_off'),
                          ('Opt high counts', '@high_counts'),
                          ('Opt low counts', '@low_counts'),
                          ],
                names=['to_inspect'],
                mode='mouse', point_policy='snap_to_data'))

        if mode in ['volcano', 'fishtail']:
            self.plt.toolbar.active_inspect = None

            self.plt.select(name='to_inspect').glyph.size = 8
            self.plt.select(name='opt_line').glyph.line_alpha = 0

            self.plt.select(name='to_inspect').glyph.fill_alpha = 0.2
            self.plt.select(name='to_inspect').glyph.line_alpha = 0
            self.plt.select(name='to_inspect').selection_glyph.fill_alpha = 1
            # self.plt.select(
            #     name='to_inspect').selection_glyph.line_color = None
            self.plt.select(
                name='to_inspect').nonselection_glyph.fill_alpha = 0.2
            self.plt.select(
                name='to_inspect').nonselection_glyph.line_alpha = 0

            self.plt.xgrid.visible = True

            self.plt.tools = [BoxZoomTool(), PanTool(), ResetTool(), TapTool()]
            self.plt.toolbar.active_multi = self.plt.select_one(BoxZoomTool)

    def for_range_selection(self, other_plot: figure) -> None:

        (a, b) = other_plot.x_range.bounds
        self.plt.x_range = Range1d(a, b, bounds='auto')

        self.plt.toolbar.active_inspect = None
        self.plt.toolbar_location = None
        self.plt.plot_height = 50
        self.plt.xgrid.grid_line_color = None
        self.plt.ygrid.grid_line_color = None
        self.plt.xaxis.visible = False
        self.plt.yaxis.visible = False

        range_tool = RangeTool(x_range=other_plot.x_range)
        range_tool.overlay.fill_color = 'navy'
        range_tool.overlay.fill_alpha = 0.1
        self.plt.add_tools(range_tool)
        self.plt.toolbar.active_multi = range_tool

        self.plt.select(name='to_inspect').glyph.size = 3
        self.plt.select(name='to_inspect').selection_glyph.fill_color = 'gray'
        self.plt.select(name='to_inspect').selection_glyph.line_color = None

        self.plt.select(name='opt_line').glyph.line_color = '#E9E9EC'

    def add_title(self) -> None:
        self.plt.height = self.plt.plot_height + 20
        self.plt.title.text = (f'Flagged genes for screen: {self.screen} -'
                               f' Assembly: {self.assembly}')

    def add_background(self, background: pd.DataFrame) -> None:

        background['log_counts'] = np.log10(background.high_counts
                                            + background.low_counts)
        background['log_p'] = [-np.log10(x) if x > 0 else 50
                               for x in background.p]

        src = ColumnDataSource(background[['log2_mi', 'log_p', 'log_counts']])

        if self.mode == 'volcano':
            x = 'log2_mi'
            y = 'log_p'
        elif self.mode == 'fishtail':
            x = 'log_counts'
            y = 'log2_mi'

        bg = self.plt.circle(x=x, y=y, source=src,
                             fill_color='gray', fill_alpha=0.2,
                             line_color=None, size=6)
        bg.level = 'underlay'

        if self.mode == 'fishtail':
            max_x = max(background.log_counts)
            self.plt.x_range.end = max_x
