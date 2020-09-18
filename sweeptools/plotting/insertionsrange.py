
from math import pi
from typing import Optional
from bokeh.plotting import figure
from bokeh.models import (ColumnDataSource, HoverTool, CrosshairTool,
                          RangeTool, Range1d, LinearAxis, NumeralTickFormatter)
from bokeh.palettes import PiYG8

from ..analyzeinsertions import (read_insertions_region)


def plot_insertions(data_dir: str, screen_name: str, assembly: str,
                    trim_length: int, chrom: str, start: int,
                    end: Optional[int] = None,
                    padd: Optional[int] = None) -> figure:

    padd = padd or 5000
    insertions = read_insertions_region(data_dir, screen_name, assembly,
                                        trim_length, chrom, start, end, padd)

    ylim = (0, 6)
    ins = figure(title=f'Insertions of screen {screen_name} at '
                 f'chr{chrom}:{start:,} - {end:,}',
                 plot_width=1000,
                 plot_height=300,
                 x_range=(start, end),
                 y_range=ylim,
                 tools='reset, save')
    ins.ygrid.grid_line_color = None
    ins.xgrid.grid_line_color = None
    ins.yaxis.major_tick_line_color = None
    ins.yaxis.minor_tick_line_color = None
    ins.yaxis.axis_line_color = None
    ins.outline_line_color = None
    ins.xaxis.formatter = NumeralTickFormatter(format='0,0')
    ins.yaxis.ticker = [1, 2, 4, 5, 10, 16]
    ins.yaxis.major_label_overrides = {1: 'Low - strand',
                                       2: 'Low + strand',
                                       4: 'High - strand',
                                       5: 'High + strand'}

    # Transform data source for insertions
    ins_colors = {'h+': PiYG8[0],
                  'l+': PiYG8[-1],
                  'h-': PiYG8[2],
                  'l-': PiYG8[-3]}
    ins_pos = {'h+': 5,
               'l+': 2,
               'h-': 4,
               'l-': 1}

    insertions['color'] = insertions.apply(lambda row:
                                           ins_colors[row['chan']
                                                      [0]+row['strand'][0]],
                                           axis=1)
    insertions['ypos'] = insertions.apply(lambda row:
                                          ins_pos[row['chan'][0]
                                                  + row['strand'][0]],
                                          axis=1)
    insertions['xpos'] = insertions.apply(lambda row: row['pos'],
                                          axis=1)
    source_ins = ColumnDataSource(insertions)

    # Plot insertions
    x_line = (start - padd, end + padd)
    ins.line(x=x_line, y=[1, 1], color='#BCBCBF', line_width=2, name='line')
    ins.line(x=x_line, y=[2, 2], color='#BCBCBF', line_width=2)
    ins.line(x=x_line, y=[4, 4], color='#BCBCBF', line_width=2)
    ins.line(x=x_line, y=[5, 5], color='#BCBCBF', line_width=2)

    ins.dash(x='xpos', y='ypos', color='color', source=source_ins,
             angle=pi/2, line_width=1, size=15, name='insertions')

    ins.line(x=(start, start), y=ylim, color='#E2E2E4', line_dash=[5, 2])
    ins.line(x=(end, end), y=ylim, color='#E2E2E4', line_dash=[5, 2])

    # Add dummy line for hover tool
    x = list(range(start - padd, end + padd + 1))
    y = [3 for i in x]
    ins.line(x=x, y=y, line_color='white', name='needshover')

    ins.add_tools(HoverTool(tooltips=[('Position', '$x{0,0}')], mode='vline',
                            line_policy='nearest',
                            names=['needshover']
                            ),
                  CrosshairTool(dimensions='height',
                                line_color='#363638',
                                line_width=1))

    return ins


def ins_select_range(ins: figure) -> figure:

    line = [x for x in ins.renderers if x.name == 'line'][0]
    x_range = line.data_source.data['x']

    padd = ins.x_range.start - x_range[0]

    select = figure(plot_height=150, plot_width=ins.frame_width,
                    y_range=ins.y_range,
                    x_range=x_range,
                    margin=(40, 0, 0, 0),
                    x_axis_location='above',
                    x_axis_label='Absolute position',
                    tools='', toolbar_location=None)
    select.extra_x_ranges = {'relative_pos': Range1d(start=-padd,
                                                     end=ins.x_range.end
                                                     - ins.x_range.start+padd)}
    select.add_layout(LinearAxis(x_range_name='relative_pos',
                                 axis_label='Relative position'), 'below')
    select.xaxis.formatter = NumeralTickFormatter(format='0,0')
    select.yaxis.visible = False
    select.xgrid.grid_line_color = None

    # Plot insertions
    select.line(x=x_range, y=[1, 1], color='#BCBCBF', line_width=1)
    select.line(x=x_range, y=[2, 2], color='#BCBCBF', line_width=1)
    select.line(x=x_range, y=[4, 4], color='#BCBCBF', line_width=1)
    select.line(x=x_range, y=[5, 5], color='#BCBCBF', line_width=1)

    source_ins = [x.data_source for x in ins.renderers
                  if x.name == 'insertions'][0]
    select.dash(x='xpos', y='ypos', color='color', source=source_ins,
                angle=pi/2, line_width=1, size=5,)

    select.line(x=(ins.x_range.start, ins.x_range.start),
                y=(ins.y_range.start, ins.y_range.end),
                color='#E2E2E4', line_dash=[5, 2])
    select.line(x=(ins.x_range.end, ins.x_range.end),
                y=(ins.y_range.start, ins.y_range.end),
                color='#E2E2E4', line_dash=[5, 2])

    select.ygrid.grid_line_color = None

    range_tool = RangeTool(x_range=ins.x_range)
    range_tool.overlay.fill_color = "navy"
    range_tool.overlay.fill_alpha = 0.1
    select.add_tools(range_tool)
    select.toolbar.active_multi = range_tool

    return select
