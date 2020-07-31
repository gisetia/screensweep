from math import pi, ceil
from numpy import interp
from bokeh.plotting import figure, output_file, show
from bokeh.models import ColumnDataSource, Arrow, OpenHead, HoverTool, Text
from bokeh.palettes import PiYG8
from bokeh.transform import linear_cmap
from bokeh.layouts import row, column
from bokeh.models.annotations import Title

from .analyzesweep import (read_analyzed_sweep, get_gene_info,
                           get_flags_for_gene)
from .analyzeinsertions import (get_exon_regions, read_gene_insertions, 
                                get_gene_positions)


def plot_insertions(gene, params, data_dir):
    padd = 2000
    gene_pos = get_gene_positions(gene, params['assembly'])
    insertions = read_gene_insertions(gene, data_dir, params, 
                                    gene_pos=gene_pos, padding=padd)
    exons = get_exon_regions(gene_pos)

    # Plot setup --------------------------------------------------------------
    #region -------------------------------------------------------------------

    zero = min(gene_pos['txStart'])
    xlim = (min(gene_pos['txStart']) - padd - zero,
            max(gene_pos['txEnd']) + padd - zero)
    ylim = (0, 17 + len(gene_pos))

    plt = figure(title=f'Insertions in gene: {gene.upper()} - '
                    f'Screen: {params["screen_name"]}',
                x_axis_label='Position (bp)',
                plot_width=1000,
                plot_height=400,

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
    plt.yaxis.ticker = [1, 2, 4, 5, 9, 14]
    plt.yaxis.major_label_overrides = {1: 'Low anti-sense',
                                    2: 'High anti-sense',
                                    4: 'Low sense',
                                    5: 'High sense',
                                    9: 'Sense insertion density',
                                    14: 'Transcript(s)'}

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

    y_off = 15
    # Plot lines and arrows through transcript
    gene_pos = gene_pos.reset_index(drop=True)
    gene_pos['tx_id'] = gene_pos.index + 1
    for idx, tx in gene_pos.iterrows():

        tx_length = tx['txEnd'] - tx['txStart']

        # Position of arrow
        if tx['strand'] == '+':
            x_end = tx['txEnd'] - zero + tx_length/30
            x_srt = x_end - tx_length/30
        else:
            x_end = tx['txStart'] - zero - tx_length/30
            x_srt = x_end + tx_length/30
        ypos = tx['tx_id'] + y_off

        plt.add_layout(Arrow(end=OpenHead(size=6, line_width=2, 
                            line_color='#5f5f61'), 
                            line_color='#5f5f61', line_width=2,
                            x_start=x_srt, y_start=ypos, 
                            x_end=x_end, y_end=ypos))
        plt.line(x=(tx['txStart']-zero, tx['txEnd']-zero), y=(ypos, ypos), 
                color='#5f5f61', line_width=2)

    # Plot exons
    # Transform datasource for exons
    ex_height = {'exCds': 0.8, 'exUtr': 0.5}
    exons['height'] = exons.apply(lambda x: ex_height[x.reg_type], axis=1)

    ex_color = {'exCds': '#fdae61', 'exUtr': '#3288bd'}
    exons['color'] = exons.apply(lambda x: ex_color[x.reg_type], axis=1)

    exons['width'] = exons.apply(lambda x: x.reg_lims[1]-x.reg_lims[0], axis=1)
    exons['start_plot'] = exons.apply(lambda x: x.reg_lims[0]-zero + x.width/2, 
                                    axis=1)
    exons['ypos'] = exons.apply(lambda x: x.tx_id + y_off, axis=1)

    ex_source = ColumnDataSource(exons)

    plt.rect(x='start_plot', y='ypos', width='width', height='height', 
            fill_color='color', line_color=None, source=ex_source)

    #endregion

    # Plot insertions ---------------------------------------------------------
    #region -------------------------------------------------------------------

    # Transform data source for insertions
    ins_colors = {  'hs': PiYG8[0],
                    'ls': PiYG8[-1],
                    'ha': PiYG8[3],
                    'la': PiYG8[-3]}
    def get_ins_color(row, ins_colors):
        return ins_colors[row['chan'][0]+row['dir'][0]]

    insertions['color'] = insertions.apply(lambda row: 
                                        get_ins_color(row, ins_colors),
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

    # Plot insertion density ------------------------------------------------------
    #region -------------------------------------------------------------------

    from scipy import stats

    dens_off = 10
    x = range(xlim[0], xlim[-1] + 1, 10)
    plt.line(x=x, y=dens_off, color='#E2E2E4', line_width=1)

    dens_dict = dict()
    for name, group in insertions.groupby(['chan','dir']):

        if name[1] == 'sense':
            dens = stats.gaussian_kde(group.xpos)
            dens = dens.evaluate(x)
            dens_dict[name[0]] = 3*dens/max(dens)

            dens = dens_off + 3*dens/max(dens)

            dens_col = ins_colors[f'{name[0][0]}s']
            plt.line(x=x, y=dens, color=dens_col, line_width=2)

    dens_diff = dens_off + dens_dict['high'] - dens_dict['low']
    plt.line(x=x, y=dens_diff, color='#7450A2', line_width=2)


    #endregion

    # output_file('plots/test_insertions.html')
    # legend = column(col_leg, p_leg, flag_leg)
    # layout = row(plt, legend)
    # show(plt)

    return plt


def plot_sweep(gene, params, data_dir, grouped_sweep = None, 
               plot_flags = None):

    plot_flags = plot_flags or 1

    if not grouped_sweep:
        print('Loading sweep data.')
        grouped_sweep = read_analyzed_sweep(data_dir, params)

    gene_info = get_gene_info(gene, grouped_sweep)
    
    palette = PiYG8[::-1]
    plot_width = 600

    # Arrange data to plot ----------------------------------------------------
    #region -------------------------------------------------------------------

    src = gene_info.stack().reset_index()

    # Get p values under given threshold
    p_thr = 0.00001
    src['log2_mi_masked'] = src.log2_mi.where(src.p_fdr < p_thr)

    # Rescale log2 MI to get nice sizes in plot
    min_size = plot_width/100
    max_size = plot_width/38
    src['log2_mi_rescaled'] = src['log2_mi'].map(lambda x: interp(abs(x), (0, 8), 
                                                            (min_size, 
                                                            max_size)))

    source = ColumnDataSource(src)
    #endregion

    # Set plot properties -----------------------------------------------------
    #region -------------------------------------------------------------------

    # Axes ranges
    padd = ceil(params['step']/1.5)
    xlim = (-10000 - padd, 2000 + padd)
    ylim = (10000 + padd, -2000 - padd)

    plt = figure(title=f'Gene: {gene.upper()} - Screen: {params["screen_name"]}',
                x_axis_label='End offset (bp)',
                y_axis_label='Start offset (bp)',
                plot_width=plot_width,
                plot_height=plot_width-15,
                match_aspect=True,
                aspect_scale=1,
                x_range=xlim,
                y_range=ylim
                )

    # Rename 0 to position (tx, cds end,start)
    plt.xaxis.major_label_overrides = {0: (f'End '
                                        f'{gene_info["end_pos"].iloc[0,0]}')}
    plt.yaxis.major_label_overrides = {0: (f'Start '
                                        f'{gene_info["srt_pos"].iloc[0,0]}')}
    # Lines to mark 0,0
    plt.line(x=[0, 0], y=[ylim[0], ylim[1]], color='#5f5f61', line_width=2)
    plt.line(x=[xlim[0], xlim[1]], y=[0, 0], color='#5f5f61', line_width=2)

    # Set colormaps
    line_color = linear_cmap(field_name='log2_mi', palette=palette,
                            low=-8 ,high=8)
    fill_color = linear_cmap(field_name='log2_mi_masked', palette=palette,
                            low=-8 ,high=8, nan_color='white')

    # Hover tool
    plt.add_tools(HoverTool(tooltips=[  ('Log2 MI', '@log2_mi'),
                                        ('Start offset', '@srt_off'),
                                        ('End offset', '@end_off'),
                                        ('Low counts', '@low_counts'),
                                        ('High counts', '@high_counts'),
                                        ('P-value', '@p_fdr'),
                                        ],
                            names=['datapoints']
                                #  show_arrow=False
                                ))
    #endregion

    # Plot flags --------------------------------------------------------------
    #region -------------------------------------------------------------------

    if plot_flags:

        slope_thr = 3
        p_ratio_thr= 5
        flags = get_flags_for_gene(gene, grouped_sweep, slope_thr, p_ratio_thr)

        f_col = ['#EBECEF', '#BCBCBF']

        flg_s_src = flags[0].reset_index()
        flg_e_src = flags[1].reset_index()

        flg_s_source = ColumnDataSource(flg_s_src)
        plt.square(x='end_off', y='srt_off', source=flg_s_source, size=21,
                line_color=f_col[0], fill_color=f_col[0], line_width=2)
        
        flg_e_source = ColumnDataSource(flg_e_src)
        plt.square(x='end_off', y='srt_off', source=flg_e_source, size=21,
                line_color=f_col[1], fill_color=f_col[1], line_width=2)

    #endregion

    # Plot data ---------------------------------------------------------------
    #region -------------------------------------------------------------------

    plt.circle(x='end_off', y='srt_off', source=source, 
            size='log2_mi_rescaled', line_width=3,
            line_color=line_color, color=fill_color,
            name='datapoints')
    #endregion

    # Legend ------------------------------------------------------------------
    #region -------------------------------------------------------------------

    def remove_grid_and_ticks(plot):
        plot.xgrid.grid_line_color = None
        plot.xaxis.major_tick_line_color = None
        plot.xaxis.minor_tick_line_color = None
        plot.xaxis.major_label_text_color = None
        plot.xaxis.axis_line_color = None
        plot.ygrid.grid_line_color = None
        plot.yaxis.major_tick_line_color = None
        plot.yaxis.minor_tick_line_color = None
        plot.yaxis.major_label_text_color = None
        plot.yaxis.axis_line_color = None
        plot.outline_line_color = None

    # Color legend
    col_leg = figure(title='Log2 (MI)',
                x_range=(0, 10), y_range=(0, 100),
                plot_width=200, plot_height=250,
                toolbar_location=None)
    remove_grid_and_ticks(col_leg)

    x = [1]*8
    y = [i*12 + 7 for i in range(0,8)]
    size = [7, 5, 3, 1, -1, -3, -5, -7]
    size_rescaled = list(map(lambda x: interp(abs(x), (0, 8),
                                            (min_size, max_size)), size))
    txt = ['-7', '-5', '-3', '-1', '1', '3', '5', '7']

    col_leg_source = ColumnDataSource(dict(x=x, y=y, text=txt, size=size, 
                                    size_r=size_rescaled))

    color = linear_cmap(field_name='size', palette=palette[::-1],
                            low=-8 ,high=8)

    col_leg.circle(x='x', y='y', size='size_r', source=col_leg_source,
                fill_color=color, line_color=color,
                line_width=2)

    glyph = Text(x='x', y='y', text='text', text_font_size='10pt',
                x_offset=15, y_offset=7)
    col_leg.add_glyph(col_leg_source, glyph)

    # P value legend
    p_leg = figure(title='P-value',
                x_range=(0, 10), y_range=(0, 100),
                plot_width=200, plot_height=120,
                toolbar_location=None)
    remove_grid_and_ticks(p_leg)

    x = [1]*2
    y = [i*35 + 40 for i in range(0,2)]
    size = [5, 5]
    size_rescaled = list(map(lambda x: interp(abs(x), (0, 8), 
                                            (min_size, max_size)), size))
    txt = [f'p > {p_thr}', f'p < {p_thr}']
    fill = ['white', 'gray']

    p_leg_source = ColumnDataSource(dict(x=x, y=y, text=txt, size=size, 
                                    size_r=size_rescaled, fill=fill))
    p_leg.circle(x='x', y='y', size='size_r', source=p_leg_source,
                fill_color='fill', line_color='gray',
                line_width=3)
    glyph = Text(x='x', y='y', text='text', text_font_size='10pt',
                x_offset=15, y_offset=7)
    p_leg.add_glyph(p_leg_source, glyph)


    # flags legend
    flag_leg = figure(title=' ',
                x_range=(0, 10), y_range=(0, 100),
                plot_width=200, plot_height=120,
                toolbar_location=None)
    remove_grid_and_ticks(flag_leg)

    if plot_flags:
        t = Title()
        t.text = 'Flags legend'
        flag_leg.title = t
        x = [1]*2
        y = [i*40 + 30 for i in range(0,2)]
        size = [20, 20]
        txt = ['Flag start (y) dir', 'Flag end (x) dir']

        f_leg_source = ColumnDataSource(dict(x=x, y=y, text=txt, 
                                            size=size, fill=f_col))

        flag_leg.square(x='x', y='y', size='size', source=f_leg_source,
                fill_color='fill', line_color='fill')
        glyph = Text(x='x', y='y', text='text', text_font_size='10pt',
                x_offset=15, y_offset=7)
        flag_leg.add_glyph(f_leg_source, glyph)

    #endregion

    # Output ------------------------------------------------------------------
    #region -------------------------------------------------------------------

    # output_file('plots/test_plot.html')
    legend = column(col_leg, p_leg, flag_leg)
    layout = row(plt, legend)
    # show(layout)

    return layout

    #endregion
