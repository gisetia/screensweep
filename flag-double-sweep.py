# %%

from tools.analyzesweep import (read_analyzed_sweep, get_flagged_genes,
                                get_flags_for_gene)

# Define parameters of screen to read
params = {'screen_name': 'IkBa',
        'assembly': 'hg38',
        'trim_length': '50',
        'mode': 'collapse',
        'start': 'tx',
        'end': 'tx',
        'overlap': 'both',
        'direction': 'sense',
        'step': 500}

data_dir = 'data/analyzed-data'
# data_dir = 'sample-data/analyzed-data'

grouped_sweep = read_analyzed_sweep(data_dir, params)

# %% get all gene flags

slope_thr = 3
p_ratio_thr= 6
flagged_genes = get_flagged_genes(grouped_sweep, slope_thr, p_ratio_thr)

# %% get flags for single gene

gene = 'ALG5'
slope_thr = 3
p_ratio_thr= 5

flags = get_flags_for_gene(gene, grouped_sweep, slope_thr, p_ratio_thr)

flags
# %%
