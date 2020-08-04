# %%

from tools.analyzesweep import (read_analyzed_sweep, get_flagged_genes,
                                get_flags_for_gene)

import tools as tls
from importlib import reload
reload(tls)

# Define parameters of screen to read
params = {'screen_name': 'p-AKT',
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

from importlib import reload
reload(tls)

slope_thr = 2
p_ratio_thr = 10
p_thr = 1e-5
flagged_genes = tls.analyzesweep.get_flagged_genes(grouped_sweep, p_thr, 
                                                   slope_thr, p_ratio_thr)

# %% get flags for single gene

gene = 'ALG5'
slope_thr = 4
p_ratio_thr = 6
p_thr = 0.00001

flags = tls.analyzesweep.get_flags_for_gene(gene, grouped_sweep, 
                                            slope_thr, p_ratio_thr)

flags
# %%
