# %%
import matplotlib.pyplot as plt
import pandas as pd

filename = '../data/genes/ncbi-genes-hg38.txt'
cols = ['name', 'chrom', 'strand', 'txStart', 'txEnd', 'cdsStart', 'cdsEnd',
        'name2']

genes = pd.read_csv(filename, sep='\\t', usecols=cols)

# Get only coding entries (starting with NM or XM)
genes = genes.query('name.str.startswith("NM") | name.str.startswith("XM") ')

#%%

# Add columns
genes['txSize'] = genes.eval('txEnd - txStart')
genes['cdsSize'] = genes.eval('cdsEnd - cdsStart')


genes['startUtrSize'] = genes.eval('(cdsStart - txStart) * (strand == "+")'
                               '+ (txEnd - cdsEnd) * (strand == "-")')
genes['endUtrSize'] = genes.eval('(cdsStart - txStart) * (strand == "-")'
                               '+ (txEnd - cdsEnd) * (strand == "+")')

# %%

# Stats for longest transcript
cols = ['name2', 'txSize', 'cdsSize','startUtrSize', 'endUtrSize']

longest_tx = genes[cols].groupby('name2').max()

print(longest_tx.describe().apply(lambda s: s.apply(lambda x: f'{x:,}')))

fig, ax = plt.subplots()
longest_tx.boxplot(column=cols[1:5], ax=ax)
ax.set_ylim(0,200000)
# %%
