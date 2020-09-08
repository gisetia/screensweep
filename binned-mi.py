# %%
import pandas as pd
import numpy as np
from pyarrow import Table, parquet
from scipy.stats import rankdata, fisher_exact
from math import log2
from fisher import pvalue_npy

from tools.analyzeinsertions import read_insertions, read_refseq

screen_name = 'PDL1_IFNg'
assembly = 'hg38'
trim_length = 50

step = 1000
data_dir = 'data/analyzed-data'


ins = read_insertions(data_dir, screen_name, assembly, trim_length)
refseq = read_refseq(assembly, name_chrom=True)

# %% Remove insertions in genes
# %%time


def collapse_gene_refseq(gene_refseq):
    # assert len(gene_refseq.chrom.unique()) == 1, (
    #     f'Gene {gene_refseq.name2.iloc[0]} is in more than one chromosome.')

    gene_pos = pd.Series({'name': gene_refseq.name2.iloc[0],
                          'chrom': gene_refseq.chrom.iloc[0],
                          'txStart': gene_refseq.txStart.min(),
                          'txEnd': gene_refseq.txEnd.max()})

    return gene_pos


def drop_ins_in_genes(ins_chr, chr, coll_refseq):
    print(f'Removing insterions in genes at {chr}.')

    pos = ins_chr.pos.to_numpy()

    refseq_chr = coll_refseq.query('chrom == @chr')
    starts = refseq_chr.txStart.to_numpy()
    ends = refseq_chr.txEnd.to_numpy()

    mask = np.invert(list(map(lambda x:
                              np.any((x >= starts) & (x < ends)), pos)))
    ins_chr_out = ins_chr[mask]

    ins_chr_out = ins_chr_out.drop(columns='chr')

    return ins_chr_out


coll_refseq = refseq.groupby('name_chr').apply(collapse_gene_refseq)
# ins_out = ins.groupby('chr').apply(lambda x:
#                                    drop_ins_in_genes(x, x.name, coll_refseq))

ins_out = ins

# %% Get insertion counts for genome intervals
# %%time


def bin_insertions_slow(insertions, step):

    # Get bins to group data
    max_pos = insertions.pos.max()
    bins = pd.IntervalIndex(pd.interval_range(start=0, end=max_pos+step,
                                              freq=step), closed='left')

    # Bin data and count
    counts = insertions.groupby([pd.cut(insertions.pos, bins=bins,
                                        include_lowest=True),
                                 'chan', 'strand']).size().unstack()
    counts['both'] = counts.sum(axis=1)

    # Reshape
    counts = counts.unstack().swaplevel(i=0, j=1, axis=1).sort_index(axis=1)
    counts.reset_index(inplace=True)

    # Reset index to bin number
    counts = counts.set_index([counts.reset_index().index + 1])
    counts.index.names = ['bin']

    return counts


def bin_insertions(insertions, step):

    # Get bins to group data
    max_pos = insertions.pos.max()
    bins = range(0, max_pos+step, step)

    # Bin data and get bin end positions
    bins_idx = np.digitize(insertions.pos.to_numpy(), bins, right=False)
    insertions['pos'] = list(map(bins.__getitem__, bins_idx))

    # Count insertions per bin
    counts = insertions.groupby(['pos', 'chan', 'strand']).size().unstack()
    counts['both'] = counts.sum(axis=1)

    # Reshape
    counts = counts.unstack().swaplevel(i=0, j=1, axis=1).sort_index(axis=1)
    counts.reset_index(inplace=True)

    # Reset index to bin number
    counts = counts.set_index([counts.reset_index().index + 1])
    counts.index.names = ['bin']

    counts = counts.fillna(0).astype(int)

    return counts


def rename_bins_slow(counts):

    counts = counts.reset_index()

    # Change pos from interval to string containing start and end of interval
    counts['pos'] = counts.pos.apply(lambda x:
                                     f'{int(x.left)+1:,}-{int(x.right):,}')
    counts['pos'] = counts[['chr', 'pos']].agg(':'.join, axis=1)

    # Set index
    counts = counts.set_index(['chr', 'bin', 'pos'])
    counts.index.names = ['chr', 'bin', 'pos']

    return counts


def rename_bins(counts):

    counts = counts.reset_index()

    # Change pos from range item to string containing start and end
    counts['pos'] = counts.pos.apply(lambda x: f'{x-step+1:,}-{x:,}')
    counts['pos'] = counts[['chr', 'pos']].agg(':'.join, axis=1)

    # Set index
    counts = counts.set_index(['chr', 'bin', 'pos'])
    counts.index.names = ['chr', 'bin', 'pos']

    return counts


def get_binned_counts_slow(insertions, step):
    counts = insertions.groupby('chr').apply(func=bin_insertions_slow,
                                             step=step)
    counts = counts.pipe(rename_bins_slow)

    return counts


def get_binned_counts(insertions, step):
    counts = insertions.groupby('chr').apply(func=bin_insertions, step=step)
    counts = counts.pipe(rename_bins)

    return counts


# counts_slow = get_binned_counts_slow(ins_out, step)
counts = get_binned_counts(ins_out, step)

# %% Get MI data
# %%time


def get_log2mi(counts):

    # Replace zeros to avoid division by 0
    counts = counts.replace(0, 1)

    log2mi = (((counts['high']/(counts['high'].sum() - counts['high'])) /
               (counts['low']/(counts['low'].sum() - counts['low'])))
              .applymap(log2))

    return log2mi


def get_p(group, strand, counts):

    tot_counts = counts[['high', 'low']].sum()
    tot_high = np.ones(len(group)).astype(np.uint64)*tot_counts['high'][strand]
    tot_low = np.ones(len(group)).astype(np.uint64)*tot_counts['low'][strand]

    group_high = group['high'][strand].to_numpy().astype(np.uint64)
    group_low = group['low'][strand].to_numpy().astype(np.uint64)

    _, _, p = pvalue_npy(group_high, tot_high,
                         group_low, tot_low)

    # Vectorized function 'pvalue_npy' from 'fisher' module is fast but has
    # precision issues with low p-values (lower than approx 1e-5). An ugly fix
    # is to replace low p-values with those obtained by 'fisher_extact' from
    # 'scipy'.
    inds = np.asarray(p < 1e-5).nonzero()[0]
    for i in inds:
        p[i] = fisher_exact([[group_high[i], tot_high[i]],
                             [group_low[i], tot_low[i]]])[1]

    p = pd.Series(p, index=group.index)

    return p


def fdr(p_vals):

    ranked_p_values = rankdata(p_vals)
    fdr = p_vals * len(p_vals) / ranked_p_values
    fdr[fdr > 1] = 1

    return fdr


def get_mi_data(counts):

    # Get log2mi
    print('Getting log2mi.')
    log2mi = get_log2mi(counts)

    # Get p value
    print('Getting p-values.')
    p = (counts[['high', 'low']].groupby(axis=1, level='strand')
                                .apply(lambda x: get_p(x, x.name, counts)))

    # Get corrected p value
    print('Getting corrected p-values.')
    p_fdr = p.apply(fdr)

    # Set multiindex columns to be able to concatenate with counts dataframe
    log2mi.columns = pd.MultiIndex.from_product([['log2mi'], log2mi.columns])
    p.columns = pd.MultiIndex.from_product([['p'], p.columns])
    p_fdr.columns = pd.MultiIndex.from_product([['p_fdr'], p_fdr.columns])

    # Concatenate into single dataframe
    data = pd.concat([counts, log2mi, p, p_fdr], axis=1)
    data.columns.names = ['quantity', 'strand']

    return data


data = get_mi_data(counts)

# %% Save data

filename = (f'data/analyzed-data/{screen_name}/{assembly}/{trim_length}'
            f'/out-gene-mi_step={step}')

table = Table.from_pandas(data)
parquet.write_table(table, f'{filename}.parquet')

# df_test_read = pd.read_parquet(f'{filename}test.parquet')

# %%
