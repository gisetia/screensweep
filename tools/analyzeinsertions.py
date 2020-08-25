import os
import pandas as pd
from ctypes import Structure, c_int8, c_char, c_int32, sizeof
from typing import Optional
from itertools import groupby
from operator import itemgetter
from .utils import timer

pd.options.mode.chained_assignment = None


class CInsertion(Structure):
    _fields_ = [('c', c_int8),
                ('s', c_char),
                ('p', c_int32)]


@timer
def write_insertions(data_dir, outdata_dir, screen_name,
                     assembly, trim_length):

    data_path = (f'''{data_dir}/{screen_name}/'''
                 f'''{assembly}/{trim_length}/''')

    keys = [x for x in list(range(0, 26)) if not x == 23]
    values = [f'chr{i}' for i in range(0, 23)] + ['chrX'] + ['chrY']
    chr_dict = dict(zip(keys, values))

    channels = ['high', 'low']

    insertions = []
    for c in channels:
        filename = f'{data_path}{c}'

        with open(filename, 'rb') as file:
            c_ins = CInsertion()
            while file.readinto(c_ins) == sizeof(c_ins):
                ins = [c, chr_dict[c_ins.c],
                       str(c_ins.s, 'utf-8'), c_ins.p]
                insertions.append(ins)

    insertions = pd.DataFrame(insertions, columns=['chan', 'chr',
                                                   'strand', 'pos'])

    outdata_path = (f'''{outdata_dir}/{screen_name}/'''
                    f'''{assembly}/{trim_length}/''')

    if not os.path.exists(outdata_path):
        os.makedirs(outdata_path)
        print('Creating analyzed directory.')

    insertions.to_parquet(f'{outdata_path}insertions.parquet.snappy',
                          engine='pyarrow', compression='snappy')

    return insertions


@timer
def read_insertions(data_dir, screen_name, assembly, trim_length):
    data_path = (f'''{data_dir}/{screen_name}/'''
                 f'''{assembly}/{trim_length}/''')
    insertions = pd.read_parquet(f'{data_path}insertions.parquet.snappy')

    return insertions


@timer
def read_refseq(assembly: str) -> pd.DataFrame:
    '''Returns data from refseq file as dataframe. (Only protein coding and
    chromosomes without underscore).
    '''
    filename = f'data/genes/ncbi-genes-{assembly}.txt'
    refseq = pd.read_csv(filename, sep='\\t')

    # Get only coding entries (starting with NM or XM)
    refseq = refseq.query('name.str.startswith("NM") '
                          '| name.str.startswith("XM")')

    # Remove alternative chromosomes
    refseq = refseq[~refseq['chrom'].str.contains('_')]

    return refseq


@timer
def get_gene_positions(gene: str, refseq: pd.DataFrame) -> pd.DataFrame:

    grouped_genes = refseq.groupby('name2')
    gene_data = grouped_genes.get_group(gene)

    return gene_data


def read_insertions_region(data_dir: str, screen_name: str, assembly: str,
                           trim_length: int, chrom: str, start: int,
                           end: Optional[int] = None,
                           padd: Optional[int] = 0) -> pd.DataFrame:
    '''Returns dataframe with insertion info for the given chromosome and
    interval defined by start and end:
    chan	chr	    strand	    pos
    high	chr9	-	    	4983150
    high	chr9	+	    	4983164
    '''

    end = end or start

    start = start - padd
    end = end + padd

    data_path = (f'''{data_dir}/{screen_name}/'''
                 f'''{assembly}/{trim_length}/''')

    keys = [x for x in list(range(0, 26)) if not x == 23]
    values = [f'chr{i}' for i in range(0, 23)] + ['chrX'] + ['chrY']
    chr_dict = dict(zip(keys, values))

    channels = ['high', 'low']

    insertions = []
    for c in channels:
        filename = f'{data_path}{c}'

        with open(filename, 'rb') as file:
            c_ins = CInsertion()
            while file.readinto(c_ins) == sizeof(c_ins):
                if (chr_dict[c_ins.c] == f'chr{chrom}'
                    and c_ins.p > start
                        and c_ins.p < end):

                    ins = [c, chr_dict[c_ins.c],
                           str(c_ins.s, 'utf-8'), c_ins.p]
                    insertions.append(ins)

    insertions = pd.DataFrame(insertions, columns=['chan', 'chr',
                                                   'strand', 'pos'])

    return insertions


@timer
def read_gene_insertions(gene: str, insertions: pd.DataFrame,
                         gene_pos: Optional[pd.DataFrame] = None,
                         padding: Optional[int] = None) -> \
        pd.DataFrame:
    '''Returns dataframe with insertion info for the given gene, eg.:
    chan	chr	    strand	dir	        pos
    high	chr9	-	    antisense	4983150
    high	chr9	+	    sense	    4983164
    '''

    padding = padding or 2000

    # if not isinstance(gene_pos, pd.DataFrame):
    #     gene_pos = get_gene_positions(gene, params['assembly'])

    min_pos = min(gene_pos['txStart']) - padding
    max_pos = max(gene_pos['txEnd']) + padding

    if len(gene_pos['chrom'].unique()) > 1:
        print('Warning! Gene is located in more than one chromosome. '
              'Taking only chromosomy of first entry')
    if len(gene_pos['strand'].unique()) > 1:
        print('Warning! Gene is located in more than one strand. '
              'Taking only strand of first entry')

    chrom = gene_pos['chrom'].iloc[0]
    gene_strand = gene_pos['strand'].iloc[0]

    insertions = insertions.query(
        'chr == @chrom & pos >= @min_pos & pos <= @max_pos')
    insertions['dir'] = insertions['strand'].apply(
        lambda x: 'sense' if x == gene_strand else 'antisense')

    return insertions


def contig_list_lims(lst):
    '''Divide a list into contiguous sublists and return the lower and
    upper limits of such sublists.
    '''
    lims = []
    for k, g in groupby(enumerate(lst), lambda x: x[0] - x[1]):
        x = list(map(itemgetter(1), g))
        lims.append([x[0], x[-1]])
    return lims


def get_exon_regions(gene_pos: Optional[pd.DataFrame] = None,
                     gene: Optional[str] = None,
                     assembly: Optional[str] = None) -> pd.DataFrame:
    '''Returns dataframe with start and end positions of all gene exons and
    whether they are cds or utr. Eg.:
    name2	name	        exon_id	tx_id	reg_type	reg_lims
    JAK2	NM_001322195.1	2	    1	    exCds	    [5021987, 5022212]
    JAK2	NM_001322195.1	3	    1	    exCds	    [5029782, 5029905]
    JAK2	NM_004972.3	    126	    6	    exUtr	    [5021962, 5021986]
    '''

    # if not isinstance(gene_pos, pd.DataFrame):
    #     if not gene or not assembly:
    #         raise ValueError('Gene and assembly are required when gen_pos '
    #                          'is not given.')
    #     gene_pos = get_gene_positions(gene, assembly)

    exon = gene_pos.copy(deep=True)
    exon = exon.reset_index(drop=True)
    exon['tx_id'] = exon.index + 1

    cols = ['exonStarts', 'exonEnds']
    exon[cols] = exon[cols].applymap(lambda x: x.split(','))

    exon['cdsRange'] = exon.apply(lambda x: range(x.cdsStart, x.cdsEnd),
                                  axis=1)
    exon['exonRange'] = (exon
                         .apply(lambda t:
                                [range(int(x), int(y)) for i, x
                                 in enumerate(t['exonStarts']) for j, y
                                 in enumerate(t['exonEnds']) if i == j
                                 and x != '' and y != ''], axis=1))

    exon = exon.explode('exonRange')
    exon = exon.reset_index(drop=True)
    exon['exon_id'] = exon.index + 1
    exon = exon.drop(columns=['#bin', 'exonStarts', 'exonEnds', 'score',
                              'cdsStartStat', 'cdsEndStat', 'exonFrames',
                              'cdsStart', 'cdsEnd'])

    exon['exCds_lst'] = exon.apply(lambda t: [x for x in t.exonRange
                                              if x in t.cdsRange], axis=1)
    exon['exUtr_lst'] = exon.apply(lambda t: [x for x in t.exonRange
                                              if x not in t.cdsRange], axis=1)

    exon['exCds'] = exon.apply(lambda t: contig_list_lims(t.exCds_lst), axis=1)
    exon['exUtr'] = exon.apply(lambda t: contig_list_lims(t.exUtr_lst), axis=1)

    exon_regions = exon.melt(id_vars=['name2', 'name', 'exon_id', 'tx_id'],
                             value_vars=['exCds', 'exUtr'],
                             var_name='reg_type', value_name='reg_lims')
    exon_regions = exon_regions.explode('reg_lims').dropna()

    return exon_regions
