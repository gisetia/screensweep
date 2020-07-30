import pandas as pd
from ctypes import Structure, c_int8, c_char, c_int32, sizeof
from typing import List, Union, Optional
from dataclasses import dataclass
from .utils import timer

@dataclass
class Insertion:
    ch: str
    chr: str
    strand: str
    pos: int
    dir: Union[None, str] = None

class CInsertion(Structure):
    _fields_ = [('c', c_int8),
                ('s', c_char),
                ('p', c_int32)]

@timer
def read_all_insertions(data_dir: str, params: dict) -> List[Insertion]:

    data_path = (f'''{data_dir}/{params['screen_name']}/'''
                f'''{params['assembly']}/{params['trim_length']}/''')

    keys = [x for x in list(range(0,26)) if not x == 23]
    values = [f'chr{i}' for i in range(0,23)] + ['chrX'] + ['chrY']
    chr_dict = dict(zip(keys, values))

    channels = ['high', 'low']

    insertions = []
    for c in channels:
        filename = f'{data_path}{c}'

        with open(filename, 'rb') as file:
            c_ins = CInsertion()
            while file.readinto(c_ins) == sizeof(c_ins):
                ins = Insertion(c, chr_dict[c_ins.c], 
                                str(c_ins.s, 'utf-8'), c_ins.p)
                insertions.append(ins)

    return insertions

# @timer
def get_gene_positions(gene: str, assembly: dict) -> pd.DataFrame:

    filename = f'data/genes/ncbi-genes-{assembly}.txt'
    genes = pd.read_csv(filename, sep='\\t')

    # Get only coding entries (starting with NM or XM)
    genes = genes.query('name.str.startswith("NM") '
                        '| name.str.startswith("XM")')

    grouped_genes = genes.groupby('name2')
    gene_data = grouped_genes.get_group(gene)

    return gene_data

# noooo
@timer
def get_gene_insertions(gene: str, assembly: str, insertions: List[Insertion],
                        padding: Optional[int] = None) -> \
                        List[Insertion]:

    padding = padding or 2000

    gene_pos = get_gene_positions(gene, assembly)

    min_pos = min(gene_pos['txStart']) - padding
    max_pos = max(gene_pos['txEnd']) + padding

    if len(gene_pos['chrom'].unique()) > 1:
        print(f'Warning! Gene is located in more than one chromosome. '
              'Taking only chromosomy of first entry')
    if len(gene_pos['strand'].unique()) > 1:
        print(f'Warning! Gene is located in more than one strand. '
              'Taking only strand of first entry')

    chrom = gene_pos['chrom'].iloc[0]
    gene_strand = gene_pos['strand'].iloc[0]

    gene_ins = [x for x in insertions if (x.chr == chrom 
                                        and x.pos > min_pos 
                                        and x.pos < max_pos)]

    for i in gene_ins:
        if i.strand == gene_strand:
            i.dir = 'sense' 
        else:
            i.dir = 'antisense'
    
    return gene_ins

# noo
@timer
def read_gene_insertions_noo(gene: str, data_dir: str, params: dict,
                         gene_pos: Optional[pd.DataFrame] = None,
                         padding: Optional[int] = None) -> \
                         List[Insertion]:

    padding = padding or 2000

    if not isinstance(gene_pos, pd.DataFrame):
        gene_pos = get_gene_positions(gene, params['assembly'])

    min_pos = min(gene_pos['txStart']) - padding
    max_pos = max(gene_pos['txEnd']) + padding

    if len(gene_pos['chrom'].unique()) > 1:
        print(f'Warning! Gene is located in more than one chromosome. '
              'Taking only chromosomy of first entry')
    if len(gene_pos['strand'].unique()) > 1:
        print(f'Warning! Gene is located in more than one strand. '
              'Taking only strand of first entry')

    chrom = gene_pos['chrom'].iloc[0]
    gene_strand = gene_pos['strand'].iloc[0]

    data_path = (f'''{data_dir}/{params['screen_name']}/'''
                f'''{params['assembly']}/{params['trim_length']}/''')

    keys = [x for x in list(range(0,26)) if not x == 23]
    values = [f'chr{i}' for i in range(0,23)] + ['chrX'] + ['chrY']
    chr_dict = dict(zip(keys, values))

    channels = ['high', 'low']

    insertions = []
    for c in channels:
        filename = f'{data_path}{c}'

        with open(filename, 'rb') as file:
            c_ins = CInsertion()
            while file.readinto(c_ins) == sizeof(c_ins):
                if (chr_dict[c_ins.c] == chrom 
                    and c_ins.p > min_pos
                    and c_ins.p < max_pos):
                    ins = Insertion(c, chr_dict[c_ins.c], 
                                    str(c_ins.s, 'utf-8'), c_ins.p)

                    if ins.strand == gene_strand:
                        ins.dir = 'sense' 
                    else:
                        ins.dir = 'antisense'

                    insertions.append(ins)

    return insertions

@timer
def read_gene_insertions(gene: str, data_dir: str, params: dict,
                         gene_pos: Optional[pd.DataFrame] = None,
                         padding: Optional[int] = None) -> \
                         List[Insertion]:

    padding = padding or 2000

    if not isinstance(gene_pos, pd.DataFrame):
        gene_pos = get_gene_positions(gene, params['assembly'])

    min_pos = min(gene_pos['txStart']) - padding
    max_pos = max(gene_pos['txEnd']) + padding

    if len(gene_pos['chrom'].unique()) > 1:
        print(f'Warning! Gene is located in more than one chromosome. '
              'Taking only chromosomy of first entry')
    if len(gene_pos['strand'].unique()) > 1:
        print(f'Warning! Gene is located in more than one strand. '
              'Taking only strand of first entry')

    chrom = gene_pos['chrom'].iloc[0]
    gene_strand = gene_pos['strand'].iloc[0]

    data_path = (f'''{data_dir}/{params['screen_name']}/'''
                f'''{params['assembly']}/{params['trim_length']}/''')

    keys = [x for x in list(range(0,26)) if not x == 23]
    values = [f'chr{i}' for i in range(0,23)] + ['chrX'] + ['chrY']
    chr_dict = dict(zip(keys, values))

    channels = ['high', 'low']

    insertions = []
    for c in channels:
        filename = f'{data_path}{c}'

        with open(filename, 'rb') as file:
            c_ins = CInsertion()
            while file.readinto(c_ins) == sizeof(c_ins):
                if (chr_dict[c_ins.c] == chrom 
                    and c_ins.p > min_pos
                    and c_ins.p < max_pos):

                    ins = [c, chr_dict[c_ins.c], 
                           str(c_ins.s, 'utf-8'),
                           ('sense' if str(c_ins.s, 'utf-8') == gene_strand 
                           else 'antisense'), c_ins.p,]
                    insertions.append(ins)
    
    insertions = pd.DataFrame(insertions, columns=['chan', 'chr', 'strand',
                                                   'dir', 'pos'])

    return insertions