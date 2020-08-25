# %%

from tools.analyzeinsertions import write_insertions

screen_name = 'PDL1_IFNg'
assembly = 'hg38'
trim_length = 50

indata_dir = 'data/screen-analyzer-data'
outdata_dir = 'data/analyzed-data'

insertions = write_insertions(indata_dir, outdata_dir, screen_name,
                              assembly, trim_length)

# %%
