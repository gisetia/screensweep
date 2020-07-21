# %%
import os
# import time
import subprocess
from configparser import ConfigParser

# Define parameters for config file and analyze command
screen_name = 'PDL1_IFNg'
assembly = 'hg38'
trim_length = '50'
mode = 'collapse'
start = 'tx'
end = 'tx'
overlap = 'both'
direction = 'sense'

# Define sweep settings
sweep_param = 'end'
step = 500
limit_into_gene = 2000
limit_out_gene = 1000

# %%
# Set sweep range
if sweep_param == 'start':
    sweep_range = range(-limit_out_gene, limit_into_gene + 1, step)
elif sweep_param == 'end':
    sweep_range = range(limit_out_gene, -limit_into_gene - 1, -step)
else:
    raise ValueError('incorrect parameter to sweep')

# Output directory
out_path = (f'analyzed/{screen_name}/{assembly}/{trim_length}/mode={mode}_'
            f'direction={direction}_overlap={overlap}/single-sweep_param='
            f'{sweep_param}_step={step}/')
if not os.path.exists(out_path):
    os.makedirs(out_path)
    print('Creating analyzed directory.')


def get_config_data(config_file: str) -> dict:

    with open(config_file, 'r') as f:
        config_string = '[dummy_section]\n' + f.read()
    config = ConfigParser()
    config.read_string(config_string)
    config_dict = config._sections['dummy_section']

    return config_dict


def write_config_file(config_data: dict, config_file: str) -> None:

    config = ConfigParser()
    config['dummy_section'] = config_data
    txt = '\n'.join(['='.join(item) for item in config.items('dummy_section')])
    with open(config_file, 'w') as f:
        f.write(txt)


# Create config file
config_data = {'bowtie': '/usr/bin/bowtie',
               'assembly': assembly,
               'trim-length': trim_length,
               'threads': '30',
               'screen-dir': '/media/data/nas_scratch/maarten/screens',
               'bowtie-index-hg19': '/references/genomes/BOWTIE_HG19/hg19',
               'bowtie-index-hg38': '/references/genomes/BOWTIE_GRCh38/GCA_'
                                    '000001405.15_GRCh38_no_alt_analysis_set'}
write_config_file(config_data, 'screen-analyzer.conf')

# %% run screen-analyzer

for idx, item in enumerate(sweep_range):

    # sweep parameter setup
    if sweep_param == 'end':
        end = f'tx{item:+}'
    elif sweep_param == 'start':
        start = f'tx{item:+}'

    print(f'Running analysis {idx+1} of {len(sweep_range)}')
    print(f'{sweep_param} tx{item:+}')

    out_file = (f'start={start}_end={end}_')

    # write_config_file(config_data, f'{out_path}conf_{out_file}.conf')

    cmd = (f'./screen-analyzer analyze {screen_name} {assembly} '
           f'--output {out_path}out_{out_file}.txt '
           f'--mode {mode} '
           f'--start {start} '
           f'--end {end} '
           f'--overlap {overlap} '
           f'--direction {direction} '
           f'2> {out_path}counts_{out_file}.txt ')

    # start_time = time.time()

    subprocess.call([cmd], shell=True)

    # print('--- %s seconds ---\n' % (time.time() - start_time))
# %%
