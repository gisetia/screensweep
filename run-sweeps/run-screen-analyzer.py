import os
import subprocess
from configparser import ConfigParser
import pandas as pd

######
# Missing 2020-09-01 hg38
# ['WT_EF1a-SD_pS6', 'WT_EF1a-SD_p-p38_aniso', 'NLRP3-GFP-washout', 'MOSPD2-in-WT', 'BODIPY-in-HILPDA-KO-OAchase']
######


# Define parameters for config file and analyze command
screen_name = None
assembly = 'hg38'
trim_length = '50'
mode = 'collapse'
overlap = 'neither'
direction = 'sense'
start = 'tx'
end = 'tx'

# Load and filter screens
screens = pd.read_csv('screen-details/raw_screen-list_2020-09-01.csv')
screens = screens.query('type == "ip"')

screen_list = screens.screen_name.values


# Create config file
def write_config_file(config_data: dict, config_file: str) -> None:

    config = ConfigParser()
    config['dummy_section'] = config_data
    txt = '\n'.join(['='.join(item) for item in config.items('dummy_section')])
    with open(config_file, 'w') as f:
        f.write(txt)


config_data = {'bowtie': '/usr/bin/bowtie',
               'assembly': assembly,
               'trim-length': trim_length,
               'threads': '30',
               'screen-dir': '/media/data/nas_scratch/maarten/screens',
               'bowtie-index-hg19': '/references/genomes/BOWTIE_HG19/hg19',
               'bowtie-index-hg38': '/references/genomes/BOWTIE_GRCh38/GCA_'
                                    '000001405.15_GRCh38_no_alt_analysis_set'}
write_config_file(config_data, 'screen-analyzer.conf')

# Run screen-analyzer
missing_screens = []
analyzed_screens = []

for screen_name in screen_list:

    print(f'\nRunning screen-analyzer analyze for screen {screen_name}')

    out_path = (f'analyzed_data/{screen_name}/{assembly}/{trim_length}/'
                f'mode={mode}_direction={direction}_overlap={overlap}/')

    if not os.path.exists(out_path):
        os.makedirs(out_path)
        print(f'Creating analyzed directory for screen {screen_name}.')

    out_file = (f'start={start}_end={end}_')

    cmd = (f'./screen-analyzer analyze {screen_name} {assembly} '
           f'--output {out_path}out_{out_file}.txt '
           f'--mode {mode} '
           f'--start {start} '
           f'--end {end} '
           f'--overlap {overlap} '
           f'--direction {direction} '
           f'2>&1 | tee {out_path}counts_{out_file}.txt '
           )

    out = subprocess.check_output([cmd], shell=True)

    if 'Fatal exception' in out.decode('utf-8'):
        missing_screens.append(screen_name)
        os.remove(f'{out_path}counts_{out_file}.txt')
        os.remove(f'{out_path}out_{out_file}.txt')
    else:
        analyzed_screens.append(screen_name)

print(f'\n\nSuccessfully processed screens: {analyzed_screens}')
print(f'Warning: Missing screens: {missing_screens}')
