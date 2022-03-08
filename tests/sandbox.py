import pandas as pd
import numpy as np
from tools.dt_foo import *
#https://github.com/spacecataz/supermag


def open_superMAG_vcf(date):

    filename = date.split('-')
    filename = '%s%s%s' % (filename[0][2:], filename[1], filename[2])
    with open('C:\\Users\\ivan\\workspace\\teslaswarm\\static\\data\\obs_data\\SALU\\%s.T47'%filename, mode='r') as vcf:
        lines = vcf.readlines()

    obs_data = []
    zero_time = decode_str_dt_param(date+'T'+'00:00:00')
    for i, line in enumerate(lines):
        if i == 0:
            continue
        line = np.array(line.rstrip('\n').split(',')[0].split('    ')[1:]).astype(float)    # ['-111.600', ' -13.900', '-158.200']
        zero_time = zero_time + pd.Timedelta('%s minutes' % 1)
        obs_data.append([zero_time, line[0], line[1], line[2]])
    return obs_data

open_superMAG_vcf('2017-09-18')