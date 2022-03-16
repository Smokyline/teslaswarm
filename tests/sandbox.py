import pandas as pd
import numpy as np
from tools.dt_foo import *
from teslaswarm.settings import STATIC_OS_PATH
#https://github.com/spacecataz/supermag

with open(STATIC_OS_PATH + '/media/txt/test.txt', "r",) as f:
    for line in f.readlines():
        print(line)



