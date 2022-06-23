import pandas as pd
import numpy as np
from tools.dt_foo import *
from teslaswarm.settings import STATIC_OS_PATH, DATA_PATH
#https://github.com/spacecataz/supermag

from tools.supermag_api import SuperMAGGetIndices, supermag_testing, SuperMAGGetData, SuperMAGGetInventory

#(status,mydata2c) = SuperMAGGetIndices( 'pilipenko',  '2017-09-18',  3600, 'all,swiall,imfall', FORMAT='list')
#print(mydata2c[0:1])
with open(DATA_PATH + 'supermag-stations.txt') as f:
    for i, line in enumerate(f.readlines()):
        if i>43:
            #IAGA	GEOLON	GEOLAT	AACGMLON	AACGMLAT	STATION-NAME		OPERATOR-NUM	OPERATORS
            line = line.split()
            print(line)
#mydata2c = SuperMAGGetInventory('pilipenko',  '2017-09-18',  3600, )
#print(mydata2c)
