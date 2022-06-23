import numpy as np
import pandas as pd
import itertools
from teslaswarm.settings import STATIC_OS_PATH
from tools.dt_foo import get_timestamp
import codecs

class Data2Text():
    def __init__(self):
        self.DATA = {}

        self.SWARM_channel = None
        self.ionomodel_shp = None
        self.auroral_shp = None
        self.annotate = None

        self.columns = []
        self.columns_legend = {}

        self.id = get_timestamp()

    def append(self, data, name):
       # print(data, name)
        self.DATA[name] = np.array(data)

    def save_columns_to_txt(self, ):
        head_legend = []
        #c_legend = []
        columns = []
        for g in self.DATA.keys():
            if type(self.DATA[g]) == dict:
                #len_h = (len(g) - 1)
                for sub_k in self.DATA[g]:
                    gg = str(g)+'_'+str(sub_k)+''
                    head_legend.append(gg)
                    columns.append(self.DATA[g][sub_k])
            else:
                gg = str(g) + ' '
                head_legend.append(gg)
                columns.append(self.DATA[g])

        f = codecs.open(STATIC_OS_PATH + '/media/txt/%s.txt' % self.id, "w", "utf-8")
        f.writelines(self.annotate + '\n')
        f.writelines(' '.join(np.array(head_legend).astype(str)) + '\n')
        #f.writelines(' '.join(np.array(c_legend).astype(str)) + '\n')
        fill_columns = list(itertools.zip_longest(*columns, fillvalue=''))
        for line in fill_columns:
            # print(tuple(line))
            f.writelines(' '.join(np.array(line).astype(str)) + '\n')
        f.close()
        print('txt file saved to %s' % (STATIC_OS_PATH + '/media/txt/%s.txt' % self.id))



