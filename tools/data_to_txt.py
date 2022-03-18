import numpy as np
import pandas as pd
import itertools
from teslaswarm.settings import STATIC_OS_PATH
from tools.dt_foo import get_timestamp
import codecs

class TextData():
    def __init__(self, SWARM_liter):
        self.DATA = {}

        self.SWARM_liter = SWARM_liter
        self.SWARM_legend = ['lat', 'lon', 'r']
        self.field_legend = ['lat', 'lon', 'value']

        self.SWARM_channel = None
        self.ionomodel_shp = None
        self.auroral_shp = None
        self.annotate = None

        self.id = get_timestamp()

    def append(self, data, name):
#        print(data, name)

        self.DATA[name] = np.array(data)

    def save_proj_data_to_txt(self):
        #print(self.DATA.keys())
        columns = []
        columns_legend = {}
        #print(self.DATA)
        for k in self.DATA.keys():
            """
            каждый тип добавленных данных
            """
            array = np.array(self.DATA[k])
            try:
                _, lcol = np.shape(array)
            except:
                lcol = 1
            for c in range(lcol):
                if k == 'SWARM_pos':
                    legend_key = 'SWARM-%s_position' % self.SWARM_liter
                    if c == 0:
                        columns_legend[legend_key] = []
                    columns_legend[legend_key].append(self.SWARM_legend[c])
                if k == 'SWARM':
                    legend_key = 'SWARM-%s_value' % self.SWARM_liter
                    if c == 0:
                        columns_legend[legend_key] = []
                    if self.SWARM_channel != 'FAC2':
                        columns_legend[legend_key].append(self.SWARM_channel + ',nT')
                    else:
                        columns_legend[legend_key].append('FAC2,nA/m^2')
                if k == 'SWARM_IGRF':
                    legend_key = 'SWARM-%s_IGRF-13_diff' % self.SWARM_liter
                    if c == 0:
                        columns_legend[legend_key] = []
                    columns_legend[legend_key].append('Dd,nT')
                if k == 'measure_mu':
                    legend_key = 'SWARM-%s_anomalyLVL' % self.SWARM_liter
                    if c == 0:
                        columns_legend[legend_key] = []
                    if c != 2:
                        columns_legend[legend_key].append(self.field_legend[c])
                    else:
                        columns_legend[legend_key].append('anomaly_lvl')
                if k == 'ionomodel':
                    legend_key = 'ionomodel_%s' % self.ionomodel_shp
                    if c == 0:
                        columns_legend[legend_key] = []
                    if c != 2:
                        columns_legend[legend_key].append(self.field_legend[c])
                    else:
                        columns_legend[legend_key].append('uA/m^2')
                if k == 'auroral':
                    legend_key = 'auroral_oval_%s' % self.auroral_shp
                    if c == 0:
                        columns_legend[legend_key] = []
                    if c != 2:
                        columns_legend[legend_key].append(self.field_legend[c])
                    else:
                        columns_legend[legend_key].append('ergs/cm^2')

                if lcol == 1:
                    columns.append(array.flatten())
                else:
                    columns.append(array[:, c])
        head_legend = []
        c_legend = []
        for g in columns_legend.keys():
            gg = str(g) + ' ' * (len(columns_legend[g]) - 1)
            head_legend.append(gg)
            c_legend.extend(columns_legend[g])

        f = codecs.open(STATIC_OS_PATH + '/media/txt/%s.txt' % self.id, "w", "utf-8")
        f.writelines(self.annotate + '\n')
        f.writelines(' '.join(np.array(head_legend).astype(str))+'\n')
        f.writelines(' '.join(np.array(c_legend).astype(str))+'\n')
        fill_columns = list(itertools.zip_longest(*columns, fillvalue=''))
        for line in fill_columns:
            #print(tuple(line))
            f.writelines(' '.join(np.array(line).astype(str))+'\n')
        f.close()
        print('txt file saved to %s' % (STATIC_OS_PATH + '/media/txt/%s.txt' % self.id))



