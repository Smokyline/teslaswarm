import numpy as np
import pandas as pd
import itertools
from teslaswarm.settings import STATIC_OS_PATH
from tools.dt_foo import get_timestamp
import codecs

class Data2Text():
    def __init__(self, SWARM_liter):
        self.DATA = {}

        self.SWARM_liter = SWARM_liter
        self.SWARM_legend = ['lat', 'lon', 'r']
        self.field_legend = ['lat', 'lon', 'value']

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

    def array_to_column(self):
        #print(self.DATA.keys())

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
                    legend_key = 'SWARM-%s_GEOposition' % self.SWARM_liter
                    if c == 0:
                        self.columns_legend[legend_key] = []
                    self.columns_legend[legend_key].append(self.SWARM_legend[c])
                if k == 'SWARM':
                    legend_key = 'SWARM-%s_value' % self.SWARM_liter
                    if c == 0:
                        self.columns_legend[legend_key] = []
                    if self.SWARM_channel != 'FAC':
                        self.columns_legend[legend_key].append(self.SWARM_channel + ',nT')
                    else:
                        self.columns_legend[legend_key].append('FAC,nA/m^2')
                if k == 'SWARM_IGRF' :
                    legend_key = 'SWARM-%s_IGRF_diff' % self.SWARM_liter
                    if c == 0:
                        self.columns_legend[legend_key] = []
                    self.columns_legend[legend_key].append('dB,nT')
                if k == 'SWARM_CHAOS7':
                    legend_key = 'SWARM-%s_CHAOS7_diff' % self.SWARM_liter
                    if c == 0:
                        self.columns_legend[legend_key] = []
                    self.columns_legend[legend_key].append('dB,nT')
                if k == 'measure_mu':
                    legend_key = 'SWARM-%s_anomalyLVL' % self.SWARM_liter
                    if c == 0:
                        self.columns_legend[legend_key] = []
                    if c != 2:
                        self.columns_legend[legend_key].append(self.field_legend[c])
                    else:
                        self.columns_legend[legend_key].append('anomaly_lvl')
                if 'ionomodel' in k:
                    legend_key = 'ionomodel_%s' % self.ionomodel_shp
                    if c == 0:
                        self.columns_legend[legend_key] = []
                    if c != 2:
                        self.columns_legend[legend_key].append(self.field_legend[c])
                    else:
                        if 'sigh' in k:
                            self.columns_legend[legend_key].append('Sm')
                        elif 'pot' in k:
                            self.columns_legend[legend_key].append('kV')
                        else:
                            self.columns_legend[legend_key].append('uA/m^2')
                if 'auroral' in k:
                    #print(k, array)
                    legend_key = 'auroral_oval_%s' % self.auroral_shp
                    if c == 0:
                        self.columns_legend[legend_key] = []
                    elif c == 0 and self.auroral_shp == 'near':
                        self.columns_legend[legend_key].append(['ergs/cm^2'])
                    if c != 2:
                        self.columns_legend[legend_key].append(self.field_legend[c])
                    else:
                        self.columns_legend[legend_key].append('ergs/cm^2')

                if lcol == 1:
                    self.columns.append(array.flatten())
                    #print(array.flatten())
                else:
                    self.columns.append(array[:, c])
                    #print(array[:, c])


    def save_columns_to_txt(self, ):
        head_legend = []
        c_legend = []
        for g in self.columns_legend.keys():
            gg = str(g) + ' ' * (len(self.columns_legend[g]) - 1)
            head_legend.append(gg)
            c_legend.extend(self.columns_legend[g])

        f = codecs.open(STATIC_OS_PATH + '/media/txt/%s.txt' % self.id, "w", "utf-8")
        f.writelines(self.annotate + '\n')
        f.writelines(' '.join(np.array(head_legend).astype(str)) + '\n')
        f.writelines(' '.join(np.array(c_legend).astype(str)) + '\n')
        fill_columns = list(itertools.zip_longest(*self.columns, fillvalue=''))
        for line in fill_columns:
            # print(tuple(line))
            f.writelines(' '.join(np.array(line).astype(str)) + '\n')
        f.close()
        print('txt file saved to %s' % (STATIC_OS_PATH + '/media/txt/%s.txt' % self.id))



