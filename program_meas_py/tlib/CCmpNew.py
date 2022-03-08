import numpy as np

from program_meas_py.tlib.MTools import get_limits
from program_meas_py.tlib.CmpForm import CmpForm
from program_meas_py.tlib.MerasAlpha import MerasAlpha

class CCmpNew:
    
    kol_level = 0                  #    кол-во слоев
    dlen = 0                       #    длина массива данных
    gamma = []
    
    kol_point_plus = 0   
    hrazb = 0                       #   шаг разбиения
    razbi = []                      #   разбиение
    num_point_level = []            #   кол-во точек в слое
    number_level = []               #   номер слоя, к которому принадлежит точка
    data = []                       #   исходный массив данных
    dsize = []                      #   размерность исходного массива
    level = [] 
    value = []
    bound = []
#    hfunc1 = CmpForm
#    hfunc2 = MerasAlpha

    def __init__(self, *args):
        self.kol_level = args[0]
        self.data = args[1]
        if len(args) == 3:
            if args[2]:
                self.sfunc2 = MerasAlpha

        self.dsize = self.data.shape                         #  размерности исходного массива
        self.dlen = self.data.size
        self.data = self.data.reshape((-1,))
        ind = np.where(self.data == 0)[0]
        self.data[ind] = 10e-10
            
        amin = self.data.min()
        amax = self.data.max()
        self.hrazb = (amax - amin)/(self.kol_level - 1)
    
        self.number_level = np.array(np.round((self.data - amin)/self.hrazb) + 1).astype(int)
        self.razbi = np.linspace(amin, amax, self.kol_level)
            
        self.num_point_level = np.zeros((self.kol_level,)).astype(int)

        self.hfunc1(self.number_level, self.num_point_level)

    def hfunc1(self, number_level, num_point_level):
        self.num_point_level = CmpForm(number_level, num_point_level)

    def hfunc2(self, dlen, num_point_level, razbi1, razbi2, gamma):
        level = MerasAlpha(dlen, num_point_level, razbi1, razbi2, gamma)
        return level
                        
    def build_level(self, gamma):
            
        self.gamma = gamma
        self.level = self.hfunc2(self.dlen, self.num_point_level, self.razbi, self.razbi, self.gamma)

    def build_value(self, gamma, level):
                        
        valuel = np.zeros((len(level), len(gamma)))
            
        for ig in range(len(gamma)):
            for il in range(len(level)):
                valuel[il, ig] = get_level_one(self.hfunc2, self.dlen,
                                               self.num_point_level, self.razbi, gamma, level)

        self.gamma = gamma
        self.level = level
        self.value = valuel
       
    def get_value1(self, gamma, level):
                           
        return get_level_one(self.hfunc2, self.dlen, self.num_point_level, self.razbi, gamma, level)

    def get_value2(self, igamma, ilevel):
                            
        return self.value[ilevel, igamma]
       
    def get_level1(self, gamma, value):
            
        return self.hfunc2(self.dlen, self.num_point_level, self.razbi, value, gamma)
        
    def get_level1m(self, gamma):
            
        ret1 = self.hfunc2(self.dlen, self.num_point_level, self.razbi, self.razbi, gamma)
            
        ret = np.zeros((self.dlen, ))
        indw = np.where(self.num_point_level>0)[0]
    
        for i in range(len(indw)):
            j = indw[i]
            ind = np.where(self.number_level == j + 1)[0]
#            print(ind)
            ret[ind]= ret1[j]

        return ret
       
    def get_level2(self, gamma, value):
                            
        ret = -1
            
        ig = np.where(np.abs(self.gamma-gamma) < 0.00001)[0]
        if len(ig) == 0:
            return ret
        ig = ig[0]
            
        if value < self.razbi[0, ig]:
            ret = -1
        elif value == self.razbi[0, ig]:
            ret = self.level[0]
        elif value == self.razbi[-1, ig]:
            ret = self.level[-1]
        elif value > self.razbi[-1, ig]:
            ret = 1
        else:
#                 i = find((self.razbi(1:end-1, ig) < value) & (value <= self.razbi(2:end, ig)))
#                 skf = (value - self.razbi(i, ig))/(self.razbi(i+1, ig) - self.razbi(i, ig))
#                 ret = self.level(i) + skf*(self.level(i+1) - self.level(i))
            i = np.argmin(np.abs(self.razbi[:, ig] - value))
            ret = self.level[i]
                        
        return ret
        
    def build_bound(self):
            
        self.bound = [[[] for _ in range(len(self.gamma))] for _ in range(len(self.level))]
            
        for il in range(len(self.level)):
            for ig in range(len(self.gamma)):
                    
                valuet = self.value[il, ig]
                be = get_limits(np.where(self.data > valuet)[0])
                self.bound[il][ig] = be
       
    def get_bound1(self, *args):
                            
        ret = []
            
        if len(args) != 2:
            return ret

        gamma, level = args
            
        ig = np.where(np.abs(self.gamma-gamma) < 0.00001)[0]
        if len(ig) == 0:
            return ret
           
        il = np.where(np.abs(self.level-level) < 0.00001)[0]
        if len(il) == 0:
            return ret

        return self.bound[il[0]][ig[0]]
       
    def get_bound2(self, igamma, ilevel):
            
        return self.bound[ilevel][igamma]
        
    def get_data(self, *args):
            
        if len(args) == 0:
            ret = self.data
        else:
            ret = self.data(args[0])

        return ret

def get_level_one(hfunc, kol_point, num_point_level, razbi, gamma, level):  #   rlevel
            
    rmin = razbi.min()
    lmin = hfunc(kol_point, num_point_level, razbi, rmin, gamma)
    while lmin > level:
        rmin = rmin/2
        lmin = hfunc(kol_point, num_point_level, razbi, rmin, gamma)
        if rmin <10e-10:
            rlevel = 0
            return rlevel
            
    rmax = razbi.max()
    lmax = hfunc(kol_point, num_point_level, razbi, rmax, gamma)
    while lmax < level:
        rmax = 2*rmax
        lmax = hfunc(kol_point, num_point_level, razbi, rmax, gamma)
            
    rtmp = (rmax+rmin)/2
    ltmp = hfunc(kol_point, num_point_level, razbi, rtmp, gamma)
            
    while abs(ltmp - level) > 0.00001:
        if ltmp < level:
            rmin = rtmp
        elif ltmp > level:
            rmax = rtmp
        else:
            break

        rtmp = (rmax+rmin)/2
        ltmp = hfunc(kol_point, num_point_level, razbi, rtmp, gamma)
            
    rlevel  = rtmp

    return rlevel