import numpy as np
import subprocess
import time
from teslaswarm.settings import BASE_DIR
from tools.data_foo import xyz_to_cast, xyz_to_spher, geo_to_mag_coord
import platform


def write_param(param_dict):
    keys = ['bz', 'by', 'doy', 'kp', 'f107', 'ut']
    format_param = [' Bz=    v', ' By=    v', ' DOY=    v', ' Kp=    v', ' F107=  v', ' UT=    v']
    f = open(BASE_DIR + '/E_FIELD_FAC_MODEL_FCP/input data', 'w')
    for i, p in enumerate(format_param):
        line = p.replace('v', str(np.around(float(param_dict[keys[i]]), 1)))
        f.write(line + '\n')
    f.write('   BT=   %s\n' % (np.sqrt(float(param_dict['bz'])**2 + float(param_dict['by'])**2).round(1)))
    f.write('   IMF clock angle=  45.0')
    f.close()

def run_exe_program():

    if platform.system() == 'Linux' or platform.system() == 'Darwin':
            out = subprocess.call(
                'cd E_FIELD_FAC_MODEL_FCP \nwine coef_ut.exe \nwine fac_bt_season.exe \nwine tph.exe \nwine ris_surfer_sigma.exe \nwine ris_surfer_fac.exe',
                shell=True)
    elif platform.system() == 'Windows':
        #path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '\\E_FIELD_FAC_MODEL_FCP')
        path = BASE_DIR + "/E_FIELD_FAC_MODEL_FCP"
        subprocess.call('.\\coef_ut.exe', shell=True, cwd=path)
        subprocess.call('.\\fac_bt_season.exe', shell=True, cwd=path)
        subprocess.call('.\\tph.exe', shell=True, cwd=path)
        subprocess.call('.\\ris_surfer_sigma.exe', shell=True, cwd=path)
        subprocess.call('.\\ris_surfer_fac.exe', shell=True, cwd=path)



def find_th_name(param_dict):
    type_hem = np.array(['n_fac', 's_fac', 'n_pot', 's_pot', 'n_sigh', 's_sigh'])
    map_types = ['n_fac_surf', 's_fac_surf', 'n_potential_surf', 's_potential_surf',
                 'n_sigma_surf', 's_sigma_surf']

    input_th = '%s_%s' % (param_dict['hem'], param_dict['type'])
    idx = int(np.where(type_hem == input_th)[0])
    return map_types[idx]


def array_to_string(xyz):
    str_array = ''
    for x, y, z, in xyz:
        s = '%.2f %.2f %.2f\n' % (x, y, z)
        str_array += s
    return str_array

def get_ionomodel_surf_file(param_dict):

    write_param(param_dict)     # dict to file
    run_exe_program()   # run wine
    time.sleep(3)
    file_name = find_th_name(param_dict)

    """получение иономодели из файла"""
    #s_fac_surf
    file = open(BASE_DIR + '/E_FIELD_FAC_MODEL_FCP/' + file_name + '.dat')
    xyz_array = np.empty((0, 3))
    for line in file.readlines():
        xyz_array = np.append(xyz_array, [np.array(line.split()).astype(float)], axis=0)
    X = xyz_array[:, 1]
    Y = xyz_array[:, 0] * -1
    Z = xyz_array[:, 2]

    if file_name[0] == 'n':
        PHI_P, THETA_P = 252.6, 80.4
    else:
        PHI_P, THETA_P = 107.4, -80.4

    desc_mag = xyz_to_cast(X, Y)
    desc_geo = geo_to_mag_coord(desc_mag, PHI_P, THETA_P)
    x, y = xyz_to_spher(desc_geo)

    # print('\nMEAN: %s\n--------------' % np.mean(Z))
    return np.array([x, y, Z]).T
