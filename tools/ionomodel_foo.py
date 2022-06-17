import numpy as np
import subprocess
import time
from teslaswarm.settings import BASE_DIR
from tools.data_foo import MAG2GEO, cart_to_cartGeomag_coord, sph_to_cart, cart_to_geo, get_lon_from_LT
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

def get_ionomodel_surf_file(param_dict, dt):
    # https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2005JA011465
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
    mglat = xyz_array[:, 0]
    mglon = xyz_array[:, 1]
    Z = xyz_array[:, 2]


    #GEO_latlon = np.array([MAG2GEO(geolat[i], geolon[i], dt) for i in range(len(geolat))])
    #lat = GEO_latlon[:, 1]
    #lon = GEO_latlon[:, 0]

    if file_name[0] == 'n':
        THETA_P, PHI_P  = 213.00, 85.54
    else:
        THETA_P, PHI_P = 223., -64.24
    lon_to360 = lambda x: (x - 180) % 360 - 180  # -180 180 to 0 360
    lon_to180 = lambda x: (x + 180) % 360 - 180  # 0 360 to -180 180

    desc_mag = sph_to_cart(mglon, mglat)
    desc_mag = cart_to_cartGeomag_coord(desc_mag, THETA_P, PHI_P)
    lat, lon = cart_to_geo(desc_mag)

    seconds_per_day = 24 * 60 * 60.0
    ut_sec = float(param_dict['ut'])*(60 * 60.0)
    midnight_lon = (ut_sec / seconds_per_day) * 360
    print(midnight_lon, 'midn lon')
    lon = lon_to360(lon) + midnight_lon
    for i in range(len(lon)):
        if lon[i] < 0:
            lon[i] = 360 + lon[i]
        elif lon[i] > 360:
            lon[i] = lon[i] - 360
    #lon = lon_to180(lon)
    print(lon.min(), lon.max(), np.mean(lon), 'X')
    print(lat.min(), lat.max(), np.mean(lat), 'Y')
    # print('\nMEAN: %s\n--------------' % np.mean(Z))
    return np.array([lon, lat, Z]).T
