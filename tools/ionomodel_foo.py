import datetime

import numpy as np
import subprocess
import time
from teslaswarm.settings import BASE_DIR
from tools.data_foo import *
import platform
import math
from spacepy import coordinates as coord
from spacepy.time import Ticktock
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


def to_cast(sph_x, sph_y):
    r = 6371000
    desc_data = np.empty((0, 3))
    for i in range(len(sph_x)):

        theta, phi = math.radians(sph_x[i]), math.radians(sph_y[i])
        x = r * np.sin(phi) * np.cos(theta)
        y = r * np.sin(phi) * np.sin(theta)
        z = r * np.cos(phi)
        desc_data = np.append(desc_data, [[x, y, z]], axis=0)
    return desc_data


"""def to_spher(data):
    spher_data = np.empty((0, 2))
    for x, y, z in data:
        hxy = np.hypot(x, y)
        # r = np.hypot(hxy, z)
        el = np.arctan2(z, hxy)
        az = np.arctan2(y, x)
        spher_data = np.append(spher_data, [[math.degrees(az), math.degrees(el)]], axis=0)
    return spher_data[:, 0], spher_data[:, 1]"""


def to_spher(data):
    spher_data = np.empty((0, 2))
    for x, y, z in data:
        radius = np.sqrt(x ** 2 + y ** 2 + z ** 2)
        latitude = np.arcsin(z / radius)
        longitude = np.arctan2(y, x)
        spher_data = np.append(spher_data, [[math.degrees(latitude), math.degrees(longitude)]], axis=0)
    return spher_data[:, 0], spher_data[:, 1]   # lat, lon

def mag_coord_to_geo(mag_xyz, phi, theta):
    PHI_P, THETA_P = math.radians(phi), math.radians(theta)

    R = np.array([
        [np.cos(PHI_P) * np.cos(THETA_P), -1 * np.sin(PHI_P), np.cos(PHI_P) * np.sin(THETA_P)],
        [np.sin(PHI_P) * np.cos(THETA_P), np.cos(PHI_P), np.sin(PHI_P) * np.sin(THETA_P)],
        [-1 * np.sin(THETA_P), 0, np.cos(THETA_P)]])

    geo_data = np.empty((0, 3))
    for xyz in mag_xyz:
        geo_xyz = np.sum(np.array([xyz]).T * R.T, axis=0)
        geo_data = np.append(geo_data, [geo_xyz], axis=0)
    return geo_data

def mag_crt_to_geo_shp(r, lat, lon, dt):
    cvals = coord.Coords([np.float(r), np.float(lat), np.float(lon)], 'MAG', 'car', ['Re', 'deg', 'deg'])
    cvals.ticks = Ticktock(dt)
    new_coord = np.array(cvals.convert('GEO', 'sph').data)[0]
    return [new_coord[2], new_coord[1]]  # mag lon, lat


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
    x = xyz_array[:, 0]
    y = xyz_array[:, 1]
    value = xyz_array[:, 2]

    if file_name[0] == 'n':
        THETA_P, PHI_P  = 213.00, 85.54
    else:
        THETA_P, PHI_P = 223., -64.24
    lon_to360 = lambda x: (x - 180) % 360 - 180  # -180 180 to 0 360
    lon_to180 = lambda x: (x + 180) % 360 - 180  # 0 360 to -180 180

    pollonlat = np.array([cart2polar(x, y) for x, y in zip(x, y)])
    lat, lon = pollonlat[:, 0], pollonlat[:, 1]
    lat = 90. - lat
    geo_latlon = geomag2geo(latlon=np.array([lat, lon]).T, THETA=THETA_P, PHI=PHI_P)
    lat, lon = geo_latlon[:, 0], geo_latlon[:, 1]

    #seconds_per_day = 24 * 60 * 60.0
    #ut_sec = float(param_dict['ut'])*(60 * 60.0)
    ut_h = int(float(param_dict['ut']))
    if float(param_dict['ut']) == 0.:
        ut_m = int((float(param_dict['ut'])*10)*60)
    else:
        ut_m = int((float(param_dict['ut'])%ut_h)*60)
    dt_model = datetime.datetime(dt.year, dt.month, dt.day, ut_h, ut_m)
    midnight_lat, midnight_lon = sun_pos(dt_model)

    print(midnight_lon, 'midn lon')
    lon = midnight_lon + lon_to360(lon)
    for i in range(len(lon)):
        if lon[i] < 0:
            lon[i] = 360 + lon[i]
        elif lon[i] > 360:
            lon[i] = lon[i] - 360
    return np.array([lon, lat, value]).T
