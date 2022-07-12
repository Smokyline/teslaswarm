import numpy as np
from teslaswarm.settings import BASE_DIR, STATIC_OS_PATH, DATA_PATH
from program_meas_py.meas_extr import get_measure_mu
import math
from numpy import degrees, radians
import pandas as pd
import shapefile
import subprocess
import time
from ovationpyme import ovation_prime, ovation_utilities
from geospacepy import satplottools, special_datetime
import igrf
from tools.dt_foo import *
from tools.coordinates_convert import geomag2geo, geo2GSM
from tools.sun_position import sun_pos
from tools.supermag_api import SuperMAGGetData, sm_grabme, SuperMAGGetIndices
from scipy.interpolate import griddata
import shapely
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
import json
import geog

import warnings
from pyproj import Proj, transform, CRS



def data_reduction(respond, delta, fac2_mod=False):
    warnings.simplefilter("ignore", category=RuntimeWarning)
    """сжимание секундных данных до delta шага"""
    # fac2 = (Y, X, R), dt, (fac2)
    # vector, measure mu, fac, chaos = (Y, X, R), dt, (N, E, C)
    N, M = respond.shape
    if fac2_mod:
        idx999 = np.where(respond[:, 4] >= 999)[0]
        # respond = respond[idx999]
        respond[idx999, 4] = np.nan
        redu_resp = np.empty((0, 5))
    else:
        redu_resp = np.empty((0, M))

    if delta <= 1:
        window = 1
    else:
        window = int(delta / 2)
    st_idx = 0
    while st_idx < N:
        if st_idx != 0:
            left_idx = st_idx - window
        else:
            left_idx = 0
        right_idx = st_idx + window
        delta_resp = respond[left_idx:right_idx]
        if fac2_mod:
            dt, y, x, r, fac2 = delta_resp[-1, (0, 1, 2, 3, 4)]
            fac2 = np.nanmean(delta_resp[:, 4])
            #if np.isnan(fac2):
            #    fac2 = 0.
            redu_resp = np.append(redu_resp, [[dt, y, x, r, fac2]], axis=0)

        else:
            dt, y, x, r = delta_resp[-1, (0, 1, 2, 3)]
            n = np.mean(delta_resp[:, 4])
            e = np.mean(delta_resp[:, 5])
            c = np.mean(delta_resp[:, 6])
            f = np.mean(delta_resp[:, 7])
            redu_resp = np.append(redu_resp, [[dt, y, x, r, n, e, c, f]], axis=0)
        st_idx += delta

        """if fac2_mod:
        #https://pandas.pydata.org/docs/user_guide/missing_data.html
        redu_resp[:, 4] = pd.DataFrame(np.array(redu_resp[:, 4]), dtype='float32').interpolate().to_numpy().T[0]
        """

        """
        miss_values_pd = pd.DataFrame(redu_resp[:, 4])
        fac2_miss_values_pd = miss_values_pd.fillna(miss_values_pd.mean())
        #fac2_miss_values_pd = miss_values_pd.fillna(value=miss_values_pd)
        redu_resp[:, 4] = fac2_miss_values_pd.T.to_numpy()
        """

            #for i, resp in enumerate(fac2_miss_values_pd):
             #       if np.isnan(respond[i, 3]):
             #                   print(respond[i], 'is nan')
    return redu_resp



def calc_ACvector(sw_a_cd, sw_c_cd, sw_a_values, sw_c_values, channel):
    """вычисление разницы векторов между Swarm A и Swarm C"""
    coord = np.empty((0, 3))
    values = np.empty((0, 1))

    for cd_a, cd_c, av, cv in zip(sw_a_cd, sw_c_cd, sw_a_values, sw_c_values):
        #x = np.mean([cd_a[1], cd_c[1]])
        #y = np.mean([cd_a[0], cd_c[0]])
        #r = np.mean([cd_a[2], cd_c[2]])
        x, y, r = cd_a[1], cd_a[0], cd_a[2]
        coord = np.append(coord, [[y, x, r]], axis=0)

        if channel != None:
            v = np.sqrt(np.power(av[channel] - cv[channel], 2))
        else:
            v = np.sqrt(np.power(av - cv, 2))
        values = np.append(values, v)

    return coord, np.array(values)


# shapely polygons
def points_in_poly(poly_points, swarm_pos, shapely_convert):
    if shapely_convert:
        polygon = Polygon([(lon, lat) for lon, lat in poly_points])
    else:
        polygon = poly_points
    p_containts = [polygon.contains(Point(lat, lon)) for lon, lat in swarm_pos]
    return p_containts, polygon

def get_points_in_poly(swarm_pos, swarm_poly_loc, proj_type, shapely_convert_to_poly=True):
    if proj_type == 'ortho_n' or proj_type == 'ortho_s':
        p_in_p, poly = points_in_poly(poly_points=swarm_poly_loc, swarm_pos=swarm_pos, shapely_convert=shapely_convert_to_poly)
    elif proj_type == 'miller':
        # idk but work
        p_in_p1, poly = points_in_poly(poly_points=swarm_poly_loc, swarm_pos=swarm_pos, shapely_convert=shapely_convert_to_poly)
        p_in_p2, poly = points_in_poly(poly_points=swarm_poly_loc, swarm_pos=swarm_pos, shapely_convert=shapely_convert_to_poly)
        for i, point in enumerate(p_in_p2):
            if p_in_p1[i] == False:
                p_in_p1[i] = p_in_p2[i]
        p_in_p = p_in_p1
    return p_in_p, poly   # bool array


def get_poly_around_point(xy, n_points, radius):
    print('create poly around point...')
    p = shapely.geometry.Point([xy[0], xy[1]])
    #print(p)
    d = radius * 1000  # meters
    angles = np.linspace(0, 360, n_points)
    polygon = geog.propagate(p, angles, d)
    print(json.dumps(shapely.geometry.mapping(shapely.geometry.Polygon(polygon))))
    return shapely.geometry.Polygon(polygon)

def get_swarm_value_near_obs(swarm_pos, code, proj_type):
    # https://gis.stackexchange.com/questions/268250/generating-polygon-representing-rough-100km-circle-around-latitude-longitude-poi
    obs_location = get_INTERMAGNET_observ_loc(code)[:2]
    poly = get_poly_around_point(obs_location, n_points=35, radius=1000)
    p_in_p = get_points_in_poly(swarm_pos, poly, proj_type, shapely_convert_to_poly=False)
    return p_in_p  # bool array


def get_position_near_point(swarm_pos, obs_location, degr_radius):
    km = 111.134861111
    eucl = eucl_distance(obs_location, swarm_pos)
    near_swarm_pos = np.full(len(swarm_pos), False)
    near_swarm_pos[np.where(eucl < degr_radius)] = True
    #print('near pos|', len(near_swarm_pos), 'where True:', len(np.where(near_swarm_pos == True)[0]))

    return near_swarm_pos

def get_swarm_poly_loc(point, deg_radius):
    # poly = [-x, +x, -y, +y]
    #delta = deg_radius * 2
    x1, x2 = point[0] - deg_radius * 2, point[0] + deg_radius * 2
    y1, y2 = point[1] - deg_radius * 2, point[1] + deg_radius * 2
    return [x1, x2, y1, y2]

#############


def get_mag_coordinate_system_lines(date, geomag_pole):
    # geomag == True -> pole is Geomagnetic, else pole is geocentric Magnetospheric
    lat_coord = [-90.0, 91.]
    lon_coord = [-180.0, 181.]
    lat_step = 15   # step of legend
    lon_step = 15   # step of legend
    r = 6371.2
    mag_lat_lines = []
    mag_lon_lines = []
    annotate_points = []
    for lat_lines_start in np.arange(lat_coord[0], lat_coord[1], lat_step):
        lat_lines = []  # coord lines in geo
        for lon in np.arange(-180, 181, 1):
            lat_lines.append([lat_lines_start, lon, r])
        #   convert to mag\geomag coords
        if geomag_pole:
            #mag_lat_lines.append(apex_coord_system_convert(np.array(lat_lines)[:, :2], source_sys='geo', dest_sys='apex'))
            #mag_lat_lines.append(geo2geomag(np.array(lat_lines)[:, :2]))
            mag_lat_lines.append(geo2GSM(np.array(lat_lines), np.full((1, len(lat_lines)), date)[0], to_coord='MAG'))

        else:
            mag_lat_lines.append(geo2GSM(np.array(lat_lines), np.full((1, len(lat_lines)), date)[0], to_coord='GSM'))


    for lon_lines_start in np.arange(lon_coord[0], lon_coord[1], lon_step):
        lon_lines = []
        for lat in np.arange(-90, 91, 1):
            lon_lines.append([lat, lon_lines_start, r])

        #   convert to mag\geomag coords
        if geomag_pole:
            #mag_lon_lines.append(apex_coord_system_convert(np.array(lon_lines)[:, :2], source_sys='geo', dest_sys='apex'))
            #mag_lat_lines.append(geo2geomag(np.array(lon_lines)[:, :2]))
            mag_lon_lines.append(geo2GSM(np.array(lon_lines), np.full((1, len(lon_lines)), date)[0], to_coord='MAG'))

        else:
            mag_lon_lines.append(geo2GSM(np.array(lon_lines), np.full((1, len(lon_lines)), date)[0], to_coord='GSM'))

    for lat in np.arange(-90+lat_step, 91-lat_step, lat_step):
        for lon in np.arange(-180, 181., lon_step):
            if geomag_pole:
                ap = geo2GSM(np.array([[lat, lon, r]]), np.full((1, 1), date)[0], to_coord='MAG')[0]
            else:
                ap = geo2GSM(np.array([[lat, lon, r]]), np.full((1, 1), date)[0], to_coord='GSM')[0]

            annotate_points.append([ap[0], ap[1], [lat, lon]])
    for lat in [-90, 90]:
        if geomag_pole:
            ap = geo2GSM(np.array([[lat, 0, r]]), np.full((1, 1), date)[0], to_coord='MAG')[0]
        else:
            ap = geo2GSM(np.array([[lat, 0, r]]), np.full((1, 1), date)[0], to_coord='GSM')[0]
        annotate_points.append([ap[0], ap[1], [lat, 0]])



    return mag_lat_lines, mag_lon_lines, annotate_points

def get_igrf_value(date, coords):
    """
     Parameters
     ----------

     date: datetime.date or decimal year yyyy.dddd
     glat, glon: geographic Latitude, Longitude
     alt_km: altitude [km] above sea level for itype==1
     isv: 0 for main geomagnetic field
     itype: 1: altitude is above sea level
     """
    igrf_model = np.empty((0, 4))
    for i in range(len(coords)):
        glat, glon, alt_km = coords[i, 0], coords[i, 1], coords[i, 2]
        alt_km = alt_km - earth_radius_in_meters(glat)/1000
        model = igrf.igrf(date[i], glat, glon, alt_km)
        igrf_model = np.append(igrf_model, [[float(model.north), float(model.east), float(model.down),
                                             float(model.total)]], axis=0)
    return igrf_model

def swarm_egrf_vector_subtraction(swarm_pos, swarm_values_full, swarm_date):
    #swarm_pos[:, 2] = swarm_pos[:, 2]
    igrf_value = get_igrf_value(swarm_date, swarm_pos)
    sw_N, sw_E, sw_C, sw_F = swarm_values_full[:, 0], swarm_values_full[:, 1], swarm_values_full[:, 2], swarm_values_full[:, 3],
    dN = sw_N - igrf_value[:, 0]
    dE = sw_E - igrf_value[:, 1]
    dC = sw_C - igrf_value[:, 2]
    dF = sw_F - igrf_value[:, 3]
    sat_H = np.sqrt(sw_F ** 2 - sw_C ** 2)
    model_H = np.sqrt(igrf_value[:, 3] ** 2 - igrf_value[:, 2] ** 2)
    dH = sat_H - model_H

    #print(np.array([sw_N, sw_E, sw_C, sw_F]).T, 'swarm')
    #print(np.array([igrf_value[:, 0],igrf_value[:, 1], igrf_value[:, 2], igrf_value[:, 3]]).T, 'chaos')
    #print(sat_H - model_H, np.min(sat_H - model_H), np.max(sat_H - model_H), 'dH')
    return np.array([dN, dE, dC, dH]).T

def magfield_variation(n_swarm, e_swarm, x_igrf, y_igrf, ):
    """
         n same as x
         e same as y
    """
    FACT = 180. / np.pi
    x, y,  = (n_swarm + x_igrf) / 2, (e_swarm + y_igrf) / 2,
    dx, dy, = (x_igrf - n_swarm) / 2, (y_igrf - e_swarm) / 2,
    h = np.sqrt(x * x + y * y)
    dd = (FACT * (x * dy - y * dx)) / (h * h)
    return dd, dx, dy

######################################################################################################

"""def calculate_magnetic_field_intensity(lon, lat, alt, year):
    #   https://pypi.org/project/pyIGRF/
    igrf_value = pyIGRF.igrf_value(lat, lon, alt, year)
    return igrf_value"""

def get_INTERMAGNET_observ_loc(codes):
    """получение локации обсерватории_code"""
    observ_data = pd.read_csv(BASE_DIR + '/static/data/observ_INT2.txt', sep='\t', decimal=',', comment='#',
                              encoding='unicode_escape')


    # observ_data.columns = ["LAT", "LON", "NAME", "INSTITUTION",	"FULL_TITLE", "NAME_ENG", "INST_ENG", "GIN", "COMMENT"]
    observ_data = observ_data[['LAT', 'LON', 'CODE']]

    observ_data = observ_data.reindex(columns=['LON', 'LAT', 'CODE'])

    x, y, code = observ_data.loc[observ_data['CODE'].isin([codes])].to_numpy()[0]
    x, y = (x + 180) % 360 - 180, y  # 360 to 180 -180
    # observ_loc = observ_raw[['LAT', 'LON']].to_numpy(dtype=float)
    #x, y = apex_convert(lat=observ_raw[1], lon=observ_raw[0], source='geo', dest='mlt', height=0, datetime=None)

    return [x, y, code]    # x y liter
    #return observ_raw  # x y liter

def get_superMAG_observ_loc(code):
    lat, lon = 0, 0
    with open(DATA_PATH + 'supermag-stations.txt') as f:
        for i, line in enumerate(f.readlines()):
            if i > 43:
                # IAGA	GEOLON	GEOLAT	AACGMLON	AACGMLAT	STATION-NAME		OPERATOR-NUM	OPERATORS
                line = line.split()
                if line[0] == code:
                    lat = float(line[2])
                    lon = (float(line[1]) + 180) % 360 - 180
    return [lon, lat, code]

def get_superMAG_value_from_web(date, station):
    #print(date, station)
    answer = SuperMAGGetData(logon='pilipenko',start='%sT00:00:00'%date, extent=86400,station=station, flagstring='',)
    #answer = SuperMAGGetData(logon='pilipenko',start='%sT00:00:00'%date, 86400, 'all', station, FORMAT = 'list')

    print(answer)
    station_data = []
    i, zero_time = 0, decode_str_dt_param(date + 'T' + '00:00:00')
    #https://supermag.jhuapl.edu/mag/?fidelity=low&start=2017-01-01T00%3A00%3A00.000Z&interval=23%3A59&tab=api#pythonClientDocumentationSection
    #for n, e, c in zip(sm_grabme(answer[1], 'N', 'nez'), sm_grabme(answer[1], 'E', 'nez'), sm_grabme(answer[1], 'Z', 'nez')):
    for n, e, c in zip(sm_grabme(answer[1], 'N', 'geo'), sm_grabme(answer[1], 'E', 'geo'), sm_grabme(answer[1], 'Z', 'geo')):
        station_data.append([zero_time, n, e, c])
        zero_time = zero_time + pd.Timedelta('%s minutes' % 1)
        i += 1
    return np.array(station_data)

def get_sMAGstation_value_by_time(date_array, time_array, delta, station):


    def fill_day_value(this_day_time_array, station_times, station_values):
        station_sec_values = []
        for sec in this_day_time_array:
            this_sec_in_station_times = np.where(station_times == sec)[0]
            if len(this_sec_in_station_times) == 0:
                station_sec_values.append(np.nan)
            else:
                station_sec_values.append(station_values[this_sec_in_station_times])
        df = pd.DataFrame({'value': station_sec_values}, dtype=float)
        df = df.interpolate('linear')
        station_sec_values = np.array([v for v in df.to_numpy()]).astype(float)
        return np.ravel(station_sec_values)
    ##########
    print('start import data from station %s date: %s' % (station, date_array[0]))
    station_data = get_superMAG_value_from_web(date_array[0], station)
    print('superMAG ground station %s %s len %s' % (station, date_array[0], len(station_data)))
    this_day_time_array = np.arange(station_data[0, 0], station_data[0, 0] + datetime.timedelta(days=1) + datetime.timedelta(seconds=1),
                                    datetime.timedelta(seconds=1)).astype(datetime.datetime)
    station_sec_time = []
    station_sec_value_N = []
    station_sec_value_E = []
    station_sec_value_Z = []
    for i, date in enumerate(date_array):
        if i == 0:
            station_times = station_data[:, 0]
            station_values = station_data[:, 1:]
            day_data_per_sec_N = fill_day_value(this_day_time_array, station_times, station_values[:, 0])
            station_sec_time.extend(this_day_time_array)
            station_sec_value_N.extend(day_data_per_sec_N)
            day_data_per_sec_E = fill_day_value(this_day_time_array, station_times, station_values[:, 1])
            station_sec_value_E.extend(day_data_per_sec_E)
            day_data_per_sec_Z = fill_day_value(this_day_time_array, station_times, station_values[:, 2])
            station_sec_value_Z.extend(day_data_per_sec_Z)

        if i != 0 and date_array[i-1] != date_array[i]:
            station_data = get_superMAG_value_from_web(date, station)
            this_day_time_array = np.arange(station_data[0, 0], station_data[0, 0] + datetime.timedelta(days=1),
                                            datetime.timedelta(seconds=1)).astype(datetime.datetime)
            station_times = station_data[:, 0]
            station_values = station_data[:, 1:]
            day_data_per_sec_N = fill_day_value(this_day_time_array, station_times, station_values[:, 0])
            station_sec_time.extend(this_day_time_array)
            station_sec_value_N.extend(day_data_per_sec_N)
            day_data_per_sec_E = fill_day_value(this_day_time_array, station_times, station_values[:, 1])
            station_sec_value_E.extend(day_data_per_sec_E)
            day_data_per_sec_Z = fill_day_value(this_day_time_array, station_times, station_values[:, 2])
            station_sec_value_Z.extend(day_data_per_sec_Z)



    #print(station_sec_time)
    #sw_dt = np.array(decode_str_dt_param(date_array[i] + 'T' + time_array[i]) for i in range(len(date_array)))

    station_delta_value_N = []
    station_delta_value_E = []
    station_delta_value_Z = []
    for i, cd in enumerate(date_array):
        sw_dt = decode_str_dt_param(cd + 'T' + time_array[i])
        for k, station_dt in enumerate(station_sec_time):
            if station_dt == sw_dt:
                station_delta_value_N.append(station_sec_value_N[k])
                station_delta_value_E.append(station_sec_value_E[k])
                station_delta_value_Z.append(station_sec_value_Z[k])
    station_delta_value_F = np.sqrt(np.power(station_delta_value_N, 2) + np.power(station_delta_value_E, 2)
                                    + np.power(station_delta_value_Z, 2))

    return np.array([station_delta_value_N, station_delta_value_E, station_delta_value_Z, station_delta_value_F]).T




def get_avroral_oval_file():
    """получение аврорального овала из файла"""
    # avroral_oval_data = pd.read_csv(BASE_DIR + '/static/data/20161030_0005_north_forecast_diff_aacgm.txt', skiprows=1, header=None, decimal='.', sep='\t')[0]
    # avroral_oval_data = avroral_oval_data.str.split('   ',expand=True).to_numpy(dtype=float)

    avroral_oval_data = pd.read_csv(BASE_DIR + '/static/data/20150101_0015_north_forecast_aacgm.txt', skiprows=1,
                                    header=None, decimal='.', sep='\t')
    row = avroral_oval_data[avroral_oval_data[0] == 'Totals:'].index.tolist()[0]
    avroral_oval_data = avroral_oval_data.iloc[:row - 1][0]
    avroral_oval_data = avroral_oval_data.str.split('   ', expand=True)

    avroral_oval_data.columns = ['LON', 'LAT', 'VALUE']
    avroral_oval_data = avroral_oval_data.to_numpy(dtype=float)
    avroral_oval_data[:, 0] = avroral_oval_data[:, 0] / 23.75 * 360  # x 24h setup to 360 degrees

    return avroral_oval_data  # lat, lon, value



def get_shapefile():
    """загрузка shapefile"""
    shpFilePath = BASE_DIR + '/static/data/perm_ice_shp/PERM_ICE_poly.shp'
    sf = shapefile.Reader(shpFilePath)
    return sf


"""def mera_mu(swarm_resp, led):
    # TODO  wrong import
    pd.DataFrame(swarm_resp).to_csv(
        STATIC_OS_PATH + '/swarm/' + '%s.csv' %
        led, index=False, header=False, )
    cmd = [
        MATLAB_PATH,
        'matlab -nodesktop -nosplash -r "meas_extr("%s"); exit"' %
        led]
    p = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        cwd=BASE_DIR + '/swarm')

    time.sleep(10)
    meas = np.array(
        pd.read_csv(
            BASE_DIR +
            '/swarm/swarm_mera/' +
            '%s.csv' %
            led,
            header=None)).astype(float)
    print('mera success saves and read ')
    return meas
"""

def mera_mu(swarm_resp):
    meas = get_measure_mu(swarm_resp)
    return meas


def cart_to_cartGeomag_coord(mag_xyz, phi, theta):
    PHI_P, THETA_P = math.radians(phi), math.radians(theta)

    R = np.array([
        [np.cos(PHI_P) * np.cos(THETA_P), -1 * np.sin(PHI_P), np.cos(PHI_P) * np.sin(THETA_P)],
        [np.sin(PHI_P) * np.cos(THETA_P), np.cos(PHI_P), np.sin(PHI_P) * np.sin(THETA_P)],
        [-1 * np.sin(THETA_P), 0, np.cos(THETA_P)]])

    geo_data = np.empty((0, 3))
    for xyz in mag_xyz:
        geo_xyz = np.sum(np.array([np.deg2rad(xyz)]).T * R.T, axis=0)
        geo_data = np.append(geo_data, [np.rad2deg(geo_xyz)], axis=0)
    return geo_data

def switch_coord_GEO_zeromodel_to_geomag(X, Y, magpol_theta, magpol_phi ):
    PHI_POL, THETA_POL = math.radians(magpol_phi), math.radians(magpol_theta)
    xy_mag = np.empty((0, 2))
    for i, x in enumerate(X):
        x, y = math.radians(x), math.radians(Y[i])
        x_ = math.cos(x-THETA_POL) + np.sin(y-PHI_POL)
        y_ = -math.sin(x-THETA_POL) + math.cos(y-PHI_POL)
        xy_mag = np.append(xy_mag, [[np.rad2deg(x_), np.rad2deg(y_)]], axis=0)
    return xy_mag

def rotate_midnight_coord(lat, lon, lt):
    def M(phi):
        return np.array([[np.cos(phi), -np.sin(phi)],
                         [np.sin(phi), np.cos(phi)]])
    phi_0 = 2*np.pi * (24/24)
    psi = np.pi/2 - np.deg2rad(lat)
    phi = np.deg2rad(lon)

    ans = np.matmul(M(phi_0), np.array([psi*np.cos(phi), psi*np.sin(phi)]).T)
    return np.rad2deg(ans.T)


def get_auroral_flux(dt, atype='diff', jtype='energy', hemishpere='N'):
    """
    Test automatic generation of omni_intervals in ovation_utilities
    also by not specifying a start and end time for the FluxEstimator
    """
    """
    fluxtupleN = estimatorN.get_gridded_flux(dF, combined_N_and_S=False)
    (mlatgridN, mltgridN, fluxgridN) = fluxtupleN[:3]

    fluxtupleS = estimatorS.get_gridded_flux(dF, combined_N_and_S=False)
    (mlatgridS, mltgridS, fluxgridS) = fluxtupleS[3:]
    """

    print('get auroral oval from date: %s' % dt)
    estimator = ovation_prime.FluxEstimator(atype, jtype)

    if hemishpere == 'N':
        THETA_P, PHI_P  = 213.00, 85.54
        mlatgrid, mltgrid, fluxgrid = estimator.get_flux_for_time(dt, hemi='N')

    else:
        THETA_P, PHI_P = 223., -64.24
        mlatgrid, mltgrid, fluxgrid = estimator.get_flux_for_time(dt, hemi='S')

    Y = mlatgrid.flatten()
    X = mltgrid.flatten() / 24 * 360
    geo_latlon = geomag2geo(latlon=np.array([Y, X]).T, THETA=THETA_P, PHI=PHI_P)
    Y, X = geo_latlon[:, 0], geo_latlon[:, 1]

    #midnight_lon = get_lon_from_LT(dt)
    midnight_lat, midnight_lon = sun_pos(dt)
    X = midnight_lon + X
    for i in range(len(X)):
        if X[i] < 0:
            X[i] = 360 + X[i]
        elif X[i] >= 360:
            X[i] = X[i] - 360

    Xmag, Ymag = X, Y

    #XYmag = np.array([rotate_midnight_coord(lat=Y[i], lon=X[i], lt=mltgrid.flatten()[i]) for i in range(len(X))])
    #Y, X = satplottools.latlt2cart(mlatgrid.flatten(), mltgrid.flatten(), hemishpere)
    #desc_mag = [gc_latlon_to_xyz(lat=lat, lon=lon, R=6371) for lat, lon in zip(Y, X)]
    #desc_mag = cart_to_cartGeomag_coord(desc_mag, THETA_P, PHI_P)
    #Y, X = satplottools.latlt2polar(mlatgrid.flatten(), mltgrid.flatten(), hemishpere)


    """X, Y = satplottools.latlt2cart(mlatgridS.flatten(), mltgridS.flatten(), 'S')
    desc_mag = xyz_to_cast(X, Y)
    desc_geo = cart_to_geomag_coord(desc_mag, PHI_P, THETA_P)
    Xmag, Ymag = xyz_to_spher(desc_geo)"""

    print(X.min(), X.max(), np.mean(X), 'X')
    print(Y.min(), Y.max(), np.mean(Y), 'Y')
    #xi = np.linspace(Xmag.min(), Xmag.max() + 1, 125)
    xi = np.linspace(0, 360., 125)
    yi = np.linspace(Ymag.min(), Ymag.max() + 1, 125)
    value = griddata((Xmag, Ymag), fluxgrid.flatten(), (xi[None, :], yi[:, None]), method='nearest')  # create a uniform spaced grid
    XX, YY = np.meshgrid(xi, yi)
    x_2d, y_2d, value_2d = X, Y, fluxgrid.flatten()

    #   lat, lon, fluxgridN  north
    #   lat, lon, fluxgridS  south

    return XX, YY, value, np.array([y_2d, x_2d, value_2d]).T

def eucl_distance(p1, p2_array):
    distY = (p1[0] - p2_array[:, 1]) ** 2  # check if y > deg
    distX = (p1[1] - p2_array[:, 0]) ** 2  # check if x > deg
    dist = np.sum(np.array([distX, distY]).T, axis=1)
    return np.sqrt(dist)

def get_nearest_auroral_point_to_swarm(swarm_set, atype):
    swarm_liter, swarm_pos, swarm_date, swarm_time, swarm_values = swarm_set
    print('nearest auroral point to swarm will get use data from swarm %s dt from %s to %s' % (
        swarm_liter, swarm_date[0] + ' ' + swarm_time[0],
        swarm_date[-1] + ' ' + swarm_time[-1]))

    #   auroral без reshape
    estimator = ovation_prime.FluxEstimator(atype=atype, energy_or_number='energy')
    def get_auraral_xyz(datetime):
        midnight_lat, midnight_lon = sun_pos(datetime)
        THETA_N, PHI_N = 213.00, 85.54
        THETA_S, PHI_S = 223., -64.24

        mlatgridN, mltgridN, fluxgridN = estimator.get_flux_for_time(datetime, hemi='N')
        mlatgridS, mltgridS, fluxgridS = estimator.get_flux_for_time(datetime, hemi='S')
        YN = mlatgridN.flatten()
        XN = mltgridN.flatten() / 24 * 360
        geo_latlon = geomag2geo(latlon=np.array([YN, XN]).T, THETA=THETA_N, PHI=PHI_N)
        YN, XN = geo_latlon[:, 0], geo_latlon[:, 1]
        XN = midnight_lon + XN
        for i in range(len(XN)):
            if XN[i] < 0:
                XN[i] = 360 + XN[i]
            elif XN[i] >= 360:
                XN[i] = XN[i] - 360

        YS = mlatgridS.flatten()
        XS = mltgridS.flatten() / 24 * 360
        geo_latlon = geomag2geo(latlon=np.array([YS, XS]).T, THETA=THETA_S, PHI=PHI_S)
        YS, XS = geo_latlon[:, 0], geo_latlon[:, 1]
        XS = midnight_lon + XS
        for i in range(len(XS)):
            if XS[i] < 0:
                XS[i] = 360 + XS[i]
            elif XS[i] >= 360:
                XS[i] = XS[i] - 360

        X = np.append(XN, XS)
        X = (X + 180) % 360 - 180  # 0 360 to -180 180
        Y = np.append(YN, YS)
        XY = np.vstack((X, Y)).T
        Z = np.append(fluxgridN, fluxgridS)
        return XY, Z

    near_auroral_array = []
    zero_time = decode_str_dt_param(swarm_date[0] + 'T' + swarm_time[0])
    XY, Z = get_auraral_xyz(decode_str_dt_param(swarm_date[0] + 'T' + swarm_time[0]))
    for j, pos in enumerate(swarm_pos):
        time_delta = decode_str_dt_param(swarm_date[j] + 'T' + swarm_time[j]) - pd.Timedelta('%s minutes' % 5)
        if time_delta > zero_time:
            #print('new auroral at time:', decode_str_dt_param(swarm_date[j] + 'T' + swarm_time[j]))
            XY, Z = get_auraral_xyz(decode_str_dt_param(swarm_date[j] + 'T' + swarm_time[j]))
            zero_time = decode_str_dt_param(swarm_date[j] + 'T' + swarm_time[j])
        eucl = eucl_distance(pos, XY)
        #print(eucl[np.argmin(eucl)], decode_str_dt_param(swarm_date[j] + 'T' + swarm_time[j]), pos[0], XY[j])
        if eucl[np.argmin(eucl)] < 5:
            near_auroral_array.append(Z[np.argmin(eucl)])
        else:
            near_auroral_array.append(np.nan)

    return np.array(near_auroral_array)


def open_superMAG_vcf(date):

    filename = date.split('-')
    filename = '%s%s%s' % (filename[0][2:], filename[1], filename[2])
    with open(DATA_PATH + '%s.T47' % filename, mode='r') as vcf:
        lines = vcf.readlines()

    station_data = []
    zero_time = decode_str_dt_param(date+'T'+'00:00:00')
    for i, line in enumerate(lines):
        if i == 0:
            continue
        line = np.array(line.rstrip('\n').split(',')[0].split('    ')[1:]).astype(float)    # ['-111.600', ' -13.900', '-158.200']
        zero_time = zero_time + pd.Timedelta('%s minutes' % 1)
        station_data.append([zero_time, line[0], line[1], line[2]])
    return np.array(station_data)



"""
########################################################################################

изменения системы координат
"""
# lat lon r geocentric to lat lon geodetic

def geocentric_to_geodetic(sat_data, geoid_alt=False):
    """
        # geocentric_to_geodetic
        :param sat_data[:, 0]: geocentric lat -90... 90
        :param sat_data[:, 1]: geocentric lon -180... 180
        :param sat_data[:, 2]: alt above center of Earth, km
        :geoid_alt==True -> alt above Earth sea lvl (geoid), km
        :return: [geodetic lat, geodetic lon, alt]
        """
    geocentric_xyz = np.empty((0, 3))
    for i in range(len(sat_data)):
        # lat, lon, r
        geocentric_xyz = np.append(geocentric_xyz, [gc_latlon_to_xyz(sat_data[i, 0], sat_data[i, 1], sat_data[i, 2], )],
                                   axis=0)
    geodetic_latlonR = gc2gd(geocentric_xyz)
    if geoid_alt:
        geocentric_latlonAlt = np.append(sat_data[:, :2], np.array([geodetic_latlonR[:, 2]]).T, axis=1).astype(float)
        return geocentric_latlonAlt
    else:
        return geodetic_latlonR

def gc2gd(gc_xyz):

    wgs84 = Proj('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')
    geocentric = Proj('+proj=geocent +datum=WGS84 +units=m +no_defs')

    # same
    #wgs84 = CRS("EPSG:4326")   # (deg)  # Horizontal component of 3D system. Used by the GPS satellite navigation system and for NATO military geodetic surveying.
    #geocentric = CRS("EPSG:4978")    # WGS 84 (geocentric) # X OTHER, Y EAST, Z NORTH,
    x, y, z = gc_xyz[:, 0], gc_xyz[:, 1], gc_xyz[:, 2]
    lon, lat, alt = transform(geocentric, wgs84, x, y, z)
    alt = alt/1000  # m to km
    return np.array([lat, lon, alt]).T

def gc_latlon_to_xyz(lat, lon, R):

    # Convert (lat, lon, elv) to (x, y, z).
    lat = lat * math.pi / 180.0
    lon = lon * math.pi / 180.0
    radius = R*1000     # km to m
    geo_centric_lat = geocentric_latitude(lat)

    cos_lon = math.cos(lon)
    sin_lon = math.sin(lon)
    cos_lat = math.cos(geo_centric_lat)
    sin_lat = math.sin(geo_centric_lat)
    x = radius * cos_lon * cos_lat
    y = radius * sin_lon * cos_lat
    z = radius * sin_lat

    return x, y, z

def geocentric_latitude(lat):
    # Convert geodetic latitude 'lat' to a geocentric latitude 'clat'.
    # Geodetic latitude is the latitude as given by GPS.
    # Geocentric latitude is the angle measured from center of Earth between a point and the equator.
    # https:#en.wikipedia.org/wiki/Latitude#Geocentric_latitude
    e2 = 0.00669437999014
    return math.atan((1.0 - e2) * math.tan(lat))


############


def get_local_time(dt, lon):
    UTtime_h, UTtime_m, UTtime_s = dt.hour, dt.minute, dt.second
    #loct = np.rint(lon / 15) + UTtime_h + (UTtime_m / 60) + (UTtime_s / 3600)
    loct = int(lon / 15) + UTtime_h + (UTtime_m / 60) + (UTtime_s / 3600)
    if loct < 0:
        loct = 24 + loct
    if loct >= 24:
        loct = loct - 24
    local_time = datetime.datetime.strptime('%s %s' %
                                                ('%i-%i-%i' % (dt.year, dt.month, dt.day),
                                                 '%i:%i:%i' % (loct, UTtime_m, UTtime_s)),
                                            '%Y-%m-%d %H:%M:%S')
    local_time = local_time.time()
    return local_time

def get_lon_from_LT(dt):
    """print('time:', dt.time())
    lon_lt = (dt.hour*15)
    try:
        lon_lt += (dt.minute/60*15)
    except:
        print('minute is 0')
    try:
        lon_lt += (dt.second/60*15)
    except:
        print('second is 0')
    #return 360-lon_lt
    return -lon_lt
    #return lon_lt
    #return 2*np.pi *(dt.hour/1440)*360"""

    seconds_per_day = 24 * 60 * 60.0
    sec_since_midnight = (dt - datetime.datetime(dt.year, dt.month, dt.day)).seconds
    lng = -(sec_since_midnight / seconds_per_day) * 360
    #lng = (lng + 180) % 360 - 180  # 0 360 to -180 180
    lng = (lng + 180) % 360 - 180  # 0 360 to -180 180

    return lng


def get_solar_coord(date, lon):
    lon_to360 = lambda x: (x - 180) % 360 - 180  # -180 180 to 0 360
    lon_to180 = lambda x: (x + 180) % 360 - 180  # 0 360 to -180 180
    midnight_lat, midnight_lon = sun_pos(date)

    #midnight_lon = get_lon_from_LT(date)
    solar_coord_lon = midnight_lon + lon_to360(lon)
    if solar_coord_lon < 0:
        solar_coord_lon = 360 + solar_coord_lon
    elif solar_coord_lon >= 360:
        solar_coord_lon = solar_coord_lon - 360
    return solar_coord_lon


def earth_radius_in_meters(latitude_radians):
    # latitudeRadians is geodetic, i.e. that reported by GPS.
    # http:#en.wikipedia.org/wiki/Earth_radius
    a = 6378137.0 # equatorial radius in meters
    b = 6356752.3 # polar radius in meters
    cos = math.cos(latitude_radians)
    sin = math.sin(latitude_radians)
    t1 = a * a * cos
    t2 = b * b * sin
    t3 = a * cos
    t4 = b * sin
    return math.sqrt((t1*t1 + t2*t2) / (t3*t3 + t4*t4))








