import numpy as np
from teslaswarm.settings import BASE_DIR, STATIC_OS_PATH, OBS_DATA_PATH
from program_meas_py.meas_extr import get_measure_mu
import math
from numpy import degrees, radians
import pandas as pd
import shapefile
import subprocess
import time
from ovationpyme import ovation_prime, ovation_utilities
from geospacepy import satplottools, special_datetime
import pyIGRF
from pyIGRF import calculate
from tools.dt_foo import *
from tools.supermag_api import SuperMAGGetData, sm_grabme
from scipy.interpolate import griddata
import shapely
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
import json
import geog
from spacepy import coordinates as coord
from spacepy.time import Ticktock
import aacgmv2
import warnings
from pyproj import Proj, transform, CRS




def data_reduction(respond, delta, fac2_mod=False):
    warnings.simplefilter("ignore", category=RuntimeWarning)
    """сжимание секундных данных до delta шага"""
    # fac2 = (Y, X, R), dt, (fac2)
    # vector, measure mu, fac, chaos = (Y, X, R), dt, (N, E, C)
    N, M = respond.shape
    if fac2_mod:
        idx999 = np.where(respond[:, 4] == 999)[0]
        # respond = respond[idx999]
        respond[idx999, 4] = np.nan
        redu_resp = np.empty((0, 5))
    else:
        redu_resp = np.empty((0, M))

    window = int(delta / 2)
    if window == 0:
        window = 1
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

    if fac2_mod:
        #https://pandas.pydata.org/docs/user_guide/missing_data.html
        redu_resp[:, 4] = pd.DataFrame(np.array(redu_resp[:, 4]), dtype='float32').interpolate().to_numpy().T[0]


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
        x = np.mean([cd_a[1], cd_c[1]])
        y = np.mean([cd_a[0], cd_c[0]])
        r = np.mean([cd_a[2], cd_c[2]])
        coord = np.append(coord, [[y, x, r]], axis=0)

        if channel == 3:
            v = np.sqrt(np.power(av[0] - cv[0], 2) +
                        np.power(av[1] - cv[1], 2) +
                        np.power(av[2] - cv[2], 2))
        else:
            v = np.sqrt(np.power(av[channel] - cv[channel], 2))
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
            mag_lat_lines.append(geo2mag(np.array(lat_lines), np.full((1, len(lat_lines)), date)[0], to_coord='MAG'))

        else:
            mag_lat_lines.append(geo2mag(np.array(lat_lines), np.full((1, len(lat_lines)), date)[0], to_coord='GSM'))


    for lon_lines_start in np.arange(lon_coord[0], lon_coord[1], lon_step):
        lon_lines = []
        for lat in np.arange(-90, 91, 1):
            lon_lines.append([lat, lon_lines_start, r])

        #   convert to mag\geomag coords
        if geomag_pole:
            #mag_lon_lines.append(apex_coord_system_convert(np.array(lon_lines)[:, :2], source_sys='geo', dest_sys='apex'))
            #mag_lat_lines.append(geo2geomag(np.array(lon_lines)[:, :2]))
            mag_lon_lines.append(geo2mag(np.array(lon_lines), np.full((1, len(lon_lines)), date)[0], to_coord='MAG'))

        else:
            mag_lon_lines.append(geo2mag(np.array(lon_lines), np.full((1, len(lon_lines)), date)[0], to_coord='GSM'))

    for lat in np.arange(-90+lat_step, 91-lat_step, lat_step):
        for lon in np.arange(-180, 181., lon_step):
            if geomag_pole:
                ap = geo2mag(np.array([[lat, lon, r]]), np.full((1, 1), date)[0], to_coord='MAG')[0]
            else:
                ap = geo2mag(np.array([[lat, lon, r]]), np.full((1, 1), date)[0], to_coord='GSM')[0]

            annotate_points.append([ap[0], ap[1], [lat, lon]])
    for lat in [-90, 90]:
        if geomag_pole:
            ap = geo2mag(np.array([[lat, 0, r]]), np.full((1, 1), date)[0], to_coord='MAG')[0]
        else:
            ap = geo2mag(np.array([[lat, 0, r]]), np.full((1, 1), date)[0], to_coord='GSM')[0]
        annotate_points.append([ap[0], ap[1], [lat, 0]])



    return mag_lat_lines, mag_lon_lines, annotate_points

def swarm_egrf_vector_subtraction(swarm_pos, swarm_values_full, swarm_date):
    #TODO igrf12 to 13
    year = int(swarm_date[0].split("-")[0])
    #switched_swarm_pos = convert_coord_system(swarm_pos[:, :2], dest_sys='apex')  # convert geo to mag coords
    #switched_swarm_pos = geo2mag(swarm_pos, swarm_date, to_coord='GSM')
    B = []
    idx = 0
    for n, e in swarm_values_full:
        """  
                d, i, h, x, y, z, f = calculate_magnetic_field_intensity(lon=switched_swarm_pos[idx, 1], lat=switched_swarm_pos[idx, 0], alt=switched_swarm_pos[idx, 2], date=date)
                d, i, h, x, y, z, f = variation_of_magnetic_filed_intensity(lon=swarm_pos[idx, 1], lat=swarm_pos[idx, 0], alt=swarm_pos[idx, 2], date=date)
                print('n:%s e:%s' % (n, e))
                print('x:%s y:%s' % (x, y))
                Bb = np.sqrt((n-x)**2+(e-y)**2)
                Bb = np.sqrt((n-x)**2+(e-y)**2)
                print('d:%s dx:%s dy%s' % (dd, dx, dy))
                Bd = dd - d
                print(dd, end='======\n')
                b1, b2 = np.sqrt(n-x), np.sqrt(e-y)

                """
        d, i, h, x, y, z, f = igrf_value(lat=swarm_pos[idx, 0], lon=180.-swarm_pos[idx, 1],  alt=swarm_pos[idx, 2], year=year)
        #d, i, h, x, y, z, f = igrf_value(lat=switched_swarm_pos[idx, 0], lon=switched_swarm_pos[idx, 1], alt=switched_swarm_pos[idx, 2], year=year)
        #print('sw_n:%.2f sw_e:%.2f x:%.2f y:%.2f' % (n, e, x, y))
        dd, dx, dy = magfield_variation(n, e, x, y)
        B.append([dd, dx, dy])
        idx += 1
    return np.array(B)


def igrf_value(lat, lon, alt=0., year=2005.):

    FACT = 180. / np.pi
    """
    calculate.igrf12syn
         This is a synthesis routine for the 12th generation IGRF as agreed
         in December 2014 by IAGA Working Group V-MOD. It is valid 1900.0 to
         2020.0 inclusive. Values for dates from 1945.0 to 2010.0 inclusive are
         definitive, otherwise they are non-definitive.
       INPUT
         date  = year A.D. Must be greater than or equal to 1900.0 and
                 less than or equal to 2025.0. Warning message is given
                 for dates greater than 2020.0. Must be double precision.
         itype = 1 if geodetic (spheroid)
         itype = 2 if geocentric (sphere)
         alt   = height in km above sea level if itype = 1
               = distance from centre of Earth in km if itype = 2 (>3485 km)
         lat = latitude (-90~90)
         elong = east-longitude (0-360)
         alt, colat and elong must be double precision.
       OUTPUT
         x     = north component (nT) if isv = 0, nT/year if isv = 1
         y     = east component (nT) if isv = 0, nT/year if isv = 1
         z     = vertical component (nT) if isv = 0, nT/year if isv = 1
         f     = total intensity (nT) if isv = 0, rubbish if isv = 1
         
        """
    x, y, z, f = calculate.igrf12syn(year, 2, alt, lat, lon)
    d = FACT * np.arctan2(y, x)
    h = np.sqrt(x * x + y * y)
    i = FACT * np.arctan2(z, h)
    """
        :return
             D is declination (+ve east)
             I is inclination (+ve down)
             H is horizontal intensity
             X is north component
             Y is east component
             Z is vertical component (+ve down)
             F is total intensity
    """
    return d, i, h, x, y, z, f

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


def get_superMAG_value_from_web(date, station):
    #print(date, station)
    answer = SuperMAGGetData(logon='pilipenko',start='%sT00:00:00'%date, extent=86400,
                             station=station, flagstring='',)
    #print(answer)
    station_data = []
    i, zero_time = 0, decode_str_dt_param(date + 'T' + '00:00:00')
    for n, e, c in zip(sm_grabme(answer[1], 'N', 'geo'), sm_grabme(answer[1], 'E', 'geo'), sm_grabme(answer[1], 'Z', 'geo')):
        zero_time = zero_time + pd.Timedelta('%s minutes' % 1)
        station_data.append([zero_time, n, e, c])
        i += 1
    return np.array(station_data)

def get_sMAGstation_value_by_time(date_array, time_array, channel, delta, station):
    #print(channel, station)

    """
        def nearest(items, pivot):
        min_idx = 0
        min_value = None
        for k, x in enumerate(items):
            value = abs(x - pivot)
            if min_value is not None:
                if value <= min_value:
                    min_value = value
                    min_idx = k
            else:
                min_value = value
                min_idx = k
        return min_idx
    """

    def fill_day_value(this_day_time_array, station_times, station_values):
        def get_future_value(v1, v2, step):
            if v1 > v2:
                return np.flip(np.arange(v2, v1, step))
            else:
                return np.arange(v1, v2, step)

        extend_station_data = []
        station_last_value_idx = 0
        future_value_idx = 0
        station_next_dt = station_times[station_last_value_idx + 1]

        for i, ct in enumerate(this_day_time_array):
            if i == 0:
                step = np.abs(station_values[station_last_value_idx]-station_values[station_last_value_idx + 1])/60
                future_values = get_future_value(station_values[station_last_value_idx],
                                                 station_values[station_last_value_idx+1], step)
            if ct != station_next_dt:

                extend_station_data.append([ct, future_values[future_value_idx]])
                future_value_idx += 1
            else:
                future_value_idx = 0
                station_last_value_idx += 1
                if station_last_value_idx+1 < len(station_times):
                    station_next_dt = station_times[station_last_value_idx+1]
                else:
                    continue

                extend_station_data.append([ct, station_values[station_last_value_idx]])
                step = np.abs(station_values[station_last_value_idx]-station_values[station_last_value_idx + 1])/60
                future_values = get_future_value(station_values[station_last_value_idx],
                                                 station_values[station_last_value_idx + 1], step)
        return np.array(extend_station_data)

    station_data = get_superMAG_value_from_web(date_array[0], station)
    print('superMAG ground station %s %s len %s' % (station, date_array[0], len(station_data)))
    this_day_time_array = np.arange(station_data[0, 0], station_data[0, 0] + datetime.timedelta(days=1),
                                    datetime.timedelta(seconds=1)).astype(datetime.datetime)
    station_sec_time = []
    station_sec_value = []
    for i, date in enumerate(date_array):
        if i == 0:
            station_times = station_data[:, 0]
            station_values = station_data[:, 1:]
            day_data_per_sec = fill_day_value(this_day_time_array, station_times, station_values[:, channel])
            station_sec_time.extend(day_data_per_sec[:, 0])
            station_sec_value.extend(day_data_per_sec[:, 1])
        if i != 0 and date_array[i-1] != date_array[i]:
            station_data = get_superMAG_value_from_web(date, station)
            this_day_time_array = np.arange(station_data[0, 0], station_data[0, 0] + datetime.timedelta(days=1),
                                            datetime.timedelta(seconds=1)).astype(datetime.datetime)
            station_times = station_data[:, 0]
            station_values = station_data[:, 1:]
            day_data_per_sec = fill_day_value(this_day_time_array, station_times, station_values[:, channel])
            station_sec_time.extend(day_data_per_sec[:, 0])
            station_sec_value.extend(day_data_per_sec[:, 1])
    #print(station_sec_value)
    #print(station_sec_time)


    station_delta_value = []
    for i, cd in enumerate(date_array):
        sw_dt = decode_str_dt_param(cd+'T'+time_array[i])
        for k, station_dt in enumerate(station_sec_time):
            if station_dt == sw_dt:
                station_delta_value.append(station_sec_value[k])

    #print(np.array(station_delta_value))
    return np.array(station_delta_value)




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
        mlatgridN, mltgridN, fluxgridN = estimator.get_flux_for_time(dt, hemi='N')
    else:
        mlatgridS, mltgridS, fluxgridS = estimator.get_flux_for_time(dt, hemi='S')
    if hemishpere == 'N':
        PHI_P, THETA_P = 252.68, 80.65

        #XN, YN = satplottools.latlt2cart(mlatgridN.flatten(), mltgridN.flatten(), 'N')
        XN, YN = satplottools.latlt2cart(mlatgridN.flatten(), mltgridN.flatten(), 'N')
        desc_mag = xyz_to_cast(XN, YN)
        desc_geo = geo_to_mag_coord(desc_mag, PHI_P, THETA_P)
        XN, YN = xyz_to_spher(desc_geo)
    else:
        PHI_P, THETA_P = 107.62, -80.37
        #XS, YS = satplottools.latlt2cart(mlatgridS.flatten(), mltgridS.flatten(), 'S')
        XS, YS = satplottools.latlt2cart(mlatgridS.flatten(), mltgridS.flatten(), 'S')
        desc_mag = xyz_to_cast(XS, YS)
        desc_geo = geo_to_mag_coord(desc_mag, PHI_P, THETA_P)
        XS, YS = xyz_to_spher(desc_geo)

    if hemishpere == 'N':
        xi = np.linspace(XN.min(), XN.max() + 1, 125)
        yi = np.linspace(YN.min(), YN.max() + 1, 125)
        value = griddata((XN, YN), fluxgridN.flatten(), (xi[None, :], yi[:, None]), method='linear')  # create a uniform spaced grid
        XX, YY = np.meshgrid(xi, yi)
        x_2d, y_2d, value_2d = XN, YN, fluxgridN.flatten()

    else:
        xi = np.linspace(XS.min(), XS.max() + 1, 125)
        yi = np.linspace(YS.min(), YS.max() + 1, 125)
        value = griddata((XS, YS), fluxgridS.flatten(), (xi[None, :], yi[:, None]), method='linear')  # create a uniform spaced grid
        XX, YY = np.meshgrid(xi, yi)
        x_2d, y_2d, value_2d = XS, YS, fluxgridS.flatten()

    #   lat, lon, fluxgridN  north
    #   lat, lon, fluxgridS  south

    return XX, YY, value, np.array([y_2d, x_2d, value_2d]).T

def eucl_distance(p1, p2_array):
    distY = (p1[0] - p2_array[:, 1]) ** 2  # check if y > deg
    distX = (p1[1] - p2_array[:, 0]) ** 2  # check if x > deg
    dist = np.sum(np.array([distX, distY]).T, axis=1)
    return np.sqrt(dist)

def get_nearest_auroral_point_to_swarm(swarm_set):
    swarm_liter, swarm_pos, swarm_date, swarm_time, swarm_values = swarm_set
    print('nearest auroral point to swarm will get use data from swarm %s dt from %s to %s' % (
        swarm_liter, swarm_date[0] + ' ' + swarm_time[0],
        swarm_date[-1] + ' ' + swarm_time[-1]))

    #   auroral без reshape
    estimator = ovation_prime.FluxEstimator('diff', 'energy')
    def get_auraral_xyz(datetime):
        #middle_datetime = decode_dt_param( swarm_date[int((len(swarm_time) / 2))] + 'T' + swarm_time[int((len(swarm_time) / 2))])
        mlatgridN, mltgridN, fluxgridN = estimator.get_flux_for_time(datetime, hemi='N')
        mlatgridS, mltgridS, fluxgridS = estimator.get_flux_for_time(datetime, hemi='S')
        XN, YN = satplottools.latlt2cart(mlatgridN.flatten(), mltgridN.flatten(), 'N')
        desc_mag = xyz_to_cast(XN, YN)
        desc_geo = geo_to_mag_coord(desc_mag, 252.68, 80.65)
        XN, YN = xyz_to_spher(desc_geo)

        XS, YS = satplottools.latlt2cart(mlatgridS.flatten(), mltgridS.flatten(), 'S')
        desc_mag = xyz_to_cast(XS, YS)
        desc_geo = geo_to_mag_coord(desc_mag, 107.62, -80.37)
        XS, YS = xyz_to_spher(desc_geo)

        X = np.append(XN, XS)
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
        #print(eucl[np.argmin(eucl)])
        if eucl[np.argmin(eucl)] < 5:
            near_auroral_array.append(Z[np.argmin(eucl)])
        else:
            near_auroral_array.append(np.nan)

    return np.array(near_auroral_array)



def open_superMAG_vcf(date):

    filename = date.split('-')
    filename = '%s%s%s' % (filename[0][2:], filename[1], filename[2])
    with open(OBS_DATA_PATH+'%s.T47' % filename, mode='r') as vcf:
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


def geo2mag(sw_pos, swarm_date, to_coord='MAG'):
    #input sw_pos format: lat, lon, rad

    # coordinate init_coords = [Z, Y, X] or [rad, lat, lon]
    init_coords = np.array([sw_pos[:, 2] * 1000, sw_pos[:, 0], sw_pos[:, 1]]).T
    cvals = coord.Coords(init_coords, 'GEO', 'sph')


    dt_array = []
    for k, sw_date in enumerate(swarm_date):
        dt_array.append(str(sw_date)+'T'+'00:00:00')
    #print(dt_array)
    cvals.ticks = Ticktock(dt_array, 'ISO')

    if to_coord == 'MAG':
        # геомагнитные координаты
        newcoord = np.array(cvals.convert('MAG', 'sph').data)  # return radius, latitude, longitude

    elif to_coord == 'GSM':
        # Geocentric Solar Magnetospheric
        newcoord = np.array(cvals.convert('GSM', 'sph').data)  # return radius, latitude, longitude

    #lan, lon, r
    return np.array([newcoord[:, 1], newcoord[:, 2], newcoord[:, 0] / 1000]).T

def mag2mlt(sw_pos, swarm_date, swarm_time):
    #   magnetic longitude to MLT
    """
             When referencing this package, please cite both the package DOI and the AACGM-v2 journal article:

            Shepherd, S. G. (2014), Altitude‐adjusted corrected geomagnetic coordinates: Definition and functional approximations, Journal of Geophysical Research: Space Physics, 119, 7501–7521, doi:10.1002/2014JA020264.
            """
    convert_mlt = []
    for pos, date, time in zip(sw_pos, swarm_date, swarm_time):
        mlt_coords = np.array(aacgmv2.convert_mlt(pos, datetime.datetime.strptime(date+' '+time, '%Y-%m-%d %H:%M:%S'), m2a=False))
        convert_mlt.append(mlt_coords)

    return np.array(convert_mlt)


def geo2geomag(swarm_pos):
    switched_swarm_pos = []
    for lat, lon in swarm_pos:
        lat_a, lon_a = geo2geomag_foo(lat, lon)
        switched_swarm_pos.append([lat_a, lon_a])
    return np.array(switched_swarm_pos)

def geo2geomag_foo(glat: float, glon: float) :
    """
    Converts GEOGRAPHIC (latitude,longitude) to GEOMAGNETIC (latitude, longitude).
    Ground-level altitude
    This is just a rotation, so the accuracy is low compared to actual
    geomagnetic models like WMM, IGRF, etc.
    Latitudes and longitudes (east, 0..360) are expressed in degrees.
    They may be SCALAR, LIST or Numpy.ndarray (any shape and rank)
    > geo2mag(79.3,288.59) == pytest.approx([89.999992, -173.02325])
    Written by Pascal Saint-Hilaire (Saint-Hilaire@astro.phys.ethz.ch), May 2002
    https://github.com/wlandsman/IDLAstro
    converted to Python by Michael Hirsch June 2020
    http://wdc.kugi.kyoto-u.ac.jp/poles/polesexp.html "geomagnetic south pole"
    """

    # longitude (in degrees east) of Earth's magnetic south pole
    # 2015
    Dlong = 360 - (180 - 107.4)
    Dlat = 80.4

    # 1995
    # Dlong=288.59
    # Dlat=79.30

    Dlong = np.radians(Dlong)
    Dlat = np.radians(Dlat)

    R = 1

    glat = np.radians(glat)
    glon = np.radians(glon)
    galt = R

    # %% handle array shape
    glat_shape = glat.shape
    glon_shape = glon.shape

    glat = glat.ravel()
    glon = glon.ravel()

    # %% convert to rectangular coordinates
    #       X-axis: defined by the vector going from Earth's center towards
    #            the intersection of the equator and Greenwitch's meridian.
    #       Z-axis: axis of the geographic poles
    #       Y-axis: defined by Y=Z^X
    x = galt * np.cos(glat) * np.cos(glon)
    y = galt * np.cos(glat) * np.sin(glon)
    z = galt * np.sin(glat)

    # %% Compute 1st rotation matrix
    # rotation around plane of the equator,
    # from the Greenwich meridian to the meridian containing the magnetic
    # dipole pole.
    geolong2maglong = np.zeros((3, 3))
    geolong2maglong[0, 0] = np.cos(Dlong)
    geolong2maglong[1, 0] = np.sin(Dlong)
    geolong2maglong[0, 1] = -np.sin(Dlong)
    geolong2maglong[1, 1] = np.cos(Dlong)
    geolong2maglong[2, 2] = 1.0
    out = geolong2maglong.T @ np.array([x, y, z])

    # %% Second rotation
    # in the plane of the current meridian from geographic
    #                  pole to magnetic dipole pole.
    tomaglat = np.zeros((3, 3))
    tomaglat[0, 0] = np.cos(np.pi / 2 - Dlat)
    tomaglat[2, 0] = -np.sin(np.pi / 2 - Dlat)
    tomaglat[0, 2] = np.sin(np.pi / 2 - Dlat)
    tomaglat[2, 2] = np.cos(np.pi / 2 - Dlat)
    tomaglat[1, 1] = 1.0
    out = tomaglat.T @ out

    # %% convert back to latitude, longitude and altitude
    mlat = np.arctan2(out[2], np.sqrt(out[0] ** 2 + out[1] ** 2))
    mlat = np.degrees(mlat)
    mlon = np.arctan2(out[1], out[0])
    mlon = np.degrees(mlon)
    # malt=sqrt(out[0,*]^2+out[1,*]^2+out[2,*]^2)-R

    return mlat.reshape(glat_shape), mlon.reshape(glon_shape)


def latlt2cart(lat, lt, hemisphere):
    """
    Latitude and local time to cartesian for a top-down dialplot
    """
    r, theta = latlt2polar(lat, lt, hemisphere)

    return r * np.cos(theta), r * np.sin(theta)
    #return r, theta


def latlt2polar(lat, lt, hemisphere):
    """
    Converts an array of latitude and lt points to polar for a top-down dialplot (latitude in degrees, LT in hours)
    i.e. makes latitude the radial quantity and MLT the azimuthal

    get the radial displacement (referenced to down from northern pole if we want to do a top down on the north,
        or up from south pole if visa-versa)
    """
    if hemisphere == 'N':
        r = 90. - lat
    elif hemisphere == 'S':
        r = 90. - (-1 * lat)
    else:
        raise ValueError('%s is not a valid hemisphere, N or S, please!' % (hemisphere))
    # convert lt to theta (azimuthal angle) in radians
    theta = lt / 24. * 2 * np.pi - np.pi / 2

    # the pi/2 rotates the coordinate system from
    # theta=0 at negative y-axis (local time) to
    # theta=0 at positive x axis (traditional polar coordinates)
    return r, theta


def cart2polar(x, y):
    rho = np.sqrt(x ** 2 + y ** 2)
    phi = np.arctan2(y, x)
    return (rho, phi)


def pol2cart(rho, phi):
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return (x, y)

def to_cast(sph_x, sph_y):
    # geodetic_to_geocentric
    r = 6370
    desc_data = np.empty((0, 3))
    for i in range(len(sph_x)):
        # az, el = sph_x[i] * math.pi / 180, sph_y[i] * math.pi / 180
        az, el = math.radians(sph_x[i]), math.radians(sph_y[i])
        rcos_theta = r * np.cos(el)
        x = rcos_theta * np.cos(az)
        y = rcos_theta * np.sin(az)
        z = r * np.sin(el)
        desc_data = np.append(desc_data, [[x, y, z]], axis=0)
    return desc_data


def to_spher(data):
    spher_data = np.empty((0, 2))
    for x, y, z in data:
        hxy = np.hypot(x, y)
        # r = np.hypot(hxy, z)
        el = np.arctan2(z, hxy)
        az = np.arctan2(y, x)
        spher_data = np.append(spher_data, [[math.degrees(az), math.degrees(el)]], axis=0)
    return spher_data[:, 0], spher_data[:, 1]


def geo_to_mag_coord(mag_xyz, phi, theta):
    PHI_P, THETA_P = math.radians(phi), math.radians(-1 * theta)

    R = np.array([
        [np.cos(PHI_P) * np.cos(THETA_P), -1 * np.sin(PHI_P), np.cos(PHI_P) * np.sin(THETA_P)],
        [np.sin(PHI_P) * np.cos(THETA_P), np.cos(PHI_P), np.sin(PHI_P) * np.sin(THETA_P)],
        [-1 * np.sin(THETA_P), 0, np.cos(THETA_P)]])

    geo_data = np.empty((0, 3))
    for xyz in mag_xyz:
        geo_xyz = np.sum(np.array([xyz]).T * R.T, axis=0)
        geo_data = np.append(geo_data, [geo_xyz], axis=0)
    return geo_data

def xyz_to_cast(sph_x, sph_y):
    r = 6370
    desc_data = np.empty((0, 3))
    for i in range(len(sph_x)):
        # az, el = sph_x[i] * math.pi / 180, sph_y[i] * math.pi / 180
        az, el = math.radians(sph_x[i]), math.radians(sph_y[i])
        rcos_theta = r * np.cos(el)
        x = rcos_theta * np.cos(az)
        y = rcos_theta * np.sin(az)
        z = r * np.sin(el)
        desc_data = np.append(desc_data, [[x, y, z]], axis=0)
    return desc_data


def xyz_to_spher(data):
    spher_data = np.empty((0, 2))
    for x, y, z in data:
        hxy = np.hypot(x, y)
        # r = np.hypot(hxy, z)
        el = np.arctan2(z, hxy)
        az = np.arctan2(y, x)
        spher_data = np.append(spher_data, [[math.degrees(az), math.degrees(el)]], axis=0)
    return spher_data[:, 0], spher_data[:, 1]


def shiftlon(lon, lat):
    xy = np.empty((0, 2))
    for i, x, y in zip(range(len(lon)), lon, lat):
        try:
            if abs(x - lon[i + 1]) > 180:
                if np.sign(y) == -1:
                    y = -90
                    lat[i + 1] = -89
                else:
                    y = 90
                    lat[i + 1] = 89
        except BaseException:
            pass
        xy = np.append(xy, [[x, y]], axis=0)
    return xy


def geo_to_gg(radius, theta):
    """
    Compute geodetic colatitude and vertical height above the ellipsoid from
    geocentric radius and colatitude.

    Parameters
    ----------
    radius : ndarray, shape (...)
        Geocentric radius in kilometers.
    theta : ndarray, shape (...)
        Geocentric colatitude in degrees.

    Returns
    -------
    height : ndarray, shape (...)
        Altitude in kilometers.
    beta : ndarray, shape (...)
        Geodetic colatitude

    Notes
    -----
    Round-off errors might lead to a failure of the algorithm especially but
    not exclusively for points close to the geographic poles. Corresponding
    geodetic coordinates are returned as NaN.

    References
    ----------
    Function uses Heikkinen's algorithm taken from:

    Zhu, J., "Conversion of Earth-centered Earth-fixed coordinates to geodetic
    coordinates", IEEE Transactions on Aerospace and Electronic Systems}, 1994,
    vol. 30, num. 3, pp. 957-961

    """

    # Use WGS-84 ellipsoid parameters
    a = 6378.137  # equatorial radius
    b = 6356.752  # polar radius

    a2 = a ** 2
    b2 = b ** 2

    e2 = (a2 - b2) / a2  # squared eccentricity
    e4 = e2 * e2
    ep2 = (a2 - b2) / b2  # squared primed eccentricity

    r = radius * np.sin(radians(theta))
    z = radius * np.cos(radians(theta))

    r2 = r ** 2
    z2 = z ** 2

    F = 54 * b2 * z2

    G = r2 + (1. - e2) * z2 - e2 * (a2 - b2)

    c = e4 * F * r2 / G ** 3

    s = (1. + c + np.sqrt(c ** 2 + 2 * c)) ** (1. / 3)

    P = F / (3 * (s + 1. / s + 1.) ** 2 * G ** 2)

    Q = np.sqrt(1. + 2 * e4 * P)

    r0 = -P * e2 * r / (1. + Q) + np.sqrt(
        0.5 * a2 * (1. + 1. / Q) - P * (1. - e2) * z2 / (Q * (1. + Q)) - 0.5 * P * r2)

    U = np.sqrt((r - e2 * r0) ** 2 + z2)

    V = np.sqrt((r - e2 * r0) ** 2 + (1. - e2) * z2)

    z0 = b2 * z / (a * V)

    height = U * (1. - b2 / (a * V))

    beta = 90. - degrees(np.arctan2(z + ep2 * z0, r))

    return height, beta


"""
def sphericalToVector(inp_ar):
    ar = np.array([0.0, 0.0, 0.0])
    ar[0] = -sin(inp_ar[1]) 
    ar[1] = sin(inp_ar[0]) * cos(inp_ar[1]) 
    ar[2] = cos(inp_ar[0]) * cos(inp_ar[1])
    return ar
    
def vectorToGeogr(vect):
    ar = np.array([0.0, 0.0])
    ar[0] = np.arctan2(vect[1], vect[2])
    ar[1] = np.arcsin(-vect[0] / np.linalg.norm(vect))
    return ar

"""

"""
def apex_coord_system_convert(swarm_pos, source_sys, dest_sys):
    switched_swarm_pos = []
    for lat, lon in swarm_pos:
        lat_a, lon_a = apex_convert(lat, lon, source=source_sys, dest=dest_sys)
        switched_swarm_pos.append([lat_a, lon_a])
    return np.array(switched_swarm_pos)

def apex_convert(lat, lon, source='geo', dest='qd', height=0, datetime=None,):
    a = Apex()
    lat_a, lon_a = a.convert(lat, lon, source, dest, height, datetime, precision=1e-10, ssheight=50*6371)
    return lat_a, lon_a

def geodetic2geocentric(theta, alt):
    ###https://github.com/zzyztyy/pyIGRF
    geocentric_colatitude, gccolat_minus_theta, geocentric_radius = calculate.geodetic2geocentric(theta, alt)
    return geocentric_colatitude, gccolat_minus_theta

def mlng2mlt(ut, mlng):
    '''
     Relationship between MLT, geomagnetic longitude, and UT
        MLT=UT+(MLNG-70)/15

    Relationship between LT, geographic longitude, and UT
        LT=UT+LNG/15
    '''
    MLT = ut + (mlng - 70) / 15
"""