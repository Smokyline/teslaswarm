import numpy as np
import pymysql
import warnings
import datetime
from datetime import timezone

def get_sql_response(swarm_liter, from_date, to_date, fac2_mod=False):
    swarm_liter = 'SW' + swarm_liter.upper()
    from_date = ut_dt_to_unix(from_date, out_type='str')
    to_date = ut_dt_to_unix(to_date, out_type='str')
    #   подключение к sql таблице на сервере imagdb.gcras.ru
    sql_connect = pymysql.connect(host='imagdb.gcras.ru', port=3306,
                                  user='data_reader',
                                  passwd='V9e2OeroufrluwluN2u88leP9lAPhluzL7dlU67oAb5afROub9iUv7unLEhiEs9o',
                                  db='intermagnet',
                                  charset='utf8mb4',
                                  )
    #   получение данных swarm_type, from_date, to_date
    if fac2_mod:
        request = "SELECT date, latitude, longitude, radius, fac FROM sat_sec_fac_plain WHERE code='%s' " \
                  "AND date BETWEEN %s AND %s" % (
                      swarm_liter, from_date, to_date)
        # a % c SW_
        """
        SELECT `date`, latitude, longitude, radius, fac 
        FROM sat_sec_fac_plain WHERE code='SW_' 
        AND date BETWEEN UNIX_TIMESTAMP('2015-1-10 00:00:59') AND UNIX_TIMESTAMP('2015-1-10 05:59:59')
        """
    else:
        #request = "SELECT date, latitude, longitude, radius, n, e, c FROM sat_sec_plain WHERE code='%s' " \
        request = "SELECT date, latitude, longitude, radius, n, e, c, f FROM sat_sec_plain WHERE code='%s' " \
                  "AND date BETWEEN %s AND %s" % (
                      swarm_liter, from_date, to_date)
    print('SQL REQUEST: ', request)
    cur = sql_connect.cursor()
    cur.execute(request)
    respond = cur.fetchall()
    sql_connect.close()
    print('The SWORD render will use data from swarm-%s, datetime from %s to %s' % (
        swarm_liter, str(unix_to_ut_dt(int(respond[0][0]))).split(' '),  str(unix_to_ut_dt(int(respond[-1][0]))).split(' ')))
    print('-------------------------------------------------------------------------------------')
    return np.array(respond)


def get_swarm_set(respond, swarm_liter, delta=None, fac2_mod=False):
    """Парсинг ответа SQL в [[swarm liter, XY, dt, values], ...]"""


    # -> dt, y, x, r, n, e, c, f
    if delta is not None:
        respond = data_averaging(np.array(respond), delta, fac2_mod)

    swarm_date = []
    swarm_time = []
    swarm_position = []

    if fac2_mod:
        # fac2 = (Y, X, R), dt, (fac2)
        swarm_dt_unix, Y, X, R, FAC = respond.T
        for i, dt in enumerate(swarm_dt_unix):
            sw_date, sw_time = str(unix_to_ut_dt(int(dt))).split(' ')
            swarm_date.append(sw_date)
            swarm_time.append(sw_time)
            swarm_position.append([Y[i], X[i], R[i] / 1000])  # Y X
        swarm_values = np.empty((len(swarm_date), 0))
        for v in [FAC]:
            swarm_values = np.append(swarm_values, np.array([v]).T, axis=1)
    else:
        # chaos = (Y, X, R), dt, (N, E, C, F)
        # vector, measure mu, fac = (Y, X), dt, (N, E, C, F)
        swarm_dt_unix, Y, X, R, N, E, C, F = respond.T
        for i, dt in enumerate(swarm_dt_unix):
            sw_date, sw_time = str(unix_to_ut_dt(int(dt))).split(' ')
            swarm_date.append(sw_date)
            swarm_time.append(sw_time)
            swarm_position.append([Y[i], X[i], (R[i] / 1000)])  # km => r/1000  Y X R

        swarm_values = np.empty((len(swarm_date), 0))
        for v in [N, E, C, F]:
            swarm_values = np.append(swarm_values, np.array([v]).T, axis=1)
    if swarm_liter == '_':
        swarm_liter = 'AC'

    return [swarm_liter, np.array(swarm_position),
            np.array(swarm_date), np.array(swarm_time), swarm_values]


def data_averaging(respond, delta, fac2_mod=False):
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

def decode_str_dt_param(str_dt):
    date, time = str_dt.split('T')
    dt = datetime.datetime.strptime(
        '%s %s' %
        (date, time), '%Y-%m-%d %H:%M:%S')
    return dt

def ut_dt_to_unix(dt, out_type='str'):
    """Конвертирование UT datetime в unix time"""
    if out_type == 'str':
        #unix_time = int(time.mktime(dt.timetuple()))
        #unix_time = int(dt.timestamp())
        #unix_time = int(dt.replace(tzinfo=timezone.utc).timestamp())
        unix_time = dt.replace(tzinfo=timezone.utc).astimezone(tz=None).timestamp()

    elif out_type == 'float':
        unix_time = dt.timestamp()
    else:
        print('Error input unix type')
        unix_time = None
    return unix_time  # int ot float


def unix_to_ut_dt(unix_time):
    """Конвертирование unix int time в UT datetime"""

    if not isinstance(unix_time, int):
        print('Error input unix type')
        ut_time_dt = None
    else:
        #ut_time_dt = datetime.datetime.fromtimestamp(unix_time)
        ut_time_dt = datetime.datetime.utcfromtimestamp(unix_time).strftime('%Y-%m-%d %H:%M:%S')
    return ut_time_dt  # dt


def get_sql_swarm_data(sw_liter='A',
                       from_date='2015-02-17T03:23:00', to_date='2015-02-17T03:27:00',
                       delta=1, fac2_mod=False):
    """
    input
    sw_liter: A, B, C, _ (A&C)
    from_date, to_date: YY-MM-DDTh:m:s
    delta: int, sec average swarm values
    fac2_mod: bool

    return
    swarm_liter,
    swarm_latGEO, swarm_lonGEO,  swarm_r_km,
    swarm_datetime (datetime) UTC
    swarm_values (X,Y,Z,F) or (FAC) if fac2_mod==True)
    """
    from_date, to_date = decode_str_dt_param(from_date), decode_str_dt_param(to_date)



    sql_response = get_sql_response(sw_liter, from_date, to_date, fac2_mod)
    swarm_set = get_swarm_set(sql_response, sw_liter, delta, fac2_mod)
    sw_liter = swarm_set[0]
    sw_lat = swarm_set[1][:, 0]
    sw_lon = swarm_set[1][:, 1]
    sw_rad = swarm_set[1][:, 2]
    sw_dt = np.array([decode_str_dt_param(d+'T'+t) for d, t in zip(swarm_set[2], swarm_set[3])])
    sw_vector = swarm_set[4]    # [x,y,z,f] for [fac] if fac2_mod==True
    return sw_liter, sw_lat, sw_lon, sw_rad, sw_dt, sw_vector

print(get_sql_swarm_data())