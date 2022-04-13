import pymysql
import datetime
from tools.data_foo import *
from tools.dt_foo import ut_dt_to_unix, unix_to_ut_dt

### !!!!!
import pandas as pd
from teslaswarm.settings import STATIC_OS_PATH

#  unix time to ut
# ut_time_dt = datetime.fromtimestamp(int_value)

#  ut time to unix
# unix_time = ut_time_dt.strftime('%s') # '1485896400'
# or
# unix_time = ut_time_dt..timestamp() # 1485896400.0

def get_sql_response(swarm_type, from_date, to_date, fac2_mod=False):
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
                      swarm_type, from_date, to_date)
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
                      swarm_type, from_date, to_date)
    print(request)
    cur = sql_connect.cursor()
    cur.execute(request)
    respond = cur.fetchall()
    sql_connect.close()
    return respond

    #df = pd.read_csv(STATIC_OS_PATH + '/data/dataset_fac2_3days.csv', header=None)
    #return df


def get_swarm_set(swarm_liter, date_from, date_to, delta, fac2_mod=False):
    """Парсинг ответа SQL в [[swarm liter, XY, dt, values], ...]"""
    from_date = ut_dt_to_unix(date_from, out_type='str')
    to_date = ut_dt_to_unix(date_to, out_type='str')

    respond = get_sql_response('SW' + swarm_liter.upper(), from_date, to_date, fac2_mod)
    if len(np.array(respond)) == 0:
        print('MISSING VALUES IN SQL TABLE!!!')
    # -> dt, y, x, r, n, e, c, f
    respond = data_reduction(np.array(respond), delta, fac2_mod)

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

    print('-------------------------------------------------------------------------------------')
    print('The SWORD render will use data from swarm-%s, datetime from %s to %s' % (
    swarm_liter, swarm_date[0] + ' ' + swarm_time[0],
    swarm_date[-1] + ' ' + swarm_time[-1]))
    return [swarm_liter, np.array(swarm_position),
            np.array(swarm_date), np.array(swarm_time), swarm_values]
