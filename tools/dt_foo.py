import datetime
import time
import numpy as np
from datetime import timezone, timedelta
import calendar

def decode_str_dt_param(str_dt):
    date, time = str_dt.split('T')
    dt = datetime.datetime.strptime(
        '%s %s' %
        (date, time), '%Y-%m-%d %H:%M:%S')
    return dt

def decode_str_time_param(str_time):
    dt_time = datetime.datetime.strptime(str_time, '%H:%M:%S')
    return dt_time

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

def get_timestamp():
    # unix time ms
    return int(round(time.time() * 1000))

def get_local_time(utdate, uttime, lon):
    UTtime_h, UTtime_m, UTtime_s = np.array(uttime.split(':')).astype(int)
    loct = lon / 15 + UTtime_h + UTtime_m / 60 + UTtime_s / 3600
    if loct < 0:
        loct = 24 + loct
    if loct >= 24:
        loct = loct - 24

    local_time = datetime.datetime.strptime(
        '%s %s' %
        (utdate, '%i:%i:%i' %
         (loct, UTtime_m, UTtime_s)), '%Y-%m-%d %H:%M:%S')

    local_time = local_time.time()
    return local_time
