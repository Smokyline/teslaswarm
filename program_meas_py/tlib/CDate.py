import datetime
import numpy as np


def datenum(year=2014, month=1, day=1, hour=0, minute=0, seconds=0, smicro = 0):

    dt = datetime.date(year, month, day)

    num = dt.toordinal() + 366
    num += (3600.0*hour + 60.0*minute + seconds + smicro/1000000)/86400

    return num


def datevec(dt):

    dtt = datetime.date.fromordinal(int(dt-366))
    tm_year = dtt.year
    tm_mon = dtt.month
    tm_day = dtt.day
    ost = 86400*np.mod(dt, 1)
    tm_hour = int(np.floor(ost/3600))
    ost = np.mod(ost, 3600)
    tm_min = int(np.floor(ost/60))
    tm_sec = np.mod(ost, 60)
#    ostint = np.round(ost)
#    tm_hour = int(ostint/3600)
#    ostint = np.mod(ostint, 3600)
#    tm_min = int(ostint/60)
#    ostint = np.mod(ostint, 60)
#    tm_sec = float(ostint)
#    tm_micro = np.mod(ost, 1)*1000000

    return tm_year, tm_mon, tm_day, tm_hour, tm_min, tm_sec

def datevecfull(dt):

    dtt = datetime.date.fromordinal(int(dt-366))
    tm_year = dtt.year
    tm_mon = dtt.month
    tm_day = dtt.day
    ost = 86400*np.mod(dt, 1)
    ostint = np.round(ost)
    tm_hour = int(ostint/3600)
    ostint = np.mod(ostint, 3600)
    tm_min = int(ostint/60)
    ostint = np.mod(ostint, 60)
    tm_sec = int(ostint)
    tm_micro = np.mod(ost, 1)*1000000

    return tm_year, tm_mon, tm_day, tm_hour, tm_min, tm_sec, tm_micro