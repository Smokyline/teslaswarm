import numpy as np
from collections import OrderedDict

from tlib.CDate import datenum, datevecfull


def get_IAGA2002(name):   #   [header, kol]


    datag = []

    fid = open(name, 'rt')

    header = OrderedDict()
    header['Format'] = fid.readline()[24:-2].rstrip()
    header['SourceofData'] = fid.readline()[24:-2].rstrip()
    header['StationName'] = fid.readline()[24:-2].rstrip()
    header['IAGACODE'] = fid.readline()[24:-2].rstrip()
    header['GeodticLatitude'] = fid.readline()[24:-2].rstrip()
    header['GeodticLongitude'] = fid.readline()[24:-2].rstrip()
    header['Elevation'] = fid.readline()[24:-2].rstrip()
    header['Reported'] = fid.readline()[24:-2].rstrip()
    header['SensorOrientation'] = fid.readline()[24:-2].rstrip()
    header['DigitalSampling'] = fid.readline()[24:-2].rstrip()
    header['DataIntervalType'] = fid.readline()[24:-2].rstrip()
    header['DataType'] = fid.readline()[24:-2].rstrip()
    header['DataComment'] = np.zeros((20,), dtype=np.object)
    header['Date_begin'] = []
    header['Frequency'] = 0

    i = 0
    for str in fid:
        if str.find('DATE') >= 0:
            break
        str = str[:-2].rstrip()
        header['DataComment'][i] = str
        i += 1
    for j in range(i, 20):
        header['DataComment'][j] = ' '


    i = 0
    date0 = -1
    for str in fid:
        data = np.zeros((5, ))
        dt = str.split()
        if date0 == -1:
            dtmp = np.array(dt[0].split('-')).astype(int)
            date0 = datenum(dtmp[0], dtmp[1], dtmp[2])
        time = np.array(dt[1].split(':')).astype(float)
        time = (3600.0 * time[0] + 60.0 * time[1] + time[2]) / 86400
        data[0] = date0 + time
        data[1:] = np.array(dt[3:]).astype(float)
        datag.append(data)

    fid.close()
    datag = np.array(datag)

    return header, datag