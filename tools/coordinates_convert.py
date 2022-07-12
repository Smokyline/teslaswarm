import math
import numpy as np
from tools.sun_position import sun_pos
from tools.rotate2 import rotate
from spacepy import coordinates as coord
from spacepy.time import Ticktock
import datetime

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



def rotate_GEO_vector_to_MFA(sw_lat, sw_lon, sw_rad, sw_Bvector, sw_dtime):
    from geopack import geopack, t89

    '''
    convert from GEO b vector to mean field alignet (MFA) coord system
    https://geo.phys.spbu.ru/~tsyganenko/empirical-models/coordinate_systems/gsw_system
    https://github.com/tsssss/geopack
    '''
    t0 = datetime.datetime(1970, 1, 1)
    rzero = 1+(100/6371.0)
    #xyzGEO = to_cast(sw_lon, sw_lat, sw_rad/6371.0)
    ut = (sw_dtime[0] - t0).total_seconds()
    ps = geopack.recalc(ut)

    coord_GSM = np.empty(shape=(0, 3))
    field_MFA = np.empty(shape=(0, 3))
    field_MFA_lines = []
    for i in range(len(sw_lat)):


        cvals = coord.Coords([sw_rad[i]/6371.0, np.float(sw_lat[i]), np.float(sw_lon[i])], 'GEO', 'sph', ['Re', 'deg', 'deg'])
        cvals.ticks = Ticktock(sw_dtime[i])
        xgeo,ygeo,zgeo = np.array(cvals.convert('GEO', 'car').data)[0]
        #xgeo,ygeo,zgeo = geopack.sphcar(sw_rad[i]/6371.2, np.deg2rad(sw_lat[i]), np.deg2rad(sw_lon[i]), 1)


        xgsm, ygsm, zgsm = geopack.geogsm(xgeo,ygeo,zgeo,  1)

        # TODO 1 north, -1 south
        x1gsm, y1gsm, z1gsm, xx,yy,zz = geopack.trace(xgsm, ygsm, zgsm,
                                            dir=1, rlim=10, r0=rzero,
                                            parmod=2,exname='t89',inname='igrf', maxloop=1)
        dxgsm, dygsm, dzgsm = x1gsm - xgsm, y1gsm-y1gsm, z1gsm-zgsm

        xnew = rotate(v1=[xgsm, ygsm, zgsm], v2=[dxgsm, dygsm, dzgsm], xold=sw_Bvector[i, :3])
        field_MFA = np.append(field_MFA, [xnew], axis=0)

        #print(xgsm, ygsm, zgsm)
        #print(dxgsm, dygsm, dzgsm)
        #rgeo, latgeo, longeo = geopack.sphcar(xgsm, ygsm, zgsm, -1)
        cvals = coord.Coords([xgsm, ygsm, zgsm], 'GSM', 'car')
        cvals.ticks = Ticktock(sw_dtime[i])
        rgeo, latgeo, longeo = np.array(cvals.convert('GSM', 'sph').data)[0]

        coord_GSM = np.append(coord_GSM, [[rgeo, latgeo, longeo]], axis=0)
        #coord_GSM = np.append(coord_GSM, [[xgsm, ygsm, zgsm]], axis=0)
        field_MFA_lines.append([xx,yy,zz])
        """

        b0xgsm,b0ygsm,b0zgsm = geopack.dip(xgsm,ygsm,zgsm)    		# calc dipole B in GSM.
        dbxgsm,dbygsm,dbzgsm = t89.t89(2, ps, xgsm,ygsm,zgsm)       # calc T89 dB in GSM.
        #bxgsm,bygsm,bzgsm = [b0xgsm+dbxgsm,b0ygsm+dbygsm,b0zgsm+dbzgsm]
        
        bcvals = coord.Coords([bxgsm,bygsm,bzgsm], 'GSM', 'car')
        bcvals.ticks = Ticktock(dt)
        r, lat, lon = np.array(cvals.convert('GEO', 'shp').data)[0]
        field_var_geocoord = np.append(field_var_geocoord, [[lat, lon, r]], axis=0)"""
    return np.array(field_MFA),field_MFA_lines, coord_GSM

def GEO2MAG(lat, lon, dt):
    cvals = coord.Coords([np.float(6371.0), np.float(lat), np.float(lon)], 'GEO', 'sph', ['Re', 'deg', 'deg'])
    cvals.ticks = Ticktock(dt)
    new_coord = np.array(cvals.convert('MAG', 'sph').data)[0]
    return [new_coord[2], new_coord[1]]     # mag lon, lat

def MAG2GEO(lat, lon, dt):
    cvals = coord.Coords([np.float(6371.0), np.float(lat), np.float(lon)], 'MAG', 'sph', ['Re', 'deg', 'deg'])
    cvals.ticks = Ticktock(dt)
    new_coord = np.array(cvals.convert('GEO', 'sph').data)[0]
    return [new_coord[2], new_coord[1]]     # mag lon, lat

def geo2GSM(sw_pos, swarm_date, to_coord='MAG'):
    # input sw_pos format: lat, lon, rad

    # coordinate init_coords = [Z, Y, X] or [km from center of earth, lat, lon]
    init_coords = np.array([sw_pos[:, 2], sw_pos[:, 0], sw_pos[:, 1]]).T
    # alt above sea lvl to rad
    for i, r in enumerate(init_coords[:, 0]):
        earth_r = earth_radius_in_meters(init_coords[i, 1])
        init_coords[i, 0] = (init_coords[i, 0] + earth_r) * 1000 / earth_r
    cvals = coord.Coords(init_coords, 'GEO', 'sph')

    dt_array = []
    for k, sw_date in enumerate(swarm_date):
        dt_array.append(str(sw_date) + 'T' + '00:00:00')
    # print(dt_array)
    cvals.ticks = Ticktock(dt_array, 'UTC')

    newcoord = np.array(cvals.convert(to_coord, 'sph').data)  # return radius, latitude, longitude

    for i, r in enumerate(newcoord[:, 0]):
        earth_r = earth_radius_in_meters(init_coords[i, 1])
        newcoord[i, 0] = newcoord[i, 0] * earth_r / 1000
    # lan, lon, r
    return np.array([newcoord[:, 1], newcoord[:, 2], newcoord[:, 0]]).T



def geomag2geo(latlon, PHI, THETA):
    switched_swarm_pos = []
    for lat, lon in latlon:
        lat_a, lon_a = geomag2geo_foo(lat, lon, PHI, THETA)
        switched_swarm_pos.append([lat_a, lon_a])
    return np.array(switched_swarm_pos)

def geomag2geo_foo(glat: float, glon: float, PHI, THETA) :
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
    #Dlong = 360 - (180 - 107.4)
    Dlong = THETA
    Dlat = PHI

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



def cart2polar(x, y):
    rho = np.sqrt(x ** 2 + y ** 2)
    phi = np.arctan2(y, x)
    return (rho, np.rad2deg(phi))








"""

def sph_to_cart(lon, lat):
    r = 6371.2
    desc_data = np.empty((0, 3))
    for i in range(len(lon)):
        # az, el = sph_x[i] * math.pi / 180, sph_y[i] * math.pi / 180
        az, el = math.radians(lon[i]), math.radians(lat[i])     # lon, lat
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

def cart_to_geo(data):
    spher_data = np.empty((0, 2))
    for x, y, z in data:
        radius = np.sqrt(x ** 2 + y ** 2 + z ** 2)
        latitude = np.arcsin(z / radius)
        longitude = np.arctan2(y, x)
        spher_data = np.append(spher_data, [[math.degrees(latitude), math.degrees(longitude)]], axis=0)
    return spher_data[:, 0], spher_data[:, 1]   # lat, lon


def pol2cart(rho, phi):
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return [x, y]

def to_cast(sph_x, sph_y):
    # geodetic_to_geocentric
    r = 6371.2
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




def latlt2cart(lat, lt, hemisphere):
    ''''''
    Latitude and local time to cartesian for a top-down dialplot
    ''''''
    r, theta = latlt2polar(lat, lt, hemisphere)

    return r * np.cos(theta), r * np.sin(theta)
    #return r, theta


def latlt2polar(lat, lt, hemisphere):
    ''''''
    Converts an array of latitude and lt points to polar for a top-down dialplot (latitude in degrees, LT in hours)
    i.e. makes latitude the radial quantity and MLT the azimuthal

    get the radial displacement (referenced to down from northern pole if we want to do a top down on the north,
        or up from south pole if visa-versa)
    ''''''
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


def mag2mlt(sw_pos, swarm_date, swarm_time):
    ''''''
        magnetic longitude to MLT
        When referencing this package, please cite both the package DOI and the AACGM-v2 journal article:

        Shepherd, S. G. (2014), Altitude‐adjusted corrected geomagnetic coordinates: Definition and functional approximations,
        Journal of Geophysical Research: Space Physics, 119, 7501–7521, doi:10.1002/2014JA020264.

    ''''''
    convert_mlt = []
    for pos, date, time in zip(sw_pos, swarm_date, swarm_time):
        mlt_coords = np.array(aacgmv2.convert_mlt(pos, datetime.datetime.strptime(date+' '+time, '%Y-%m-%d %H:%M:%S'), m2a=False))
        convert_mlt.append(mlt_coords)

    return np.array(convert_mlt)



def geo_to_gg(radius, theta):


    ''''''

from geocentric radius and colatitude.
    to geodetic colatitude and vertical height above the ellipsoid

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
    ''''''

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




Relationship between MLT, geomagnetic longitude, and UT
MLT=UT+(MLNG-70)/15
Relationship between LT, geographic longitude, and UT
LT=UT+LNG/15


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