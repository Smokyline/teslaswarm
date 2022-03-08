from io import StringIO

import numpy as np
from datetime import datetime
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap


def sun_pos(dt=None):
    """This function computes a rough estimate of the coordinates for
    the point on the surface of the Earth where the Sun is directly
    overhead at the time dt. Precision is down to a few degrees. This
    means that the equinoxes (when the sign of the latitude changes)
    will be off by a few days.
    The function is intended only for visualization. For more precise
    calculations consider for example the PyEphem package.
    Parameters
    ----------
        dt: datetime
            Defaults to datetime.utcnow()
    Returns
    -------
        lat, lng: tuple of floats
            Approximate coordinates of the point where the sun is
            in zenith at the time dt.
    """
    if dt is None:
        dt = datetime.utcnow()

    axial_tilt = 23.4
    ref_solstice = datetime(2020, 1, 1, 0, 0)
    days_per_year = 365.2425
    seconds_per_day = 24 * 60 * 60.0

    days_since_ref = (dt - ref_solstice).total_seconds() / seconds_per_day
    lat = axial_tilt * np.cos(2 * np.pi * days_since_ref / days_per_year)
    sec_since_midnight = (dt - datetime(dt.year, dt.month, dt.day)).seconds
    lng = -(sec_since_midnight / seconds_per_day - 0.5) * 360
    return lat, lng


def fill_dark_side(sw_date, sw_time):
    """
    Plot a fill on the dark side of the planet (without refraction).
    Parameters
    ----------
        ax : matplotlib axes
            The axes to plot on.
        time : datetime
            The time to calculate terminator for. Defaults to datetime.utcnow()
        **kwargs :
            Passed on to Matplotlib's ax.fill()
    """
    time = datetime.strptime('%s %s' % (sw_date, sw_time), '%Y-%m-%d %H:%M:%S')
    lat, lng = sun_pos(time)
    pole_lng = lng
    if lat > 0:
        pole_lat = -90 + lat
        central_rot_lng = 180
    else:
        pole_lat = 90 + lat
        central_rot_lng = 0

    rotated_pole = ccrs.RotatedPole(pole_latitude=pole_lat,
                                    pole_longitude=pole_lng,
                                    central_rotated_longitude=central_rot_lng)

    x = np.empty(360)
    y = np.empty(360)
    x[:180] = -90
    y[:180] = np.arange(-90, 90.)
    x[180:] = 90
    y[180:] = np.arange(90, -90., -1)

    return x, y, rotated_pole
