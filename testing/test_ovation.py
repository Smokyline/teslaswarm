import datetime
import numpy as np
import matplotlib.pyplot as plt
import math
from ovationpyme import ovation_prime, ovation_utilities
from geospacepy import satplottools, special_datetime
from teslaswarm.settings import STATIC_OS_PATH
from tools.data_foo import *
from engine.sword import SWORD
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.ticker import (LongitudeFormatter, LatitudeFormatter,
                                LatitudeLocator)
import matplotlib.ticker as mticker

def draw_weighted_flux(dt, atype='diff', jtype='energy'):
    """
    Test automatic generation of omni_intervals in ovation_utilities
    also by not specifying a start and end time for the FluxEstimator
    """
    estimator = ovation_prime.FluxEstimator(atype,jtype)

    mlatgridN, mltgridN, fluxgridN = estimator.get_flux_for_time(dt, hemi='N')
    #mlatgridS, mltgridS, fluxgridS = estimator.get_flux_for_time(dt, hemi='S')


    sword = SWORD('ortho_n')


    #XN, YN = latlt2cart(mlatgridN.flatten(), mltgridN.flatten(), 'N')
    XN, YN = satplottools.latlt2cart(mlatgridN.flatten(), mltgridN.flatten(), 'N')


    #XN, YN = satplottools.latlon2polar(XN, YN, 'N')


    #XN, YN = satplottools.latlt2cart(mlatgridN.flatten(), mltgridN.flatten(), 'N')

    #XS, YS = satplottools.latlt2cart(mlatgridS.flatten(), mltgridS.flatten(), 'S')
    #print('min Xn: {0} max Xn: {1}'.format(min(YN), max(YN)))
    #print('min Xs: {0} max Xs: {1}'.format( min(YS), max(YS)))

    #YN, XN = cart2pol(YN, XN)


    XN = np.array(XN).reshape(mlatgridN.shape)
    YN = np.array(YN).reshape(mltgridN.shape)
    #XS = XS.reshape(mlatgridS.shape)
    #YS = YS.reshape(mltgridS.shape)


    cmap = plt.get_cmap('seismic')
    sword.draw_colorbar(fluxgridN, cmap, 'auroral')

    rotated_pole = ccrs.RotatedPole(pole_longitude=0, pole_latitude=0)
    sword.ax.contourf(YN, XN, fluxgridN, cmap=cmap, transform=rotated_pole)

    #mappableN = ax.pcolormesh(XN, YN, fluxgridN, transform=ccrs.PlateCarree())
    #mappableS = ax.pcolormesh(XS, YS, fluxgridS, vmin=0, vmax=2, transform=ccrs.PlateCarree())


    im = sword.fig_to_PIL_image()
    im.save(STATIC_OS_PATH + '/media/images/ova_test.png')
    return im

dt1 = datetime.datetime(2017, 9, 8, 0, 00, 0)
draw_weighted_flux(dt1)