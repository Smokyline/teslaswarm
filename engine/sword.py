from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.gridspec as gridspec
import matplotlib.patches as patches
from mpl_toolkits.axes_grid1 import make_axes_locatable

import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.ticker import (LongitudeFormatter, LatitudeFormatter,
                                LatitudeLocator)
from cartopy.io.shapereader import Reader
from cartopy.feature import ShapelyFeature
from pyproj import Proj, transform
from mpl_toolkits import axisartist
from mpl_toolkits.axes_grid1 import host_subplot
import matplotlib.ticker as mticker
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import io
from PIL import Image

from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from scipy.interpolate import griddata
from teslaswarm.settings import STATIC_OS_PATH
from tools.dt_foo import decode_str_dt_param, decode_str_time_param
from tools.data_foo import get_mag_coordinate_system_lines, geo2mag, mag2mlt, get_obs_value_by_time

class SWORD():
    """
    SWarm Ordinary RenDer
    """
    def __init__(self, proj_type='', extend=None):
        #   custom projection
        #   https://matplotlib.org/stable/gallery/misc/custom_projection.html#sphx-glr-gallery-misc-custom-projection-py
        self.proj_type = proj_type

        # Orthographic or Miller image projection
        if proj_type == 'ortho_n':
            projection = ccrs.Orthographic(central_latitude=90., central_longitude=0)
            # projection = ccrs.NorthPolarStereo(central_longitude=0)
            fig = plt.figure(figsize=(22, 22))  # a*dpi x b*dpi aka 3000px x 3000px
            ax = fig.add_subplot(1, 1, 1, projection=projection)
            ax.set_facecolor((1.0, 1.0, 1.0))
            ax.coastlines(resolution='110m', color='k', zorder=7)
            ax.add_feature(cfeature.LAND, facecolor='0.75', zorder=0)
            gl = ax.gridlines(draw_labels=True)
            gl.xlocator = mticker.FixedLocator([-180, -150, -120, -90, -60, -30, 0,
                                                30, 60, 90, 120, 150, 180])
            gl.ylocator = mticker.FixedLocator(np.arange(-90, 91, 10))
            gl.xformatter = LongitudeFormatter()
            gl.yformatter = LatitudeFormatter()
            gl.ylabel_style = {'size': 12, 'color': 'gray'}
            gl.xlabel_style = {'size': 12, 'color': 'gray'}
            ax.set_global()
            global_extent = ax.get_extent(crs=ccrs.PlateCarree())
            ax.set_extent(global_extent[:2] + (50, 90), crs=ccrs.PlateCarree())
            # ax.set_extent([xlims[0], xlims[1], ylims[0], ylims[1]], crs=ccrs.Geodetic())
            #ax.set_extent([-180, 180, 60, 90], projection)     #  lon_min, lon_max, lat_min, lat_max,
            self.ax = ax

        if proj_type == 'ortho_s':
            projection = ccrs.Orthographic(central_latitude=-90., central_longitude=0)
            # projection = ccrs.NorthPolarStereo(central_longitude=0)
            fig = plt.figure(figsize=(20, 20))  # a*dpi x b*dpi aka 3000px x 3000px
            ax = fig.add_subplot(1, 1, 1, projection=projection)
            ax.set_facecolor((1.0, 1.0, 1.0))
            ax.coastlines(resolution='110m', color='k', zorder=7)
            ax.add_feature(cfeature.LAND, facecolor='0.75', zorder=0)
            gl = ax.gridlines(draw_labels=True)
            gl.xlocator = mticker.FixedLocator([-180, -150, -120, -90, -60, -30, 0,
                                                30, 60, 90, 120, 150, 180])
            gl.ylocator = mticker.FixedLocator(np.arange(-90, 91, 10))
            gl.xformatter = LongitudeFormatter()
            gl.yformatter = LatitudeFormatter()
            gl.ylabel_style = {'size': 12, 'color': 'gray'}
            gl.xlabel_style = {'size': 12, 'color': 'gray'}
            ax.set_global()
            global_extent = ax.get_extent(crs=ccrs.PlateCarree())
            ax.set_extent(global_extent[:2] + (-50, -90), crs=ccrs.PlateCarree())
            # ax.set_extent([xlims[0], xlims[1], ylims[0], ylims[1]], crs=ccrs.Geodetic())
            self.ax = ax

        if proj_type == 'miller':
            projection = ccrs.Miller()
            #projection = ccrs.PlateCarree()
            # projection = ccrs.NorthPolarStereo(central_longitude=0)
            fig = plt.figure(figsize=(20, 20))  # a*dpi x b*dpi aka 3000px x 3000px
            ax = fig.add_subplot(1, 1, 1, projection=projection)
            ax.set_facecolor((1.0, 1.0, 1.0))
            ax.coastlines(resolution='50m', color='k', zorder=7)
            ax.add_feature(cfeature.LAND, facecolor='0.75', zorder=0)
            #ax.stock_img()
            gl = ax.gridlines(draw_labels=True)
            gl.xlocator = mticker.FixedLocator([-180, -150, -120, -90, -60, -30, 0,
                                                30, 60, 90, 120, 150, 180])
            gl.ylocator = mticker.FixedLocator(np.arange(-90, 91, 10))
            gl.xformatter = LongitudeFormatter()
            gl.yformatter = LatitudeFormatter()
            gl.ylabel_style = {'size': 12, 'color': 'gray'}
            gl.xlabel_style = {'size': 12, 'color': 'gray'}
            ax.set_global()
            self.ax = ax

        # увеличение конкретной области
        if extend is not None:
            self.ax.set_extent(extend, ccrs.PlateCarree())

        else:
            """no projection"""
            pass

        #plt.gca().set_axis_off()
        plt.subplots_adjust(top=1, bottom=0.025, right=1, left=0.025,
                            hspace=.025, wspace=.025)
        plt.margins(0, 0)
        #plt.gca().xaxis.set_major_locator(plt.NullLocator())
        #plt.gca().yaxis.set_major_locator(plt.NullLocator())

    def draw_plot(self, draw_list, include=None, draw_obs=None):
        # long_set
        short_set, long_set = False, False
        if len(draw_list) == 1:
            short_set = True
        elif len(draw_list) > 1:
            long_set = True

        #   отрисовка графика
        if short_set:
            fig = plt.figure(figsize=(20, 10))
            #axs = host_subplot(111, axes_class=axisartist.Axes)
            axs = fig.add_subplot(111)
            axs.set_aspect('auto')
        elif long_set:
            fig, axs = plt.subplots(len(draw_list), figsize=(20, 20))

        for i, sw_set in enumerate(draw_list):
            label, date_list, time_list, value, pos = sw_set
            if short_set:
                ax1 = axs
                ax2 = ax1.twiny()
                ax3 = ax2.twinx()

            elif long_set:
                ax1 = axs[i]
                ax2 = ax1.twiny()
                ax3 = ax2.twinx()

            """if len(np.array(sw_set[4]).shape) > 1:
                Z = np.array(sw_set[4])[:, i]
            else:
                Z = np.array(sw_set[4])"""
            #   time ticks
            #str_ticks = [decode_dt_param(d + 'T' + t).strftime('%c') for d, t in zip(date_list, time_list)]
            str_ticks = time_list
            date_ticks = np.array([[x, label] for x, label in zip(np.arange(len(time_list)), str_ticks)])

            magcoord_pos = np.array(np.round(mag2mlt(geo2mag(pos, date_list, to_coord='GSM')[:, :2], date_list, time_list), 2)).astype(str)
            #magcoord_pos = np.array(np.round(mag2mlt(pos[:, :2], date_list, time_list), 2)).astype(str)
            magcoord_pos = np.array([lat + ' ' + lon for lat, lon in magcoord_pos]).astype(str)
            coord_pos = np.array([str(lat) + ' ' + str(lon) for lat, lon in np.round(pos[:, :2], 2)])

            #pos_ticks = np.array([[x, label] for x, label in zip(np.arange(len(time_list)), coord_pos)])
            #magpos_ticks = np.array([[x, label] for x, label in zip(np.arange(len(time_list)), magcoord_pos)])

            coord_ticks = []
            for mag, geo in zip(magcoord_pos, coord_pos):
                coord_ticks.append(str(mag))
                coord_ticks.append('\n\n'+str(geo))
            coord_ticks = np.array(coord_ticks).astype(str)
            coord_ticks = np.array([''.join(x) for x in zip(coord_ticks[0::2], coord_ticks[1::2])])
            coordtick_locations = np.arange(len(pos))


            # decode from str to datetime format
            points_time = np.array([decode_str_dt_param(d + 'T' + t) for d, t in zip(date_list, time_list)], dtype='datetime64')
            last_point = points_time[0]
            line_idx = []
            for k, point_time in enumerate(points_time):
                minute_before = point_time - np.timedelta64(1, 'm')
                if minute_before == last_point:
                    last_point = point_time
                    line_idx.append(k)
            if len(date_ticks) >= 12:
                date_ticks = date_ticks[line_idx]
                coord_ticks_arange = np.arange(0, len(coord_ticks), int(len(coord_ticks) / 14))
            else:
                coord_ticks_arange = np.arange(0, len(coord_ticks))
            ax1.vlines(x=line_idx, ymin=np.min(value), ymax=np.max(value), color='b', linestyle='--', alpha=.5)


            ax1.plot(np.arange(len(time_list)), value, c='r', label=label)
            if include[i] is not None:
                ax1.plot(np.arange(len(include[i])), include[i], c='b', label='auroral')
            if draw_obs is not None:
                ax3.set_ylabel("%s ground station, nT" % draw_obs)
                ch_str = label.split(' ')[1]
                channels = ['n', 'e', 'c']
                try:
                    ch_where_nec = channels.index(ch_str)
                    obs_value = get_obs_value_by_time(date_list, time_list, channel=ch_where_nec, obs=draw_obs)
                    ax3.plot(np.arange(len(obs_value)), obs_value, label='%s %s' % (draw_obs, ch_str))

                except:
                    if ch_str == 'fac2':
                        for chnl in [0, 1, 2]:     # n, e, c
                            obs_value = get_obs_value_by_time(date_list, time_list, channel=chnl, obs=draw_obs)
                            ax3.plot(np.arange(len(obs_value)), obs_value, label='%s %s' % (draw_obs, channels[chnl]))


            ax1.annotate("MLat MLT\n\nLat Lon", xy=(-.05, 1.03), xycoords="axes fraction")
            ax1.set_ylabel(label)
            ax1.set_xticks(np.array(date_ticks[:, 0]).astype(int))
            ax1.set_xticklabels(date_ticks[:, 1])

            ax2.set_xlim(ax1.get_xlim())
            ax2.set_xticks(coordtick_locations[coord_ticks_arange])
            ax2.set_xticklabels(coord_ticks[coord_ticks_arange])

            #https://stackoverflow.com/questions/20532614/multiple-lines-of-x-tick-labels-in-matplotlib

            ax1.grid()
            ax1.xaxis.grid(False)
            ax2.grid()
            ax1.legend()
            ax3.legend(loc=4)

        dt_from, dt_to = decode_str_dt_param(draw_list[0][1][0]+'T'+draw_list[0][2][0]), decode_str_dt_param(draw_list[0][1][-1]+'T'+draw_list[0][2][-1])
        fig.suptitle('from %s to %s' % (dt_from, dt_to), fontsize=16)
        plt.grid(True)


    def draw_mag_coord_lines(self, date, geomag_pole=False):
        mag_lat_lines, mag_lon_lines, annotate_points = get_mag_coordinate_system_lines(date, geomag_pole)

        for lat_lines in mag_lat_lines:
            self.ax.plot(lat_lines[:, 1], lat_lines[:, 0], 'r:', zorder=1, lw=1, alpha=0.75,
                         transform=ccrs.Geodetic())
        for lon_lines in mag_lon_lines:
            self.ax.plot(lon_lines[:, 1], lon_lines[:, 0], 'r:', zorder=1, lw=1, alpha=0.75,
                         transform=ccrs.Geodetic())
        for i, xy in enumerate(annotate_points):
            if self.proj_type == 'ortho_n':
                if xy[2][0] > 44:
                    self.ax.text(xy[1], xy[0], '%i,%i' % (xy[2][0], xy[2][1]), transform=ccrs.PlateCarree(), fontsize=6,
                                 color='r')
            if self.proj_type == 'ortho_s':
                if xy[2][0] < -29:
                    self.ax.text(xy[1], xy[0], '%i,%i' % (xy[2][0], xy[2][1]), transform=ccrs.PlateCarree(), fontsize=6,
                                 color='r')
            else:
                self.ax.text(xy[1], xy[0], '%i,%i' % (xy[2][0], xy[2][1]), transform=ccrs.PlateCarree(), fontsize=6, color='r')


    def draw_swarm_scatter(self, points, values, channel=0, annotate=False, custom_value_name=None):
        """отрисовка """
        # prepare colorbar, based on values
        if custom_value_name is not None:
            label = custom_value_name
            if custom_value_name == "Dd":
                cmap = plt.get_cmap('jet')
                #cmap = plt.get_cmap('seismic')
            else:
                cmap = plt.get_cmap('seismic')

        elif channel is None:
            label = 'SWARM ' + "fac2"
            if len(np.array(values).shape) > 1:
                values = values[:, channel]
            cmap = plt.get_cmap('seismic')

        else:
            label_set = ['dBn', 'dBe', 'dBd', 'dB']
            label = 'SWARM ' + label_set[channel]
            cmap = plt.get_cmap('seismic')

        self.draw_colorbar(values=values, cmap=cmap, label=label, vcenter=0.)


        # draw scatter point, color based on values and legend
        levels = MaxNLocator(nbins=10).tick_values(values.min(), values.max())
        norm = BoundaryNorm(levels, cmap.N, clip=True)
        s = np.power(len(values), -1/16) * 300
        self.ax.scatter(points[:, 1], points[:, 0], s=s, c=values, cmap=cmap, norm=norm, transform=ccrs.PlateCarree(), edgecolors='k', lw=0.3)
        #self.ax.scatter(points[:, 1], points[:, 0], c=values, cmap=cmap, norm=norm, transform=ccrs.PlateCarree())

        #str_values = np.format_float_positional(values, precision=3)
        str_values =[str('%.2f' % x) for x in values]
        #self.ax.text(points[:, 1], points[:, 0], str_values, )
        #self.ax.scatter(points[:, 1], points[:, 0])

        if annotate:
            for i, txt in enumerate(str_values):
                # self.ax.annotate(txt, (points[i, 1], points[i, 0]), transform=ccrs.PlateCarree())
                self.ax.text(points[i, 1], points[i, 0], txt, transform=ccrs.PlateCarree(), fontsize=8)

    def draw_swarm_path(self, points, points_time=None):
        x, y = points[:, 1], points[:, 0]
        self.ax.plot(x, y,'k--',zorder=1,lw=1,alpha=0.4,transform=ccrs.Geodetic())
        line_length = 1

        if points_time is not None:
            # decode from str to datetime format
            points_time = np.array(points_time, dtype='datetime64')

            last_point = points_time[0]
            last_annotate_point = points_time[0]
            for k, point_time in enumerate(points_time):
                minute_before = point_time - np.timedelta64(1, 'm')
                if minute_before == last_point:
                    try:
                        x1, y1 = x[k], y[k]
                        x2, y2 = x[k + 1], y[k + 1]
                        # (x1, y1) + a(y2-y1, x1-x2)
                        p1 = np.array((x1, y1))
                        a1 = 1
                        a2 = -1
                        l1 = np.array([y2 - y1, x1 - x2])
                        l1_length = np.sqrt(l1[0]**2+l1[1]**2)
                        l1 = l1/l1_length
                        l2 = np.array((a1 * l1[0], a1 * l1[1]))
                        l3 = np.array((a2 * l1[0], a2 * l1[1]))
                        c1 = p1 + l2
                        c2 = p1 + l3
                        line = Line2D([c1[0], x1, c2[0]],
                                      [c1[1], y1, c2[1]], transform=ccrs.Geodetic())
                        self.ax.add_line(line, )
                        self.ax.text(c1[0], c1[1], str(point_time).split('T')[1], fontsize=8,  transform=ccrs.Geodetic())

                    except Exception as e:
                        pass
                    last_point = point_time

                five_minute_before = point_time - np.timedelta64(5, 'm')
                if five_minute_before == last_annotate_point:
                    try:
                        x1, y1 = x[k], y[k]
                        x2, y2 = x[k + 1], y[k + 1]

                        # (x1, y1) + a(y2-y1, x1-x2)
                        p1 = np.array((x1, y1))
                        a1 = 1
                        a2 = -1

                        l1 = np.array([y2 - y1, x1 - x2])
                        l1_length = np.sqrt(l1[0] ** 2 + l1[1] ** 2)
                        l1 = l1 / l1_length

                        l2 = np.array((a1 * l1[0], a1 * l1[1]))
                        l3 = np.array((a2 * l1[0], a2 * l1[1]))

                        c1 = p1 + l2
                        c2 = p1 + l3

                        self.ax.text(c1[0], c1[1], str(point_time).split('T')[1], fontsize=8,  transform=ccrs.Geodetic())
                    except Exception as e:
                        pass
                    last_annotate_point = point_time

        """for i, xy in enumerate(points):
            try:
                pass
                #self.ax.annotate(s='', xy=(xy[i + 1], xy[i + 1]), xytext=(xy[i], xy[i]), arrowprops=dict(arrowstyle='->'), transform=ccrs.PlateCarree())
                #self.ax.annotate(s='', xy=(xy[i + 1], xy[i + 1]), xytext=(xy[i], xy[i]), transform=ccrs.PlateCarree())
            except Exception as e:
                pass"""

    def draw_vector(self, sw_pos, B):
        X, Y = sw_pos[:, 1], sw_pos[:, 0]
        #B = sw_value[:, :2]
        u = 5000
        q = self.ax.quiver(X, Y, B[:, 0], B[:, 1], transform=ccrs.PlateCarree(), width=0.0008)
        self.ax.quiverkey(q, X=0.7, Y=0.99, U=u, labelpos='E', label='{N, E} vector length = %s' % u, transform=ccrs.PlateCarree())


    def draw_ionomodel(self, surf_data, type=['north','seismic']):
        x, y, z = surf_data[:, 0], surf_data[:, 1], surf_data[:, 2]
        # grid the data.
        xi = np.linspace(x.min()-1, x.max()+1, 200)
        yi = np.linspace(y.min(), y.max(), 200)
        Vi = griddata((x, y), z, (xi[None, :], yi[:, None]), method='linear')  # create a uniform spaced grid
        X, Y = np.meshgrid(xi, yi)

        lin = np.max(np.abs(z)) + (np.max(np.abs(z)) / 10)
        if type[1] =='sigma':
            cmap = 'cool'
            bounds = np.linspace(-0.5, lin, 50)
        else:
            cmap = 'PiYG'
            bounds = np.linspace(-1 * lin, lin, 80)

        self.draw_colorbar(Vi, cmap, 'ionomodel ' + str(type[0]) + ' μА/m²')
        self.ax.contourf(X, Y, Vi, transform=ccrs.PlateCarree(), cmap=cmap)

    def draw_avroral_oval(self, surf_data, hemisphere='north'):
        """x, y, z = surf_data[:, 0], surf_data[:, 1], surf_data[:, 2]
        # grid the data.
        xi = np.linspace(x.min()-1, x.max()+1, 200)
        yi = np.linspace(y.min(), y.max(), 200)
        Vi = griddata((x, y), z, (xi[None, :], yi[:, None]), method='nearest')  # create a uniform spaced grid
        X, Y = np.meshgrid(xi, yi)

        self.draw_colorbar(z, cmap='coolwarm', label='avroral egrs/cm2s')
        self.ax.contourf(X, Y, Vi, transform=ccrs.PlateCarree(), cmap='coolwarm')"""


        X, Y, flux_value = surf_data
        cmap_args = self.draw_colorbar(flux_value, cmap='plasma', label='auroral %s egrs/cm2s' % hemisphere, vmin=0.5, )
        #cmap_args = self.draw_colorbar(flux_value, cmap=aurora_cmap, label='auroral %s egrs/cm2s' % hemisphere, vmin=0.5, )
        """if hemisphere == 'north':
            rotated_pole = ccrs.RotatedPole(pole_longitude=0, pole_latitude=90)
        else:
            rotated_pole = ccrs.RotatedPole(pole_longitude=0, pole_latitude=90)"""
            #rotated_pole = ccrs.RotatedPole(pole_longitude=0, pole_latitude=0)
        rotated_pole = ccrs.PlateCarree()
        self.ax.contourf(X, Y, flux_value, **cmap_args, transform=rotated_pole, alpha=0.85)


    def draw_shapefile(self, shapefile):
        #inProj = Proj(init='laea', preserve_units=True)
        #outProj = Proj(init='epsg:4326')
        lccproj  = ccrs.LambertConformal(central_latitude=45,
                              central_longitude=100,
                              standard_parallels=(25, 25))
        for shape in shapefile.shapeRecords():
            for i in range(len(shape.shape.parts)):
                i_start = shape.shape.parts[i]
                if i == len(shape.shape.parts) - 1:
                    i_end = len(shape.shape.points)
                else:
                    i_end = shape.shape.parts[i + 1]
                x = [i[0] for i in shape.shape.points[i_start:i_end]]
                y = [i[1] for i in shape.shape.points[i_start:i_end]]
                #print(x, y)
                #x, y = transform(inProj,outProj,x,y)
                #print(x)
                self.ax.plot(x, y, transform=lccproj, zorder=1)


    def draw_point_with_annotate(self, points):
        #   for intermag observ for example
        for point in points:
            x, y, annotate = float(point[0]), float(point[1]), str(point[2])
            self.ax.scatter(x, y, c='r', marker='*', transform=ccrs.PlateCarree())
            self.ax.text(x, y, annotate, transform=ccrs.PlateCarree(), fontsize=10)

    def draw_colorbar(self, values, cmap, label, vcenter=None, vmin=None):
        if vmin is not None:
            cmap_edit = plt.get_cmap(cmap)
            cmap_edit.set_under('white', alpha="0")
            #cmap_edit.set_under((1, 1, 1, 0))
            sm = plt.cm.ScalarMappable(cmap=cmap_edit, norm=plt.Normalize(vmin=vmin, vmax=values.max()))
            cmap_args = dict(cmap=cmap_edit, vmin=vmin, vmax=values.max())
        elif vcenter is not None:
            if np.abs(values.min()) > np.abs(values.max()):
                vmin, vmax = -1 * np.abs(values.min()), np.abs(values.min())
            else:
                vmin, vmax = -1 * np.abs(values.max()), np.abs(values.max())
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=vmin, vmax=vmax))
            cmap_args = None
        else:
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=values.min(), vmax=values.max()))
            cmap_args = dict(cmap=cmap, vmin=values.min(), vmax=values.max())


        sm._A = []
        if self.proj_type == 'miller':
            rotation, orientation = 0, 'horizontal'
        else:
            rotation, orientation = 270, 'vertical'
        #   https://stackoverflow.com/questions/15908371/matplotlib-colorbars-and-its-text-labels
        cb = plt.colorbar(sm, ax=self.ax, shrink=0.4, orientation=orientation, fraction=0.046, pad=0.015)
        cb.set_label(label, labelpad=-40, rotation=rotation)
        return cmap_args

    def draw_polygon(self, poly):
        x, y = poly.exterior.xy
        self.ax.plot(x, y, transform=ccrs.PlateCarree())

    def set_axis_label(self, label, zoom_axis=False):
        print('axis label:   ', label)
        if zoom_axis:
            self.ax.set_title(label, fontsize=20, y=-0.075)
        else:
            self.ax.set_title(label, fontsize=20)

    def fig_to_PIL_image(self, name='1'):

        buf = io.BytesIO()
        plt.savefig(buf, format='jpeg', dpi=300)
        buf.seek(0)
        im = Image.open(buf)
        im.save(STATIC_OS_PATH +'/media/tmp/%s.jpg' % name)
        buf.close()
        return im

