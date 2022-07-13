import datetime

import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.gridspec as gridspec
import matplotlib.patches as patches
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.cm as cm
import matplotlib.colors as mcolors
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
from tools.sun_position import sun_pos
from tools.data_foo import *

class SWORD():
    """
    SWarm Ordinary RenDer
    """
    def __init__(self, proj_type='', extend=None):
        #   custom projection
        #   https://matplotlib.org/stable/gallery/misc/custom_projection.html#sphx-glr-gallery-misc-custom-projection-py
        self.proj_type = proj_type
        self.transform = ccrs.PlateCarree()

        # Orthographic or Miller image projection
        if proj_type == 'ortho_n':
            self.projection = ccrs.Orthographic(central_latitude=90., central_longitude=0)
            # projection = ccrs.NorthPolarStereo(central_longitude=0)
            self.fig = plt.figure(figsize=(15, 15))  # a*dpi x b*dpi aka 3000px x 3000px
            ax = self.fig.add_subplot(1, 1, 1, projection=self.projection)
            ax.set_facecolor((1.0, 1.0, 1.0))
            ax.coastlines(resolution='110m', color='k', zorder=1)
            ax.add_feature(cfeature.LAND, facecolor='0.75', zorder=0)
            gl = ax.gridlines(draw_labels=True)
            gl.xlocator = mticker.FixedLocator(np.arange(-180, 180, 15))
            gl.ylocator = mticker.FixedLocator(np.arange(-90, 91, 10))
            gl.xformatter = LongitudeFormatter()
            gl.yformatter = LatitudeFormatter()
            gl.ylabel_style = {'size': 11, 'color': 'gray'}
            gl.xlabel_style = {'size': 11, 'color': 'gray'}
            ax.set_global()
            global_extent = ax.get_extent(crs=ccrs.PlateCarree())
            ax.set_extent(global_extent[:2] + (50, 90), crs=ccrs.PlateCarree())
            # ax.set_extent([xlims[0], xlims[1], ylims[0], ylims[1]], crs=ccrs.Geodetic())
            #ax.set_extent([-180, 180, 60, 90], projection)     #  lon_min, lon_max, lat_min, lat_max,
            self.ax = ax

        if proj_type == 'ortho_s':
            self.projection = ccrs.Orthographic(central_latitude=-90., central_longitude=0)
            # projection = ccrs.NorthPolarStereo(central_longitude=0)
            self.fig = plt.figure(figsize=(15, 15))  # a*dpi x b*dpi aka 3000px x 3000px
            ax = self.fig.add_subplot(1, 1, 1, projection=self.projection)
            ax.set_facecolor((1.0, 1.0, 1.0))
            ax.coastlines(resolution='110m', color='k', zorder=1)
            ax.add_feature(cfeature.LAND, facecolor='0.75', zorder=0)
            gl = ax.gridlines(draw_labels=True)
            gl.xlocator = mticker.FixedLocator(np.arange(-180, 180, 15))
            gl.ylocator = mticker.FixedLocator(np.arange(-90, 91, 10))
            gl.xformatter = LongitudeFormatter()
            gl.yformatter = LatitudeFormatter()
            gl.ylabel_style = {'size': 11, 'color': 'gray'}
            gl.xlabel_style = {'size': 11, 'color': 'gray'}
            ax.set_global()
            global_extent = ax.get_extent(crs=ccrs.PlateCarree())
            ax.set_extent(global_extent[:2] + (-50, -90), crs=ccrs.PlateCarree())
            # ax.set_extent([xlims[0], xlims[1], ylims[0], ylims[1]], crs=ccrs.Geodetic())
            self.ax = ax

        if proj_type == 'miller':
            self.projection = ccrs.Miller()
            #projection = ccrs.PlateCarree()
            # projection = ccrs.NorthPolarStereo(central_longitude=0)
            self.fig = plt.figure(figsize=(15, 15))  # a*dpi x b*dpi aka 3000px x 3000px
            ax = self.fig.add_subplot(1, 1, 1, projection=self.projection)
            ax.set_facecolor((1.0, 1.0, 1.0))
            ax.coastlines(resolution='50m', color='k', zorder=1)
            ax.add_feature(cfeature.LAND, facecolor='0.75', zorder=0)
            #ax.stock_img()
            gl = ax.gridlines(draw_labels=True)
            gl.xlocator = mticker.FixedLocator(np.arange(-180, 181, 15))
            gl.ylocator = mticker.FixedLocator(np.arange(-90, 91, 10))
            gl.xformatter = LongitudeFormatter()
            gl.yformatter = LatitudeFormatter()
            gl.ylabel_style = {'size': 11, 'color': 'gray'}
            gl.xlabel_style = {'size': 11, 'color': 'gray'}
            ax.set_global()
            self.ax = ax

        # увеличение конкретной области
        if (extend is None and proj_type == 'miller') or proj_type == 'plot':
            print('proj')
            self.extend = None
        else:
            if extend is None:
                if proj_type == 'ortho_s':
                    self.extend = [-179, 179, -40, -90]
                elif proj_type == 'ortho_n':
                    self.extend = [-179, 179, 40, 90]
            else:
                self.extend = extend
            self.ax.set_extent(self.extend, self.transform)
            self.ax.set_aspect('equal')

        #plt.gca().set_axis_off()
        """plt.subplots_adjust(top=1, bottom=0.025, right=1, left=0.025,
                            hspace=.025, wspace=.025)
        plt.margins(0, 0)"""
        #plt.gca().xaxis.set_major_locator(plt.NullLocator())
        #plt.gca().yaxis.set_major_locator(plt.NullLocator())

    def draw_plot(self, draw_list, auroral_list, delta=1, draw_station=None):
        station_name = [k for k in draw_station.keys()][0]
        obs_value = draw_station[station_name]

        sat_plot_color = {'A': 'g', 'B': 'k', 'C': 'b', 'AC': 'm', '|AC|': 'c'}
        plot_count = len(draw_list)
        if len(obs_value) > 0:
            plot_count += 3
        fig, axs = plt.subplots(ncols=1, nrows=plot_count, figsize=(15, 5*plot_count),
                                constrained_layout=True)
        liters = [sw_set[0].split('-')[1].split()[0] for i, sw_set in enumerate(draw_list)]
        #j2_list_pos = np.where(liters == 'AC')
        j2_list_pos = [ i for i in range(len(liters)) if liters[i]=='AC']
        if len(j2_list_pos)> 0:
            j2_list_pos = j2_list_pos[0]
        else:
            j2_list_pos = None
        for plot_num, sw_set in enumerate(draw_list):
            label, date_list, time_list, value, pos = sw_set
            sw_liter = label.split('-')[1].split()[0]
            if len(draw_list) == 1 and len(obs_value) == 0:
                host = axs
            elif len(draw_list) == 1 and len(obs_value) > 0:
                host = axs[0]
            elif len(draw_list) > 1:
                host = axs[plot_num]
            ax2 = host.twiny()
            #   main value plot
            if sw_liter != 'AC':
                host.plot(np.arange(len(time_list)), value, c=sat_plot_color[sw_liter], label=label, lw=2)
            if ((j2_list_pos != None and sw_liter in ['A', 'C']) or sw_liter == 'AC') and label.split('-')[1].split(' ')[1] == 'fac':
                AC_value, AC_label = draw_list[j2_list_pos][3],  'swarm-AC (J₂) ' + draw_list[j2_list_pos][0].split('-')[1].split(' ')[1]
                AC_value_upper_zero = AC_value.copy()
                AC_value_upper_zero[np.where(AC_value < 0)] = np.nan
                host.plot(np.arange(len(time_list)), AC_value, c='k', label=AC_label + ' (< 0)', lw=2)
                host.plot(np.arange(len(time_list)), AC_value_upper_zero, c='m', label=AC_label + ' (>= 0)', lw=2)

            # prepare dt ticks
            dt_tick = []
            for i, d in enumerate(date_list):
                # solar_lat, solar_lon = sun_pos(dt)
                # solarcoord_pos.append('%s  %.1f' % (lt.strftime('%H:%M:%S'), 90. - solar_lat))
                dt = decode_str_dt_param(date_list[i] + 'T' + time_list[i])
                lon = pos[i, 1]
                lt = get_local_time(dt, lon)
                dt_tick.append(dt.strftime('%H:%M:%S') + ' UT')
                dt_tick.append('\n' +lt.strftime('%H:%M:%S') + ' LT')
            dt_tick = np.array(dt_tick).astype(str)
            dt_tick = np.array([''.join(x) for x in zip(dt_tick[0::2], dt_tick[1::2])])
            date_ticks = np.array([[x, label] for x, label in zip(np.arange(len(time_list)), dt_tick)])

            # prepare GEO coord ticks (top)
            # pos_ticks = np.array([[x, y_label] for x, y_label in zip(np.arange(len(time_list)), coord_pos)])
            # magpos_ticks = np.array([[x, y_label] for x, y_label in zip(np.arange(len(time_list)), magcoord_pos)])
            coord_pos = np.array([[lat, lon] for lat, lon in np.round(np.array(pos[:, :2]).astype(float), 2)])
            coord_ticks = []
            for lat, lon in coord_pos:
                coord_ticks.append(str(lat))
                coord_ticks.append('\n' + str(lon))
            coord_ticks = np.array(coord_ticks).astype(str)
            coord_ticks = np.array([''.join(x) for x in zip(coord_ticks[0::2], coord_ticks[1::2])])

            # decode from str to datetime format
            points_time = np.array([decode_str_dt_param(d + 'T' + t) for d, t in zip(date_list, time_list)],
                                   dtype='datetime64')
            if len(date_ticks) >= 12:
                tick_count = 9
                if points_time[-1] - points_time[0] < np.timedelta64(1, 'h'):  # отрисовка vlines раз в минуту
                    coord_ticks_arange = np.arange(0, len(coord_ticks),
                                                   int(len(coord_ticks) / tick_count))  # top xlabels tick
                    last_point = points_time[0]
                    line_idx = []
                    for k, point_time in enumerate(points_time):
                        if points_time[-1] - points_time[0] > np.timedelta64(30, 'm'):
                            minute_before = point_time - np.timedelta64(5, 'm')
                        else:
                            minute_before = point_time - np.timedelta64(1, 'm')
                        if minute_before == last_point:
                            last_point = point_time
                            line_idx.append(k)
                    date_ticks = date_ticks[line_idx]
                    #nonnan_value = [x for x in value if not np.isnan(x)]
                    #host.vlines(x=line_idx, ymin=np.min(nonnan_value), ymax=np.max(nonnan_value), color='b',linestyle='--', alpha=.5, linewidths=1)
                else:
                    coord_ticks_arange = np.arange(0, len(coord_ticks),
                                                   int(len(coord_ticks) / tick_count))  # top xlabels tick
                    date_ticks = date_ticks[coord_ticks_arange]

            else:
                coord_ticks_arange = np.arange(0, len(coord_ticks))
            coord_ticks[0] = ''

            right_axes = []

            if len(auroral_list) > 0:
                auroral_diff = auroral_list[plot_num][0]
                auroral_mono = auroral_list[plot_num][1]
                right_axes.append(host.twinx)
                #print(len(auroral_list), 'i=', plot_num)
                ax4 = host.twinx()
                ax4.set_ylabel("OVATION auroral oval near %s, ergs/cm²" % label)
                ax4.plot(np.arange(len(auroral_diff)), auroral_diff, c='b',
                         label='auroral oval diff near ' + label, lw=2)
                ax4.plot(np.arange(len(auroral_mono)), auroral_mono, c='y',
                         label='auroral oval mono near ' + label, lw=2)
                ax4.legend(loc=2)
                ax4.grid(linestyle='-')
                ax4.set_zorder(host.get_zorder() + 1)

            if 'fac' in label:
                label += ' %sA/m²' % chr(956)
            else:
                label += ' nT'
            host.set_ylabel(label)
            host.set_xticks(np.array(date_ticks[:, 0]).astype(int))
            host.set_xticklabels(date_ticks[:, 1], rotation=40)
            ax2.grid(False)
            ax2.set_xlim(host.get_xlim())
            ax2.set_xticks(np.arange(len(pos))[coord_ticks_arange])
            ax2.set_xticklabels(coord_ticks[coord_ticks_arange])
            ax2.tick_params(axis='x', labelsize=11, grid_alpha=0)
            ax2.annotate("GEO lat\nGEO lon", xy=(-0.04, 1.037), xycoords=ax2.transAxes, size=11)

            #https://stackoverflow.com/questions/20532614/multiple-lines-of-x-tick-labels-in-matplotlib
            #ax2.xaxis.grid(False)
            #ax2.grid(False)
            host.grid(linestyle='--')
            host.legend()

        # plot superMag station value
        if len(obs_value) > 0:
            station_plot_color = {'X': 'g', 'Y': 'k', 'Z': 'b'}
            label, date_list, time_list, value, pos = draw_list[0]

            try:
                obs_loc = get_superMAG_observ_loc(station_name)[:2]
                p_in_p = get_position_near_point(pos, obs_loc, degr_radius=10)
                switch_idxs = []
                last_bool_value = p_in_p[0]
                for i, b in enumerate(p_in_p):
                    if b != last_bool_value:
                        switch_idxs.append(i)
                        last_bool_value = b
            except :
                switch_idxs = []
            for col, vect_comp in enumerate(['X', 'Y', 'Z']):
                ax = axs[col - 3]
                ax.set_ylabel("SuperMAG ground station %s, nT" % station_name)
                ax.plot(np.arange(len(obs_value)), obs_value[:, col], lw=2, c=station_plot_color[vect_comp], label='%s %s, nT' % (station_name, vect_comp))

                if len(switch_idxs)>=2:
                    ax.vlines(x=switch_idxs, ymin=np.min(obs_value[:, col]), ymax=np.max(obs_value[:, col]), color='k',
                              linestyle='--', alpha=.7, linewidths=1)
                    window_idx = switch_idxs[0] + int((switch_idxs[1] - switch_idxs[0]) / 2)
                    center_y_axis = np.min(obs_value[:, col]) +  (np.max(obs_value[:, col])-np.min(obs_value[:, col]))/2
                    style = dict(size=10, color='gray')
                    ax.text(window_idx, center_y_axis, "SWARM flight window", ha='center', **style)
                ax.set_xlim(axs[0].get_xlim())
                ax.set_xticks(np.array(date_ticks[:, 0]).astype(int))
                ax.set_xticklabels(date_ticks[:, 1], rotation=40)
                ax.legend(loc=4)
                ax.grid(linestyle='--')

        if draw_list[0][1][0] != draw_list[0][1][-1]:  # first date != second date
            dt_from, dt_to = decode_str_dt_param(draw_list[0][1][0] + 'T' + draw_list[0][2][0]), decode_str_dt_param(
                draw_list[0][1][-1] + 'T' + draw_list[0][2][-1])
            fig.suptitle('from %s to %s' % (dt_from, dt_to), fontsize=16)
        else:
            fig.suptitle(draw_list[0][1][0], fontsize=16)

        plt.grid(True)

    def draw_swarm_path(self, points, points_time=None, annotate=False):
        ax_proj_xyz = self.ax.projection.transform_points(self.transform, points[:, 1], points[:, 0], points[:, 2] )
        X, Y, Z = ax_proj_xyz[:, 0], ax_proj_xyz[:, 1], ax_proj_xyz[:, 1]
        x_deg, y_deg = points[:, 1], points[:, 0]

        # swarm path -----------------
        if self.proj_type != 'miller':
            self.ax.plot(points[:, 1], points[:, 0], 'k--', zorder=8, lw=1, alpha=0.4, transform=ccrs.Geodetic())
        else:
            #self.ax.plot(X, Y, 'k--', zorder=8, lw=1, alpha=0.4, transform=self.ax.projection)
            self.ax.plot(points[:, 1], points[:, 0], 'k--', zorder=8, lw=1, alpha=0.4, transform=ccrs.Geodetic())

        #self.ax.plot(X, Y, 'k--', zorder=8, lw=1, alpha=0.4, transform=ccrs.Geodetic())

        if points_time is not None and len(points_time) > 0:
            print('draw vline')
            # decode from str to datetime format
            points_time = np.array(points_time, dtype='datetime64')

            last_point = points_time[0]
            for k, point_time in enumerate(points_time):
                minute_before = point_time - np.timedelta64(1, 'm')
                if minute_before >= last_point:
                    # if minute_before == last_point:
                    try:
                        if k < len(X) - 1:
                            x1, y1, z1 = X[k], Y[k], Z[k]
                            x2, y2, z2 = X[k + 1], Y[k + 1], Z[k + 1]
                            target_x, target_y = (x2 - x1), (y2 - y1)

                        else:
                            x1, y1, z1 = X[k], Y[k], Z[k]
                            x2, y2, z2 = X[k - 1], Y[k - 1], Z[k - 1]
                            target_x, target_y = (x1 - x2), (y1 - y2)
                        if self.proj_type != 'miller':
                            m = 70000/np.sqrt(target_x ** 2 + target_y ** 2)
                        else:
                            m = 1 / np.sqrt(target_x ** 2 + target_y ** 2)
                        target_x, target_y = target_x * m, target_y * m
                        #target_x, target_y = (x2 - x1) * 7, (y2 - y1) * 7


                        text = point_time.astype(datetime.datetime).strftime('%H:%M') + ' UT'

                        #self.ax.arrow(x1, y1, (x2 - x1) * 7, (y2 - y1) * 7, head_width=0.75, head_length=0.25, zorder=5, alpha=.65)
                        #self.ax.arrow(x1, y1, target_x, target_y, head_width=0.75, head_length=0.25, zorder=5,
                        self.ax.arrow(x1, y1, target_x, target_y, head_width=0.75, head_length=0.25, zorder=5, alpha=.65)

                        extend_bool = self.extend is not None
                        annotate_and_long_range = annotate == False and points_time[-1] - points_time[0] < np.timedelta64(1, 'D')
                        if annotate_and_long_range or annotate:
                            if extend_bool:
                                if self.extend[0] < x_deg[k] < self.extend[1] and self.extend[2] < y_deg[k] < self.extend[3]:
                                    self.ax.text(x1, y1, text, fontsize=12, clip_on=True, zorder=5, transform=self.ax.projection)
                                # self.ax.annotate(text, xytext=(x2, y2), xy=(x1, y1), size=10, fontsize=12)
                            else:
                                self.ax.text(x1, y1, text, fontsize=6, clip_on=True, zorder=5, transform=self.ax.projection)
                        elif extend_bool:
                            # D M  H:M  LT
                            text = point_time.astype(datetime.datetime).strftime('%d-%m %H:%M') + '\nLT:%s' % get_local_time(point_time.astype(datetime.datetime),points[k, 1]).strftime('%H:%M')
                            self.ax.text(x1, y1, text, fontsize=5, clip_on=True, zorder=5, transform=self.ax.projection)

                    #else:
                        #    self.ax.text(x1, y1, text, fontsize=5, clip_on=True, zorder=5, transform=self.ax.projection)

                    except Exception as e:
                        print(e)
                    last_point = point_time


    def draw_swarm_scatter(self, points, values, custom_label=None, annotate=False):
        values = np.array(values).astype(np.float).flatten()
        print(custom_label, 'custom label')
        X, Y, Z = points[:, 1], points[:, 0], points[:, 2]
        """отрисовка """
        # prepare colorbar, based on values
        cb_type = 'zero_center'
        unit = 'no unit'
        cmap = plt.get_cmap('jet')
        if str(custom_label).rsplit(' ')[1] in ['N', 'E', 'C', 'F']:
            unit = 'SWARM-' + custom_label + ', nT'
            cb_type = 'vmin_vmax'
        if '|AC|' in custom_label:
            cb_type = 'zero_min'
        if 'FAC' in custom_label:
            unit = 'SWARM-%s, %sA/m²' % (custom_label, chr(956))
            cmap = plt.get_cmap('seismic')
        if 'IGRF' in custom_label:
            unit = custom_label + ' (nT)'
        if 'CHAOS7' in custom_label:
            unit = custom_label + ' (nT)'
        if 'mu' in custom_label:
            unit += ' DMA anomaly lvl'
            cb_type = 'DMA_anomaly'
            #cmap = plt.get_cmap('RdYlBu')
            c = mcolors.ColorConverter().to_rgb
            cmap = mpl.colors.ListedColormap([c('#1F2FFF'), c('#00FF4E'), c('#FFCD0D'), c('#FF1C19')])



        print('scatter max:%s scatter min:%s' % (np.max(values), np.min(values)), )
        nonnan_values = np.array([x for x in values if not pd.isnull(x)])
        levels = [np.min(nonnan_values), -100.0, -75., -50., -25., 0, 25., 50., 75., 100., np.max(nonnan_values)]
        # draw scatter point, color based on values and legend
        cmap_args = self.draw_colorbar(values=nonnan_values, cmap=cmap, label=unit, cb_type=cb_type)
        m = cm.ScalarMappable(norm=cmap_args['norm'], cmap=cmap_args['cmap'])
        m.set_array(levels)
        if 'FAC' in custom_label:
            if self.proj_type == 'miller':
                line_mp_lenght = 20
            else:
                line_mp_lenght = 1e6
        elif 'd' in custom_label:
            if self.proj_type == 'miller':
                line_mp_lenght = 10
            else:
                line_mp_lenght = 75e4
        else:
            if self.proj_type == 'miller':
                line_mp_lenght = 3
            else:
                line_mp_lenght = 1e5

        if self.proj_type == 'miller':
            line_mp_lenght_max = 1500
        else:
            line_mp_lenght_max = 750000

        value_max = np.nanmax(values)
        a_length = []
        for k in range(len(X)):
            if k < len(X) - 1:
                x1, y1, z1 = X[k], Y[k], Z[k]
                x2, y2, z2 = X[k + 1], Y[k + 1], Z[k + 1]
            else:
                x1, y1, z1 = X[k], Y[k], Z[k]
                x2, y2, z2 = X[k - 1], Y[k - 1], Z[k - 1]

            x1, y1, z1 = self.ax.projection.transform_points(self.transform, np.array([x1]), np.array([y1]), )[0]
            x2, y2, z2 = self.ax.projection.transform_points(self.transform, np.array([x2]), np.array([y2]), )[0]
            p1 = np.array((x1, y1))
            a1 = values[k] / value_max * line_mp_lenght
            a_length.append(a1)
            a2 = -a1
            l1 = np.array([y2 - y1, x1 - x2])
            l1_length = np.sqrt(l1[0] ** 2 + l1[1] ** 2)
            l1 = l1 / l1_length

            l2 = np.array((a1 * l1[0], a1 * l1[1]))
            l3 = np.array((a2 * l1[0], a2 * l1[1]))
            l2[np.where(l2>line_mp_lenght_max)] = line_mp_lenght_max
            l3[np.where(l3<-line_mp_lenght_max)] = -line_mp_lenght_max

            c1 = p1 + l2
            c2 = p1 + l3

            line = Line2D([c1[0], x1, c2[0]],
                          [c1[1], y1, c2[1]], color=m.to_rgba(values[k]), alpha=.7, zorder=4)
            # [c1[1], y1, c2[1]], transform=transf, color='k', alpha=.9)
            self.ax.add_line(line)
        print('line length min:%s max:%s ' % (np.min(a_length), np.max(a_length)) )
        # solar coords
        if annotate:
            str_values = [str('%.2f' % x) for x in values]
            for i, txt in enumerate(str_values):
                # self.ax.annotate(txt, (points[i, 1], points[i, 0]), transform=self.transform)
                if self.extend is not None:
                    if self.extend[0] < points[i, 1] < self.extend[1] and\
                            self.extend[2] < points[i, 0] < self.extend[3]:
                        self.ax.text(points[i, 1], points[i, 0], txt, transform=self.ax.projection, fontsize=8)
                else:
                    self.ax.text(points[i, 1], points[i, 0], txt, transform=self.ax.projection, fontsize=8)


    def draw_ionomodel(self, surf_data, type=['north','seismic']):
        if type[1] =='sigh':
            cmap = plt.get_cmap('plasma')
            cb_type = 'zero_min'
            grid_method = 'linear'
            unit = ' Sm'
        else:
            #cmap = plt.get_cmap('PiYG')
            cmap = plt.get_cmap('coolwarm')
            cb_type = 'zero_center'
            grid_method = 'nearest'
            if type[1] == 'pot':
                unit = ' kV'
                type[1] = 'potential'
            else:
                unit = ' μА/m²'
        x, y, z = surf_data[:, 0], surf_data[:, 1], surf_data[:, 2]
        # grid the data.
        print(y.min(), y.max(), 'yyyy')
        print(x.min(), x.max(), 'xxx')
        xi = np.linspace(0, 360.3, 100)
        """if type[0] == 'northern':
            yi = np.linspace(40.4, 89.59, 125)
        else:
            yi = np.linspace(-89.59, -40.4, 125)"""
        #xi = np.linspace(x.min(), x.max() + 1, 100)
        if y.max()<0:
            yi = np.linspace(y.max(), -90., 75)
        else:
            yi = np.linspace(y.min(), 90., 75)
        yi = np.linspace(y.min(), y.max(), 100)

        Vi = griddata((x, y), z, (xi[None, :], yi[:, None]), method='nearest')  # create a uniform spaced grid
        X, Y = np.meshgrid(xi, yi)
        print('Vi max', np.max(z))
        print('Vi min', np.min(z))

        lin = np.max(np.abs(z)) + (np.max(np.abs(z)) / 10)

        cmap_args = self.draw_colorbar(z, cmap, 'ionomodel %s %s' % (type[0], type[1]) +unit, cb_type=cb_type)
        self.ax.contourf(X, Y, Vi, cmap=cmap, norm=cmap_args['norm'], transform=ccrs.PlateCarree(), alpha=0.75)

    def draw_avroral_oval(self, surf_data, hemisphere='north', atype='diff'):
        """x, y, z = surf_data[:, 0], surf_data[:, 1], surf_data[:, 2]
        # grid the data.
        xi = np.linspace(x.min()-1, x.max()+1, 200)
        yi = np.linspace(y.min(), y.max(), 200)
        Vi = griddata((x, y), z, (xi[None, :], yi[:, None]), method='nearest')  # create a uniform spaced grid
        X, Y = np.meshgrid(xi, yi)

        self.draw_colorbar(z, cmap='coolwarm', label='avroral egrs/cm2s')
        self.ax.contourf(X, Y, Vi, transform=ccrs.PlateCarree(), cmap='coolwarm')"""
        cb_type = 'zero_min_oval'

        X, Y, flux_value = surf_data
        #print(flux_value[np.isnan(flux_value)])
        #nonnan_values = flux_value[not np.isnan(flux_value)]

        flux_value[np.isnan(flux_value)] = 0
        cmap_args = self.draw_colorbar(flux_value, cmap=plt.get_cmap('cool'),
                                       label='OVATION auroral oval %s %s_hemi ergs/cm²' % (atype ,hemisphere), cb_type=cb_type)
        #cmap.set_under('white', 0)

        levels = np.linspace(0.25, 4, 20)
        #cmap_args = self.draw_colorbar(flux_value, cmap=aurora_cmap, label='auroral %s egrs/cm2s' % hemisphere, vmin=0.5, )
        """if hemisphere == 'north':
            rotated_pole = ccrs.RotatedPole(pole_longitude=0, pole_latitude=90)
        else:
            rotated_pole = ccrs.RotatedPole(pole_longitude=0, pole_latitude=90)"""
            #rotated_pole = ccrs.RotatedPole(pole_longitude=0, pole_latitude=0)
        rotated_pole = ccrs.PlateCarree()
        #self.ax.contourf(X, Y, flux_value, **cmap_args, transform=rotated_pole, alpha=0.85)
        self.ax.contourf(X, Y, flux_value, cmap=plt.get_cmap('cool'), levels=levels, transform=rotated_pole, alpha=0.85)
        #self.ax.contourf(X, Y, flux_value, cmap=cmap, levels=levels, transform=rotated_pole, alpha=0.85)


    def draw_colorbar(self, values, cmap, label, cb_type='zero_center'):
        print(label, cb_type, values.max(), values.min())

        if cb_type == 'zero_center':
            # scatter
            #levels = MaxNLocator(nbins=10).tick_values(values.min(), values.max())
            #norm = BoundaryNorm(levels, cmap.N, clip=True)

            if np.abs(values.min()) > np.abs(values.max()):
                vmin, vmax = -1 * np.abs(values.min()), np.abs(values.min())
            else:
                vmin, vmax = -1 * np.abs(values.max()), np.abs(values.max())
            if ('FAC' in label):
                if vmax < 5:
                    bounds = [-6, -5.0, -4., -3., -2., -1., 0, 1., 2., 3., 4., 5., 6]
                else:
                    bounds = [vmin, -5.0, -4., -3., -2., -1., 0, 1., 2., 3., 4., 5., vmax]
                #bounds = [vmin, 0, vmax]
            elif ('d' in label):
                if vmax < 50:
                    bounds = [-100, -50, -40, -30, -20, -10, 0, 10, 20, 30, 40, 50, 100]
                else:
                    bounds = [vmin, -50, -40, -30, -20, -10, 0, 10, 20, 30, 40, 50, vmax]
            else:
                delta = vmax/5
                bounds = [x for x in np.arange(vmin, vmax+delta, delta)]
            norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
            cm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
            #cm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=-10, vmax=10))
            #cm = plt.cm.ScalarMappable(cmap=cmap)
            cmap_args = dict(cmap=cmap, norm=norm)
        elif cb_type == 'DMA_anomaly':
            #bounds = [-1.0, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5, 0.75, 1.0]
            cmap.set_under('white', 0)
            bounds = [-1.0, -0.5, 0, 0.5, 1]
            norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
            cm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
            cmap_args = dict(cmap=cmap, norm=norm)
        elif 'zero_min' in cb_type:
            """
            cmap_edit = plt.get_cmap(cmap)
            cmap_edit.set_under('white', alpha="0")
            cm = plt.cm.ScalarMappable(cmap=cmap_edit, norm=plt.Normalize(vmin=0.25, vmax=np.max(values)))
            cmap_args = dict(cmap=cmap_edit, vmin=0.25, vmax=np.max(values))"""
            #colors = [(0, 1, 0), (1, 0, 0)]
            #lincm = LinearSegmentedColormap.from_list(
            #    "Custom", colors, N=20)
            if 'oval' in cb_type:
                cm = LinearSegmentedColormap.from_list(
                    'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=0.45, b=1.), cmap(np.linspace(0.45, 1., 100)))
                cm.set_under('white', 0)
            else:
                cmap = plt.get_cmap('jet')
            norm = plt.Normalize(vmin=0, vmax=values.max())
            cm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
            #cm.set_under('white', 0)
            cmap_args = dict(cmap=cmap, norm=norm)

        else:
            vmin, v_max = values.min(), values.max()
            #vmin, v_max = values[np.where(values>vmin)].min(), values[np.where(values<v_max)].max()
            norm = plt.Normalize(vmin=vmin, vmax=v_max)
            cm = plt.cm.ScalarMappable(cmap='jet', norm=norm)
            # cm = plt.cm.ScalarMappable(cmap=cmap)
            cmap_args = dict(cmap='jet', norm=norm)

        cm._A = []
        if self.proj_type == 'miller':
            rotation, orientation, pad, label_pad = 0, 'horizontal', 0.01, None
        else:
            #rotation, orientation, pad, label_pad = 270, 'vertical', 0.015, -0.1
            rotation, orientation, pad, label_pad = 0, 'horizontal', 0.01, None
        #   https://stackoverflow.com/questions/15908371/matplotlib-colorbars-and-its-text-labels
        if cb_type == 'zero_center' and ('d' in label or 'FAC' in label):
            cb = plt.colorbar(cm, ax=self.ax, orientation=orientation, shrink=0.35, fraction=0.05, pad=pad)
            cb.ax.set_yticks([int(v) for v in bounds])
            cb.ax.set_yticklabels([int(v) for v in bounds])
        else:
            cb = plt.colorbar(cm, ax=self.ax, orientation=orientation, shrink=0.35, fraction=0.05, pad=pad)
        cb.set_label(label, labelpad=label_pad, rotation=rotation)
        return cmap_args

    def draw_mag_coord_lines(self, date, geomag_pole=False):
        mag_lat_lines, mag_lon_lines, annotate_points = get_mag_coordinate_system_lines(date, geomag_pole)

        for lat_lines in mag_lat_lines:
            self.ax.plot(lat_lines[:, 1], lat_lines[:, 0], 'r:', zorder=1, lw=1, alpha=0.75,
                         transform=ccrs.Geodetic())
        for lon_lines in mag_lon_lines:
            self.ax.plot(lon_lines[:, 1], lon_lines[:, 0], 'r:', zorder=1, lw=1, alpha=0.75,
                         transform=ccrs.Geodetic())
        for i, xy in enumerate(annotate_points):
            if self.extend is not None:
                if self.extend[0] < xy[1] < self.extend[1] and self.extend[2] < xy[0] < self.extend[
                    3]:
                    self.ax.text(xy[1], xy[0], '%i,%i' % (xy[2][0], xy[2][1]), transform=self.transform, fontsize=6,
                                 color='r')
            else:
                self.ax.text(xy[1], xy[0], '%i,%i' % (xy[2][0], xy[2][1]), transform=self.transform, fontsize=6,
                             color='r')

    def draw_vector(self, sw_pos, B):
        X, Y = sw_pos[:, 1], sw_pos[:, 0]
        #B = sw_value[:, :2]
        u = 100
        q = self.ax.quiver(X, Y, B[:, 0], B[:, 1], transform=ccrs.PlateCarree(),
                           width=0.0006, color='m', zorder=8, alpha=0.75)
        self.ax.quiverkey(q, X=0.35, Y=-0.02, U=u,  labelpos='E', label='SWARM-model vector length {dN, dE} = %s' % u, transform=ccrs.PlateCarree())

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

    def fill_dark_side(self, datetime):
        lat, lng = sun_pos(datetime)
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

    def draw_sun_terminator(self, swarm_date, swarm_time):
        date_array = [decode_str_dt_param(d+'T'+t) for d, t in zip(swarm_date, swarm_time)]
        date = date_array[int((len(date_array)-1)/2)]
        x_term, y_term, rot_pol_term = self.fill_dark_side(date)
        self.ax.fill(x_term, y_term, transform=rot_pol_term, zorder=1, color='k', alpha=0.25)

    def draw_point_with_annotate(self, pos, annotate):
        #   for intermag observ for example
        #for point in points:
        x, y = float(pos[0]), float(pos[1])
        self.ax.scatter(x, y, c='r', marker='*', s=150, transform=self.transform, zorder=10)
        self.ax.text(x, y, annotate, transform=self.transform, fontsize=10,zorder=10)

    def draw_polygon(self, poly):
        x, y = poly.exterior.xy
        self.ax.plot(x, y, transform=ccrs.PlateCarree())

    def set_axis_label(self, label, zoom_axis=False):
        print('axis label:   ', label)
        if len(label)>100:
            fs = 12
        else:
            fs = 17
        if zoom_axis:
            #self.ax.set_title(label, fontsize=fs, y=-0.1)
            #self.ax.set_title(label, fontsize=fs, pad=1)
            self.fig.suptitle(label, fontsize=fs)
        else:
            self.ax.set_title(label, fontsize=fs)

    def fig_to_PIL_image(self, name='1'):

        buf = io.BytesIO()
        plt.savefig(buf, format='jpeg', dpi=400)
        buf.seek(0)
        im = Image.open(buf)
        im.save(STATIC_OS_PATH +'/media/tmp/%s.jpg' % name)
        buf.close()
        return im

