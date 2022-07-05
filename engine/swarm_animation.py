import matplotlib
matplotlib.use('Agg')
import time
import tqdm
import copy

from teslaswarm.settings import BASE_DIR, FFMPEG_PATH
from tools.data_foo import *
from tools.sun_position import fill_dark_side
from tools.sql_table import *
import cartopy.feature as cfeature
import cartopy.crs as ccrs

import matplotlib as mpl
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.collections import LineCollection
from matplotlib import animation
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt




class SwarmAnimRender():
    def __init__(self, dt_from, dt_to, LT_interval, delta,
                 mod, swarm_char, channel, frame_step):
        self.led = int(round(time.time() * 1000))

        print('led:%s' % self.led)
        LT_interval = [datetime.time(0, 0, 0), datetime.time(23, 59, 59)]
        self.dt_from, self.dt_to, self.LT_interval, self.delta, self.mod, self.swarm_char, self.channel, self.frame_step \
            = dt_from, dt_to, LT_interval, delta, mod, swarm_char, channel, frame_step

        self.fig, self.ax1, self.ax2, self.ax3 = self.get_axis()



    def get_axis(self):
        """создание fig и axis"""
        fig = plt.figure(figsize=(17, 9))  # x, y
        gs = gridspec.GridSpec(6, 7)  # y, x
        ax1 = plt.subplot(gs[:3, 5:], projection=ccrs.Orthographic(central_latitude=90.,
                                                                   central_longitude=0))
        ax2 = plt.subplot(gs[:, :5], projection=ccrs.Miller())
        ax3 = plt.subplot(gs[3:, 5:], projection=ccrs.Orthographic(central_latitude=-90.,
                                                                   central_longitude=0))
        for ax in [ax1, ax2, ax3]:
            ax.coastlines(resolution='110m', color='k', zorder=7)
            ax.add_feature(cfeature.LAND, facecolor='0.75', zorder=0)
            ax.gridlines()
        return fig, ax1, ax2, ax3

    def im_swarm_satellite(self):
        swarm_char, swarm_pos, swarm_date, swarm_time, Z = self.swarm_set
        sw_LAT, sw_LON, sw_R = swarm_pos[:, 0], swarm_pos[:, 1], swarm_pos[:, 2]

        if self.mod == 'vector':

            label_set = ['dBn', 'dBe', 'dBd', 'dB']
            label = label_set[int(self.channel)]
            # cmap = mpl.colors.LinearSegmentedColormap.from_list('mycolors', ['blue', 'red'])

            cmap = plt.get_cmap('cool')
            sm = plt.cm.ScalarMappable(
                cmap=cmap, norm=plt.Normalize(
                    vmin=Z.min(), vmax=Z.max()))
        elif self.mod == 'mera':
            label_set = ['MAn', 'MAe', 'MAd', 'MAf']
            label = label_set[int(self.channel)]
            cmap = mpl.colors.ListedColormap(['#794044', '#8a2be2', 'b', '#0099cc',
                                              '#00ffff', '#008000', '#00ff00', '#ccff00',
                                              '#ffff00', '#ff7f50', '#ff00ff', 'r'])
            bounds = [-1.0, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0,
                      0.1, 0.2, 0.3, 0.4, 0.5, 1]
            norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        elif self.mod == 'fac2':
            label = 'FAC\nμА/m^2'
            # cmap = mpl.colors.LinearSegmentedColormap.from_list('mycolors', ['blue', 'red'])
            cmap = plt.get_cmap('seismic')

            sm = plt.cm.ScalarMappable(
                cmap=cmap, norm=plt.Normalize(
                    vmin=-np.abs(Z).max(), vmax=np.abs(Z).max()))
        sm._A = []
        cb = plt.colorbar(sm, ax=self.ax2, shrink=0.7)
        cb.ax.text(0.05, 1.025, label, rotation=0)
        try:
            cb.set_ticks(bounds)
        except BaseException:
            pass


        self.ax1.plot(sw_LON, sw_LAT,'k--', zorder=1,lw=1, alpha=0.4,transform=ccrs.Geodetic())
        self.ax2.plot(sw_LON, sw_LAT,'k--',zorder=1,lw=1,alpha=0.4,transform=ccrs.Geodetic())
        self.ax3.plot(sw_LON, sw_LAT,'k--',zorder=1,lw=1,alpha=0.4,transform=ccrs.Geodetic())

        ims_swarm = []
        swarm_lt = []
        for i in np.arange(len(swarm_pos)):
            lat, lon, r = sw_LAT[i], sw_LON[i], sw_R[i]
            if self.mod == 'vector':
                lt = get_local_time(decode_str_dt_param(swarm_date[i]+'T'+swarm_time[i]), lon)
                satel_frame = [
                    (self.ax2.scatter(lon, lat, c='k', marker='^', zorder=4,
                                      transform=self.ax2.projection)),
                    (self.ax2.scatter(lon, lat, c='k', marker='^', zorder=4,
                                      transform=self.ax2.projection)),
                    (self.ax2.text(lon - 10,
                                   lat - 0,
                                   'A',
                                   color='k',
                                   transform=self.ax2.projection)),
                    (self.ax2.text(lon + 4,
                                   lat + 0,
                                   'C',
                                   color='k',
                                   transform=self.ax2.projection)),
                    (self.ax2.annotate('LT: ' + str(lt),
                                       xy=(0.8, 1), xycoords='axes fraction', fontsize=16,
                                       xytext=(0, +20), textcoords='offset points', ha='left', va='top')),
                ]
            else:  # mera or fac2
                lt = get_local_time(decode_str_dt_param(swarm_date[i]+'T'+swarm_time[i]), lon)
                satel_frame = [
                    (self.ax2.scatter(lon, lat, c='k', marker='^',
                                      zorder=4, transform=self.ax2.projection)),
                    (self.ax2.text(lon + 10, lat + 10, 'Swarm-%s' % swarm_char.upper(), color='k',
                                   transform=self.ax2.projection)),
                    (self.ax2.annotate('LT: ' + str(lt),
                                       xy=(0.8, 1), xycoords='axes fraction', fontsize=16,
                                       xytext=(0, +20), textcoords='offset points', ha='left', va='top')), ]

            x_term, y_term, rot_pol_term = fill_dark_side(
                swarm_date[i], swarm_time[i])

            satel_frame.extend(
                [(self.ax2.annotate('UTC: ' + swarm_time[i], xy=(0.6, 1), xycoords='axes fraction', fontsize=16,
                                    xytext=(0, +20), textcoords='offset points', ha='right', va='top')),
                 (self.ax2.annotate(swarm_date[i], xy=(0, 1), xycoords='axes fraction', fontsize=16,
                                    xytext=(0, +20), textcoords='offset points', ha='left', va='top')),

                 (self.ax1.fill(x_term, y_term, transform=rot_pol_term,
                                zorder=1, color='k', alpha=0.65)[0]),
                 (self.ax2.fill(x_term, y_term, transform=rot_pol_term,
                                zorder=1, color='k', alpha=0.65)[0]),
                 (self.ax3.fill(x_term, y_term, transform=rot_pol_term,
                                zorder=1, color='k', alpha=0.65)[0]),
                 ])

            ims_swarm.append(satel_frame)
            swarm_lt.append(lt)
        self.swarm_set.append(swarm_lt)
        return ims_swarm

    def im_swarm_response_plot(self, ims_swarm):
        swarm_char, swarm_pos, swarm_date, swarm_time, Z, swarm_lt = self.swarm_set

        if self.mod == 'mera':
            cmap = mpl.colors.ListedColormap(['#794044', '#8a2be2', 'b', '#0099cc',
                                              '#00ffff', '#008000', '#00ff00', '#ccff00',
                                              '#ffff00', '#ff7f50', '#ff00ff', 'r'])
            levels = [-1.0, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0,
                      0.1, 0.2, 0.3, 0.4, 0.5, 1]
        elif self.mod == 'fac2':  # if swarm.mod == 'fac2' or 'vector'
            cmap = plt.get_cmap('seismic')
            # cmap = mpl.colors.LinearSegmentedColormap.from_list('mycolors', ['blue', 'red'])
            levels = MaxNLocator(nbins=10).tick_values(Z.min(), Z.max())

        elif self.mod == 'vector':
            cmap = plt.get_cmap('cool')
            levels = MaxNLocator(nbins=10).tick_values(Z.min(), Z.max())


        norm = BoundaryNorm(levels, cmap.N, clip=True)
        lat, lon, r = swarm_pos.T

        Z = Z.ravel()
        IMS = []
        del_segm_idx = []
        print('rendered')

        for i in tqdm.tqdm(range(len(ims_swarm))):
            points = np.array([lon[:i], lat[:i]]).T.reshape(-1, 1, 2)
            segments = np.concatenate([points[:-1], points[1:]], axis=1)

            time_bool = (
                self.LT_interval[0] < swarm_lt[i]) and (
                swarm_lt[i] < self.LT_interval[1])

            if time_bool is not True:
                del_segm_idx.append(i)
            segments = np.delete(segments, del_segm_idx, axis=0)

            lc1 = LineCollection(segments, cmap=cmap,norm=norm, lw=2,alpha=0.98,
                                 transform=ccrs.Geodetic())
            lc2 = LineCollection( segments, cmap=cmap,norm=norm,lw=4,alpha=0.90,
                                  transform=ccrs.Geodetic())
            lc3 = LineCollection(segments,cmap=cmap,norm=norm,lw=2,alpha=0.98,
                                 transform=ccrs.Geodetic())
            #self.ax1.projection
            lc1.set_array(Z)
            lc2.set_array(Z)
            lc3.set_array(Z)

            #IMS.append(np.append(ims_swarm[i], [self.ax1.add_collection(lc1),self.ax2.add_collection( copy.copy(lc2)),self.ax3.add_collection(copy.copy(lc3)),]))
            IMS.append(np.append(ims_swarm[i], [self.ax1.add_collection(lc1),self.ax2.add_collection( copy.copy(lc2)),self.ax3.add_collection(copy.copy(lc3)),]))
        return IMS

    def render_swarm_video(self):
        """вызывается из views.py"""

        if self.mod == 'vector':
            sws_a = get_swarm_set('A', self.dt_from, self.dt_to, self.delta, fac2_mod=False)
            sws_c = get_swarm_set('C', self.dt_from, self.dt_to, self.delta, fac2_mod=False)

            # swarm_char, swarm_position, swarm_date, swarm_time, swarm_resp
            sw_coord, sw_respond = calc_ACvector(sw_a_cd=sws_a[1], sw_c_cd=sws_c[1],
                                                 sw_a_values=sws_a[4],
                                                 sw_c_values=sws_c[4], channel=self.channel)
            self.swarm_set = ["|AC|", sw_coord, sws_a[2], sws_a[3], sw_respond]

        if self.mod == 'mera':

            self.swarm_set = get_swarm_set(str(self.swarm_char).upper(), self.dt_from, self.dt_to,self.delta,
                                           fac2_mod=False)

            # swarm_char, swarm_position, swarm_date, swarm_time, swarm_resp
            sw_respond = get_measure_mu(self.swarm_set[4][:, self.channel])
            self.swarm_set[4] = sw_respond

        if self.mod == 'fac2':
            # swarm_char, swarm_position, swarm_date, swarm_time, swarm_resp
            self.swarm_set = get_swarm_set(str(self.swarm_char).upper(), self.dt_from, self.dt_to,self.delta,
                                           fac2_mod=True)


        ims_swarm_satel = self.im_swarm_satellite()
        IMS = self.im_swarm_response_plot(ims_swarm_satel)

        im_ani = animation.ArtistAnimation(self.fig, IMS, blit=True)

        plt.rcParams['animation.ffmpeg_path'] = FFMPEG_PATH
        # plt.rcParams['animation.bitrate'] = 2000
        Writer = animation.writers['ffmpeg']
        writer = Writer(fps=25, metadata=dict(artist='gcras'),
                        extra_args=['-vcodec', 'libvpx-vp9',
                                    # extra_args=['-vcodec', 'libx264',
                                    # '-preset','ultrafast',
                                    '-b:v', '2000k', '-vsync', '0',
                                    ])
        path = STATIC_OS_PATH + '/media/videos/'
        print('figure animation complete!')
        im_ani.save(path + '%s.webm' % self.led, writer=writer)
        print('video saved')
        return self.led