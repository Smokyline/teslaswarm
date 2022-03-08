from tools.sql_table import get_sql_response, get_swarm_set
from tools.dt_foo import decode_str_dt_param
import os
import sys
import pymysql
import datetime
from tools.data_foo import *
from tools.dt_foo import ut_dt_to_unix, unix_to_ut_dt
import matplotlib.pyplot as plt
from teslaswarm.settings import STATIC_OS_PATH


def plot_data(sw_set):
    sw_liter, sw_position, sw_date, sw_time, sw_values = sw_set
    fig = plt.figure(figsize=(16, 9))
    ax = fig.add_subplot(111)
    #ax.set_aspect('auto')
    str_ticks = [decode_str_dt_param(d + 'T' + t).strftime('%c') for d, t in zip(sw_date, sw_time)]
    ax_ticks = np.array([[x, label] for x, label in zip(np.arange(len(sw_time)), str_ticks)])
    ax_ticks = ax_ticks[np.arange(0, len(ax_ticks), int(len(ax_ticks) / 6))]

    ax.vlines(x=np.array(ax_ticks[:, 0]).astype(int), ymin=np.min(sw_values), ymax=np.max(sw_values), color='purple',linestyle='--', alpha=.5)
    ax.plot(np.arange(len(sw_values)), sw_values, c='r', label='swb')



    ax.set_ylabel('fac')
    ax.set_xticks(np.array(ax_ticks[:, 0]).astype(int))
    ax.set_xticklabels(ax_ticks[:, 1])
    ax.grid()
    ax.legend()
    plt.grid(True)
    plt.savefig(STATIC_OS_PATH +'/media/images/swb_fac.png')


sw_liter = 'b'
from_date = '2017-9-7T23:39:59'
to_date = '2017-9-8T00:20:59'
from_date, to_date = decode_str_dt_param(from_date), decode_str_dt_param(to_date)
delta = 1
fac2_mod = True
chaos_mod = False




#swarm_set = [sw_liter, np.array(swarm_position), swarm_date, swarm_time, swarm_values]
swarm_set = get_swarm_set(sw_liter, from_date, to_date, delta=delta, fac2_mod=True)
plot_data(swarm_set)
