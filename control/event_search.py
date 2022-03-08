import os
import sys
import glob

import pandas as pd

from engine.stacker import stack_images, single_image
import matplotlib.pyplot as plt
from tools.dt_foo import decode_str_dt_param
from tools.sql_table import get_swarm_set
import numpy as np
from tools.data_foo import *
from control.proj_im_creator import *

# Find code directory relative to our directory
sys.path.append(os.getcwd())
from teslaswarm.settings import STATIC_OS_PATH





def read_obs_events_txt():
    "чтение ВСЕХ событий"
    input_events = []

    with open('C:\\Users\\ivan\\workspace\\teslaswarm\\static\\data\\MIE_Events_Salluit__2015_and_2017.txt') as f:
        lines = f.readlines()
    for i, line in enumerate(lines):
        if i == 0:
            continue
        st_line = line.strip("\t")
        if len(st_line) < 20:
            continue
        st_line = st_line.split("	")
        title = st_line[0]
        y = '20' + st_line[0][:2]
        m = st_line[1]
        if len(m) > 2:
            m = m[:2]
        d = st_line[2]
        time = st_line[3]
        str_dt = '%s-%s-%sT%s:00' % (y, m, d, time)
        try:
            event_dt = decode_str_dt_param(str_dt)
        except Exception as e:
            print(e)
        input_events.append([title, str_dt])

    return np.array(input_events)   # title, str_dt

def read_selected_events_txt(target=[]):
    "чтение событий, выделенных main на 1 этапе"
    array_event_liter = []
    with open('C:\\Users\\ivan\\workspace\\teslaswarm\\static\\media\\images\\SALU\\SALU_events.txt', mode="r",  encoding="utf-8") as f:
        lines = f.readlines()
    for line in lines:
        event_liter = line.rsplit('U', 1)[1][:7]
        event, liter = event_liter.split('_')
        if len(target) == 0:
            array_event_liter.append([event, liter])
        else:
            if int(event) in target:
                array_event_liter.append([event, liter])
    return np.array(array_event_liter)  # event, liter


def main(input_events, selected=None, delta_minute=10):


        if selected is not None:
            selected_titles, selected_liters = selected[:, 0], selected[:, 1]  # from read_selected_events_txt()

            for i, event in enumerate(input_events):    # from read_obs_events_txt()
                event_title, event_dt = event

                idx = np.where(selected_titles == event_title)[0]
                if len(idx) > 0:
                    selected_liter = selected_liters[idx][0]
                    #selected_liter = "A"
                    event_dt = decode_str_dt_param(event_dt)
                    from_date = event_dt - pd.Timedelta('%s minutes' % delta_minute)
                    to_date = event_dt + pd.Timedelta('%s minutes' % delta_minute)
                    try:
                        """swarm_setA = get_swarm_set('A', from_date, to_date, delta=delta,
                                                                              fac2_mod=fac2)
                                                    swarm_setC = get_swarm_set('C', from_date, to_date, delta=delta,
                                                                              fac2_mod=fac2)
                                                    sw_coord, sw_respond = calc_ACvector(sw_a_cd=swarm_setA[1], sw_c_cd=swarm_setC[1], sw_a_values=swarm_setA[4], sw_c_values=swarm_setC[4], channel=sw_channel)
                                                    swarm_set = ["|AC|", sw_coord, swarm_setA[2], swarm_setA[3], sw_respond]
                                                    """
                        swarm_set = get_swarm_set(str(selected_liter), from_date, to_date, delta=delta,
                                                      fac2_mod=fac2)
                    except:
                        continue
                    im1 = get_proj_image(swarm_info=[swarm_set, from_date, to_date, sw_channel], proj_type='miller',
                                         proj_extend_loc=proj_extend_loc, intermag_observ_code=intermag_observ_code, measure_mu=False, draw_auroral_n=from_date,)

                    #auroral_to_swarm = get_nearest_auroral_point_to_swarm(swarm_set)
                    #im1 = get_plot_im(swarm_sets=[swarm_set], labels=['fac sw%s'%selected_liter, None, None ], include=[auroral_to_swarm, None, None], channel=None)
                    #im1 = get_plot_im(swarm_sets=[swarm_set], labels=['fac sw%s'%selected_liter, None, None ], include=[None, None, None], channel=None)

                    filename = 'SALU' + event_title + '_%s-anomaly' % str(selected_liter)
                    single_image(im1).save(STATIC_OS_PATH + '/media/images/%s.jpg' % filename)
                    print('saved as %s' % STATIC_OS_PATH + '/media/images/%s.jpg' % filename)
        else:
            for i, event in enumerate(input_events):    # from read_obs_events_txt()
                event_title, event_dt = event
                event_dt = decode_str_dt_param(event_dt)
                from_date = event_dt - pd.Timedelta('%s minutes' % 20)
                to_date = event_dt + pd.Timedelta('%s minutes' % 20)

                # swarm_liter, swarm_position, swarm_date, swarm_time, swarm_values
                try:
                    swarm_set_A = get_swarm_set('A', from_date, to_date, delta=delta, fac2_mod=fac2)
                    swarm_set_B = get_swarm_set('B', from_date, to_date, delta=delta, fac2_mod=fac2)
                except:
                    continue

                im1 = get_proj_image(swarm_info=[swarm_set_A, from_date, to_date, sw_channel], proj_type='miller',
                                     proj_extend_loc=proj_extend_loc, intermag_observ_code=intermag_observ_code)
                im2 = get_proj_image(swarm_info=[swarm_set_B, from_date, to_date, sw_channel], proj_type='miller',
                                     proj_extend_loc=proj_extend_loc, iintermag_observ_code=intermag_observ_code)

                filename = 'SALU' + event_title + '_A'
                single_image(im1).save(STATIC_OS_PATH + '/media/images/%s.jpg' % filename)
                filename = 'SALU' + event_title + '_B'
                single_image(im2).save(STATIC_OS_PATH + '/media/images/%s.jpg' % filename)
                print('saved as %s' % STATIC_OS_PATH + '/media/images/%s.jpg' % filename)




"""jpeg_in_folder = glob.glob(glob.escape('C:\\Users\\ivan\\workspace\\teslaswarm\\static\\media\\images\\SALU\\') + "*.jpg")
for jpeg in jpeg_in_folder:
    elements = jpeg.split('\\')
    print(elements[9].split('.')[0])"""
intermag_observ_code = ['RBY', 'CDR', 'IQA', 'SALU', 'KJPK']
delta = 4

fac2=True
sw_channel = None   # None equal fac2

proj_extend_loc = [-100, -60, 50, 75]    #  lon_min, lon_max, lat_min, lat_max,

selected_event_liter = read_selected_events_txt(target=[17261])
all_input_events = read_obs_events_txt()

# main(input_events=all_input_events)    # 1st stage
main(input_events=all_input_events, selected=selected_event_liter, delta_minute=10)  # 2nd stage
