import os
import sys
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

intermag_observ_code = ['RBY', 'CDR', 'IQA', 'SALU', 'KJPK']
cut_obs_swarm_value_bool = False
sw_liter = 'B'
sw_channel = 0
from_date = '2017-09-18T02:25:00'
to_date = '2017-09-18T02:45:00'
auroral_date = '2017-9-18T02:35:00'
delta = 3

auroral_date = decode_str_dt_param(auroral_date)
from_date, to_date = decode_str_dt_param(from_date), decode_str_dt_param(to_date)

# swarm_liter, swarm_position, swarm_date, swarm_time, swarm_values
swarm_set_A = get_swarm_set('A', from_date, to_date, delta=delta, fac2_mod=False)
#swarm_set_B = get_swarm_set('B', from_date, to_date, delta=delta, fac2_mod=True)
#swarm_set_C = get_swarm_set('C', from_date, to_date, delta=delta, fac2_mod=False)

# swarm |AC|
#sw_coord, sw_respond = calc_ACvector(sw_a_cd=swarm_set_A[1], sw_c_cd=swarm_set_C[1], sw_a_values=swarm_set_A[4], sw_c_values=swarm_set_C[4], channel=sw_channel)
#swarm_set_AC = ["|AC|", sw_coord, swarm_set_A[2], swarm_set_A[3], sw_respond]



###   idk
"""swarm_liter, swarm_pos, swarm_date, swarm_time, swarm_values = swarm_set
X, Y = swarm_pos[:, 0], swarm_pos[:, 1]
Z = swarm_values.ravel()"""

fac_set_A = get_swarm_set('A', from_date, to_date, delta=delta, fac2_mod=True)
#fac_set_B = get_swarm_set('B', from_date, to_date, delta=delta, fac2_mod=True)
#fac_set_C = get_swarm_set('C', from_date, to_date, delta=delta, fac2_mod=True)
#   if use in plot -> sw_channel = None

#auroral_to_swarm_A = get_nearest_auroral_point_to_swarm(fac_set_A)
#auroral_to_swarm_B = get_nearest_auroral_point_to_swarm(fac_set_B)
#auroral_to_swarm_C = get_nearest_auroral_point_to_swarm(fac_set_C)
#nearest_swarm_value_to_auroral = None
#print(nearest_swarm_value_to_auroral)

"""
    poly_points = [i1, i2, i3, i4]
    i1 = [lat, lon]
    [i1]--------[i2]
      |         |
      |  p_in_p |
      |         |
    [i4]--------[i3]
"""
deg_radius = 5
cut_swarm_value_bool = True
#swarm_poly_loc = [[5.216440319940263, 36.14439512806877], [5.036148385331252, 39.089590556206], [3.583236376927949, 41.82329674087288], [0.8822804435796463, 44.01053854353759], [-2.8151035614305653, 45.346561204742905], [-7.00262518765703, 45.62511937772345], [-11.022005347427989, 44.801017981065726], [-14.269177732767334, 43.00571102110124], [-16.361639709284102, 40.503845423467006], [-17.17028028589505, 37.622771177745385], [-16.753427449254605, 34.693377424056], [-15.275963973227991, 32.01527990046319], [-12.95557003962975, 29.839575863666585], [-10.036853955099946, 28.35919842581214], [-6.780384987028307, 27.701494104710893], [-3.4558840366958283, 27.922292893516445], [-0.33499161959349766, 29.00289553518439], [2.3163077717404232, 30.851114902687968], [4.244494587473613, 33.305954911117304], [5.216440319940263, 36.14439512806877]]
#swarm_poly_loc = [[-15.864257812499998,36.756490329505176],[-46.494140625,19.72534224805787 ],[ -3.1640625,17.476432197195518],[-1.60400390625,39.13006024213511 ],[  -3.80126953125,46.13417004624326],[  -15.402832031250002,45.336701909968134],[ -15.864257812499998,36.756490329505176]]
#swarm_poly_loc = [[-15.864257812499998,76.756490329505176],[-46.494140625,69.72534224805787 ],[ -3.1640625,67.476432197195518],[-1.60400390625,79.13006024213511 ],[  -3.80126953125,86.13417004624326],[  -15.402832031250002,85.336701909968134],[ -15.864257812499998,76.756490329505176]]
#swarm_poly_loc = [[-90, 50], [-90, 90], [25, 90], [25,50], [-90, 50]]
#proj_extend_loc = None    #  [-x, +x, -y, +y]
proj_extend_loc = [-100, -60, 50, 75]    #  [-x, +x, -y, +y]
#proj_extend_loc = [-30, 30, -90, 90]    #  [-x, +x, -y, +y]
"""
    [im1]--------[im3]
      |          |
      |          |
      |          |
    [im2]--------[im4]
"""

#   swarm_info=[fac_set_A, from_date, to_date, sw_channel]
im1 = get_proj_image(swarm_info=[swarm_set_A, from_date, to_date, sw_channel],
                     #swarm_info=[fac_set_A, from_date, to_date, None],
                     #proj_type='ortho_n', draw_auroral_n=auroral_date,)
                     #proj_type='ortho_n',
                     proj_type='miller',
                     #draw_vector_diff=True,
                     #proj_extend_loc=proj_extend_loc,
                     draw_auroral_n=auroral_date,
                    #draw_ionomodel_s=True,

                        txt_out=True
                     )
"""im2 = projection_im(swarm_info=[swarm_set_A, from_date, to_date, sw_channel],
                    proj_type='ortho_s', draw_ionomodel_s=True, proj_extend_loc=[1, 360, -45, -90])

im3 = projection_im(swarm_info=[swarm_set_B, from_date, to_date, sw_channel],
                    proj_type='miller', intermag_observ_code=intermag_observ_code, cut_swarm_value_bool=cut_swarm_value_bool, swarm_poly_loc=swarm_poly_loc)

#im3 = projection_im(swarm_info=[swarm_set_B, from_date, to_date, sw_channel], proj_type='miller')
#im4 = plot_vector(swarm_sets=[fac_set_A, fac_set_B, fac_set_C], labels=['fac swA', 'fac swB', 'fac swC'], include=[auroral_to_swarm_A, auroral_to_swarm_B, None], channel=None)

#im4 = plot_vector(swarm_sets=[swarm_set_B,], labels=['dBn', 'dBe', 'dBd', 'dB'], include=[None, None, None, None], channel=None)
#im4 = plot_vector(swarm_sets=[fac_set_A,], labels=['fac sw_a'], channel=None)

"""
#im4 = get_plot_im(swarm_sets=[fac_set_A, fac_set_B, fac_set_C], labels=['fac swA', 'fac swB', 'fac swC'], include=[None, None, None], channel=None)
#im4 = get_plot_im(swarm_sets=[fac_set_A, fac_set_B], labels=['fac swA', 'fac swB', None], include=[None, None, None], channel=None)
#im4 = get_plot_im(swarm_sets=[fac_set_A, fac_set_B], labels=['fac swA', 'fac swB', None], include=[auroral_to_swarm_A, auroral_to_swarm_B, None], channel=None)
#im4 = get_plot_im(swarm_sets=[fac_set_A, fac_set_B, fac_set_C], labels=['fac swA', 'fac swB', 'fac swC'], include=[auroral_to_swarm_A, auroral_to_swarm_B, None], channel=None)
#im4 = get_plot_im(swarm_sets=[swarm_set_A], labels=['swarm-A n', None, None ], include=[None, None, None], channel=0, obs='SALU')

#im4 = get_plot_im(swarm_sets=[swarm_set_A], labels=['swarm-A n', 'swarm-A e', 'swarm-A c' ], include=[None, None, None], channel=None, delta=delta, ground_station='T47')
#im4 = get_plot_im(swarm_sets=[fac_set_A], labels=['swarm-A fac2', None, None ], include=[None, None, None], channel=None, delta=delta, ground_station='T47')

#   fourths
#out_image = stack_images(2000, im1, im2, im3, im4)
#out_image.save(STATIC_OS_PATH + '/media/images/test_v8.png')

#   single
out_image = single_image(im1)
#out_image = single_image(im4)
filename = 'single_test_v16_test'





print('saved as %s' % STATIC_OS_PATH + '/media/images/%s.jpg' % filename)
out_image.save(STATIC_OS_PATH + '/media/images/%s.jpg' % filename)


# get_shapefile()

# sword = SWORD()
# sword.draw_plot([swarm_set_A[4][:, 0],swarm_set_A[4][:, 1], swarm_set_A[4][:, 2]], swarm_time)
# im = sword.fig_to_PIL_image()
# plt.show(im)

# plt.show(projection_im(type='ortho_n'))
