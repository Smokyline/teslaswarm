import os
import sys
import matplotlib.pyplot as plt
from tools.dt_foo import decode_str_dt_param
from tools.sql_table import get_swarm_set
import numpy as np
from tools.data_foo import *
from tools.dt_foo import get_timestamp
from control.proj_im_creator import *
from engine.stacker import stack_images, single_image
# Find code directory relative to our directory
sys.path.append(os.getcwd())
from teslaswarm.settings import STATIC_OS_PATH


class teslaswarmControl():
    def __init__(self, from_date, to_date):
        # swarm

        self.from_date = from_date    # <class 'datetime.datetime'> 2017-09-03 00:00:59
        self.to_date = to_date

        self.ionomodel_param = None

        self.auroral_date = None
        self.draw_auroral_n = None  # dt
        self.draw_auroral_s = None  # dt

        self.proj_extend_loc = None     # [-lon, +lon, -lat, +lat]
        self.swarm_poly_loc = None      # [[-90, 50], [-90, 90], [25, 90], [25,50], [-90, 50]]
        self.intermag_observ_code = None    # ['SFS']
        self.supermag_observ_code = None    # ['SFS']

        self.swarm_value_delta = 300
        self.deg_radius = 5
        self.cut_swarm_value_bool = False
        self.cut_obs_swarm_value_bool = False
        self.draw_vector_diff = False
        self.measure_mu = False
        self.mag_grid_coord = False
        self.draw_shape = False
        self.txt_out = False

    def get_id(self):
        self.id = get_timestamp()
        return self.id


    def save_single_image(self, im):
        out_image = single_image(im)
        out_image.save(STATIC_OS_PATH + '/media/images/request/%s.jpg'%self.id)

    def save_fourths_image(self, im1, im2, im3, im4):
        out_image = stack_images(2000, im1, im2, im3, im4)
        out_image.save(STATIC_OS_PATH + '/media/images/request/%s.jpg'%self.id)

    def set_swarm_value_delta(self, value):
        self.swarm_value_delta = value

    def set_ionomodel(self, param=None):
        # ionomodel?type=fac&hem=n&bz=1&f107=100&by=1&doy=1&kp=1&ut=1&img_w=800&img_h=600&out=plot
        if param is not None:
            param_dict = {
                'type': 'fac',
                'hem': 'n',
                'bz': '6',
                'f107': '100',
                'by': '6',
                'doy': '1',
                'kp': '1',
                'ut': '1'
            }
            for i, p in enumerate(param):
                if i == 0:
                    if p == 'FAC':
                        param_dict['type'] = 'fac'
                    elif p == 'Potential':
                        param_dict['type'] = 'pot'
                    elif p == 'sigH':
                        param_dict['type'] = 'sigh'
                if i == 1:
                    if p == 'Northern':
                        param_dict['hem'] = 'n'
                    elif p == 'Southern':
                        param_dict['hem'] = 's'
                if i == 2:
                    param_dict['bz'] = p
                if i == 3:
                    param_dict['by'] = p
                if i == 4:
                    param_dict['doy'] = p
                if i == 5:
                    param_dict['kp'] = p
                if i == 6:
                    param_dict['ut'] = p
                if i == 7:
                    param_dict['f107'] = p
        else:
            param_dict = None
        self.ionomodel_param = param_dict

    def set_auroral_oval(self, north=False, south=False, date=None):
        if north:
            self.draw_auroral_n = date
        if south:
            self.draw_auroral_s = date

    def set_proj_extend_loc(self, loc=None):
        self.proj_extend_loc = loc

    def set_swarm_cutpoly_loc(self, loc=None):
        self.swarm_poly_loc = loc

    def set_obs_codes(self, codes=None, deg_radius=5, cut_obs_swarm_value=False):
        self.intermag_observ_code = codes
        self.deg_radius = deg_radius
        self.cut_obs_swarm_value_bool = cut_obs_swarm_value

    def set_supermag_obs_codes(self, codes=None):
        self.supermag_observ_code = codes
        #TODO im creator

    def set_swarm_igrf_vectorDiff(self, b=False):
        self.draw_vector_diff = b

    def set_measure_mu(self, b=False):
        self.measure_mu = b

    def set_mag_grid_coord(self, b=False):
        self.mag_grid_coord = b

    def set_shape_file(self, file):
        # TODO import file
        self.draw_shape = file

    def set_txt_out(self, b=False):
        self.txt_out = b

    def get_swarm_set(self, sw_liter='A', fac2_mod=False):
        #   получение датасета данных SWARM_X
        return get_swarm_set(sw_liter, self.from_date, self.to_date, self.swarm_value_delta, fac2_mod)

    def get_swarmAC_diff(self, swarm_set_A, swarm_set_C, sw_channel):
        #   получение датасета данных разности значений поля между SWARM_A и SWARM_C
        sw_coord, sw_respond = calc_ACvector(sw_a_cd=swarm_set_A[1], sw_c_cd=swarm_set_C[1], sw_a_values=swarm_set_A[4],
                                             sw_c_values=swarm_set_C[4], channel=sw_channel)
        swarm_set_AC = ["|AC|", sw_coord, swarm_set_A[2], swarm_set_A[3], sw_respond]
        return swarm_set_AC

    def get_strParam2arrayParam(self, param):
        paramlist = param.split(":")
        floatlist = []
        for p in paramlist:
            if p == 'true':
                floatlist.append(True)
            elif p =='false':
                floatlist.append(False)
            else:
                try:
                    floatlist.append(float(p))
                except:
                    floatlist.append(str(p))
        return np.array(floatlist)


    def get_projection_image(self, swarm_set, swarm_channel, proj_type, annotate_sw_value_bool):
        """
        (swarm_info, proj_type,
        draw_ionomodel_n=False, draw_ionomodel_s=False, draw_auroral_s=None,
        draw_auroral_n=None, draw_shape=None, draw_vector_diff=False,
        intermag_observ_code=None, measure_mu=False,
        cut_swarm_value_bool=False, swarm_poly_loc=None, proj_extend_loc=None, annotate_sw_value_bool=False, convert_coord=None, cut_deg_radius=5):

        """
        swarm_info = [swarm_set, self.from_date, self.to_date, swarm_channel]

        (status, out) = get_proj_image(swarm_info=swarm_info, proj_type=proj_type,
                             ionomodel_param=self.ionomodel_param, draw_auroral_s=self.draw_auroral_s,
                             draw_auroral_n=self.draw_auroral_n, draw_shape=self.draw_shape, draw_vector_diff=self.draw_vector_diff,
                             intermag_observ_code=self.intermag_observ_code, measure_mu=self.measure_mu, mag_grid_coord=self.mag_grid_coord,
                             cut_swarm_value_bool=self.cut_swarm_value_bool, cut_obs_swarm_value_bool=self.cut_obs_swarm_value_bool,
                             swarm_poly_loc=self.swarm_poly_loc, proj_extend_loc=self.proj_extend_loc,
                             annotate_sw_value_bool=annotate_sw_value_bool, cut_deg_radius=self.deg_radius, txt_out=self.txt_out)
        if status == 1 and self.txt_out:
            message = out
        if status == 1 and self.txt_out == False:
            self.save_single_image(out)
            message = self.id
        if status == 0:
            message = out
        return (status, message)   # Image.open(buf)

    def get_plotted_image(self, swarm_sets, labels, auroral, sw_channel, delta, station):
        """
        swarm_sets, labels, include, channel
        swarm_sets=[fac_set_A, fac_set_B, C]
        swarm_sets=[fac_set_A, fac_set_B]
        """
        (status, out) = get_plot_im(swarm_sets, labels, auroral, sw_channel, delta, self.measure_mu, station, self.txt_out)

        if status == 1 and self.txt_out:
            message = out
        if status == 1 and self.txt_out == False:
            self.save_single_image(out)
            message = self.id
        return (status, message)


