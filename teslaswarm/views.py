from django.core.mail import send_mail
from django.shortcuts import render
from teslaswarm import settings

import os
from tools.dt_foo import decode_str_dt_param
from engine.stacker import stack_images, single_image
from teslaswarm.settings import STATIC_OS_PATH, STATICFILES_DIRS
from control.request_control import teslaswarmControl

def get_image_page(request):
    # main backend foo

    def save_single_image(im, id):
        out_image = single_image(im)
        out_image.save(STATIC_OS_PATH + '/media/images/request/%s.jpg'%id)

    def save_fourths_image(im1, im2, im3, im4, id):
        out_image = stack_images(2000, im1, im2, im3, im4)
        out_image.save(STATIC_OS_PATH + '/media/images/request/%s.jpg'%id)



    original_umask = os.umask(0)

    param_name = [
        'proj_type',    # str
        'plot_liter',    # str []
        'sw_liter',     # str
        'sw_delta',     # int
        'from_date',    # str   2019-5-13T22:59:59
        'to_date',      # str   2019-5-13T22:59:59
        'sw_channel',   # str
        # 'supermag_obs',   # str
        'annotate_sw_value',    # bool
        'ionomodel_date',  # str
        'ionomodel_n',  # bool
        'ionomodel_s',  # bool
        'auroral_date', # str
        'auroral_n',    # bool
        'auroral_s',    # bool
        'proj_extend_loc',  # str []
        'swarm_poly_loc',   # str [[]]
        'cut_obs_swarm_value',     # bool
        'obs_code',     # str []
        'deg_radius',   # float
        'igrf_vector_diff',  # bool
        'measure_mu',   # bool
        #'mag_grid',   # bool
        #'shape',        # file
        'near_auroral_points',  # bool
        ]
    param_dict = {}
    for p in param_name:
        param_dict[p] = request.GET[p]
    print(param_dict)
    """
    {'proj_type': 'miller', 'plot_liter': 'true:false:false:false', 'sw_liter': 'A', 'sw_delta': '300', 
    'from_date': '2017-09-07T00:00:00', 'to_date': '2017-09-09T00:00:00', 'sw_channel': '1', 'annotate_sw_value': 'true', 
    'ionomodel_date': 'undefinedTundefined:00', 'ionomodel_n': 'false', 'ionomodel_s': 'false', 
    'auroral_date': '2017-09-08T00:00:00', 'auroral_n': 'true', 'auroral_s': 'false', 
    'proj_extend_loc': 'undefined:undefined:undefined:undefined',  'swarm_poly_loc': 'undefined:undefined:undefined:undefined', 
    'cut_obs_swarm_value': 'false', 'obs_code': 'undefined', 'deg_radius': '5', 'igrf_vector_diff': 'true', 
    'measure_mu': 'false', 'near_auroral_points': 'false'}
    """


    dt_from, dt_to = decode_str_dt_param(param_dict['from_date']), decode_str_dt_param(param_dict['to_date'])
    sm = teslaswarmControl(from_date=dt_from, to_date=dt_to)
    id = sm.get_id()  # name as timestamp in ms

    fac2_mod = False
    if param_dict['sw_channel'] == 'n':
        sw_channel = 0
    elif param_dict['sw_channel'] == 'e':
        sw_channel = 1
    elif param_dict['sw_channel'] == 'c':
        sw_channel = 2
    elif param_dict['sw_channel'] == 'fac2':
        sw_channel = None
        fac2_mod = True
    else:
        sw_channel = 0

    if not ('undefined' in str(param_dict['sw_delta'])):
        delta = int(param_dict['sw_delta'])
        if delta >= 2:
            sm.set_swarm_value_delta(delta)
    else:
        delta = 1


    if param_dict['annotate_sw_value'] == 'true':
        annotate_sw_value_bool = True
    else:
        annotate_sw_value_bool = False

    if param_dict['proj_type'] != 'plot':
        if param_dict['sw_liter'] == '|AC|':
            swarm_set_A = sm.get_swarm_set(sw_liter='A', fac2_mod=fac2_mod)
            swarm_set_C = sm.get_swarm_set(sw_liter='C', fac2_mod=fac2_mod)
            SWARM_SET = sm.get_swarmAC_diff(swarm_set_A=swarm_set_A, swarm_set_C=swarm_set_C, sw_channel=sw_channel)
        else:
            SWARM_SET = sm.get_swarm_set(sw_liter=param_dict['sw_liter'], fac2_mod=fac2_mod)

    #   модель ионосферы
    if param_dict['ionomodel_n'] == 'true' or param_dict['ionomodel_s'] == 'true':
        #TODO ionomodel date
        ionomodel_date = decode_str_dt_param(param_dict['ionomodel_date'])
        if param_dict['ionomodel_n'] == 'true':
            ionomodel_n = True
        else:
            ionomodel_n = False
        if  param_dict['ionomodel_s'] == 'true':
            ionomodel_s = True
        else:
            ionomodel_s = False
        sm.set_ionomodel(request=None)

    #   модель аврорального овала
    if param_dict['auroral_n'] == 'true' or param_dict['auroral_s'] == 'true':
        auroral_date = decode_str_dt_param(param_dict['auroral_date'])
        if param_dict['auroral_n'] == 'true':
            auroral_n = True
        else:
            auroral_n = False
        if param_dict['auroral_s'] == 'true':
            auroral_s = True
        else:
            auroral_s = False
        sm.set_auroral_oval(north=auroral_n, south=auroral_s, date=auroral_date)


    if not ('undefined' in str(param_dict['proj_extend_loc'])):
        sm.set_proj_extend_loc(loc=sm.get_strParam2arrayParam(param_dict['proj_extend_loc']))
    if not ('undefined' in str(param_dict['swarm_poly_loc'])):

        sm.set_swarm_cutpoly_loc(loc=sm.get_strParam2arrayParam(param_dict['swarm_poly_loc']))

    # INTERMAGNET observatories
    if not ('undefined' in str(param_dict['obs_code'])):
        sm.set_obs_codes(codes=param_dict['obs_code'], deg_radius=float(param_dict['deg_radius']), cut_obs_swarm_value=True)

    # vector difference between swarm and IGRF-13 model
    if param_dict['igrf_vector_diff'] == 'true':
        igrf_vector_diff = True
    else:
        igrf_vector_diff = False
    sm.set_swarm_igrf_vectorDiff(b=igrf_vector_diff)

    # convert swarm value to DMA anomaly measure mu
    if param_dict['measure_mu'] == 'true':
        measure_mu = True
    else:
        measure_mu = False
        sm.set_measure_mu(b=measure_mu)

    # геомагнитная сетка координат mag_coordinate_system_lines
    if param_dict['mag_grid'] == 'true':
        mag_grid = True
    else:
        mag_grid = False
        sm.mag_grid(b=mag_grid)

    #TODO add file import
    #if param_dict['shape'] == 'true':
    sm.set_shape_file(file=None)


    ###############################################
    if param_dict['proj_type'] == 'plot':
        # append auroral value to plot
        # [x for x in bool_set1 if x==True]
        swarm_sets = []
        labels = []
        include_data = []
        bool_set = sm.get_strParam2arrayParam(param_dict['plot_liter'])
        station = None
        for i, selected in enumerate(bool_set):
            if selected:
                if i == 0:   # if A, B, C
                    swarm_set = sm.get_swarm_set(sw_liter='A', fac2_mod=fac2_mod)
                elif i == 1:   # if A, B, C
                    swarm_set = sm.get_swarm_set(sw_liter='B', fac2_mod=fac2_mod)
                elif i == 2:   # if A, B, C
                    swarm_set = sm.get_swarm_set(sw_liter='C', fac2_mod=fac2_mod)
                elif i == 3:   # if |AC|
                    swarm_set_A = sm.get_swarm_set(sw_liter='A', fac2_mod=fac2_mod)
                    swarm_set_C = sm.get_swarm_set(sw_liter='C', fac2_mod=fac2_mod)
                    swarm_set = sm.get_swarmAC_diff(swarm_set_A=swarm_set_A, swarm_set_C=swarm_set_C, sw_channel=sw_channel)

                if param_dict['near_auroral_points'] == 'true':
                    include_data.append(sm.get_near_auroral_points_to_swarm(swarm_set=swarm_set))
                else:
                    include_data.append(None)
                if fac2_mod:
                    labels.append('swarm-%s %s' % (swarm_set[0], 'fac2'))
                else:
                    lbs = ['n', 'e', 'c']
                    labels.append('swarm-%s %s' % (swarm_set[0], lbs[sw_channel]))
                swarm_sets.append(swarm_set)

                if param_dict['supermag_obs'] != '-':
                    station = param_dict['supermag_obs']


        print(bool_set, '=> plot sets:', len(swarm_sets), labels, 'include:')
        im = sm.get_plotted_image(swarm_sets=swarm_sets, labels=labels, include=include_data, sw_channel=sw_channel, delta=delta, station=station)
    else:
        im = sm.get_projection_image(swarm_set=SWARM_SET, swarm_channel=sw_channel, proj_type=param_dict['proj_type'],
                                     annotate_sw_value_bool=annotate_sw_value_bool)

    save_single_image(im, id)
    os.umask(original_umask)

    return render(request, 'show_image.html', {'STATIC_URL': settings.STATIC_URL, 'IMAGE_NAME': id})

def test(request):
    conda_env = os.environ['CONDA_DEFAULT_ENV']
    return render(request, 'test_page.html', {'env_name': conda_env})

def get_homepage(request):
    return render(request, 'index.html')