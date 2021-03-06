import copy
import datetime

from django.core.mail import send_mail
from django.shortcuts import render
from teslaswarm import settings
from django.http import HttpResponse
import os
from tools.dt_foo import decode_str_dt_param
from teslaswarm.settings import STATIC_OS_PATH, ALLOWED_HOSTS
from control.request_control import teslaswarmControl
from engine.swarm_animation import SwarmAnimRender
from tools.ionomodel_foo import get_ionomodel_surf_file, draw_map

TeslaSwarm = teslaswarmControl()
Swarm_anim = SwarmAnimRender

def get_teslaswarm(request):
    # main backend foo

    original_umask = os.umask(0)

    param_name = [
        'proj_type',    # str
        'plot_liter',    # str []
        'sw_liter',     # str
        'sw_delta',     # int
        'from_date',    # str   2019-5-13T22:59:59
        'to_date',      # str   2019-5-13T22:59:59
        'sw_channel',   # str
        'supermag_obs',   # str
        'annotate_sw_value',    # bool
        'dataType',  # str
        'hemisphere',  # str
        'iono_bz',  # str
        'iono_by',  # str
        'iono_kp',  # str
        'iono_f107',  # str
        'iono_doy',  # str
        'iono_ut',  # str
        'auroral_date', # str
        'auroral_n',    # bool
        'auroral_s',    # bool
        #'auroral_type',    # str
        'proj_extend_loc',  # str []
        'swarm_poly_loc',   # str [[]]
        'cut_obs_swarm_value',     # bool
        'obs_code',     # str []
        'supermag_obs',     # str []
        'deg_radius',   # float
        'igrf_vector_diff',  # bool
        'igrf_vector_diffChaos',  # bool
        'measure_mu',   # bool
        'mag_grid',   # bool
        #'geomag_proj',        # bool
        'near_auroral_points',  # bool
        'txt_out',  # bool
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
    sm = copy.deepcopy(TeslaSwarm)
    id = sm.get_id()  # name as timestamp in ms

    fac_mod = False
    if param_dict['sw_channel'] == 'n':
        sw_channel = 0
    elif param_dict['sw_channel'] == 'e':
        sw_channel = 1
    elif param_dict['sw_channel'] == 'c':
        sw_channel = 2
    elif param_dict['sw_channel'] == 'f':
        sw_channel = 3
    elif param_dict['sw_channel'] == 'fac':
        sw_channel = None
        fac_mod = True
    else:
        sw_channel = 0

    if not ('undefined' in str(param_dict['sw_delta'])):
        delta = int(param_dict['sw_delta'])
        sm.set_swarm_value_delta(delta)


    if param_dict['annotate_sw_value'] == 'true':
        annotate_sw_time_bool = True
    else:
        annotate_sw_time_bool = False

    if param_dict['igrf_vector_diffChaos'] == 'true':
        sm.set_swarm_chaos_vectorDiff(b=True)
    elif param_dict['igrf_vector_diff'] == 'true':
        sm.set_swarm_igrf_vectorDiff(b=True)

    if param_dict['sw_liter'] == 'AC' and not fac_mod:
        return render(request, 'error.html', {'ERROR_MESSAGE': 'only FAC is available for SWARM-AC '})

    if param_dict['proj_type'] != 'plot':
        try:
            # vector difference between swarm and IGRF-13 model
            if param_dict['sw_liter'] == '|A-C|':  # |A-C|
                swarm_set_A = sm.get_swarm_set(sw_liter='A', from_date=dt_from, to_date=dt_to, fac_mod=fac_mod)
                swarm_set_C = sm.get_swarm_set(sw_liter='C', from_date=dt_from, to_date=dt_to, fac_mod=fac_mod)
                SWARM_SET = sm.get_swarmAC_diff(swarm_set_A=swarm_set_A, swarm_set_C=swarm_set_C, sw_channel=sw_channel)
            elif param_dict['sw_liter'] == 'AC':  # J2
                SWARM_SET = sm.get_swarm_set(sw_liter='_', from_date=dt_from, to_date=dt_to, fac_mod=True)
            else:  # J1
                SWARM_SET = sm.get_swarm_set(sw_liter=param_dict['sw_liter'], from_date=dt_from, to_date=dt_to,
                                             fac_mod=fac_mod)
        except:
            return render(request, 'error.html', {'ERROR_MESSAGE': 'NO SWARM VALUE'})

    #   ???????????? ??????????????????
    if not ('undefined' in str(param_dict['iono_bz'])):
        iono_param = [param_dict['dataType'], param_dict['hemisphere'], param_dict['iono_bz'], param_dict['iono_by'],
                      param_dict['iono_doy'], param_dict['iono_kp'], param_dict['iono_ut'], param_dict['iono_f107']]

    else:
        iono_param = None
    sm.set_ionomodel(param=iono_param)

    #   ???????????? ???????????????????????? ??????????
    if param_dict['auroral_n'] == 'true' or param_dict['auroral_s'] == 'true':
        auroral_date = decode_str_dt_param(param_dict['auroral_date'])
        # TODO
        #auroral_type = param_dict['auroral_type']
        auroral_type = 'mono'
        if param_dict['auroral_n'] == 'true':
            auroral_n = True
        else:
            auroral_n = False
        if param_dict['auroral_s'] == 'true':
            auroral_s = True
        else:
            auroral_s = False
        sm.set_auroral_oval(north=auroral_n, south=auroral_s, date=auroral_date, type=auroral_type)


    if not ('undefined' in str(param_dict['proj_extend_loc'])):
        sm.set_proj_extend_loc(loc=sm.get_strParam2arrayParam(param_dict['proj_extend_loc']))
        #sm.set_swarm_cutpoly_loc(loc=sm.get_strParam2arrayParam(param_dict['swarm_poly_loc']))


    # INTERMAGNET observatories
    if not ('-' in str(param_dict['obs_code'])):
        sm.set_intermag_obs_codes(codes=param_dict['obs_code'], deg_radius=float(param_dict['deg_radius']))

    # superMAG stations
    if not ('-' in str(param_dict['supermag_obs'])):
        sm.set_supermag_obs_codes(codes=param_dict['supermag_obs'], deg_radius=float(param_dict['deg_radius']))

    # convert swarm value to DMA anomaly measure mu
    if param_dict['measure_mu'] == 'true':
        measure_mu = True
    else:
        measure_mu = False
    sm.set_measure_mu(b=measure_mu)

    # ???????????????????????? ?????????? ?????????????????? mag_coordinate_system_lines
    if param_dict['mag_grid'] == 'true':
        mag_grid = True
    else:
        mag_grid = False
    sm.set_mag_grid_coord(b=mag_grid)

    #TODO
    """if param_dict['geomag_proj'] == 'true':
        geomag_proj = True
    else:
        geomag_proj = False
    sm.set_geomag_proj(b=geomag_proj)"""

    if param_dict['txt_out'] == 'true':
        sm.set_txt_out(b=True)

    print(vars(sm), 'request control')
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
                try:
                    if i == 0:   # if A, B, C
                        swarm_set = sm.get_swarm_set(sw_liter='A', from_date=dt_from, to_date=dt_to, fac_mod=fac_mod)
                    elif i == 1:   # if A, B, C
                        swarm_set = sm.get_swarm_set(sw_liter='B', from_date=dt_from, to_date=dt_to, fac_mod=fac_mod)
                    elif i == 2:   # if A, B, C
                        swarm_set = sm.get_swarm_set(sw_liter='C', from_date=dt_from, to_date=dt_to, fac_mod=fac_mod)
                    elif i == 3:   # if A&C
                        # A&C
                        swarm_set = sm.get_swarm_set(sw_liter='_', from_date=dt_from, to_date=dt_to, fac_mod=True)
                    elif i == 4:    # if |A-C|
                        swarm_set_A = sm.get_swarm_set(sw_liter='A', from_date=dt_from, to_date=dt_to, fac_mod=fac_mod)
                        swarm_set_C = sm.get_swarm_set(sw_liter='C', from_date=dt_from, to_date=dt_to, fac_mod=fac_mod)
                        swarm_set = sm.get_swarmAC_diff(swarm_set_A=swarm_set_A, swarm_set_C=swarm_set_C,
                                                        sw_channel=sw_channel)
                except:
                    return render(request, 'error.html', {'ERROR_MESSAGE': 'NO SWARM VALUE'})
                if fac_mod:
                    labels.append('swarm-%s %s' % (swarm_set[0], 'fac'))
                else:
                    lbs = ['n', 'e', 'c', 'f']
                    labels.append('swarm-%s %s' % (swarm_set[0], lbs[sw_channel]))
                swarm_sets.append(swarm_set)

                if param_dict['supermag_obs'] != '-':
                    station = param_dict['supermag_obs']

        if param_dict['near_auroral_points'] == 'true':
            # include_data.append(sm.get_near_auroral_points_to_swarm(swarm_set=swarm_set))
            get_auroral = True
        else:
            get_auroral = False
            # include_data.append(None)
        print(bool_set, '=> plot sets:', len(swarm_sets), labels, 'include:', include_data, 'superMAG station:', station)
        (status, message) = sm.get_plotted_image(swarm_sets=swarm_sets, labels=labels, auroral=get_auroral,
                                                 sw_channel=sw_channel, delta=delta, station=station)
    else:
        (status, message) = sm.get_projection_image(swarm_set=SWARM_SET, from_date=dt_from, to_date=dt_to, swarm_channel=sw_channel,
                                                    proj_type=param_dict['proj_type'], annotate_sw_time_bool=annotate_sw_time_bool)
    os.umask(original_umask)

    if status in [0, 2]:
        return render(request, 'error.html', {'ERROR_MESSAGE': message})

    id = message
    if sm.txt_out:
        data_file = open(STATIC_OS_PATH + '/media/txt/%s.txt' % id, 'r', encoding="utf-8")
        data = data_file.readlines()
        """context = {'context': data}
        return render(request, "show_txt.html", context)   #{{ file_content }}."""
        #response_content = '\n'.join(lines)
        return HttpResponse(data, content_type="text/plain")
    else:
        return render(request, 'show_image.html', {'STATIC_URL': settings.STATIC_URL, 'IMAGE_NAME': id})


def test(request):
    conda_env = os.environ['CONDA_DEFAULT_ENV']
    return render(request, 'test_page.html', {'env_name': conda_env})

def get_teslaswarm_homepage(request):
    return render(request, 'index.html')


def get_swarmAnim(request):

    # yield render_to_response('send_mail.html',)

    original_umask = os.umask(0)
    param_name = [
        'swchar',
        'mod',
        'timeFrom',
        'timeTo',
        'channel',
        'delta',
        'mailTo']
    param_dict = {}
    for p in param_name:
        param_dict[p] = request.GET[p]
    dt_from, dt_to = decode_str_dt_param(
        param_dict['timeFrom']), decode_str_dt_param(
        param_dict['timeTo'])
    print(param_dict)
    r = Swarm_anim(dt_from=dt_from, dt_to=dt_to, LT_interval=[None], delta=int(param_dict['delta']),
               mod=param_dict['mod'], swarm_char=param_dict['swchar'],
               channel=int(param_dict['channel']), frame_step=1)
    print('data downloaded from server...')
    led = r.render_swarm_video()
    link = 'http://%s/static/media/videos/%s.webm' % (ALLOWED_HOSTS[0], led)
    #link = 'http://%s/static/videos/%s.webm' % ('http://127.0.0.1:8000', led)

    channels = ['North', 'East', 'Down', 'Vector modulus']
    if param_dict['mod'] == 'vector':
        swarm_sat = 'Swarm A & Swarm C'
        mod = 'Difference vector magnitude'
    elif param_dict['mod'] == 'mera':
        swarm_sat = 'Swarm %s' % str(param_dict['swchar']).upper()
        mod = 'Measure of anomalies'
    elif param_dict['mod'] == 'fac2':
        swarm_sat = 'Swarm %s' % str(param_dict['swchar']).upper()
        mod = 'Time series of field-aligned currents'
    channel = channels[int(param_dict['channel'])]
    delta = param_dict['delta']

    mail_text = "%s %s\nComponent: %s\nData period from: %s\nto: %s\nAveraging interval, seconds: %s\n\nlink to video: %s\nPlease, save the video If you want to use it in the future." % (
        mod, swarm_sat, channel, param_dict['timeFrom'], param_dict['timeTo'], delta, link
    )

    send_mail('Swarm map', mail_text, 'swarm@gcras.ru', [
              param_dict['mailTo']], fail_silently=False)
    os.umask(original_umask)

    return render(request, 'send_mail.html')

def swarm_anim_form(request):
    return render(request, 'dataserv-swarmtracks.html')

def get_ionimodel(request):
    param_name = ['type', 'hem', 'bz', 'f107', 'by', 'doy', 'kp', 'ut', 'out', 'img_w', 'img_h']
    param_dict = {}
    for p in param_name:
        param_dict[p] = request.GET[p]

    iono_xyz = get_ionomodel_surf_file(param_dict, datetime.date(2020, 1, 1), polar=True)
    if param_dict['out'] == 'plot':
        # show map
        #img_name = draw_geo_map(data_type)
        img_name = draw_map(iono_xyz, param_dict['hem'], param_dict['hem'])
        path = 'http://%s/static/media/images/request' % ALLOWED_HOSTS[0]
        return render(request, 'ionoimage_view.html', {'MEDIA_URL': path, 'IMAGE_NAME': img_name,
                                                  'IMG_W': param_dict['img_w'], 'IMG_H': param_dict['img_h']})
    else:
        # ascii
        str_array = ''
        for x, y, z, in iono_xyz:
            s = '%.2f %.2f %.2f\n' % (x, y, z)
            str_array += s
        return render(request,'print_ascii.html', {'XYZ': str_array})

def ionomodel_form(request):
    return render(request, 'dataserv-ionomodel-ru.html')

def show_mail_page(request):

    return render(request, 'send_mail.html')

def show_video(request, led):
    led = '/static/media/videos/%s.webm' % led
    return render(request, 'show_video.html', {'led': led})