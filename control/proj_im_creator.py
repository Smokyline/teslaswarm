from engine.sword import SWORD
from tools.data_to_txt import Data2Text
from tools.data_foo import *
from tools.ionomodel_foo import get_ionomodel_surf_file
from chaos7_model.chaos_model import CHAOS7
from tools.sql_table import get_swarm_set

def get_proj_image(swarm_info, proj_type,
                   ionomodel_param=None, draw_auroral_s=None,
                   draw_auroral_n=None, draw_shape=None, draw_IGRFvector_diff=False, draw_CHAOSvector_diff=False,
                   observ_code_value=None, measure_mu=False, mag_grid_coord=False,
                   cut_swarm_value_bool=False, proj_extend_loc=None,
                   annotate_sw_time_bool=False, cut_deg_radius=5, txt_out=False, delta=1):
    swarm_liter, swarm_pos, swarm_date, swarm_time, swarm_values_necf = swarm_info[0]   # swarm_set
    from_date, to_date, swarm_channel = swarm_info[1], swarm_info[2], swarm_info[3]
    d2txt = Data2Text()
    STATUS = 1

    legend_label = swarm_liter
    if len(np.array(swarm_values_necf).shape) > 1 and swarm_channel is not None:
        swarm_values = swarm_values_necf[:, swarm_channel]
        channel_str = ['N', 'E', 'C', 'F'][swarm_channel]
    else:
        swarm_values = swarm_values_necf  # fac
        channel_str = 'FAC'
    legend_label += ' ' + channel_str

    #   приближение до OBS области
    #   extend_loc = [-x, +x, -y, +y]
    #if cut_obs_swarm_value_bool == True and proj_type=='miller':
    #if cut_obs_swarm_value_bool == True and proj_type != 'plot':
    if observ_code_value is not None:
        obs_org, obs_code = observ_code_value.split('_')
        obs_location = [0, 0]
        if obs_org == 'intermag':
            obs_location = get_INTERMAGNET_observ_loc(obs_code)[:2]
        if obs_org == 'supermag':
            obs_location = get_superMAG_observ_loc(obs_code)[:2]
        if proj_extend_loc is None:
            proj_extend_loc = get_swarm_poly_loc(obs_location, deg_radius=cut_deg_radius)

    #   установка границ проекции на основе координат области значений swarm
    if proj_extend_loc is not None:
        #swarm_poly_loc = np.array(swarm_poly_loc)
        lat_min, lat_max = np.min(proj_extend_loc[2:]), np.max(proj_extend_loc[2:])
        lon_min, lon_max = np.min(proj_extend_loc[:2]), np.max(proj_extend_loc[:2])
        if lat_max >= 90:
            lat_max = 89.99
        if lat_min <= -90:
            lat_min = -89.99
        if lon_max >= 180:
            lon_max = 179.99
        if lon_min <= -180:
            lon_min = -179.99
        proj_extend_loc = [lon_min, lon_max, lat_min, lat_max]

    sword = SWORD(proj_type=proj_type, extend=proj_extend_loc)

    dt_from = decode_str_dt_param(swarm_date[0]+'T'+swarm_time[0])
    dt_to = decode_str_dt_param(swarm_date[-1]+'T'+swarm_time[-1])
    if (dt_to - dt_from) < datetime.timedelta(minutes=30):
        sword.draw_sun_terminator(swarm_date, swarm_time)

    #   подпись к проеции
    ax_label = fr'swarm_{swarm_liter} {from_date} -- {to_date} '
    if sword.extend is not None or observ_code_value is not None:
        ax_label += 'loc:lat{%0.2f:%0.2f},lon{%0.2f:%0.2f}' % (sword.extend[2], sword.extend[3],
                                                   sword.extend[0], sword.extend[1],)
        cut_swarm_value_bool = True


    ######################################################################################################
    #   инициация рендера, указание проекции, и если требуется extend_loc - приближение конкретной области

    #   отрисовка на проекции иономодели
    if ionomodel_param is not None:
        iono_lonlatZ = get_ionomodel_surf_file(ionomodel_param, dt_from)
        if ionomodel_param['hem'] == 'n':
            sword.draw_ionomodel(iono_lonlatZ, type=['northern', ionomodel_param['type']])
            #d2txt.ionomodel_shp = 'n'
        elif ionomodel_param['hem'] == 's':
            sword.draw_ionomodel(iono_lonlatZ, type=['southern', ionomodel_param['type']])
            #d2txt.ionomodel_shp = 's'
        key = 'ionomodel'+str(ionomodel_param['hem']).upper()+'_'+ionomodel_param['type']
        d2txt.DATA[key] = {}
        d2txt.DATA[key]['lon'] = iono_lonlatZ[:, 0]
        d2txt.DATA[key]['lat'] = iono_lonlatZ[:, 1]
        if ionomodel_param['type'] == 'sigh':
            value_key = 'Sm'
        elif ionomodel_param['type'] == 'pot':
            value_key = 'kV'
        else:
            value_key = 'uA/m^2'
        d2txt.DATA[key][value_key] = iono_lonlatZ[:, 2]
        ax_label += '\n iono' + str(ionomodel_param)


    #   отрисовка на проекции аврорального овала
    if draw_auroral_s is not None:
        x, y, z, yxz = get_auroral_flux(draw_auroral_s['date'], hemishpere='S', atype=draw_auroral_s['type'])  # f(dt)
        sword.draw_avroral_oval([x,y,z], hemisphere='south', atype=draw_auroral_s['type'])
        #sword.draw_mag_coord_lines(swarm_date[0])
        ax_label += ' auroralS'
        d2txt.DATA['auroral_s'] = {}
        d2txt.DATA['auroral_s']['lat'] = yxz[:, 0]
        d2txt.DATA['auroral_s']['lon'] = yxz[:, 1]
        d2txt.DATA['auroral_s']['%s_ergs/cm^2'%draw_auroral_s['type']] = yxz[:, 2]
    if draw_auroral_n is not None:
        x, y, z, yxz = get_auroral_flux(draw_auroral_n['date'], hemishpere='N', atype=draw_auroral_n['type'])  # f(dt)
        sword.draw_avroral_oval([x,y,z], hemisphere='north', atype=draw_auroral_n['type'])
        #sword.draw_mag_coord_lines(swarm_date[0])
        ax_label += ' auroralN'
        d2txt.DATA['auroral_n'] = {}
        d2txt.DATA['auroral_n']['lat'] = yxz[:, 0]
        d2txt.DATA['auroral_n']['lon'] = yxz[:, 1]
        d2txt.DATA['auroral_n']['%s_ergs/cm^2'%draw_auroral_n['type']] = yxz[:, 2]

    #   отрисовка на проекции shape файла
    if draw_shape is not None:
        sword.draw_shapefile(get_shapefile())
        ax_label += ' shapefile '

    if mag_grid_coord:
        sword.draw_mag_coord_lines(swarm_date[0], geomag_pole=True)    # True to_coord='MAG' /  False to_coord='GSM'

    d2txt.annotate = ax_label

    # выбор точек измерений SWARM в указанном полигоне
    if cut_swarm_value_bool:
        if sword.extend is not None:
            print('cut swarm_pos manual by proj_extend_loc')
            poly_loc = [[sword.extend[1],sword.extend[2],], [sword.extend[1],sword.extend[3],],
                        [sword.extend[0],sword.extend[3],], [sword.extend[0],sword.extend[2],]]
            p_in_p, poly = get_points_in_poly(swarm_pos[:, :2], poly_loc, proj_type)
        if observ_code_value is not None:

            print('cut swarm_pos near obs')
            #p_in_p, poly = get_swarm_value_near_obs(swarm_pos, intermag_observ_code, proj_type)
            p_in_p = get_position_near_point(swarm_pos[:, :2], obs_location, degr_radius=cut_deg_radius)

        swarm_pos_in_poly = swarm_pos[p_in_p]
        swarm_values_in_poly = swarm_values[p_in_p]
        swarm_values_full_in_poly = swarm_values_necf[p_in_p]
        swarm_date_in_poly, swarm_time_in_poly = swarm_date[p_in_p], swarm_time[p_in_p]
        print('len cut pos', len(swarm_pos_in_poly))
        vline_min_path = True
        if len(swarm_pos_in_poly) == 0:
            STATUS = 0
            try:
                out = 'NO SWARM DATA IN TIME INTERVAL FROM %s TO %s IN LOCATION lat:%s... %s lon:%s... %s\n\n' \
                      'hint: expand specified location, date range or decrease delta' % (
                          swarm_date[0] + 'T' + swarm_time[0], swarm_date[-1] + 'T' + swarm_time[-1],
                          proj_extend_loc[0],
                          proj_extend_loc[1], proj_extend_loc[2], proj_extend_loc[3],)
            except:
                out = 'NO SWARM DATA IN TIME INTERVAL FROM %s TO %s \n\n' \
                      'hint: expand specified location, date range or decrease delta' % (
                          swarm_date[0] + 'T' + swarm_time[0], swarm_date[-1] + 'T' + swarm_time[-1])
            return (STATUS, out)

    else:

        swarm_pos_in_poly = swarm_pos
        swarm_values_in_poly = swarm_values
        swarm_values_full_in_poly = swarm_values_necf
        swarm_date_in_poly, swarm_time_in_poly = swarm_date, swarm_time
        vline_min_path = False
    d2txt.DATA['SWARM_dt'] = [date + 'T' + time for date, time in zip(swarm_date_in_poly, swarm_time_in_poly)]
    d2txt.DATA['SWARM_GEOpos'] = {}
    d2txt.DATA['SWARM_GEOpos']['lat'] = swarm_pos_in_poly[:, 0]
    d2txt.DATA['SWARM_GEOpos']['lon'] = swarm_pos_in_poly[:, 1]
    d2txt.DATA['SWARM_GEOpos']['r'] = swarm_pos_in_poly[:, 2]
    # отрисовка временных линий на пролете спутника
    if vline_min_path:
        swarm_dt = []
        for date, time in zip(swarm_date, swarm_time):
            swarm_dt.append(date + 'T' + time)
        vline_dt = np.array(swarm_dt)
    else:
        vline_dt = None
    sword.draw_swarm_path(swarm_pos, points_time=vline_dt, annotate=annotate_sw_time_bool)

    # отрисовка вектора (X, Y, n, e) или (X, Y, |n-x|, |e-y|)
    if draw_IGRFvector_diff or draw_CHAOSvector_diff:
        if swarm_liter not in ['A', 'B', 'C'] or swarm_channel is None:
            STATUS = 2  # ERROR: |A-C|, A&C or FAC
            out = 'If you need a difference between the satellite data and the magnetic field model,' \
                  ' please select only the satellite SWARM-A, SWARM-B or SWARM-C, and also do not choose FAC'
            return (STATUS, out)
        else:
            vector_components = swarm_values_full_in_poly  # n,e,c,f component
            if draw_IGRFvector_diff:
                swarm_dt_in_poly = [decode_str_dt_param(d + 'T' + t) for d, t in
                                    zip(swarm_date_in_poly, swarm_time_in_poly)]
                # dN, dE, dC, dH
                vector_subtraction = swarm_egrf_vector_subtraction(swarm_pos_in_poly, vector_components,
                                                                   swarm_dt_in_poly)
                model_name = 'IGRF'
            elif draw_CHAOSvector_diff:
                chaosm = CHAOS7(
                    swarm_set=[swarm_liter, swarm_pos_in_poly, swarm_date_in_poly, swarm_time_in_poly,
                               vector_components])
                # dN, dE, dC, dH
                vector_subtraction = chaosm.get_swarm_chaos_vector_subtraction()
                model_name = 'CHAOS7'
            # print(vector_subtraction, 'vector substr')
            swarm_values_in_poly = vector_subtraction[:, 3]  # dH
            vector_components = vector_subtraction[:, (0, 1)]  # dx, dy

            # legend_label = 'SWARM-%s_%s' % (model_name, ['dN', 'dE', 'dC', 'dF'][swarm_channel])
            legend_label = 'SWARM_%s-%s $\\mathrm{d} H$' % (legend_label.rstrip(' ')[0], model_name)
            key = 'SWARM_%s-%s_dH' % (legend_label.rstrip(' ')[0], model_name)
            d2txt.DATA[key] = {}
            d2txt.DATA[key]['dH_nT'] = swarm_values_in_poly
            d2txt.DATA[key]['dN_nT'] = vector_components[:, 0]
            d2txt.DATA[key]['dE_nT'] = vector_components[:, 1]
            ax_label += ' %s ' % legend_label
    if measure_mu:
        swarm_values_in_poly = get_measure_mu(swarm_values_in_poly)
        ax_label += ' measure_mu '
        legend_label += ' mu'
        d2txt.DATA['DMA_anomaly'] = swarm_values_in_poly
    else:
        if channel_str != 'FAC':
            value_key = 'SWARM_%s_' % swarm_liter + channel_str + '_nT'
        else:
            value_key = 'SWARM_%s_' % swarm_liter + channel_str + '_uA/m^2'
        d2txt.DATA[value_key] = np.ravel(swarm_values_in_poly)

    """    
        if convert_coord is not None:
            # конвертация географических координат в dest_sys
            swarm_pos_in_poly = convert_coord_system(swarm_pos_in_poly, dest_sys=convert_coord)
            ax_label += 'coord_sys:%s ' % convert_coord
        vector_components = swarm_egrf_vector_subtraction(swarm_pos_in_poly, vector_components, swarm_date)
        sword.draw_vector(swarm_pos_in_poly, B=vector_components)"""
    #sword.draw_mag_coord_lines(swarm_date_in_poly[0], geomag_pole=False)
    #sword.draw_swarm_path(swarm_pos_in_poly[:, :2])

    # отрисовка точек измерений swarm в в указанном полигоне
    # если полигона нет - отрисовка всех (swarm_pos_in_poly = swarm_pos)

    sword.draw_swarm_scatter(swarm_pos_in_poly, swarm_values_in_poly, custom_label=legend_label,
                                 annotate=False)  # отрисовка значение точек на орбите
    if draw_CHAOSvector_diff or draw_IGRFvector_diff:
        sword.draw_vector(swarm_pos_in_poly, B=vector_components)
    if observ_code_value is not None:
        sword.draw_point_with_annotate(obs_location, annotate=obs_code)
    # конвертация figure matplotlib в PIL image для stacker.py
    if cut_swarm_value_bool == True or proj_extend_loc is not None:
        sword.set_axis_label(ax_label, zoom_axis=True)
    else:
        sword.set_axis_label(ax_label)


    if txt_out:
        d2txt.save_columns_to_txt()
        out = d2txt.id
    else:
        out = sword.fig_to_PIL_image()

    return (STATUS, out)


def get_plot_im(swarm_sets, labels, auroral, channel, delta, measure_mu, ground_station=None, txt_out=False):
    #swarm_liter, swarm_pos, swarm_date, swarm_time, swarm_values = swarm_info[0]
    #from_date, to_date = swarm_info[1], swarm_info[2]
    # инициация рендера
    sword = SWORD(proj_type='plot')

    #   отрисовка графика
    draw_list = []
    auroral_list = []
    d2txt = Data2Text()

    obs_value = []
    if ground_station is not None:
        position_list, date_list, time_list = swarm_sets[0][1], swarm_sets[0][2], swarm_sets[0][3]
        # N, E, C, F
        obs_value = get_sMAGstation_value_by_time(date_list, time_list, delta=delta, station=ground_station)
        d2txt.DATA[ground_station] = {}
        for col, vect_comp in enumerate(['N', 'E', 'Z', 'F']):
            d2txt.DATA[ground_station][vect_comp + '_nT'] = obs_value[:, col]

    for i, swarm_set in enumerate(swarm_sets):
        position_list, date_list, time_list = swarm_set[1], swarm_set[2], swarm_set[3]
        d2txt.annotate = 'from %s to %s' % (date_list[0]+'T'+time_list[0], date_list[-1]+'T'+time_list[-1])
        #   find swarm value shape size
        if len(np.array(swarm_set[4]).shape) == 1:
            shape_size = 1
        else:
            shape_size = np.array(swarm_set[4]).shape[1]



        d2txt.DATA['SWARM_dt'] = [date+'T'+time for date, time in zip(date_list, time_list)]
        #   append value[channel] or value
        if channel is not None:
            label = labels[i]
            try:
                value = np.array(swarm_set[4])[:, channel]
            except:
                value = swarm_set[4]
            if measure_mu:
                value = get_measure_mu(value)
                value_key = 'SWARM_%s_' % label.rsplit('-')[1][0] + '_DMA_anomaly'
            else:
                value_key = 'SWARM_%s_' % label.rsplit('-')[1][0] + ['N', 'E', 'C', 'F'][channel] + '_nT'

            draw_list.append([label, date_list, time_list, value, position_list])

            d2txt.DATA['SWARM_GEOpos'] = {}
            d2txt.DATA['SWARM_GEOpos']['lat'] = position_list[:, 0]
            d2txt.DATA['SWARM_GEOpos']['lon'] = position_list[:, 1]
            d2txt.DATA['SWARM_GEOpos']['r'] = position_list[:, 2]
            d2txt.DATA[value_key] = np.ravel(value)


        if channel is None:
            if shape_size == 1:
                label = labels[i]
                value = swarm_set[4]
                if measure_mu:
                    value = get_measure_mu(value)
                    value_key = 'SWARM_%s_' % label.rsplit('-')[1][0] + '_FAC_DMA_anomaly'
                else:
                    value_key = 'SWARM_%s_' % label.rsplit('-')[1][0] + 'FAC' + '_nT'
                draw_list.append([label, date_list, time_list, value, position_list])
                d2txt.DATA['SWARM_GEOpos'] = {}
                d2txt.DATA['SWARM_GEOpos']['lat'] = position_list[:, 0]
                d2txt.DATA['SWARM_GEOpos']['lon'] = position_list[:, 1]
                d2txt.DATA['SWARM_GEOpos']['r'] = position_list[:, 2]
                d2txt.DATA[value_key] = np.ravel(value)

            else:
                for ch in range(4):
                    label = labels[ch]
                    value = np.array(swarm_set[4])[:, ch]
                    if measure_mu:
                        value = get_measure_mu(value)
                        value_key = 'SWARM_%s_' % label.rsplit('-')[1][0] + '_DMA_anomaly'
                    else:
                        value_key = 'SWARM_%s_' % label.rsplit('-')[1][0] + ['N', 'E', 'C', 'F'][ch] + '_nT'
                    draw_list.append([label, date_list, time_list, value, position_list])
                    if len(obs_value) > 0:
                        draw_list.append([label, date_list, time_list, value, position_list])
                    d2txt.DATA['SWARM_GEOpos'] = {}
                    d2txt.DATA['SWARM_GEOpos']['lat'] = position_list[:, 0]
                    d2txt.DATA['SWARM_GEOpos']['lon'] = position_list[:, 1]
                    d2txt.DATA['SWARM_GEOpos']['r'] = position_list[:, 2]
                    d2txt.DATA[value_key] = np.ravel(value)

        if auroral:
            auroral_near_swarm_diff = get_nearest_auroral_point_to_swarm(swarm_set, atype='diff')
            auroral_near_swarm_mono = get_nearest_auroral_point_to_swarm(swarm_set, atype='mono')
            auroral_list.append([auroral_near_swarm_diff, auroral_near_swarm_mono])
            d2txt.DATA['%s_auroral_near_diff_ergs/cm^2' % label] = auroral_near_swarm_diff
            d2txt.DATA['%s_auroral_near_mono_ergs/cm^2' % label] = auroral_near_swarm_mono



    status = 1
    if status == 1:
        if txt_out:
            d2txt.save_columns_to_txt()
            out = d2txt.id
        else:
            sword.draw_plot(draw_list, auroral_list, delta, draw_station={ground_station:obs_value})
            out = sword.fig_to_PIL_image()
    return (status, out)

