from engine.sword import SWORD
from tools.data_to_txt import Data2Text
from tools.data_foo import *
from tools.ionomodel_foo import get_ionomodel_surf_file
from chaos7_model.chaos_model import CHAOS7

def get_proj_image(swarm_info, proj_type,
                   ionomodel_param=None, draw_auroral_s=None,
                   draw_auroral_n=None, draw_shape=None, draw_vector_diff=False,
                   intermag_observ_code=None, measure_mu=False, mag_grid_coord=False,
                   cut_swarm_value_bool=False, cut_obs_swarm_value_bool=False, swarm_poly_loc=None, proj_extend_loc=None,
                   annotate_sw_value_bool=False, cut_deg_radius=5, txt_out=False):
    swarm_liter, swarm_pos, swarm_date, swarm_time, swarm_values_nec = swarm_info[0]   # swarm_set
    from_date, to_date, swarm_channel = swarm_info[1], swarm_info[2], swarm_info[3]
    d2txt = Data2Text(SWARM_liter=swarm_liter)
    STATUS = 1

    if len(np.array(swarm_values_nec).shape) > 1 and swarm_channel is not None:
        FAC2 = True
    else:
        FAC2 = False

    #   приближение до OBS области
    #   extend_loc = [-x, +x, -y, +y]
    #if cut_obs_swarm_value_bool == True and proj_type=='miller':
    if cut_obs_swarm_value_bool == True and proj_type != 'plot':
        obs_location = get_INTERMAGNET_observ_loc(intermag_observ_code)[:2]
        proj_extend_loc = get_swarm_poly_loc(obs_location, deg_radius=cut_deg_radius)
        cut_swarm_value_bool = True

    #   установка границ проекции на основе координат области значений swarm
    if swarm_poly_loc is not None:
        swarm_poly_loc = np.array(swarm_poly_loc)
        lat_min, lat_max = np.min(swarm_poly_loc[:, 1]) - 5, np.max(swarm_poly_loc[:, 1]) + 5
        lon_min, lon_max = np.min(swarm_poly_loc[:, 0]) - 5, np.max(swarm_poly_loc[:, 0]) + 5
        if lat_max >= 90:
            lat_max = 89.9
        if lat_min <= -90:
            lat_min = -89.9
        if lon_max >= 180:
            lon_max = 179.9
        if lon_min <= -180:
            lon_min = -179.9
        proj_extend_loc = [lon_min, lon_max, lat_min, lat_max, ]
        cut_swarm_value_bool = True

    #   подпись к проеции
    ax_label = fr'swarm_{swarm_liter} {from_date} -- {to_date} '
    if proj_extend_loc is not None:
        ax_label += 'loc:lat{%0.2f:%0.2f},lon{%0.2f:%0.2f}' % (proj_extend_loc[2], proj_extend_loc[3],
                                                   proj_extend_loc[0], proj_extend_loc[1],)
        cut_swarm_value_bool = True

    ######################################################################################################
    #   инициация рендера, указание проекции, и если требуется extend_loc - приближение конкретной области
    sword = SWORD(proj_type=proj_type, extend=proj_extend_loc)

    #   отрисовка на проекции иономодели
    if ionomodel_param is not None:
        fac_surf = get_ionomodel_surf_file(ionomodel_param)
        if ionomodel_param['hem'] == 'n':
            sword.draw_ionomodel(fac_surf, type=['northern', ionomodel_param['type']])
            d2txt.ionomodel_shp = 'n'
        elif ionomodel_param['hem'] == 's':
            sword.draw_ionomodel(fac_surf, type=['southern', ionomodel_param['type']])
            d2txt.ionomodel_shp = 's'
        d2txt.append(np.array([fac_surf[:,1], fac_surf[:,0], fac_surf[:,2]]).T,
                     name='ionomodel'+str(ionomodel_param['hem']).upper()+'_'+ionomodel_param['type'])
        ax_label += ' iono' + str(ionomodel_param)


    #   отрисовка на проекции аврорального овала
    if draw_auroral_s is not None:
        x, y, z, yxz = get_auroral_flux(draw_auroral_s, hemishpere='S')  # f(dt)
        sword.draw_avroral_oval([x,y,z], hemisphere='south')
        #sword.draw_mag_coord_lines(swarm_date[0])
        ax_label += ' auroralS'
        d2txt.append(yxz, name='auroralS')
        d2txt.auroral_shp = 's'
    if draw_auroral_n is not None:
        x, y, z, yxz = get_auroral_flux(draw_auroral_n, hemishpere='N')  # f(dt)
        sword.draw_avroral_oval([x,y,z], hemisphere='north')
        #sword.draw_mag_coord_lines(swarm_date[0])
        ax_label += ' auroralN'
        d2txt.append(yxz, name='auroralN')
        d2txt.auroral_shp = 'n'

    #   отрисовка на проекции shape файла
    if draw_shape is not None:
        sword.draw_shapefile(get_shapefile())
        ax_label += ' shapefile '


    if mag_grid_coord:
        sword.draw_mag_coord_lines(swarm_date[0], geomag_pole=True)    # True to_coord='MAG' /  False to_coord='GSM'


    #   отрисовка на проеции местоположения обсерваторий INTERMAGNET
    if intermag_observ_code is not None:
        obs_location = get_INTERMAGNET_observ_loc(intermag_observ_code)  # lon, lat, code
        sword.draw_point_with_annotate(obs_location)
        lon, lat = obs_location[:2]

    #   выбор канала n, e, c
    if FAC2:
        swarm_values = swarm_values_nec[:, swarm_channel]
        legend_label = ['dBn', 'dBe', 'dBd'][swarm_channel]
        d2txt.SWARM_channel = ['N', 'E', 'C'][swarm_channel]
    else:
        swarm_values = swarm_values_nec  # fac2
        legend_label = 'fac2'
        d2txt.SWARM_channel = 'FAC2'

    d2txt.annotate = ax_label

    # выбор точек измерений SWARM в указанном полигоне
    # lat, lon
    """
    poly_points = [i1, i2, i3, i4] # lat... lat, lon... lon
    i1 = [lat, lon]
    [i1]--------[i2]
      |         |
      |  p_in_p |
      |         |
    [i4]--------[i3]
    """

    if cut_swarm_value_bool:
        if swarm_poly_loc is not None:
            print('cut swarm_pos manual')
            p_in_p, poly = get_points_in_poly(swarm_pos[:, :2], swarm_poly_loc, proj_type)
            sword.draw_polygon(poly)
        elif proj_extend_loc is not None:
            print('cut swarm_pos manual by proj_extend_loc')
            poly_loc = [[proj_extend_loc[1],proj_extend_loc[2],], [proj_extend_loc[1],proj_extend_loc[3],],
                        [proj_extend_loc[0],proj_extend_loc[3],], [proj_extend_loc[0],proj_extend_loc[2],]]
            p_in_p, poly = get_points_in_poly(swarm_pos[:, :2], poly_loc, proj_type)
        elif intermag_observ_code is not None:
            print('cut swarm_pos near obs')
            #p_in_p, poly = get_swarm_value_near_obs(swarm_pos, intermag_observ_code, proj_type)
            p_in_p = get_position_near_point(swarm_pos[:, :2], intermag_observ_code, degr_radius=cut_deg_radius)
        #print(len(swarm_pos))
        #print(len(p_in_p))
        #print(p_in_p)
        swarm_pos_in_poly = swarm_pos[p_in_p]
        swarm_values_in_poly = swarm_values[p_in_p]
        swarm_values_full_in_poly = swarm_values_nec[p_in_p]
        swarm_date_in_poly, swarm_time_in_poly = swarm_date[p_in_p], swarm_time[p_in_p]
        print('len cut pos', len(swarm_pos_in_poly))
        vline_min_path = True
        if len(swarm_pos_in_poly) == 0:
            STATUS = 0
    else:
        """if proj_type == 'ortho_n':
            p_in_p = data_lat_up(swarm_pos, lat=0, hemisphere='N')    # bool
            swarm_pos_in_poly = swarm_pos[p_in_p]
            swarm_values_in_poly = swarm_values[p_in_p]
        elif proj_type == 'ortho_s':
            p_in_p = data_lat_up(swarm_pos, lat=0, hemisphere='S')     # bool
            swarm_pos_in_poly = swarm_pos[p_in_p]
            swarm_values_in_poly = swarm_values[p_in_p]"""

        swarm_pos_in_poly = swarm_pos
        swarm_values_in_poly = swarm_values
        swarm_values_full_in_poly = swarm_values_nec
        swarm_date_in_poly, swarm_time_in_poly = swarm_date, swarm_time
        vline_min_path = False
    d2txt.append(swarm_pos_in_poly, name='SWARM_pos')
    # отрисовка временных линий на пролете спутника
    if vline_min_path:
        swarm_dt = []
        for date, time in zip(swarm_date, swarm_time):
            swarm_dt.append(date + 'T' + time)
        vline_dt = np.array(swarm_dt)
    else:
        vline_dt = None
    sword.draw_swarm_path(swarm_pos[:, :2], points_time=vline_dt)

    # отрисовка вектора (X, Y, n, e) или (X, Y, |n-x|, |e-y|)
    if draw_vector_diff:
        #vector_components = swarm_values_nec[:, :2]   # n,e component
        vector_components = swarm_values_full_in_poly[:, :2]   # n,e component
        chaosm = CHAOS7(swarm_set=[swarm_liter, swarm_pos_in_poly, swarm_date_in_poly, swarm_time_in_poly, vector_components])
        vector_subtraction = chaosm.get_swarm_chaos_vector_subtraction()
        #vector_subtraction = swarm_egrf_vector_subtraction(swarm_pos_in_poly, vector_components, swarm_date_in_poly)

        swarm_values_in_poly = vector_subtraction[:, 0]     # dd
        vector_components = vector_subtraction[:, (1, 2)]   # dx, dy
        legend_label = 'Dd'

        # convert to geomagnetic coords
        #swarm_pos_in_poly = geo2mag(swarm_pos_in_poly, swarm_date_in_poly)

        d2txt.append(swarm_values_in_poly, name='SWARM_CHAOS7')
        ax_label += '|SWARM-CHAOS7|'
    elif measure_mu:
        swarm_values_in_poly = get_measure_mu(swarm_values_in_poly)
        ax_label += ' measure_mu '
        legend_label = 'mu'
        d2txt.append(swarm_values_in_poly, name='measure_mu')
    else:
        d2txt.append(swarm_values_in_poly, name='SWARM')

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
    if STATUS == 1:

        sword.draw_swarm_scatter(swarm_pos_in_poly[:, :2], swarm_values_in_poly, custom_label=legend_label,
                                 annotate=annotate_sw_value_bool)  # отрисовка значение точек на орбите
        if draw_vector_diff:
            sword.draw_vector(swarm_pos_in_poly, B=vector_components)

    # конвертация figure matplotlib в PIL image для stacker.py
    if cut_swarm_value_bool == True or proj_extend_loc is not None:
        sword.set_axis_label(ax_label, zoom_axis=True)
    else:
        sword.set_axis_label(ax_label)

    if STATUS == 1:
        if txt_out:
            d2txt.array_to_column()
            d2txt.save_columns_to_txt()
            out = d2txt.id
        else:
            out = sword.fig_to_PIL_image()
    if STATUS == 0:
        out = 'NO SWARM DATA IN TIME INTERVAL FROM %s TO %s IN LOCATION lat:%s... %s lon:%s... %s\n\n' \
              'hint: expand specified location or decrease delta' % (
           swarm_date[0]+'T'+swarm_time[0], swarm_date[-1]+'T'+swarm_time[-1], proj_extend_loc[0], proj_extend_loc[1],
                                                                       proj_extend_loc[2], proj_extend_loc[3],)
    return (STATUS, out)


def get_plot_im(swarm_sets, labels, auroral, channel, delta, measure_mu, ground_station=None, txt_out=False):
    #swarm_liter, swarm_pos, swarm_date, swarm_time, swarm_values = swarm_info[0]
    #from_date, to_date = swarm_info[1], swarm_info[2]
    # инициация рендера
    sword = SWORD()

    #   отрисовка графика
    draw_list = []
    auroral_list = []
    d2txt = Data2Text(SWARM_liter='A')

    for i, swarm_set in enumerate(swarm_sets):
        position_list, date_list, time_list = swarm_set[1], swarm_set[2], swarm_set[3]
        d2txt.annotate = 'from %s to %s' % (date_list[0]+'T'+time_list[0], date_list[-1]+'T'+time_list[-1])
        #   find swarm value shape size
        if len(np.array(swarm_set[4]).shape) == 1:
            shape_size = 1
        else:
            shape_size = np.array(swarm_set[4]).shape[1]

        #   append value[channel] or value
        if channel is not None:
            try:
                value = np.array(swarm_set[4])[:, channel]
            except:
                value = swarm_set[4]
            if measure_mu:
                value = get_measure_mu(value)
            label = labels[i]
            draw_list.append([label, date_list, time_list, value, position_list])
            d2txt.SWARM_liter = label.rsplit('-')[1][0]
            d2txt.SWARM_channel = ['N', 'E', 'C'][channel]
            d2txt.append(position_list, 'SWARM_pos')
            d2txt.append(value, 'SWARM')

        if channel is None:
            if shape_size == 1:
                value = swarm_set[4]
                if measure_mu:
                    value = get_measure_mu(value)
                label = labels[i]
                draw_list.append([label, date_list, time_list, value, position_list])
                d2txt.SWARM_liter = label.rsplit('-')[1][0]
                d2txt.SWARM_channel = 'FAC2'
                d2txt.append(position_list, 'SWARM_pos')
                d2txt.append(value, 'SWARM')
            else:
                for ch in range(3):
                    label = labels[ch]
                    value = np.array(swarm_set[4])[:, ch]
                    if measure_mu:
                        value = get_measure_mu(value)
                    draw_list.append([label, date_list, time_list, value, position_list])
                    d2txt.SWARM_liter = label.rsplit('-')[1][0]
                    d2txt.SWARM_channel = str(channel)
                    d2txt.append(position_list, 'SWARM_pos')
                    d2txt.append(value, 'SWARM')
        if auroral:
            auroral_near_swarm = get_nearest_auroral_point_to_swarm(swarm_set)
            auroral_list.append(auroral_near_swarm)
            d2txt.auroral_shp = 'near'
            d2txt.append(auroral_near_swarm, 'auroral')



    status = 1
    if status == 1:
        if txt_out:
            d2txt.array_to_column()
            d2txt.save_columns_to_txt()
            out = d2txt.id
        else:
            sword.draw_plot(draw_list, auroral_list, delta, draw_station=ground_station)
            out = sword.fig_to_PIL_image()
    return (status, out)

