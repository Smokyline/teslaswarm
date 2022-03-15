import numpy as np

from engine.sword import SWORD
from tools.data_foo import *
from tools.ionomodel_foo import get_ionomodel_surf_file

def get_proj_image(swarm_info, proj_type,
                   ionomodel_param=None, draw_auroral_s=None,
                   draw_auroral_n=None, draw_shape=None, draw_vector_diff=False,
                   intermag_observ_code=None, measure_mu=False, mag_grid_coord=False,
                   cut_swarm_value_bool=False, cut_obs_swarm_value_bool=False, swarm_poly_loc=None, proj_extend_loc=None, annotate_sw_value_bool=False, cut_deg_radius=5):
    swarm_liter, swarm_pos, swarm_date, swarm_time, swarm_values_nec = swarm_info[0]   # swarm_set
    from_date, to_date, swarm_channel = swarm_info[1], swarm_info[2], swarm_info[3]
    custom_value_name = None


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
        ax_label += 'loc:lat{%0.2f:%0.2f},lon{%0.2f:%0.2f}' % (proj_extend_loc[0], proj_extend_loc[1],
                                                   proj_extend_loc[2], proj_extend_loc[3],)

    ######################################################################################################
    #   инициация рендера, указание проекции, и если требуется extend_loc - приближение конкретной области
    sword = SWORD(proj_type=proj_type, extend=proj_extend_loc)

    #   отрисовка на проекции иономодели
    if ionomodel_param is not None:
        s_fac_surf = get_ionomodel_surf_file(ionomodel_param)
        sword.draw_ionomodel(s_fac_surf, type=['south', 'seismic'])
        #TODO iono
        ax_label += ' iono_s+'
        #        ax_label += 'iono_n+'

    #   отрисовка на проекции аврорального овала
    if draw_auroral_s is not None:
        x, y, z = get_auroral_flux(draw_auroral_s, hemishpere='S')  # f(dt)
        sword.draw_avroral_oval([x,y,z], hemisphere='south')
        #sword.draw_mag_coord_lines(swarm_date[0])
        ax_label += ' auroral+'
    if draw_auroral_n is not None:
        x, y, z = get_auroral_flux(draw_auroral_n, hemishpere='N')  # f(dt)
        sword.draw_avroral_oval([x,y,z], hemisphere='north')
        #sword.draw_mag_coord_lines(swarm_date[0])
        ax_label += ' auroral+'


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
    if len(np.array(swarm_values_nec).shape) > 1 and swarm_channel is not None:
        swarm_values = swarm_values_nec[:, swarm_channel]
        legend_label = ['dBn', 'dBe', 'dBd'][swarm_channel]
    else:
        swarm_values = swarm_values_nec  # fac2
        legend_label = 'fac2'
    if measure_mu:
        swarm_values = get_measure_mu(swarm_values)
        ax_label += ' measure_mu '
        legend_label = 'mu'



    # выбор точек измерений SWARM в указанном полигоне
    # lat, lon
    """
    poly_points = [i1, i2, i3, i4]
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
        elif intermag_observ_code is not None:
            print('cut swarm_pos near obs')
            #p_in_p, poly = get_swarm_value_near_obs(swarm_pos, intermag_observ_code, proj_type)
            p_in_p = get_position_near_point(swarm_pos[:, :2], intermag_observ_code, degr_radius=cut_deg_radius)
        print(len(swarm_pos))
        print(len(p_in_p))
        print(p_in_p)
        swarm_pos_in_poly = swarm_pos[p_in_p]
        swarm_values_in_poly = swarm_values[p_in_p]
        swarm_values_full_in_poly = swarm_values_nec[p_in_p]
        swarm_date_in_poly = swarm_date[p_in_p]
        print('len cut pos', len(swarm_pos_in_poly))
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
        swarm_date_in_poly = swarm_date




    # отрисовка вектора (X, Y, n, e) или (X, Y, |n-x|, |e-y|)
    if draw_vector_diff:
        #vector_components = swarm_values_nec[:, :2]   # n,e component
        vector_components = swarm_values_full_in_poly[:, :2]   # n,e component
        vector_subtraction = swarm_egrf_vector_subtraction(swarm_pos_in_poly, vector_components, swarm_date_in_poly)
        swarm_values_in_poly = vector_subtraction[:, 0]     # dd
        vector_components = vector_subtraction[:, (1, 2)]   # dx, dy
        legend_label = 'Dd'

        # convert to geomagnetic coords
        #swarm_pos_in_poly = geo2mag(swarm_pos_in_poly, swarm_date_in_poly)
        sword.draw_vector(swarm_pos_in_poly, B=vector_components)



    """    
        if convert_coord is not None:
            # конвертация географических координат в dest_sys
            swarm_pos_in_poly = convert_coord_system(swarm_pos_in_poly, dest_sys=convert_coord)
            ax_label += 'coord_sys:%s ' % convert_coord
        vector_components = swarm_egrf_vector_subtraction(swarm_pos_in_poly, vector_components, swarm_date)
        sword.draw_vector(swarm_pos_in_poly, B=vector_components)"""
    #sword.draw_mag_coord_lines(swarm_date_in_poly[0], geomag_pole=False)
    #sword.draw_swarm_path(swarm_pos_in_poly[:, :2])
    vline_min_path = True
    if vline_min_path:
        swarm_dt = []
        for date, time in zip(swarm_date, swarm_time):
            swarm_dt.append(date+'T'+time)
        sword.draw_swarm_path(swarm_pos[:, :2], points_time=np.array(swarm_dt))
    else:
        sword.draw_swarm_path(swarm_pos[:, :2])
    # отрисовка точек измерений swarm в в указанном полигоне
    # если полигона нет - отрисовка всех (swarm_pos_in_poly = swarm_pos)
    sword.draw_swarm_scatter(swarm_pos_in_poly[:, :2], swarm_values_in_poly, custom_label=legend_label,
                                 annotate=annotate_sw_value_bool)  # отрисовка значение точек на орбите


    # конвертация figure matplotlib в PIL image для stacker.py
    if cut_swarm_value_bool == True or proj_extend_loc is not None:
        sword.set_axis_label(ax_label, zoom_axis=True)
    else:
        sword.set_axis_label(ax_label)
    im = sword.fig_to_PIL_image()
    return im


def get_plot_im(swarm_sets, labels, include, channel, delta, ground_station=None):
    #swarm_liter, swarm_pos, swarm_date, swarm_time, swarm_values = swarm_info[0]
    #from_date, to_date = swarm_info[1], swarm_info[2]
    # инициация рендера
    sword = SWORD()

    #   отрисовка графика
    draw_list = []

    for i, swarm_set in enumerate(swarm_sets):
        position_list, date_list, time_list = swarm_set[1], swarm_set[2], swarm_set[3]

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
            label = labels[i]
            draw_list.append([label, date_list, time_list, value, position_list])

        if channel is None and shape_size == 1:
            value = swarm_set[4]
            label = labels[i]
            draw_list.append([label, date_list, time_list, value, position_list])

        elif channel is None and shape_size > 1:
            for ch in range(3):
                label = labels[ch]
                value = np.array(swarm_set[4])[:, ch]
                draw_list.append([label, date_list, time_list, value, position_list])

    sword.draw_plot(draw_list, include, delta, draw_station=ground_station)

    im = sword.fig_to_PIL_image()
    # im.save(STATIC_OS_PATH + '/media/images/2.png')
    return im
