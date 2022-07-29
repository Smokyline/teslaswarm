import pandas as pd
import numpy as np
from tools.dt_foo import *
from tools.coordinates_convert import rotate_GEO_vector_to_MFA, GEO_to_GDZgeomag_projection
from teslaswarm.settings import STATIC_OS_PATH
from spacepy import coordinates as coord
from spacepy.time import Ticktock
from tools.sql_table import get_swarm_set, get_sql_response
from tools.dt_foo import decode_str_dt_param
sw_liter = 'A'
from_date = '2021-4-17T08:00:59'
to_date = '2021-4-17T09:45:59'
from_date, to_date = decode_str_dt_param(from_date), decode_str_dt_param(to_date)
delta = 15
fac2_mod = False

sql_response = get_sql_response(sw_liter,from_date, to_date, fac2_mod)
swarm_set_full = get_swarm_set(sql_response, sw_liter, None, fac2_mod)
swarm_set = get_swarm_set(sql_response, sw_liter, delta, fac2_mod)



sw_lat = swarm_set[1][:, 0]
sw_lon = swarm_set[1][:, 1]
sw_rad = swarm_set[1][:, 2]
sw_vector = swarm_set[4]
sw_dt = [decode_str_dt_param(d+'T'+t) for d, t in zip(swarm_set[2], swarm_set[3])]
MFA_field, field_MFA_lines, coord_GSM, _ = rotate_GEO_vector_to_MFA(sw_lat, sw_lon, sw_rad, sw_vector, sw_dt,
                                                  mfield_min=1,dir=-1)
geomag_trace_coord = GEO_to_GDZgeomag_projection(sw_lat, sw_lon, sw_rad, sw_dt)
print(MFA_field)
print(coord_GSM)
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation


"""
fig = plt.figure()
ax = fig.add_subplot(111)

xrange = np.arange(0, len(sw_vector))
plt.plot(xrange, sw_vector[:, 0], c='r', label='X_NEC')
plt.plot(xrange, MFA_field[:, 0], c='r', linestyle='--', label='X_MFA')

plt.plot(xrange, sw_vector[:, 1], c='g', label='Y_NEC')
plt.plot(xrange, MFA_field[:, 1], c='g', linestyle='--', label='Y_MFA')

plt.plot(xrange, sw_vector[:, 2], c='b', label='Z_NEC')
plt.plot(xrange, MFA_field[:, 2], c='b', linestyle='--', label='Z_MFA')

plt.plot(xrange, np.sqrt(sw_vector[:, 0]**2 + sw_vector[:, 1]**2 + sw_vector[:, 2]**2), c='k', label='|B|_NEC')
plt.plot(xrange, np.sqrt(MFA_field[:, 0]**2 + MFA_field[:, 1]**2 + MFA_field[:, 2]**2), c='k', linestyle='--', label='|B|_MFA')


plt.legend()
plt.grid(True)
ax.set_xticks(np.arange(len(sw_lat))[::225])
ax.set_xticklabels(np.round(coord_GSM[::225, 1], 2))
ax.set_ylabel('nT')
ax.set_xlabel('Lat GSM')
plt.show()

"""
fig = plt.figure(figsize=(10, 10))
ax = Axes3D(fig)
ax.scatter3D(0, 0, 0, c='k', marker='X')
def init():

    for i in range(len(coord_GSM)):
        """
                satgeo = coord.Coords([sw_rad[i]/6371.0, np.float(sw_lat[i]), np.float(sw_lon[i])], 'GEO', 'sph', ['Re', 'deg', 'deg'])
        satgeo.ticks = Ticktock(sw_dt[i])
        xsat, ysat, ysat = np.array(satgeo.convert('GEO', 'car').data)[0]

        ax.scatter3D(xsat, ysat, ysat, c='r', s=2, marker='o')

        gdz = coord.Coords([gdz_coord[i,2]/6371.0, np.float(gdz_coord[i,0]), np.float(gdz_coord[i, 1])], 'GDZ', 'sph', ['Re', 'deg', 'deg'])
        satgeo.ticks = Ticktock(sw_dt[i])
        xgdz, ygdz, ygdz = np.array(satgeo.convert('GEO', 'car').data)[0]
        ax.scatter3D(xgdz, ygdz, ygdz, c='k', s=2, marker='x')

        """



        xgsm, ygsm, zgsm = coord_GSM[i]
        xx, yy, zz = field_MFA_lines[i]
        if i == 0:
            ax.scatter3D(xgsm, ygsm, zgsm, c='b', s=2)
        else:
            ax.scatter3D(xgsm, ygsm, zgsm, c='k', s=2)
        #ax.scatter3D(xx[-1], yy[-1], zz[-1], c='g')
        #ax.plot3D(xx, yy, zz, c='b', lw=1.3)

        #xgeo,ygeo,zgeo = np.array(cgeo.convert('GEO', 'car').data)[0]
        #xgeo,ygeo,zgeo = geopack.sphcar(sw_rad[i]/6371.2, np.deg2rad(sw_lat[i]), np.deg2rad(sw_lon[i]), 1)
        #xgsm, ygsm, zgsm = geopack.geogsm(xgeo, ygeo, zgeo, 1)

        cgeo = coord.Coords([1, np.float(sw_lat[i]), np.float(sw_lon[i])], 'GEO', 'sph', ['Re', 'deg', 'deg'])
        cgeo.ticks = Ticktock(sw_dt[i])
        xgsm0, ygsm0, zgsm0 = np.array(cgeo.convert('GSM', 'car').data)[0]
        ax.scatter3D(xgsm0, ygsm0, zgsm0, c='k', s=2, marker='x')

        cgeo = coord.Coords([6471.0/6371.0, np.float(sw_lat[i]), np.float(sw_lon[i])], 'GEO', 'sph', ['Re', 'deg', 'deg'])
        cgeo.ticks = Ticktock(sw_dt[i])
        xgsm0, ygsm0, zgsm0 = np.array(cgeo.convert('GSM', 'car').data)[0]
        ax.scatter3D(xgsm0, ygsm0, zgsm0, c='b', s=2, marker='x')

        xgdz, ygdz, zgdz = geomag_trace_coord[i]
        ax.scatter3D(xgdz, ygdz, zgdz, c='r', s=2, marker='x')


        ax.plot3D(xx, yy, zz, lw=1.3)
        ax.plot3D([xgsm, 0], [ygsm, 0], [zgsm, 0], c='k', lw=0.4, alpha=0.6)
    return fig,

init()
plt.show()

#plt.clf()
#ax.scatter3D(0, 0, 0, c='k', marker='X')

def animate(i):
    ax.view_init(elev=10., azim=i)
    return fig,

# Animate
anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=360, interval=20, blit=True)
# Save
anim.save(STATIC_OS_PATH + '\\media\\videos\\basic_animation.mp4', fps=30,)


