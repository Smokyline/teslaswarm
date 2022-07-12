import pandas as pd
import numpy as np
from tools.dt_foo import *
from tools.coordinates_convert import rotate_GEO_vector_to_MFA
from teslaswarm.settings import STATIC_OS_PATH

from tools.sql_table import get_swarm_set
from tools.dt_foo import decode_str_dt_param
sw_liter = 'A'
from_date = '2021-4-17T08:11:59'
to_date = '2021-4-17T09:45:59'
from_date, to_date = decode_str_dt_param(from_date), decode_str_dt_param(to_date)
delta = 1
fac2_mod = False

swarm_set = get_swarm_set(sw_liter, from_date, to_date, delta, fac2_mod)



sw_lat = swarm_set[1][:, 0]
sw_lon = swarm_set[1][:, 1]
sw_rad = swarm_set[1][:, 2]
sw_vector = swarm_set[4]
sw_dt = [decode_str_dt_param(d+'T'+t) for d, t in zip(swarm_set[2], swarm_set[3])]
MFA_field, field_MFA_lines, coord_GSM = rotate_GEO_vector_to_MFA(sw_lat, sw_lon, sw_rad, sw_vector, sw_dt)
print(MFA_field)
print(coord_GSM)
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation



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
        xgsm, ygsm, zgsm = coord_GSM[i]
        xx, yy, zz = field_MFA_lines[i]
        if i == 0:
            ax.scatter3D(xgsm, ygsm, zgsm, c='b', s=2)
        else:
            ax.scatter3D(xgsm, ygsm, zgsm, c='k', s=2)
        #ax.scatter3D(xx[-1], yy[-1], zz[-1], c='g')
        #ax.plot3D(xx, yy, zz, c='b', lw=1.3)
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



"""