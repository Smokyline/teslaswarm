import numpy as np
from chaosmagpy import load_CHAOS_matfile
from chaosmagpy.model_utils import synth_values
from teslaswarm.settings import CHAOS_PATH
from astropy.time import Time


class CHAOS7():
    def __init__(self, swarm_set):
        self.swarm_set = swarm_set
        self.chaos7_mat_model = load_CHAOS_matfile(CHAOS_PATH)
        self.get_coef()
    def dt_unix_to_mjd(self, times):
        t = Time(times, format='iso', scale='utc')
        t.format = 'mjd'
        time_mjd_ad = np.array(t.value).astype(float) - 51544
        return time_mjd_ad

    def get_coef(self):
        swarm_liter, swarm_pos, swarm_date, swarm_time, swarm_values = self.swarm_set
        sw_n, sw_e, sw_c, sw_f = swarm_values[:, 0], swarm_values[:, 1], swarm_values[:, 2], swarm_values[:, 3]
        # B_r = -Z; B_phi = Y; B_theta = -X
        theta = 90. - swarm_pos[:, 1]  # colat deg
        phi = swarm_pos[:, 0]  # deg
        radius = swarm_pos[:, 2]  # radius in km

        time = self.dt_unix_to_mjd(
            [str(a) + " " + str(b) for a, b in zip(swarm_date, swarm_time)])  # time in modified Julian date 2000

        # computing core field
        coeffs = self.chaos7_mat_model.synth_coeffs_tdep(time)  # SV max. degree 16
        self.B_radius, self.B_theta, self.B_phi = synth_values(coeffs, radius, theta, phi)
        # B_radius, B_theta, B_phi = self.chaos7_mat_model(time, radius, theta, phi)
        self.SW_C_theta, self.SW_C_phi = (sw_n * -1) - self.B_theta, sw_e - self.B_phi

    def get_chaos_swarm_diff(self):
        return self.SW_C_theta, self.SW_C_phi
    def get_chaos_theta_phi(self):
        return self.B_radius, self.B_theta, self.B_phi
    # TODO

