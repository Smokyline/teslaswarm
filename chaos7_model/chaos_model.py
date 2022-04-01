import numpy as np
from chaosmagpy import load_CHAOS_matfile
from chaosmagpy.model_utils import synth_values
from chaosmagpy.coordinate_utils import transform_points
from teslaswarm.settings import CHAOS_PATH
from astropy.time import Time
import cdflib
from tools.dt_foo import ut_dt_to_unix, decode_str_dt_param
from chaosmagpy.data_utils import mjd2000
import datetime
from spacepy import pycdf

class CHAOS7():
    def __init__(self, swarm_set):
        #   https://chaosmagpy.readthedocs.io/_/downloads/en/latest/pdf/

        self.swarm_set = swarm_set
        self.chaos7_mat_model = load_CHAOS_matfile(CHAOS_PATH)
        self.get_coef()
        self.cdf = pycdf.Library()
    def dt_unix_to_mjd(self, times):
        t = Time(times, format='iso', scale='utc')
        t.format = 'mjd'
        time_mjd_ad = np.array(t.value).astype(float) - 51544
        return time_mjd_ad

    def datetime_to_epoch(self, dt):
        """
        Converts a Python datetime to a CDF Epoch value

        Parameters
        ==========
        dt : :class:`datetime.datetime`
            date and time to convert

        Returns
        =======
        out : float
            epoch corresponding to dt

        See Also
        ========
        v_datetime_to_epoch
        """
        if dt.tzinfo != None and dt.utcoffset() != None:
            dt = dt - dt.utcoffset()
        dt.replace(tzinfo=None)
        micro = dt.microsecond % 1000
        if micro >= 500 and dt.year < 9999:
            dt += datetime.timedelta(0, 0, 1000)
        return self._library.computeEPOCH(dt.year, dt.month, dt.day, dt.hour,
                                          dt.minute, dt.second,
                                          int(dt.microsecond / 1000))
    def magfield_variation(self, n_swarm, e_swarm, x_chaos, y_chaos, ):
        """
             n same as x
             e same as y
        """
        FACT = 180. / np.pi
        x, y, = (n_swarm + x_chaos) / 2, (e_swarm + y_chaos) / 2,
        dx, dy, = (x_chaos - n_swarm) / 2, (y_chaos - e_swarm) / 2,
        h = np.sqrt(x * x + y * y)
        dd = (FACT * (x * dy - y * dx)) / (h * h)
        return dd, dx, dy

    def get_coef(self):
        #lat=swarm_pos[idx, 0], lon=swarm_pos[idx, 1],  alt=swarm_pos[idx, 2]
        #swarm_liter, swarm_pos, swarm_date, swarm_time, swarm_values = self.swarm_set
        swarm_liter, swarm_pos, swarm_date, swarm_time, vector_components = self.swarm_set
        sw_n, sw_e = vector_components[:, 0], vector_components[:, 1]
        # B_r = -Z; B_phi = Y; B_theta = -X
        theta = 90. - swarm_pos[:, 0]  # colat deg
        phi = swarm_pos[:, 1]  # deg
        radius = swarm_pos[:, 2]  # radius in km

        time = self.dt_unix_to_mjd(
            [str(a) + " " + str(b) for a, b in zip(swarm_date, swarm_time)])  # time in modified Julian date 2000

        # computing core field
        coeffs = self.chaos7_mat_model.synth_coeffs_tdep(time)  # SV max. degree 16
        self.B_radius, self.B_theta, self.B_phi = synth_values(coeffs, radius, theta, phi)
        # B_radius, B_theta, B_phi = self.chaos7_mat_model(time, radius, theta, phi)
        self.SW_C_theta, self.SW_C_phi = (sw_n * -1) - self.B_theta, sw_e - self.B_phi

    def get_swarm_chaos_vector_subtraction(self):
        model = load_CHAOS_matfile(CHAOS_PATH)

        swarm_liter, swarm_pos, swarm_date, swarm_time, vector_components = self.swarm_set
        sw_n, sw_e = vector_components[:, 0], vector_components[:, 1]

        theta = 90. - swarm_pos[:, 0]  # colat deg
        phi = swarm_pos[:, 1]  # deg
        radius = swarm_pos[:, 2]  # radius in km
        print(theta[:10])
        print(phi[:10])
        print(radius[:10])
        sw_dt = np.array([decode_str_dt_param(a+'T'+b) for a, b in zip(swarm_date, swarm_time)])    # datetime format
        time = self.cdf.v_datetime_to_epoch(sw_dt)  # CDF epoch format
        time = time / (1e3 * 3600 * 24) - 730485    # time in modified Julian date 2000

        #theta_gsm, phi_gsm = transform_points(theta, phi,
        #                                      time=time, reference='gsm')
        #index_day = np.logical_and(phi_gsm < 90, phi_gsm > -90)
        #index_night = np.logical_not(index_day)

        # complete forward computation: pre-built not customizable (see ex. 1)
        B_radius, B_theta, B_phi = model(time, radius, theta, phi)

        B = []
        for n, e, x, y in zip(sw_n, sw_e, B_theta, B_phi):
            dd, dx, dy = self.magfield_variation(n, e, x*-1, y)
            B.append([dd, dx, dy])
        return np.array(B)
        # compute field strength and plot together with data
        #F = np.sqrt(B_radius ** 2 + B_theta ** 2 + B_phi ** 2)

        #print('RMSE of F: {:.5f} nT'.format(np.std(F - F_swarm)))

    def get_chaos_swarm_diff(self):
        return self.SW_C_theta, self.SW_C_phi
    def get_chaos_theta_phi(self):
        return self.B_radius, self.B_theta, self.B_phi


