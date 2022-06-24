import os
import sys
from teslaswarm.settings import *

# Find code directory relative to our directory
sys.path.append(os.getcwd())


def test_lib_load():
    try:
        import numpy
        import pandas
        import matplotlib
        import django
        import pymysql
        import cartopy  # request Proj 8.0.0 https://stackoverflow.com/questions/53697814/using-pip-install-to-install-cartopy-but-missing-proj-version-at-least-4-9-0
        import pyproj
        import astropy
        import shapely
        import geog     # pip
        import geospacepy   # pip below
        import spacepy  # conda install networkx h5py   pip install ffnet
        import aacgmv2
        import pyproj
        import astropy
        import chaosmagpy
        import ovationpyme  # https://github.com/lkilcommons/OvationPyme https://github.com/lkilcommons/nasaomnireader   #https://spdf.gsfc.nasa.gov/pub/software/cdf/dist/cdf38_0/
        #import pyIGRF   # https://pypi.org/project/pyIGRF/#description #12
        import igrf # 13


        print('all libraries successfully load')
    except Exception as e:
        print(e)
        print('libraries load is failed')


def test_chaos7_model_load():
    try:
        from chaos7_model.chaos_model import CHAOS7
        chaos = CHAOS7(None)
        print('chaos7 model successfully load')
    except Exception as e:
        print(e)
        print('chaos7 model load is failed')


def test_sql_connect():
    try:
        from datetime import datetime
        from tools.sql_table import get_sql_response
        from tools.dt_foo import ut_dt_to_unix
        swarm_type = ['SWA', 'SWB', 'SWC']
        from_date = ut_dt_to_unix(datetime(2018, 5, 26, 00, 00, 00), out_type='str')
        to_date = ut_dt_to_unix(datetime(2018, 5, 26, 23, 00, 00), out_type='str')
        respond = get_sql_response(swarm_type[0], from_date, to_date, fac2_mod=False)
        print('sql table successfully load')
    except Exception as e:
        print(e)
        print('sql table load is failed')


def test_matlab_status():
    try:
        from teslaswarm.settings import BASE_DIR, MATLAB_PATH
        import subprocess

        cmd = [MATLAB_PATH, 'matlab -nodesktop -nosplash -r "meas_extr_SA(); exit"']
        p = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            cwd=BASE_DIR + '/DMA')
        while p.poll() is None:
            l = p.stdout.readline()  # This blocks until it receives a newline.
            print(l)
            print(p.stdout.read())
            print('matlab successfully load')
    except Exception as e:
        print(e)
        print('matlab model load is failed')


if __name__ == '__main__':
    #test_lib_load()
    #test_chaos7_model_load()
    test_sql_connect()
    #test_matlab_status()