import os
import sys
import subprocess
import platform
import numpy as np
"""
NEED TO ADD RE_1_FAC to C:\\RE_1_FAC\\fac_data\\
"""

# Find code directory relative to our directory
def run_exe_program():
    if platform.system() == 'Linux' or platform.system() == 'Darwin':
            sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
            out = subprocess.call(
                'cd E_FIELD_FAC_MODEL_FCP \nwine coef_ut.exe \nwine fac_bt_season.exe \nwine tph.exe \nwine ris_surfer_sigma.exe \nwine ris_surfer_fac.exe',
                shell=True)
    elif platform.system() == 'Windows':
        path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '\\E_FIELD_FAC_MODEL_FCP')
        subprocess.call('.\\coef_ut.exe', shell=True, cwd=path)
        subprocess.call('.\\fac_bt_season.exe', shell=True, cwd=path)
        subprocess.call('.\\tph.exe', shell=True, cwd=path)
        subprocess.call('.\\ris_surfer_sigma.exe', shell=True, cwd=path)
        subprocess.call('.\\ris_surfer_fac.exe', shell=True, cwd=path)


def set_param():
    keys = ['bz', 'by', 'doy', 'kp', 'f107', 'ut']
    format_param = [' Bz=    v', ' By=    v', ' DOY=    v', ' Kp=    v', ' F107=  v', ' UT=    v']
    param_dict = {
        'type': 'fac',
        'hem': 'n',
        'bz': '6',
        'f107': '100',
        'by': '6',
        'doy': '300',
        'kp': '1',
        'ut': '1'
    }
    BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

    f = open(BASE_DIR + '/E_FIELD_FAC_MODEL_FCP/input data', 'w')
    for i, p in enumerate(format_param):
        line = p.replace('v', str(np.around(float(param_dict[keys[i]]), 1)))
        f.write(line + '\n')
    f.write('   BT=   %s\n' % (np.sqrt(float(param_dict['bz']) ** 2 + float(param_dict['by']) ** 2).round(1)))
    f.write('   IMF clock angle=  45.0')
    f.close()


if __name__ == '__main__':
    print('please, launch this script from project dir \n for this use $ python testing/test_iono.py')
    set_param()
    run_exe_program()