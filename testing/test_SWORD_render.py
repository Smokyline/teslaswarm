import os
import sys

# Find code directory relative to our directory
sys.path.append(os.getcwd())

from engine.sword_engine import SWORD
from tools.dt_foo import decode_str_dt_param

sw_liter = 'A'
from_date = '2019-5-13T22:59:59'
to_date = '2019-5-14T22:59:59'
from_date, to_date = decode_str_dt_param(from_date), decode_str_dt_param(to_date)
delta = 180

sword = SWORD(sw_liter, from_date, to_date, delta, swarm_channel=1, fac2_mod=True, chaos_mod=False, def_mod='')
id = sword.run_render()
