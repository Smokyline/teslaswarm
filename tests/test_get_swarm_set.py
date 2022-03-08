import os
import sys

# Find code directory relative to our directory
sys.path.append(os.getcwd())

def test_get_swarm_data(fac2_mod=False, chaos_mod=False):
        try:
            from tools.sql_table import get_swarm_set
            from tools.dt_foo import decode_str_dt_param
            sw_liter = 'A'
            from_date = '2019-5-13T22:59:59'
            to_date = '2019-5-14T22:59:59'
            from_date, to_date = decode_str_dt_param(from_date), decode_str_dt_param(to_date)
            delta = 300

            swarm_set = get_swarm_set(sw_liter, from_date, to_date, delta, fac2_mod, chaos_mod)
            print('fac2_mod=%s, chaos_mod=%s' % (fac2_mod, chaos_mod))
            print('swarm_liter:', swarm_set[0])
            print('swarm_position:', swarm_set[1])
            print('swarm_date:', swarm_set[2])
            print('swarm_time', swarm_set[3])
            print('swarm_values:', swarm_set[4])
            print('swarm values len:', len(swarm_set[4]))

            print('\nswarm set successfully created')
        except Exception as e:
            #print(e)
            print('swarm set create is failed!')


if __name__ == '__main__':
    fac2_mod = False
    chaos_mod = False
    try:
        if str(sys.argv[1]) == 'fac2':
            fac2_mod = True
        elif str(sys.argv[1]) == 'chaos':
            chaos_mod = True
    except Exception as e:
        pass
    test_get_swarm_data(fac2_mod, chaos_mod)