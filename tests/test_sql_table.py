import os
import sys

# Find code directory relative to our directory
sys.path.append(os.getcwd())

def test_sql_connect(fac2_mod):
    try:
        from datetime import datetime
        import numpy as np
        import pandas as pd
        from tools.sql_table import get_sql_response, ut_dt_to_unix
        from teslaswarm.settings import STATIC_OS_PATH


        swarm_type = ['SWA', 'SWB', 'SWC']
        from_date = ut_dt_to_unix(datetime(2017, 9, 7, 23, 00, 00), out_type='str')
        to_date = ut_dt_to_unix(datetime(2017, 9, 8, 1, 00, 00), out_type='str')

        respond = get_sql_response(swarm_type[0], from_date, to_date, fac2_mod)
        print(pd.DataFrame(respond))

        respond = np.array(respond)
        print('len before remove nines', len(respond))
        # data = data_reduction(np.array(respond), delta=1, fac2_mod=True)
        idx999 = np.where(respond[:, 3] != 999)[0]
        respond = respond[idx999]
        print('len after remove nines', len(respond))
        print('max value: %s' % np.max(respond[:, 3]))
        print('min value: %s' % np.min(respond[:, 3]))


        pd.DataFrame(respond).to_csv(STATIC_OS_PATH + '/data/dataset_fac2_1day.csv', index=False, header=False, )
        df = pd.read_csv(STATIC_OS_PATH + '/data/dataset_fac2_1day.csv', header=None)
        print(df)

        print('sql table successfully load')
    except Exception as e:
        print(e)
        print('sql table load failed!')

if __name__ == '__main__':
    fac2_mod = False
    try:
        if str(sys.argv[1]) == 'fac2':
            fac2_mod = True
    except Exception as e:
        pass
    test_sql_connect(fac2_mod)

