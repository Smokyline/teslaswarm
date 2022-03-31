import pandas as pd
import numpy as np
from tools.dt_foo import *
from teslaswarm.settings import STATIC_OS_PATH
#https://github.com/spacecataz/supermag

from tools.supermag_api import SuperMAGGetIndices, supermag_testing, SuperMAGGetData

out = SuperMAGGetIndices(logon='pilipenko', start='2017-09-18', extent=3600)
print(out)
