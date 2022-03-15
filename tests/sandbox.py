import pandas as pd
import numpy as np
from tools.dt_foo import *
#https://github.com/spacecataz/supermag

from tools.supermag_api import SuperMAGGetData, sm_grabme

#['AAA', 'ABG', 'ABK', 'AIA', 'AMD', 'AMS', 'AND', 'API', 'ARS', 'ASC', 'ASP', 'ATU', 'BBG', 'BDV', 'BEL', 'BEY', 'BFE', 'BFO', 'BJN', 'BLC', 'BMT', 'BOU', 'BOX', 'BRD', 'BRN', 'BRW', 'BSL', 'C02', 'C03', 'C04', 'C05', 'C06', 'C07', 'C10', 'C13', 'CAN', 'CBB', 'CDC', 'CKI', 'CLF', 'CMO', 'CNB', 'CPS', 'CSY', 'CTA', 'CUL', 'CYG', 'CZT', 'DAW', 'DED', 'DIK', 'DLT', 'DMC', 'DMH', 'DOB', 'DON', 'DOU', 'DRV', 'DUR', 'EBR', 'ESK', 'EYR', 'FCC', 'FHB', 'FRD', 'FRN', 'FSP', 'FUR', 'GCK', 'GDH', 'GHB', 'GHC', 'GIM', 'GNG', 'GUA', 'GUI', 'HAD', 'HAN', 'HER', 'HLP', 'HON', 'HOP', 'HRB', 'HRN', 'IGC', 'IPM', 'IQA', 'IRT', 'ISL', 'IVA', 'JAI', 'JAN', 'JCK', 'JCO', 'KAG', 'KAK', 'KAR', 'KDU', 'KEP', 'KEV', 'KHB', 'KIL', 'KIV', 'KLI', 'KMH', 'KNY', 'KOU', 'KTB', 'KZN', 'LER', 'LET', 'LON', 'LRM', 'LVV', 'LYC', 'LYR', 'LZH', 'M01', 'MAB', 'MAS', 'MAW', 'MBO', 'MCQ', 'MEA', 'MEK', 'MGD', 'MMB', 'MNK', 'MOS', 'MSR', 'MUO', 'NAD', 'NAL', 'NAQ', 'NCK', 'NEW', 'NGK', 'NOR', 'NUR', 'NVS', 'OTT', 'OUJ', 'PAF', 'PAG', 'PEG', 'PEL', 'PET', 'PGC', 'PHU', 'PIL', 'PIN', 'PKR', 'PPT', 'PST', 'RAL', 'RAN', 'RES', 'RIK', 'ROE', 'RVK', 'SBA', 'SBL', 'SFS', 'SHU', 'SIT', 'SJG', 'SMI', 'SOD', 'SOL', 'SON', 'SOR', 'SPG', 'SPT', 'STF', 'STJ', 'SUA', 'SUM', 'SVS', 'T16', 'T22', 'T23', 'T29', 'T31', 'T32', 'T36', 'T37', 'T41', 'T43', 'T44', 'T45', 'T47', 'T48', 'T49', 'T50', 'T51', 'T52', 'T54', 'T55', 'T56', 'TAB', 'TAM', 'TAR', 'TDC', 'TEO', 'THL', 'THY', 'TIK', 'TRO', 'TRW', 'TSU', 'TUC', 'TWN', 'UPN', 'UPS', 'VAL', 'VIC', 'VIZ', 'VNA', 'VOS', 'VSS', 'W01', 'W02', 'WIC', 'WNG', 'YAK', 'YKC', 'HUA']

answer = SuperMAGGetData(logon='pilipenko',start='2017-09-18T00:00:00', extent=86400,
                             station='T47', flagstring='',)

data = answer[1]
#print(answer[1].to_string())

format = 'geo'
i = 2
for n, e, c in zip(sm_grabme(data,'N', format), sm_grabme(data,'E', format), sm_grabme(data,'Z', format)):
    print(i, n, e, c)
    i+=1

#print(sm_grabme(data,'N','nez'))
#print(sm_grabme(data,'N','geo'))

