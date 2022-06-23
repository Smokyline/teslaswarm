from tools.supermag_api import SuperMAGGetIndices, supermag_testing, SuperMAGGetData

#out = SuperMAGGetIndices(logon='pilipenko', start='2019-10-15', extent=3600)
#print(out)


#supermag_testing('pilipenko')
(status,mydata1b) = SuperMAGGetData('pilipenko','2017-09-18T02:30',3600,'all','T47',FORMAT='list')
print(status)
print(mydata1b)