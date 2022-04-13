from tools.supermag_api import SuperMAGGetIndices, supermag_testing, SuperMAGGetData

#out = SuperMAGGetIndices(logon='pilipenko', start='2019-10-15', extent=3600)
#print(out)


#supermag_testing('pilipenko')
(status,mydata1b) = SuperMAGGetData('pilipenko','2019-11-15T10:40',3600,'all','HBK',FORMAT='list')
print(status)
print(mydata1b)