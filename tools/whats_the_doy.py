import datetime

dt = datetime.date(2017, 9, 28)

day_of_year = dt.timetuple().tm_yday
print(dt, 'doy:', day_of_year)