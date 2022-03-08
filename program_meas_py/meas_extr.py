import program_meas_py.tlib.Drect as Drect

from program_meas_py.tlib.CCmpNew import CCmpNew

def get_measure_mu(data):
    #   Построение выпрямления
    halfwindow = 10
    rect = Drect.length(data, halfwindow)

    #   Построение меры экстремальности
    cmp = CCmpNew(1000, rect, 0)
    meas = cmp.get_level1m(1)  # out

    return meas