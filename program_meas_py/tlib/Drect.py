import numpy as np

def energy(data, delta=5):

    midd = np.zeros_like(data)
    rect = np.zeros_like(data)
    kpoint = len(data)

    for i in range(kpoint):
        i1 = max(0, i-delta)
        i2 = min(kpoint-1, i+delta)
        midd[i] = np.mean(data[i1:i2])

    for i in range(kpoint):
        i1 = max(0, i - delta)
        i2 = min(kpoint - 1, i + delta)
        tmp = data[i1:i2] - midd[i]
        rect[i] = np.mean(tmp**2)

    return rect

def length(data, delta=5):

    rect = np.zeros_like(data)
    kpoint = len(data)

    for i in range(kpoint):
        i1 = max(0, i-delta)
        i2 = min(kpoint, i+delta+1)
        rect[i] = np.sum(abs(data[i1+1:i2]-data[i1:i2-1]))

    return rect

def minmax(data, delta=5):

    rect = np.zeros_like(data)
    kpoint = len(data)

    for i in range(kpoint):
        i1 = max(0, i-delta)
        i2 = min(kpoint-1, i+delta)
        rect[i] = np.max(data[i1:i2])-np.min(data[i1:i2])

    return rect