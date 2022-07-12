import numpy as np

def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)


def angle_between(v1, v2):
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))


def rotate(v1, v2, xold):
    r = np.cross(v1, v2)
    if np.linalg.norm(r) == 0:
        print('Входные вектора коллинеарны')
        quit(-1)
    e = unit_vector(r)

    phi = -angle_between(v1, v2)

    cphi = np.cos(phi)
    sphi = np.sin(phi)

    xnew = xold * cphi + e * (1 - cphi) * np.dot(e, xold) + np.cross(e, xold) * sphi

    return xnew



"""if __name__ == '__main__':
    v1 = np.array([0, 0, 1])    #   старая ось Z    # среднее за t замеров
    v2 = np.array([0, 1, 1])    #   направление новой оси Z     # модель Циганенко t96

    xold = np.array([0, 0, 1])  #  преобразуемый вектор     # замер спутника

    xnew = rotate(v1, v2, xold)
    print(xnew)"""