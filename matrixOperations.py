def multiply(a, b):
    h = len(a)
    mid = len(b)
    w = len(b[0])
    return [[sum([a[i][k] * b[k][j] for k in range(mid)]) for j in range(w)] for i in range(h)]


def inverse(a):
    import numpy as np
    np.array(a)
    return np.linalg.inv(a).tolist()

