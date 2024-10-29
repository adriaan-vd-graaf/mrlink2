import numpy as np

def dxx_function_v0(th, lam, cX, cY, n_x, n_y):
    a = th[0]
    tX = abs(th[1])
    tY = abs(th[2])

    return 1. / (((a ** 2 * n_y + n_x) * lam + tX) - a ** 2 * n_y ** 2 * (lam ** 2) / (n_y * lam + tY))


def dxx_function_v2(th, lam, cX, cY, n_x, n_y):
    a = th[0]
    tX = abs(th[1])
    tY = abs(th[2])

    return 1. / (np.exp(np.log(a ** 2 * n_y + n_x) + np.log(lam)) + tX - np.exp(np.log(a ** 2 * n_y ** 2 * (lam ** 2)) - np.log(n_y * lam + tY)))


def test_dXX():
    th = np.array([1.7976931348623158e+306, 1.7976931348623158e-308, 0.00])
    lam = np.array([np.nextafter(0, 1), 2, 3])
    cX = np.array([-0, 1, 1])
    cY = np.array([2, 2, 2])

    n_x = 10
    n_y = 10
    print()

    print(dxx_function_v0(th, lam, cX, cY, n_x, n_y))
    print(dxx_function_v2(th, lam, cX, cY, n_x, n_y))
