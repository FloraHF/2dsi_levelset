import numpy as np
from math import sqrt, atan2
from Config import Config
r = Config.CAP_RANGE
#R = Config.TAG_RANGE
a = Config.VD/Config.VI

def s_strategy(x):
    
    xt, _ = tangent(x)
    IT = np.concatenate((xt - x[2:4], np.zeros((1,))))
    DT = np.concatenate((xt - x[0:2], np.zeros(1,)))
    xaxis = np.array([1, 0, 0])
    phi = atan2(np.cross(xaxis, DT)[-1], np.dot(xaxis, DT))
    psi = atan2(np.cross(xaxis, IT)[-1], np.dot(xaxis, IT))
    
    return np.array([phi, psi])

#print(min_dist_to_target(np.array([0, 6, 5, 7])))