import numpy as np
from math import cos, sin, pi

def dthtr_to_phy(s):
    xd = np.array([0, s[2]])
    xi = np.array([s[0]*cos(s[1] - pi/2), s[2] - s[0]*sin(s[1] - pi/2)])
    return xd, xi

def phy_to_DI(x, s):
    delta = -(s[1] - pi/2)
    C = np.array([[cos(delta), sin(delta)], [-sin(delta), cos(delta)]])
    xs = C.dot(x - np.array([0, s[2]]))
    return xs