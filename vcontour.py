import numpy as np
from math import cos, acos, sqrt, sin, pi
from scipy.optimize import NonlinearConstraint, minimize

from Config import Config
from coords import dthtr_to_phy

r = Config.CAP_RANGE
a = Config.VD/Config.VI

def dominant_region(x, s):
    xd, xi = dthtr_to_phy(s)
    return a*np.linalg.norm(x-xi) - (np.linalg.norm(x-xd) - r)

def target(x, R=Config.TAG_RANGE):
    return np.linalg.norm(x) - R

def tangent(x):
    xd, xi = x[0:2], x[2:4]
    def dominant_region(x, xi=xi, xd=xd):
        return a*np.linalg.norm(x-xi) - (np.linalg.norm(x-xd) - r)
    def target_square(x):
        return x[0]**2 + x[1]**2
    on_dr = NonlinearConstraint(dominant_region, -np.inf, 0)
    sol = minimize(target_square, xi, constraints=(on_dr,))
    return sol.x, sqrt(target_square(sol.x)), np.linalg.norm(sol.x - x[2:4])/Config.VI

def min_dist_to_target(x):
    return tangent(x)[1] - Config.TAG_RANGE

def time_to_cap(x):
    return tangent(x)[2]
#def intersect(s, R=Config.TAG_RANGE):
#    
#    def residual(x, s=s):
#        in_dr = dominant_region(x, s)
#        in_tg = target(x, R=R)
#        return in_dr**2 + in_tg**2
#    
#    d_ = a/(1+a)*(s[0]-r)+r
#    x0 = np.array([d_*cos(s[1] - pi/2), s[2] - d_*sin(s[1] - pi/2)])
#    sol = minimize(residual, x0)
#    
#    return sol.x, residual(sol.x, s=s)
#
#def tangent(d, r1, R=Config.TAG_RANGE):
#    
#    if r1 >= R + r:
#        lb, ub = 0, 0.99*pi
#    else:
#        if d >= r1 + R:
#            lb, ub = 0, 0.99*pi
#        else:
#            lb, ub = 0, pi - acos((r1**2 + d**2 - R**2)/(2*d*r1))
#                              
#    while ub - lb > 1e-6:
#        tht = (lb + ub)/2
##        x_ub, res_ub = intersect(np.array([d, ub, r1]), R=R)
##        x_lb, res_lb = intersect(np.array([d, lb, r1]), R=R)
#        x, res = intersect(np.array([d, tht, r1]), R=R)
##        print('[%.3f, %.3f, %.3f]'%(lb, tht, ub), '[%.6f, %.6f, %.6f]'%(res_lb, res, res_ub))
#        if res < 1e-10:
#            ub = tht
#        else:
#            lb = tht
#        
#    return (lb + ub)/2
#        
#def const_v(R, r1=1.5*Config.TAG_RANGE):
#    
#    if r1 >= R + r:
#        d_min = (r1 - R) - (r1 - R - r)/a
#        d_max = (r1 + R) + (r1 + R - r)/a
#    else:
#        d_min = r
#        d_max = (r1 + R) + (r1 + R - r)/a
#    ds = np.linspace(d_min, d_max, 25)
#
#    ss = []
#    for d in ds:
#        tht = tangent(d, r1, R=R)
#        ss.append(np.array([d, tht, r1]))
#    #    print(tht)
#    return ss
    
    
    
    
    
    
    
    
    
    
    
    