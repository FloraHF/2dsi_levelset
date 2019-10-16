import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

import numpy as np
from math import cos, sin, pi
from copy import deepcopy

from Config import Config
from vcontour import tangent, target, dominant_region, min_dist_to_target, time_to_cap
from coords import dthtr_to_phy

r = Config.CAP_RANGE
R = Config.TAG_RANGE
a = Config.VD/Config.VI

def get_defender(r1):
    x, y = [], []
    for tht in np.linspace(0, 2*pi, 50):
        x.append(r*cos(tht))
        y.append(r1 + r*sin(tht))
    # ax.plot(x, y, color='r')
    # ax.plot(0, r1, 'r.')
    return np.asarray(x), np.asarray(y)

def get_target(R=Config.TAG_RANGE, color='b', linestyle='dashed'):
    
    def get_tg():
        k = 1.1
        x = np.linspace(-k*R, k*R)
        y = np.linspace(-k*R, k*R)
        X, Y = np.meshgrid(x, y)
        T = np.zeros(np.shape(X))
        for i, (xx, yy) in enumerate(zip(X, Y)):
            for j, (xxx, yyy) in enumerate(zip(xx, yy)):
        #        print(dominant_region(np.array([xxx, yyy])))
                T[i,j] = target(np.array([xxx, yyy]), R=R)
        
        return {'X': X, 'Y': Y, 'T': T}
    
    tgt = get_tg()
    return tgt
    # CT = ax.contour(tgt['X'], tgt['Y'], tgt['T'], [0], linestyles=(linestyle,))
    # plt.contour(CT, levels = [0], colors=(color,), linestyles=(linestyle,))

def get_dr(s, color='b', skip=16):
    
    k = 2
    def get_dr(s):
        x = np.linspace(-k*R, k*R)
        y = np.linspace(-k*R, k*R)
        X, Y = np.meshgrid(x, y)
        D = np.zeros(np.shape(X))
        for i, (xx, yy) in enumerate(zip(X, Y)):
            for j, (xxx, yyy) in enumerate(zip(xx, yy)):
        #        print(dominant_region(np.array([xxx, yyy])))
                D[i,j] = dominant_region(np.array([xxx, yyy]), s)
        
        return {'X': X, 'Y': Y, 'D': D}

    return get_dr(s)
    
    # drs = []
    # for s in ss:
    #     drs.append(get_dr(s))
    
    # for i, dr in enumerate(drs):
    #     if i==1 or i==len(drs)-1 or i%skip==0:
    #         CD = ax.contour(dr['X'], dr['Y'], dr['D'], [0], linestyles='dashed')
    #         plt.contour(CD, levels = [0], colors=(color,), linestyles=('dashed',))
    #         xi = deepcopy(dthtr_to_phy(ss[i])[-1])
    #         ax.plot(xi[0], xi[1], '.', color=color)

def get_constd(r1, R=Config.TAG_RANGE, color='b', linestyle='solid'):
    def get_constd(r1):
        k = 4
        x = np.linspace(-k*Config.TAG_RANGE, k*Config.TAG_RANGE)
        y = np.linspace(-k*Config.TAG_RANGE, k*Config.TAG_RANGE)
        X, Y = np.meshgrid(x, y)
        C = np.zeros(np.shape(X))
        for i, (xx, yy) in enumerate(zip(X, Y)):
            for j, (xxx, yyy) in enumerate(zip(xx, yy)):
        #        print(dominant_region(np.array([xxx, yyy])))
                C[i,j] = min_dist_to_target(np.array([0, r1, xxx, yyy]))
        
        return {'X': X, 'Y': Y, 'C': C}
    vctr = get_constd(r1)
    return vctr
#     CC = ax.contour(vctr['X'], vctr['Y'], vctr['C'], [-2, 0, 2], linestyles=(linestyle,))
#     ax.clabel(CC, inline=True, fontsize=10)
# #    ax.clabel(CC, inline=True, fontsize=10)
#     ax.contour(vctr['X'], vctr['Y'], vctr['C'], levels = [-2, 0, 2], colors=('b',), linestyles=(linestyle,))

#    CC = ax.contour(vctr['X'], vctr['Y'], vctr['C'], [R], linestyles=(linestyle,))
#    ps = CC.collections[0].get_paths()[0].vertices
#    print(R)
#    for p in ps:
#        print(min_dist_to_target(np.array([0, r1, p[0], p[1]])))
#    plt.contour(CT, levels = [R], colors=(color,), linestyles=(linestyle,))

def get_constt(r1, R=Config.TAG_RANGE, color='k', linestyle='solid'):
    def get_constt(r1):
        k = 2
        x = np.linspace(-k*Config.TAG_RANGE, k*Config.TAG_RANGE)
        y = np.linspace(-k*Config.TAG_RANGE, k*Config.TAG_RANGE)
        X, Y = np.meshgrid(x, y)
        T = np.zeros(np.shape(X))
        for i, (xx, yy) in enumerate(zip(X, Y)):
            for j, (xxx, yyy) in enumerate(zip(xx, yy)):
        #        print(dominant_region(np.array([xxx, yyy])))
                T[i,j] = time_to_cap(np.array([0, r1, xxx, yyy]))
        
        return {'X': X, 'Y': Y, 'T': T}
    tctr = get_constt(r1)
    return tctr
    # CT = ax.contour(tctr['X'], tctr['Y'], tctr['T'], linestyles=(linestyle,))
    # plt.contour(CT, colors=(color,), linestyles=(linestyle,))
    # ax.clabel(CT, inline=True, fontsize=10)
#    ax.contour(tctr['X'], tctr['Y'], tctr['T'], colors=(color,), linestyles=(linestyle,))
    
    
# def plot_vcontour(r1, Rs, R=Config.TAG_RANGE):
#     fig, ax = plt.subplots()
#     plot_defender(ax, r1)
#     plot_target(ax, R=R, linestyle='solid')
# #    for R in Rs:
# #        plot_constd(ax, r1, R=R, linestyle='solid')
#     plot_constd(ax, r1)
# #    plot_constt(ax, r1)
    
#     k = 2
#     plt.ylim((-k*R, k*R))
#     plt.xlim((-k*R, k*R))
#     ax.axis('equal')
#     plt.grid()
#     plt.savefig('test.png')    

#def plot_vcontour(ss, Rs, colors=['b']*5, drs=[False, False]):
#    fig, ax = plt.subplots()
#    for s, R, color, dr in zip(ss, Rs, colors, drs):
#        plot_defender(ax, s[-1][2])
#        if dr:
#            plot_dr(ax, s, color=color)
#            plot_target(ax, R=R, color=color, linestyle='solid')
#        else:
#            plot_target(ax, R=R, color=color)
#        plot_contour(ax, s, color=color)
#    
#    plt.ylim((-2*R, 2*R))
#    plt.xlim((-2*R, 2*R))
#    ax.axis('equal')
#    plt.grid()
#    plt.savefig('test.png')