import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np
from math import pi, acos

from Config import Config
from vcontour import tangent
from plotter import plot_vcontour

#ss = []
#colors = ['g', 'c', 'b']
#for R in Rs:
#    ss.append(const_v(R))
#plot_vcontour(ss, Rs, colors=colors, drs=[True, False, False])

r = Config.CAP_RANGE
R = Config.TAG_RANGE
a = Config.VD/Config.VI
Rs = np.linspace(0.8*R, 1.8*R, 3)
plot_vcontour(1.5*R, [Rs[0]])
