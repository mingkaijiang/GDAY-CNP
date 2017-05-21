#!/usr/bin/env python

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import plot_settings as ps
import pandas as pd
import datetime as dt

# load data
data = np.loadtxt("y")


plt.rcParams['figure.subplot.hspace'] = 0.15
plt.rcParams['figure.subplot.wspace'] = 0.15
plt.rcParams['font.size'] = 12
plt.rcParams['legend.fontsize'] = 12
plt.rcParams['xtick.labelsize'] = 12.0
plt.rcParams['ytick.labelsize'] = 12.0
plt.rcParams['axes.labelsize'] = 12.0
plt.rcParams['font.family'] = "sans-serif"
plt.rcParams['font.style'] = "normal"
plt.rcParams['font.serif'] = "Helvetica"
plt.rcParams['text.usetex'] = False


fig = plt.figure()

ax1 = fig.add_subplot(311)

ax1.plot(1.0/data[:,0], "b-")

ax1.set_ylabel("leaf N:C")
ax1.set_ylim(0., 0.05)
plt.setp(ax1.get_xticklabels(), visible=False)

ax2 = fig.add_subplot(312)
ax2.plot(data[:,1], "b-")
ax2.set_ylabel("CO$_2$ (ppm)")
plt.setp(ax2.get_xticklabels(), visible=False)

ax3 = fig.add_subplot(313)
ax3.plot(data[:,2]*365.25, "b-")
ax3.set_ylabel("N dep (t ha$^{-1}$ yr$^{-1}$)")
ax3.set_xlabel("DOY since 1760")
#ax3.legend(numpoints=1, loc="upper left")
plt.show()
