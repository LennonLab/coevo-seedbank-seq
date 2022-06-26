from __future__ import division
import numpy
import sys
from math import fabs
import glob, os, sys, re
import config
import utils
import parse_file
import pickle
import pandas
import scipy

import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from matplotlib.patches import Ellipse
import matplotlib.transforms as transforms
from mpl_toolkits.axes_grid1.inset_locator import inset_axes



fig, ax = plt.subplots(figsize=(4,4))

phage = ax.twinx()

ax.set_xlabel("Transfer")
ax.set_ylabel("Allele frequency, host")
phage.set_ylabel("Allele frequency, phage")


t = numpy.asarray([1, 4, 7, 10, 14])
host_pos = numpy.asarray([0, 0, 0.25, 0.38, 0.68])
phage_pos = numpy.asarray([0, 0.2, 0.32, 0.4, 0.82])

host_neg = numpy.asarray([0, 0.22, 0.15, 0.08, 0.01])
phage_neg = numpy.asarray([0, 0.1, 0.3, 0.32, 0.54])



ax_host, = ax.plot(t, host_pos, color='dodgerblue', ls='-')
ax_phage, = ax.plot(t, phage_pos, color='orangered', ls='-')

ax.plot(t, host_neg, color='dodgerblue', ls=':')
ax.plot(t, phage_neg, color='orangered', ls=':')




#lns = [ax_host, ax_phage]
#ax.legend(handles=lns, loc='upper left')



ax.yaxis.label.set_color(ax_host.get_color())
phage.yaxis.label.set_color(ax_phage.get_color())

phage.set_yticklabels([])
phage.set_yticks([])
ax.set_yticklabels([])
ax.set_yticks([])

ax.set_xticklabels(['1', '4', '7', '10', '14'])
ax.set_xticks(t)

ax.set_xlim([1,14])
ax.set_ylim([0,1])



from matplotlib.lines import Line2D
custom_lines = [Line2D([0], [0], color='k', lw=2, ls='-'),
                Line2D([0], [0], color='k', lw=2, ls='--')]

ax.legend(custom_lines, [r'$\rho >0$', r'$\rho < 0$'], loc ='upper left')


fig_name = '%sconceptual_fig.png' % (config.analysis_directory)
fig.savefig(fig_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
