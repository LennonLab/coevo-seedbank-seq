from __future__ import division
import numpy
import sys
from math import fabs
import glob, os, sys, re
import config
import utils
import parse_file

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

fig, ax = plt.subplots(figsize=(4,4))

count = 0

for phage_or_host_type in utils.phage_or_host_types:
    for seed_bank_type in  utils.seed_bank_types:
        for phage_treatment_type in utils.phage_treatment_types:
            for subpop_type in utils.subpop_types:
                for replicate in utils.replicates:

                    coverage_all = []

                    exists, annotated_mapgd_dict = utils.load_annotated_breseq(phage_or_host_type, seed_bank_type, phage_treatment_type, subpop_type, replicate)

                    if exists == False:
                        continue

                    for key, value in annotated_mapgd_dict.items():
                        coverage_all.append(value['coverage_trajectory'].tolist())

                    median_coverage_trajectory = []

                    coverage_all_zip = list(map(list, zip(*coverage_all)))
                    for coverage_all_zip_t in coverage_all_zip:

                        coverage_all_zip_t = numpy.asarray(coverage_all_zip_t)
                        coverage_all_zip_t = coverage_all_zip_t[~numpy.isnan(coverage_all_zip_t)]

                        median_coverage_trajectory.append(numpy.median(coverage_all_zip_t))

                    if phage_or_host_type == 'phage':
                        c = 'firebrick'
                        continue
                    else:
                        c = 'dodgerblue'

                    ax.plot(utils.transfers, median_coverage_trajectory, c=c, alpha=0.7)


                    count += 1


ax.set_xlabel('Transfer', fontsize=12)
ax.set_ylabel('Median coverage of\nsites that passed breseq', fontsize=12)
#ax.set_yscale('log', basey=10)

#custom_lines = [Line2D([0], [0], color='firebrick'),
#                Line2D([0], [0], color='dodgerblue')]
#ax.legend(custom_lines, ['Phage', 'Host'], loc="upper left")


fig_name = "%s%s" % (config.analysis_directory, 'coverage.png')
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()














                    #list(map(list, zip(*l)))
