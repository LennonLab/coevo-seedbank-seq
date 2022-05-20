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

from sklearn.metrics import pairwise_distances
from skbio.stats.ordination import pcoa
import scipy.stats as stats

numpy.random.seed(123456789)

iter =10000



gene_intersection = utils.get_gene_intersection()
phage_treatment_type = 'SPO1'


def filter_trajectoties(phage_or_host_type, annotated_mapgd_dict):

    frequency_trajectory_all = []

    for position, position_dict in annotated_mapgd_dict.items():

        gene = position_dict['gene']

        if (gene == 'intergenic') or (gene == 'repeat'):
            continue

        if phage_or_host_type == 'host':
            # ignore genes that aren't in both strains
            if gene not in gene_intersection:
                continue

        if (position_dict['annotation'] == 'synonymous') or (position_dict['annotation'] == 'noncoding'):
            continue

        frequency_trajectory = position_dict['frequency_trajectory']

        if sum(frequency_trajectory>0) < utils.min_n_non_zero_freqs:
            continue

        frequency_trajectory_all.append(frequency_trajectory)

    return frequency_trajectory_all



def get_frequency_trajectory_dict():

    frequency_trajectory_dict = {}

    for seed_bank_type in utils.seed_bank_types:

        rho_all = []

        for replicate in utils.replicates:
            exists_host, annotated_mapgd_dict_host = utils.load_annotated_breseq('host', seed_bank_type, phage_treatment_type, 'revived_total', replicate)
            exists_phage, annotated_mapgd_dict_phage = utils.load_annotated_breseq('phage', seed_bank_type, phage_treatment_type, 'filtered_phage', replicate)

            if (exists_host == False) or (exists_phage == False):
                continue

            frequency_trajectory_host = filter_trajectoties('host', annotated_mapgd_dict_host)
            frequency_trajectory_phage = filter_trajectoties('phage', annotated_mapgd_dict_phage)


            for frequency_trajectory_host_i_idx, frequency_trajectory_host_i in enumerate(frequency_trajectory_host):

                for frequency_trajectory_phage_i_idx, frequency_trajectory_phage_i in enumerate(frequency_trajectory_phage[frequency_trajectory_host_i_idx:]):

                    idx_to_keep = (frequency_trajectory_host_i>0) & (frequency_trajectory_phage_i>0)

                    frequency_trajectory_host_i_to_keep = frequency_trajectory_host_i[idx_to_keep]
                    frequency_trajectory_phage_i_to_keep = frequency_trajectory_phage_i[idx_to_keep]

                    if len(frequency_trajectory_host_i_to_keep) < 3:
                        continue

                    rho = numpy.corrcoef(frequency_trajectory_host_i_to_keep, frequency_trajectory_phage_i_to_keep)[1,0]
                    if numpy.isnan(rho) == True:
                        continue

                    rho_all.append(rho)

        frequency_trajectory_dict[seed_bank_type] = {}
        frequency_trajectory_dict[seed_bank_type]['rho'] = numpy.asarray(rho_all)

    return frequency_trajectory_dict



frequency_trajectory_dict = get_frequency_trajectory_dict()

#fig = plt.figure(figsize = (5, 8)) #
#fig.subplots_adjust(bottom= 0.15,  wspace=0.25)

#ax = plt.subplot2grid((1, 0), (0, 0), colspan=1)
fig, ax = plt.subplots(figsize=(4, 4))
for seed_bank_type_idx, seed_bank_type in enumerate(utils.seed_bank_types):

    rho_all_sort = numpy.sort(frequency_trajectory_dict[seed_bank_type]['rho'])
    cdf = 1-numpy.arange(len(rho_all_sort))/float(len(rho_all_sort))

    label = utils.seed_bank_types_format_dict[seed_bank_type] + ', ' + utils.phage_treatment_types_format_dict[phage_treatment_type]
    ax.plot(rho_all_sort, cdf, c=utils.color_dict[seed_bank_type], ls='-', label=label, lw=2, alpha=1)
    ax.axvline(x=numpy.median(rho_all_sort), color=utils.color_dict[seed_bank_type], linestyle=':', alpha = 0.8, zorder=1)

ax.set_xlim(-1, 1)
#ax.set_yscale('log', basey=10)
ax.set_xlabel('Correlation coefficient between host and phage trajectories', fontsize = 10)
ax.set_ylabel('Fraction of coefficients >= ', fontsize = 8)
ax.legend(loc="upper right", fontsize=6)

fig_name = '%srho.png' % (config.analysis_directory)
fig.savefig(fig_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
