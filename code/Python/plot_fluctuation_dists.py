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

import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from matplotlib.patches import Ellipse
import matplotlib.transforms as transforms

from sklearn.metrics import pairwise_distances
from skbio.stats.ordination import pcoa


gene_intersection = utils.get_gene_intersection()
phage_or_host_type = 'host'
subpop_type = 'revived_total'
transfers = numpy.asarray(utils.transfers)

measures = ['delta_f', 'log_ratio']

log_ratio_dict = {}

for seed_bank_type in utils.seed_bank_types:

    log_ratio_dict[seed_bank_type] = {}

    for phage_treatment_type in utils.phage_treatment_types:

        log_ratio_dict[seed_bank_type][phage_treatment_type] = {}
        log_ratio_dict[seed_bank_type][phage_treatment_type]['log_ratio'] = []
        log_ratio_dict[seed_bank_type][phage_treatment_type]['delta_f'] = []

        log_ratio_freqs_all = []
        delta_f_all = []

        exists_all = False
        for replicate in utils.replicates:
            exists, annotated_mapgd_dict = utils.load_annotated_breseq(phage_or_host_type, seed_bank_type, phage_treatment_type, subpop_type, replicate)

            if exists == False:
                continue

            else:
                exists_all = True

            replicate_name = '%s_%s_%d' % (seed_bank_type, phage_treatment_type, replicate)

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

                frequency_trajectory_idx = (frequency_trajectory>0) & (frequency_trajectory<1)
                transfers_i = transfers[frequency_trajectory_idx]
                frequency_trajectory = frequency_trajectory[frequency_trajectory_idx]

                delta_f_trajectory = numpy.absolute(frequency_trajectory[1:] - frequency_trajectory[:-1]) / (transfers_i[1:] - transfers_i[:-1])
                log_ratio_frequency_trajectory = numpy.log10(frequency_trajectory[1:]/frequency_trajectory[:-1])

                delta_f_all.extend(delta_f_trajectory.tolist())
                log_ratio_freqs_all.extend(log_ratio_frequency_trajectory.tolist())

        log_ratio_dict[seed_bank_type][phage_treatment_type]['delta_f'] = numpy.asarray(delta_f_all)
        log_ratio_dict[seed_bank_type][phage_treatment_type]['log_ratio'] = numpy.asarray(log_ratio_freqs_all)




fig = plt.figure(figsize = (4, 8)) #
fig.subplots_adjust(bottom= 0.15,  wspace=0.25)


import matplotlib as mpl
mpl.rcParams['xtick.labelsize'] = 7
mpl.rcParams['ytick.labelsize'] = 7

x_axis_labels = ['Absolute change in allele frequencies per-transfer, ' + r'$|\Delta f| / \Delta t$', 'Ratio of allele frequency changes, ' + r'$f(t+\delta t)/f(t)$']
y_axis_labels = ['Fraction ' + r'$\geq |\Delta f| / \Delta t$', 'Fraction ' + r'$\geq f(t+\delta t)/f(t)$']
for measure_idx, measure in enumerate(measures):

    ax = plt.subplot2grid((2, 1), (measure_idx, 0), colspan=1)

    for seed_bank_type in utils.seed_bank_types:
        for phage_treatment_type in utils.phage_treatment_types:

             log_ratio_array_sort = numpy.sort(log_ratio_dict[seed_bank_type][phage_treatment_type][measure])
             cdf = 1-numpy.arange(len(log_ratio_array_sort))/float(len(log_ratio_array_sort))

             label = utils.seed_bank_types_format_dict[seed_bank_type] + ', ' + utils.phage_treatment_types_format_dict[phage_treatment_type]

             ax.plot(10**log_ratio_array_sort, cdf, c =utils.color_dict[seed_bank_type], ls=utils.line_dict[phage_treatment_type], label=label, lw=1.5, alpha=0.8)


    #ax.set_xlabel('' , fontsize = 12)
    #ax.set_ylabel('PCo 2 (' + str(round(df_pcoa.proportion_explained[1]*100, 1)) + '%)' , fontsize = 12)

    if measure_idx == 0:
        ax.legend(loc="upper right", fontsize=6)

    if measure_idx == 1:
        ax.set_xlim(10**-2, 10**4)

    #ax.tick_params(labelsize=10)
    #ax.tick_params(axis='both', which='major', pad=1)

    #ax.yaxis.set_tick_params(labelsize=7)

    ax.set_xlabel(x_axis_labels[measure_idx], fontsize = 10)
    ax.set_ylabel(y_axis_labels[measure_idx], fontsize = 10)



    ax.set_xscale('log', basex=10)
    ax.set_yscale('log', basey=10)

fig_name = '%sfluctuation_survival.png' % (config.analysis_directory)
fig.savefig(fig_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
