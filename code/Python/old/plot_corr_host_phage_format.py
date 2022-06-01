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
from scipy.stats import gaussian_kde


numpy.random.seed(123456789)

iter =10000


ks_path = "%sks_rho_dist.pickle" % (config.data_directory)



seedbank_pairs = utils.seedbank_pairs

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




def make_ks_dist_rho_dict():

    frequency_trajectory_dict = get_frequency_trajectory_dict()

    ks_dist_dict = {}

    for seedbank_pair in seedbank_pairs:

        seed_bank_type_i = seedbank_pair[0]
        seed_bank_type_j = seedbank_pair[1]

        rho_i = frequency_trajectory_dict[seed_bank_type_i]['rho']
        rho_j = frequency_trajectory_dict[seed_bank_type_j]['rho']

        D, p = scipy.stats.ks_2samp(rho_i, rho_j)

        measure_array_merged = numpy.concatenate((rho_i, rho_j))

        D_null_all = []
        for i in range(iter):
            numpy.random.shuffle(measure_array_merged)

            null_i = measure_array_merged[:len(rho_i)]
            null_j = measure_array_merged[len(rho_i):]

            D_null, p_null = scipy.stats.ks_2samp(null_i, null_j)
            D_null_all.append(D_null)

        D_null_all = numpy.asarray(D_null_all)
        D_null_all = numpy.sort(D_null_all)
        p_perm = sum(D_null_all > D)/iter
        lower_ci = D_null_all[int(iter*0.025)]
        upper_ci = D_null_all[int(iter*0.975)]

        ks_dist_dict[seedbank_pair] = {}
        ks_dist_dict[seedbank_pair]['D'] = D
        ks_dist_dict[seedbank_pair]['p'] = p_perm
        ks_dist_dict[seedbank_pair]['null_mean'] = numpy.mean(D_null_all)
        ks_dist_dict[seedbank_pair]['null_median'] = numpy.median(D_null_all)
        ks_dist_dict[seedbank_pair]['lower_ci'] = lower_ci
        ks_dist_dict[seedbank_pair]['upper_ci'] = upper_ci

        print(seed_bank_type_i, seed_bank_type_j, D, p_perm)

    with open(ks_path, 'wb') as handle:
        pickle.dump(ks_dist_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)


def load_ks_dict():

    with open(ks_path, 'rb') as handle:
        ks_dist_dict = pickle.load(handle)
    return ks_dist_dict


def jitter_dots(dots):
    offsets = dots.get_offsets()
    jittered_offsets = offsets
    # only jitter in the x-direction
    jittered_offsets[:, 0] += numpy.random.uniform(-0.2, 0.2, offsets.shape[0])
    dots.set_offsets(jittered_offsets)


frequency_trajectory_dict = get_frequency_trajectory_dict()
ks_dist_dict = load_ks_dict()

fig = plt.figure(figsize = (5, 12)) #
fig.subplots_adjust(bottom= 0.15,  wspace=0.25)

ax_survial = plt.subplot2grid((3, 1), (0, 0))
ax_pdf = plt.subplot2grid((3, 1), (1, 0))
ax_jitter = plt.subplot2grid((3, 1), (2, 0))
#fig, ax = plt.subplots(figsize=(5, 4))

import matplotlib as mpl
mpl.rcParams['xtick.labelsize'] = 7
mpl.rcParams['ytick.labelsize'] = 7

for seed_bank_type_idx, seed_bank_type in enumerate(utils.seed_bank_types):

    rho_all_sort = numpy.sort(frequency_trajectory_dict[seed_bank_type]['rho'])
    cdf = 1-numpy.arange(len(rho_all_sort))/float(len(rho_all_sort))

    label = utils.seed_bank_types_format_dict[seed_bank_type] + ', ' + utils.phage_treatment_types_format_dict[phage_treatment_type]
    ax_survial.plot(rho_all_sort, cdf, c=utils.color_dict[seed_bank_type], ls='-', label=label, lw=2, alpha=1)
    #ax.axvline(x=numpy.median(rho_all_sort), color=utils.color_dict[seed_bank_type], linestyle=':', alpha = 0.8, zorder=1)

    #hist, bins = numpy.histogram(rho_all_sort, bins=20)
    #logbins = numpy.logspace(numpy.log10(bins[0]), numpy.log10(bins[-1]),len(bins))
    #ax_pdf.hist(rho_all_sort, bins=12, color=utils.color_dict[seed_bank_type], alpha=1, density=True, histtype='step')

    density = gaussian_kde(rho_all_sort)
    xs = numpy.linspace(-1,1,1000)
    density.covariance_factor = lambda : .25
    density._compute_covariance()

    ax_pdf.plot(xs, density(xs), ls='-', lw=2, c=utils.color_dict[seed_bank_type])


    # cmap=utils.seed_bank_type_cmap_dict[seed_bank_type]
    dots = ax_jitter.scatter(numpy.full((len(rho_all_sort), 1), seed_bank_type_idx), rho_all_sort, s=50, alpha=0.12, c=utils.color_dict[seed_bank_type])
    jitter_dots(dots)


ax_survial.set_xlim(-1, 1)
ax_survial.set_xlabel('Corr. coefficient between host and phage trajectories, ' +  r'$\rho$', fontsize = 10)
ax_survial.set_ylabel('Fraction ' + r'$\geq \rho$', fontsize = 10)

ax_pdf.set_xlim(-1, 1)
ax_pdf.set_xlabel('Corr. coefficient between host and phage trajectories, ' +  r'$\rho$', fontsize = 10)
ax_pdf.set_ylabel('Probability density', fontsize = 10)


ax_jitter.set_xlim(-0.5, 2.5)
ax_jitter.set_ylabel('Corr. coefficient between host and phage trajectories, ' +  r'$\rho$', fontsize = 8)

ax_jitter.set_xticks([0, 1, 2])
ax_jitter.set_xticklabels(['No seed bank', 'Short seed bank', 'Long seed bank'] )



for seedbank_pair_idx, seedbank_pair in enumerate(seedbank_pairs):

    seedbank_pair_D = ks_dist_dict[seedbank_pair]['D']

    marker_style = dict(color='k', marker='o',
        markerfacecoloralt=utils.color_dict[seedbank_pair[1]],
        markerfacecolor=utils.color_dict[seedbank_pair[0]],
        mew=0.5)


    delta = 0.05
    fill_between_x = numpy.asarray([seedbank_pair_idx-delta, seedbank_pair_idx+delta])
    fill_between_y_lower = numpy.asarray([ks_dist_dict[seedbank_pair]['lower_ci'], ks_dist_dict[seedbank_pair]['lower_ci']])
    fill_between_y_upper = numpy.asarray([ks_dist_dict[seedbank_pair]['upper_ci'], ks_dist_dict[seedbank_pair]['upper_ci']])




fig_name = '%srho_test.png' % (config.analysis_directory)
fig.savefig(fig_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
