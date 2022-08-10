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

iter = 10000


ks_path = "%sks_rho_dist.pickle" % (config.data_directory)

seedbank_pairs = [('no_seed_bank', 'long_seed_bank')]


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
        rho_null_all = []

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

                    if (sum(frequency_trajectory_host_i_to_keep==1) == len(frequency_trajectory_host_i_to_keep)) or (sum(frequency_trajectory_phage_i_to_keep==1) == len(frequency_trajectory_phage_i_to_keep)):
                        continue

                    rho = numpy.corrcoef(frequency_trajectory_host_i_to_keep, frequency_trajectory_phage_i_to_keep)[1,0]

                    if numpy.isnan(rho) == True:
                        continue

                    rho_all.append(rho)

                    for i in range(100):

                        numpy.random.shuffle(frequency_trajectory_host_i_to_keep)
                        numpy.random.shuffle(frequency_trajectory_phage_i_to_keep)
                        rho_null = numpy.corrcoef(frequency_trajectory_host_i_to_keep, frequency_trajectory_phage_i_to_keep)[1,0]

                        if numpy.isnan(rho_null) == True:
                            continue

                        rho_null_all.append(rho_null)


        frequency_trajectory_dict[seed_bank_type] = {}
        frequency_trajectory_dict[seed_bank_type]['rho'] = numpy.asarray(rho_all)
        frequency_trajectory_dict[seed_bank_type]['rho_null'] = numpy.asarray(rho_null_all)

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

    # same versus null distribution
    for seed_bank_type in utils.seed_bank_types:

        rho_i = frequency_trajectory_dict[seed_bank_type]['rho']
        rho_j = frequency_trajectory_dict[seed_bank_type]['rho_null']

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

        ks_dist_dict[seed_bank_type] = {}
        ks_dist_dict[seed_bank_type]['D'] = D
        ks_dist_dict[seed_bank_type]['p'] = p_perm
        ks_dist_dict[seed_bank_type]['null_mean'] = numpy.mean(D_null_all)
        ks_dist_dict[seed_bank_type]['null_median'] = numpy.median(D_null_all)
        ks_dist_dict[seed_bank_type]['lower_ci'] = lower_ci
        ks_dist_dict[seed_bank_type]['upper_ci'] = upper_ci

        print(seed_bank_type, D, p_perm)


    with open(ks_path, 'wb') as handle:
        pickle.dump(ks_dist_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)


def load_ks_dict():

    with open(ks_path, 'rb') as handle:
        ks_dist_dict = pickle.load(handle)
    return ks_dist_dict





def plot_distribution():

    frequency_trajectory_dict = get_frequency_trajectory_dict()
    ks_dist_dict = load_ks_dict()

    #fig = plt.figure(figsize = (5, 8)) #
    #fig.subplots_adjust(bottom= 0.15,  wspace=0.25)

    #ax = plt.subplot2grid((1, 0), (0, 0), colspan=1)
    fig, ax = plt.subplots(figsize=(5, 4))

    import matplotlib as mpl
    mpl.rcParams['xtick.labelsize'] = 7
    mpl.rcParams['ytick.labelsize'] = 7

    for seed_bank_type_idx, seed_bank_type in enumerate(utils.seed_bank_types):

        rho_all_sort = numpy.sort(frequency_trajectory_dict[seed_bank_type]['rho'])
        cdf = 1-numpy.arange(len(rho_all_sort))/float(len(rho_all_sort))

        label = utils.seed_bank_types_format_dict[seed_bank_type] + ', ' + utils.phage_treatment_types_format_dict[phage_treatment_type]
        #ax.plot(rho_all_sort, cdf, c=utils.color_dict[seed_bank_type], ls='-', label=label, lw=2, alpha=1)
        #ax.axvline(x=numpy.median(rho_all_sort), color=utils.color_dict[seed_bank_type], linestyle=':', alpha = 0.8, zorder=1)

        density = gaussian_kde(rho_all_sort)
        xs = numpy.linspace(-1,1,1000)
        density.covariance_factor = lambda : .25
        density._compute_covariance()

        ax.plot(xs, density(xs), ls='-', lw=2, c=utils.color_dict[seed_bank_type], label=label, alpha=1)


    ax.set_xlim(-1, 1)
    #ax.set_yscale('log', basey=10)
    ax.set_xlabel('Correlation coefficient between\nhost and phage mutation trajectories', fontsize = 10)
    ax.set_ylabel('Probability density', fontsize = 12)
    ax.legend(loc="upper right", fontsize=6)

    # plot inset
    #ins_ks = inset_axes(ax, width="100%", height="100%", loc='lower right', bbox_to_anchor=(0.14,0.07,0.4,0.38), bbox_transform=ax.transAxes)
    ins_ks = inset_axes(ax, width="100%", height="100%", loc='lower right', bbox_to_anchor=(0.16,0.05,0.4,0.38), bbox_transform=ax.transAxes)

    for seedbank_pair_idx, seedbank_pair in enumerate(seedbank_pairs):

        seedbank_pair_D = ks_dist_dict[seedbank_pair]['D']
        p_perm = ks_dist_dict[seedbank_pair]['p']

        marker_style = dict(color='k', marker='o',
            markerfacecoloralt=utils.color_dict[seedbank_pair[1]],
            markerfacecolor=utils.color_dict[seedbank_pair[0]],
            mew=0.5)

        ins_ks.plot(seedbank_pair_idx, seedbank_pair_D, markersize = 11,   \
                linewidth=0.4,  alpha=1, zorder=3, fillstyle='left', **marker_style)

        #ins_ks.plot(seed_bank_type_idx, D, markersize = 11, marker = 'o',  \
        #    linewidth=0.4,  alpha=1, color=utils.color_dict[seed_bank_type], zorder=2)

        if p_perm < 0.05:
            ins_ks.text(seedbank_pair_idx, seedbank_pair_D+0.012, '*', ha='center', fontsize=9)

        delta = 0.05
        fill_between_x = numpy.asarray([seedbank_pair_idx-delta, seedbank_pair_idx+delta])
        fill_between_y_lower = numpy.asarray([ks_dist_dict[seedbank_pair]['lower_ci'], ks_dist_dict[seedbank_pair]['lower_ci']])
        fill_between_y_upper = numpy.asarray([ks_dist_dict[seedbank_pair]['upper_ci'], ks_dist_dict[seedbank_pair]['upper_ci']])

        ins_ks.fill_between(fill_between_x, y1=fill_between_y_lower, y2=fill_between_y_upper, color='grey')


    ins_ks.tick_params(labelsize=4)
    ins_ks.tick_params(axis='both', which='major', pad=1)
    ins_ks.set_xlim([-0.5, 2.5])
    ins_ks.set_ylim([-0.005, 0.18])
    ins_ks.set_xticks([0, 1, 2])

    ins_ks.set_xticklabels(['No vs. short\nseed bank', 'No vs. long\nseed bank', 'Short vs. long\nseed bank'] )

    ins_ks.set_ylabel("KS distance", fontsize=6)
    ins_ks.tick_params(axis='x', labelsize=5, length = 0)

    ins_ks.axhline(y=0, color='k', linestyle=':', alpha = 0.8, zorder=1)
    ins_ks.yaxis.set_tick_params(labelsize=4)

    fig_name = '%srho_pdf.png' % (config.analysis_directory)
    fig.savefig(fig_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()



def plot_distribution_null():

    frequency_trajectory_dict = get_frequency_trajectory_dict()
    ks_dist_dict = load_ks_dict()

    fig = plt.figure(figsize = (12, 4)) #
    fig.subplots_adjust(bottom= 0.15,  wspace=0.25)

    #x_axis_labels = ['Maximum observed allele frequency, ' + r'$f_{max}$', 'Absolute change in allele frequencies per-transfer, ' + r'$|\Delta f| / \Delta t$', 'Ratio of allele frequency changes, ' + r'$f(t+\delta t)/f(t)$']
    #y_axis_labels = ['Fraction ' + r'$\geq f_{max}$', 'Fraction ' + r'$\geq |\Delta f| / \Delta t$', 'Fraction ' + r'$\geq f(t+\delta t)/f(t)$']
    #for measure_idx, measure in enumerate(measures):

    seed_bank_types = ['no_seed_bank', 'long_seed_bank']

    ax_concept = plt.subplot2grid((1, 3), (0, 0), colspan=1)

    phage = ax_concept.twinx()

    ax_concept.set_xlabel("Transfer", fontsize = 12)
    ax_concept.set_ylabel("Allele frequency, host", fontsize = 12)
    phage.set_ylabel("Allele frequency, phage", fontsize = 12)

    t = numpy.asarray([1, 4, 7, 10, 14])
    host_pos = numpy.asarray([0, 0, 0.25, 0.38, 0.68])
    phage_pos = numpy.asarray([0, 0.2, 0.32, 0.4, 0.82])

    host_neg = numpy.asarray([0, 0.22, 0.15, 0.08, 0.01])
    phage_neg = numpy.asarray([0, 0.1, 0.3, 0.32, 0.54])


    # dodgerblue, orangered
    ax_host, = ax_concept.plot(t, host_pos, color=utils.color_dict['no_seed_bank'], ls='-')
    ax_phage, = ax_concept.plot(t, phage_pos, color=utils.color_dict['long_seed_bank'], ls='-')

    ax_concept.plot(t, host_neg, color=utils.color_dict['no_seed_bank'], ls=':')
    ax_concept.plot(t, phage_neg, color=utils.color_dict['long_seed_bank'], ls=':')

    ax_concept.yaxis.label.set_color(ax_host.get_color())
    phage.yaxis.label.set_color(ax_phage.get_color())

    phage.set_yticklabels([])
    phage.set_yticks([])
    ax_concept.set_yticklabels([])
    ax_concept.set_yticks([])

    ax_concept.set_xticklabels(['1', '4', '7', '10', '14'])
    ax_concept.set_xticks(t)

    ax_concept.set_xlim([1,14])
    ax_concept.set_ylim([0,1])


    from matplotlib.lines import Line2D
    custom_lines = [Line2D([0], [0], color='k', lw=2, ls='-'),
                    Line2D([0], [0], color='k', lw=2, ls='--')]

    ax_concept.legend(custom_lines, [r'$\rho >0$', r'$\rho < 0$'], loc ='upper left')
    ax_concept.text(-0.1, 1.08, utils.subplot_labels[0], fontsize=11, fontweight='bold', ha='center', va='center', transform=ax_concept.transAxes)


    for seed_bank_type_idx, seed_bank_type in enumerate(seed_bank_types):

        ax = plt.subplot2grid((1, 3), (0, seed_bank_type_idx+1), colspan=1)

        rho_all_sort = numpy.sort(frequency_trajectory_dict[seed_bank_type]['rho'])
        rho_null_all_sort = numpy.sort(frequency_trajectory_dict[seed_bank_type]['rho_null'])

        xs = numpy.linspace(-1,1,1000)

        density = gaussian_kde(rho_all_sort)
        density.covariance_factor = lambda : .25
        density._compute_covariance()

        density_null = gaussian_kde(rho_null_all_sort)
        density_null.covariance_factor = lambda : .25
        density_null._compute_covariance()


        ax.plot(xs, density(xs), ls='-', lw=2, c=utils.color_dict[seed_bank_type], label='Observed', alpha=1)
        ax.plot(xs, density_null(xs), ls='-', lw=2, c='lightgray', label='Null', alpha=1)

        ax.fill_between(xs, density(xs), color=utils.color_dict[seed_bank_type], alpha=0.7)
        ax.fill_between(xs, density_null(xs), color='lightgray', alpha=0.4)

        ax.set_xlim(-1, 1)
        ax.set_ylim(0, 0.8)
        #ax.set_yscale('log', basey=10)
        ax.set_xlabel('Correlation coefficient between\nhost and phage mutation trajectories, ' + r'$\rho$', fontsize = 11)
        ax.set_ylabel('Probability density', fontsize = 12)
        ax.legend(loc="upper right", fontsize=6)
        ax.set_title(utils.seed_bank_types_format_dict[seed_bank_type], fontsize=12, color='k')

        ax.text(0.15,  0.95, r'$D =$' + str(round(ks_dist_dict[seed_bank_type]['D'], 3)), fontsize=8, color='k', ha='center', va='center', transform=ax.transAxes  )

        p_value = ks_dist_dict[seed_bank_type]['p']

        if p_value == 0:
            ax.text(0.18,  0.89, r'$P <$' + str(1/iter), fontsize=8, color='k', ha='center', va='center', transform=ax.transAxes  )
        else:
            ax.text(0.15,  0.89, r'$P =$' + str(round(p_value, 3)), fontsize=8, color='k', ha='center', va='center', transform=ax.transAxes  )


        ax.text(-0.1, 1.08, utils.subplot_labels[seed_bank_type_idx+1], fontsize=11, fontweight='bold', ha='center', va='center', transform=ax.transAxes)

    fig.subplots_adjust(wspace=0.34) #hspace=0.3, wspace=0.5
    fig_name = '%srho_null_pdf.png' % (config.analysis_directory)
    fig.savefig(fig_name, bbox_inches = "tight", pad_inches = 0.3, dpi = 600)
    plt.close()




#make_ks_dist_rho_dict()
plot_distribution_null()
