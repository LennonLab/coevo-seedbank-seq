from __future__ import division
import numpy
import sys
from math import fabs
import glob, os, sys, re
import config
import utils
import parse_file

import matplotlib.pyplot as plt


#phage_or_host_type = 'host'
#seed_bank_type = 'long_seed_bank'
#phage_treatment_type = 'noPhage'
#replicate = 1

phage_or_host_type = 'phage'
seed_bank_type = 'long_seed_bank'
phage_treatment_type = 'noPhage'


for seed_bank_type in  utils.seed_bank_types:
    for phage_treatment_type in utils.phage_treatment_types:

        # set up plot

        n_rows = 0
        n_columns = 3

        for phage_or_host_type in utils.phage_or_host_types:
            if (phage_or_host_type == 'phage') and (phage_treatment_type == 'noPhage'):
                continue

            for subpop_type in utils.subpop_types:

                for replicate in utils.replicates:
                    exists, annotated_mapgd_dict = utils.load_annotated_mapgd(phage_or_host_type, seed_bank_type, phage_treatment_type, subpop_type, replicate)

                    if exists == False:
                        continue
                    if replicate == 1:
                        n_rows += 1
                    #if exists == False:
                    #    continue

        fig = plt.figure(figsize = (n_columns*2, n_rows*2))
        row_count = 0
        for phage_or_host_type in utils.phage_or_host_types:
            if (phage_or_host_type == 'phage') and (phage_treatment_type == 'noPhage'):
                continue

            for subpop_type_idx, subpop_type in enumerate(utils.subpop_types):

                column_count = 0

                for replicate in utils.replicates:

                    exists, annotated_mapgd_dict = utils.load_annotated_breseq(phage_or_host_type, seed_bank_type, phage_treatment_type, subpop_type, replicate)

                    if exists == False:
                        continue

                    ax_i = plt.subplot2grid((n_rows, n_columns), (row_count, column_count), colspan=1)

                    for key, value in annotated_mapgd_dict.items():
                        frequency_trajectory = value['frequency_trajectory']

                        if sum(frequency_trajectory>0) < utils.min_n_non_zero_freqs:
                            continue

                        # add zero at transfer 0 for aesthetic reasons
                        frequency_trajectory = numpy.insert(frequency_trajectory, 0, 0)

                        transfers = utils.transfers
                        transfers = numpy.insert(transfers, 0, 0)

                        # draw random color
                        rgb = utils.mut_freq_colormap()
                        rgb = utils.lighten_color(rgb, amount=0.5)

                        ax_i.plot(transfers, frequency_trajectory, '.-', c=rgb, alpha=0.3)

                    ax_i.set_xlim([0, max(utils.transfers)])
                    ax_i.set_ylim([0, 1])

                    ax_i.tick_params(axis="x", labelsize=6)
                    ax_i.tick_params(axis="y", labelsize=6)

                    if replicate == 1:
                        ax_i.set_ylabel(utils.subpop_types_format_dict[subpop_type], fontsize=10)

                    if subpop_type_idx == 0:
                        ax_i.set_title("Replicate %d" % replicate, fontsize=9)

                    column_count += 1


                if column_count == 0:
                    continue

                row_count += 1

        title = '%s, %s' % (utils.seed_bank_types_format_dict[seed_bank_type], utils.phage_treatment_types_format_dict[phage_treatment_type])

        if n_rows == 2:

            fig.text(0.5, 0.04, 'Transfer, ' + r'$t$', ha='center', va='center', fontsize=12)
            fig.text(0.03, 0.5, 'Allele frequency, ' + r'$f(t)$', ha='center', va='center', rotation='vertical',  fontsize=12)
            fig.suptitle(title, fontsize=14)

        else:
            fig.text(0.5, -0.01, 'Transfer, ' + r'$t$', ha='center', va='center', fontsize=10)
            fig.text(0.03, 0.5, 'Allele frequency, ' + r'$f(t)$', ha='center', va='center', rotation='vertical',  fontsize=10)
            fig.suptitle(title, y=1.06, fontsize=12)



        fig_name = "%smutation_trajectories_breseq/%s_%s_mutation_trajectories.png" % (config.analysis_directory, seed_bank_type, phage_treatment_type)
        fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
        plt.close()
