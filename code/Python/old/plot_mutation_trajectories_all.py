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
#subpop_type = 'revived_spore'
#replicate = 1

phage_or_host_type = 'phage'
seed_bank_type = 'long_seed_bank'
phage_treatment_type = 'noPhage'



fig = plt.figure(figsize = (20, 10)) #
fig.subplots_adjust(bottom= 0.15,  wspace=0.25)


for seed_bank_type_idx, seed_bank_type in enumerate(utils.seed_bank_types):
    for phage_treatment_type_idx, phage_treatment_type in enumerate(utils.phage_treatment_types):

        # set up plot
        for replicate in utils.replicates:

            replicate_idx = (replicate-1)+(seed_bank_type_idx * len(utils.replicates))

            ax_host = plt.subplot2grid((3, 6), (phage_treatment_type_idx*2, replicate_idx))

            if phage_treatment_type == 'SPO1':
                ax_phage = plt.subplot2grid((3, 6), (1, replicate_idx))
                exists_phage, annotated_mapgd_dict_phage = utils.load_annotated_breseq('phage', seed_bank_type, phage_treatment_type, 'filtered_phage', replicate)

            exists_host, annotated_mapgd_dict_host = utils.load_annotated_breseq('host', seed_bank_type, phage_treatment_type, 'revived_total', replicate)

            for key, value in annotated_mapgd_dict_host.items():
                frequency_trajectory = value['frequency_trajectory']

                if sum(frequency_trajectory>0) < utils.min_n_non_zero_freqs:
                    continue

                # add zero at transfer 0 for aesthetic reasons
                frequency_trajectory = numpy.insert(frequency_trajectory, 0, 0)

                transfers = utils.transfers
                transfers = numpy.insert(transfers, 0, 0)

                # draw random color
                #rgb = utils.mut_freq_colormap()
                #rgb = utils.lighten_color(rgb, amount=0.5)

                #ax_host.plot(transfers, frequency_trajectory, '.-', c=rgb, alpha=0.3)

            ax_host.set_xlim([0, max(utils.transfers)])
            ax_host.set_ylim([0, 1])

            ax_host.tick_params(axis="x", labelsize=6)
            ax_host.tick_params(axis="y", labelsize=6)

            #if replicate == 1:
            #    ax_i.set_ylabel(utils.subpop_types[subpop_type], fontsize=10)

            if phage_treatment_type_idx == 0:
                ax_host.set_title("Replicate %d" % replicate, fontsize=10)



        #title = '%s, %s' % (utils.seed_bank_types_format_dict[seed_bank_type], utils.phage_treatment_types_format_dict[phage_treatment_type])
        #fig.suptitle(title, fontsize=14)


fig.text(0.5, 0.04, 'Transfer', ha='center', va='center', fontsize=12)
fig.text(0.03, 0.5, 'Allele frequency, ' + r'$f(t)$', ha='center', va='center', rotation='vertical',  fontsize=12)


fig_name = "%smutation_trajectories.png" % config.analysis_directory
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
