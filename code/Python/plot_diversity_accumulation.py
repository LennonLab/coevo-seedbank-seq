from __future__ import division
import numpy
import sys
from math import fabs
import glob, os, sys, re
import config
import utils
import parse_file

import matplotlib.pyplot as plt




phage_or_host_type = 'host'
seed_bank_type = 'long_seed_bank'
phage_treatment_type = 'noPhage'
subpop_type = 'revived_spore'
replicate = 1


fig = plt.figure(figsize = (12, 4))

for subpop_type_idx, subpop_type in enumerate(utils.subpop_types):



    ax_i = plt.subplot2grid((1, len(utils.subpop_types)), (0, subpop_type_idx), colspan=1)

    for phage_or_host_type in utils.phage_or_host_types:
        for seed_bank_type in  utils.seed_bank_types:
            for phage_treatment_type in utils.phage_treatment_types:

                log_m_trajectory_all = []

                for replicate in utils.replicates:

                    m_trajectory = numpy.zeros(len(utils.transfers))
                    # check if it exists
                    exists, annotated_mapgd_dict = utils.load_annotated_mapgd(phage_or_host_type, seed_bank_type, phage_treatment_type, subpop_type, replicate)
                    if exists == False:
                        continue

                    for key, value in annotated_mapgd_dict.items():

                        frequency_trajectory = value['frequency_trajectory']

                        if sum(frequency_trajectory>0) < utils.min_n_non_zero_freqs:
                            continue

                        m_trajectory += frequency_trajectory

                    log_m_trajectory = numpy.log10(m_trajectory)
                    log_m_trajectory_all.append(log_m_trajectory)

                if len(log_m_trajectory_all) == 0:
                    continue

                mean_log_m_trajectory = numpy.mean(log_m_trajectory_all, axis=0)
                se_log_m_trajectory = numpy.std(log_m_trajectory_all, axis=0) / numpy.sqrt(len(log_m_trajectory_all))
                #print(subpop_type, se_log_m_trajectory)


                ax_i.plot(utils.transfers, 10**mean_log_m_trajectory, '.-', c='dodgerblue', alpha=0.7)

    ax_i.set_yscale('log', basey=10)
    ax_i.set_ylabel('Sum of derived allele frequenies, ' + r'$M(t)$', fontsize = 11)
    ax_i.set_xlabel('Transfer', fontsize=11)
    ax_i.set_title(utils.subpop_types[subpop_type], fontsize=12)





fig_name = "%s%s" % (config.analysis_directory, 'diversity_accumulation.png')
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
