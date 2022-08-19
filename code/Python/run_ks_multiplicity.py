from __future__ import division
import numpy
import sys
from math import fabs
import glob, os, sys, re
import config
import utils
import parse_file
import pickle

import scipy.stats as stats



mult_dict = {}
mult_dict['phage'] = {}
mult_dict['no_phage'] = {}

mult_dict['phage']['seedbank'] = []
mult_dict['phage']['no_seedbank'] = []
mult_dict['no_phage']['seedbank'] = []
mult_dict['no_phage']['no_seedbank'] = []
for line in open("%smult_host.csv" % config.data_directory, 'r'):

    if 'A8O17_RS00100' in line:
        continue

    line = line.strip().split(',')
    mult = [float(m) for m in line[1:]]
    population = line[0]
    if 'noPhage_' in population:

        if 'no_seed_bank_' in population:
            mult_dict['no_phage']['no_seedbank'].extend(mult)
        else:
            mult_dict['no_phage']['seedbank'].extend(mult)

    else:
        if 'no_seed_bank_' in population:
            mult_dict['phage']['no_seedbank'].extend(mult)
        else:
            mult_dict['phage']['seedbank'].extend(mult)



iter = 10000
for phage_status in ['phage', 'no_phage']:

    seedbank = numpy.asarray(mult_dict[phage_status]['seedbank'])
    no_seedbank = numpy.asarray(mult_dict[phage_status]['no_seedbank'])

    D, p = stats.ks_2samp(seedbank, no_seedbank)

    measure_array_merged = numpy.concatenate((seedbank, no_seedbank))

    D_null_all = []
    for i in range(iter):
        numpy.random.shuffle(measure_array_merged)

        seedbank_null = measure_array_merged[:len(seedbank)]
        no_seedbank_null = measure_array_merged[len(no_seedbank):]

        D_null, p_null = stats.ks_2samp(seedbank_null, no_seedbank_null)
        D_null_all.append(D_null)

    D_null_all = numpy.asarray(D_null_all)
    D_null_all = numpy.sort(D_null_all)

    p_perm = sum(D_null_all > D)/iter

    print(phage_status, D, p_perm)



# same thing for phage


mult_dict = {}
seedbank_mult = []
no_seedbank_mult = []
for line in open("%smult_phage.csv" % config.data_directory, 'r'):

    if 'SPO1_8' in line:
        continue

    line = line.strip().split(',')
    mult = [float(m) for m in line[1:]]
    population = line[0]

    if 'no_seed_bank' in population:
        no_seedbank_mult.extend(mult)
    else:
        seedbank_mult.extend(mult)




seedbank = numpy.asarray(seedbank_mult)
no_seedbank = numpy.asarray(no_seedbank_mult)

D, p = stats.ks_2samp(seedbank, no_seedbank)

measure_array_merged = numpy.concatenate((seedbank, no_seedbank))

D_null_all = []
for i in range(iter):
    numpy.random.shuffle(measure_array_merged)

    seedbank_null = measure_array_merged[:len(seedbank)]
    no_seedbank_null = measure_array_merged[len(no_seedbank):]

    D_null, p_null = stats.ks_2samp(seedbank_null, no_seedbank_null)
    D_null_all.append(D_null)

D_null_all = numpy.asarray(D_null_all)
D_null_all = numpy.sort(D_null_all)

p_perm = sum(D_null_all > D)/iter

print('Phage multiplicity test', D, p_perm)



#print(mult_dict)
