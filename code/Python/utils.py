from __future__ import division
import numpy
import sys
from math import fabs
import glob, os, sys, re
import config

from Bio.Alphabet import generic_dna
from Bio import SeqIO
from BCBio import GFF


phage_or_host_types = ['phage', 'host']
seed_bank_types = ['short_seed_bank', 'long_seed_bank', 'no_seed_bank']
phage_treatment_types = ['noPhage', 'SPO1']
subpop_types = ['filtered_phage', 'revived_total', 'revived_spore']
references = ['delta6-ANC', 'dspoIIE-ANC', 'SPO1-ANC']
replicates = [1, 2, 3]


transfers = [1,4,7,10,14]


#def gff3_to_fasta():
#    for reference in references:
#        with open("%s%s.gff3" % (config.data_directory, reference), "rU") as input_handle:
#            with open("%sGCF_001660525.1_ASM166052v1_genomic.fa" % config.data_directory , "w") as output_handle:
#                sequences = SeqIO.parse(input_handle, "genbank")
#                count = SeqIO.write(sequences, output_handle, "fasta")


def gff3_to_genbank():

    for reference in references:

        gff_path = "%s%s_no_fasta.gff3" % (config.data_directory, reference)
        fasta_path = "%s%s.fa" % (config.data_directory, reference)
        gbk_path = "%s%s.gbk" % (config.data_directory, reference)

        fasta_input = SeqIO.to_dict(SeqIO.parse(fasta_path, "fasta", generic_dna))
        gff_iter = GFF.parse(gff_path, fasta_input)

        SeqIO.write(gff_iter, gbk_path, "genbank")






def get_sample_metadata_dict():

    metadata_path = '%s%s' % (config.data_directory, 'samples_evolved.csv')

    metadata_dict = {}
    phageORhost_all = []

    for line in open(metadata_path, 'r'):

        line_split = line.strip().split(',')

        if line_split[0] == 'unq.sample':
            continue

        unq_sample = line_split[0]
        phage_or_host = line_split[1]
        seed_bank = line_split[6]
        phage_treatment = line_split[7]
        subpop = line_split[8]

        reference = '%s%s.gbk' % (config.data_directory, line_split[-1].split('/')[1].split('.')[0])
        sample_name = '%s_%s_%s_%s_%s' % (unq_sample, phage_or_host, seed_bank, phage_treatment, subpop)

        metadata_dict[sample_name] = {}
        metadata_dict[sample_name]['unq_sample'] = unq_sample
        metadata_dict[sample_name]['phage_or_host'] = phage_or_host
        metadata_dict[sample_name]['seed_bank'] = seed_bank
        metadata_dict[sample_name]['phage_treatment'] = phage_treatment
        metadata_dict[sample_name]['subpop'] = subpop

        metadata_dict[sample_name]['transfer'] = int(unq_sample.split('-')[-1][1:])
        metadata_dict[sample_name]['replicate'] = int(unq_sample.split('-')[1][1:])

        metadata_dict[sample_name]['reference'] = reference


    return metadata_dict



def get_samples_from_metadata(phage_or_host_type, seed_bank_type, phage_treatment_type, subpop_type, replicate):

    metadata_dict = get_sample_metadata_dict()

    samples_all = []

    for key in metadata_dict.keys():

        if metadata_dict[key]['phage_or_host'] != phage_or_host_type:
            continue

        if metadata_dict[key]['seed_bank'] != seed_bank_type:
            continue

        if metadata_dict[key]['phage_treatment'] != phage_treatment_type:
            continue

        if metadata_dict[key]['subpop'] != subpop_type:
            continue

        if metadata_dict[key]['replicate'] != replicate:
            continue


        samples_all.append(key)

    return samples_all
