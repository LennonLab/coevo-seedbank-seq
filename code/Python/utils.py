from __future__ import division
import numpy
import sys
from math import fabs
import glob, os, sys, re, glob
import config
import parse_file

import matplotlib.colors as clr


from Bio.Alphabet import generic_dna
from Bio import SeqIO
from BCBio import GFF


phage_or_host_types = ['phage', 'host']
seed_bank_types = ['short_seed_bank', 'long_seed_bank', 'no_seed_bank']
phage_treatment_types = ['noPhage', 'SPO1']
subpop_types = ['filtered_phage', 'revived_total', 'revived_spore']
references = ['delta6-ANC', 'dspoIIE-ANC', 'SPO1-ANC']
replicates = [1, 2, 3]


color_dict = {'no_seed_bank':'#87CEEB', 'short_seed_bank': '#FFA500', 'long_seed_bank':'#FF6347'}
marker_dict = {'noPhage': 'o', 'SPO1': '^'}
line_dict = {'noPhage': '--', 'SPO1': ':'}

background = ['SNO', 'WLCt', 'WLO', 'SNCt', 'WSCt', 'WSO']
#to_ignore = [('host', 'revived_spore')]

transfers = [1,4,7,10,14]

min_n_non_zero_freqs = 3

high_coverage_idx = numpy.asarray([0,1,4])
low_coverage_idx = numpy.asarray([2,3])

seed_bank_types_format_dict = {'short_seed_bank': 'Short seed bank', 'long_seed_bank': 'Long seed bank', 'no_seed_bank': 'No seed bank'}
phage_treatment_types_format_dict = {'noPhage': 'no phage', 'SPO1': 'with phage'}
subpop_types = {'filtered_phage': 'Phage', 'revived_total': 'All cells', 'revived_spore': 'Seed bank'}




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





def get_sample_metadata_breseq_dict():

    sample_metadata_dict = get_sample_metadata_dict()
    sample_metadata_dict_keys = list(sample_metadata_dict.keys())

    subpop_breseq_to_metadata_dict = {'phage': 'filtered_phage', 'rS': 'revived_spore', 'rV': 'revived_total'}

    metadata_breseq_dict = {}

    for file in os.listdir("%sbreseq_output/" % config.data_directory):

        file_strip = file.strip()
        uniq_label, subpop_breseq = file_strip.split('.')[0].rsplit('-', 1)

        subpop_metadata = subpop_breseq_to_metadata_dict[subpop_breseq]
        sample = [e for e in sample_metadata_dict_keys if e.startswith(uniq_label+'_') and e.endswith(subpop_metadata)][0]

        metadata_breseq_dict[sample] = file_strip

    return metadata_breseq_dict





def get_samples_from_metadata(phage_or_host_type, seed_bank_type, phage_treatment_type, subpop_type, replicate):

    metadata_dict = get_sample_metadata_dict()

    samples_all = []
    #print(metadata_dict['SNCt-L3-T7_host_no_seed_bank_noPhage_revived_total'])

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




def load_annotated_mapgd(phage_or_host_type, seed_bank_type, phage_treatment_type, subpop_type, replicate):

    population = 'L%d_%s_%s_%s_%s_annotated_timecourse.txt' % (replicate, phage_or_host_type, seed_bank_type, phage_treatment_type, subpop_type)

    #print(population)

    annotated_mapgd_dict = {}

    partial_matches = glob.glob('%s/timecourse_final/*%s' % (config.data_directory, population))

    if len(partial_matches) == 0:
        exists = False
    else:
        exists = True

        filename = partial_matches[0]
        file = open(filename,"r")
        header = file.readline() # header
        header = header.strip().split('/t')
        for line in file:

            line = line.strip().split(', ')
            position = int(line[0])
            annotated_mapgd_dict[position] = {}

            annotated_mapgd_dict[position]['gene'] = line[1]
            annotated_mapgd_dict[position]['allele'] = line[2]
            annotated_mapgd_dict[position]['annotation'] = line[3]
            annotated_mapgd_dict[position]['codon'] = line[4]

            if line[5] == 'None':
                annotated_mapgd_dict[position]['position_in_codon'] = 'None'
                annotated_mapgd_dict[position]['aa_fold_count'] = 'None'
            elif line[5] == 'unknown':
                annotated_mapgd_dict[position]['position_in_codon'] = 'unknown'
                annotated_mapgd_dict[position]['aa_fold_count'] = 'unknown'
            else:
                annotated_mapgd_dict[position]['position_in_codon'] = int(line[5])
                annotated_mapgd_dict[position]['aa_fold_count'] = int(line[6])

            frequency_trajectory = numpy.asarray([float(f) for f in line[7:12]])
            annotated_mapgd_dict[position]['frequency_trajectory'] = frequency_trajectory

            coverage_trajectory = numpy.asarray([float(f) for f in line[12:]])
            annotated_mapgd_dict[position]['coverage_trajectory'] = coverage_trajectory


        file.close()


    return exists, annotated_mapgd_dict





def load_annotated_breseq(phage_or_host_type, seed_bank_type, phage_treatment_type, subpop_type, replicate):

    population = 'L%d_%s_%s_%s_%s_annotated_timecourse.txt' % (replicate, phage_or_host_type, seed_bank_type, phage_treatment_type, subpop_type)

    annotated_mapgd_dict = {}

    partial_matches = glob.glob('%s/timecourse_final_breseq/*%s' % (config.data_directory, population))

    if len(partial_matches) == 0:
        exists = False
    else:
        exists = True

        filename = partial_matches[0]
        file = open(filename,"r")
        header = file.readline() # header
        header = header.strip().split('/t')
        for line in file:

            line = line.strip().split(', ')
            position = int(line[0])
            annotated_mapgd_dict[position] = {}

            annotated_mapgd_dict[position]['mutation_type'] = line[1]
            annotated_mapgd_dict[position]['gene'] = line[2]
            annotated_mapgd_dict[position]['allele'] = line[3]
            annotated_mapgd_dict[position]['annotation'] = line[4]
            annotated_mapgd_dict[position]['codon'] = line[5]

            if line[5] == 'None':
                annotated_mapgd_dict[position]['position_in_codon'] = 'None'
                annotated_mapgd_dict[position]['aa_fold_count'] = 'None'
            elif line[5] == 'unknown':
                annotated_mapgd_dict[position]['position_in_codon'] = 'unknown'
                annotated_mapgd_dict[position]['aa_fold_count'] = 'unknown'
            else:
                annotated_mapgd_dict[position]['position_in_codon'] = int(line[6])
                annotated_mapgd_dict[position]['aa_fold_count'] = int(line[7])

            frequency_trajectory = numpy.asarray([float(f) for f in line[8:13]])
            annotated_mapgd_dict[position]['frequency_trajectory'] = frequency_trajectory

            coverage_trajectory = numpy.asarray([float(f) for f in line[13:]])
            annotated_mapgd_dict[position]['coverage_trajectory'] = coverage_trajectory

        file.close()


    return exists, annotated_mapgd_dict




def lighten_color(color, amount=0.5):
    """
    Lightens the given color by multiplying (1-luminosity) by the given amount.
    Input can be matplotlib color string, hex string, or RGB tuple.
    Examples:
    >> lighten_color('g', 0.3)
    >> lighten_color('#F034A3', 0.6)
    >> lighten_color((.3,.55,.1), 0.5)
    """
    import matplotlib.colors as mc
    import colorsys
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])



def mut_freq_colormap():
    #cmap = clr.LinearSegmentedColormap.from_list('Zissou1', ["#3B9AB2", "#78B7C5", "#EBCC2A", "#E1AF00", "#F21A00"], N=256)
    cmap = clr.LinearSegmentedColormap.from_list('Darjeeling1', ["#FF0000", "#00A08A", "#F2AD00", "#F98400", "#5BBCD6"], N=256)

    # sample from cmap using uniform dist b/w 0 and 1
    u = numpy.random.uniform()
    rgb = '#%02x%02x%02x' % tuple([int(x * 100) for x in list(cmap(u))[:-1]])
    #tuple([int(x * 100) for x in list(cmap(u))[:-1]])
    # RGB six digit code
    return rgb



def get_gene_intersection():

    reference_1 = "%sdspoIIE-ANC.gbk" % config.data_directory
    gene_data_1 = parse_file.parse_gene_list(reference_1)
    gene_names_1, gene_start_positions_1, gene_end_positions_1, promoter_start_positions_1, promoter_end_positions_1, gene_sequences_1, strands_1, genes_1, features_1, protein_ids_1 = gene_data_1

    reference_2 = "%sdelta6-ANC.gbk" % config.data_directory
    gene_data_2 = parse_file.parse_gene_list(reference_2)
    gene_names_2, gene_start_positions_2, gene_end_positions_2, promoter_start_positions_2, promoter_end_positions_2, gene_sequences_2, strands_2, genes_2, features_2, protein_ids_2 = gene_data_2
    #gene_names_1 = numpy.asarray(gene_names_1)

    gene_names_intersection = set(gene_names_1).intersection(set(gene_names_2))
    gene_names_union = set(gene_names_1).union(set(gene_names_2))

    return gene_names_intersection
