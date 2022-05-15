from __future__ import division
import numpy
import sys
from math import fabs
import glob, os, sys, re
import config
import utils
import parse_file
from Bio import SeqIO


#phage_or_host_type = 'host'
#seed_bank_type = 'long_seed_bank'
#phage_treatment_type = 'noPhage'
#subpop_type = 'revived_spore'
#replicate = 1



metadata_dict = utils.get_sample_metadata_dict()


def parse_and_annotate_pol_files_all_timepoints(samples, transfers):

    all_positions = []
    pol_dict = {}

    reference = metadata_dict[samples[0]]['reference']
    records = list(SeqIO.parse('%s.fa' % (reference.split('.')[0]), "fasta"))
    sequence = records[0].seq

    for sample_idx, sample in enumerate(samples):

        transfer = transfers[sample_idx]

        pol_path = "%spol/%s.pol" % (config.data_directory, sample)

        for line in open(pol_path, 'r'):

            if 'VERSION:0.4.40' in line:
                continue

            if '@SCFNAME' in line:
                continue

            if '@END_TABLE' in line:
                continue

            line_split = line.strip().split('\t')

            position = int(line_split[1])
            major_allele = line_split[3]
            minor_allele = line_split[4]
            coverage = float(line_split[5])
            frequency = float(line_split[7].split('/')[0])

            # mapgd starts counting at one
            reference_allele = sequence[position-1]
            # frequency in mapgd output is of major allele

            if major_allele == reference_allele:
                frequency = 1-frequency
                alt_allele = minor_allele

            else:
                alt_allele = major_allele


            if position not in pol_dict:
                pol_dict[position] = {}
                pol_dict[position]['alt_allele'] = alt_allele
                pol_dict[position]['transfers'] = []
                pol_dict[position]['frequency_trajectory'] = []
                pol_dict[position]['coverage_trajectory'] = []


            else:
                # skip new mutations at sites that were previously mutated
                if pol_dict[position]['alt_allele'] != alt_allele:
                    continue

            pol_dict[position]['transfers'].append(transfer)
            pol_dict[position]['frequency_trajectory'].append(frequency)
            pol_dict[position]['coverage_trajectory'].append(coverage)

            all_positions.append(position)

    # get annotation for these sites
    all_positions = list(set(all_positions))

    gene_data = parse_file.parse_gene_list(reference)
    gene_names, gene_start_positions, gene_end_positions, promoter_start_positions, promoter_end_positions, gene_sequences, strands, genes, features, protein_ids = gene_data

    position_gene_map, effective_gene_lengths, substitution_specific_synonymous_fraction = parse_file.create_annotation_map(reference=reference, gene_data=gene_data)


    annotated_mutations = []

    for position in all_positions:
        # get gene
        gene_name = parse_file.annotate_gene(position, position_gene_map)
        alt_allele = pol_dict[position]['alt_allele']


        if gene_name=='intergenic':
            var_type = 'noncoding'
            codon=None
            position_in_codon=None
            fold_count=None
        elif gene_name=='repeat':
            var_type = 'repeat'
            codon=None
            position_in_codon=None
            fold_count=None
        else:

            i = gene_names.index(gene_name)

            gene_start_position = gene_start_positions[i]
            gene_end_position = gene_end_positions[i]
            promoter_start_position = promoter_start_positions[i]
            promoter_end_position = promoter_end_positions[i]
            gene_sequence = gene_sequences[i]
            strand = strands[i]

            if position<gene_start_position or position>gene_end_position:
                var_type='noncoding' # (promoter)
                codon=None
                position_in_codon=None
                fold_count=None
            else:

                if gene_name.startswith('tRNA') or gene_name.startswith('rRNA'):
                    var_type='noncoding'
                    codon=None
                    position_in_codon=None
                    fold_count=None
                else:

                    # calculate position in gene
                    if strand=='forward':
                        position_in_gene = position-gene_start_position
                        oriented_gene_sequence = gene_sequence
                        new_base = alt_allele
                    else:
                        position_in_gene = gene_end_position-position
                        oriented_gene_sequence = parse_file.calculate_reverse_complement_sequence(gene_sequence)
                        new_base = parse_file.base_table[alt_allele]

                    # calculate codon start
                    codon_start = int(position_in_gene/3)*3
                    codon = oriented_gene_sequence[codon_start:codon_start+3]
                    codon_list = list(codon)
                    position_in_codon = position_in_gene%3
                    if (len(codon_list) == 0) or (len(set(codon_list) - set('ACGT')) != 0) :
                        var_type='unknown'
                        fold_count='unknown'
                        codon='unknown'
                        position_in_codon='unknown'
                        fold_count='unknown'

                    else:
                        codon_list[position_in_codon]=new_base
                        new_codon="".join(codon_list)
                        # one wrong bp in caulobacter reference genome
                        if parse_file.codon_table[codon]==parse_file.codon_table[new_codon]:
                            var_type='synonymous'
                        else:

                            if parse_file.codon_table[new_codon]=='!':
                                var_type='nonsense'
                            else:
                                var_type='missense'

                        # count fold
                        if position_in_codon >= 3:
                            fold_count = 'unknown'

                        else:
                            amino_acids = []
                            for base in ['A','C','T','G']:

                                if position_in_codon == 0:
                                    new_fold_codon = list(base) + codon_list[1:]
                                elif position_in_codon == 1:
                                    new_fold_codon = list(codon_list[0]) + list(base) + list(codon_list[2])
                                else:
                                    new_fold_codon = codon_list[:2] + list(base)

                                amino_acids.append(parse_file.codon_table["".join(new_fold_codon)])

                            fold_count=len(set(amino_acids))


        print_strings = [str(position), gene_name, alt_allele, var_type, str(codon), str(position_in_codon), str(fold_count)]
        # once for frequency
        for transfer in utils.transfers:
            if transfer not in pol_dict[position]['transfers']:
                print_strings.append(str(0))

            else:
                print_strings.append(str(pol_dict[position]['frequency_trajectory'][pol_dict[position]['transfers'].index(transfer)]))
        # once for coverage
        for transfer in utils.transfers:
            if transfer not in pol_dict[position]['transfers']:
                print_strings.append("Nan")

            else:
                print_strings.append(str(pol_dict[position]['coverage_trajectory'][pol_dict[position]['transfers'].index(transfer)]))


        annotated_mutations.append((position, ", ".join(print_strings)))



    file_name = samples[0]
    file_name = file_name.replace('-', '_')
    file_name = file_name.split('_')
    file_name.pop(2)
    file_name = '_'.join(file_name)

    output_filename = "%stimecourse_final/%s_annotated_timecourse.txt" % (config.data_directory, file_name)

    header = ['Position', 'Gene', 'Allele', 'Annotation', 'Codon', 'Position in codon', 'AA fold count', 'Freq:1', 'Freq:4', 'Freq:7', 'Freq:10', 'Freq:14', 'Cov:1', 'Cov:4', 'Cov:7', 'Cov:10', 'Cov:14']
    header = ", ".join(header)

    output_file = open(output_filename,"w")
    output_file.write(header)
    output_file.write("\n")
    for location, mutation_str in sorted(annotated_mutations,key=lambda x: x[0]):
        output_file.write(mutation_str)
        output_file.write("\n")
    output_file.close()





    #output_filename_template = config.data_directory() + "/data/timecourse_final/%s_annotated_timecourse.txt"



def annotate_all_line():

    for phage_or_host_type in utils.phage_or_host_types:
        for seed_bank_type in  utils.seed_bank_types:
            for phage_treatment_type in utils.phage_treatment_types:
                print(phage_or_host_type, seed_bank_type, phage_treatment_type)

                for subpop_type in utils.subpop_types:

                    for replicate in utils.replicates:

                        samples = utils.get_samples_from_metadata(phage_or_host_type, seed_bank_type, phage_treatment_type, subpop_type, replicate)
                        if len(samples) == 0:
                            continue

                        # get timepoints
                        transfer_all = []
                        for sample in samples:
                            transfer_all.append(metadata_dict[sample]['transfer'])

                        transfer_all, samples = zip(*sorted(zip(transfer_all, samples)))
                        samples = list(samples)
                        transfer_all = list(transfer_all)


                        parse_and_annotate_pol_files_all_timepoints(samples, transfer_all)



#annotate_all_line()
