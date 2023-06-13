from __future__ import division
import numpy
import sys
from math import fabs
import glob, os, sys, re
import config
import utils
import parse_file
from Bio import SeqIO


metadata_dict = utils.get_sample_metadata_dict()
sample_metadata_breseq_dict = utils.get_sample_metadata_breseq_dict()

breseq_types_to_keep = ['SNP', 'DEL', 'INS']
breseq_non_snp_types = ['DEL', 'INS']


def parse_and_annotate_breseq_files_all_timepoints(breseq_samples, samples, transfers):

    pol_dict = {}

    reference = metadata_dict[samples[0]]['reference']
    records = list(SeqIO.parse('%s.fa' % (reference.split('.')[0]), "fasta"))
    sequence = records[0].seq

    all_positions = []

    for breseq_sample_idx, breseq_sample in enumerate(breseq_samples):

        transfer = transfers[breseq_sample_idx]

        breseq_sample_path = "%sbreseq_output/%s" % (config.data_directory, breseq_sample)

        for line in open(breseq_sample_path, 'r'):

            line = line.strip().split('\t')
            line_type = line[0]
            if line_type not in breseq_types_to_keep:
                continue

            position = int(line[4])
            alt_allele = line[5]
            frequency_string = line[6].split('=')[1]
            if frequency_string == 'NA':
                frequency_string = '0'
            frequency = float(frequency_string)

            # breseq starts counting at one
            reference_allele = sequence[position-1]

            if position not in pol_dict:
                pol_dict[position] = {}
                pol_dict[position]['mutation_type'] = line_type
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

            all_positions.append(position)

        for line in open(breseq_sample_path, 'r'):

            line = line.strip().split('\t')
            if line[0] == 'RA':

                position = int(line[4])
                if position in all_positions:

                    coverage = [x for x in line if 'total_cov=' in x][0].split('=')[1]
                    coverage_1, coverage_2 = coverage.split('/')
                    coverage_total = int(coverage_1) + int(coverage_2)

                    pol_dict[position]['coverage_trajectory'].append(coverage_total)


    # go back through and set frequencies of observations with no coverage succeeding fixations to zero
    all_positions = list(pol_dict.keys())
    for position in all_positions:

        t = pol_dict[position]['frequency_trajectory']
        t = numpy.asarray(t)

        if sum(t) == 0:
            del pol_dict[position]
            continue

        if len(t) < 2:
            del pol_dict[position]
            continue

        if len(pol_dict[position]['coverage_trajectory']) < 2:
            del pol_dict[position]
            continue

        # not real mutations
        if t[0] > 0.8:
            del pol_dict[position]
            continue

        if sum(t==0.5) == len(t):
            del pol_dict[position]
            continue

        zip_list = list(zip(t.tolist(), pol_dict[position]['coverage_trajectory']))

        if (0.5, 0) in zip_list:
            del pol_dict[position]
            continue


        #else:
        #    # keep if theres at least one fixation and one polymophic observation
        #    if (sum(t==1)>0) and (sum((t>0.01) & (t<0.99)) >0):
        #        init_idx = numpy.argwhere(t==1)[0][0]
        #        # set all future samples to 1
        #        print(t)
        #        t[init_idx+1:] = 1
        #        print(t)


    all_positions = list(pol_dict.keys())
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
        mutation_type = pol_dict[position]['mutation_type']
        # insertion or deletion
        if mutation_type in breseq_non_snp_types:
            var_type = alt_allele
            codon=None
            position_in_codon=None
            fold_count=None
            position_in_gene=None
            codon_start=None

        elif gene_name=='intergenic':
            var_type = 'noncoding'
            codon=None
            position_in_codon=None
            fold_count=None
            position_in_gene=None
            codon_start=None
        
        elif gene_name=='repeat':
            var_type = 'repeat'
            codon=None
            position_in_codon=None
            fold_count=None
            position_in_gene=None
            codon_start=None
        
        else:

            #position_in_gene = position_in_gene_map[position]
            i = gene_names.index(gene_name)

            gene_start_position = gene_start_positions[i]
            gene_end_position = gene_end_positions[i]
            promoter_start_position = promoter_start_positions[i]
            promoter_end_position = promoter_end_positions[i]
            gene_sequence = gene_sequences[i]
            strand = strands[i]

            if position < gene_start_position or position>gene_end_position:
                var_type='noncoding' # (promoter)
                codon=None
                position_in_codon=None
                fold_count=None
                position_in_gene=None
                codon_start=None
            else:

                #position_in_gene = position-gene_start_position
                if gene_name.startswith('tRNA') or gene_name.startswith('rRNA'):
                    var_type='noncoding'
                    codon=None
                    position_in_codon=None
                    fold_count=None
                    position_in_gene=None
                    codon_start=None
                else:

                    # calculate position in gene
                    if strand=='forward':
                        position_in_gene = (position-gene_start_position) -1
                        oriented_gene_sequence = gene_sequence
                        new_base = alt_allele
                    else:
                        position_in_gene = (gene_end_position-position) -1
                        oriented_gene_sequence = parse_file.calculate_reverse_complement_sequence(gene_sequence)

                        if alt_allele not in parse_file.base_table:
                            new_base = 'NA'

                        else:
                            new_base = parse_file.base_table[alt_allele]


           
                    # for printing to file
                    alt_allele = new_base
                    
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

                        if new_codon not in parse_file.codon_table:
                            var_type='unknown'

                        elif parse_file.codon_table[codon]==parse_file.codon_table[new_codon]:
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


        print_strings = [str(position), mutation_type, gene_name, alt_allele, str(position_in_gene), str(codon_start), var_type, str(codon), str(position_in_codon), str(fold_count)]
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
                if pol_dict[position]['mutation_type'] not in breseq_non_snp_types:
                    print_strings.append(str(pol_dict[position]['coverage_trajectory'][pol_dict[position]['transfers'].index(transfer)]))

        annotated_mutations.append((position, ", ".join(print_strings)))


    file_name = samples[0]
    file_name = file_name.replace('-', '_')
    file_name = file_name.split('_')
    file_name.pop(2)
    file_name = '_'.join(file_name)

    output_filename = "%stimecourse_final_breseq/%s_annotated_timecourse.txt" % (config.data_directory, file_name)

    header = ['Position', 'Mutation type', 'Gene', 'Allele', 'Position in gene', 'Codon position in gene', 'Annotation', 'Codon', 'Position in codon', 'AA fold count', 'Freq:1', 'Freq:4', 'Freq:7', 'Freq:10', 'Freq:14', 'Cov:1', 'Cov:4', 'Cov:7', 'Cov:10', 'Cov:14']
    #header = ['Position', 'Gene', 'Allele', 'Annotation', 'Codon', 'Position in codon', 'AA fold count', 'Freq:1', 'Freq:4', 'Freq:7', 'Freq:10', 'Freq:14']
    header = ", ".join(header)

    output_file = open(output_filename,"w")
    output_file.write(header)
    output_file.write("\n")
    for location, mutation_str in sorted(annotated_mutations,key=lambda x: x[0]):
        output_file.write(mutation_str)
        output_file.write("\n")
    output_file.close()






def annotate_all_line():

    for phage_or_host_type in utils.phage_or_host_types:
        for seed_bank_type in  utils.seed_bank_types:
            for phage_treatment_type in utils.phage_treatment_types:
                for subpop_type in utils.subpop_types:
                    for replicate in utils.replicates:

                        samples = utils.get_samples_from_metadata(phage_or_host_type, seed_bank_type, phage_treatment_type, subpop_type, replicate)
                        if len(samples) == 0:
                            continue

                        # Daniel needs to rerun this sample
                        #if ('WLO-L3' in samples[0]) and ('_filtered_phage' in samples[0]):
                        #    continue

                        breseq_samples = [sample_metadata_breseq_dict[s] for s in samples]
                        # get timepoints
                        transfer_all = []
                        for sample in samples:
                            transfer_all.append(metadata_dict[sample]['transfer'])

                        transfer_all, breseq_samples, samples = zip(*sorted(zip(transfer_all, breseq_samples, samples)))
                        transfer_all = list(transfer_all)
                        breseq_samples = list(breseq_samples)
                        samples = list(samples)

                        parse_and_annotate_breseq_files_all_timepoints(breseq_samples, samples, transfer_all)





# Daniel needs to rerun this sample
#if ('WLO-L3' in samples[0]) and ('_filtered_phage' in samples[0]):
#    continue

#breseq_samples = [sample_metadata_breseq_dict[s] for s in samples]
# get timepoints
#transfer_all = []
#for sample in samples:
#    transfer_all.append(metadata_dict[sample]['transfer'])

#transfer_all, breseq_samples, samples = zip(*sorted(zip(transfer_all, breseq_samples, samples)))
#transfer_all = list(transfer_all)
#breseq_samples = list(breseq_samples)
#samples = list(samples)

#parse_and_annotate_breseq_files_all_timepoints(breseq_samples, samples, transfer_all)

annotate_all_line()
