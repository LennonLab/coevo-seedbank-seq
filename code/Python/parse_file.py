from __future__ import division
import numpy
import sys
from math import fabs
import glob, os, sys, re
from bz2 import BZ2File
from Bio import SeqIO
import config

#import phylo_tools as pt
#import timecourse_utils
#data_directory='data_files/'
#figure_directory='manuscript/figures/'

#default_min_depth=5

#Y = pYrimidines
#R = puRines
#S = strong ineractions, C or G
#W = weak interactions, A or T
#K = Ketones, G or T
#M = aMino groups, C or A
bases_to_skip = ['K', 'S', 'R', 'N', 'Y', 'M', 'W']
base_table = {'A':'T','T':'A','G':'C','C':'G',
            'Y':'R', 'R':'Y', 'S':'W', 'W':'S', 'K':'M', 'M':'K', 'N':'N'}

codon_table = { 'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'CGT': 'R', 'CGC': 'R', 'CGA':'R',
'CGG':'R', 'AGA':'R', 'AGG':'R', 'AAT':'N', 'AAC':'N', 'GAT':'D', 'GAC':'D', 'TGT':'C', 'TGC':'D',
'CAA':'Q', 'CAG':'Q', 'GAA':'E', 'GAG':'E', 'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G', 'CAT':'H',
'CAC':'H', 'ATT':'I', 'ATC':'I', 'ATA':'I', 'TTA':'L', 'TTG':'L', 'CTT':'L', 'CTC':'L', 'CTA':'L',
'CTG':'L', 'AAA':'K', 'AAG':'K', 'ATG':'M', 'TTT':'F', 'TTC':'F', 'CCT':'P', 'CCC':'P', 'CCA':'P',
'CCG':'P', 'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S', 'AGT':'S', 'AGC':'S', 'ACT':'T', 'ACC':'T',
'ACA':'T', 'ACG':'T', 'TGG':'W', 'TAT':'Y', 'TAC':'Y', 'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V',
'TAA':'!', 'TGA':'!', 'TAG':'!'}#, 'KTC':'F', 'KAC':'Y', 'KCC':'A', 'KGC':'D'}


var_types = ['synonymous','missense','nonsense','noncoding','indel','sv']


# calculate number of synonymous opportunities for each codon
codon_synonymous_opportunity_table = {}
for codon in codon_table.keys():
    codon_synonymous_opportunity_table[codon] = {}
    for i in range(0,3):
        codon_synonymous_opportunity_table[codon][i] = -1 # since G->G is by definition synonymous, but we don't want to count it
        codon_list = list(codon)
        for base in ['A','C','T','G']:
            codon_list[i]=base
            new_codon = "".join(codon_list)
            if 'K' in new_codon:
                continue
            if codon_table[codon]==codon_table[new_codon]:
                # synonymous!
                codon_synonymous_opportunity_table[codon][i]+=1

bases = set(['A','C','T','G'])
substitutions = []
for b1 in bases:
    for b2 in bases:
        if b2==b1:
            continue

        substitutions.append( '%s->%s' % (b1,b2) )

codon_synonymous_substitution_table = {}
codon_nonsynonymous_substitution_table = {}
for codon in codon_table.keys():
    codon_synonymous_substitution_table[codon] = [[],[],[]]
    codon_nonsynonymous_substitution_table[codon] = [[],[],[]]

    for i in range(0,3):
        reference_base = codon[i]

        codon_list = list(codon)
        for derived_base in ['A','C','T','G']:
            if derived_base==reference_base:
                continue
            substitution = '%s->%s' % (reference_base, derived_base)
            codon_list[i]=derived_base
            new_codon = "".join(codon_list)
            if codon_table[codon]==codon_table[new_codon]:
                # synonymous!
                codon_synonymous_substitution_table[codon][i].append(substitution)
            else:
                codon_nonsynonymous_substitution_table[codon][i].append(substitution)





def get_genome_size(reference):
    genome_size_dict = {"delta6-ANC": 3876919,
                        "dspoIIE-ANC": 3874590,
                        "SPO1-ANC": 132564}


    r = reference.split('/')[-1].split('.')[0]

    if 'dspoIIE' in reference:
        size =  3874590

    elif 'delta6' in reference:
        size = 3876919

    elif ('SPO1' in reference):
        size = 132564

    else:
        print("don't know the reference!")


    #return genome_size_dict[r]
    return size





def calculate_synonymous_nonsynonymous_target_sizes(reference):
    position_gene_map, effective_gene_lengths, substitution_specific_synonymous_fraction  = create_annotation_map(reference=reference)
    return effective_gene_lengths['synonymous'], effective_gene_lengths['nonsynonymous'], substitution_specific_synonymous_fraction



def calculate_reverse_complement_sequence(dna_sequence):
    return "".join(base_table[base] for base in dna_sequence[::-1])




#####################################################################
#
# Reads through the Genbank file for the reference and
# compiles a list of genes, tRNAs, etc.
#
#####################################################################
def parse_gene_list(reference, reference_sequence=None):
    gene_names = []
    start_positions = []
    end_positions = []
    promoter_start_positions = []
    promoter_end_positions = []
    gene_sequences = []
    strands = []
    genes = []
    features = []
    protein_ids = []

    filename = reference
    gene_features = ['CDS', 'tRNA', 'rRNA', 'ncRNA', 'tmRNA']
    recs = [rec for rec in SeqIO.parse(filename, "gb")]
    count_riboswitch = 0
    for rec in recs:
        reference_sequence = rec.seq
        contig = rec.annotations['accessions'][0]
        for feat in rec.features:
            if 'pseudo' in list((feat.qualifiers.keys())):
                continue
            if (feat.type == "source") or (feat.type == "gene"):
                continue

            locations = re.findall(r"[\w']+", str(feat.location))
            if feat.type in gene_features:

                if 'locus_tag' not in feat.qualifiers:
                    locus_tag = feat.qualifiers['Alias'][0]
                else:
                    locus_tag = feat.qualifiers['locus_tag'][0]



            elif (feat.type=="regulatory"):
                if 'regulatory_class' not in feat.qualifiers:
                    locus_tag = 'regulatory_region'
                else:
                    locus_tag = feat.qualifiers["regulatory_class"][0] + '_' + str(count_riboswitch)
                count_riboswitch += 1
            else:
                continue
            # for frameshifts, split each CDS seperately and merge later
            # Fix this for Deinococcus, it has a frameshift in three pieces
            split_list = []
            if 'join' in locations:
                location_str = str(feat.location)
                minus_position = []
                if '-' in location_str:
                    minus_position = [r.start() for r in re.finditer('-', location_str)]
                pos_position = []

                if '+' in location_str:
                    #if taxon == 'D':
                    #    pos_position = [pos for pos, char in enumerate(location_str) if char == '+']
                    #elif taxon == 'J':
                    #    pos_position = [pos for pos, char in enumerate(location_str) if char == '+']
                    #else:
                    #    pos_position = [r.start() for r in re.finditer('+', location_str)]
                    pos_position = [r.start() for r in re.finditer('+', location_str)]


                if len(minus_position) + len(pos_position) == 2:
                    if len(minus_position) == 2:
                        strand_symbol_one = '-'
                        strand_symbol_two = '-'
                    elif len(pos_position) == 2:
                        strand_symbol_one = '+'
                        strand_symbol_two = '+'
                    else:
                        # I don't think this is possible, but might as well code it up
                        if minus_position[0] < pos_position[0]:
                            strand_symbol_one = '-'
                            strand_symbol_two = '+'
                        else:
                            strand_symbol_one = '+'
                            strand_symbol_two = '-'

                    start_one = int(locations[1])
                    stop_one = int(locations[2])
                    start_two = int(locations[3])
                    stop_two = int(locations[4])

                    locus_tag1 = locus_tag + '_1'
                    locus_tag2 = locus_tag + '_2'

                    split_list.append([locus_tag1, start_one, stop_one, strand_symbol_one])
                    split_list.append([locus_tag2, start_two, stop_two, strand_symbol_two])

                else:
                    if len(pos_position) == 3:
                        strand_symbol_one = '+'
                        strand_symbol_two = '+'
                        strand_symbol_three = '+'
                    start_one = int(locations[1])
                    stop_one = int(locations[2])
                    start_two = int(locations[3])
                    stop_two = int(locations[4])
                    start_three = int(locations[5])
                    stop_three = int(locations[6])

                    locus_tag1 = locus_tag + '_1'
                    locus_tag2 = locus_tag + '_2'
                    locus_tag3 = locus_tag + '_3'

                    split_list.append([locus_tag1, start_one, stop_one, strand_symbol_one])
                    split_list.append([locus_tag2, start_two, stop_two, strand_symbol_two])
                    split_list.append([locus_tag3, start_three, stop_three, strand_symbol_three])


            else:
                strand_symbol = str(feat.location)[-2]
                start = int(locations[0])
                stop = int(locations[1])
                split_list.append([locus_tag, start, stop, strand_symbol])

            for split_item in split_list:
                locus_tag = split_item[0]
                start = split_item[1]
                stop = split_item[2]
                strand_symbol = split_item[3]


                if feat.type == 'CDS':
                    # biopython accounts for the -1
                    #gene_sequence = reference_sequence[start-1:stop]
                    gene_sequence = str(reference_sequence[start:stop])

                    # gene_sequence and biopython_gene_sequence are equivalent
                    gene_sequence_ = feat.location.extract(rec).seq
                    #biopython_gene_sequence = biopython_gene_sequence.strip()

                else:
                    gene_sequence = ""


                if 'gene' in list((feat.qualifiers.keys())):
                    gene = feat.qualifiers['gene'][0]
                else:
                    gene = ""

                if 'protein_id' in list((feat.qualifiers.keys())):
                    protein_id = feat.qualifiers['protein_id'][0]
                else:
                    protein_id = ""


                if strand_symbol == '+':
                    promoter_start = start - 100 # by arbitrary definition, we treat the 100bp upstream as promoters
                    promoter_end = start - 1
                    strand = 'forward'
                else:
                    promoter_start = stop+1
                    promoter_end = stop+100
                    strand = 'reverse'


                if gene_sequence!="" and (not len(gene_sequence)%3==0):
                    print(locus_tag, start, "Not a multiple of 3")
                    continue

                # dont need to check if gene names are unique because we're using
                # locus tags

                start_positions.append(start)
                end_positions.append(stop)
                promoter_start_positions.append(promoter_start)
                promoter_end_positions.append(promoter_end)
                gene_names.append(locus_tag)
                gene_sequences.append(gene_sequence)
                strands.append(strand)
                genes.append(gene)
                features.append(feat.type)
                protein_ids.append(protein_id)

    gene_names, start_positions, end_positions, promoter_start_positions, promoter_end_positions, gene_sequences, strands, genes, features, protein_ids = (list(x) for x in zip(*sorted(zip(gene_names, start_positions, end_positions, promoter_start_positions, promoter_end_positions, gene_sequences, strands, genes, features, protein_ids), key=lambda pair: pair[1])))

    return gene_names, numpy.array(start_positions), numpy.array(end_positions), numpy.array(promoter_start_positions), numpy.array(promoter_end_positions), gene_sequences, strands, genes, features, protein_ids




# get the annotations from genbank files
def get_gene_annotation_dict(reference, reference_sequence=None):

    annotation_dict = {}

    filename = reference
    gene_features = ['CDS', 'tRNA', 'rRNA', 'ncRNA', 'tmRNA']
    recs = [rec for rec in SeqIO.parse(filename, "genbank")]
    count_riboswitch = 0
    for rec in recs:
        reference_sequence = rec.seq
        contig = rec.annotations['accessions'][0]
        for feat in rec.features:
            if 'pseudo' in list((feat.qualifiers.keys())):
                continue
            if (feat.type == "source") or (feat.type == "gene"):
                continue

            locations = re.findall(r"[\w']+", str(feat.location))
            if feat.type in gene_features:
                locus_tag = feat.qualifiers['locus_tag'][0]
            elif (feat.type=="regulatory"):
                locus_tag = feat.qualifiers["regulatory_class"][0] + '_' + str(count_riboswitch)
                count_riboswitch += 1
            else:
                continue
            # for frameshifts, split each CDS seperately and merge later
            # Fix this for Deinococcus, it has a frameshift in three pieces
            split_list = []
            if 'join' in locations:
                location_str = str(feat.location)
                minus_position = []
                if '-' in location_str:
                    minus_position = [r.start() for r in re.finditer('-', location_str)]
                pos_position = []

                if '+' in location_str:
                    if taxon == 'D':
                        pos_position = [pos for pos, char in enumerate(location_str) if char == '+']
                    elif taxon == 'J':
                        pos_position = [pos for pos, char in enumerate(location_str) if char == '+']
                    else:
                        pos_position = [r.start() for r in re.finditer('+', location_str)]


                if len(minus_position) + len(pos_position) == 2:
                    if len(minus_position) == 2:
                        strand_symbol_one = '-'
                        strand_symbol_two = '-'
                    elif len(pos_position) == 2:
                        strand_symbol_one = '+'
                        strand_symbol_two = '+'
                    else:
                        # I don't think this is possible, but might as well code it up
                        if minus_position[0] < pos_position[0]:
                            strand_symbol_one = '-'
                            strand_symbol_two = '+'
                        else:
                            strand_symbol_one = '+'
                            strand_symbol_two = '-'

                    start_one = int(locations[1])
                    stop_one = int(locations[2])
                    start_two = int(locations[3])
                    stop_two = int(locations[4])

                    locus_tag1 = locus_tag + '_1'
                    locus_tag2 = locus_tag + '_2'

                    split_list.append([locus_tag1, start_one, stop_one, strand_symbol_one])
                    split_list.append([locus_tag2, start_two, stop_two, strand_symbol_two])

                else:
                    if len(pos_position) == 3:
                        strand_symbol_one = '+'
                        strand_symbol_two = '+'
                        strand_symbol_three = '+'
                    start_one = int(locations[1])
                    stop_one = int(locations[2])
                    start_two = int(locations[3])
                    stop_two = int(locations[4])
                    start_three = int(locations[5])
                    stop_three = int(locations[6])

                    locus_tag1 = locus_tag + '_1'
                    locus_tag2 = locus_tag + '_2'
                    locus_tag3 = locus_tag + '_3'

                    split_list.append([locus_tag1, start_one, stop_one, strand_symbol_one])
                    split_list.append([locus_tag2, start_two, stop_two, strand_symbol_two])
                    split_list.append([locus_tag3, start_three, stop_three, strand_symbol_three])


            else:
                strand_symbol = str(feat.location)[-2]
                start = int(locations[0])
                stop = int(locations[1])
                split_list.append([locus_tag, start, stop, strand_symbol])

            for split_item in split_list:
                locus_tag = split_item[0]
                start = split_item[1]
                stop = split_item[2]
                strand_symbol = split_item[3]


                if feat.type == 'CDS':
                    #gene_sequence = reference_sequence[start-1:stop]
                    gene_sequence = str(reference_sequence[start:stop])
                else:
                    gene_sequence = ""


                if 'gene' in list((feat.qualifiers.keys())):
                    gene = feat.qualifiers['gene'][0]
                else:
                    gene = ""

                if 'protein_id' in list((feat.qualifiers.keys())):
                    protein_id = feat.qualifiers['protein_id'][0]
                else:
                    protein_id = ""


                if strand_symbol == '+':
                    promoter_start = start - 100 # by arbitrary definition, we treat the 100bp upstream as promoters
                    promoter_end = start - 1
                    strand = 'forward'
                else:
                    promoter_start = stop+1
                    promoter_end = stop+100
                    strand = 'reverse'


                if gene_sequence!="" and (not len(gene_sequence)%3==0):
                    print(locus_tag, start, "Not a multiple of 3")
                    continue

                if 'product' in feat.qualifiers:

                    function = feat.qualifiers['product'][0]

                else:
                    function = 'hypothetical protein'

                annotation_dict[locus_tag] = function

    return annotation_dict









def annotate_gene(position, position_gene_map):

    if position in position_gene_map:
        gene_name = position_gene_map[position]
    else:
        gene_name = 'intergenic'

    return gene_name



def annotate_variant(position, allele, gene_data, position_gene_map):

    gene_names, gene_start_positions, gene_end_positions, promoter_start_positions, promoter_end_positions, gene_sequences, strands, genes, features, protein_ids = gene_data

    # get gene
    gene_name = annotate_gene(position, position_gene_map)

    if allele.startswith('Depth'):
        var_type = 'unknown'
        codon=None
        position_in_codon=None
        fold_count=None
    elif allele.startswith('MOB') or allele.startswith('junction'):
        var_type = 'sv'
        codon=None
        position_in_codon=None
        fold_count=None
    elif allele.startswith('indel'):
        var_type = 'indel'
        codon=None
        position_in_codon=None
        fold_count=None
    elif allele[1:3]=='->':
        # a SNP, so annotate it
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
            # must be in a real gene
            # so get it
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
                        new_base = allele[3]
                    else:
                        position_in_gene = gene_end_position-position
                        oriented_gene_sequence = calculate_reverse_complement_sequence(gene_sequence)
                        new_base = base_table[allele[3]]

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
                        if codon_table[codon]==codon_table[new_codon]:
                            var_type='synonymous'
                        else:

                            if codon_table[new_codon]=='!':
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

                                amino_acids.append(codon_table["".join(new_fold_codon)])

                            fold_count=len(set(amino_acids))

    else:
        sys.stderr.write("Unknown: %s\n" % allele)
        var_type='unknown'
        codon='unknown'
        position_in_codon='unknown'
        fold_count='unknown'

    return gene_name, var_type, codon, position_in_codon, fold_count








def create_annotation_map(reference, gene_data=None):

    if gene_data==None:
        gene_data = parse_gene_list(reference)

    gene_names, gene_start_positions, gene_end_positions, promoter_start_positions, promoter_end_positions, gene_sequences, strands, genes, features, protein_ids = gene_data
    position_gene_map = {}
    gene_position_map = {}
    # new
    gene_feature_map = {}

    #position_in_gene_map = {}

    # then greedily annotate genes at remaining sites
    overlapping_genes = []
    for gene_name, feature, start, end in zip(gene_names, features, gene_start_positions, gene_end_positions):
        gene_feature_map[gene_name] = feature

        if start>end:
            print(gene_name, feature, start, end)
        
        for position in range(start,end+1):
            if position not in position_gene_map:
                position_gene_map[position] = gene_name

                if gene_name not in gene_position_map:
                    gene_position_map[gene_name]=[]
                
                gene_position_map[gene_name].append(position)

            #else:
            #    overlapping_genes.append(gene_name)
            #    overlapping_genes.append(position_gene_map[position])

            #if position not in position_in_gene_map:
            #    position_in_gene_map[position] = position-start
    
    #print(len(set(overlapping_genes)), len(gene_names))
    #print(list(set(overlapping_genes)))

    # remove 'partial' genes that have < 10bp unmasked sites
    for gene_name in list(sorted(gene_position_map.keys())):
        if len(gene_position_map[gene_name]) < 10:
            for position in gene_position_map[gene_name]:
                position_gene_map[position] = 'repeat'
            del gene_position_map[gene_name]

    # count up number of synonymous opportunities
    effective_gene_synonymous_sites = {}
    effective_gene_nonsynonymous_sites = {}

    substitution_specific_synonymous_sites = {substitution: 0 for substitution in substitutions}
    substitution_specific_nonsynonymous_sites = {substitution: 0 for substitution in substitutions}

    for gene_name, start, end, gene_sequence, strand in zip(gene_names, gene_start_positions, gene_end_positions, gene_sequences, strands):

        if gene_name not in gene_position_map:
            continue

        if strand=='forward':
            oriented_gene_sequence = gene_sequence
        else:
            oriented_gene_sequence = calculate_reverse_complement_sequence(gene_sequence)

        for position in gene_position_map[gene_name]:

            if gene_name not in effective_gene_synonymous_sites:
                effective_gene_synonymous_sites[gene_name]=0
                effective_gene_nonsynonymous_sites[gene_name]=0

            if 'CDS' not in gene_feature_map[gene_name]:
                continue

            else:
                # calculate position in gene
                if strand=='forward':
                    position_in_gene = position-start
                else:
                    position_in_gene = end-position

                position_in_gene = position_in_gene-1


                if position_in_gene < 0:
                    continue

                # calculate codon start
                codon_start = int((position_in_gene)/3)*3
                if codon_start+3 > len(gene_sequence):
                    continue


                #codon = gene_sequence[codon_start:codon_start+3]
                codon = oriented_gene_sequence[codon_start:codon_start+3]
                if any(codon_i in codon for codon_i in bases_to_skip):
                    continue
                position_in_codon = position_in_gene%3


                if codon == '':
                    continue

                #if position == 2093:
                #    print(position, codon, strand, start, end, codon_start, position_in_gene, position_in_gene%3)

                effective_gene_synonymous_sites[gene_name] += codon_synonymous_opportunity_table[codon][position_in_codon]/3.0
                effective_gene_nonsynonymous_sites[gene_name] += 1-codon_synonymous_opportunity_table[codon][position_in_codon]/3.0

                for substitution in codon_synonymous_substitution_table[codon][position_in_codon]:
                    substitution_specific_synonymous_sites[substitution] += 1

                for substitution in codon_nonsynonymous_substitution_table[codon][position_in_codon]:
                    substitution_specific_nonsynonymous_sites[substitution] += 1

    substitution_specific_synonymous_fraction = {substitution: substitution_specific_synonymous_sites[substitution]*1.0/(substitution_specific_synonymous_sites[substitution]+substitution_specific_nonsynonymous_sites[substitution]) for substitution in substitution_specific_synonymous_sites.keys()}
    # then annotate promoter regions at remaining sites
    for gene_name,start,end in zip(gene_names,promoter_start_positions,promoter_end_positions):
        for position in range(start,end+1):
            if position not in position_gene_map:
                # position hasn't been annotated yet

                if gene_name not in gene_position_map:
                    # the gene itself has not been annotated
                    # so don't annotate the promoter
                    continue
                
                else:
                    position_gene_map[position] = gene_name
                    gene_position_map[gene_name].append(position)

    # calculate effective gene lengths
    effective_gene_lengths = {gene_name: len(gene_position_map[gene_name])-effective_gene_synonymous_sites[gene_name] for gene_name in gene_position_map.keys()}
    effective_gene_lengths['synonymous'] = sum([effective_gene_synonymous_sites[gene_name] for gene_name in gene_position_map.keys()])
    effective_gene_lengths['nonsynonymous'] = sum([effective_gene_nonsynonymous_sites[gene_name] for gene_name in gene_position_map.keys()])
    effective_gene_lengths['noncoding'] = get_genome_size(reference=reference)-effective_gene_lengths['synonymous']-effective_gene_lengths['nonsynonymous']


    return position_gene_map, effective_gene_lengths, substitution_specific_synonymous_fraction
