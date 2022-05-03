import numpy
import sys
import bz2
import parse_file

# Load this information about the reference genome so that we know if a snp is in a repeat-masked region
input_filename = sys.argv[1]
depth_filename = sys.argv[2]
snp_filename = sys.argv[3]
strain = sys.argv[4]

input_file = bz2.open(input_filename,"rt")
snp_file = bz2.open(snp_filename,"wt")
depth_file = bz2.open(depth_filename,"wt")

taxon = input_filename.split('/')[-1].split('_')[0][1]


avg_depths = None
times = None
alts = None

depth_records = []


for line in input_file:
    items = line.split(",")
    position = int(items[1])
    allele = items[2].strip()

    if allele[1:3]!='->':
        continue # not a snp!

    snp_file.write(line)

    # calculate depths and add them
    times = numpy.array([float(subitem) for subitem in items[3].split()])
    depths = [float(subitem) for subitem in items[5].split()]

    depth_records.append(depths)


depths = numpy.array(depth_records)
# Could do median or mean
#avg_depths = depths.mean(axis=0)
avg_depths = numpy.median(depths, axis=0)

alts = numpy.array([0 for t in times])

depth_line = ", ".join([strain, "0", "Depth", " ".join([str(t) for t in times]), " ".join([str(alt) for alt in alts]), " ".join([str(avg_depth) for avg_depth in avg_depths])])

depth_file.write(depth_line)
depth_file.write("\n")

input_file.close()
snp_file.close()
depth_file.close()
