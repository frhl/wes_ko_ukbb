#!/usr/bin/python

## 21-02-19
## this script returns individual that are heterozygous for the variants

import pandas as pd
import numpy as np
import scipy
import sys
import re
from pysnptools.snpreader import Bed, Pheno
from pandas import DataFrame
from scipy import stats

# check command line args
print 'Number of arguments:', len(sys.argv), 'arguments.'
print 'Argument List:', str(sys.argv)

# args[0] python file
# args[1] input
# args[2] output

# infile
infile=sys.argv[1]
chrom=str(re.findall("\d+",infile)[0])
infile_name= re.search(r"^[^.]*", infile).group(0)

# outfile
outfile_name = sys.argv[2]
outfile = open(outfile_name, "w")

# Read file
bedfile=Bed(infile, count_A1 = True) # A1 = Ref
d=bedfile.read()
variant = d.sid
individual = d.iid
position = d.pos
d_val=np.array(d.val)
maf = np.nanmean(d_val, axis = 0)/2
print d_val.shape

# set homozygous to zero
# and keep heterozygous as one
#d_val[d_val==1]=1
d_val[d_val==2]=0
n_homo_ind = np.nansum(d_val, axis = 1)
n_homo_snp = np.nansum(d_val, axis = 0)
pairs = np.where(d_val==1)
entries = len(pairs[0])

# write some info to outline
info_line1 = '# ' + str(len(n_homo_snp)) + ' SNPS ' + str(len(n_homo_ind)) + ' individuals\n'
info_line2 = '# ' + str(entries) + " heterozygous (for alt allele) individuals.\n"
outfile.write(info_line1)
outfile.write(info_line2)

# get pairs of positions and snps where individuals are homozygous
outfile.write("#SID1\tSID2\tCHR\tPOS\tSNP\tMAF\n")
for i in range(entries):
	cur_var_index = pairs[1][i]
	cur_ind_index = pairs[0][i]
	cur_var = str(variant[cur_var_index].tolist())
	cur_maf = str(maf[cur_var_index].tolist())
	cur_pos_chrom = str(position[cur_var_index].tolist()[0])
	cur_pos_bp = str(position[cur_var_index].tolist()[2])
	cur_ind = individual[cur_ind_index].tolist()
	line = cur_ind
	line.append(cur_pos_chrom)
	line.append(cur_pos_bp)
	line.append(cur_var)
	line.append(cur_maf)
	line.append("\n")
	result = "\t".join(line)
	outfile.write(result)

outfile.close()


