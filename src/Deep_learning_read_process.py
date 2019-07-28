 #! /usr/bin/env python
# Copyright 2018 Terry Jones, Wei Li

# This file is part of CHALM.

# CHALM is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# CHALM is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with CHALM.  If not, see <https://www.gnu.org/licenses/>.

import sys, time, os, array, optparse, random
import pandas as pd
from collections import defaultdict
import numpy as np
from operator import itemgetter
from optparse import OptionParser

usage = "usage: %prog [options] BSMAP_MAPPING_FILES"
parser = optparse.OptionParser(usage=usage)

parser.add_option("-o", "--output-dir", dest="out_dir", metavar="FILE", help="output methylation ratio file name. [default: STDOUT]", default="")
parser.add_option("-n", "--sample-name", dest="name", metavar="FILE", help="sample name for output files", default="Sample")
parser.add_option("-d", "--ref", dest="reffile", metavar="FILE", help="reference genome fasta file. (required)", default="")
parser.add_option("-c", "--chr", dest="chroms", metavar="CHR", help="process only specified chromosomes, separated by ','. [default: all]\nexample: --chroms=chr1,chr2", default=[])
parser.add_option("-s", "--sam-path", dest="sam_path", metavar="PATH", help="path to samtools. [default: none]", default='')
parser.add_option("-u", "--unique", action="store_true", dest="unique", help="process only unique mappings/pairs.", default=False)                                    
parser.add_option("-p", "--pair", action="store_true", dest="pair", help="process only properly paired mappings.", default=False)
parser.add_option("-q", "--quiet", action="store_true", dest="quiet", help="don't print progress on stderr.", default=False)
parser.add_option("-r", "--remove-duplicate", action="store_true", dest="rm_dup", help="remove duplicated reads.", default=False)
parser.add_option("-t", "--trim-fillin", dest="trim_fillin", type="int", metavar='N', help="trim N end-repairing fill-in nucleotides. [default: 0]", default=0)
parser.add_option("-x", "--context", dest="context", metavar='TYPE', help="methylation pattern type [CG|CHG|CHH], multiple types separated by ','. [default: all]", default='')
parser.add_option("-S", "--shuffle", action="store_true", dest="shuffle", help="shuffle the mCpG to destory clonal information", default=False) 
parser.add_option("--region", action='store', type = 'string', dest="region", metavar="FILE", help="the modified bed file of interested regions", default="") 
parser.add_option("--depth_cut", dest="depth_cut", type="int", metavar='N', help="mininal number of reads mapped to the interested region. [default: 50]", default=50)
parser.add_option("--read_bins", dest="read_bins", type="int", metavar='N', help="number of bins for containing mapped reads. [default: 200]", default=200)

options, infiles = parser.parse_args()

##### label read for methylation level: mCpG==3, CpG==2, non-CpG==1
if len(options.reffile) == 0: parser.error("Missing reference file, use -d or --ref option.")
if len(infiles) == 0: parser.error("Require at least one BSMAP_MAPPING_FILE.") 
if len(options.chroms) > 0: options.chroms = set(options.chroms.split(','))
if options.trim_fillin < 0: parser.error('Invalid -t value, must >= 0')
seq_context_str, CG, CHG, CHH = ['CG','CHG','CHH'], 1, 2, 3
if len(options.context) > 0: 
    options.context = set(options.context.upper().split(','))
    try: seq_context = set([seq_context_str.index(cc)+1 for cc in options.context])
    except ValueError: parser.error('Invalid -x value, must be one or more of "CG", "CHG", or "CHH"')
else: seq_context = set([1, 2, 3])

if len(options.sam_path) > 0: 
    if options.sam_path[-1] != '/': options.sam_path += '/'
if len(options.out_dir) > 0: 
    if options.out_dir[-1] != '/': options.out_dir += '/'

def disp(txt, nt=0):
    if not options.quiet: print >> sys.stderr, '[methratio] @%s \t%s' %(time.asctime(), txt)

def get_alignment(line):
    col = line.split('\t')
    if sam_format:
        if line[0] == '@': return []
        flag = col[1] 
        if 'u' in flag: return []
        if options.unique and 's' in flag: return []
        if options.pair and 'P' not in flag: return []
        cr, pos, cigar, seq, strand, insert = col[2], int(col[3])-1, col[5], col[9], '', int(col[8])
        if cr not in options.chroms: return []
        strand_index = line.find('ZS:Z:')
        assert strand_index >= 0, 'missing strand information "ZS:Z:xx"'
        strand = line[strand_index+5:strand_index+7]
        gap_pos, gap_size = 0, 0
        while 'I' in cigar or 'D' in cigar:
            for sep in 'MID':
                try: gap_size = int(cigar.split(sep, 1)[0])
                except ValueError: continue
                break
            if sep == 'M': gap_pos += gap_size
            elif sep == 'I': seq = seq[:gap_pos] + seq[gap_pos+gap_size:]
            elif sep == 'D': 
                seq = seq[:gap_pos] + '-' * gap_size + seq[gap_pos:]                        
                gap_pos += gap_size
            cigar = cigar[cigar.index(sep)+1:]
    else:
        flag = col[3][:2]
        if flag == 'NM' or flag == 'QC': return []
        if options.unique and flag != 'UM': return []
        if options.pair and col[7] == '0': return []
        seq, strand, cr, pos, insert, mm = col[1], col[6], col[4], int(col[5])-1, int(col[7]), col[9]
        if cr not in options.chroms: return []
        if ':' in mm:
            tmp = mm.split(':')
            gap_pos, gap_size = int(tmp[1]), int(tmp[2])
            if gap_size < 0: seq = seq[:gap_pos] + seq[gap_pos-gap_size:]  # insertion on reference
            else: seq = seq[:gap_pos] + '-' * gap_size + seq[gap_pos:]
    if pos + len(seq) >= len(ref[cr]): return []
    if options.rm_dup:  # remove duplicate hits
        if strand == '+-' or strand == '-+': frag_end, direction = pos+len(seq), 2
        else: frag_end, direction = pos, 1
        if coverage[cr][frag_end] & direction: return []
        coverage[cr][frag_end] |= direction
    if options.trim_fillin > 0: # trim fill in nucleotides
        if strand == '+-' or strand == '-+': seq = seq[:-options.trim_fillin]
        elif strand == '++' or strand == '--': seq, pos = seq[options.trim_fillin:], pos+options.trim_fillin
    if sam_format and insert > 0: seq = seq[:int(col[7])-1-pos] # remove overlapped regions in paired hits, SAM format only
    return (col, seq, strand[0], cr, pos)   

# open pipes to alignment files
pipes = []
for infile in infiles:
    nline = 0
    if infile.strip() == '-': sam_format, fin, infile = True, os.popen('%ssamtools view -XSh -' % options.sam_path), 'STDIN'
    elif infile[-4:].upper() == '.SAM': sam_format, fin = True, os.popen('%ssamtools view -XS %s' % (options.sam_path, infile)) 
    elif infile[-4:].upper() == '.BAM': sam_format, fin = True, os.popen('%ssamtools view -X %s' % (options.sam_path, infile))
    else: sam_format, fin = False, open(infile)
    pipes.append((infile,sam_format,fin))

# Read in chromosomes
ref, cr, seq = {}, '', ''
disp('loading reference genome file: %s ...' % options.reffile)
for line in open(options.reffile):
    if line[0] == '>': 
        if len(cr) > 0: 
            if len(options.chroms) == 0 or cr in options.chroms: ref[cr] = seq.upper()
        cr, seq = line[1:-1].split()[0], ''
    else: seq += line.strip()

if len(options.chroms) == 0 or cr in options.chroms: ref[cr] = seq.upper()
del seq

coverage, refmark = {}, {}
for cr in ref:
    refmark[cr] = array.array('b', [0]) * len(ref[cr])
    if options.rm_dup: coverage[cr] = array.array('B', [0]) * len(ref[cr])

options.chroms = set(ref.keys())

disp('marking reference genome ...')
for cr in ref:
    refcr, refmarkcr = ref[cr], refmark[cr]
    index = refcr.find('C', 0, len(refcr)-2)
    while index >= 0:
        if refcr[index+1] == 'G': refmarkcr[index] = CG
        elif refcr[index+2] == 'G': refmarkcr[index] = CHG
        else: refmarkcr[index] = CHH
        index = refcr.find('C', index+1, len(refcr)-2)
    index = refcr.find('G', 2, len(refcr))
    while index >= 0:
        if refcr[index-1] == 'C': refmarkcr[index] = CG
        elif refcr[index-2] == 'C': refmarkcr[index] = CHG
        else: refmarkcr[index] = CHH
        index = refcr.find('G', index+1, len(refcr))

BS_conversion = {'+': ('C','T','G','A'), '-': ('G','A','C','T')}

fout = open(options.out_dir+options.name+'_methlabel.txt', 'w')
disp('writing %s ...' % (options.out_dir+options.name+'_methlabel.txt'))

for infile, sam_format, fin in pipes:
    disp('reading alignment file: %s ...' % infile)
    nline = 0
    for line in fin:
        nline += 1
        if nline % 10000000 == 0: disp('read %d lines' % nline, nt=1)
        map_info = get_alignment(line)
        if len(map_info) <= 1: continue
        ori_line, seq, strand, cr, pos = map_info  
        if len(seq) == 0: continue
        fout.write('%s\t' % ori_line[2])
        fout.write('%s\t' % ori_line[3])
        fout.write('%d\t' % (pos + len(seq)))
        refseq, refmarkcr = ref[cr], refmark[cr]
        match, convert, rc_match, rc_convert = BS_conversion[strand]
        for i in range(0, len(seq)-1):
            index = pos + i
            if refmarkcr[index] in seq_context:
                if seq[i] == match:
                    fout.write('%s' % '3')  ##mCpG == 3
                elif seq[i] == convert:
                    fout.write('%s' % '2')  ##CpG == 2
                else:
                    fout.write('%s' % '1')  ##non-C == 1
            else:
                fout.write('%s' % '1')      ##non-C == 1
        if refmarkcr[index+1] in seq_context:
            if seq[-1] == match:
                fout.write('%s\n' % '3')
            elif seq[-1] == convert:
                fout.write('%s\n' % '2')
            else:
                fout.write('%s\n' % '1')
        else: 
            fout.write('%s\n' % '1')
                           
fout.close()

##### converting methylation label into methylation 2D code
def read_process(TSS, read_start, read_end, strand, read_seq, total_C):
    out_seq=[]
    if strand == '+':
        tss_dis = TSS - (read_end + read_start)/2  ## upstream will be assigned a positive value
    else:
        tss_dis = (read_end + read_start)/2 - TSS
    if tss_dis > 0:
        tss_dis = np.log2(tss_dis)
        tss_dis = float("{0:.1f}".format(tss_dis))
    elif tss_dis <0:
        tss_dis = np.log2(-tss_dis)
        tss_dis = -float("{0:.1f}".format(tss_dis))
    read_seq = read_seq.replace('1', '')
    if total_C < 10:
        for i in range(0, total_C):
            out_seq.append(int(read_seq[i])-2)
            if (i < 10 - total_C):
                out_seq.append(int(read_seq[i])-2)
    else:
        size_factor = total_C/10
        remainder = total_C - 10*size_factor
        for i in range(0,10):
            if(i < remainder):
                read_fraction = read_seq[i*(size_factor + 1):(i + 1)*(size_factor + 1)]
                n_mC = read_fraction.count('3')
                out_seq.append(float(n_mC)/(size_factor + 1))
            else:
                read_fraction = read_seq[(remainder*(size_factor+1)+(i-remainder)*size_factor):(remainder*(size_factor+1)+(i-remainder+1)*size_factor)]
                n_mC = read_fraction.count('3')
                out_seq.append(float(n_mC)/size_factor)
    total_meth = sum(out_seq)
    out_seq.append(total_meth)
    out_seq.append(tss_dis)
    return out_seq

def read_shuffle(matrix_fraction):
    matrix_fraction = np.array(matrix_fraction)
    n_read = len(matrix_fraction)
    read_pos = list(matrix_fraction[:,11])
    list_mC = []
    for i in range(0,n_read):
        list_mC += list(matrix_fraction[i][0:10])
    np.random.shuffle(list_mC)
    list_mC = np.reshape(list_mC, (n_read,10))
    shuffled_fraction=[]
    for i in range(0,n_read):
        element=list(list_mC[i]) + [sum(list_mC[i])] + [read_pos[i]]
        shuffled_fraction.append(element)
    return shuffled_fraction

os.system('sort -k1,1 -k2,2n %s > %s' %(options.out_dir + options.name + '_methlabel.txt', options.out_dir + options.name + '_methlabel_sorted.txt'))
os.system('bedtools intersect -wa -wb -a %s -b %s > %s' % (options.region, options.out_dir + options.name + '_methlabel_sorted.txt', options.out_dir + options.name + '_mapped_reads.txt'))

fin = open(options.out_dir + options.name + '_mapped_reads.txt', 'r')
dict_promoter = defaultdict()
for line in fin:
    flag = 0
    try: dict_promoter[line.split('\t')[1]]
    except KeyError: flag = 1
    read_seq = line.split('\t')[9]
    total_C = read_seq.count('2') + read_seq.count('3') 
    if total_C < 5:
        continue
    if flag == 1:
        dict_promoter[line.split('\t')[1]] = [line.split('\t')[0], line.split('\t')[1], line.split('\t')[2], line.split('\t')[3], line.split('\t')[4],line.split('\t')[5],[]]
        read_meth = read_process(int(line.split('\t')[3]), int(line.split('\t')[7]), int(line.split('\t')[8]), line.split('\t')[5], line.split('\t')[9], total_C)
        (dict_promoter[line.split('\t')[1]])[6].append(read_meth)

    else:
        read_meth = read_process(int(line.split('\t')[3]), int(line.split('\t')[7]), int(line.split('\t')[8]), line.split('\t')[5], line.split('\t')[9], total_C)
        (dict_promoter[line.split('\t')[1]])[6].append(read_meth)

## downsampling and sort 
dict_dp = defaultdict()
depth_cut = options.depth_cut #default: 50
read_bins = options.read_bins #default: 200
for i in dict_promoter.keys():
    if len(dict_promoter[i][6]) < depth_cut: continue
    ## upsampling
    if len(dict_promoter[i][6]) < read_bins and len(dict_promoter[i][6]) >= depth_cut:
        list_element = [dict_promoter[i][0],dict_promoter[i][1],dict_promoter[i][2],dict_promoter[i][3],dict_promoter[i][4],dict_promoter[i][5]]
        fraction = dict_promoter[i][6]
        if options.shuffle: fraction=read_shuffle(fraction)
        fraction_size = int(read_bins/len(dict_promoter[i][6]))
        fraction_extended = fraction*fraction_size
        sample_n = read_bins - len(fraction_extended)
        fraction_upsampled = fraction_extended + random.sample(fraction, sample_n)
        fraction_upsampled_sorted = sorted(fraction_upsampled, key=itemgetter(10)) ## sort the matrix by the total meth
        list_element.append(fraction_upsampled_sorted)
        dict_dp[i] = list_element
    ## devide by size factor and then downsampling
    if len(dict_promoter[i][6]) >= read_bins:
        list_element = [dict_promoter[i][0],dict_promoter[i][1],dict_promoter[i][2],dict_promoter[i][3],dict_promoter[i][4],dict_promoter[i][5]]
        fraction = dict_promoter[i][6]
        if options.shuffle: fraction=read_shuffle(fraction)
        fraction_size = int(len(dict_promoter[i][6]))/read_bins
        fraction_downsampled = random.sample(fraction, read_bins*fraction_size)
        fraction_downsampled_sorted = sorted(fraction_downsampled, key=itemgetter(10))
        fraction_new = []
        for k in range(0,read_bins):
            fraction_slice = fraction_downsampled_sorted[k*fraction_size:(k+1)*fraction_size]
            fraction_slice_mean = [float(sum(l))/len(l) for l in zip(*fraction_slice)]
            fraction_new.append(fraction_slice_mean)
        list_element.append(fraction_new)
        dict_dp[i] = list_element
## output methylation and position files
keys_intersect = list(set(dict_dp.keys()))
array_out_meth = []
array_out_pos = []
for i in keys_intersect:
    gene_ID, X_vector, X_vector_meth, X_vector_pos = [dict_dp[i][4]], dict_dp[i][6], [], [] 
    for i in range(0,read_bins):
        X_vector_meth += X_vector[i][0:10]
        X_vector_pos += X_vector[i][11:12]*10
    array_element_meth = gene_ID + X_vector_meth
    array_element_pos = gene_ID + X_vector_pos
    array_out_meth.append(array_element_meth)
    array_out_pos.append(array_element_pos)  
df_output_meth = pd.DataFrame(array_out_meth)
df_output_pos = pd.DataFrame(array_out_pos)
if options.shuffle:
    df_output_meth.to_csv(options.out_dir + options.name + '_meth_2D_code_control.txt', sep='\t', index=False, header=False)
    df_output_pos.to_csv(options.out_dir + options.name + '_distance_2_TSS_control.txt', sep='\t', index=False, header=False)
else:
    df_output_meth.to_csv(options.out_dir + options.name + '_meth_2D_code.txt', sep='\t', index=False, header=False)
    df_output_pos.to_csv(options.out_dir + options.name + '_distance_2_TSS.txt', sep='\t', index=False, header=False)
















