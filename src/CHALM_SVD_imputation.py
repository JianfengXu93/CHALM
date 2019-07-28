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

import sys, time, os, array, optparse, json
import numpy as np
from collections import defaultdict
import readline
import rpy2 
import rpy2.robjects
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
pandas2ri.activate()
import pandas as pd
bcv = importr('bcv')
NA = rpy2.rinterface.NA_Real
usage = "usage: %prog [options] BSMAP_MAPPING_FILES"
parser = optparse.OptionParser(usage=usage)

parser.add_option("-o", "--out", dest="outfile", metavar="FILE", help="output methylation ratio file name. [default: STDOUT]", default="")
parser.add_option("-O", "--alignment-copy", dest="alignfile", metavar="FILE", help="save a copy of input alignment for BSMAP pipe input. (in BAM format) [default: none]", default="")
parser.add_option("-b", "--wig-bin", dest="wigbin", type="int", metavar='BIN', help="wiggle file bin size. [default: 25]", default=25)
parser.add_option("-d", "--ref", dest="reffile", metavar="FILE", help="reference genome fasta file. (required)", default="")
parser.add_option("-c", "--chr", dest="chroms", metavar="CHR", help="process only specified chromosomes, separated by ','. [default: all]\nexample: --chroms=chr1,chr2", default=[])
parser.add_option("-s", "--sam-path", dest="sam_path", metavar="PATH", help="path to samtools. [default: none]", default='')
parser.add_option("-u", "--unique", action="store_true", dest="unique", help="process only unique mappings/pairs.", default=False)                                    
parser.add_option("-p", "--pair", action="store_true", dest="pair", help="process only properly paired mappings.", default=False)
parser.add_option("-z", "--zero-meth", action="store_true", dest="meth0", help="report loci with zero methylation ratios. (depreciated, -z will be always enabled)", default=True)
parser.add_option("-q", "--quiet", action="store_true", dest="quiet", help="don't print progress on stderr.", default=False)
parser.add_option("-r", "--remove-duplicate", action="store_true", dest="rm_dup", help="remove duplicated reads.", default=False)
parser.add_option("-t", "--trim-fillin", dest="trim_fillin", type="int", metavar='N', help="trim N end-repairing fill-in nucleotides. [default: 0]", default=0)
parser.add_option("-n", "--no-header", action="store_true", dest="no_header", help="don't print a header line", default=False)
parser.add_option("-x", "--context", dest="context", metavar='TYPE', help="methylation pattern type [CG|CHG|CHH], multiple types separated by ','. [default: all]", default='')
parser.add_option("-l", "--read-meth", metavar="N", type="int", dest="readmeth", help="calculate clone level methylation for reads containg >= N methyl-C. [default: 1]", default=1)
parser.add_option("-R", "--region", dest="region", metavar="FILE", help="bedfile of interested region. (required)", default='')
parser.add_option("-L", "--read-length", dest="read_length", metavar="N", type="int", help="the length of read to define the border of region. [default: 0]", default=0)
parser.add_option("-e", "--read-extend", metavar="N", type="int", dest="extend", help="the length to extend the raw read. [default: 0]", default=0)


options, infiles = parser.parse_args()


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
options.extend = options.extend/2
if len(options.sam_path) > 0: 
    if options.sam_path[-1] != '/': options.sam_path += '/'

def disp(txt, nt=0):
    if not options.quiet: print >> sys.stderr, '[methratio] @%s \t%s' %(time.asctime(), txt)

if len(options.outfile) == 0: disp("Missing output file name, write to STDOUT.")
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
    return (seq, strand[0], cr, pos)

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

meth_read, depth_read, coverage, refmark, readseq, readseq_pos, readseq_neg = {}, {}, {}, {}, {}, defaultdict(), defaultdict()
for cr in ref:
    refmark[cr] = array.array('b', [0]) * len(ref[cr])
    readseq_pos[cr] = defaultdict()
    readseq_neg[cr] = defaultdict()
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

BS_conversion, readseq = {'+': ('C','T','G','A'), '-': ('G','A','C','T')}, {'+': readseq_pos, '-': readseq_neg}
nmap = 0
for infile, sam_format, fin in pipes:
    disp('reading alignment file: %s ...' % infile)
    nline = 0
    if len(options.alignfile) > 0: pout = os.popen('%ssamtools view -bS - > %s' % (options.sam_path, options.alignfile), 'w')
    for line in fin:
        if len(options.alignfile) > 0: pout.write(line)
        nline += 1
        if nline % 10000000 == 0: disp('read %d lines' % nline, nt=1)
        map_info = get_alignment(line)
        if len(map_info) == 0: continue
        seq, strand, cr, pos = map_info         
        nmap += 1
        refseq, refmarkcr, readseq_cr = ref[cr], refmark[cr], readseq[strand][cr]
        match, convert, rc_match, rc_convert = BS_conversion[strand]
        try: readseq_cr[pos].append(seq)
        except KeyError: readseq_cr[pos] = [seq]


    fin.close()
    if len(options.alignfile) > 0: pout.close()

disp('read %d lines' % nline, nt=1)

def CpG_map(cr, pos, strand, reads, CGI_start, CGI_end, extend):
    refseq, refmarkcr, mapped_read, i, read_p, end_p = ref[cr], refmark[cr], [], 0, [], []
    match, convert, rc_match, rc_convert = BS_conversion[strand]
    pos = pos - extend
    for read in reads:
        read = 'N'*extend + read + 'N'*extend
    	if len(read) <=2 : continue
    	d=0
        found, last = False, False
        mapped_read.append([])
        index = refseq.find(match, CGI_start, CGI_end)
        while index >= 0:
            if (index-pos >= len(read)  or index-pos < 0) and refmarkcr[index] in seq_context: 
                if index-pos >= len(read) and last == False and found == True:
                    end_p.append(d)
                    last = True
                mapped_read[i].append(NA)
                d +=1
            elif index-pos >= len(read)  or index-pos < 0:
            	pass
            elif read[index-pos] == match and refmarkcr[index] in seq_context: 
                mapped_read[i].append(1)
                if found == False: 
                    read_p.append(d)
                    found = True
                d +=1
            elif read[index-pos] == convert and refmarkcr[index] in seq_context: 
                mapped_read[i].append(-1)
                if found == False: 
                    read_p.append(d)
                    found = True
                d +=1
            elif refmarkcr[index] in seq_context:
                mapped_read[i].append(NA)
                if found == False:
                    read_p.append(d)
                    found = True
                d +=1
            index = refseq.find(match, index+1, CGI_end)
        if found == False: 
            del mapped_read[i]
            found = 'deleted'
            i -= 1
        elif last == False: end_p.append(len(mapped_read[i]))
        if found == 'deleted':
        	i += 1
        	continue
        if ((1 in mapped_read[i]) == False) and ((-1 in mapped_read[i]) == False) and ((NA in mapped_read[i]) == True):
            del mapped_read[i]
            del read_p[i]
            del end_p[i]
            i -= 1

        i += 1
    return (mapped_read, read_p, end_p)


z95, z95sq = 1.96, 1.96 * 1.96
def count_mc_imp(read):
    return (1 - reduce(lambda x,y: x*y, read))
def count_mc(read):
    if 0 in read: return 1
    else: return 0

fout = open(options.outfile, 'w')
fout.write('chr\tCGI_start\tCGI_end\ttotal_count\tmeth_count\tratio\tmeth_count_imp\tratio_imp\n')
fout.close()
for line in open(options.region):
    read_p_pos, read_p_neg, end_p_pos, end_p_neg = [], [], [], []
    cr, start, end, CGI_start, CGI_end = line.split('\t')[0], int((int(line.split('\t')[1])-options.read_length*0.5)), int((int(line.split('\t')[2])-options.read_length*0.5)), int(line.split('\t')[1]), int(line.split('\t')[2])
    meth_matrix_pos, meth_matrix_neg , readseq_pos_cr, readseq_neg_cr = [], [], readseq['+'][cr], readseq['-'][cr]
    for i in range(start, end+1):
        try: mapped_read_pos, p_pos, end_pos = CpG_map(cr, i, '+', readseq_pos_cr[i], CGI_start, CGI_end, options.extend)
        except KeyError: mapped_read_pos = []
        try: mapped_read_neg, p_neg, end_neg = CpG_map(cr, i, '-', readseq_neg_cr[i], CGI_start, CGI_end, options.extend)
        except KeyError: mapped_read_neg = []
        if mapped_read_pos != []: 
            meth_matrix_pos.extend(mapped_read_pos)
            read_p_pos.extend(p_pos)
            end_p_pos.extend(end_pos)
        if mapped_read_neg != []: 
            meth_matrix_neg.extend(mapped_read_neg)
            read_p_neg.extend(p_neg)
            end_p_neg.extend(end_neg)
    ## do SVD to the meth_matrix
    if len(meth_matrix_pos) <=20 or len(meth_matrix_neg) <=20: continue 
    df_pos, df_neg = pd.DataFrame(meth_matrix_pos), pd.DataFrame(meth_matrix_neg)
    k_pos = min([50,len(df_pos.index), len(df_pos.columns)])
    k_neg = min([50,len(df_neg.index), len(df_neg.columns)])
    df_pos, df_neg = pandas2ri.py2ri(df_pos), pandas2ri.py2ri(df_neg)
    try: df_pos_svd = bcv.impute_svd(df_pos,k=k_pos,maxiter=10)
    except rpy2.rinterface.RRuntimeError: continue
    try: df_neg_svd = bcv.impute_svd(df_neg,k=k_neg,maxiter=10)
    except rpy2.rinterface.RRuntimeError: continue
    mx_pos_filled, mx_neg_filled = np.matrix(df_pos_svd[0]), np.matrix(df_neg_svd[0])
    mx_pos_filled[mx_pos_filled > 1] = 1
    mx_neg_filled[mx_neg_filled < -1] = -1
    mx_pos_filled, mx_neg_filled = 1 - (mx_pos_filled*0.5 + 0.5), 1 - (mx_neg_filled*0.5 + 0.5)
    mx_pos_filled, mx_neg_filled = mx_pos_filled.tolist(), mx_neg_filled.tolist()
    extend_read_pos, extend_read_neg = [], []
    for i in range(0, len(mx_pos_filled)):
        read_start = max(read_p_pos[i], 0)
        read_end = min(end_p_pos[i], len(mx_pos_filled[i]))
        extend_read_pos.append(mx_pos_filled[i][read_start:read_end])
    for i in range(0, len(mx_neg_filled)):
        read_start = max(read_p_neg[i], 0)
        read_end = min(end_p_neg[i], len(mx_neg_filled[i]))
        extend_read_neg.append(mx_neg_filled[i][read_start:read_end])
    extend_read_pos, extend_read_neg = [np.matrix(x) for x in extend_read_pos], [np.matrix(x) for x in extend_read_neg]    
    totalC = len(extend_read_pos) + len(extend_read_neg)
    methC_imp = sum([float(np.apply_along_axis(count_mc_imp, axis=1, arr=x.tolist())) for x in extend_read_pos]) + sum([float(np.apply_along_axis(count_mc_imp, axis=1, arr=x.tolist())) for x in extend_read_neg])
    methC = sum([float(np.apply_along_axis(count_mc, axis=1, arr=x)) for x in extend_read_pos]) + sum([float(np.apply_along_axis(count_mc, axis=1, arr=x)) for x in extend_read_neg])
    ratio, ratio_imp = methC/float(totalC), methC_imp/float(totalC)
    fout = open(options.outfile, 'a')
    fout.write('%s\t%d\t%d\t%d\t%d\t%.3f\t%.3f\t%.3f\n' % (cr, CGI_start, CGI_end, totalC, methC, ratio, methC_imp, ratio_imp))
    fout.close()

