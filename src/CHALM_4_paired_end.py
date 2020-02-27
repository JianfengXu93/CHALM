#! /usr/bin/env python
import sys, time, os, array, json
import numpy as np
from collections import defaultdict
import readline
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='quantifying methylation level from aligned files')
subparsers = parser.add_subparsers(dest='mode', help='quantification mode')
parser.add_argument('infiles',help="input aligned files") ## positional argument
### traditional mode subparser
sp_trad = subparsers.add_parser('trad', help='traditional mode')
sp_trad.add_argument("-o", "--out", dest="outfile", metavar="FILE", help="output methylation ratio file name. [default: STDOUT]", default="")
sp_trad.add_argument("-O", "--alignment-copy", dest="alignfile", metavar="FILE", help="save a copy of input alignment for BSMAP pipe input. (in BAM format) [default: none]", default="")
sp_trad.add_argument("-w", "--wig", dest="wigfile", metavar="FILE", help="output methylation ratio wiggle file. [default: none]", default="")
sp_trad.add_argument("-b", "--wig-bin", dest="wigbin", type=int, metavar='BIN', help="wiggle file bin size. [default: 25]", default=25)
sp_trad.add_argument("-d", "--ref", dest="reffile", metavar="FILE", help="reference genome fasta file. (required)", default="")
sp_trad.add_argument("-c", "--chr", dest="chroms", metavar="CHR", help="process only specified chromosomes, separated by ','. [default: all]\nexample: --chroms=chr1,chr2", default=[])
sp_trad.add_argument("-s", "--sam-path", dest="sam_path", metavar="PATH", help="path to samtools. [default: none]", default='')
sp_trad.add_argument("-u", "--unique", action="store_true", dest="unique", help="process only unique mappings/pairs.", default=False)                                    
sp_trad.add_argument("-p", "--pair", action="store_true", dest="pair", help="process only properly paired mappings.", default=False)
sp_trad.add_argument("-z", "--zero-meth", action="store_true", dest="meth0", help="report loci with zero methylation ratios. (depreciated, -z will be always enabled)", default=True)
sp_trad.add_argument("-q", "--quiet", action="store_true", dest="quiet", help="don't print progress on stderr.", default=False)
sp_trad.add_argument("-r", "--remove-duplicate", action="store_true", dest="rm_dup", help="remove duplicated reads.", default=False)
sp_trad.add_argument("-t", "--trim-fillin", dest="trim_fillin", type=int, metavar='N', help="trim N end-repairing fill-in nucleotides. [default: 0]", default=0)
sp_trad.add_argument("-g", "--combine-CpG", action="store_true", dest="combine_CpG", help="combine CpG methylaion ratios on both strands.", default=False)
sp_trad.add_argument("-m", "--min-depth", dest="min_depth", type=int, metavar='FOLD', help="report loci with sequencing depth>=FOLD. [default: 1]", default=1)
sp_trad.add_argument("-n", "--no-header", action="store_true", dest="no_header", help="don't print a header line", default=False)
sp_trad.add_argument("-i", "--ct-snp", dest="CT_SNP", help='how to handle CT SNP ("no-action", "correct", "skip"), default: "correct".', default="correct")
sp_trad.add_argument("-x", "--context", dest="context", metavar='TYPE', help="methylation pattern type [CG|CHG|CHH], multiple types separated by ','. [default: all]", default='')
sp_trad.add_argument("-l", "--read-meth", metavar="N", type=int, dest="readmeth", help="calculate clone level methylation for reads containg >= N methyl-C, default: N=1000 (turn off clone-level ratio)", default=1000)
### CHALM mode subparser
sp_CHALM = subparsers.add_parser('CHALM', help='CHALM mode')
sp_CHALM.add_argument("-o", "--out", dest="outfile", metavar="FILE", help="output methylation ratio file name. [default: STDOUT]", default="")
sp_CHALM.add_argument("-O", "--alignment-copy", dest="alignfile", metavar="FILE", help="save a copy of input alignment for BSMAP pipe input. (in BAM format) [default: none]", default="")
sp_CHALM.add_argument("-w", "--wig", dest="wigfile", metavar="FILE", help="output methylation ratio wiggle file. [default: none]", default="")
sp_CHALM.add_argument("-b", "--wig-bin", dest="wigbin", type=int, metavar='BIN', help="wiggle file bin size. [default: 25]", default=25)
sp_CHALM.add_argument("-d", "--ref", dest="reffile", metavar="FILE", help="reference genome fasta file. (required)", default="")
sp_CHALM.add_argument("-c", "--chr", dest="chroms", metavar="CHR", help="process only specified chromosomes, separated by ','. [default: all]\nexample: --chroms=chr1,chr2", default=[])
sp_CHALM.add_argument("-s", "--sam-path", dest="sam_path", metavar="PATH", help="path to samtools. [default: none]", default='')
sp_CHALM.add_argument("-u", "--unique", action="store_true", dest="unique", help="process only unique mappings/pairs.", default=False)                                    
sp_CHALM.add_argument("-p", "--pair", action="store_true", dest="pair", help="process only properly paired mappings.", default=False)
sp_CHALM.add_argument("-z", "--zero-meth", action="store_true", dest="meth0", help="report loci with zero methylation ratios. (depreciated, -z will be always enabled)", default=True)
sp_CHALM.add_argument("-q", "--quiet", action="store_true", dest="quiet", help="don't print progress on stderr.", default=False)
sp_CHALM.add_argument("-r", "--remove-duplicate", action="store_true", dest="rm_dup", help="remove duplicated reads.", default=False)
sp_CHALM.add_argument("-t", "--trim-fillin", dest="trim_fillin", type=int, metavar='N', help="trim N end-repairing fill-in nucleotides. [default: 0]", default=0)
sp_CHALM.add_argument("-g", "--combine-CpG", action="store_true", dest="combine_CpG", help="combine CpG methylaion ratios on both strands.", default=False)
sp_CHALM.add_argument("-n", "--no-header", action="store_true", dest="no_header", help="don't print a header line", default=False)
sp_CHALM.add_argument("-x", "--context", dest="context", metavar='TYPE', help="methylation pattern type [CG|CHG|CHH], multiple types separated by ','. [default: all]", default='')
sp_CHALM.add_argument("-l", "--read-meth", metavar="N", type=int, dest="readmeth", help="calculate clone level methylation for reads containg >= N methyl-C, default: N=1000 (turn off clone-level ratio)", default=1000)
sp_CHALM.add_argument("-R", "--region", dest="region", metavar="FILE", help="bedfile of interested region. (required)", default='')  
sp_CHALM.add_argument("-L", "--read-length", dest="read_length", metavar="N", type=int, help="the length of sequencing read. [default: 100]", default=100) 

args = parser.parse_args()

args.infiles = [args.infiles]
if len(args.reffile) == 0: parser.error("Missing reference file, use -d or --ref option.")
if len(args.infiles) == 0: parser.error("Require at least one BSMAP_MAPPING_FILE.") 
if len(args.chroms) > 0: args.chroms = set(args.chroms.split(','))
if args.trim_fillin < 0: parser.error('Invalid -t value, must >= 0')
seq_context_str, CG, CHG, CHH = ['CG','CHG','CHH'], 1, 2, 3
if len(args.context) > 0: 
    args.context = set(args.context.upper().split(','))
    try: seq_context = set([seq_context_str.index(cc)+1 for cc in args.context])
    except ValueError: parser.error('Invalid -x value, must be one or more of "CG", "CHG", or "CHH"')
else: seq_context = set([1, 2, 3])

if len(args.sam_path) > 0: 
    if args.sam_path[-1] != '/': args.sam_path += '/'

def disp(txt, nt=0):
    if not args.quiet: print >> sys.stderr, '[methratio] @%s \t%s' %(time.asctime(), txt)

if len(args.outfile) == 0: disp("Missing output file name, write to STDOUT.")
def get_alignment(line):
    col = line.split('\t')
    if sam_format:
        if line[0] == '@': return []
        flag = col[1] 
        if 'u' in flag: return []
        if args.unique and 's' in flag: return []
        if args.pair and 'P' not in flag: return []
        rn, cr, pos, cigar ,seq, strand, insert = col[0] ,col[2], int(col[3])-1, col[5], col[9], '', int(col[8]) ## get read name 
        if cr not in args.chroms: return []
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
        if args.unique and flag != 'UM': return []
        if args.pair and col[7] == '0': return []
        seq, strand, cr, pos, insert, mm = col[1], col[6], col[4], int(col[5])-1, int(col[7]), col[9]
        if cr not in args.chroms: return []
        if ':' in mm:
            tmp = mm.split(':')
            gap_pos, gap_size = int(tmp[1]), int(tmp[2])
            if gap_size < 0: seq = seq[:gap_pos] + seq[gap_pos-gap_size:]  # insertion on reference
            else: seq = seq[:gap_pos] + '-' * gap_size + seq[gap_pos:]
    if pos + len(seq) >= len(ref[cr]): return []
    if args.rm_dup:  # remove duplicate hits
        if strand == '+-' or strand == '-+': frag_end, direction = pos+len(seq), 2
        else: frag_end, direction = pos, 1
        if coverage[cr][frag_end] & direction: return []
        coverage[cr][frag_end] |= direction
    if args.trim_fillin > 0: # trim fill in nucleotides
        if strand == '+-' or strand == '-+': seq = seq[:-args.trim_fillin]
        elif strand == '++' or strand == '--': seq, pos = seq[args.trim_fillin:], pos+args.trim_fillin
    if sam_format and insert > 0: seq = seq[:int(col[7])-1-pos] # remove overlapped regions in paired hits, SAM format only
    return (rn, seq, strand[0], cr, pos) ## return read name a

# open pipes to alignment files
pipes = []
for infile in args.infiles:
    nline = 0
    if infile.strip() == '-': sam_format, fin, infile = True, os.popen('%ssamtools view -XSh -' % args.sam_path), 'STDIN'
    elif infile[-4:].upper() == '.SAM': sam_format, fin = True, os.popen('%ssamtools view -XS %s' % (args.sam_path, infile)) 
    elif infile[-4:].upper() == '.BAM': sam_format, fin = True, os.popen('%ssamtools view -X %s' % (args.sam_path, infile))
    else: sam_format, fin = False, open(infile)
    pipes.append((infile,sam_format,fin))

# Read in chromosomes
ref, cr, seq = {}, '', ''
disp('loading reference genome file: %s ...' % args.reffile)
for line in open(args.reffile):
    if line[0] == '>': 
        if len(cr) > 0: 
            if len(args.chroms) == 0 or cr in args.chroms: ref[cr] = seq.upper()
        cr, seq = line[1:-1].split()[0], ''
    else: seq += line.strip()

if len(args.chroms) == 0 or cr in args.chroms: ref[cr] = seq.upper()
del seq


if args.mode == 'trad':
    if args.min_depth <= 0: parser.error('Invalid -m value, must >= 1')
    CT_SNP_val = {"no-action": 0, "correct": 1, "skip": 2}
    try: args.CT_SNP = CT_SNP_val[args.CT_SNP.lower()]
    except KeyError: parser.error('Invalid -i value, select "no-action", "correct" or "skip"')
    meth, depth, coverage, meth1, depth1, refmark = {}, {}, {}, {}, {}, {}
    for cr in ref:
        meth[cr] = array.array('H', [0]) * len(ref[cr])
        depth[cr] = array.array('H', [0]) * len(ref[cr])
        refmark[cr] = array.array('b', [0]) * len(ref[cr])
        if args.rm_dup: coverage[cr] = array.array('B', [0]) * len(ref[cr])
        if args.CT_SNP > 0:
            meth1[cr] = array.array('H', [0]) * len(ref[cr])
            depth1[cr] = array.array('H', [0]) * len(ref[cr])
elif args.mode == 'CHALM':
    meth_read, depth_read, coverage, refmark = {}, {}, {}, {} #CHALM
    for cr in ref:
        refmark[cr] = array.array('b', [0]) * len(ref[cr]) #CHALM
        meth_read[cr] = array.array('H', [0]) * len(ref[cr]) #CHALM
        depth_read[cr] = array.array('H', [0]) * len(ref[cr]) #CHALM
        if args.rm_dup: coverage[cr] = array.array('B', [0]) * len(ref[cr])

args.chroms = set(ref.keys())
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

if args.mode == "trad":
    BS_conversion = {'+': ('C','T','G','A'), '-': ('G','A','C','T')}
    nmap = 0
    for infile, sam_format, fin in pipes:
        disp('reading alignment file: %s ...' % infile)
        nline = 0
        if len(args.alignfile) > 0: pout = os.popen('%ssamtools view -bS - > %s' % (args.sam_path, args.alignfile), 'w')
        map_info_pre=[]
        for line in fin:
            if len(args.alignfile) > 0: pout.write(line)
            nline += 1
            # if nline % 10000000 == 0: disp('read %d lines' % nline, nt=1)
            map_info = get_alignment(line)
            if len(map_info) == 0: 
                map_info_pre = map_info
                continue
            rn, seq, strand, cr, pos = map_info
            if len(map_info_pre) == 0:
                map_info_pre = map_info
                continue
            rn_pre, seq_pre, strand_pre, cr_pre, pos_pre = map_info_pre
            rn_list = rn.split('.')
            rn_pre_list = rn_pre.split('.')
            if rn_list[0] != rn_pre_list[0] or rn_list[1] != rn_pre_list[1]: ## if not paired, then continue
                map_info_pre = map_info
                continue

            ## for current read, find mCpG #
            depthcr = depth[cr]
            pos2 = pos + len(seq)
            nmap += 1
            methcr = meth[cr]
            refseq, refmarkcr = ref[cr], refmark[cr]
            match, convert, rc_match, rc_convert = BS_conversion[strand]
            index, readmeth = refseq.find(match, pos, pos2), 0
            while index >=0 and readmeth < args.readmeth:
                if seq[index-pos] == match and refmarkcr[index] in seq_context: readmeth += 1
                index = refseq.find(match, index+1, pos2)
            ## for previous paired read, find mCpG #
            depthcr_pre = depth[cr_pre]
            pos2_pre = pos_pre + len(seq_pre)
            nmap += 1
            methcr_pre = meth[cr_pre]
            refseq_pre, refmarkcr_pre = ref[cr_pre], refmark[cr_pre]
            match_pre, convert_pre, rc_match_pre, rc_convert_pre = BS_conversion[strand_pre]
            index_pre = refseq_pre.find(match_pre, pos_pre, pos2_pre)
            while index_pre >=0 and readmeth < args.readmeth:
                if seq_pre[index_pre - pos_pre] == match_pre and refmarkcr_pre[index_pre] in seq_context: readmeth += 1
                index_pre = refseq_pre.find(match_pre, index_pre + 1, pos2_pre)

            ## for current read, count mCpG
            index = refseq.find(match, pos, pos2)
            while index >= 0:
                if refmarkcr[index] in seq_context and depthcr[index] < 65535:
                    if seq[index-pos] == convert: 
                        depthcr[index] += 1
                        if readmeth >= args.readmeth: methcr[index] += 1
                    elif seq[index-pos] == match:
                        depthcr[index] += 1
                        methcr[index] += 1
                index = refseq.find(match, index+1, pos2)
            ## for previous paired read, count mCpG
            index_pre = refseq_pre.find(match_pre, pos_pre, pos2_pre)
            while index_pre >= 0:
                if refmarkcr_pre[index_pre] in seq_context and depthcr_pre[index_pre] < 65535:
                    if seq_pre[index_pre - pos_pre] == convert_pre: 
                        depthcr_pre[index_pre] += 1
                        if readmeth >= args.readmeth: methcr_pre[index_pre] += 1
                    elif seq_pre[index_pre - pos_pre] == match_pre:
                        depthcr_pre[index_pre] += 1
                        methcr_pre[index_pre] += 1
                index_pre = refseq_pre.find(match_pre, index_pre + 1, pos2_pre)

            ## for current read
            if args.CT_SNP == 0: 
                map_info_pre = []
                continue
            methcr1 = meth1[cr]
            depthcr1 = depth1[cr]
            index = refseq.find(rc_match, pos, pos2)
            while index >= 0:
                if refmarkcr[index] in seq_context and depthcr1[index] < 65535:
                    if seq[index-pos] == rc_convert: 
                        depthcr1[index] += 1
                    elif seq[index-pos] == rc_match: 
                        depthcr1[index] += 1
                        methcr1[index] += 1
                index = refseq.find(rc_match, index+1, pos2)
            ## for previous paired read
            methcr1_pre = meth1[cr_pre]
            depthcr1_pre = depth1[cr_pre]
            index_pre = refseq_pre.find(rc_match_pre, pos_pre, pos2_pre)
            while index_pre >= 0:
                if refmarkcr_pre[index_pre] in seq_context and depthcr1_pre[index_pre] < 65535:
                    if seq_pre[index_pre - pos_pre] == rc_convert_pre: 
                        depthcr1_pre[index_pre] += 1
                    elif seq_pre[index_pre - pos_pre] == rc_match_pre: 
                        depthcr1_pre[index_pre] += 1
                        methcr1_pre[index_pre] += 1
                index_pre = refseq_pre.find(rc_match_pre, index_pre+1, pos2_pre)

            ## reset map_info_pre
            map_info_pre = []
        fin.close()
        if len(args.alignfile) > 0: pout.close()

    disp('read %d lines' % nline, nt=1)
    if args.combine_CpG:
        disp('combining CpG methylation from both strands ...')
        for cr in depth:
            depthcr, methcr, refcr = depth[cr], meth[cr], ref[cr]
            if args.CT_SNP > 0: depthcr1, methcr1 = depth1[cr], meth1[cr]
            pos = refcr.find('CG')
            while pos >= 0:
                try: 
                    depthcr[pos] += depthcr[pos+1]
                    methcr[pos] += methcr[pos+1]
                except OverflowError: 
                    depthcr[pos] = (depthcr[pos] + depthcr[pos+1]) / 2
                    methcr[pos] = (methcr[pos] + methcr[pos+1]) / 2
                depthcr[pos+1] = 0
                methcr[pos+1] = 0
                if args.CT_SNP > 0:
                    try: 
                        depthcr1[pos] += depthcr1[pos+1]
                        methcr1[pos] += methcr1[pos+1]
                    except OverflowError: 
                        depthcr1[pos] = (depthcr1[pos] + depthcr1[pos+1]) / 2
                        methcr1[pos] = (methcr1[pos] + methcr1[pos+1]) / 2
                pos = refcr.find('CG', pos+2)
    if len(args.outfile) == 0: fout, outfile = sys.stdout, 'STDOUT'
    else: fout = open(args.outfile, 'w')
    disp('writing %s ...' % args.outfile)
    if args.wigfile: 
        fwig = open(args.wigfile, 'w')
        fwig.write('track type=wiggle_0\n')
    if not args.no_header: 
        fout.write('chr\tpos\tstrand\tcontext\tratio\teff_CT_count\tC_count\tCT_count\trev_G_count\trev_GA_count\tCI_lower\tCI_upper\n')
    z95, z95sq = 1.96, 1.96 * 1.96
    nc, nd, dep0 = 0, 0, args.min_depth
    for cr in sorted(depth.keys()):
        depthcr, methcr, refcr, refmarkcr = depth[cr], meth[cr], ref[cr], refmark[cr]
        if args.CT_SNP > 0: depthcr1, methcr1 = depth1[cr], meth1[cr]
        if args.wigfile:
            fwig.write('variableStep chrom=%s span=%d\n' % (cr, args.wigbin))
            bin = wigd = wigm = 0
        for i, dd in enumerate(depthcr):
            if dd < dep0: continue
            if args.CT_SNP > 0: 
                m1, d1 = methcr1[i], depthcr1[i]
                if m1 != d1:
                    if args.CT_SNP == 2: continue
                    d = float(dd) * m1 / d1
                else: d = float(dd)
            else: d = float(dd)
            if refmarkcr[i] not in seq_context: continue
            else: seq = seq_context_str[refmarkcr[i]-1]
            if refcr[i] == 'C': strand = '+'
            else: strand = '-'
            m = methcr[i]
            try: ratio = min(m, d) / d
            except ZeroDivisionError: continue
            nc += 1
            nd += d
            if args.wigfile:
                if i / args.wigbin != bin:
                    if wigd > 0: fwig.write('%d\t%.3f\n' % (bin*args.wigbin+1, min(wigm/wigd,1)))
                    bin = i / args.wigbin
                    wigd = wigm = 0.  
                wigd += d
                wigm += m
            pmid = ratio + z95sq / (2 * d)
            sd = z95 * ((ratio*(1-ratio)/d + z95sq/(4*d*d)) ** 0.5)
            denorminator = 1 + z95sq / d
            CIl, CIu = (pmid - sd) / denorminator, (pmid + sd) / denorminator
            if args.CT_SNP: fout.write('%s\t%d\t%c\t%s\t%.3f\t%.2f\t%d\t%d\t%d\t%d\t%.3f\t%.3f\n' % (cr, i+1, strand, seq, ratio, d, m, dd, m1, d1, CIl, CIu))
            else: fout.write('%s\t%d\t%c\t%s\t%.3f\t%.2f\t%d\t%d\tNA\tNA\t%.3f\t%.3f\n' % (cr, i+1, strand, seq, ratio, d, m, dd, CIl, CIu))                             
    if args.outfile != 'STDOUT': fout.close()
    if args.wigfile: fwig.close()
    disp('total %d valid mappings, %d covered cytosines, average coverage: %.2f fold.' % (nmap, nc, float(nd)/nc))
  
