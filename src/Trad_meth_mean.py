#! /usr/bin/env python

import pandas as pd
from collections import defaultdict

## for unix use
from optparse import OptionParser
usage = "usage: %prog [options] arg"
parser = OptionParser(usage=usage)
parser.add_option("-r", "--region", action='store', type = 'string', dest="region", metavar="FILE", help="input analyzed region bedfile", default="")
parser.add_option("-m", "--dnameth", action='store', type = 'string', dest="DNAmeth", metavar="FILE", help="input DNA methylation file", default="")
parser.add_option("-o", "--output", action='store', type = 'string', dest="output", metavar="FILE", help="output", default="")

(options, args) = parser.parse_args()

## input CGI bedfile
df_region = pd.read_table(options.region)
df_region = df_region.sort_values(by = 'start', ascending =True)
#df_region.drop(['#bin','cpgNum','gcNum','perCpg','perGc','obsExp'], axis=1, inplace=True) ##drop multiple columns

## input meth ratio file
df_meth = pd.read_table(options.DNAmeth)
#df_meth.drop(['context','eff_CT_count', 'C_count', 'CT_count', 'rev_G_count', 'rev_GA_count', 'CI_lower', 'CI_upper'], axis =1, inplace = True)

## label for each chromosome
list_chr = list(df_region.chr.unique())
dict_chr = defaultdict()
for item in list_chr:
    dict_chr[item]=0
    
## functions to calculate CGI methyaltion
def CGI_meth(CGI_chrom, CGI_start, CGI_end, strand, Gene_symbol):
    global dict_chr
    df_meth_chr = df_meth[df_meth['chr'] == CGI_chrom]
    df_meth_chr = df_meth_chr.sort_values(by='pos', ascending = True)
    index_CpG = len(df_meth_chr.index)
    sums, array_meth = [], []
    for i in range(dict_chr[CGI_chrom], index_CpG):
        if(df_meth_chr.iloc[i]['pos'] >= CGI_start) and (df_meth_chr.iloc[i]['pos'] <= CGI_end):
            sums.append(df_meth_chr.iloc[i]['ratio'])
            if len(sums) == 1:
                dict_chr[CGI_chrom] = i
        if(df_meth_chr.iloc[i]['pos'] > CGI_end):
            break 
    series_sums = pd.Series(sums)
    average = series_sums.mean()
    if len(series_sums.index) ==0:
        array_meth.append(CGI_chrom)
        array_meth.append(CGI_start)
        array_meth.append(CGI_end)
        array_meth.append(strand)
        array_meth.append(Gene_symbol)
        array_meth.append(-1)
        array_meth.append(len(series_sums.index))
    else:
        array_meth.append(CGI_chrom)
        array_meth.append(CGI_start)
        array_meth.append(CGI_end)
        array_meth.append(strand)
        array_meth.append(Gene_symbol)
        array_meth.append(average)
        array_meth.append(len(series_sums.index))
    return array_meth

index_GE = len(df_region.index)
array_GE_meth = []

for i in range(0, index_GE):
    array_GE_meth.append(CGI_meth(df_region.iloc[i]['chr'], df_region.iloc[i]['start'], df_region.iloc[i]['end'], df_region.iloc[i]['strand'], df_region.iloc[i]['Gene_symbol']))


df_GE_meth = pd.DataFrame(array_GE_meth, columns = ['chr','start', 'end', 'strand', 'Gene_symbol', 'meth_mean', 'number_CpG'])
df_GE_meth.to_csv(options.output, sep='\t', index = False)



    
    
