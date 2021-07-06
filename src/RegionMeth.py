#!/usr/bin/env python2
# -*- coding: utf-8 -*-
import os, re, sys, time, argparse, array, random
import pandas as pd
import numpy as np
from copy import deepcopy

ScriptName=re.split('/', sys.argv[0])[-1]
cmds={"RegionMeth":"Calculate traditional mean methylation ratio or CHALM of given regions in Bed format."};

def disp(text):
    print >> sys.stderr, '[%s] %s' %(time.asctime(), text)

def read_methy_files(ifile, cols=[0,1,2,6,7]):
    names = ['chr', 'pos', 'strand', 'methy', 'total']
    disp('Loading CpG_ratio: %s' % ifile)
    meth_file = pd.read_csv(ifile, sep='\t', header=0, usecols=cols, names=names, compression='infer')
    meth_file.index = meth_file['pos']
    meth_file.drop(['pos'], axis=1, inplace=True)
    return meth_file

def merge_strand_each_chr(df):
    df_p = df[df['strand']=='+']
    df_n = df[df['strand']=='-']
    df_n.index =  df_n.index.values - 1
    merge_index = np.sort(np.unique(np.append(df_n.index.values, df_p.index.values)))
    df_merge = pd.DataFrame(np.zeros(shape=[len(merge_index),2]), index=merge_index)
    df_merge.loc[df_p.index,:] = df_merge.loc[df_p.index,:] + df_p.loc[:,['methy','total']].values
    df_merge.loc[df_n.index,:] = df_merge.loc[df_n.index,:] + df_n.loc[:,['methy','total']].values
    df_merge.columns = ['methy','total']
    df_merge = df_merge.loc[0:,:]
    return df_merge

def merge_strand(df):
    chs = df["chr"].unique().tolist()
    df_merge = pd.DataFrame()
    for ch in chs:
        chr_sub = df[df["chr"] == ch]
        if chr_sub.shape[0] > 0:
            chr_sub = merge_strand_each_chr(chr_sub)
            chr_sub['chr']=pd.Series([ch] * chr_sub.shape[0], index=chr_sub.index)
            df_merge=df_merge.append(chr_sub)
    return df_merge

def Region_Meth_Ratio(methratio_sub,start=0,end=0):
    #methratio_sub: methratio df of one chrom
    region_meth=methratio_sub.loc[start:(end+1),:]
    count_C=region_meth.shape[0]
    methy_C=region_meth["methy"].sum()
    total_C=region_meth["total"].sum()
    if total_C==0:
        region_methratio=np.nan
    else:
        region_methratio=methy_C*1.0/total_C
    return region_methratio

def printHelp():
    print "For help information of each function, try:\n"
    print "  python "+ScriptName+" <Function> -h\n"
    print "Availible Functions:\n"
    for i in cmds.keys():
        print '  '+i+'\t'+cmds[i]+'\n'

def main(cmd=''):
    if cmd not in cmds.keys():
        print "Function "+cmd+" not found, please see the help below:\n"
        printHelp()
        return False

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,\
                                     usage='python '+ScriptName+' '+cmd+' <positional arguments> [optional arguments]\n',\
                                     description=cmds[cmd]+"\n",\
                                     epilog='')
    parser.add_argument(' ',default=None,help='')

    # parse input parameters of each function
    if cmd=="RegionMeth":
        parser.add_argument('Bed',default=None,\
                            help="bed file with at least first 3 columns")
        parser.add_argument('CpG_ratio',default=None,\
                            help="Mean methylation ratio or CHALM of each CpGs generated from BSMAP alignments")
        parser.add_argument('-s', '--usestrand',action="store_true",dest="usestrand",default=False,\
                            help="Use the strand information(6th column) in bed file, and methylation on both strands will not be merged.")
        parser.add_argument('-o', '--output',dest="OUT",metavar='',default='RegionMeth.tsv',\
                            help="Output file name")

    if '-h' in sys.argv or '--help' in sys.argv:
        parser.print_help()
        return False
    else:
        if len(sys.argv)==2:
            print "Too few arguments. Try:\npython "+ScriptName+" "+cmd+" -h"
            return False

    args=parser.parse_args()

    # run function
    if cmd=='RegionMeth':
        disp("RegionMeth Started")
        MethRatio_df=read_methy_files(ifile=args.CpG_ratio, cols=[0,1,2,6,7])
        o1=open(args.OUT,'w')
        if args.usestrand==True:
            Bed=pd.read_csv(args.Bed, sep="\t",usecols=[0,1,2,5],header=None)
            Bed.columns=["chr","start","end","strand"]
            Bed.sort_values(['chr','strand','start','end'], inplace=True, ascending=True)
            disp('Generating methylation level for %d Regions' % Bed.shape[0])
            o1.write('\t'.join(["chr","start","end","strand","methratio"])+'\n')
            chr0="";strand0="";
            for row in Bed.iterrows():
                chr1=row[1]['chr'];
                start1=int(row[1]['start']);
                end1=int(row[1]['end']);
                strand1=row[1]['strand']
                if chr1 != chr0 or strand1!=strand0:
                    MethRatio_sub=MethRatio_df[((MethRatio_df['chr']==chr1) & (MethRatio_df['strand']==strand1))]
                ratio=Region_Meth_Ratio(methratio_sub=MethRatio_sub,start=start1,end=end1)
                aline=[chr1,start1,end1,strand1,ratio];
                o1.write('\t'.join(map(str,aline))+'\n')
                chr0=chr1;strand0=strand1;
        else:
            MethRatio_df=merge_strand(df=MethRatio_df)
            Bed=pd.read_csv(args.Bed, sep="\t",usecols=[0,1,2],header=None)
            Bed.columns=["chr","start","end"]
            Bed.sort_values(['chr','start','end'], inplace=True, ascending=True)
            disp('Generating methylation level for %d Regions' % Bed.shape[0])
            o1.write('\t'.join(["chr","start","end","methratio"])+'\n')
            chr0="";
            for row in Bed.iterrows():
                chr1=row[1]['chr'];
                start1=row[1]['start'];
                end1=row[1]['end'];
                if chr1 != chr0:
                    MethRatio_sub=MethRatio_df[MethRatio_df['chr']==chr1]
                ratio=Region_Meth_Ratio(methratio_sub=MethRatio_sub,start=start1,end=end1)
                aline=[chr1,start1,end1,ratio];
                o1.write('\t'.join(map(str,aline))+'\n')
                chr0=chr1;
        o1.close()
        disp("RegionMeth Finished")

if len(sys.argv)>1:
    main(cmd=sys.argv[1])
else:
    printHelp()
