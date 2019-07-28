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

import sys, time, os
import pandas as pd
import argparse
from scipy.stats import chi2_contingency

parser = argparse.ArgumentParser(description='differential CHALM methylation level between two conditions')
parser.add_argument("-o", "--out", dest="outfile", metavar="FILE", help="output methylation ratio file name. [default: STDOUT]", default="")
parser.add_argument("-f1", "--file1", dest="file1", metavar="FILE", help="methylation files from the first condition [default: none]", default="")
parser.add_argument("-f2", "--file2", dest="file2", metavar="FILE", help="methylation files from the second condition [default: none]", default="")
args = parser.parse_args()

files1 = args.file1.split(',')
files2 = args.file2.split(',')

def file_loader(files):
    for i in range(0,len(files)):
        if i == 0:
            df_combine = pd.read_table(files[i])
            df_combine.drop(['ratio', 'CI_lower', 'CI_upper'], axis =1, inplace = True)
        else:
            df_file = pd.read_table(files[i])
            df_file.drop(['ratio', 'CI_lower', 'CI_upper'], axis =1, inplace = True)
            df_combine = pd.merge(df_combine, df_file, how='inner', left_on=['chr', 'start', 'end'], right_on=['chr', 'start', 'end'])
            df_combine['total_read_count'] = df_combine.total_read_count_x + df_combine.total_read_count_y
            df_combine['methylated_read_count'] = df_combine.methylated_read_count_x + df_combine.methylated_read_count_y
            df_combine = df_combine[['chr', 'start', 'end', 'total_read_count', 'methylated_read_count']]
    return df_combine

def CHALM_dif(df):
    pvalues=[]
    for i in range(0, len(df)):
        pvalue = chi2_contingency([[(df.total_read_count_x[i] - df.methylated_read_count_x[i]),df.methylated_read_count_x[i]],[(df.total_read_count_y[i] - df.methylated_read_count_y[i]),df.methylated_read_count_y[i]]])[1]
        pvalues.append(pvalue)
    df['CHALM_dif'] = (df.methylated_read_count_y)/df.total_read_count_y - (df.methylated_read_count_x)/df.total_read_count_x
    df['p_value'] = pvalues
    return df

df_condition_1 = file_loader(files1)
df_condition_2 = file_loader(files2)
df_merge = pd.merge(df_condition_1, df_condition_2, how='inner', left_on=['chr', 'start', 'end'], right_on=['chr', 'start', 'end'])
df_output = CHALM_dif(df_merge)
df_output = [['chr', 'start', 'end', 'CHALM_dif', 'p_value']]
df_output.to_csv(args.outfile)









