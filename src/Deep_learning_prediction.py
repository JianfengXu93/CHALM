#! /usr/bin/env python

import torch
import numpy as np
import os
import random
import time
from torch.utils.data import TensorDataset
from torch.autograd import Variable
import pandas as pd
import argparse
from collections import defaultdict
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description='quantifying methylation level from aligned files')
## positional argument
## optional argument
parser.add_argument("-o", "--output-dir", dest="out_dir", metavar="FILE", help="output methylation ratio file name. [default: STDOUT]", default="")
parser.add_argument("-m", "--meth", dest="trad_meth", metavar="FILE", help="traditional mean methylation level", default="")
parser.add_argument("-e", "--exp", dest="expression", metavar="FILE", help="gene expression file for prediction", default="")
parser.add_argument("-s", "--sample-name", dest="sample", metavar="SAMPLE", help="", default="input_sample")
parser.add_argument("-d", "--down-sample", dest="down_sample", action="store_true", default=False)
parser.add_argument("-f1", "--2D-code", dest="meth_code", metavar="FILE", help="input methylation 2D code of regions of interest, output file from Deep_learning_read_rpocess.py", default="")
parser.add_argument("-f2", "--TSS-distance", dest="TSS_distance", metavar="FILE", help="input the distance of aligned read to TSS, output file from Deep_learning_read_rpocess.py", default="")
args = parser.parse_args()

#df_meth = pd.read_table('/Users/jianfengxu/Desktop/lab_work_temporary/clone_meth/output_GE_Meth/multiple_samples/all_meth/K562_GE_Me_match.txt')
df_meth = pd.read_table(args.trad_meth)
data = df_meth
if args.down_sample:
    bins = []
    bin_start = 0
    for i in range(0,200):
        bin_start, bin_end = 0.005*i, 0.005*(i+1)
        if(i != 199):
            data_bin = data[(data.meth_mean >=bin_start) & (data.meth_mean < bin_end)]
            if len(data_bin) > 60:
                data_bin = data_bin.sample(n=60)
        else:
            data_bin = data[(data.meth_mean >=bin_start) & (data.meth_mean <= bin_end)]
            if len(data_bin) > 60:
                data_bin = data_bin.sample(n=60)
        bins.append(data_bin)        
    data = pd.concat(bins)
included_gene = list(data.Gene_symbol)


#1.preprocess the data
## expression
dict_exp = defaultdict()
with open(args.expression, 'r') as f:
    lines = f.readlines()
for line in lines:
    eles = line.rstrip().split('\t')
    if (eles[0] in included_gene) == False or eles[0] == 'gene_id':
        continue
    dict_exp[eles[0]] = np.float(eles[5])

## 2D-code
with open(args.meth_code, 'r') as f:
    lines = f.readlines()
meth_list = []
gene_list = []
exp_list = []
for line in lines:
    eles = line.rstrip().split('\t')
    try: dict_exp[eles[0]]
    except KeyError: continue
    if (eles[0] in included_gene) == False :
        continue
    gene_list.append(eles[0])
    exp_list.append(dict_exp[eles[0]])
    eles = [np.float(ele) for ele in eles[1:]]
    meth_list.append(eles)
meth_set = np.array(meth_list)
exp_set = np.log2(np.array(exp_list).reshape(-1,) + 0.01)

## TSS distance
with open(args.TSS_distance,'r') as f:
    lines = f.readlines()
pos_list = []
for line in lines:
    eles = line.rstrip().split('\t')
    try: dict_exp[eles[0]]
    except KeyError: continue
    if (eles[0] in included_gene) == False :
        continue
    eles = [np.float(ele) for ele in eles[1:]]
    pos_list.append(eles)
pos_set = np.array(pos_list)


# 1.2 combine the pos_set and meth_set together
def norm_meth(meth_set):
    meth_set = (meth_set - 0.5) / 0.5
    return meth_set

meth_set = norm_meth(meth_set)

def norm_pos(pos_set, val=10):
    pos_set = pos_set.clip(min=-val, max=val)
    pos_set = pos_set / pos_set.max()
    return pos_set

pos_set = norm_pos(pos_set)

# 1.3 concatenate the numpy into a big numpy array
def merge_meth_pos(meth_set, pos_set, read_length=10):
    reads_num = pos_set.shape[1] / read_length
    if reads_num * read_length != pos_set.shape[1] or meth_set.shape[0] != pos_set.shape[0]:
        raise ValueError('alignment errors for reads')
    pos_set = pos_set.reshape(pos_set.shape[0],1,reads_num,read_length)
    meth_set = meth_set.reshape(meth_set.shape[0],1,reads_num,read_length)
    data = np.concatenate((meth_set,pos_set),axis=1)
    return data

data = torch.from_numpy(merge_meth_pos(meth_set, pos_set))
exp_set = torch.from_numpy(exp_set)

split_ratio = 0.5
split_index = int(data.shape[0]*split_ratio)
data_train = data[:split_index]
data_test = data[split_index:]
exp_train = exp_set[:split_index]
exp_test = exp_set[split_index:]

data_exp_train = TensorDataset(data_train,exp_train)
data_exp_test = TensorDataset(data_test,exp_test)
train_loader = torch.utils.data.DataLoader(data_exp_train, shuffle=True, batch_size=32)
test_loader = torch.utils.data.DataLoader(data_exp_test, shuffle=True, batch_size=32)

# let the data to be torch.FloatTensor
#2. write the cnn model parts
class cnn_meth(torch.nn.Module):
    def __init__(self):
        super(cnn_meth,self).__init__()
        #cnn 32, 4, 100, 10 -> pool 32, 4, 50, 10 -> cnn 32,16,25,10 -> pool 32,16,12,9 -> cnn 32,32,5,4 -> pool 32,32, 2,2
        self.cnn_layer = torch.nn.Sequential(torch.nn.Conv2d(2,8,kernel_size=(1,5),stride=(1,1), padding=(1,1)),
            torch.nn.BatchNorm2d(8), torch.nn.ReLU(True), torch.nn.Conv2d(8,32,kernel_size=(1,4),stride=(1,1)), torch.nn.BatchNorm2d(32), torch.nn.ReLU(True), torch.nn.Conv2d(32,128,kernel_size=(1,3),stride=(1,1)), torch.nn.BatchNorm2d(128), torch.nn.ReLU(True))
        self.reg_layer = torch.nn.Sequential(torch.nn.Linear(77568,1))
        #self.reg_layer = torch.nn.Linear(8000,1)
        self.dropout_layer = torch.nn.Dropout(p=0.2)

    def forward(self, x):
        x = self.cnn_layer(x)
        x = x.view(x.shape[0],-1)
        x = self.dropout_layer(x)
        x = self.reg_layer(x)
        return x
    

# 3. training progress
net = cnn_meth()
net = net.double() # will be removed
optimizer = torch.optim.Adam(net.parameters(), weight_decay=1.5) ## best value for weight decay??
criterion = torch.nn.MSELoss()

train_loss = []
test_loss = []
for epoch in range(20):
    loss_epoch = []
    for batch_data in train_loader:
        batch_x = Variable(batch_data[0])
        batch_y = Variable(batch_data[1])
        batch_y_pred = net(batch_x)
        # cacluate loss
        loss = criterion(batch_y_pred, batch_y)
        loss_epoch.append(loss.data[0])
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()

    # store the averaged loss per epoch
    train_loss.append(np.array(loss_epoch).mean())
    print('[epoch of {} --> Train:loss: {}]'.format(epoch, train_loss[-1]))
    
    
# Use test mode
net.eval()
loss_epoch = []
for batch_data_test in test_loader:
    batch_x_test = Variable(batch_data_test[0])
    batch_y_test = Variable(batch_data_test[1])
    batch_y_test_pred = net(batch_x_test)
    loss = criterion(batch_y_test_pred, batch_y_test)
    loss_epoch.append(loss)

test_loss.append(np.array(loss_epoch).mean())
print('[epoch of {} --> Test:loss: {}]'.format(epoch, test_loss[-1]))


import scipy.stats
test_loader = torch.utils.data.DataLoader(data_exp_test, shuffle=True, batch_size=10000)
for batch_data_test in test_loader:
    batch_x_test = Variable(batch_data_test[0])
    batch_y_test = Variable(batch_data_test[1])
    batch_y_test_pred = net(batch_x_test)

    y = batch_y_test.data.numpy().reshape(-1,)
    y_pred = batch_y_test_pred.data.numpy().reshape(-1,)
    print('Spearman correlation: {}'.format(scipy.stats.spearmanr(y,y_pred)))

plt.scatter(y_pred, y, alpha =0.5)
plt.title('Spearman correlation: {0:.2f}'.format(scipy.stats.spearmanr(y,y_pred)[0]), horizontalalignment='center', color = 'red')
plt.xlabel("Predicted Expression")
plt.ylabel("Measured Expression")
plt.savefig('{}/{}_deep_learning_prediction.png'.format(args.out_dir, args.sample), dpi=300)


