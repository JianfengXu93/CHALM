import torch
import numpy as np
import os
import random
import time
import scipy
import scipy.stats as ss
from torch.utils.data import TensorDataset
from torch.autograd import Variable
import pandas as pd
import glob
import argparse
from collections import defaultdict
import matplotlib.pyplot as plt

# normalize the position to let it between -1 to 1
def norm_pos(pos_set, val=10):
    pos_set = pos_set.clip(min=-val, max=val)
    pos_set = pos_set / pos_set.max()
    return pos_set

# normalize the methylation to let it between -1 to 1
def norm_meth(meth_set):
    meth_set = (meth_set - 0.5) / 0.5
    return meth_set

def merge_meth_pos(meth_set, pos_set, read_length=10):
    reads_num = pos_set.shape[1] // read_length
    if reads_num * read_length != pos_set.shape[1] or meth_set.shape[0] != pos_set.shape[0]:
        raise ValueError('alignment errors for reads')
    pos_set = pos_set.reshape(pos_set.shape[0],1,reads_num,read_length)
    meth_set = meth_set.reshape(meth_set.shape[0],1,reads_num,read_length)
    data = np.concatenate((meth_set,pos_set),axis=1)
    return data

def predict_exp_corr(net, data_loader):
    y_pred = []
    y = []
    for batch_data in data_loader:
        batch_x_test = Variable(batch_data[0], volatile=True)
        batch_y_test = batch_data[1].float()
        batch_y_test_pred = net(batch_x_test)
        y_pred.append(batch_y_test_pred.data.numpy())
        y.append(batch_y_test.numpy())

    y_pred_np = np.concatenate(y_pred).squeeze()
    y_np = np.concatenate(y)
    return ss.pearsonr(y_np,y_pred_np), ss.spearmanr(y_np,y_pred_np), y_np, y_pred_np

def file_process(trad_meth, expression, meth_code, TSS_distance):
    df_meth = pd.read_table(trad_meth)
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
    with open(expression, 'r') as f:
        lines = f.readlines()
    for line in lines:
        eles = line.rstrip().split('\t')
        if (eles[0] in included_gene) == False or eles[0] == 'gene_id':
            continue
        dict_exp[eles[0]] = np.float(eles[5])

    ## 2D-code
    with open(meth_code, 'r') as f:
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
    meth_data = np.array(meth_list)
    exp_data = np.log2(np.array(exp_list).reshape(-1,) + 0.01)

    ## TSS distance
    with open(TSS_distance,'r') as f:
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
    pos_data = np.array(pos_list)
    return meth_data, exp_data, pos_data

class cnn_meth(torch.nn.Module):
    def __init__(self):
        super(cnn_meth,self).__init__()
        #cnn 32, 4, 100, 10 -> pool 32, 4, 50, 10 -> cnn 32,16,25,10 -> pool 32,16,12,9 -> cnn 32,32,5,4 -> pool 32,32, 2,2
        self.cnn_layer = torch.nn.Sequential(torch.nn.Conv2d(2,8,kernel_size=(1,5),stride=(1,1), padding=(1,1)),
            torch.nn.BatchNorm2d(8), torch.nn.ReLU(True), torch.nn.Conv2d(8,32,kernel_size=(1,4),stride=(1,1)), torch.nn.BatchNorm2d(32), torch.nn.ReLU(True), torch.nn.Conv2d(32,128,kernel_size=(1,3),stride=(1,1)), torch.nn.BatchNorm2d(128), torch.nn.ReLU(True))
        self.reg_layer = torch.nn.Sequential(torch.nn.Linear(77568,128), torch.nn.Linear(128,1))
        #self.reg_layer = torch.nn.Linear(8000,1)
        self.dropout_layer = torch.nn.Dropout(p=0.2)

    def forward(self, x):
        x = self.cnn_layer(x)
        x = x.view(x.size(0),-1)
        x = self.dropout_layer(x)
        x = self.reg_layer(x)
        return x
        

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Use pretained Deep Neural Network to predict the justified methylation ratio')
    parser.add_argument('--model', dest='model', help='Pretained model to be loaded for prediction', default="")
    parser.add_argument("-o", "--output-dir", dest="out_dir", metavar="FILE", help="output methylation ratio file name. [default: STDOUT]", default="")
    parser.add_argument("-m", "--meth", dest="trad_meth", metavar="FILE", help="traditional mean methylation level", default="")
    parser.add_argument("-e", "--exp", dest="expression", metavar="FILE", help="gene expression file for prediction", default="")
    parser.add_argument("-s", "--sample-name", dest="sample", metavar="SAMPLE", help="", default="input_sample")
    parser.add_argument("-d", "--down-sample", dest="down_sample", action="store_true", default=False)
    parser.add_argument("-f1", "--2D-code", dest="meth_code", metavar="FILE", help="input methylation 2D code of regions of interest, output file from Deep_learning_read_rpocess.py", default="")
    parser.add_argument("-f2", "--TSS-distance", dest="TSS_distance", metavar="FILE", help="input the distance of aligned read to TSS, output file from Deep_learning_read_rpocess.py", default="")
    args = parser.parse_args()

    print('Reading file:' + args.meth_code)
    print('Reading file:' + args.TSS_distance)
    print('Reading file:' + args.expression)
    print('Reading file:' + args.trad_meth)
    print('Reading file:' + args.model)
    
    # preprocessing the input data
    meth_data, exp_data, pos_data = file_process(args.trad_meth, args.expression, args.meth_code, args.TSS_distance)
    meth_data = norm_meth(meth_data)
    pos_data = norm_pos(pos_data, val=12)
    data = torch.from_numpy(merge_meth_pos(meth_data, pos_data)).float()
    exp_data = torch.from_numpy(exp_data).float()
    data_exp_train = TensorDataset(data,exp_data)
    data_loader = torch.utils.data.DataLoader(data_exp_train, batch_size=32)

    # reload the pretained module
    net = cnn_meth()
    net = torch.load(args.model)
    net.cpu().eval() # put the model into cpu mode

    # calculate the correlation
    pearson, spearman, exp, exp_pred = predict_exp_corr(net, data_loader)
    print(pearson[0], spearman[0])
    plt.scatter(exp_pred, exp, alpha =0.5)
    plt.title('Spearman correlation: {0:.2f}'.format(scipy.stats.spearmanr(exp,exp_pred)[0]), horizontalalignment='center', color = 'red')
    plt.xlabel("Predicted Expression")
    plt.ylabel("Measured Expression")
    plt.savefig('{}/{}_deep_learning_prediction.png'.format(args.out_dir, args.sample), dpi=300)





    


