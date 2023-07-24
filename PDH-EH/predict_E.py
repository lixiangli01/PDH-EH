import numpy as np
from sklearn import preprocessing
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GridSearchCV

from sklearn import metrics
from sklearn.metrics import roc_curve
from sklearn.metrics import confusion_matrix
from sklearn.metrics import precision_recall_curve,auc
from sklearn.metrics import make_scorer
from sklearn.model_selection import cross_val_score
#import xgboost as xgb
import pandas as pd
import math
import pickle


dictaa = {'GLY': 'G', 'ALA': 'A', 'VAL': 'V', 'LEU': 'L',
          'ILE': 'I', 'PHE': 'F', 'TRP': 'W', 'TYR': 'Y',
          'ASP': 'D', 'ASN': 'N', 'GLU': 'E', 'LYS': 'K',
          'GLN': 'Q', 'MET': 'M', 'SER': 'S', 'THR': 'T',
          'CYS': 'C', 'PRO': 'P', 'HIS': 'H', 'ARG': 'R','MSE':'M'}

aalist = ['GLY', 'ALA', 'VAL', 'LEU', 'ILE',
          'PHE', 'TRP', 'TYR', 'ASP', 'ASN',
          'GLU', 'LYS', 'GLN', 'MET', 'SER',
          'THR', 'CYS', 'PRO', 'HIS', 'ARG','MSE']

ls=['G','A','V','L','I',
    'P','F','Y','W','S',
    'T','C','M','N','Q',
    'D','E','K','R','H']


def ReadPDB(dirfile,n):
    dirstr = dirfile + '.pdb'
    with  open(dirstr, "r") as myfile:
        num = 0
        sum = 0
        ls_sort = []
        ls_sort_new = []
        ls_seq_pdb = []
        for line in myfile:
            if line[0:14] == 'MODEL        2':
                break
            elif line[0:4] == 'ATOM' and line[17:20] in aalist and num != int(line[22:26].replace(' ', ''))and line[21]==n or (line[0:4] == 'ANIS' and line[17:20]=='MSE' and num != int(line[22:26].replace(' ', ''))and line[21]==n):
                
                sum += 1
                ls_sort.append(line[22:26].replace(' ', ''))
                ls_sort_new.append(str(sum))
                ls_seq_pdb.append(line[17:20])
                
                num = int(line[22:26].replace(' ', ''))
        
    myfile.close()
    strAA = ''
    for i in ls_seq_pdb:
        strAA += dictaa[i]
    return ls_sort, ls_sort_new, strAA



def findint(s1,s2):
    len1 = len(s1)
    len2 = len(s2)
    xmax = 0
    xindex = 0
    
    matrix = [[0]*len1 for i in range(len2)]
    for i,x in enumerate(s2):
        for j ,y in enumerate(s1):  
            if x==y: 
                if i==0 or j ==0:       
                    matrix[i][j]=1     
                else:
                    matrix[i][j] = matrix[i-1][j-1]+1   
                if matrix[i][j] > xmax and xmax<=30:  
                    xmax = matrix[i][j]  
                    xindex = j  
    return s1[xindex-xmax+1:xindex+1]

def getfastaallline(filename):
    dirname=filename+'.fasta'
    readfile = open(dirname, 'r')
    lsline = []
    for line in readfile:
        lsline.append(line.replace('\n', ''))
    return lsline

def oneline(mystr):
    for i in mystr:
        if i not in ls:
            break
        else:
            strword = ''
            for word in mystr:
                if word not in strword:
                    strword += word
            if len(strword) > 5 : 
                return mystr
            break

def getallfasta(name):
    tmpprolist = []
    for i in getfastaallline(name):
        if oneline(i) != None:
            
            tmpprolist.append(oneline(i))
    return tmpprolist

def getOneChainPro(name,chain):#name,chain
    PDBseq=ReadPDB(name,chain)[2]
    tempro=getallfasta(name)
    temi = ''
    num=0
    for i in tempro:
        if len(findint(PDBseq, i))>num:
            num=len(findint(PDBseq, i))
            temi=i
    return temi
        

def getfastasite(name,chain,Res,pdbsite):
    pdbname = name
    Res = Res
    pdbsite = pdbsite
    chain = chain
    
    reslist_pdb = ReadPDB(pdbname, chain)[2]
    sortlist_pdb = ReadPDB(pdbname, chain)[0]
    res_index = sortlist_pdb.index(pdbsite)


    if res_index - 4 in range(len(reslist_pdb)):
        a=reslist_pdb[res_index - 4]
    else:
        a='-4'

    if res_index - 3 in range(len(reslist_pdb)):
        b=reslist_pdb[res_index - 3]
    else:
        b='-3'

    if res_index - 2 in range(len(reslist_pdb)):
        c=reslist_pdb[res_index - 2]
    else:
        c='-2'

    if res_index - 1 in range(len(reslist_pdb)):
        d=reslist_pdb[res_index - 1]
    else:
        d='-1'

    if res_index  in range(len(reslist_pdb)):
        e=reslist_pdb[res_index]
    else:
        e='0'

    if res_index +1 in range(len(reslist_pdb)):
        f=reslist_pdb[res_index+1]
    else:
        f='1'

    if res_index +2 in range(len(reslist_pdb)):
        g=reslist_pdb[res_index+2]
    else:
        g='2'

    if res_index +3 in range(len(reslist_pdb)):
        h=reslist_pdb[res_index+3]
    else:
        h='3'

    if res_index +4 in range(len(reslist_pdb)):
        i=reslist_pdb[res_index+4]
    else:
        i='4'
    
    strtmp1 = c+d+e+f+g
    strtmp2 = e+f+g+h+i
    strtmp3 = a+b+c+d+e

    seq = getOneChainPro(pdbname, chain)

    if seq[seq.find(strtmp1) + 2] == Res:
        return seq.find(strtmp1) + 2, seq
    elif seq[seq.find(strtmp2)] == Res:
        return seq.find(strtmp2), seq
    elif seq[seq.find(strtmp3)+4] == Res:
        return seq.find(strtmp3)+4, seq
    else:
        print(pdbname,chain,Res,pdbsite)
        raise 'check again'


import urllib.request
import time

def downOnePdb(filename):
    def getHtml(url):
        html = urllib.request.urlopen(url).read()
        return html

    def saveHtml(file_name, file_content):
        
        with open(file_name.replace('/', '_') + ".pdb", "wb") as f:
            
            f.write(file_content)

    aurl = "https://files.rcsb.org/view/"+filename+'.pdb'
    
    html = getHtml(aurl)
    
    saveHtml(filename, html)
    time.sleep(2)
    #print("Successfully download the pdb file"+':'+filename)

def downOneFASTA(filename):
    def getHtml(url):
        html = urllib.request.urlopen(url).read()
        return html

    def saveHtml(file_name, file_content):
        
        with open(file_name.replace('/', '_') + ".fasta", "wb") as f:
            
            f.write(file_content)
    aurl="https://www.rcsb.org/fasta/entry/"+filename+'/display'
    html = getHtml(aurl)

    saveHtml(filename, html)
    time.sleep(2)
    #print("Successfully download the fasta file"+':'+filename)

import argparse
parser = argparse.ArgumentParser(description='Calculate volume of a cylinder')
parser.add_argument('PDB_name', type=str, help='PDB uppercase')
parser.add_argument('CHAIN', type=str, help='only one chain')
parser.add_argument('RES', type=str, help='Abbreviation for amino acid residues, such as Y')
parser.add_argument('PDBsite', type=str, help='The location of the residue in the PDB file')
args = parser.parse_args()


PDB_name=args.PDB_name#'1TN9'
CHAIN=args.CHAIN
RES=args.RES
PDBsite=args.PDBsite


# PDB_name='1TN9'
# CHAIN='A'
# RES='Y'
# PDBsite='40'

downOnePdb(PDB_name)
downOneFASTA(PDB_name)
siteinfasta,sequ=getfastasite(PDB_name,CHAIN,RES,PDBsite)
#print(a,b)



#1AAY,D120A,D19A,0,MERPYACPVESCDRRFSRSDELTRHIRIHTGQKPFQCRICMRNFSRSDHLTTHIRTHTGEKPFACDICGRKFARSDERKRHTKIHLRQKD
from transformers import BertModel, BertTokenizer
import re
import numpy


tokenizer = BertTokenizer.from_pretrained('/PDH-EH/set_file', do_lower_case=False )
model = BertModel.from_pretrained("/PDH-EH/set_file")
##################################################



def genstr(str1):
    c = ''
    for j in str1:
        c += j + ' '
    mystr = c[0:len(c) - 1]
    return mystr


def genbertfeature(seq):
    
    mystr=genstr(seq)
    #sequence_Example = "C T G C C T C T G T G T C T T G T C A C C A G G C C T T A C C C T G G G G A C C C C T G C T C C C A G C G G A G C C A G T A G T G A T G A C A G G C G C A G C T G G G A G C A G C T T G G T A G A C C T A G G G G G T C T T T C T A G A A G C C A A G G G G G C C C T T G G C A C A C A C A T G T G G A T G C A G G G C T G C C C A C C C A A C A C T G C T G A G C C C A C A C A G G C C C A A G A G A A A A GC T G C C T C T G T G T C T T G T C A C C A G G C C T T A C C C T G G G G A C C C C T G C T C C C A G C G G A G C C A G T A G T G A T G A C A G G C G C A G C T G G G A G C A G C T T G G T A G A C C T A G G G G G T C T T T C T A G A A G C C A A G G G G G C C C T T G G C A C A C A C A T G T G G A T G C A G G G C T G C C C A C C C A A C A C T G C T G A G C C C A C A C A G G C C C A A G A G A A A A GC T G C C T C T G T G T C T T G T C A C C A G G C C T T A C C C T G G G G A C C C C T G C T C C C A G C G G A G C C A G T A G T G A T G A C A G G C G C A G C T G G G A G C A G C T T G G T A G A C C T A G G G G G T C T T T C T A G A A G C C A A G G G G G C C C T T G G C A C A C A C A T G T G G A T G C A G G G C T G C C C A C C C A A C A C T G C T G A G C C C A C A C A G G C C C A A G A G A A A A GC T G C C T C T G T G T C T T G T C A C C A G G C C T T A C C C T G G G G A C C C C T G C T C C C A G C G G A G C C A G T A G T G A T G A C A G G C G C A G C T G G G A G C A G C T T G G T A G A C C T A G G G G G T C T T T C T A G A A G C C A A G G G G G C C C T T G G C A C A C A C A T G T G G A T G C A G G G C T G C C C A C C C A A C A C T G C T G A G C C C A C A C A G G C C C A A G A G A A A A G"
    sequence_Example=mystr
    sequence_Example = re.sub(r"[UZOB]", "X", sequence_Example)
    encoded_input = tokenizer(sequence_Example, return_tensors='pt')
    output = model(**encoded_input)

    feaoutput = output[0][0].detach().numpy()

    return feaoutput


def extrafun(featuredata,i):

    data_cut = featuredata[i+1]
    data_cut=data_cut.tolist()

    return data_cut


featuredata=genbertfeature(sequ)


testfeature=extrafun(featuredata,siteinfasta)

X1 = np.array(testfeature).reshape(1, -1)


def readtxt(filename):
    label_train=[]
    l_train = []
    with open(filename,'r') as f:
        for line in f:
            line = line.split( )
            l_train.append([float(x) for x in line[3:]])
            label_train.append(int(line[2]))
    return l_train,label_train


(X,Y) = readtxt('train_BFD_1024.txt')
X = np.array(X)
le = preprocessing.LabelEncoder()
le.fit(Y)
y_train = le.transform(Y)
min_max_scaler = preprocessing.MinMaxScaler()
X_train = min_max_scaler.fit_transform(X)

X_test = min_max_scaler.transform(X1)

with open('sort_BFD_1024_mrmr.pickle', 'rb') as f2:
    idx_sorted = pickle.load(f2)

i=945#1024-79

X_test_tmp = X_test[:,idx_sorted[i:]]
f1 = open('rf_mrmr_1024_79.pickle','rb')
rlf1 = pickle.load(f1)
x_test = np.array(X_test_tmp)
y_pred = rlf1.predict(x_test)
y_scores = rlf1.predict_proba(x_test)[:,1]
print(y_scores)

import os
pdbfilename=PDB_name+'.pdb'
fastafilename=PDB_name+'.fasta'
os.remove(pdbfilename)
os.remove(fastafilename)

testset=False
if  testset == True:
    (X2,Y2) = readtxt('test_BFD_1024.txt')#######
    X2 = np.array(X2)

    y_test2 = le.transform(Y2)

    X_test2 = min_max_scaler.transform(X2)
    X_test_tmp2 = X_test2[:,idx_sorted[i:]]

    y_pred2 = rlf1.predict(X_test_tmp2)
    y_scores2 = rlf1.predict_proba(X_test_tmp2)[:,1]
    TN, FP, FN, TP = confusion_matrix(y_test2, y_pred2).ravel()
    print('TN, FP, FN, TP:',TN, FP, FN, TP)

    Specificity=TN/(TN+FP)
    ACC = (TP + TN) / (TP + FP + FN + TN)
    Precision=TP/(TP+FP)
    Recall=TP/(TP+FN)
    F1Score=2*TP/(2*TP+FP+FN)
    if math.sqrt(float(TP + FP) * float(TP + FN) * float(TN + FP) * float(TN + FN))==0:
        MCC=0
    else:
        MCC=float(TP*TN-FP*FN)/ math.sqrt(float(TP + FP) * float(TP + FN) * float(TN + FP) * float(TN + FN))
    p, r, thresh = metrics.precision_recall_curve(y_test2, y_scores2)
    pr_auc = metrics.auc(r, p)
    ro_auc = metrics.roc_auc_score(y_test2, y_scores2)
    print('Specificity_test,ACC_test,Precision_test,Recall_test,F1Score_test,MCC_test,auprc_test,auroc_test:',
        Specificity,ACC,Precision,Recall,F1Score,MCC,pr_auc,ro_auc)
