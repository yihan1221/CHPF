import pandas as pd
import sys
import lightgbm as lgb
import numpy as np
from sklearn.model_selection import KFold
from sklearn.metrics import recall_score
from scipy import interpolate
import random
import os
import argparse
import subprocess

parser = argparse.ArgumentParser(description= 'Cellular hypoxia status predicting.')
parser.add_argument('-n',  '--n_splits',  action="store", dest="n_splits", default='5', type=int,  help='The fold of cross-validation.Default value is 5.')
parser.add_argument('-t',  '--tree',  action="store", dest="tree", default='100',  type=int,help='The number of decision tree bulit in one cross-validation.Default value is 100.')
parser.add_argument('-T',  '--THEO',  action="store", dest="THEO", default='0.5', type=float,help='The threshold value of the predicting outcome.Default value is 0.5.')
parser.add_argument('-g',  '--gmt', action="store", dest="gmt", help='The input file of gene sets.',required=True)
parser.add_argument('-e',  '--expr',  action="store", dest="expr", help='The input file of expression profile from scRNA-seq.Columns correspond to cells.Rows correspond to genes.',required=True)
parser.add_argument('-o',  '--outpath', action="store", dest="outpath",  help='The output path.',required=True)
args = parser.parse_args()
gmt=os.path.abspath(args.gmt)
expr=os.path.abspath(args.expr)

outputFile=args.outpath
n_splits=args.n_splits
T=args.tree
THEO=args.THEO

os.makedirs(outputFile+"/temp", exist_ok=True)
os.makedirs(outputFile+"/feature", exist_ok=True)
temp=outputFile+"/temp"
feature=outputFile+"/feature/"

code1=subprocess.call(['Rscript', 'src/hypoxia_score_calculation.r',gmt,expr,temp])
if code1 != 0 :
    print("code: ", code1)
    sys.exit(1)
code2=subprocess.call(['Rscript', 'src/Gaussian_mixture_model.r',temp])
if code2 != 0 :
    rint("code: ", code2)
    sys.exit(1)
code3=subprocess.call(['Rscript', 'src/high_confidence_cells.r',temp])
if code3 != 0 :
    rint("code: ", code3)
    sys.exit(1)


X_train_value = pd.read_csv(temp+'/expr_highconfi.csv', sep=',', index_col=0)
X_test_value = pd.read_csv(temp+'/expr_others.csv', sep=',', index_col=0)
Y_train = pd.read_csv(temp+'/label_highconfi.csv', sep=',', index_col=0)
N_t = len(X_test_value)


params = {
    'task': 'train',
    'boosting_type': 'gbdt',  
    'objective': 'multiclass',
    'num_class': 2,  
    'num_leaves': 120, 
    'max_depth': 7, 
    'min_data_in_leaf': 50,
    'learning_rate': 0.01,  
    "min_split_gain": 0.1,
    'feature_fraction': 0.8,  
    'bagging_fraction': 0.8,  
    'bagging_freq': 5,  
    'verbose': -1,  
}

num_round = 10  

seed = 12234  
cv_pred = []  
kfold = KFold(n_splits=n_splits, random_state=seed, shuffle=True)
index = kfold.split(X=X_train_value, y=Y_train)
k = 0 
total_re = []

for train_index, test_index in index:
    k = k + 1
    train_value = X_train_value.iloc[train_index] 
    train_target = Y_train.iloc[train_index]
    test_value = X_train_value.iloc[test_index]
    test_target = Y_train.iloc[test_index]
    N = len(train_target.loc[train_target['group'] == 1])  
    h_train_index_all = train_index[0:N]  
    M = len(train_target.loc[train_target['group'] == 0])  
    n_train_index_all = train_index[N:]  
    Px = np.zeros(len(test_index))  
    Px_test = np.zeros(N_t)  
    sum_Wi = 0  


    for i in range(0, T):
        if N<=M :
            ran = np.array(sorted(random.sample(range(0, M), N)))
            n_train_index = n_train_index_all[ran]
            h_train_index = h_train_index_all
        else:
            ran = np.array(sorted(random.sample(range(0, N), M)))
            h_train_index = h_train_index_all[ran]
            n_train_index = n_train_index_all

        tree_train_index = np.union1d(h_train_index, n_train_index)
        tree_test_index = test_index
        x_train = X_train_value.iloc[tree_train_index]
        y_train = Y_train.iloc[tree_train_index]
        x_test = X_train_value.iloc[tree_test_index]
        y_test = Y_train.iloc[tree_test_index]
        train_data = lgb.Dataset(x_train, label=y_train)  
        test_data = lgb.Dataset(x_test, label=y_test)  

        from lightgbm import log_evaluation,early_stopping
        callbacks = [log_evaluation(period=100),early_stopping(stopping_rounds=30)]

        model = lgb.train(params, train_data, num_boost_round=100000, valid_sets=[test_data], callbacks=callbacks)  
        x_pred = model.predict(x_test, num_iteration=model.best_iteration)
        x_pred = [np.argmax(x) for x in x_pred]
        re = recall_score(y_test, x_pred, average='binary') 

        sum_Wi += re
        total_re.append(re)
        P = np.array(x_pred) * re
        Px = Px + P

        pd.DataFrame({
            'column': model.feature_name(),
            'importance': model.feature_importance(),
        }).sort_values(by='importance', ascending=False).to_csv(feature + str((k-1) * T + i) + '.txt')

        y_test = model.predict(X_test_value, num_iteration=model.best_iteration)  
        y_test = [np.argmax(x) for x in y_test]
        Px_test += np.array(y_test) * re 
        print(i)
    Px = Px / sum_Wi


    Px_test = Px_test / sum_Wi  
    target = np.zeros(N_t)
    target[np.where(Px_test > THEO)] = 1  
    target = [int(x) for x in target]
    if k == 1:
        cv_pred = np.array(target).reshape(-1, 1)
    else:
        cv_pred = np.hstack((cv_pred, np.array(target).reshape(-1, 1)))



print(cv_pred)

submit = []
for line in cv_pred:
    submit.append(np.argmax(np.bincount(line)))

df_test = pd.DataFrame()
df_test['id'] = list(X_test_value.index)
df_test['predict'] = submit
df_test.to_csv(temp+'/label_others.csv', index=False)

os.remove(temp+'/high_confi_cellID.RData')
os.remove(temp+'/hypoxia_score.RData')
os.remove(temp+'/score_cluster.RData')
os.remove(temp+'/scRNA_matrix.RData')

ofile = open(outputFile+'/cell_status.csv', 'w')
ofile.write(','.join(['cell', 'status'])+'\n')
file1 = temp+'/label_highconfi.csv'
file2 = temp+'/label_others.csv'
with open(file1) as f:
    line = f.readline()
    while True:
        line = f.readline()
        if not line:
            break
        ofile.write(line)
with open(file2) as f:
    line = f.readline()
    while True:
        line = f.readline()
        if not line:
            break
        ofile.write(line)
ofile.close()


