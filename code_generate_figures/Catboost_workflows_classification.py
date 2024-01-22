import pandas as pd
import numpy as np
from sklearn.model_selection import KFold
from sklearn.metrics import accuracy_score as acc, recall_score as recall, precision_score as precision, f1_score as f1, matthews_corrcoef as mcc
import random
import os
import csv
from sklearn import preprocessing
from fpgrowth_py import fpgrowth
from mlxtend.preprocessing import TransactionEncoder
from mlxtend.frequent_patterns import fpgrowth as fpgrowth_ml
import catboost as cbt
from catboost import Pool


def seed_torch(seed=1029):
    random.seed(seed)
    os.environ['PYTHONHASHSEED'] = str(seed)
    np.random.seed(seed)

seed_torch()

def read_file(name, sheet, dty, save_root):
    ranks_raw = pd.read_excel(name, header=0,
                                     index_col=None,
                                     sheet_name=sheet)

    ranks_raw.sort_values('avg_rank_mean', inplace=True)
    ranks = ranks_raw.values
    workflows = ranks[:, 0]

    #avg_med = []
    feas = []

    #peps = ['msqrob2', 'ProteoMM']
    pros = ['limma', 'ttest', 'DEP', 'DEqMS','ANOVA', 'SAM', 'ROTS', 'proDA', 'MSstats']#, 'siggenes'
    scs = ['edgeR', 'plgem', 'beta_binomial']

    for i in range(len(ranks[:, 0])):
        fea = ['None' if j =='' else j for j in workflows[i].split('|')]
        if fea[3] == 'blank':
            fea[3] = 'No_MVI'
        if fea[4] == 'blank':
            fea[4] = 'No_norm'

        feas.append(fea)

    rank_wf = np.array([i for i in range(1, len(ranks[:, 0])+1)])

    rank_wf = np.array([rank_wf], dtype='int').T
    feas = np.array([feas]).squeeze()

    labs = np.zeros((len(feas[:, 0]), 1))
    Qs = [np.percentile(rank_wf, 5), np.percentile(rank_wf,25), np.percentile(rank_wf,50)]
    labs[np.where(rank_wf <= Qs[0])[0], 0] = 1
    labs[np.where((rank_wf >= Qs[0]) & (rank_wf < Qs[1]))[0], 0] = 2
    labs[np.where((rank_wf >= Qs[1]) & (rank_wf < Qs[2]))[0], 0] = 3

    labs[np.where(rank_wf >= Qs[2])[0], 0] = 0

    feas_fp = feas[:, [0,2,3,4]]
    itms_H = feas_fp[np.where((labs==1))[0], :]
    save_path = save_root + platform + '_' + dty + '_item_H1.csv'

    with open(save_path, 'a+', newline='') as ws:
        writer = csv.writer(ws)
        writer.writerows(np.array(np.column_stack((itms_H, rank_wf[np.where((labs==1))[0]]))))
    freqItemSet_H, rules_H = fpgrowth(itms_H, minSupRatio=0.1, minConf=0.5)

    itms_L = feas_fp[np.where((labs == 0))[0], :]
    freqItemSet_L, rules_L = fpgrowth(itms_L, minSupRatio=0.1, minConf=0.5)

    te_H = TransactionEncoder()
    itms_H = te_H.fit(itms_H).transform(itms_H)
    itms_H = pd.DataFrame(itms_H, columns=te_H.columns_)
    res_H = fpgrowth_ml(itms_H, min_support=0.1, use_colnames=True)

    te_L = TransactionEncoder()
    itms_L = te_L.fit(itms_L).transform(itms_L)
    itms_L = pd.DataFrame(itms_L, columns=te_L.columns_)
    res_L = fpgrowth_ml(itms_L, min_support=0.1, use_colnames=True)

    save_path = save_root + platform + '_' + dty + 'FP_freitem_H.csv'

    with open(save_path, 'a+', newline='') as ws:
        writer = csv.writer(ws)
        writer.writerows(np.array(res_H))

    save_path = save_root + platform + '_' + dty + 'FP_freitem_L.csv'

    with open(save_path, 'a+', newline='') as ws:
        writer = csv.writer(ws)
        writer.writerows(np.array(res_L))

    OHfea = []

    uni_dea = np.array([np.unique(feas[:, 0])])
    uni_ext = np.array([np.unique(feas[:, 2])])
    uni_mvi = np.array([np.unique(feas[:, 3])])
    uni_norm = np.array([np.unique(feas[:, 4])])
    header = np.column_stack((uni_dea, uni_ext, uni_mvi, uni_norm))
    for i in range(len(feas[:, 0])):
        dea = np.zeros((1, len(uni_dea[0])))
        dea[0, np.where((uni_dea[0] == feas[i, 0]))[0][0]] = 1

        ext = np.zeros((1, len(uni_ext[0])))
        ext[0, np.where((uni_ext[0] == feas[i, 2]))[0][0]] = 1

        mvi = np.zeros((1, len(uni_mvi[0])))
        mvi[0, np.where((uni_mvi[0] == feas[i, 3]))[0][0]] = 1

        norm = np.zeros((1, len(uni_norm[0])))
        norm[0, np.where((uni_norm[0] == feas[i, 4]))[0][0]] = 1

        OH = np.column_stack((dea, ext, mvi, norm))

        OHfea.append(OH)

    OHfea = np.array(OHfea).squeeze()

    out_fea_lab = np.column_stack((feas_fp, labs))

    save_path = save_root + platform + '_' + dty + '_string_Fea_lab_DT1.csv'

    with open(save_path, 'a+', newline='') as ws:
        writer = csv.writer(ws)
        writer.writerow(['DEA', 'expression_type', 'MVI', 'Normalization', 'label'])
        writer.writerows(np.array(out_fea_lab))

    return ranks, OHfea, labs, header, feas_fp

def cv_cbt_cls(feas, labs, seed, platform, dty, save_root):

    seed_torch(seed)

    logo = KFold(n_splits=10, shuffle=True)

    pred_labs = np.zeros((len(labs)))
    Accs = []
    F1s = []
    recs = []
    pres = []
    Mccs = []
    for train_index, test_index in logo.split(labs):
        train_x = feas[train_index]
        train_lab = labs[train_index]

        test_x = feas[test_index]
        test_lab = labs[test_index]

        train_X,train_y=train_x,train_lab
        test_X,test_y=test_x,test_lab

        params = {'iterations': 1000,
                'learning_rate':0.3,
                'l2_leaf_reg':3,
                'bagging_temperature':1,
                'random_strength':1,
                'depth':6,
                'rsm':1,
                'one_hot_max_size':2,
                'leaf_estimation_method':'Gradient',
                'fold_len_multiplier':2,
                'border_count':128,}
        bst=[]
        bst = cbt.CatBoostClassifier(**params, loss_function='MultiClass')
        bst.fit(train_X, train_y, cat_features=[0, 1, 2, 3])
        ypred = bst.predict(test_X)
        ypred = np.array(ypred).astype(int)

        pred_labs[test_index] = ypred[:, 0]

        Accs.append(acc(test_y, ypred))
        F1s.append(f1(test_y, ypred, average='macro'))
        recs.append(recall(test_y, ypred, average='macro'))
        pres.append(precision(test_y, ypred, average='macro'))
        Mccs.append(mcc(test_y, ypred))

    Acc = acc(labs, pred_labs)
    F1 = f1(labs, pred_labs, average='macro')
    rec = recall(labs, pred_labs, average='macro')
    pre = precision(labs, pred_labs, average='macro')
    Mcc = mcc(labs, pred_labs)

    out = [Acc, F1, rec, pre, Mcc]
    print(out)

    out_res = np.column_stack((np.array(Accs), np.array(F1s), np.array(recs), np.array(pres), np.array(Mccs)))

    save_path = save_root + platform + '_' + dty + '_cv_cbt_cls_res_labs_str'
    with open(save_path, 'a+', newline='') as ws:
        writer = csv.writer(ws)
        writer.writerows(np.column_stack((labs,pred_labs)))

    save_path = save_root + platform + '_' + dty + '_cv_cbt_cls_res_folds_str'
    with open(save_path, 'a+', newline='') as ws:
        writer = csv.writer(ws)
        writer.writerows(out_res)
    return Mcc

def feature_importance_catboost(model):
    result=pd.DataFrame(model.get_feature_importance(),index=model.feature_names_,columns=['FeatureImportance'])
    return result.sort_values('FeatureImportance',ascending=False)


def test_CBT(fea_train, labs_train, fea_test, labs_test, seed, platform, header, dty,save_root):
    seed_torch(seed)

    params = {'iterations': 1000,
              'learning_rate': 0.3,
              'l2_leaf_reg': 3,
              'bagging_temperature': 1,
              'random_strength': 1,
              'depth': 6,
              'rsm': 1,
              'one_hot_max_size': 2,
              'leaf_estimation_method': 'Gradient',
              'fold_len_multiplier': 2,
              'border_count': 128, }

    train_X, train_y = fea_train, labs_train
    test_X, test_y = fea_test, labs_test

    pool = Pool(train_X, train_y, cat_features=[0, 1, 2, 3], feature_names=['DEA', 'exp_ty', 'MVI', 'norm'])
    bst = cbt.CatBoostClassifier(**params, loss_function='MultiClass')

    bst.fit(pool)
    ypred = bst.predict(test_X)

    ypred = np.array(ypred).astype(int)

    Acc = acc(test_y, ypred)
    F1 = f1(test_y, ypred, average='micro')
    rec = recall(test_y, ypred, average='micro')
    pre = precision(test_y, ypred, average='micro')
    Mcc = mcc(test_y, ypred)

    out = [Acc, F1, rec, pre, Mcc]
    print(out)
    save_path = save_root + platform + '_' + dty + '_test_DT_cls_res'
    with open(save_path, 'a+', newline='') as ws:
        writer = csv.writer(ws)
        writer.writerow(out)
    #return Mcc
    fea_imp = feature_importance_catboost(bst)
    fea_imp.to_csv(save_root + platform + '_' + dty + '_fea_importance.csv')

    return Acc, F1, Mcc

if __name__ == '__main__':
    root = 'data/'
    save_root = 'data/'
    filename = root + 'ranks_all_Fragpipe_TMT.xlsx'
    sheet = 'ranking_all'
    platform = 'Fragpipe'
    dty = 'TMT'
    seed = 100
    ranks, OHfea, labs, header, feas_str = read_file(filename, sheet, dty,save_root)
    cv_cbt_cls(feas_str, labs, seed, platform, dty,save_root)
    test_CBT(feas_str, labs, feas_str, labs, seed, platform, header, dty, save_root)

    filename = root + 'ranks_all_Fragpipe_DDA.xlsx'
    sheet = 'ranking_all'
    platform = 'Fragpipe'
    dty = 'DDA'
    seed = 100
    ranks, OHfea, labs, header, feas_str = read_file(filename, sheet, dty,save_root)
    cv_cbt_cls(feas_str, labs, seed, platform, dty,save_root)
    test_CBT(feas_str, labs, feas_str, labs, seed, platform, header, dty,save_root)

    filename = root + 'ranks_all_Maxquant_DDA.xlsx'
    sheet = 'ranking_all'
    platform = 'Maxquant'
    dty = 'DDA'

    ranks, OHfea, labs, header, feas_str = read_file(filename, sheet, dty,save_root)
    cv_cbt_cls(feas_str, labs, seed, platform, dty,save_root)
    test_CBT(feas_str, labs, feas_str, labs, seed, platform, header, dty,save_root)

    filename = root + 'ranks_all_DIANN_DIA.xlsx'
    sheet = 'ranking_all'
    platform = 'DIANN'
    dty = 'DIA'
    seed = 100
    ranks, OHfea, labs, header, feas_str = read_file(filename, sheet, dty,save_root)
    cv_cbt_cls(feas_str, labs, seed, platform, dty,save_root)
    test_CBT(feas_str, labs, feas_str, labs, seed, platform, header, dty,save_root)

    filename = root + 'ranks_all_spt_DIA.xlsx'
    sheet = 'ranking_all'
    platform = 'spt'
    dty = 'DIA'

    ranks, OHfea, labs, header, feas_str = read_file(filename, sheet, dty,save_root)
    cv_cbt_cls(feas_str, labs, seed, platform, dty,save_root)
    test_CBT(feas_str, labs, feas_str, labs, seed, platform, header, dty,save_root)


