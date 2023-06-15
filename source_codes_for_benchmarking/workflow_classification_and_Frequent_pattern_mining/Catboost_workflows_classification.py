import pandas as pd
import numpy as np
from sklearn.model_selection import KFold
from sklearn.metrics import accuracy_score as acc, recall_score as recall, precision_score as precision, f1_score as f1, matthews_corrcoef as mcc
import random
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
import os
import optuna
import csv
from sklearn import preprocessing
from sklearn.tree import DecisionTreeClassifier as DT
from sklearn import tree
import matplotlib.pyplot as plt
from fpgrowth_py import fpgrowth
from mlxtend.preprocessing import TransactionEncoder
from mlxtend.frequent_patterns import fpgrowth as fpgrowth_ml
import xgboost as xg
from xgboost import XGBClassifier
import matplotlib
import catboost as cbt
from catboost import Pool
#import shap

def seed_torch(seed=1029):
    random.seed(seed)
    os.environ['PYTHONHASHSEED'] = str(seed)
    np.random.seed(seed)
    #torch.manual_seed(seed)
    #torch.cuda.manual_seed(seed)
    #torch.backends.cudnn.deterministic = True
    #torch.backends.cudnn.benchmark = False

seed_torch()

def read_file(name, sheet):
    ranks = pd.read_excel(name, header=0,
                                     index_col=None,
                                     sheet_name=sheet).values

    workflows = ranks[:, 0]

    avg_med = []
    feas = []

    peps = ['msqrob2', 'ProteoMM']
    pros = ['limma', 'siggenes', 'ttest', 'DEP', 'DEqMS','ANOVA', 'SAM', 'ROTS', 'proDA', 'MSstats']
    scs = ['edgeR', 'plgem', 'beta_binomial']

    for i in range(len(ranks[:, 0])):
        fea = ['None' if j =='' else j for j in workflows[i].split('|')]
        if fea[3] == 'None':
            fea[3] = 'No_MVI'
        if fea[4] == 'None':
            fea[4] = 'No_norm'
        if fea[0] in scs:
            fea[2] = 'sc'
        elif (fea[0] in peps) & (fea[2] == 'None'):
            fea[2] = 'pep_inten'
        elif (fea[0] in peps) & ((fea[2] == 'MaxLFQ') | (fea[2] == 'LFQ')) :
            fea[2] = 'pep_LFQ'
        elif (fea[0] in pros) & (fea[2] == 'None'):
            fea[2] = 'pro_inten'
        elif (fea[0] in pros) & ((fea[2] == 'MaxLFQ') | (fea[2] == 'LFQ')):
            fea[2] = 'pro_LFQ'

        feas.append(fea)
        avg_med.append((ranks[i, 3] + ranks[i, 7] + ranks[i, 11] + ranks[i, 15] + ranks[i, 19])/5)

    rank_wf = np.array([i for i in range(1, len(ranks[:, 0])+1)])

    rank_wf = np.array([rank_wf], dtype='int').T
    feas = np.array([feas]).squeeze()

    labs = np.zeros((len(avg_med), 1))
    Qs = [np.percentile(rank_wf, 5), np.percentile(rank_wf,25), np.percentile(rank_wf,50)]
    labs[np.where(rank_wf <= Qs[0])[0], 0] = 1
    labs[np.where((rank_wf >= Qs[0]) & (rank_wf < Qs[1]))[0], 0] = 2
    labs[np.where((rank_wf >= Qs[1]) & (rank_wf < Qs[2]))[0], 0] = 3

    labs[np.where(rank_wf >= Qs[2])[0], 0] = 0

    feas_fp = feas[:, [0,2,3,4]]
    itms_H = feas_fp[np.where((labs==1))[0], :]
    save_path = root + platform + '_item_H1.csv'
    #out = np.column_stack((hurdle_real, hurdle_real_ctrs))
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

    save_path = root + platform + 'FP_freitem_H.csv'
    #out = np.column_stack((hurdle_real, hurdle_real_ctrs))
    with open(save_path, 'a+', newline='') as ws:
        writer = csv.writer(ws)
        writer.writerows(np.array(res_H))

    save_path = root + platform + 'FP_freitem_L.csv'
    #out = np.column_stack((hurdle_real, hurdle_real_ctrs))
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

    save_path = root + platform + '_string_Fea_lab_DT1.csv'
    #out = np.column_stack((hurdle_real, hurdle_real_ctrs))
    with open(save_path, 'a+', newline='') as ws:
        writer = csv.writer(ws)
        writer.writerow(['DEA', 'expression_type', 'MVI', 'Normalization', 'label'])
        writer.writerows(np.array(out_fea_lab))

    return ranks, OHfea, labs, header, feas_fp

def cv_RF_cls(feas, labs, seed, platform, N, md):

    seed_torch(seed)
    logo = KFold(n_splits=10, shuffle=True)

    pred_labs = np.zeros((len(labs)))

    for train_index, test_index in logo.split(labs):
        train_fea = feas[train_index]
        train_lab = labs[train_index]

        test_fea = feas[test_index]
        test_lab = labs[test_index]

        rf = RandomForestClassifier(n_estimators=N, max_depth=md, class_weight='balanced', random_state=seed)
        rf.fit(train_fea, train_lab)

        pred_l = rf.predict(test_fea)

        pred_labs[test_index] = pred_l

    Acc = acc(labs, pred_labs)
    F1 = f1(labs, pred_labs, average='micro')
    Mcc = mcc(labs, pred_labs)

    out = [N, md, Acc, F1, Mcc]
    print(out)
    save_path = root + platform + '_' + 'cv_RF_cls_res'
    with open(save_path, 'a+', newline='') as ws:
        writer = csv.writer(ws)
        writer.writerow(out)
    return Mcc

def cv_DT_cls(feas, labs, seed, platform):

    seed_torch(seed)
    logo = KFold(n_splits=10, shuffle=True)

    pred_labs = np.zeros((len(labs)))
    Accs = []
    F1s = []
    recs = []
    pres = []
    Mccs = []
    for train_index, test_index in logo.split(labs):
        train_fea = feas[train_index]
        train_lab = labs[train_index]

        test_fea = feas[test_index]
        test_lab = labs[test_index]

        rf = DT(random_state=seed)
        rf.fit(train_fea, train_lab)

        pred_l = rf.predict(test_fea)

        pred_labs[test_index] = pred_l

        Accs.append(acc(test_lab, pred_l))
        F1s.append(f1(test_lab, pred_l, average='macro'))
        recs.append(recall(test_lab, pred_l, average='macro'))
        pres.append(precision(test_lab, pred_l, average='macro'))
        Mccs.append(mcc(test_lab, pred_l))

    Acc = acc(labs, pred_labs)
    F1 = f1(labs, pred_labs, average='macro')
    rec = recall(labs, pred_labs, average='macro')
    pre = precision(labs, pred_labs, average='macro')
    Mcc = mcc(labs, pred_labs)

    out = [Acc, F1, rec, pre, Mcc]
    print(out)

    out_res = np.column_stack((np.array(Accs), np.array(F1s), np.array(recs), np.array(pres), np.array(Mccs)))

    save_path = root + platform + '_' + 'cv_DT_cls_res_labs_str'
    with open(save_path, 'a+', newline='') as ws:
        writer = csv.writer(ws)
        writer.writerows(np.column_stack((labs,pred_labs)))

    save_path = root + platform + '_' + 'cv_DT_cls_res_folds_str'
    with open(save_path, 'a+', newline='') as ws:
        writer = csv.writer(ws)
        writer.writerows(out_res)
    return out

def format_data_for_xg(fea, lab):
    df_data = pd.DataFrame(fea, columns=['DEA', 'exp_ty', 'MVI', 'norm'])
    labs = pd.DataFrame(lab, columns=['lab'])
    X, y = df_data, \
    labs[['lab']]
    cats = X.select_dtypes(exclude=np.number).columns.tolist()
    for col in cats:
        X[col] = X[col].astype('category')

    y = y.fillna(0)

    ncats = X.select_dtypes(exclude='category').columns.tolist()
    for col in ncats:
        X[col] = X[col].fillna(0)

    return X, y

def cv_XG_cls(feas, labs, seed, platform):

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

        dtrain = xg.DMatrix(train_X, label=train_y, enable_categorical=True)
        dtest = xg.DMatrix(test_X, enable_categorical=True)

        params = {'booster': 'gbtree',
                  'objective': 'multi:softmax',
                  'num_class':4,
                  'seed': 2023,
                  'nthread': 3,
                  'silent': 1,
                  'enable_categorical': True}

        watchlist = [(dtrain, 'train')]

        bst = xg.train(params, dtrain, num_boost_round=1000, evals=watchlist)
        ypred = bst.predict(dtest)
        ypred = [round(value) for value in ypred]
        ypred = np.array(ypred).astype(int)

        pred_labs[test_index] = ypred

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

    return Mcc

def cv_cbt_cls(feas, labs, seed, platform):

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

    save_path = root + platform + '_' + 'cv_cbt_cls_res_labs_str'
    with open(save_path, 'a+', newline='') as ws:
        writer = csv.writer(ws)
        writer.writerows(np.column_stack((labs,pred_labs)))

    save_path = root + platform + '_' + 'cv_cbt_cls_res_folds_str'
    with open(save_path, 'a+', newline='') as ws:
        writer = csv.writer(ws)
        writer.writerows(out_res)
    return Mcc

def feature_importance_catboost(model):
    result=pd.DataFrame(model.get_feature_importance(),index=model.feature_names_,columns=['FeatureImportance'])
    return result.sort_values('FeatureImportance',ascending=False)


def test_CBT(fea_train, labs_train, fea_test, labs_test, seed, platform, header):
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
    save_path = root + platform + '_' + 'test_DT_cls_res'
    with open(save_path, 'a+', newline='') as ws:
        writer = csv.writer(ws)
        writer.writerow(out)
    #return Mcc
    fea_imp = feature_importance_catboost(bst)
    fea_imp.to_csv(root + platform + 'fea_importance.csv')
    # shap.initjs()

    # explainer = shap.TreeExplainer(bst)
    # shap_values = explainer.shap_values(Pool(train_X, train_y, cat_features=[0, 1, 2, 3]))
    # shap.summary_plot(shap_values[0], train_X)
    # shap.summary_plot(shap_values[1], train_X)
    # shap.summary_plot(shap_values[2], train_X)
    # shap.summary_plot(shap_values[3], train_X)
    # shap.plots.scatter(shap_values[0][:, 0], color=shap_values[0][:, 1])


    return Acc, F1, Mcc

def test_XG(fea_train, labs_train, fea_test, labs_test, param, seed, platform, header):
    seed_torch(seed)
    if param['objective'] == 0:
        object = 'binary:logistic'
    elif param['objective'] == 1:
        object = 'binary:hinge'

    params = {'booster': 'gbtree',
              'objective': object,
              'max_depth': param['md'],
              'lambda': param['lam'],
              'subsample': param['sub'],
              'colsample_bytree': param['cb'],
              'min_child_weight': param['mcw'],
              'eta': param['eta'],
              'seed': param['seed'],
              'nthread': 3,
              'silent': 1,
              'tree_method': 'auto',
              'enable_categorical': True}

    train_X, train_y = format_data_for_xg(fea_train, labs_train)
    test_X, test_y = format_data_for_xg(fea_test, labs_test)

    dtrain = xg.DMatrix(train_X, label=train_y, enable_categorical=True)
    dtest = xg.DMatrix(test_X, enable_categorical=True)

    watchlist = [(dtrain, 'train')]

    bst = xg.train(params, dtrain, num_boost_round=param['nbr'], evals=watchlist)
    ypred = bst.predict(dtest)
    ypred = [round(value) for value in ypred]
    ypred = np.array(ypred).astype(int)
    #print_decision_rules(rf, header[0], ['0','1','2','3'])
    #pred_l = rf.predict(fea_test)
    Acc = acc(test_y, ypred)
    F1 = f1(test_y, ypred, average='micro')
    rec = recall(test_y, ypred, average='micro')
    pre = precision(test_y, ypred, average='micro')
    Mcc = mcc(test_y, ypred)

    out = [Acc, F1, rec, pre, Mcc]
    print(out)
    # save_path = root + platform + '_' + 'test_DT_cls_res'
    # with open(save_path, 'a+', newline='') as ws:
    #     writer = csv.writer(ws)
    #     writer.writerow(out)
    #return Mcc
    ax = xg.plot_importance(bst, importance_type='weight')
    ax.figure.tight_layout()
    ax.figure.savefig(platform + '_feature_importance_weight.png')

    ax = xg.plot_importance(bst, importance_type='gain')
    ax.figure.tight_layout()
    ax.figure.savefig(platform + '_feature_importance_gain.png')

    ax = xg.plot_importance(bst, importance_type='cover')
    ax.figure.tight_layout()
    ax.figure.savefig(platform + '_feature_importance_cover.png')

    xg.plot_tree(bst)
    fig = matplotlib.pyplot.gcf()
    fig.set_size_inches(150, 100)
    fig.savefig(platform + '_learned_tree.png')
    return Acc, F1, Mcc

def print_decision_rules(rf, fea_name, class_name):
    from sklearn.tree import export_graphviz
    # Export as dot file
    for i in range(len(rf.estimators_)):
        estimator = rf.estimators_[i]
        export_graphviz(estimator, out_file='tree.dot',
                        feature_names=fea_name,
                        class_names=class_name,
                        rounded=True, proportion=False,
                        precision=2, filled=True)

        # Convert to png using system command (requires Graphviz)
        from subprocess import call
        call(['dot', '-Tpng', 'tree.dot', '-o', root + platform + '/tree' + str(i) + '.png', '-Gdpi=600'])

        # Display in jupyter notebook
        from IPython.display import Image
        Image(filename=root + platform + '/tree' + str(i) + '.png')


def test_RF(fea_train, labs_train, fea_test, labs_test, best_params, seed, platform, header):
    rf = RandomForestClassifier(n_estimators=best_params['nestimator'], max_depth=best_params['max_depth'], class_weight='balanced', random_state=seed)
    rf.fit(fea_train, labs_train)

    print(rf.feature_importances_)
    save_path = root + platform + 'DT_feature_importance.csv'
    # out = np.column_stack((hurdle_real, hurdle_real_ctrs))
    with open(save_path, 'a+', newline='') as ws:
        writer = csv.writer(ws)
        writer.writerow(header[0])
        writer.writerow(rf.feature_importances_)

    #print_decision_rules(rf, header[0], ['0','1','2','3'])
    pred_l = rf.predict(fea_test)
    Acc = acc(labs_test, pred_l)
    F1 = f1(labs_test, pred_l, average='micro')
    rec = recall(labs_test, pred_l, average='micro')
    pre = precision(labs_test, pred_l, average='micro')
    Mcc = mcc(labs_test, pred_l)

    out = [Acc, F1, rec, pre, Mcc]
    print(out)
    save_path = root + platform + '_' + 'test_DT_cls_res'
    with open(save_path, 'a+', newline='') as ws:
        writer = csv.writer(ws)
        writer.writerow(out)
    #return Mcc

    return Acc, F1, Mcc

def test_DT(fea_train, labs_train, fea_test, labs_test, seed, platform, header):
    #rf = RandomForestClassifier(n_estimators=best_params['nestimator'], max_depth=best_params['max_depth'], class_weight='balanced', random_state=seed)
    rf = DT(random_state=seed)
    rf.fit(fea_train, labs_train)
    tree_to_pseudo(rf, list(header[0]))
    print(rf.feature_importances_)
    save_path = root + platform + '_DT_feature_importance_str.csv'
    # out = np.column_stack((hurdle_real, hurdle_real_ctrs))
    with open(save_path, 'a+', newline='') as ws:
        writer = csv.writer(ws)
        writer.writerow(header[0])
        writer.writerow(rf.feature_importances_)
    #print_decision_rules(rf, header[0], ['0','1','2','3'])

    with open(root + platform + '_DT_tree_str.txt', "w") as fout:
        fout.write(tree.export_text(rf, feature_names=list(header[0])))

    plt.figure(figsize=(80,120))
    tree.plot_tree(rf, feature_names=list(header[0]), fontsize=10, filled=True)

    plt.title('Decision tree trained on all ' + platform + ' workflow')
    #plt.show()
    plt.savefig(root + platform + '_DT_tree_str.pdf', dpi=360)

    #save_path = root + platform + 'DT_trr_path.csv'
    # out = np.column_stack((hurdle_real, hurdle_real_ctrs))

    pred_l = rf.predict(fea_test)
    Acc = acc(labs_test, pred_l)
    F1 = f1(labs_test, pred_l, average='micro')
    rec = recall(labs_test, pred_l, average='micro')
    pre = precision(labs_test, pred_l, average='micro')
    Mcc = mcc(labs_test, pred_l)

    out = [Acc, F1, rec, pre, Mcc]
    print(out)
    save_path = root + platform + '_' + 'test_DT_cls_res_str'
    with open(save_path, 'a+', newline='') as ws:
        writer = csv.writer(ws)
        writer.writerow(out)
    #return Mcc

    return Acc, F1, Mcc

def objective(trail, feas, labs, seed, platform):
    seed_torch(seed)

    objective = trail.suggest_int('objective', 0, 0, step=1)
    md = trail.suggest_int('md', 1, 20)
    lam = trail.suggest_int('lam', 1, 800, step=10)

    sub = trail.suggest_float('sub', 0.01, 1, step=0.01)
    cb = trail.suggest_float('cb', 0.01, 1, step=0.01)
    mcw = trail.suggest_int('mcw', 0, 200, step=5)

    eta = trail.suggest_float('eta', 0, 1, step=0.01)

    nbr = trail.suggest_int('nbr', 50, 500, step=50)

    seed = trail.suggest_int('seed', 2023, 2023, step=1)

    param = {
        'objective': objective, 'md': md, 'lam': lam, 'sub': sub, 'cb': cb, 'mcw': mcw,
        'eta': eta, 'seed': 2023, 'nbr': nbr, 'seed': seed
    }
    mcc_test = cv_XG_cls(feas, labs, seed, platform, param)
    return mcc_test

def optuna_optimize(feas, labs, seed, platform):

    #save_sty = 'sqlite:///' + root_fold + 'op.db'
    study = optuna.create_study(study_name='test', direction='maximize', sampler=optuna.samplers.TPESampler(),
                                pruner=optuna.pruners.HyperbandPruner())#, storage=save_sty)

    func = lambda trial: objective(trial, feas, labs, seed, platform)
    study.optimize(func, n_trials=100)
    print(study.best_params)
    print(study.best_trial)
    print(study.best_trial.value)
    return study.best_params


def tree_to_pseudo(tree, feature_names):
    """
    Outputs a decision tree model as if/then pseudocode

    Parameters:
    -----------
    tree: decision tree model
        The decision tree to represent as pseudocode
    feature_names: list
        The feature names of the dataset used for building the decision tree
    """

    left = tree.tree_.children_left
    right = tree.tree_.children_right
    threshold = tree.tree_.threshold
    features = [feature_names[i] for i in tree.tree_.feature]
    value = tree.tree_.value

    def recurse(left, right, threshold, features, node, depth=0):
        indent = "  " * depth
        if (threshold[node] != -2):
            print(indent, "if ( " + features[node] + " <= " + str(threshold[node]) + " ) {")
            if left[node] != -1:
                recurse(left, right, threshold, features, left[node], depth + 1)
                print(indent, "} else {")
                if right[node] != -1:
                    recurse(left, right, threshold, features, right[node], depth + 1)
                print(indent, "}")
        else:
            print(indent, "return " + str(value[node]))

    recurse(left, right, threshold, features, 0)

if __name__ == '__main__':
    root = ''
    filename = root + 'ranking_fragpipe_mean.xlsx'
    sheet = 'ranking_all'
    platform = 'fragpipe'
    seed = 100
    ranks, OHfea, labs, header, feas_str = read_file(filename, sheet)
    #cv_XG_cls(OHfea, labs, seed, platform)
    cv_cbt_cls(feas_str, labs, seed, platform)
    test_CBT(feas_str, labs, feas_str, labs, seed, platform, header)

    filename = root + 'ranking_maxquant_mean.xlsx'
    sheet = 'ranking_all'
    platform = 'maxquant'

    ranks, OHfea, labs, header, feas_str = read_file(filename, sheet)
    #cv_XG_cls(feas_str, labs, seed, platform)
    cv_cbt_cls(feas_str, labs, seed, platform)
    test_CBT(feas_str, labs, feas_str, labs, seed, platform, header)

    filename = root + 'ranking_diann_mean.xlsx'
    sheet = 'ranking_all'
    platform = 'diann'
    ranks, OHfea, labs, header, feas_str = read_file(filename, sheet)
    cv_cbt_cls(feas_str, labs, seed, platform)
    test_CBT(feas_str, labs, feas_str, labs, seed, platform, header)


