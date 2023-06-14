import sys
from sklearn.metrics import roc_auc_score
import pandas as pd
import numpy as np
import os
import random
import csv
import time
from sklearn.metrics import mean_squared_error
from scipy import stats
from sklearn.metrics import accuracy_score as acc, recall_score as recall, precision_score as precision, f1_score as f1, matthews_corrcoef as mcc
from sklearn.metrics import confusion_matrix

def seed_torch(seed=1029):
    random.seed(seed)
    os.environ['PYTHONHASHSEED'] = str(seed)
    np.random.seed(seed)
    # torch.manual_seed(seed)
    # torch.cuda.manual_seed(seed)
    # torch.backends.cudnn.deterministic = True
    # torch.backends.cudnn.benchmark = False

seed_torch()
R_fold = 'benchmark_R_scripts/'
save_folder = 'ensemble_metrics_topK_mean_DDA/'

def cls_metric(adjp, label, logFC, q, lfc):
    y_true = label
    y_pred = (adjp <= q) & (abs(logFC) >= lfc)
    Acc = acc(y_true, y_pred)
    Prec = precision(y_true, y_pred)
    Rec = recall(y_true, y_pred)
    F1 = f1(y_true, y_pred)
    F1w = f1(y_true, y_pred, average='weighted')
    Mcc = mcc(y_true, y_pred)
    nMcc = (Mcc + 1) / 2
    if len(np.unique(label)) == 1:
        spec = 0
    else:
        tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()
        spec = tn / (tn + fp)
    #spec = tn / (tn + fp)
    geomean = np.sqrt(Rec * spec)
    return Acc, Prec, Rec, F1, F1w, Mcc, nMcc, geomean


def cal_auc(score, label, const, p, logFC, rea_logFC):
    label = np.array(label)
    score = np.array(score)
    score[np.where((np.isnan(score)))[0]] = 0
    const = np.array(const)

    aucs = []
    paucs1 = []
    paucs2 = []
    paucs3 = []
    sprs = []
    rmses = []
    lens = []

    Cls_met_1_001 = []
    Cls_met_1_005 = []

    cls_met_1_001 = cls_metric(1 - score, label, logFC, 0.01, 1)
    cls_met_1_005 = cls_metric(1 - score, label, logFC, 0.05, 1)

    Cls_met_1_001.extend(cls_met_1_001)
    Cls_met_1_005.extend(cls_met_1_005)

    if len(np.unique(label)) == 1:
        auc = 0
        pauc1 = 0
        pauc2 = 0
        pauc3 = 0
    else:
        auc = roc_auc_score(label, score)
        pauc1 = roc_auc_score(label, score, max_fpr=0.01)
        pauc2 = roc_auc_score(label, score, max_fpr=0.05)
        pauc3 = roc_auc_score(label, score, max_fpr=0.1)

    spr = stats.spearmanr(logFC, rea_logFC)[0]
    rmse = np.sqrt(mean_squared_error(logFC, rea_logFC))
    aucs.append(auc)
    paucs1.append(pauc1)
    paucs2.append(pauc2)
    paucs3.append(pauc3)
    sprs.append(spr)
    rmses.append(rmse)
    lens.append(len(label))

    uni_cons = np.unique(const)
    for con in uni_cons:
        idx = np.where((const == con))[0]

        cls_met_1_001 = cls_metric(1 - score[idx], label[idx], logFC[idx], 0.01, 1)
        cls_met_1_005 = cls_metric(1 - score[idx], label[idx], logFC[idx], 0.05, 1)

        Cls_met_1_001.extend(cls_met_1_001)
        Cls_met_1_005.extend(cls_met_1_005)

        if len(np.unique(label[idx])) == 1:
            auc = 0
            pauc1 = 0
            pauc2 = 0
            pauc3 = 0
        else:
            auc = roc_auc_score(label[idx], score[idx])
            pauc1 = roc_auc_score(label[idx], score[idx], max_fpr=0.01)
            pauc2 = roc_auc_score(label[idx], score[idx], max_fpr=0.05)
            pauc3 = roc_auc_score(label[idx], score[idx], max_fpr=0.1)
        aucs.append(auc)
        paucs1.append(pauc1)
        paucs2.append(pauc2)
        paucs3.append(pauc3)

        sprs.append(stats.spearmanr(logFC[idx], rea_logFC[idx])[0])
        rmses.append(np.sqrt(mean_squared_error(logFC[idx], rea_logFC[idx])))
        lens.append(len(idx))

    return aucs, paucs1, paucs2, paucs3, sprs, rmses, lens, Cls_met_1_001, Cls_met_1_005

def cal_auc1(score, label, const, p, logFC):
    label = np.array(label)
    score = np.array(score)
    score[np.where((np.isnan(score)))[0]] = 0
    const = np.array(const)

    aucs = []
    paucs = []
    sprs = []
    rmses = []
    lens = []

    Cls_met_1_001 = []
    Cls_met_1_005 = []

    cls_met_1_001 = cls_metric(1 - score, label, logFC, 0.01, 1)
    cls_met_1_005 = cls_metric(1 - score, label, logFC, 0.05, 1)

    Cls_met_1_001.extend(cls_met_1_001)
    Cls_met_1_005.extend(cls_met_1_005)

    if len(np.unique(label)) <= 1:
        auc = 0
        pauc = 0
    else:
        auc = roc_auc_score(label, score)
        pauc = roc_auc_score(label, score, max_fpr=p)

    aucs.append(auc)
    paucs.append(pauc)
    lens.append(len(label))

    uni_cons = np.unique(const)
    for con in uni_cons:
        idx = np.where((const == con))[0]

        cls_met_1_001 = cls_metric(1 - score[idx], label[idx], logFC[idx], 0.01, 1)
        cls_met_1_005 = cls_metric(1 - score[idx], label[idx], logFC[idx], 0.05, 1)

        Cls_met_1_001.extend(cls_met_1_001)
        Cls_met_1_005.extend(cls_met_1_005)

        if len(np.unique(label[idx])) == 1:
            auc = 0
            pauc = 0
        else:
            auc = roc_auc_score(label[idx], score[idx])
            pauc = roc_auc_score(label[idx], score[idx], max_fpr=p)
        aucs.append(auc)
        paucs.append(pauc)

        lens.append(len(idx))

    return aucs, paucs, sprs, rmses, lens, Cls_met_1_001, Cls_met_1_005

def get_true_logFC(label, contrast, dataset):
    Label = np.zeros((len(label), 1))
    Label[np.where((label))[0]] = 1
    label = Label
    true_logFC=np.zeros((len(label), 1))
    if dataset == 'CPTAC_DDA':
        true_logFC[np.where((label == 1) & (contrast == 'conditionB-conditionA'))[0], 0] = np.log2(0.74 / 0.25)
        true_logFC[np.where((label == 1) & (contrast == 'conditionC-conditionA'))[0], 0] = np.log2(2.2 / 0.25)
        true_logFC[np.where((label == 1) & (contrast == 'conditionD-conditionA'))[0], 0] = np.log2(6.7 / 0.25)
        true_logFC[np.where((label == 1) & (contrast == 'conditionE-conditionA'))[0], 0] = np.log2(20 / 0.25)
        true_logFC[np.where((label == 1) & (contrast == 'conditionC-conditionB'))[0], 0] = np.log2(2.2 / 0.74)
        true_logFC[np.where((label == 1) & (contrast == 'conditionD-conditionB'))[0], 0] = np.log2(6.7 / 0.74)
        true_logFC[np.where((label == 1) & (contrast == 'conditionE-conditionB'))[0], 0] = np.log2(20 / 0.74)
        true_logFC[np.where((label == 1) & (contrast == 'conditionD-conditionC'))[0], 0] = np.log2(6.7 / 2.2)
        true_logFC[np.where((label == 1) & (contrast == 'conditionE-conditionC'))[0], 0] = np.log2(20 / 2.2)
        true_logFC[np.where((label == 1) & (contrast == 'conditionE-conditionD'))[0], 0] = np.log2(20 / 6.7)
    elif dataset == 'yeast_DDA':
        true_logFC[np.where((label == 1) & (contrast == 'conditionB-conditionA'))[0], 0] = np.log2(4 / 2)
        true_logFC[np.where((label == 1) & (contrast == 'conditionC-conditionA'))[0], 0] = np.log2(10 / 2)
        true_logFC[np.where((label == 1) & (contrast == 'conditionD-conditionA'))[0], 0] = np.log2(25 / 2)
        true_logFC[np.where((label == 1) & (contrast == 'conditionE-conditionA'))[0], 0] = np.log2(50 / 2)
        true_logFC[np.where((label == 1) & (contrast == 'conditionC-conditionB'))[0], 0] = np.log2(10 / 4)
        true_logFC[np.where((label == 1) & (contrast == 'conditionD-conditionB'))[0], 0] = np.log2(25 / 4)
        true_logFC[np.where((label == 1) & (contrast == 'conditionE-conditionB'))[0], 0] = np.log2(50 / 4)
        true_logFC[np.where((label == 1) & (contrast == 'conditionD-conditionC'))[0], 0] = np.log2(25 / 10)
        true_logFC[np.where((label == 1) & (contrast == 'conditionE-conditionC'))[0], 0] = np.log2(50 / 10)
        true_logFC[np.where((label == 1) & (contrast == 'conditionE-conditionD'))[0], 0] = np.log2(50 / 25)
    elif dataset == 'yeast1819_DDA':
        true_logFC[np.where((label == 1) & (contrast == 'conditionB-conditionA'))[0], 0] = np.log2(125 / 50)
        true_logFC[np.where((label == 1) & (contrast == 'conditionC-conditionA'))[0], 0] = np.log2(250 / 50)
        true_logFC[np.where((label == 1) & (contrast == 'conditionD-conditionA'))[0], 0] = np.log2(500 / 50)
        true_logFC[np.where((label == 1) & (contrast == 'conditionE-conditionA'))[0], 0] = np.log2(2500 / 50)
        true_logFC[np.where((label == 1) & (contrast == 'conditionF-conditionA'))[0], 0] = np.log2(5000 / 50)
        true_logFC[np.where((label == 1) & (contrast == 'conditionG-conditionA'))[0], 0] = np.log2(12500 / 50)
        true_logFC[np.where((label == 1) & (contrast == 'conditionH-conditionA'))[0], 0] = np.log2(25000 / 50)
        true_logFC[np.where((label == 1) & (contrast == 'conditionI-conditionA'))[0], 0] = np.log2(50000 / 50)
        true_logFC[np.where((label == 1) & (contrast == 'conditionC-conditionB'))[0], 0] = np.log2(250 / 125)
        true_logFC[np.where((label == 1) & (contrast == 'conditionD-conditionB'))[0], 0] = np.log2(500 / 125)
        true_logFC[np.where((label == 1) & (contrast == 'conditionE-conditionB'))[0], 0] = np.log2(2500 / 125)
        true_logFC[np.where((label == 1) & (contrast == 'conditionF-conditionB'))[0], 0] = np.log2(5000 / 125)
        true_logFC[np.where((label == 1) & (contrast == 'conditionG-conditionB'))[0], 0] = np.log2(12500 / 125)
        true_logFC[np.where((label == 1) & (contrast == 'conditionH-conditionB'))[0], 0] = np.log2(25000 / 125)
        true_logFC[np.where((label == 1) & (contrast == 'conditionI-conditionB'))[0], 0] = np.log2(50000 / 125)
        true_logFC[np.where((label == 1) & (contrast == 'conditionD-conditionC'))[0], 0] = np.log2(500 / 250)
        true_logFC[np.where((label == 1) & (contrast == 'conditionE-conditionC'))[0], 0] = np.log2(2500 / 250)
        true_logFC[np.where((label == 1) & (contrast == 'conditionF-conditionC'))[0], 0] = np.log2(5000 / 250)
        true_logFC[np.where((label == 1) & (contrast == 'conditionG-conditionC'))[0], 0] = np.log2(12500 / 250)
        true_logFC[np.where((label == 1) & (contrast == 'conditionH-conditionC'))[0], 0] = np.log2(25000 / 250)
        true_logFC[np.where((label == 1) & (contrast == 'conditionI-conditionC'))[0], 0] = np.log2(50000 / 250)
        true_logFC[np.where((label == 1) & (contrast == 'conditionE-conditionD'))[0], 0] = np.log2(2500 / 500)
        true_logFC[np.where((label == 1) & (contrast == 'conditionF-conditionD'))[0], 0] = np.log2(5000 / 500)
        true_logFC[np.where((label == 1) & (contrast == 'conditionG-conditionD'))[0], 0] = np.log2(12500 / 500)
        true_logFC[np.where((label == 1) & (contrast == 'conditionH-conditionD'))[0], 0] = np.log2(25000 / 500)
        true_logFC[np.where((label == 1) & (contrast == 'conditionI-conditionD'))[0], 0] = np.log2(50000 / 500)
        true_logFC[np.where((label == 1) & (contrast == 'conditionF-conditionE'))[0], 0] = np.log2(5000 / 2500)
        true_logFC[np.where((label == 1) & (contrast == 'conditionG-conditionE'))[0], 0] = np.log2(12500 / 2500)
        true_logFC[np.where((label == 1) & (contrast == 'conditionH-conditionE'))[0], 0] = np.log2(25000 / 2500)
        true_logFC[np.where((label == 1) & (contrast == 'conditionI-conditionE'))[0], 0] = np.log2(50000 / 2500)
        true_logFC[np.where((label == 1) & (contrast == 'conditionG-conditionF'))[0], 0] = np.log2(12500 / 5000)
        true_logFC[np.where((label == 1) & (contrast == 'conditionH-conditionF'))[0], 0] = np.log2(25000 / 5000)
        true_logFC[np.where((label == 1) & (contrast == 'conditionI-conditionF'))[0], 0] = np.log2(50000 / 5000)
        true_logFC[np.where((label == 1) & (contrast == 'conditionH-conditionG'))[0], 0] = np.log2(25000 / 12500)
        true_logFC[np.where((label == 1) & (contrast == 'conditionI-conditionG'))[0], 0] = np.log2(50000 / 12500)
        true_logFC[np.where((label == 1) & (contrast == 'conditionI-conditionH'))[0], 0] = np.log2(50000 / 25000)
    elif dataset == 'human_ecoli':
        true_logFC[np.where((label == 1) & (contrast == 'conditionB-conditionA'))[0], 0] = np.log2(2 / 1)
    return true_logFC

def get_uniport_id(protein_list):
    pro_uniport = []
    for i in protein_list:
        pros = i.split(';')
        pro =  ''
        for j in pros:

            proj = j.split('|')
            if proj[0] == 'sp':
                pro = pro + proj[1]
            else:
                pro = pro + proj[0]
        pro_uniport.append(pro)
    return pro_uniport

def get_candidate_proteins(path, acq, in_type, dataset, flag):
    save_path = save_folder + acq + '_' + in_type + '_' + dataset + '_' + str(flag) + '.csv'
    os.system('Rscript ' + R_fold + 'get_candidate_proteins.R ' + in_type + ' ' + path + ' ' + save_path)
    candidate_pros = pd.read_csv(save_path, header=None, sep=',').values
    return candidate_pros

def get_label_ups(proteins, key):
    label = [key in proteins[i] for i in range(len(proteins))]
    label = np.array(label)

    return label

def get_label_AK(proteins, in_files):
    labeling = pd.read_csv(in_files['labels'], sep='\t', header=0)
    trues = labeling['Uniprot_Acc'][0:100]

    all_DEs = []

    for true in trues:
        DE = true.split(';')
        for de in DE:
            de0 = de.split('-')
            all_DEs.append(de0[0])

    label = []
    for pro in proteins:
        pros = pro.split(';')
        l = 0
        for L in pros:
            if L in all_DEs:
                l += 1
            elif '|' in L:
                if L.split('|')[1] in all_DEs:
                    l += 1
        if l == len(pros):
            label.append('TRUE')
        else:
            label.append('FALSE')

    Label = ['TRUE' in label[i] for i in range(len(label))]
    Label = np.array(Label)
    return Label

def get_comparable_res(proteins, score, label, const, logFC, real_logFC, in_files, dataset):

    proteins_all = in_files['all_candidate'][:, 0]
    pro_uniport_all = get_uniport_id(proteins_all)
    uni_cons = np.unique(const)
    proteins_c = []
    score_c = []
    label_c = []
    const_c = []
    logFC_c = []
    real_logFC_c = []
    for con in uni_cons:
        idx_con = np.where((const==con))[0]
        pro_con = proteins[idx_con]
        pro_con_uniport = get_uniport_id(pro_con)
        non_pro = np.setdiff1d(pro_uniport_all, pro_con_uniport)

        comm, idx1, idx2 = np.intersect1d(non_pro, pro_uniport_all, return_indices=True)

        proteins_c.extend(proteins[idx_con])
        proteins_c.extend(proteins_all[idx2])
        score_c.extend(score[idx_con])
        score_c.extend([0 for i in range(len(non_pro))])
        label_c.extend(label[idx_con])
        label_c.extend(get_label_ups(proteins_all[idx2], '_HUMAN'))
        const_c.extend(const[idx_con])
        const_c.extend([con for i in range(len(non_pro))])
        logFC_c.extend(logFC[idx_con])
        logFC_c.extend([0 for i in range(len(non_pro))])

    real_logFC_c = get_true_logFC(label_c, const_c, dataset)
    return np.array(proteins_c), np.array(score_c), np.array(label_c), np.array(const_c), np.array(logFC_c), np.array(real_logFC_c)

def cal_mse_auc(tool, acq, in_type, inten_ty, imput, normal, res_file, dataset, run_time, in_files, mode, k):
    if(os.path.exists(res_file)):
        res = pd.read_csv(res_file, sep=',', header=0)
        idx = np.where(((~np.isnan(np.array(res['logFC_avg']))) & (~np.isnan(np.array(res['pensemble']))) & (~np.isnan(np.array(res['qensemble']))) &
                        (~np.isinf(np.array(res['logFC_avg']))) & (~np.isinf(np.array(res['pensemble']))) & (~np.isinf(np.array(res['qensemble'])))))[0]

        proteins = np.array(res['protein'])[idx]
        score = 1 - np.array(res['qensemble'])[idx]

        label = get_label_ups(proteins, '_HUMAN')
        const = np.array(res['contrast'])[idx]
        logFC = np.array(res['logFC_avg'])[idx]
        real_logFC = get_true_logFC(label, const, dataset)

        protein_c, score_c, label_c, const_c, logFC_c, real_logFC_c = get_comparable_res(proteins, score, label, const,
                                                                                         logFC, real_logFC, in_files, dataset)

        p = 0.01
        aucs, paucs1, paucs2, paucs3, sprs, rmses, lens, Cls_met_1_001, Cls_met_1_005 = cal_auc(score, label, const, p,
                                                                                                logFC, real_logFC)
        out1 = [tool, in_type, inten_ty, imput, normal, run_time, k]
        out2 = [tool, in_type, inten_ty, imput, normal, k]
        out3 = [tool, in_type, inten_ty, imput, normal, run_time, k]
        out4 = [tool, in_type, inten_ty, imput, normal, k]

        out1.extend(lens)
        out1.extend(aucs)
        out1.extend(sprs)
        out1.extend(rmses)
        out2.extend(paucs1)
        out2.extend(paucs2)
        out2.extend(paucs3)
        out3.extend(Cls_met_1_001)
        out4.extend(Cls_met_1_005)

        aucs, paucs1, paucs2, paucs3, sprs, rmses, lens, Cls_met_1_001, Cls_met_1_005 = cal_auc(score_c, label_c,
                                                                                                const_c, p, logFC_c,
                                                                                                real_logFC_c)
        out5 = [tool, in_type, inten_ty, imput, normal, run_time, k]
        out6 = [tool, in_type, inten_ty, imput, normal, k]
        out7 = [tool, in_type, inten_ty, imput, normal, run_time, k]
        out8 = [tool, in_type, inten_ty, imput, normal, k]
        out5.extend(lens)
        out5.extend(aucs)
        out5.extend(sprs)
        out5.extend(rmses)
        out6.extend(paucs1)
        out6.extend(paucs2)
        out6.extend(paucs3)
        out7.extend(Cls_met_1_001)
        out8.extend(Cls_met_1_005)
    else:
        out1 = [tool, in_type, inten_ty, imput, normal, run_time, k]
        out2 = [tool, in_type, inten_ty, imput, normal, k]
        out3 = [tool, in_type, inten_ty, imput, normal, run_time, k]
        out4 = [tool, in_type, inten_ty, imput, normal, k]
        out1.extend(['bug'])
        out2.extend(['bug'])
        out3.extend(['bug'])
        out4.extend(['bug'])

        out5 = [tool, in_type, inten_ty, imput, normal, run_time, k]
        out6 = [tool, in_type, inten_ty, imput, normal, k]
        out7 = [tool, in_type, inten_ty, imput, normal, run_time, k]
        out8 = [tool, in_type, inten_ty, imput, normal, k]

        out5.extend(['bug'])
        out6.extend(['bug'])
        out7.extend(['bug'])
        out8.extend(['bug'])

    save_path = save_folder + dataset + '_' + acq + '_' + mode + '_' + in_type + '_aucs1.csv'
    with open(save_path, 'a+', newline='') as ws:
        writer = csv.writer(ws)
        writer.writerow(out1)

    save_path = save_folder + dataset + '_' + acq + '_' + mode + '_' + in_type + '_paucs1.csv'
    with open(save_path, 'a+', newline='') as ws:
        writer = csv.writer(ws)
        writer.writerow(out2)

    save_path = save_folder + dataset + '_' + acq + '_' + mode + '_' + in_type + '_cls_001_1.csv'
    with open(save_path, 'a+', newline='') as ws:
        writer = csv.writer(ws)
        writer.writerow(out3)

    save_path = save_folder + dataset + '_' + acq + '_' + mode + '_' + in_type + '_cls_005_1.csv'
    with open(save_path, 'a+', newline='') as ws:
        writer = csv.writer(ws)
        writer.writerow(out4)

    save_path = save_folder + dataset + '_' + acq + '_' + mode + '_' + in_type + '_aucs2.csv'
    with open(save_path, 'a+', newline='') as ws:
        writer = csv.writer(ws)
        writer.writerow(out5)

    save_path = save_folder + dataset + '_' + acq + '_' + mode + '_' + in_type + '_paucs2.csv'
    with open(save_path, 'a+', newline='') as ws:
        writer = csv.writer(ws)
        writer.writerow(out6)

    save_path = save_folder + dataset + '_' + acq  +  '_' + mode + '_' + in_type + '_cls_001_2.csv'
    with open(save_path, 'a+', newline='') as ws:
        writer = csv.writer(ws)
        writer.writerow(out7)

    save_path = save_folder + dataset + '_' + acq + '_' + mode + '_' + in_type + '_cls_005_2.csv'
    with open(save_path, 'a+', newline='') as ws:
        writer = csv.writer(ws)
        writer.writerow(out8)

def cal_mse_auc1(tool, acq, in_type, inten_ty, imput, normal, res_file, dataset, run_time, in_files, comb_idx, k):
    if(os.path.exists(res_file)):
        res = pd.read_csv(res_file, sep=',', header=0)
        idx = np.where(((~np.isnan(np.array(res['logFC_avg']))) & (~np.isnan(np.array(res['pensemble']))) & (~np.isnan(np.array(res['qensemble']))) &
                        (~np.isinf(np.array(res['logFC_avg']))) & (~np.isinf(np.array(res['pensemble']))) & (~np.isinf(np.array(res['qensemble'])))))[0]

        proteins = np.array(res['protein'])[idx]
        score = 1 - np.array(res['qensemble'])[idx]
        #label = get_label_AK(proteins, in_files)
        label = get_label_ups(proteins, '_ECOLI')
        const = np.array(res['contrast'])[idx]
        logFC = np.array(res['logFC_avg'])[idx]
        real_logFC = get_true_logFC(label, const, dataset)

        protein_c, score_c, label_c, const_c, logFC_c, real_logFC_c = get_comparable_res(proteins, score, label, const, logFC,real_logFC, in_files, dataset)

        p = 0.01
        aucs, paucs1, paucs2, paucs3, sprs, rmses, lens, Cls_met_1_001, Cls_met_1_005 = cal_auc(score, label, const, p, logFC, real_logFC)
        out1 = [tool, in_type, inten_ty, imput, normal, comb_idx, run_time, k]
        out2 = [tool, in_type, inten_ty, imput, normal, comb_idx, k]
        out3 = [tool, in_type, inten_ty, imput, normal, comb_idx, run_time, k]
        out4 = [tool, in_type, inten_ty, imput, normal, comb_idx, k]

        out1.extend(lens)
        out1.extend(aucs)
        out1.extend(sprs)
        out1.extend(rmses)
        out2.extend(paucs1)
        out2.extend(paucs2)
        out2.extend(paucs3)
        out3.extend(Cls_met_1_001)
        out4.extend(Cls_met_1_005)

        aucs, paucs1,  paucs2, paucs3, sprs, rmses, lens, Cls_met_1_001, Cls_met_1_005 = cal_auc(score_c, label_c, const_c, p, logFC_c,
                                                                               real_logFC_c)
        out5 = [tool, in_type, inten_ty, imput, normal, comb_idx, run_time, k]
        out6 = [tool, in_type, inten_ty, imput, normal, comb_idx, k]
        out7 = [tool, in_type, inten_ty, imput, normal, comb_idx, run_time, k]
        out8 = [tool, in_type, inten_ty, imput, normal, comb_idx, k]
        out5.extend(lens)
        out5.extend(aucs)
        out5.extend(sprs)
        out5.extend(rmses)
        out6.extend(paucs1)
        out6.extend(paucs2)
        out6.extend(paucs3)
        out7.extend(Cls_met_1_001)
        out8.extend(Cls_met_1_005)
    else:
        out1 = [tool, in_type, inten_ty, imput, normal, run_time, k]
        out2 = [tool, in_type, inten_ty, imput, normal, k]
        out3 = [tool, in_type, inten_ty, imput, normal, run_time, k]
        out4 = [tool, in_type, inten_ty, imput, normal, k]
        out1.extend(['bug'])
        out2.extend(['bug'])
        out3.extend(['bug'])
        out4.extend(['bug'])

        out5 = [tool, in_type, inten_ty, imput, normal, run_time, k]
        out6 = [tool, in_type, inten_ty, imput, normal, k]
        out7 = [tool, in_type, inten_ty, imput, normal, run_time, k]
        out8 = [tool, in_type, inten_ty, imput, normal, k]

        out5.extend(['bug'])
        out6.extend(['bug'])
        out7.extend(['bug'])
        out8.extend(['bug'])

    save_path = save_folder + dataset + '_' + acq + '_' + in_type + '_aucs1.csv'
    with open(save_path, 'a+', newline='') as ws:
        writer = csv.writer(ws)
        writer.writerow(out1)

    save_path = save_folder + dataset + '_' + acq + '_' + in_type + '_paucs1.csv'
    with open(save_path, 'a+', newline='') as ws:
        writer = csv.writer(ws)
        writer.writerow(out2)

    save_path = save_folder + dataset + '_' + acq + '_' + in_type + '_cls_001_1.csv'
    with open(save_path, 'a+', newline='') as ws:
        writer = csv.writer(ws)
        writer.writerow(out3)

    save_path = save_folder + dataset + '_' + acq + '_' + in_type + '_cls_005_1.csv'
    with open(save_path, 'a+', newline='') as ws:
        writer = csv.writer(ws)
        writer.writerow(out4)

    save_path = save_folder + dataset + '_' + acq + '_' + in_type + '_aucs2.csv'
    with open(save_path, 'a+', newline='') as ws:
        writer = csv.writer(ws)
        writer.writerow(out5)

    save_path = save_folder + dataset + '_' + acq + '_' + in_type + '_paucs2.csv'
    with open(save_path, 'a+', newline='') as ws:
        writer = csv.writer(ws)
        writer.writerow(out6)

    save_path = save_folder + dataset + '_' + acq + '_' + in_type + '_cls_001_2.csv'
    with open(save_path, 'a+', newline='') as ws:
        writer = csv.writer(ws)
        writer.writerow(out7)

    save_path = save_folder + dataset + '_' + acq + '_' + in_type + '_cls_005_2.csv'
    with open(save_path, 'a+', newline='') as ws:
        writer = csv.writer(ws)
        writer.writerow(out8)

def ensemble_DEA_threeinputs(dataset, platform, method, mode, in_files):
    '''
    :param dataset: CPTAC_DDA, yeast_DDA, yeast1819_DDA, human_ecoli_DDA,
    :param platform:
    :param method:
    :param save_name:
    :param mode:
    :return:
    '''
    file_root = '/data/res_DE_files/'
    if dataset == 'human_ecoli_DDA':
        if platform == 'fragpipe':
            wf1 = file_root + dataset + '/plgem/' + str(mode) + '/plgem_fragpipe__QRILC_div.mean.csv '  # plgem_maxquant__knn_sum
            wf2 = file_root + dataset + '/DEP/' + str(mode) + '/DEP_fragpipe_MaxLFQ_MinProb_.csv '
            wf3 = file_root + dataset + '/ProteoMM/' + str(mode) + '/ProteoMM_fragpipe_MaxLFQ_MinDet_max.csv '
        elif platform == 'maxquant':
            wf1 = file_root + dataset + '/plgem/' + str(mode) + '/plgem_maxquant__MinDet_.csv '  # plgem_maxquant__knn_sum
            wf2 = file_root + dataset + '/DEqMS/' + str(mode) + '/DEqMS_maxquant__MLE_center.median.csv '
            wf3 = file_root + dataset + '/ProteoMM/' + str(mode) + '/ProteoMM_maxquant__bpca_vsn.csv '
    else:
        if platform == 'fragpipe':
            wf1 = file_root + dataset + '/' + mode + '/plgem/plgem_fragpipe__QRILC_div.mean.csv ' #plgem_maxquant__knn_sum
            wf2 = file_root + dataset + '/' + mode + '/DEP/DEP_fragpipe_MaxLFQ_MinProb_.csv '
            wf3 = file_root + dataset + '/' + mode + '/ProteoMM/ProteoMM_fragpipe_MaxLFQ_MinDet_max.csv '
        elif platform == 'maxquant':
            wf1 = file_root + dataset + '/' + mode + '/plgem/plgem_maxquant__MinDet_.csv '  # plgem_maxquant__knn_sum
            wf2 = file_root + dataset + '/' + mode + '/DEqMS/DEqMS_maxquant__MLE_center.median.csv '
            wf3 = file_root + dataset + '/' + mode + '/ProteoMM/ProteoMM_maxquant__bpca_vsn.csv '

    if method == 'set':
        for operation in ['min', 'max', 'median']:
            save_fold = file_root + '/ensemle_res/' + dataset + '_' + str(mode) + '_' + platform + '_' + method + '_' + operation + '.csv'
            DEA_res_files = save_fold + ' ' + operation + ' ' + wf1 + wf2 + wf3
            print('Rscript ' + R_fold + 'ensemble_set.R ' + DEA_res_files)
            os.system('Rscript ' + R_fold + 'ensemble_set.R ' + DEA_res_files)
            if dataset == 'human_ecoli_DDA':
                cal_mse_auc1('ensemble_' + method + '_' + operation, 'DDA', platform, '_', '_', '_', save_fold, dataset,
                            0, in_files, mode, 3)
            else:
                cal_mse_auc('ensemble_'+method + '_' + operation, 'DDA', platform, '_', '_', '_', save_fold, dataset, 0, in_files,mode, 3)
    elif method == 'fisher':
        save_fold = file_root + '/ensemle_res/' + dataset + '_' + str(mode) + '_' + platform + '_' + method + '_' + '.csv'
        DEA_res_files = save_fold + ' ' + wf1 + wf2 + wf3
        print('Rscript ' + R_fold + 'ensemble_fisher.R ' + DEA_res_files)
        os.system('Rscript ' + R_fold + 'ensemble_fisher.R ' + DEA_res_files)
        if dataset == 'human_ecoli_DDA':
            cal_mse_auc1('ensemble_' + method, 'DDA', platform, '_', '_', '_', save_fold, dataset, 0,
                        in_files, mode, 3)
        else:
            cal_mse_auc('ensemble_' + method, 'DDA', platform, '_', '_', '_', save_fold, dataset, 0,
                    in_files, mode, 3)
    elif method == 'hurdle':
        save_fold = file_root + '/ensemle_res/' + dataset + '_' + str(mode) + '_' + platform + '_' + method + '_' + '.csv'
        DEA_res_files = save_fold + ' ' + wf1 + wf2 + wf3
        print('Rscript ' + R_fold + 'ensemble_hurdle.R ' + DEA_res_files)
        os.system('Rscript ' + R_fold + 'ensemble_hurdle.R ' + DEA_res_files)
        if dataset == 'human_ecoli_DDA':
            cal_mse_auc1('ensemble_' + method, 'DDA', platform, '_', '_', '_', save_fold, dataset, 0,
                        in_files, mode, 3)
        else:
            cal_mse_auc('ensemble_' + method, 'DDA', platform, '_', '_', '_', save_fold, dataset, 0,
                    in_files, mode, 3)

def getwfs(res_file, K, dataset, mode):
    file_root = '/data/res_DE_files/'
    res = pd.read_csv(res_file, sep=',', header=0).values
    wfs = ''
    for i in range(K):
        wf_s = res[i, 1].split('|')
        wf_r = res[i, 1].replace('|', '_')
        if (dataset == 'human_ecoli_DDA') :
            wf = file_root + dataset + '/' + wf_s[0] + '/' + str(mode) + '/' + wf_r + '.csv ' #
        elif (dataset =='human_ecoli_DIA'):
            wf = file_root + 'benchmark_res/' + dataset + '/' + wf_s[0] + '/' + str(mode) + '/' + wf_r + '.csv '  # /data/res_DE_files/benchmark_res/human_ecoli_DIA/ANOVA/1/ANOVA_DIANN_MaxLFQ__.csv
        elif dataset == 'ecoli_DIA':
            wf = file_root + 'benchmark_res/' + dataset + '/' + str(mode) + '/' + wf_s[0] + '/' + wf_r + '.csv '  # /data/res_DE_files/benchmark_res/ecoli_DIA/Narrow/ANOVA/ANOVA_DIANN_MaxLFQ__.csv
        else:
            wf = file_root + dataset + '/' + str(mode) + '/' + wf_s[0] + '/' + wf_r + '.csv ' #/data/res_DE_files/CPTAC_DDA/LTQ86/ANOVA/ANOVA_maxquant_LFQ_bpca_.csv

        wfs = wfs + wf
    return wfs

def ensemble_DEA_topKwf(dataset, platform, method, mode, in_files, wfs, acq, k):
    '''
    :param dataset: CPTAC_DDA, yeast_DDA, yeast1819_DDA, human_ecoli_DDA,
    :param platform:
    :param method:
    :param save_name:
    :param mode:
    :return:
    '''
    file_root = '/data/res_DE_files/'

    if method == 'set':
        for operation in ['min', 'max', 'median']:
            save_fold = file_root + '/ensemble_res_DIA_mean/' + dataset + '_' + str(mode) + '_' + platform + '_' + method + '_' + operation + '.csv'
            DEA_res_files = save_fold + ' ' + operation + ' ' + wfs
            print('Rscript ' + R_fold + 'ensemble_set.R ' + DEA_res_files)
            os.system('Rscript ' + R_fold + 'ensemble_set.R ' + DEA_res_files)
            if (dataset == 'human_ecoli_DDA') | (dataset == 'human_ecoli_DIA'):
                cal_mse_auc1('ensemble_' + method + '_' + operation, acq, platform, '_', '_', '_', save_fold, dataset,
                            0, in_files, mode, k)
            else:
                cal_mse_auc('ensemble_'+method + '_' + operation, acq, platform, '_', '_', '_', save_fold, dataset, 0, in_files, mode, k)
    elif method == 'fisher':
        save_fold = file_root + '/ensemble_res_DIA_mean/' + dataset + '_' + str(mode) + '_' + platform + '_' + method + '_' + '.csv'
        DEA_res_files = save_fold + ' ' + wfs
        print('Rscript ' + R_fold + 'ensemble_fisher.R ' + DEA_res_files)
        os.system('Rscript ' + R_fold + 'ensemble_fisher.R ' + DEA_res_files)
        if (dataset == 'human_ecoli_DDA') | (dataset == 'human_ecoli_DIA'):
            cal_mse_auc1('ensemble_' + method, acq, platform, '_', '_', '_', save_fold, dataset, 0,
                        in_files, mode, k)
        else:
            cal_mse_auc('ensemble_' + method, acq, platform, '_', '_', '_', save_fold, dataset, 0,
                    in_files, mode, k)
    elif method == 'hurdle':
        save_fold = file_root + '/ensemble_res_DIA_mean/' + dataset + '_' + str(mode) + '_' + platform + '_' + method + '_' + '.csv'
        DEA_res_files = save_fold + ' ' + wfs
        print('Rscript ' + R_fold + 'ensemble_hurdle.R ' + DEA_res_files)
        os.system('Rscript ' + R_fold + 'ensemble_hurdle.R ' + DEA_res_files)
        if (dataset == 'human_ecoli_DDA') | (dataset == 'human_ecoli_DIA'):
            cal_mse_auc1('ensemble_' + method, acq, platform, '_', '_', '_', save_fold, dataset, 0,
                        in_files, mode, k)
        else:
            cal_mse_auc('ensemble_' + method, acq, platform, '_', '_', '_', save_fold, dataset, 0,
                    in_files, mode, k)

if __name__ == '__main__':

    dataset = 'CPTAC_DDA'  ## A549-K562
    modes1 = ['LTQ86', 'LTQO65', 'LTQP65', 'LTQW56']
    modes2 = ['LTQ86', 'LTQO65', 'LTQW56']  # , 'LTQW56']
    #in_type = 'fragpipe'
    for in_type in ['fragpipe', 'maxquant']:
        acq = 'DDA'
        if in_type == 'fragpipe':
            modes = modes1
        else:
            modes = modes2
        for mode in modes:
            for k in [2, 3, 5, 7, 9, 11, 15, 20]:
                if acq == 'DIA':
                    if in_type == 'maxquant':
                        in_files = {'protein': '/home/hui/Documents/VM_share/DIA_data/A549_K562/DIA/maxquant/txt/proteinGroups.txt',
                                    'peptide': '/home/hui/Documents/VM_share/DIA_data/A549_K562/DIA/maxquant/txt/peptides.txt',
                                    'design': '/home/hui/Documents/VM_share/DIA_data/A549_K562/design_cell_line.txt',
                                    'msstats': '/home/hui/Documents/VM_share/DIA_data/A549_K562/DIA/maxquant/txt/evidence.txt',
                                    'save': '/home/hui/PycharmProjects/DL_proteomics/benchmark_res/',
                                    'design_stat': '/home/hui/Documents/VM_share/DIA_data/A549_K562/design_cell_line_msstats_DIA.txt',
                                    'labels': '/home/hui/Documents/VM_share/DIA_data/A549_K562/CL_labeling.txt'}
                    elif in_type == 'DIANN':
                        in_files = {
                            'protein': '/home/hui/Documents/VM_share/DIA_data/ecoli/DIANN/' + mode + '/report.pg_matrix.tsv',
                            'protein1': '/home/hui/Documents/VM_share/DIA_data/ecoli/DIANN/' + mode + '/report.tsv',
                            'peptide': '/home/hui/Documents/VM_share/DIA_data/yeast/frager/combined_peptide.tsv',
                            'design1': '/home/hui/Documents/VM_share/DIA_data/ecoli/DIANN/design_ecoli_' + mode + '.txt',
                            'design': '/home/hui/Documents/VM_share/DIA_data/A549_K562/design_cell_line.txt',
                            'msstats': '/home/hui/Documents/VM_share/DIA_data/A549_K562/MSstats.csv',
                            'save': '/home/hui/PycharmProjects/DL_proteomics/benchmark_res/',
                            'design_stat': '/home/hui/Documents/VM_share/DIA_data/yeast/design_yeast_msstats.txt',
                            'labels': '/home/hui/Documents/VM_share/DIA_data/A549_K562/CL_labeling.txt'
                            }
                elif acq == 'DDA':
                    if in_type == 'maxquant':
                        in_files = {
                            'protein': '/home/hui/Documents/VM_share/DIA_data/CPTAC/' + mode + '/combined/txt/proteinGroups.txt',
                            'peptide': '/home/hui/Documents/VM_share/DIA_data/CPTAC/' + mode + '/combined/txt/peptides.txt',
                            'design': '/home/hui/Documents/VM_share/DIA_data/CPTAC/design_CPTAC_mq.txt',
                            'msstats': '/home/hui/Documents/VM_share/DIA_data/CPTAC/' + mode + '/combined/txt/evidence.txt',
                            'save': '/home/hui/PycharmProjects/DL_proteomics/benchmark_res/',
                            'design_stat': '/home/hui/Documents/VM_share/DIA_data/CPTAC/design_' + mode + '_msstats.txt',
                            'labels': '/home/hui/Documents/VM_share/DIA_data/A549_K562/CL_labeling.txt',
                            'database': '/home/hui/Documents/VM_share/DIA_data/yeast_UPS_2017-09-14.fasta',
                            'pro_path': '/home/hui/Documents/VM_share/DIA_data/CPTAC/' + mode + '/combined/txt/'}
                    elif in_type == 'fragpipe':
                        in_files = {
                            'protein': '/home/hui/Documents/VM_share/DIA_data/CPTAC/fragpipe/' + mode + '/combined_protein.tsv',
                            'peptide': '/home/hui/Documents/VM_share/DIA_data/CPTAC/fragpipe/' + mode + '/combined_peptide.tsv',
                            'design': '/home/hui/Documents/VM_share/DIA_data/CPTAC/design_CPTAC_frag.txt',
                            'msstats': '/home/hui/Documents/VM_share/DIA_data/CPTAC/fragpipe/' + mode + '/MSstats.csv',
                            'save': '/home/hui/PycharmProjects/DL_proteomics/benchmark_res/',
                            'design_stat': '/home/hui/Documents/VM_share/DIA_data/CPTAC/design_' + mode + '_msstats.txt',
                            'labels': '/home/hui/Documents/VM_share/DIA_data/A549_K562/CL_labeling.txt',
                            'database': '/home/hui/Documents/VM_share/DIA_data/yeast_UPS_2017-09-14.fasta',
                            'pro_path': '/home/hui/Documents/VM_share/DIA_data/CPTAC/fragpipe/' + mode + '/'}

                flag = 1
                candidate_pros = get_candidate_proteins(in_files['pro_path'], acq, in_type, dataset, flag)
                in_files['all_candidate'] = candidate_pros

                if in_type == 'fragpipe':
                    wfs = getwfs('/data/res_DE_files/DDA_frag_rank_mean', k, dataset, mode)
                elif in_type == 'maxquant':
                    wfs = getwfs('/data/res_DE_files/DDA_mq_rank_mean', k, dataset, mode)

                for method in ['set', 'fisher', 'hurdle']:
                    ensemble_DEA_topKwf(dataset, in_type, method, mode, in_files, wfs, acq, k)

################################################################################
############################################################################
    dataset = 'yeast_DDA'  ## A549-K562
    mode = ''
    for in_type in ['fragpipe', 'maxquant']:
        acq = 'DDA'
        for k in [2, 3, 5, 7, 9, 11, 15, 20]:

            if acq == 'DIA':
                if in_type == 'maxquant':
                    in_files = {'protein': '/home/hui/Documents/VM_share/DIA_data/A549_K562/DIA/maxquant/txt/proteinGroups.txt',
                                'peptide': '/home/hui/Documents/VM_share/DIA_data/A549_K562/DIA/maxquant/txt/peptides.txt',
                                'design': '/home/hui/Documents/VM_share/DIA_data/A549_K562/design_cell_line.txt',
                                'msstats': '/home/hui/Documents/VM_share/DIA_data/A549_K562/DIA/maxquant/txt/evidence.txt',
                                'save': '/home/hui/PycharmProjects/DL_proteomics/benchmark_res/',
                                'design_stat': '/home/hui/Documents/VM_share/DIA_data/A549_K562/design_cell_line_msstats_DIA.txt',
                                'labels': '/home/hui/Documents/VM_share/DIA_data/A549_K562/CL_labeling.txt'}
                elif in_type == 'DIANN':
                    in_files = {'protein': '/home/hui/Documents/VM_share/DIA_data/A549_K562/DIA/DIANN_res/report.pg_matrix.tsv',
                                'protein1': '/home/hui/Documents/VM_share/DIA_data/A549_K562/DIA/DIANN_res/report.tsv',
                                'peptide': '/home/hui/Documents/VM_share/DIA_data/yeast/frager/combined_peptide.tsv',
                                'design1': '/home/hui/Documents/VM_share/DIA_data/A549_K562/design_cell_line_DIANN.txt',
                                'design': '/home/hui/Documents/VM_share/DIA_data/A549_K562/design_cell_line.txt',
                                'msstats': '/home/hui/Documents/VM_share/DIA_data/A549_K562/MSstats.csv',
                                'save': '/home/hui/PycharmProjects/DL_proteomics/benchmark_res/',
                                'design_stat': '/home/hui/Documents/VM_share/DIA_data/yeast/design_yeast_msstats.txt',
                                'labels': '/home/hui/Documents/VM_share/DIA_data/A549_K562/CL_labeling.txt'
                                }
            elif acq == 'DDA':
                if in_type == 'maxquant':
                    in_files = {'protein': '/home/hui/Documents/VM_share/DIA_data/yeast/maxquant/txt/proteinGroups.txt',
                                'peptide': '/home/hui/Documents/VM_share/DIA_data/yeast/maxquant/txt/peptides.txt',
                                'design': '/home/hui/Documents/VM_share/DIA_data/yeast/design_yeast.txt',
                                'msstats': '/home/hui/Documents/VM_share/DIA_data/yeast/maxquant/txt/evidence.txt',
                                'save': '/home/hui/PycharmProjects/DL_proteomics/benchmark_res/',
                                'design_stat': '/home/hui/Documents/VM_share/DIA_data/yeast/design_yeast_msstats.txt',
                                'labels': '/home/hui/Documents/VM_share/DIA_data/A549_K562/CL_labeling.txt',
                                'database': '/home/hui/Documents/VM_share/DIA_data/test.fasta',
                                'pro_path': '/home/hui/Documents/VM_share/DIA_data/yeast/maxquant/txt/'}
                elif in_type == 'fragpipe':
                    in_files = {'protein': '/home/hui/Documents/VM_share/DIA_data/yeast/frager/combined_protein.tsv',
                                'peptide': '/home/hui/Documents/VM_share/DIA_data/yeast/frager/combined_peptide.tsv',
                                'design': '/home/hui/Documents/VM_share/DIA_data/yeast/design_yeast.txt',
                                'msstats': '/home/hui/Documents/VM_share/DIA_data/yeast/frager/MSstats.csv',
                                'save': '/home/hui/PycharmProjects/DL_proteomics/benchmark_res/',
                                'design_stat': '/home/hui/Documents/VM_share/DIA_data/yeast/design_yeast_msstats.txt',
                                'labels': '/home/hui/Documents/VM_share/DIA_data/A549_K562/CL_labeling.txt',
                                'database': '/home/hui/Documents/VM_share/DIA_data/test.fasta',
                                'pro_path': '/home/hui/Documents/VM_share/DIA_data/yeast/frager/'}
            flag = 1
            candidate_pros = get_candidate_proteins(in_files['pro_path'], acq, in_type, dataset, flag)
            in_files['all_candidate'] = candidate_pros

            if in_type == 'fragpipe':
                wfs = getwfs('/data/res_DE_files/DDA_frag_rank_mean', k, dataset, mode)
            elif in_type == 'maxquant':
                wfs = getwfs('/data/res_DE_files/DDA_mq_rank_mean', k, dataset, mode)

            for method in ['set', 'fisher', 'hurdle']:
                # ensemble_DEA_threeinputs(dataset, in_type, method, mode, in_files)
                ensemble_DEA_topKwf(dataset, in_type, method, mode, in_files, wfs, acq, k)

#################################################################################
#############################################################################
    dataset = 'yeast1819_DDA'  ## A549-K562
    mode = ''
    for in_type in ['fragpipe', 'maxquant']:
        acq = 'DDA'
        for k in [2, 3, 5, 7, 9, 11, 15, 20]:

            if acq == 'DIA':
                if in_type == 'maxquant':
                    in_files = {
                        'protein': '/home/hui/Documents/VM_share/DIA_data/A549_K562/DIA/maxquant/txt/proteinGroups.txt',
                        'peptide': '/home/hui/Documents/VM_share/DIA_data/A549_K562/DIA/maxquant/txt/peptides.txt',
                        'design': '/home/hui/Documents/VM_share/DIA_data/A549_K562/design_cell_line.txt',
                        'msstats': '/home/hui/Documents/VM_share/DIA_data/A549_K562/DIA/maxquant/txt/evidence.txt',
                        'save': '/home/hui/PycharmProjects/DL_proteomics/benchmark_res/',
                        'design_stat': '/home/hui/Documents/VM_share/DIA_data/A549_K562/design_cell_line_msstats_DIA.txt',
                        'labels': '/home/hui/Documents/VM_share/DIA_data/A549_K562/CL_labeling.txt'}
                elif in_type == 'DIANN':
                    in_files = {
                        'protein': '/home/hui/Documents/VM_share/DIA_data/A549_K562/DIA/DIANN_res/report.pg_matrix.tsv',
                        'protein1': '/home/hui/Documents/VM_share/DIA_data/A549_K562/DIA/DIANN_res/report.tsv',
                        'peptide': '/home/hui/Documents/VM_share/DIA_data/yeast/frager/combined_peptide.tsv',
                        'design1': '/home/hui/Documents/VM_share/DIA_data/A549_K562/design_cell_line_DIANN.txt',
                        'design': '/home/hui/Documents/VM_share/DIA_data/A549_K562/design_cell_line.txt',
                        'msstats': '/home/hui/Documents/VM_share/DIA_data/A549_K562/MSstats.csv',
                        'save': '/home/hui/PycharmProjects/DL_proteomics/benchmark_res/',
                        'design_stat': '/home/hui/Documents/VM_share/DIA_data/yeast/design_yeast_msstats.txt',
                        'labels': '/home/hui/Documents/VM_share/DIA_data/A549_K562/CL_labeling.txt'
                        }
            elif acq == 'DDA':
                if in_type == 'maxquant':
                    in_files = {
                        'protein': '/home/hui/Documents/VM_share/DIA_data/yeast_1819/maxquant/combined/txt/proteinGroups.txt',
                        'peptide': '/home/hui/Documents/VM_share/DIA_data/yeast_1819/maxquant/combined/txt/peptides.txt',
                        'design': '/home/hui/Documents/VM_share/DIA_data/yeast_1819/design_yeast1819_mq.txt',
                        'msstats': '/home/hui/Documents/VM_share/DIA_data/yeast_1819/maxquant/combined/txt/evidence.txt',
                        'save': '/home/hui/PycharmProjects/DL_proteomics/benchmark_res/',
                        'design_stat': '/home/hui/Documents/VM_share/DIA_data/yeast_1819/design_yeast1819_msstats.txt',
                        'labels': '/home/hui/Documents/VM_share/DIA_data/A549_K562/CL_labeling.txt',
                        'database': '/home/hui/Documents/VM_share/DIA_data/test.fasta',
                        'pro_path': '/home/hui/Documents/VM_share/DIA_data/yeast_1819/maxquant/combined/txt/'}
                elif in_type == 'fragpipe':
                    in_files = {'protein': '/home/hui/Documents/VM_share/DIA_data/yeast_1819/fragpipe/combined_protein.tsv',
                                'peptide': '/home/hui/Documents/VM_share/DIA_data/yeast_1819/fragpipe/combined_peptide.tsv',
                                'design': '/home/hui/Documents/VM_share/DIA_data/yeast_1819/design_yeast1819_frag.txt',
                                'msstats': '/home/hui/Documents/VM_share/DIA_data/yeast_1819/fragpipe/MSstats.csv',
                                'save': '/home/hui/PycharmProjects/DL_proteomics/benchmark_res/',
                                'design_stat': '/home/hui/Documents/VM_share/DIA_data/yeast_1819/design_yeast1819_msstats.txt',
                                'labels': '/home/hui/Documents/VM_share/DIA_data/A549_K562/CL_labeling.txt',
                                'database': '/home/hui/Documents/VM_share/DIA_data/test.fasta',
                                'pro_path': '/home/hui/Documents/VM_share/DIA_data/yeast_1819/fragpipe/'}
            flag = 1
            candidate_pros = get_candidate_proteins(in_files['pro_path'], acq, in_type, dataset, flag)
            in_files['all_candidate'] = candidate_pros

            if in_type == 'fragpipe':
                wfs = getwfs('/data/res_DE_files/DDA_frag_rank_mean', k, dataset, mode)
            elif in_type == 'maxquant':
                wfs = getwfs('/data/res_DE_files/DDA_mq_rank_mean', k, dataset, mode)

            for method in ['set', 'fisher', 'hurdle']:
                # ensemble_DEA_threeinputs(dataset, in_type, method, mode, in_files)
                ensemble_DEA_topKwf(dataset, in_type, method, mode, in_files, wfs, acq, k)

#################################################################################
#############################################################################

    dataset = 'human_ecoli_DDA'  ## A549-K562
    comb_idxs = [i for i in range(1, 35)]

    for in_type in ['fragpipe', 'maxquant']:
        acq = 'DDA'
        for comb_idx in comb_idxs:
            for k in [2, 3, 5, 7, 9, 11, 15, 20]:
                if acq == 'DIA':
                    if in_type == 'maxquant':
                        in_files = {
                            'protein': '/home/hui/Documents/VM_share/DIA_data/A549_K562/DIA/maxquant/txt/proteinGroups.txt',
                            'peptide': '/home/hui/Documents/VM_share/DIA_data/A549_K562/DIA/maxquant/txt/peptides.txt',
                            'design': '/home/hui/Documents/VM_share/DIA_data/A549_K562/design_cell_line.txt',
                            'msstats': '/home/hui/Documents/VM_share/DIA_data/A549_K562/DIA/maxquant/txt/evidence.txt',
                            'save': '/home/hui/PycharmProjects/DL_proteomics/benchmark_res/',
                            'design_stat': '/home/hui/Documents/VM_share/DIA_data/A549_K562/design_cell_line_msstats_DIA.txt',
                            'labels': '/home/hui/Documents/VM_share/DIA_data/A549_K562/CL_labeling.txt'}
                    elif in_type == 'DIANN':
                        in_files = {
                            'protein': '/home/hui/Documents/VM_share/DIA_data/A549_K562/DIA/DIANN_res/report.pg_matrix.tsv',
                            'protein1': '/home/hui/Documents/VM_share/DIA_data/A549_K562/DIA/DIANN_res/report.tsv',
                            'peptide': '/home/hui/Documents/VM_share/DIA_data/yeast/frager/combined_peptide.tsv',
                            'design1': '/home/hui/Documents/VM_share/DIA_data/A549_K562/design_cell_line_DIANN.txt',
                            'design': '/home/hui/Documents/VM_share/DIA_data/A549_K562/design_cell_line.txt',
                            'msstats': '/home/hui/Documents/VM_share/DIA_data/A549_K562/MSstats.csv',
                            'save': '/home/hui/PycharmProjects/DL_proteomics/benchmark_res/',
                            'design_stat': '/home/hui/Documents/VM_share/DIA_data/yeast/design_yeast_msstats.txt',
                            'labels': '/home/hui/Documents/VM_share/DIA_data/A549_K562/CL_labeling.txt'
                            }
                elif acq == 'DDA':
                    if in_type == 'maxquant':
                        in_files = {
                            'protein': '/home/hui/Documents/VM_share/DIA_data/human_ecoli/combined/txt/proteinGroups' + str(
                                comb_idx) + '.txt',
                            'peptide': '/home/hui/Documents/VM_share/DIA_data/human_ecoli/combined/txt/peptides' + str(
                                comb_idx) + '.txt',
                            'design': '/home/hui/Documents/VM_share/DIA_data/human_ecoli/design_he_mq' + str(comb_idx) + '.txt',
                            'msstats': '/home/hui/Documents/VM_share/DIA_data/human_ecoli/combined/txt/evidence' + str(
                                comb_idx) + '.txt',
                            'save': '/home/hui/PycharmProjects/DL_proteomics/benchmark_res/',
                            'design_stat': '/home/hui/Documents/VM_share/DIA_data/human_ecoli/design_he' + str(
                                comb_idx) + '_msstats.txt',
                            'labels': '/home/hui/Documents/VM_share/DIA_data/A549_K562/CL_labeling.txt',
                            'database': '/home/hui/Documents/VM_share/DIA_data/ecoli_human.fasta',
                            'pro_path': '/home/hui/Documents/VM_share/DIA_data/human_ecoli/combined/txt/'}
                    elif in_type == 'fragpipe':
                        in_files = {'protein': '/home/hui/Documents/VM_share/DIA_data/human_ecoli/frag/combined_protein' + str(
                            comb_idx) + '.tsv',
                                    'peptide': '/home/hui/Documents/VM_share/DIA_data/human_ecoli/frag/combined_peptide' + str(
                                        comb_idx) + '.tsv',
                                    'design': '/home/hui/Documents/VM_share/DIA_data/human_ecoli/design_he_frag' + str(
                                        comb_idx) + '.txt',
                                    'msstats': '/home/hui/Documents/VM_share/DIA_data/human_ecoli/frag/MSstats' + str(
                                        comb_idx) + '.csv',
                                    'save': '/home/hui/PycharmProjects/DL_proteomics/benchmark_res/',
                                    'design_stat': '/home/hui/Documents/VM_share/DIA_data/human_ecoli/design_he' + str(
                                        comb_idx) + '_msstats.txt',
                                    'labels': '/home/hui/Documents/VM_share/DIA_data/A549_K562/CL_labeling.txt',
                                    'database': '/home/hui/Documents/VM_share/DIA_data/ecoli_human.fasta',
                                    'pro_path': '/home/hui/Documents/VM_share/DIA_data/human_ecoli/frag/'}
                flag = 1
                candidate_pros = get_candidate_proteins(in_files['pro_path'], acq, in_type, dataset, flag)
                in_files['all_candidate'] = candidate_pros

                if in_type == 'fragpipe':
                    wfs = getwfs('/data/res_DE_files/DDA_frag_rank_mean', k, dataset, str(comb_idx))
                elif in_type == 'maxquant':
                    wfs = getwfs('/data/res_DE_files/DDA_mq_rank_mean', k, dataset, str(comb_idx))

                for method in ['set', 'fisher', 'hurdle']:
                    ensemble_DEA_topKwf(dataset, in_type, method, str(comb_idx), in_files, wfs, acq, k)