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
save_folder = 'benchmark_res/'

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

def get_all_proteins(fasta_name):
    filepath = fasta_name
    fasta = {}
    with open(filepath) as file_one:
        for line in file_one:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                active_sequence_name = line[1:]
                if active_sequence_name not in fasta:
                    fasta[active_sequence_name] = ''
                continue
            sequence = line
            fasta[active_sequence_name] += sequence

    protein_names = []
    protein_seqs = []
    protein_len = []
    for key in fasta.keys():
        strs = key.split(' ')
        protein_names.append(strs[0])
        protein_len.append(len(fasta[key]))
        protein_seqs.append(fasta[key])

    out = np.column_stack((np.array([protein_names]).T, np.array([protein_len]).T))
    return fasta, out

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

def get_true_logFC(label, contrast):
    true_logFC=np.zeros((len(label), 1))
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

def get_comparable_res(proteins, score, label, const, logFC,real_logFC, in_files):

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

    real_logFC_c = get_true_logFC(label_c, const_c)
    return np.array(proteins_c), np.array(score_c), np.array(label_c), np.array(const_c), np.array(logFC_c), np.array(real_logFC_c)

def cal_mse_auc(tool, acq, in_type, inten_ty, imput, normal, res_file, dataset, run_time, in_files):
    if(os.path.exists(res_file)):
        res = pd.read_csv(res_file, sep=',', header=0)
        idx = np.where(((~np.isnan(np.array(res['logFC']))) & (~np.isnan(np.array(res['pvalue']))) & (~np.isnan(np.array(res['adj.pvalue']))) &
                        (~np.isinf(np.array(res['logFC']))) & (~np.isinf(np.array(res['pvalue']))) & (~np.isinf(np.array(res['adj.pvalue'])))))[0]

        proteins = np.array(res['protein'])[idx]
        score = 1 - np.array(res['adj.pvalue'])[idx]

        label = get_label_ups(proteins, '_HUMAN')
        const = np.array(res['contrast'])[idx]
        logFC = np.array(res['logFC'])[idx]
        real_logFC = get_true_logFC(label, const)

        protein_c, score_c, label_c, const_c, logFC_c, real_logFC_c = get_comparable_res(proteins, score, label, const, logFC,real_logFC, in_files)

        p = 0.01
        aucs, paucs1, paucs2, paucs3, sprs, rmses, lens, Cls_met_1_001, Cls_met_1_005 = cal_auc(score, label, const, p, logFC, real_logFC)
        out1 = [tool, in_type, inten_ty, imput, normal, run_time]
        out2 = [tool, in_type, inten_ty, imput, normal]
        out3 = [tool, in_type, inten_ty, imput, normal, run_time]
        out4 = [tool, in_type, inten_ty, imput, normal]

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
        out5 = [tool, in_type, inten_ty, imput, normal, run_time]
        out6 = [tool, in_type, inten_ty, imput, normal]
        out7 = [tool, in_type, inten_ty, imput, normal, run_time]
        out8 = [tool, in_type, inten_ty, imput, normal]
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
        out1 = [tool, in_type, inten_ty, imput, normal, run_time]
        out2 = [tool, in_type, inten_ty, imput, normal]
        out3 = [tool, in_type, inten_ty, imput, normal, run_time]
        out4 = [tool, in_type, inten_ty, imput, normal]
        out1.extend(['bug'])
        out2.extend(['bug'])
        out3.extend(['bug'])
        out4.extend(['bug'])

        out5 = [tool, in_type, inten_ty, imput, normal, run_time]
        out6 = [tool, in_type, inten_ty, imput, normal]
        out7 = [tool, in_type, inten_ty, imput, normal, run_time]
        out8 = [tool, in_type, inten_ty, imput, normal]

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

def run_Rscript(acq, tool, in_type, inten_ty, imput, normal, in_files, dataset):
    '''
    :param tool: DE tool
    :param in_type: maxquant/fragpipe
    :param inten_ty: ''/MaxLFQ for pipe or ''/'LFQ' for maxquant
    :param imput: imputation method
    :param normal: normalization method
    :param in_files: full file names for input files and save fold
    :return: DE results
    '''
    save_fold = in_files['save'] + dataset + '_' + acq + '/' + tool + '/'
    t1 = time.time()

    if in_type == 'DIANN':
        if inten_ty == 'norm':
            print('Rscript ' + R_fold + tool + '_DE.R' + ' ' + in_files['protein'] + ' ' + in_files['peptide']
                      + ' ' + in_files['design'] + ' ' + in_files['msstats'] + ' ' + in_type + ' ' + inten_ty + ' '
                      + str(imput) + ' ' + str(normal) + ' ' + save_fold)
            os.system('Rscript ' + R_fold + tool + '_DE.R' + ' ' + in_files['protein'] + ' ' + in_files['peptide']
                      + ' ' + in_files['design1'] + ' ' + in_files['msstats'] + ' ' + in_type + ' ' + inten_ty + ' '
                      + str(imput) + ' ' + str(normal) + ' ' + save_fold)
        else:
            print('Rscript ' + R_fold + tool + '_DE.R' + ' ' + in_files['protein1'] + ' ' + in_files['peptide']
                  + ' ' + in_files['design1'] + ' ' + in_files['msstats'] + ' ' + in_type + ' ' + inten_ty + ' '
                  + str(imput) + ' ' + str(normal) + ' ' + save_fold)
            os.system('Rscript ' + R_fold + tool + '_DE.R' + ' ' + in_files['protein1'] + ' ' + in_files['peptide']
                      + ' ' + in_files['design1'] + ' ' + in_files['msstats'] + ' ' + in_type + ' ' + inten_ty + ' '
                      + str(imput) + ' ' + str(normal) + ' ' + save_fold)


    else:
        print('Rscript ' + R_fold + tool + '_DE.R' + ' ' + in_files['protein'] + ' ' + in_files['peptide']
                  + ' ' + in_files['design'] + ' ' + in_files['msstats'] + ' ' + in_type + ' ' + inten_ty + ' '
                  + str(imput) + ' ' + str(normal) + ' ' + save_fold)
        os.system('Rscript ' + R_fold + tool + '_DE.R' + ' ' + in_files['protein'] + ' ' + in_files['peptide']
                  + ' ' + in_files['design'] + ' ' + in_files['msstats'] + ' ' + in_type + ' ' + inten_ty + ' '
                  + str(imput) + ' ' + str(normal) + ' ' + save_fold)
    t2 = time.time()
    run_time = t2-t1

    if (inten_ty == 'basic') | (inten_ty == 'blank'):
        inten_ty = ''

    if imput == 'blank':
        imput = ''

    if normal == 'blank':
        normal = ''
    res_file_name = save_fold + tool + '_' + in_type + '_' + inten_ty + '_' + imput + '_' + normal + '.csv'
    cal_mse_auc(tool, acq, in_type, inten_ty, imput, normal, res_file_name, dataset, run_time, in_files)


def run_all_group1(acq, tools, in_type, imputs, normals, in_files, dataset):

    for tool in tools:
        if in_type == 'maxquant':
            inten_types = ['basic', 'LFQ']
        elif in_type == 'fragpipe':
            inten_types = ['basic', 'MaxLFQ']
        elif in_type == 'DIANN':
            inten_types = ['raw', 'norm_raw', 'MaxLFQ', 'MaxLFQ_raw', 'norm']

        for inten_ty in inten_types:
            for imput in imputs:
                for normal in normals:
                    run_Rscript(acq, tool, in_type, inten_ty, imput, normal, in_files, dataset)

def run_all_group2(acq, tools, in_type, imputs, normals, in_files, dataset):
    for tool in tools:
        for imput in imputs:
            for normal in normals:
                run_Rscript(acq, tool, in_type, 'basic', imput, normal, in_files, dataset)

def run_Mstats(acq, in_files, in_type, dataset):
    in_files.update({'design':in_files['design_stat']})
    for normal in ['equalizeMedians', 'quantile', 'globalStandards', 'FALSE']:
        run_Rscript(acq, 'MSstats', in_type, 'basic', 'blank', normal, in_files, dataset)

def run_ProPCA(acq, in_files, in_type, dataset):
    tool = 'ProPCA'
    imputs = ['blank', 'bpca', 'knn', 'QRILC', 'MLE', 'MinDet', 'MinProb', 'min', 'zero', 'nbavg', 'mice', 'missForest']
    normals = ['blank', "sum", "max", "center.mean", "center.median", "div.mean", "div.median", "quantiles",
               "quantiles.robust", "vsn"]
    for inten_ty in ['basic', 'LFQ']:
        for imput1 in imputs:
            for imput2 in imputs:
                for normal1 in normals:
                    for normal2 in normals:
                        save_fold = in_files['save'] + dataset + '_' + acq + '/' + tool + '/'

                        t1 = time.time()

                        os.system(
                            'Rscript ' + R_fold + 'ProPCA' + '_DE.R' + ' ' + in_files['protein'] + ' ' + in_files['peptide']
                            + ' ' + in_files['design'] + ' ' + in_files['msstats'] + ' ' + in_type + ' ' + inten_ty + ' '
                            + imput1 + ' ' + imput2 + ' ' + normal1 + ' ' + normal2 + ' ' + save_fold)
                        t2 = time.time()
                        run_time = t2 - t1

                        if inten_ty == 'basic':
                            inten_ty0 = ''
                        else:
                            inten_ty0 = inten_ty

                        if imput1 == 'blank':
                            imput10 = ''
                        else:
                            imput10 = imput1

                        if imput2 == 'blank':
                            imput20 = ''
                        else:
                            imput20 = imput2

                        if normal1 == 'blank':
                            normal10 = ''
                        else:
                            normal10 = normal1

                        if normal2 == 'blank':
                            normal20 = ''
                        else:
                            normal20 = normal2

                        res_file_name = save_fold + tool + '_' + in_type + '_' + inten_ty0 + '_' + imput10 + '_' + imput20 + '_' + normal10 + '_' + normal20 + '.csv'
                        cal_mse_auc(tool, acq, in_type, inten_ty0, imput10 + '|' + imput20, normal10 + '|' + normal20, res_file_name, dataset, run_time, in_files)


def run_all(acq, in_files, dataset, in_type, flag):

    tool_group1 = ['limma', 'siggenes', 'ttest']
    tool_group4 = ['DEP', 'DEqMS','ANOVA']

    tool_group2 = ['edgeR', 'plgem', 'beta_binomial']

    tool_group3 = ['msqrob2']
    tool_group5 = ['ProteoMM']
    tool_group6 = ['SAM', 'ROTS']
    tool_group7 = ['proDA']

    imputs1 = ['blank', 'bpca', 'knn', 'QRILC', 'MLE', 'MinDet', 'MinProb', 'min', 'zero', 'nbavg', 'mice', 'missForest']

    normals = ['blank', "sum", "max", "center.mean", "center.median", "div.mean", "div.median", "quantiles",
               "quantiles.robust", "vsn"]

    if flag == 1:
        run_all_group1(acq, tool_group1, in_type, imputs1, normals, in_files, dataset)
    elif flag == 2:
        run_all_group2(acq, tool_group2, in_type, imputs1, normals, in_files, dataset)
    elif flag == 3:
        run_all_group1(acq, tool_group3, in_type, imputs1, normals, in_files, dataset)

    elif flag == 4:
        run_Mstats(acq, in_files, in_type, dataset)
    elif flag == 5:
        run_all_group1(acq, tool_group7, in_type, imputs1, normals, in_files, dataset)
    elif flag == 6:
        run_all_group1(acq, tool_group4, in_type, imputs1, normals, in_files, dataset)
    elif flag == 7:
        run_all_group1(acq, tool_group5, in_type, imputs1, normals, in_files, dataset)
    elif flag == 8:
        run_all_group1(acq, tool_group6, in_type, imputs1, normals, in_files, dataset)

if __name__ == '__main__':
    '''
    quantifcation tools:
    1. limma (intensity)
    2. edgeR (count)
    3. siggenes (intensity)
    4. proDA (intensity)
    5. DEP (intensity)
    6. DEqMS (intensity+count)
    7. plgem (count)
    8. msqrob2 (peptide intenisty)
    9. ProteoMM (peptide intensity)
    10. Msstats (intensity)
    11. ProPCA (count + peptide intensity) (remove)
    12. ANOVA (intensity)
    13. t-test (intensity)
    14. beta bionomial (count)
    15. perseus (intensity) (remove)
    16. SAM (intensity) (add)
    17. ROTS (intensity) (add)
    '''

    flag = int(sys.argv[1])
    acq = sys.argv[2]
    in_type = sys.argv[3]

    dataset = 'yeast1819'

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
            in_files = {'protein': '/home/hui/Documents/VM_share/DIA_data/yeast_1819/maxquant/combined/txt/proteinGroups.txt',
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
                        'database':'/home/hui/Documents/VM_share/DIA_data/test.fasta',
                        'pro_path':'/home/hui/Documents/VM_share/DIA_data/yeast_1819/fragpipe/'}

    candidate_pros = get_candidate_proteins(in_files['pro_path'], acq, in_type, dataset, flag)
    in_files['all_candidate'] = candidate_pros

    run_all(acq, in_files, dataset, in_type, flag)
