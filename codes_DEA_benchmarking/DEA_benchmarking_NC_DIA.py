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
from sklearn.metrics import accuracy_score as acc, recall_score as recall, precision_score as precision, f1_score as f1, \
    matthews_corrcoef as mcc
from sklearn.metrics import confusion_matrix


def seed_torch(seed=1029):
    random.seed(seed)
    os.environ['PYTHONHASHSEED'] = str(seed)
    np.random.seed(seed)


seed_torch()


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
    #tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()
    if len(np.unique(label)) == 1:
        tn, fp, fn, tp = 0, 0, 0, 0
        spec = 0
    else:
        tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()
        spec = tn / (tn + fp)

    geomean = np.sqrt(Rec * spec)
    return Acc, Prec, Rec, F1, F1w, Mcc, nMcc, geomean, tn, fp, fn, tp


def cal_auc(score, score_lfq, label, const, p, logFC, rea_logFC):
    label = np.array(label)
    score = np.array(score)
    score_lfq = np.array(score_lfq)
    score_lfq[np.where((np.isnan(score_lfq)))[0]] = 0
    score[np.where((np.isnan(score)))[0]] = 0
    const = np.array(const)

    aucs = []
    paucs1 = []
    paucs2 = []
    paucs3 = []
    paucs1_lfq = []
    paucs2_lfq = []
    paucs3_lfq = []

    sprs = []
    rmses = []
    lens = []

    Cls_met_1_001 = []
    Cls_met_1_005 = []

    cls_met_1_001 = cls_metric(1 - score, label, logFC, 0.01, np.log2(1.5))
    cls_met_1_005 = cls_metric(1 - score, label, logFC, 0.05, np.log2(1.5))

    Cls_met_1_001.extend(cls_met_1_001)
    Cls_met_1_005.extend(cls_met_1_005)

    if len(np.unique(label)) == 1:
        auc = 0
        pauc1 = 0
        pauc2 = 0
        pauc3 = 0

        pauc1_lfq = 0
        pauc2_lfq = 0
        pauc3_lfq = 0
    else:
        auc = roc_auc_score(label, score)
        pauc1 = roc_auc_score(label, score, max_fpr=0.01)
        pauc2 = roc_auc_score(label, score, max_fpr=0.05)
        pauc3 = roc_auc_score(label, score, max_fpr=0.1)

        pauc1_lfq = roc_auc_score(label, score_lfq, max_fpr=0.01)
        pauc2_lfq = roc_auc_score(label, score_lfq, max_fpr=0.05)
        pauc3_lfq = roc_auc_score(label, score_lfq, max_fpr=0.1)

    spr = stats.spearmanr(logFC, rea_logFC)[0]
    rmse = np.sqrt(mean_squared_error(logFC, rea_logFC))

    aucs.append(auc)
    paucs1.append(pauc1)
    paucs2.append(pauc2)
    paucs3.append(pauc3)

    paucs1_lfq.append(pauc1_lfq)
    paucs2_lfq.append(pauc2_lfq)
    paucs3_lfq.append(pauc3_lfq)

    sprs.append(spr)
    rmses.append(rmse)
    lens.append(len(label))

    uni_cons = np.unique(const)
    for con in uni_cons:
        idx = np.where((const == con))[0]

        cls_met_1_001 = cls_metric(1 - score[idx], label[idx], logFC[idx], 0.01, np.log2(1.5))
        cls_met_1_005 = cls_metric(1 - score[idx], label[idx], logFC[idx], 0.05, np.log2(1.5))

        Cls_met_1_001.extend(cls_met_1_001)
        Cls_met_1_005.extend(cls_met_1_005)

        if len(np.unique(label[idx])) == 1:
            auc = 0
            pauc1 = 0
            pauc2 = 0
            pauc3 = 0

            pauc1_lfq = 0
            pauc2_lfq = 0
            pauc3_lfq = 0

        else:
            auc = roc_auc_score(label[idx], score[idx])
            pauc1 = roc_auc_score(label[idx], score[idx], max_fpr=0.01)
            pauc2 = roc_auc_score(label[idx], score[idx], max_fpr=0.05)
            pauc3 = roc_auc_score(label[idx], score[idx], max_fpr=0.1)

            pauc1_lfq = roc_auc_score(label[idx], score_lfq[idx], max_fpr=0.01)
            pauc2_lfq = roc_auc_score(label[idx], score_lfq[idx], max_fpr=0.05)
            pauc3_lfq = roc_auc_score(label[idx], score_lfq[idx], max_fpr=0.1)

        aucs.append(auc)
        paucs1.append(pauc1)
        paucs2.append(pauc2)
        paucs3.append(pauc3)

        paucs1_lfq.append(pauc1_lfq)
        paucs2_lfq.append(pauc2_lfq)
        paucs3_lfq.append(pauc3_lfq)

        sprs.append(stats.spearmanr(logFC[idx], rea_logFC[idx])[0])
        rmses.append(np.sqrt(mean_squared_error(logFC[idx], rea_logFC[idx])))
        lens.append(len(idx))

    return aucs, paucs1, paucs2, paucs3, paucs1_lfq, paucs2_lfq, paucs3_lfq, sprs, rmses, lens, Cls_met_1_001, Cls_met_1_005


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

    cls_met_1_001 = cls_metric(1 - score, label, logFC, 0.01, np.log2(1.5))
    cls_met_1_005 = cls_metric(1 - score, label, logFC, 0.05, np.log2(1.5))

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

        cls_met_1_001 = cls_metric(1 - score[idx], label[idx], logFC[idx], 0.01, np.log2(1.5))
        cls_met_1_005 = cls_metric(1 - score[idx], label[idx], logFC[idx], 0.05, np.log2(1.5))

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
    true_logFC = np.zeros((len(label), 1))
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
    return true_logFC


def get_uniport_id(protein_list):
    pro_uniport = []
    for i in protein_list:
        pros = i.split(';')
        pro = ''
        for j in pros:

            proj = j.split('|')
            if proj[0] == 'sp':
                pro = pro + proj[1]
            else:
                pro = pro + proj[0]
        pro_uniport.append(pro)
    return pro_uniport


def get_comparable_res(proteins, score, score_lfq, label, const, logFC, real_logFC, organisms, in_files):
    proteins_all = pd.read_csv(in_files['all_candidate'], header=0, sep='\t')  # ['Protein'].values
    pro_uniport_all = get_uniport_id(proteins_all['Protein'].values)
    uni_cons = np.unique(const)
    proteins_c = []
    score_c = []
    score_lfq_c = []
    label_c = []
    const_c = []
    logFC_c = []
    real_logFC_c = []
    for con in uni_cons:
        idx_con = np.where((const == con))[0]
        pro_con = proteins[idx_con]
        pro_con_uniport = get_uniport_id(pro_con)
        non_pro = np.setdiff1d(pro_uniport_all, pro_con_uniport)

        comm, idx1, idx2 = np.intersect1d(non_pro, pro_uniport_all, return_indices=True)

        proteins_c.extend(proteins[idx_con])
        proteins_c.extend(proteins_all['Protein'].values[idx2])
        score_c.extend(score[idx_con])
        score_c.extend([0 for i in range(len(non_pro))])
        score_lfq_c.extend(score_lfq[idx_con])
        score_lfq_c.extend([0 for i in range(len(non_pro))])
        label_c.extend(label[idx_con])
        label_c.extend(proteins_all['DEP'].values[idx2])
        const_c.extend(const[idx_con])
        const_c.extend([con for i in range(len(non_pro))])
        logFC_c.extend(logFC[idx_con])
        logFC_c.extend([0 for i in range(len(non_pro))])
        real_logFC_c.extend(real_logFC[idx_con])
        ext_rlfc = []
        for i in range(len(idx2)):
            organ = proteins_all['Organism'].values[idx2[i]]
            if(len(np.where(((const == con) & (organisms == organ)))[0]) > 0):
                ext_rlfc.append(real_logFC[np.where(((const == con) & (organisms == organ)))[0][0]])
            else:
                ext_rlfc.append(0)
        real_logFC_c.extend(ext_rlfc)

    # real_logFC_c = get_true_logFC(label_c, const_c)
    return np.array(proteins_c), np.array(score_c), np.array(score_lfq_c), np.array(label_c), np.array(
        const_c), np.array(logFC_c), np.array(real_logFC_c)


def cal_mse_auc(tool, acq, in_type, inten_ty, imput, normal, res_file, dataset, run_time, in_files):
    if (os.path.exists(res_file)):
        res = pd.read_csv(res_file, sep=',', header=0)
        idx = np.where(((~np.isnan(np.array(res['logFC']))) & (~np.isnan(np.array(res['pvalue']))) & (
            ~np.isnan(np.array(res['adj.pvalue']))) &
                        (~np.isinf(np.array(res['logFC']))) & (~np.isinf(np.array(res['pvalue']))) & (
                            ~np.isinf(np.array(res['adj.pvalue'])))))[0]

        proteins = np.array(res['protein'])[idx]
        qval = np.array(res['adj.pvalue'])[idx]
        p_value = np.array(res['pvalue'])[idx]

        label = np.array(res['label'])[idx]
        const = np.array(res['contrast'])[idx]
        logFC = np.array(res['logFC'])[idx]

        lgq = np.array([-1 if qval[i] == 0 else -10 * np.log10(qval[i]) for i in range(0, len(qval))], dtype='float')
        max_q = 160
        if max(lgq) == -1:
            lgp = np.array([-1 if p_value[i] == 0 else -10 * np.log10(p_value[i]) for i in range(0, len(p_value))],
                           dtype='float')
            if max(lgp) == -1:
                max_q = 160
            else:
                max_q = max(lgp[np.where((lgp != -1))[0]])
        else:
            max_q = max(lgq[np.where((lgq != -1))[0]])

        lgq[np.where((lgq == -1))[0]] = max_q

        score = 1 - np.array(res['adj.pvalue'])[idx]
        score_lgq = np.array([abs(logFC[i] * lgq[i]) for i in range(len(lgq))], dtype='float')

        real_logFC = np.array(res['TlogFC'])[idx]
        organisms = np.array(res['organism'])[idx]

        protein_c, score_c, score_lgq_c, label_c, const_c, logFC_c, real_logFC_c = get_comparable_res(proteins, score,
                                                                                                      score_lgq, label,
                                                                                                      const,
                                                                                                      logFC, real_logFC,
                                                                                                      organisms,
                                                                                                      in_files)

        p = 0.01
        aucs, paucs1, paucs2, paucs3, paucs1_lfq, paucs2_lfq, paucs3_lfq, sprs, rmses, lens, Cls_met_1_001, Cls_met_1_005 = cal_auc(
            score, score_lgq, label, const, p,
            logFC, real_logFC)
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
        out2.extend(paucs1_lfq)
        out2.extend(paucs2_lfq)
        out2.extend(paucs3_lfq)
        out3.extend(Cls_met_1_001)
        out4.extend(Cls_met_1_005)

        aucs, paucs1, paucs2, paucs3, paucs1_lfq, paucs2_lfq, paucs3_lfq, sprs, rmses, lens, Cls_met_1_001, Cls_met_1_005 = cal_auc(
            score_c, score_lgq_c, label_c,
            const_c, p, logFC_c,
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
        out6.extend(paucs1_lfq)
        out6.extend(paucs2_lfq)
        out6.extend(paucs3_lfq)
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

    save_path = in_files['save_folder'] + dataset + '_' + acq + '_' + in_type + '_aucs1.csv'
    with open(save_path, 'a+', newline='') as ws:
        writer = csv.writer(ws)
        writer.writerow(out1)

    save_path = in_files['save_folder'] + dataset + '_' + acq + '_' + in_type + '_paucs1.csv'
    with open(save_path, 'a+', newline='') as ws:
        writer = csv.writer(ws)
        writer.writerow(out2)

    save_path = in_files['save_folder'] + dataset + '_' + acq + '_' + in_type + '_cls_001_1.csv'
    with open(save_path, 'a+', newline='') as ws:
        writer = csv.writer(ws)
        writer.writerow(out3)

    save_path = in_files['save_folder'] + dataset + '_' + acq + '_' + in_type + '_cls_005_1.csv'
    with open(save_path, 'a+', newline='') as ws:
        writer = csv.writer(ws)
        writer.writerow(out4)

    save_path = in_files['save_folder'] + dataset + '_' + acq + '_' + in_type + '_aucs2.csv'
    with open(save_path, 'a+', newline='') as ws:
        writer = csv.writer(ws)
        writer.writerow(out5)

    save_path = in_files['save_folder'] + dataset + '_' + acq + '_' + in_type + '_paucs2.csv'
    with open(save_path, 'a+', newline='') as ws:
        writer = csv.writer(ws)
        writer.writerow(out6)

    save_path = in_files['save_folder'] + dataset + '_' + acq + '_' + in_type + '_cls_001_2.csv'
    with open(save_path, 'a+', newline='') as ws:
        writer = csv.writer(ws)
        writer.writerow(out7)

    save_path = in_files['save_folder'] + dataset + '_' + acq + '_' + in_type + '_cls_005_2.csv'
    with open(save_path, 'a+', newline='') as ws:
        writer = csv.writer(ws)
        writer.writerow(out8)


# def get_candidate_proteins(in_files, acq, platform, dataset):
#     save_path = in_files['save_folder'] + dataset + '_' + acq + '_' + platform + '_all_proteins.csv'
#     os.system('Rscript ' + R_fold + 'get_candidate_proteins_NC.R ' + platform + ' ' + path + ' ' + save_path)
#     candidate_pros = pd.read_csv(save_path, header=None, sep=',').values
#     return candidate_pros


def get_label_ups(proteins, key, direct):
    if direct == 'positive':
        label = np.array([key in proteins[i] for i in range(len(proteins))])
    else:
        label = np.array([key not in proteins[i] for i in range(len(proteins))])

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


def run_Rscript(acq, tool, in_type, inten_ty, imput, normal, in_files, dataset, true_organism, DE_organism):
    '''
    :param tool: DE tool
    :param in_type: platform maxquant/fragpipe
    :param inten_ty: ''/MaxLFQ for pipe or ''/'LFQ' for maxquant
    :param imput: imputation method
    :param normal: normalization method
    :param in_files: full file names for input files and save fold
    :return: DE results
    '''
    save_fold = in_files['save_folder'] + dataset + '/'

    if not os.path.exists(save_fold):
        os.makedirs(save_fold)


    #print_lab = 'T'
    print_lab = 'T'
    if inten_ty == 'dlfq':
        logT = 'T'
    else:
        logT = 'F'

    print(in_files['Rscript_path'] + 'Rscript ' + R_fold + tool + '_DE.R' + ' ' + in_files['maxtrix_folder'] + ' ' +
          in_files['mstats_file'] +
          ' ' + in_files['evidence_file'] + ' ' + in_files['msstats_design_file'] + ' ' +
          in_files['main_output_protein_file'] + ' ' + in_files['main_output_peptide_file'] + ' ' +
          in_type + ' ' + inten_ty + ' ' +
          str(imput) + ' ' + str(normal) + ' ' + dataset + ' ' + save_fold + ' ' + print_lab + ' ' + true_organism +
          ' ' + DE_organism + ' ' + in_files['true_lgfc'] + ' ' + logT + ' ' + in_files['contrasts'])

    if (inten_ty == 'basic') | (inten_ty == 'blank'):
        inten_ty0 = ''
    else:
        inten_ty0 = inten_ty

    if imput == 'blank':
        imput0 = ''
    else:
        imput0 = imput

    if normal == 'blank':
        normal0 = ''
    else:
        normal0 = normal

    res_file_name = save_fold + dataset + '_' + tool + '_' + in_type + '_' + inten_ty0 + '_' + imput0 + '_' + normal0 + '.csv'

    t1 = time.time()
    os.system(in_files['Rscript_path'] + 'Rscript ' + R_fold + tool + '_DE.R' + ' ' + in_files['maxtrix_folder'] + ' ' +
              in_files['mstats_file'] +
              ' ' + in_files['evidence_file'] + ' ' + in_files['msstats_design_file'] + ' ' +
              in_files['main_output_protein_file'] + ' ' + in_files['main_output_peptide_file'] + ' ' +
              in_type + ' ' + inten_ty + ' ' +
              str(imput) + ' ' + str(normal) + ' ' + dataset + ' ' + save_fold + ' ' + print_lab + ' ' + true_organism +
              ' ' + DE_organism + ' ' + in_files['true_lgfc'] + ' ' + logT + ' ' + in_files['contrasts'])

    t2 = time.time()
    run_time = t2 - t1
    #in_files['save_folder'] = in_files['save_folder'] + 'new3/'
    cal_mse_auc(tool, acq, in_type, inten_ty, imput, normal, res_file_name, dataset, run_time, in_files)

    if os.path.exists(res_file_name):
        run_time = 0
        if not os.path.exists(in_files['save_folder'] + 'new3/'):
            os.makedirs(in_files['save_folder'] + 'new3/')
        in_files['save_folder'] = in_files['save_folder'] + 'new3/'
        cal_mse_auc(tool, acq, in_type, inten_ty, imput, normal, res_file_name, dataset, run_time, in_files)
        in_files['save_folder'] = in_files['save_folder'].replace('new3/', '')
    else:
        t1 = time.time()
        os.system(in_files['Rscript_path'] + 'Rscript ' + R_fold + tool + '_DE.R' + ' ' + in_files['maxtrix_folder'] + ' ' +
                  in_files['mstats_file'] +
                  ' ' + in_files['evidence_file'] + ' ' + in_files['msstats_design_file'] + ' ' +
                  in_files['main_output_protein_file'] + ' ' + in_files['main_output_peptide_file'] + ' ' +
                  in_type + ' ' + inten_ty + ' ' +
                  str(imput) + ' ' + str(normal) + ' ' + dataset + ' ' + save_fold + ' ' + print_lab + ' ' + true_organism +
                  ' ' + DE_organism + ' ' + in_files['true_lgfc'] + ' ' + logT + ' ' + in_files['contrasts'])

        t2 = time.time()
        run_time = t2 - t1

        in_files['save_folder'] = in_files['save_folder'] + 'new3/'
        cal_mse_auc(tool, acq, in_type, inten_ty, imput, normal, res_file_name, dataset, run_time, in_files)
        in_files['save_folder'] = in_files['save_folder'].replace('new3/', '')


def run_all_group1(acq, tools, platform, imputs, normals, in_files, dataset, true_organism,
                   DE_organism):  # protein intensity

    for tool in tools:
        inten_types = ['top1','top3', 'LFQ', 'dlfq']  # ['LFQ']#, 'dlfq']#

        for inten_ty in inten_types:
            for imput in imputs:
                for normal in normals:
                    run_Rscript(acq, tool, platform, inten_ty, imput, normal, in_files, dataset, true_organism,
                                DE_organism)


def run_all_group2(acq, tools, in_type, imputs, normals, in_files, dataset, true_organism,
                   DE_organism):  # protein count
    for tool in tools:
        for imput in imputs:
            for normal in normals:
                run_Rscript(acq, tool, in_type, 'count', imput, normal, in_files, dataset, true_organism, DE_organism)


def run_Mstats(acq, in_files, in_type, dataset, true_organism, DE_organism):
    # in_files.update({'design':in_files['design_stat']})
    for inten_ty in ['all', 'top3']:
        for normal in ['equalizeMedians', 'quantile', 'FALSE']: #, 'globalStandards'
            for imput in ['TRUE', 'FALSE']:
                run_Rscript(acq, 'MSstats', in_type, inten_ty, imput, normal, in_files, dataset, true_organism,
                            DE_organism)


def run_all_group3(acq, tools, platform, imputs, normals, in_files, dataset, true_organism,
                   DE_organism):  # pep intensity

    for tool in tools:
        inten_types = ['top0', 'LFQ']  # ['LFQ']#, 'dlfq']#

        for inten_ty in inten_types:
            for imput in imputs:
                for normal in normals:
                    run_Rscript(acq, tool, platform, inten_ty, imput, normal, in_files, dataset, true_organism,
                                DE_organism)


def run_all(acq, in_files, dataset, platform, flag, true_organism, DE_organism):

    tool_group1 = ['limma', 'ttest'] # 'siggenes', the same as SAM
    tool_group4 = ['DEP', 'ANOVA']# 'DEqMS',

    tool_group2 = ['edgeR', 'plgem', 'beta_binomial']

    tool_group3 = ['msqrob2']
    tool_group5 = ['ProteoMM']
    tool_group6 = ['SAM', 'ROTS']
    tool_group7 = ['proDA']

    imputs1 = ['MLE', 'missForest', 'blank', 'bpca', 'knn', 'QRILC', 'MinDet', 'MinProb', 'min', 'zero', 'nbavg', 'mice',
               "Impseq", "Impseqrob", "GMS", 'SeqKNN']

    normals = ['blank', "sum", "max", "center.mean", "center.median", "div.mean", "div.median", "quantiles",
               "quantiles.robust", "vsn", 'lossf', 'TIC', "Rlr", 'MBQN']

    if flag == 1:
        run_all_group1(acq, tool_group1, platform, imputs1, normals, in_files, dataset, true_organism, DE_organism)
    elif flag == 2:
        run_all_group2(acq, tool_group2, platform, imputs1, normals, in_files, dataset, true_organism, DE_organism)
    elif flag == 3:
        run_all_group3(acq, tool_group3, platform, imputs1, normals, in_files, dataset, true_organism, DE_organism)

    elif flag == 4:
        run_Mstats(acq, in_files, platform, dataset, true_organism, DE_organism)
    elif flag == 5:
        run_all_group1(acq, tool_group7, platform, imputs1, normals, in_files, dataset, true_organism, DE_organism)
    elif flag == 6:
        run_all_group1(acq, tool_group4, platform, imputs1, normals, in_files, dataset, true_organism, DE_organism)
    elif flag == 7:
        run_all_group3(acq, tool_group5, platform, imputs1, normals, in_files, dataset, true_organism, DE_organism)
    elif flag == 8:
        run_all_group1(acq, tool_group6, platform, imputs1, normals, in_files, dataset, true_organism, DE_organism)


if __name__ == '__main__':

    platform = sys.argv[1]  # FragPipe Maxquant DIANN spectronaut
    Rscript_pth = sys.argv[2]
    #user_save_fold = sys.argv[3]

    acq = 'DIA'
    root_fold = ''

    root_fold = 'D:/data/benchmark/'
    R_fold = root_fold + 'codes/benchmark_R_scripts_DIA/'
    save_folder = root_fold + 'benchmark_res/' + acq + '/' + platform + '/'
    maxtrix_folder = root_fold + 'data/' + acq + '/' + platform + '/'
    Rscript_path = Rscript_pth

    for data_idx in range(7):
        for flag in [1, 4, 5, 6, 8]:

            if platform == 'FragPipe':
                dataset_info_file = root_fold + 'data/dataset_info/' + acq + '_' + 'Frag.txt'
            elif platform == 'Maxquant':
                dataset_info_file = root_fold + 'data/dataset_info/' + acq + '_' + 'mq.txt'
            elif platform == 'DIANN':
                dataset_info_file = root_fold + 'data/dataset_info/' + acq + '_' + 'diann.txt'
            elif platform == 'spt':
                dataset_info_file = root_fold + 'data/dataset_info/' + acq + '_' + 'spt.txt'

            dataset_info = pd.read_csv(dataset_info_file, sep='\t', header=0)
            dataset = dataset_info['dataset'][data_idx]
            true_organism = dataset_info['organism'][data_idx]
            DE_organism = dataset_info['trueDE'][data_idx]
            true_lgfc = dataset_info['true_lgfc'][data_idx]
            conss = dataset_info['test_contrs'][data_idx]

            main_output_protein_file = '_'#dataset_info['output_folder'][data_idx] + dataset_info['main_output_protein_file'][
                #data_idx]
            main_output_peptide_file = '_' #dataset_info['output_folder'][data_idx] + dataset_info['main_output_peptide_file'][
                #data_idx]

            all_reported_protein_path = maxtrix_folder + dataset_info['all_protein_file'][
                data_idx]

            msstats_design_file = dataset_info['msstat_design_path'][
                data_idx]
            if platform != 'FragPipe':
                evidence_file = dataset_info['output_folder'][data_idx] + dataset_info['main_output_evidence'][
                    data_idx]
                mstats_file = 'NULL'
            else:
                evidence_file = 'NULL'
                mstats_file = dataset_info['output_folder'][data_idx] + 'TOP0/' + dataset_info['mstats_output'][
                    data_idx]

            in_files = {'mstats_file': mstats_file, 'evidence_file': evidence_file,
                        'msstats_design_file': msstats_design_file, 'main_output_protein_file': main_output_protein_file,
                        'main_output_peptide_file': main_output_peptide_file, 'root_fold': root_fold, 'R_fold': R_fold,
                        'save_folder': save_folder, 'maxtrix_folder': maxtrix_folder,
                        'all_candidate': all_reported_protein_path,
                        'Rscript_path': Rscript_path, 'dataset_info_file': dataset_info_file, 'true_lgfc': true_lgfc,
                        'contrasts':conss}

            run_all(acq, in_files, dataset, platform, flag, true_organism, DE_organism)
