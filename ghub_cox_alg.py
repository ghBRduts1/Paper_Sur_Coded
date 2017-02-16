#!/usr/bin/env python
import numpy as np
import scipy
from scipy import stats
from scipy.stats import kstest
from rpy2 import robjects as ro
from subprocess import call
import os
ro.r('library(survival)')

path = "####Path to directory where are the intermediate files are located"

###File with counts of number of times a given lncRNA is significantly associated with survival in the survival algorithm
conc_cox = open(path+"/sig_cox_coefficients_COUNTS_lncRNA_09:10:38.txt")

###Set threshold for how many times a given lncRNA is associated with survival
cox_lncs = [x.split()[0] for x in conc_cox if int(x.strip().split()[1]) <= -80 or int(x.strip().split()[1]) >= 80]

###File with median cox coefficient for each lncRNA
coeff_cox = open(path+"/median_cox_coefficients_lncRNA_09:10:38.txt")

coeff_cox_dic = {x.split()[0]:float(x.strip().split()[1]) for x in coeff_cox if x.split()[0] in cox_lncs}


test_expr = path+"/training_expression_matrix_09:10:38.txt"
valid_expr = path+"/validation_expression_matrix_09:10:38.txt"

clinical = open("###PATH to Clincical Data###/LGG_clinical_data.txt")
clinical.readline()

cox_clin = {}

###Converts noncontinuous pt variables to binary. Later used from cox regression
death_dic = {'Alive':0,'Dead':1}
sex_dic = {'FEMALE':0,'MALE':1}
grade_dic = {'G2':0,'G3':1}
idh = {'wt':0,'mut':1}

###Cox_clin value order:    Vital Status, Time, Grade, Gender, Age, IDH_mut Status
for i in clinical:
    x = i.strip().split("\t")
    if '[Not Available]' not in x[1:-2]:
        cox_clin[x[0]] = [death_dic[x[2]],int(x[3]),grade_dic[x[4]],sex_dic[x[5]],int(x[6]),idh[x[1]]]

printable = ['score','grade','sex','age','idh']

def cox_algorithm(data,name):
    sep_set = open(data)
    header = sep_set.readline()
    test_expr_dict = {}
    pt_id_holder = {}
    for i in sep_set:
        x = i.strip().split("\t")
        if x[0] in cox_lncs:
            for j in x[1:]:
                test_expr_dict[x[0]] = test_expr_dict.get(x[0],[]) + [float(j)]
    tr_sig = {i:0 for i in range(0,len(header.strip().split("\t")[1:]))}
    for i in test_expr_dict:
        mean = np.mean(test_expr_dict[i])
        median = np.median(test_expr_dict[i])
        std = np.std(test_expr_dict[i])
        up = median + std
        down = median - std
        for c,j in enumerate(test_expr_dict[i]):
            if j >= up:
                tr_sig[c] = tr_sig[c]+(coeff_cox_dic[i])
            elif j <= down:
                tr_sig[c] = tr_sig[c]-(coeff_cox_dic[i])
    track_clin_dict = {}
    for c,i in enumerate(header.strip().split("\t")[1:]):
        if i in cox_clin:
            track_clin_dict[c] = cox_clin[i]
            pt_id_holder[c] = i
    cox_reg = []
    for index in tr_sig:
        cox_reg.append([tr_sig[index],track_clin_dict[index]])
    ro.globalenv['score']=ro.FloatVector([ii[0] for ii in cox_reg])
    ro.globalenv['died']=ro.IntVector([ii[1][0] for ii in cox_reg])
    ro.globalenv['times']=ro.IntVector([int(ii[1][1]) for ii in cox_reg])
    ro.globalenv['grade']=ro.IntVector([ii[1][2] for ii in cox_reg])
    ro.globalenv['sex']=ro.IntVector([ii[1][3] for ii in cox_reg])
    ro.globalenv['age']=ro.IntVector([int(ii[1][4]) for ii in cox_reg])
    ro.globalenv['idh']=ro.IntVector([ii[1][5] for ii in cox_reg])
    res=ro.r('coxph(Surv(times,died) ~ score + grade + sex + age + idh)')
    print "\n\n\n"+name+"\n"
    for entry in str(res).split('\n'):
        try:
            if entry.split()[0] in printable:
                print entry
        except:
            pass
    follow = [[i,tr_sig[i]] for i in tr_sig]
    follow.sort(key=lambda x: x[1])
    wr = open(path+"/"+name+"_survival_sep.txt","w")
    wr.write("PtName\tVstatus\tDays\tGroup\n")
    for i in follow:
        if i[1] > 0:
            wr.write(pt_id_holder[i[0]]+"\t"+str(track_clin_dict[i[0]][0])+"\t"+str(track_clin_dict[i[0]][1])+"\tGroupNegProg\n")
        elif i[1] < 0:
            wr.write(pt_id_holder[i[0]]+"\t"+str(track_clin_dict[i[0]][0])+"\t"+str(track_clin_dict[i[0]][1])+"\tGroupPosProg\n")
    print follow
    wr.close()


cox_algorithm(test_expr,"TEST")
cox_algorithm(valid_expr,"VALIDATION")


