#!/usr/bin/env python
import numpy as np
import scipy
from scipy import stats
from scipy.stats import kstest
from rpy2 import robjects as ro
from subprocess import call
import random
from random import shuffle
import datetime
import os
ro.r('library(survival)')

path1 = "###Path to input files###"

expr = open(path1+"/filtered_lncMatrix.txt")
header = expr.readline()

clinical = open(path1+"/LGG_clinical_data.txt")
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
        if int(x[3]) > 0:
            cox_clin[x[0]] = [death_dic[x[2]],int(x[3]),grade_dic[x[4]],sex_dic[x[5]],int(x[6]),idh[x[1]]]
        elif int(x[3]) <= 0:
            print i,

track_clin_dict = {}

tracker = []

for c,i in enumerate(header.strip().split("\t")[1:]):
    if i in cox_clin:
        tracker.append(c)
        track_clin_dict[c] = cox_clin[i]

expr_dict = {}
temp_dict = {}
expr_final = {}

###Sets expression thresholds for downstream analysis
###Range of variables is to remove non-primary tumors from the survival analysis, consistant with the previous step
temp_med = []
for i in expr:
    x = i.strip().split("\t")    
    n = 0
    for c,j in enumerate(x[1:]):
        expr_dict[x[0]] = expr_dict.get(x[0],[]) + [float(j)]
        if c in tracker:
            if float(j) == 0:
                n+=1
            temp_dict[x[0]] = temp_dict.get(x[0],[]) + [float(j)]
            if n > 130: ###Maximum number of zeros allowed
                expr_dict.pop(x[0], None)
                break
        else:
            pass
    temp_med.append(len(expr_dict[x[0]]))

fup = int(np.median(temp_med))
remove = []
cox_coeff = {}
c_pos = {}
c_neg = {}
###Sets expression thresholds for downstream analysis
for i in expr_dict:
    if np.median(expr_dict[i]) >= 1.0:          
        if np.mean(expr_dict[i]) >= 1.0:
            if len(expr_dict[i]) == fup:
                expr_final[i]= expr_dict[i]
                cox_coeff[i] = []
                c_pos[i] = 0
                c_neg[i] = 0
            else:
                remove.append(i)


###Creates list of column indexes used to assign 60% of LGG patients to the test group
testing_num = int(len(tracker)*0.6)
sub_testing_num = int(len(tracker)*0.4)

temp = tracker
shuffle(temp)
random_assign = temp[:testing_num]
random_assign.sort()


###WRITE OUTPUT FILES: TESTING AND VALIDATION EXPRESSION MATRICES
###WRITE OUTPUT FILES: TESTING AND VALIDATION PATIENT INFORMATION
###create timestamp for output files
x = str(datetime.datetime.now())
path = "###PATH to the desired directory for the outhput files###"+x.split()[0]+"/"
y = x.split()[1].split('.')[0]

try: 
    os.makedirs(path)
except OSError:
    if not os.path.isdir(path):
        raise
    
    
print "Iterating"

wcount = open(path+"sig_cox_iteration_numbers.txt","w")

iternum = 0
while iternum < 100:
    count=0
    shuffle(random_assign)
    sub_random_assign = random_assign[:sub_testing_num]
    cox_num = 0
    for i in expr_final:
        cox_reg = []
        for j in sub_random_assign:
                cox_reg.append([i,expr_final[i][j],track_clin_dict[j]])
        data=[ii[1] for ii in cox_reg]
        ro.globalenv['expression']=ro.FloatVector(data)
        ## Perform inverse normal transformation
        res=ro.r('round(qnorm((rank(expression, na.last="keep")-0.5)/sum(!is.na(expression))), digit=5)')
        ## Parse the result
        inverse_norm=[float(ii) for ii in str(res).split() if '[' not in ii]
        ## Prepare the variables for rpy2
        ro.globalenv['gene']=ro.FloatVector(inverse_norm)
        ro.globalenv['died']=ro.IntVector([ii[2][0] for ii in cox_reg])
        ro.globalenv['times']=ro.IntVector([ii[2][1] for ii in cox_reg])
        ro.globalenv['grade']=ro.IntVector([ii[2][2] for ii in cox_reg])
        ro.globalenv['sex']=ro.IntVector([ii[2][3] for ii in cox_reg])
        ro.globalenv['age']=ro.IntVector([ii[2][4] for ii in cox_reg])
        ro.globalenv['idh']=ro.IntVector([ii[2][5] for ii in cox_reg])
        res=ro.r('coxph(Surv(times,died) ~ gene + grade + sex + age + idh)')
        ## Parse results
        for entry in str(res).split('\n'):
            try:
                if entry.split()[0]=='gene':
                    if float(entry.split()[-1]) <= 0.05:
                        count+=1
                        if float(entry.split()[1]) > 0:
                            c_pos[i] = c_pos[i]+1
                            cox_coeff[i].append(float(entry.split()[1]))
                        else:
                            c_neg[i] = c_neg[i]-1
                            cox_coeff[i].append(float(entry.split()[1]))
                    break
            except:
                pass
    iternum+=1
    wcount.write(str(count)+"\n")
    print "Iteration: "+str(iternum)

wcount.close()

c_concensus = {}
for i in c_pos:
    c_concensus[i] = c_pos[i] + c_neg[i]


wc = open(path+"sig_cox_coefficients_COUNTS_lncRNA_"+y+".txt","w")

for i in c_concensus:
    wc.write(i+"\t"+str(c_concensus[i])+"\n")

wc.close()


wm = open(path+"median_cox_coefficients_lncRNA_"+y+".txt","w")

for i in cox_coeff:
    if len(cox_coeff[i]) > 1:
        wm.write(i+"\t"+str(np.median(cox_coeff[i]))+"\n")

wm.close()


wmn = open(path+"mean_cox_coefficients_lncRNA_"+y+".txt","w")

for i in cox_coeff:
    if len(cox_coeff[i]) > 1:
        wmn.write(i+"\t"+str(np.mean(cox_coeff[i]))+"\n")

wmn.close()


test = path+"training_expression_matrix_"+y+".txt"


w = open(test,"w")
w.write("lncRNA_ID")

random_assign.sort()
for i in random_assign:
    w.write("\t"+header.strip().split("\t")[int(i)+1])

w.write("\n")
for i in expr_final:
    w.write(i)
    for j in random_assign:
        w.write("\t"+str(expr_final[i][j]))
    w.write("\n")

w.close()


    
##MAKES LIST OF ALL SAMPLES NOT INCLUDED IN THE TESTING SET i.e. VALIDATION PT SET
validation = []
tracker.sort()
for i in tracker:
    if i not in random_assign:
        validation.append(i)

valid = path+"validation_expression_matrix_"+y+".txt"

v = open(valid,"w")
v.write("lncRNA_ID")
for i in validation:
    v.write("\t"+header.strip().split("\t")[int(i)+1])

v.write("\n")
for i in expr_final:
    v.write(i)
    for j in validation:
        v.write("\t"+str(expr_final[i][j]))
    v.write("\n")

v.close()
