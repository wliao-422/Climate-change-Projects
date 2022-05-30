import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import random
from sklearn import preprocessing
from sklearn.metrics import r2_score
from sklearn.linear_model import LinearRegression
from sklearn.linear_model import Ridge
from sklearn.linear_model import Lasso
from sklearn.model_selection import cross_val_score
from sklearn.kernel_ridge import KernelRidge
from sklearn.svm import SVR
from matplotlib.mlab import PCA as mlabPCA
import time

import os
# os.chdir('')
os.getcwd()

dat2 = pd.read_csv("Liao2017_TC_data2_May2017.csv")
# dat3 = pd.read_csv("Liao2017_TC_data3.csv")
# dat4 = pd.read_csv("Liao2017_TC_data4.csv")
# dat5 = pd.read_csv("Liao2017_TC_data5.csv")
dat6 = pd.read_csv("Liao2017_TC_data6_May2017.csv")

###
### Step1: pre-process data for fit
###

# datfit = [dat2, dat3, dat4, dat5]
datfit = [dat2, dat6]
datfit = pd.concat(datfit)
print("datfit shape is",datfit.shape)

random.seed(1)
rowdatfit = list(range(datfit.shape[0]))
print("length of rows is",len(rowdatfit))
tr_row = random.sample(rowdatfit, round((datfit.shape[0]+1)*0.7)) # columns for training.
sort_tr_row = sorted(tr_row, reverse = True)
print("length of training data is",len(sort_tr_row),"which is",len(sort_tr_row)/len(rowdatfit))

# get the row numbers of training data
# get the row numbers of testing data
for i in sort_tr_row:
    rowdatfit.pop(i)

ts_row = rowdatfit
print("train row #, ",len(tr_row), "; test row #, ", len(ts_row),"which is ",len(tr_row)/(len(tr_row)+len(ts_row)))

###
### Step2: perform regression
###

	# 1) linear regression
def linear_reg_pred(datax, datay):
    linearreg = LinearRegression(fit_intercept=True, normalize=True)
    datax_train = datax.iloc[tr_row,:] # select the training x data.
    datay_train = datay.iloc[tr_row] # select the training y data.
    datax_test = datax.iloc[ts_row,:] # select the testing x data.
    datay_test = datay.iloc[ts_row] # select the testing y data
    xval_train = datax_train
    yval_train = datay_train
    xval_test = datax_test
    yval_test = datay_test
	# perform regression
    t0 = time.time()
    linearreg.fit(xval_train,yval_train)
    fit_time = time.time() - t0
    coefs_linear = linearreg.coef_
    y_pred_test=linearreg.predict(xval_test)
    y_pred_train=linearreg.predict(xval_train)
	# calculate validation scores
    rsq_test = r2_score(yval_test, y_pred_test)
    rsq_train = r2_score(yval_train, y_pred_train)
    return rsq_test,rsq_train,coefs_linear, fit_time

datay = datfit['agr']
# datax = datfit[['dbh','LFN_SUN',"SG100C_AVG",'BA','intraBA','FinterBA','NinterBA','F_LFN_BAc','N_LFN_BAc','F_WD_BAc','N_WD_BAc']]
datax = datfit[['dbh','LFN_SUN', 'LFP_SUN', 'LFC_SUN', 'LFC13_SUN', 'SG100C_AVG', 'SEED_DRY', 'DBH_AVG', 'HEIGHT_AVG', 'DIAM_AVG', 'LEAFAREA_AVD', 'LEAFTHCK_AVD', 'LMADISC_AVD', 'LDDISC_AVD', 'AVG_LAMTUF', 'fix', 'BA', 'intraBA', 'interBA', 'LN_BAc', 'LP_BAc', 'LC_BAc', 'd13C_BAc', 'WD_BAc', 'SM_BAc', 'DBH_BAc', 'Hmax_BAc', 'CD_BAc', 'LA_BAc', 'LT_BAc', 'LMA_BAc', 'LD_BAc', 'LTUF_BAc']]
rsq_test_linear,rsq_train_linear,coefs_linear, linear_fit_time = linear_reg_pred(datax,datay)

n_fold = 8
np.mean(cross_val_score(LinearRegression(fit_intercept=True, normalize=True), datax, datay, cv=n_fold))

    # 2) ridge regression
def ridge_reg_pred(datax, datay, alpha):
	# set ridge regression
    ridgereg = Ridge(alpha=alpha, normalize=True)
	# get training and testing data
    datax_train = datax.iloc[tr_row,:] # select the training x data.
    datay_train = datay.iloc[tr_row] # select the training y data.
    datax_train_std = datax_train
    datax_test = datax.iloc[ts_row,:] # select the testing x data.
    datay_test = datay.iloc[ts_row] # select the testing y data.
    datax_test_std = datax_test
    xval_train = datax_train
    yval_train = datay_train
    xval_test = datax_test
    yval_test = datay_test
    t0 = time.time()
    ridgereg.fit(xval_train,yval_train)
    fit_time = time.time() - t0
    coefs_ridge = ridgereg.coef_
    y_pred_test=ridgereg.predict(xval_test)
    y_pred_train=ridgereg.predict(xval_train)
	# calculate validation scores
    rsq_test = r2_score(yval_test, y_pred_test)
    rsq_train = r2_score(yval_train, y_pred_train)
    return rsq_test,rsq_train,coefs_ridge, fit_time

datay = datfit['agr']
# datax = datfit[['dbh','LFN_SUN',"SG100C_AVG",'BA','intraBA','FinterBA','NinterBA','F_LFN_BAc','N_LFN_BAc','F_WD_BAc','N_WD_BAc']]
datax = datfit[['dbh','LFN_SUN', 'LFP_SUN', 'LFC_SUN', 'LFC13_SUN', 'SG100C_AVG', 'SEED_DRY', 'DBH_AVG', 'HEIGHT_AVG', 'DIAM_AVG', 'LEAFAREA_AVD', 'LEAFTHCK_AVD', 'LMADISC_AVD', 'LDDISC_AVD', 'AVG_LAMTUF', 'fix', 'BA', 'intraBA', 'interBA', 'LN_BAc', 'LP_BAc', 'LC_BAc', 'd13C_BAc', 'WD_BAc', 'SM_BAc', 'DBH_BAc', 'Hmax_BAc', 'CD_BAc', 'LA_BAc', 'LT_BAc', 'LMA_BAc', 'LD_BAc', 'LTUF_BAc']]
rsq_test_ridge,rsq_train_ridge,coefs_ridge, ridge_fit_time = ridge_reg_pred(datax,datay,.1)

	# 3) Lasso regression
def lasso_reg_pred(datax, datay, alpha):
    lassoreg = Lasso(alpha=alpha, normalize=True)
	# get training and testing data
    datax_train = datax.iloc[tr_row,:] # select the training x data.
    datay_train = datay.iloc[tr_row] # select the training y data.
    datax_test = datax.iloc[ts_row,:] # select the testing x data.
    datay_test = datay.iloc[ts_row] # select the testing y data.
    xval_train = datax_train
    yval_train = datay_train
    xval_test = datax_test
    yval_test = datay_test
	# perform regression
    t0 = time.time()
    lassoreg.fit(xval_train,yval_train)
    fit_time = time.time() - t0
    coefs_lasso = lassoreg.coef_
    y_pred_test= lassoreg.predict(xval_test)
    y_pred_train = lassoreg.predict(xval_train)
	# get validation scores
    rsq_test = r2_score(yval_test, y_pred_test)
    rsq_train = r2_score(yval_train, y_pred_train)
    return rsq_test,rsq_train,coefs_lasso, fit_time

datay = datfit['agr']
# datax = datfit[['dbh','LFN_SUN',"SG100C_AVG",'BA','intraBA','FinterBA','NinterBA','F_LFN_BAc','N_LFN_BAc','F_WD_BAc','N_WD_BAc']]
datax = datfit[['dbh','LFN_SUN', 'LFP_SUN', 'LFC_SUN', 'LFC13_SUN', 'SG100C_AVG', 'SEED_DRY', 'DBH_AVG', 'HEIGHT_AVG', 'DIAM_AVG', 'LEAFAREA_AVD', 'LEAFTHCK_AVD', 'LMADISC_AVD', 'LDDISC_AVD', 'AVG_LAMTUF', 'fix', 'BA', 'intraBA', 'interBA', 'LN_BAc', 'LP_BAc', 'LC_BAc', 'd13C_BAc', 'WD_BAc', 'SM_BAc', 'DBH_BAc', 'Hmax_BAc', 'CD_BAc', 'LA_BAc', 'LT_BAc', 'LMA_BAc', 'LD_BAc', 'LTUF_BAc']]
rsq_test_lasso,rsq_train_lasso,coefs_lasso,lasso_fit_time = lasso_reg_pred(datax,datay,.05)

# print out results
print("linear r2 for test data = ",rsq_test_linear,", linear r2 for training data = ",rsq_train_linear)
print("ridge r2 for test data = ",rsq_test_ridge,", ridge r2 for training data = ",rsq_train_ridge)
print("lasso r2 for test data = ",rsq_test_lasso,", lasso r2 for training data = ",rsq_train_lasso)

###
### Step3: K-fold validation
###

	# 1) linear regression
datay = datfit['agr']
# datax = datfit[['dbh','LFN_SUN',"SG100C_AVG",'BA','intraBA','FinterBA','NinterBA','F_LFN_BAc','N_LFN_BAc','F_WD_BAc','N_WD_BAc']]
datax = datfit[['dbh','LFN_SUN', 'LFP_SUN', 'LFC_SUN', 'LFC13_SUN', 'SG100C_AVG', 'SEED_DRY', 'DBH_AVG', 'HEIGHT_AVG', 'DIAM_AVG', 'LEAFAREA_AVD', 'LEAFTHCK_AVD', 'LMADISC_AVD', 'LDDISC_AVD', 'AVG_LAMTUF', 'fix', 'BA', 'intraBA', 'interBA', 'LN_BAc', 'LP_BAc', 'LC_BAc', 'd13C_BAc', 'WD_BAc', 'SM_BAc', 'DBH_BAc', 'Hmax_BAc', 'CD_BAc', 'LA_BAc', 'LT_BAc', 'LMA_BAc', 'LD_BAc', 'LTUF_BAc']]

#linearreg = LinearRegression(fit_intercept=True, normalize=True)
#ridgereg = Ridge(alpha=alpha, normalize=True)
#lassoreg = Lasso(alpha=alpha, normalize=True)

alphas = np.logspace(-4, 1, 30)
n_fold = 8

scores_ridge = list()
scores_std_ridge = list()
for alpha in alphas:
	print(alpha)
	this_score = cross_val_score(Ridge(alpha=alpha, normalize=True), datax, datay, cv=n_fold)
	scores_ridge.append(np.mean(this_score))
	scores_std_ridge.append(np.std(this_score))

scores_ridge, scores_std_ridge = np.array(scores_ridge),np.array(scores_std_ridge)

scores_lasso = list()
scores_std_lasso = list()
for alpha in alphas:
	print(alpha)
	this_score = cross_val_score(Lasso(alpha=alpha, normalize=True), datax, datay, cv=n_fold)
	scores_lasso.append(np.mean(this_score))
	scores_std_lasso.append(np.std(this_score))

scores_lasso,scores_std_lasso = np.array(scores_lasso),np.array(scores_std_lasso)

# plotting
plt.figure()
plt.semilogx(alphas, scores_ridge,color='k')
std_error_ridge = scores_std_ridge/np.sqrt(n_fold)
plt.semilogx(alphas,scores_ridge+std_error_ridge,'k:')
plt.semilogx(alphas,scores_ridge-std_error_ridge,'k:')
plt.ylabel('Cross-validation score +/- std error')
plt.xlabel('Alpha')
plt.suptitle('Ridge regression')
plt.axhline(np.max(scores_ridge),linestyle='--',color='r')
plt.xlim(alphas[0],alphas[-1])
plt.ylim(0,0.5)
plt.savefig('Ridge regression CV.png')

plt.figure()
plt.semilogx(alphas, scores_lasso,color='k')
std_error_lasso = scores_std_lasso/np.sqrt(n_fold)
plt.semilogx(alphas,scores_lasso+std_error_lasso,'k:')
plt.semilogx(alphas,scores_lasso-std_error_lasso,'k:')
plt.ylabel('Cross-validation score +/- std error')
plt.xlabel('Alpha')
plt.suptitle('Lasso regression')
plt.axhline(np.max(scores_lasso),linestyle='--',color='r')
plt.xlim(alphas[0],alphas[-1])
plt.ylim(0,0.5)
plt.savefig('Lasso regression CV.png')

###
### Step4: Ridge regression with SVD
###

alphas = np.logspace(-4, 1, 30)
n_fold = 8

scores_ridge_svd = list()
scores_std_ridge_svd = list()
for alpha in alphas:
	print(alpha)
	this_score = cross_val_score(Ridge(alpha=alpha, normalize=True, solver="svd"), datax, datay, cv=n_fold)
	scores_ridge_svd.append(np.mean(this_score))
	scores_std_ridge_svd.append(np.std(this_score))

scores_ridge_svd, scores_std_ridge_svd = np.array(scores_ridge_svd),np.array(scores_std_ridge_svd)

def ridge_reg_pred(datax, datay, alpha):
	# set ridge regression
    ridgereg = Ridge(alpha=alpha, normalize=True, solver='svd')
	# get training and testing data
    datax_train = datax.iloc[tr_row,:] # select the training x data.
    datay_train = datay.iloc[tr_row] # select the training y data.
    datax_test = datax.iloc[ts_row,:] # select the testing x data.
    datay_test = datay.iloc[ts_row] # select the testing y data.
    xval_train = datax_train
    yval_train = datay_train
    xval_test = datax_test
    yval_test = datay_test
	# perform regression
    t0 = time.time()
    ridgereg.fit(xval_train,yval_train)
    fit_time = time.time() - t0
    coefs_ridge = ridgereg.coef_
    y_pred_test=ridgereg.predict(xval_test)
    y_pred_train=ridgereg.predict(xval_train)
	# calculate validation scores
    rsq_test = r2_score(yval_test, y_pred_test)
    rsq_train = r2_score(yval_train, y_pred_train)
    return rsq_test,rsq_train,coefs_ridge, fit_time

datay = datfit['agr']
# datax = datfit[['dbh','LFN_SUN',"SG100C_AVG",'BA','intraBA','FinterBA','NinterBA','F_LFN_BAc','N_LFN_BAc','F_WD_BAc','N_WD_BAc']]
datax = datfit[['dbh','LFN_SUN', 'LFP_SUN', 'LFC_SUN', 'LFC13_SUN', 'SG100C_AVG', 'SEED_DRY', 'DBH_AVG', 'HEIGHT_AVG', 'DIAM_AVG', 'LEAFAREA_AVD', 'LEAFTHCK_AVD', 'LMADISC_AVD', 'LDDISC_AVD', 'AVG_LAMTUF', 'fix', 'BA', 'intraBA', 'interBA', 'LN_BAc', 'LP_BAc', 'LC_BAc', 'd13C_BAc', 'WD_BAc', 'SM_BAc', 'DBH_BAc', 'Hmax_BAc', 'CD_BAc', 'LA_BAc', 'LT_BAc', 'LMA_BAc', 'LD_BAc', 'LTUF_BAc']]
rsq_test_ridge_svd,rsq_train_ridge_svd,coefs_ridge_svd,ridge_svd_fit_time = ridge_reg_pred(datax,datay,.1)

plt.figure()
plt.semilogx(alphas, scores_ridge_svd,color='k')
std_error_ridge_svd = scores_std_ridge_svd/np.sqrt(n_fold)
plt.semilogx(alphas,scores_ridge_svd+std_error_ridge_svd,'k:')
plt.semilogx(alphas,scores_ridge_svd-std_error_ridge_svd,'k:')
plt.ylabel('Cross-validation score +/- std error')
plt.xlabel('Alpha')
plt.suptitle('Ridge regression with SVD')
plt.axhline(np.max(scores_ridge_svd),linestyle='--',color='r')
plt.xlim(alphas[0],alphas[-1])
plt.ylim(0,0.5)
plt.savefig('Ridge regression CV with SVD.png')

###
### Step5: Kernal ridge regression # This did not work.
###

def kridge_reg_pred(datax, datay, alpha, gamma):
	# set ridge regression
	kridgereg = KernelRidge(alpha=alpha, gamma=gamma, kernel='rbf')
	# get training and testing data
	datax_std = datax
	for i in datax.columns.values:
		datax_std.loc[:,i] = preprocessing.scale(datax.loc[:,i])
	datax_train = datax_std.iloc[tr_row,:] # select the training x data.
	datay_train = datay.iloc[tr_row] # select the training y data.
	datax_test = datax_std.iloc[ts_row,:] # select the testing x data.
	datay_test = datay.iloc[ts_row] # select the testing y data.
	xval_train = datax_train
	yval_train = datay_train
	xval_test = datax_test
	yval_test = datay_test
	# perform regression
	kridgereg.fit(xval_train,yval_train)
	coefs_kridge = kridgereg.coef_
	y_pred_test=kridgereg.predict(xval_test)
	y_pred_train=kridgereg.predict(xval_train)
	# calculate validation scores
	rsq_test = r2_score(yval_test, y_pred_test)
	rsq_train = r2_score(yval_train, y_pred_train)
	return rsq_test,rsq_train,coefs_kridge

datay = datfit['agr']
# datax = datfit[['dbh','LFN_SUN',"SG100C_AVG",'BA','intraBA','FinterBA','NinterBA','F_LFN_BAc','N_LFN_BAc','F_WD_BAc','N_WD_BAc']]
datax = datfit[['dbh','LFN_SUN', 'LFP_SUN', 'LFC_SUN', 'LFC13_SUN', 'SG100C_AVG', 'SEED_DRY', 'DBH_AVG', 'HEIGHT_AVG', 'DIAM_AVG', 'LEAFAREA_AVD', 'LEAFTHCK_AVD', 'LMADISC_AVD', 'LDDISC_AVD', 'AVG_LAMTUF', 'fix', 'BA', 'intraBA', 'interBA', 'LN_BAc', 'LP_BAc', 'LC_BAc', 'd13C_BAc', 'WD_BAc', 'SM_BAc', 'DBH_BAc', 'Hmax_BAc', 'CD_BAc', 'LA_BAc', 'LT_BAc', 'LMA_BAc', 'LD_BAc', 'LTUF_BAc']]
rsq_test_ridge,rsq_train_ridge,coefs_ridge = kridge_reg_pred(datax,datay,.1,.1)


alphas = np.logspace(-4, 1, 30)
n_fold = 8

scores_kridge = list()
scores_std_kridge = list()

datay = datfit['agr']
# datax = datfit[['dbh','LFN_SUN',"SG100C_AVG",'BA','intraBA','FinterBA','NinterBA','F_LFN_BAc','N_LFN_BAc','F_WD_BAc','N_WD_BAc']]
datax = datfit[['dbh','LFN_SUN', 'LFP_SUN', 'LFC_SUN', 'LFC13_SUN', 'SG100C_AVG', 'SEED_DRY', 'DBH_AVG', 'HEIGHT_AVG', 'DIAM_AVG', 'LEAFAREA_AVD', 'LEAFTHCK_AVD', 'LMADISC_AVD', 'LDDISC_AVD', 'AVG_LAMTUF', 'fix', 'BA', 'intraBA', 'interBA', 'LN_BAc', 'LP_BAc', 'LC_BAc', 'd13C_BAc', 'WD_BAc', 'SM_BAc', 'DBH_BAc', 'Hmax_BAc', 'CD_BAc', 'LA_BAc', 'LT_BAc', 'LMA_BAc', 'LD_BAc', 'LTUF_BAc']]
datax_nm = datax
for i in datax.columns.values:
	datax_nm.loc[:,i] = preprocessing.scale(datax.loc[:,i])

for alpha in alphas:
	print(alpha)
	this_score = cross_val_score(KernelRidge(alpha=alpha), datax_nm, datay, cv=n_fold)
	scores_kridge.append(np.mean(this_score))
	scores_std_kridge.append(np.std(this_score))

scores_kridge, scores_std_kridge = np.array(scores_kridge),np.array(scores_std_kridge)

###
### Step6: radial-basis function support vector regression
###
datay = datfit['agr']
# datax = datfit[['dbh','LFN_SUN',"SG100C_AVG",'BA','intraBA','FinterBA','NinterBA','F_LFN_BAc','N_LFN_BAc','F_WD_BAc','N_WD_BAc']]
datax = datfit[['dbh','LFN_SUN', 'LFP_SUN', 'LFC_SUN', 'LFC13_SUN', 'SG100C_AVG', 'SEED_DRY', 'DBH_AVG', 'HEIGHT_AVG', 'DIAM_AVG', 'LEAFAREA_AVD', 'LEAFTHCK_AVD', 'LMADISC_AVD', 'LDDISC_AVD', 'AVG_LAMTUF', 'fix', 'BA', 'intraBA', 'interBA', 'LN_BAc', 'LP_BAc', 'LC_BAc', 'd13C_BAc', 'WD_BAc', 'SM_BAc', 'DBH_BAc', 'Hmax_BAc', 'CD_BAc', 'LA_BAc', 'LT_BAc', 'LMA_BAc', 'LD_BAc', 'LTUF_BAc']]
datax_nm = datax
for i in datax.columns.values:
	datax_nm.loc[:,i] = preprocessing.scale(datax.loc[:,i])

random.seed(2)
size = 1000 # len(tr_row) + len(ts_row)
tr_row_sel = random.sample(tr_row, round(size*.7))
ts_row_sel = random.sample(ts_row, round(size*.3))

svr = SVR(C=300, epsilon=.1, cache_size = 100, kernel='linear', verbose=True)
t0 = time.time()
svr.fit(datax_nm.iloc[tr_row_sel,:],datay.iloc[tr_row_sel])
svr_fit_time = time.time() - t0
coefs_svr = svr.coef_
svr_fit_time

y_pred = svr.predict(datax_nm.iloc[tr_row_sel,:])
r2_score(datay.iloc[tr_row_sel],y_pred)
y_pred_test = svr.predict(datax_nm.iloc[ts_row_sel,:])
r2_score(datay.iloc[ts_row_sel],y_pred_test)

n_fold = 8
np.mean(cross_val_score(SVR(C=300, kernel='rbf'), datax_nm.iloc[tr_row_sel,:], datay.iloc[tr_row_sel], cv=n_fold)

scores_svr=list()
scores_std_svr=list()
Cs = np.logspace(-4,4,30)
for C in Cs:
	print(C)
	this_score = cross_val_score(SVR(C=C, kernel='rbf'), datax_nm.iloc[tr_row_sel,:], datay.iloc[tr_row_sel], cv=n_fold)
	scores_svr.append(np.mean(this_score))
	scores_std_svr.append(np.std(this_score))

scores_svr, scores_std_svr = np.array(scores_svr),np.array(scores_std_svr)

n_fold = 8
scores_svr=list()
scores_std_svr=list()
Cs = np.logspace(-4,4,30)
for C in Cs:
	print(C)
	this_score = cross_val_score(SVR(C=C, kernel='linear'), datax_nm.iloc[tr_row_sel,:], datay.iloc[tr_row_sel], cv=n_fold)
	scores_svr.append(np.mean(this_score))
	scores_std_svr.append(np.std(this_score))

scores_svr, scores_std_svr = np.array(scores_svr),np.array(scores_std_svr)

plt.figure()
plt.semilogx(Cs, scores_svr,color='k')
std_error_svr = scores_std_svr/np.sqrt(n_fold)
plt.semilogx(Cs,scores_svr+std_error_svr,'k:')
plt.semilogx(Cs,scores_svr-std_error_svr,'k:')
plt.ylabel('Cross-validation score +/- std error')
plt.xlabel('C')
plt.suptitle('rbf SVR')
plt.axhline(np.max(scores_svr),linestyle='--',color='r')
plt.xlim(Cs[0],Cs[-1])
plt.ylim(0,0.5)
plt.savefig('rbf SVR.png')

###
### Step7: principal component analysis
###

datax_nm_matrix = datax_nm.as_matrix()
mlab_pca = mlabPCA(datax_nm_matrix)

plt.plot(mlab_pca.Y[:,0],mlab_pca.Y[:,1],'o', markersize=7, color='blue', alpha=0.5)
plt.xlabel('First component')
plt.ylabel('Second component')
plt.title('Principal component analysis for 14 functional traits')

###
### Step8: print out the coefficients
###

coefs = [coefs_linear, coefs_lasso, coefs_ridge, coefs_ridge_svd, coefs_svr[0,:]]
dat_coefs = pd.DataFrame(coefs)
dat_coefs.columns = datax_nm.columns.values
dat_coefs.rows = ['linear','lasso','ridge','ridge with svd','svr']
