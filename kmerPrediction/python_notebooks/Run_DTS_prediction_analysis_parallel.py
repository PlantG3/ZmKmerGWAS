#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#import needed modules
import os
import pandas as pd
pd.set_option('display.max_rows', 200)
import numpy as np
import matplotlib.pyplot as plt
from statsmodels.graphics.gofplots import qqplot
from scipy.stats import boxcox
from sklearn.linear_model import LinearRegression, RidgeCV, Ridge, LassoCV
from datetime import datetime
from sklearn.feature_extraction.text import CountVectorizer, TfidfVectorizer
from scipy.stats import pearsonr
import pickle
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
import multiprocessing
import time

import random
random.seed(12345)

import mkl
mkl.set_num_threads(1)

#required for keras model
import os

os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID" 
os.environ["CUDA_VISIBLE_DEVICES"] = "0"
from tensorflow.keras.models import Sequential
from tensorflow.keras import layers
from tensorflow.keras import backend as K
import gc
import tensorflow as tf
gpus = tf.config.experimental.list_physical_devices('GPU')
if gpus:
  # Restrict TensorFlow to only allocate 1GB of memory on the first GPU
    try:
        tf.config.experimental.set_virtual_device_configuration(
            gpus[0],
            [tf.config.experimental.VirtualDeviceConfiguration(memory_limit=3900)])
        logical_gpus = tf.config.experimental.list_logical_devices('GPU')
        print(len(gpus), "Physical GPUs,", len(logical_gpus), "Logical GPUs")
    except RuntimeError as e:
        print(e)


# In[ ]:


'''
#explore and/or transform phenotype data
def check_normality(numpy_array):
    fig, ax = plt.subplots(1, 2, figsize=(9,4))
    ax[0].hist(numpy_array)
    qqplot(numpy_array, line='s', ax=ax[1])
    plt.show()
    
def max_sqrt(x):
    return np.sqrt(np.max(x+1) - x)
def max_log10(x):
    return np.log10(np.max(x+1) - x)
def max_inverse(x):
    return 1/(np.max(x+1) - x)

all_phenos = pd.read_csv("../data/maize282_phenotype_summary11oct2021.txt", sep="\t", index_col=[0])
for pheno in all_phenos.columns:
    print(pheno)
    sng_pheno = all_phenos[pheno].dropna().copy()
    # for oil, transform with boxcox
    if pheno == "Oil":
        sng_pheno = pd.Series(boxcox(sng_pheno)[0], index=sng_pheno.index)
        all_phenos[pheno] = sng_pheno
    if pheno == "ULA":
        sng_pheno = pd.Series(boxcox(max_sqrt(np.sqrt(all_phenos[pheno].dropna())))[0], index=sng_pheno.index)
        all_phenos[pheno] = sng_pheno
    check_normality(all_phenos[pheno].dropna().copy())
all_phenos.to_csv("../data/maize282_phenotype_transformed11oct2021.txt", sep="\t")
'''
all_phenos = pd.read_csv("../data/maize282_phenotype_transformed11oct2021.txt", sep="\t", index_col=[0])


# In[ ]:


print(str(multiprocessing.cpu_count()) + " CPUs available.")


# In[ ]:


n_cpus_vectorizer = 22
n_cpus_models = 22
print("Using "+str(n_cpus_vectorizer)+ " CPUs for vectorizer.")
print("Using "+str(n_cpus_models)+ " CPUs for models.")

#select which models, vectorizors, and testing schemes
vectorizors = ["CountVectorizer", "TfidfVectorizer"]
models = ["LinearRegression", "Ridge", "Lasso", "NeuralNetwork"]
#models = ["NeuralNetwork"] # run NN part of models not run on cluster

###IN CURRENT SCRIPT splitting can only have one value
splitting = ['Cluster_282_11k']
#splitting = ["Random_Kfolds_11k"]

#phenotype = "Oil"
phenotype = "ULA"
#phenotype = "FT"

#dataSets = ["1M_random_121320_raw_"+phenotype] #for all analyses with randomly chosen k-mers

#dataSets = ["1M_associated_kmeans_102421_FT"]
#dataSets = ["1M_associated_randomFolds_102521_FT"]

#dataSets = ["1M_associated_randomFolds_110121_ULA"]
#dataSets = ["1M_associated_kmeans_110121_ULA"]

#dataSets = ["1M_associated_kmeans_110421_Oil"]
#dataSets = ["1M_associated_randomFolds_110421_Oil"]

include_abandance = [False]
#include_abandance = [True] #, False]

#specify file name for run results to be save
res_folds = datetime.now().strftime('%m%d%Y_%H%M%S')+"_results_folds.csv" #results from each fold

print(res_folds)

#put selections into table of required runs
run_params = [] #will contain run perameters from above for all runs to be performed
for vec in vectorizors:
    for mod in models:
        for spl in splitting:
            for dat in dataSets:
                for inc in include_abandance:
                    run_params.append([vec,mod,spl,dat,inc])
run_params = pd.DataFrame(run_params, columns=["vectorizor", "model", "splitting", "dataSet", "include_abandance"])
run_params = run_params.sort_values(["dataSet","include_abandance","splitting","vectorizor"]).reset_index(drop=True) #change order for most efficient runs


# In[ ]:


print(run_params)


# In[ ]:


def vectorize(run, train_set, test_set282, test_setNAM):
    #vectorize.
    if run["vectorizor"] == "CountVectorizer":
            #count vectorizer
        vectorizer = CountVectorizer()
    elif run["vectorizor"] == "TfidfVectorizer":
        #TFIDF vectorizor
        vectorizer = TfidfVectorizer(min_df = 1 , max_df = 1.0, sublinear_tf=True,use_idf=True)

        #setup vectorizors based on all sentences in dataset (for genomic prediction context we would have these)
    all_sentences = pd.concat([train_set, test_set282, test_setNAM])
    vectorizer.fit(all_sentences["sentence"].values)

    #create vectorized training and testing sets
    X_train = vectorizer.transform(train_set["sentence"].values)
    y_train = train_set["value"]

    X_test282 = vectorizer.transform(test_set282["sentence"].values)
    y_test282 = test_set282["value"]

    if len(test_setNAM)>0:
        X_testNAM = vectorizer.transform(test_setNAM["sentence"].values)
        y_testNAM = test_setNAM["value"]
    else:
        X_testNAM = pd.DataFrame([])
        y_testNAM = pd.Series([])

    #print(X_train.shape, len(y_train), X_test282.shape, len(y_test282), X_testNAM.shape, len(y_testNAM))
    return X_train, y_train, X_test282, y_test282, X_testNAM, y_testNAM


# In[ ]:


def test_regression_continuous_data(run, fold, regression, X_train, y_train, setname):
    pred = regression.predict(X_train)
    obs = y_train.values
    res=pd.DataFrame([y_train.index.tolist(),pred,obs], index=["Taxa","Pred","Obs"]).T
    for col in run.index:
        #print(col, run.loc[col])
        res[col]=run.loc[col]
    res["Fold"]=fold
    res["Set"] = setname
    res = res[['Fold', 'Set', 'vectorizor', 'model', 'splitting', 'dataSet', 'include_abandance', 'Taxa', 'Pred', 'Obs']]
    return res


# In[ ]:


def run_eval_regression_reps(LR, run, fold, X_train, y_train, X_test282, y_test282, X_testNAM, y_testNAM, reps):
    reg_reps=[]
    for rep in range(0, reps):
        print(rep)
        LR.fit(X_train, np.asarray(y_train.values))
        #test model on all datasets
        res=[]
        res.append(test_regression_continuous_data(run, fold, LR, X_train, y_train, setname="Train"))
        res.append(test_regression_continuous_data(run, fold, LR, X_test282, y_test282, setname="Test282"))
        if len(y_testNAM)>0:
            res.append(test_regression_continuous_data(run, fold, LR, X_testNAM, y_testNAM, setname="TestNAM"))
        res = pd.concat(res)
        res["rep"]=rep
        reg_reps.append(res)
    reg_reps = pd.concat(reg_reps)
    return reg_reps


# In[ ]:


def run_eval_keras_in_replicate(run, fold, X_train, y_train, X_test282, y_test282, X_testNAM, y_testNAM, reps):
    keras_reps=[]
    for rep in range(0, reps):
        print(rep)
        model, history = run_keras_sngl_model(X_train, y_train)
        res=[]
        res.append(eval_keras_model(run, fold, model, X_train, y_train, setname="Train"))
        res.append(eval_keras_model(run, fold, model, X_test282, y_test282, setname="Test282"))
        if len(y_testNAM)>0:
            res.append(eval_keras_model(run, fold, model, X_testNAM, y_testNAM, setname="TestNAM"))
        res = pd.concat(res)
        res["rep"]=rep
        keras_reps.append(res)
    keras_reps = pd.concat(keras_reps)
    return keras_reps


# In[ ]:


def run_analysis(run, fold, X_train, y_train, X_test282, y_test282, X_testNAM, y_testNAM):
    #run analysis
    if run["model"] in ["LinearRegression", "Ridge", "Lasso"]:
        #Fit regression
        if run["model"]=="LinearRegression":
            LR = LinearRegression(fit_intercept=True, normalize=False, copy_X=True, n_jobs=1)
        elif run["model"]=="Ridge":
            alphas = [1e-15, 1e-10, 1e-8, 1e-4, 1e-3, 1e-2, 1, 5, 10, 20]
            LR= RidgeCV(alphas, cv=5)
        elif run["model"]=="Lasso":
            LR = LassoCV(cv=5, n_jobs=1, max_iter=2000) #, random_state=12345

        res = run_eval_regression_reps(LR, run, fold, X_train, y_train, X_test282, y_test282, X_testNAM, y_testNAM, reps=1)
        
    elif run["model"]=="NeuralNetwork":
        #print("Neural Net not yet implemented")
        res = run_eval_keras_in_replicate(run, fold, X_train, y_train, X_test282, y_test282, X_testNAM, y_testNAM, reps=10)
    
    return res


# In[ ]:


def run_keras_sngl_model(X_train, y_train):
    #setup model
    K.clear_session()
    gc.collect()
    input_dim = X_train.shape[1]  # Number of features

    model = Sequential()
    #model.add(layers.Dense(100, input_dim=input_dim, activation='relu'))
    #model.add(layers.Dense(50, input_dim=input_dim, activation='relu'))
    model.add(layers.Dense(10, input_dim=input_dim, activation='relu'))
    model.add(layers.Dense(1, activation='linear'))
    model.compile(loss='mse', optimizer='adam')
    #model.summary()

    #train model
    callback = tf.keras.callbacks.EarlyStopping(monitor='val_loss', patience=3)
    history = model.fit(X_train.toarray(), y_train.values, epochs=200,
                        validation_split=0.1,
                        batch_size=10, callbacks=[callback], verbose=0)
    return model, history

def eval_keras_model(run, fold, model, X_train, y_train, setname):
    pred = model.predict(X_train.toarray(), batch_size=10, verbose=False).flatten()
    obs = y_train.values
    res=pd.DataFrame([y_train.index.tolist(), pred,obs], index=["Taxa","Pred","Obs"]).T
    for col in run.index:
        #print(col, run.loc[col])
        res[col]=run.loc[col]
    res["Fold"]=fold
    res["Set"] = setname
    res = res[['Fold', 'Set', 'vectorizor', 'model', 'splitting', 'dataSet', 'include_abandance', 'Taxa','Pred', 'Obs']]
    return res


# In[ ]:


def import_data_by_fold(fold, fold_sets, with_abnd, phe_file, abnd_file):
    #import set data by fold
    sets_in_fold = fold_sets[fold_sets["Fold"]==fold].copy()
    train_set_cltvrs = sets_in_fold[sets_in_fold["Set"]=="282train"]["Taxa"].tolist()
    test_set282_cltvrs = sets_in_fold[sets_in_fold["Set"]=="282test"]["Taxa"].tolist()
    print(abnd_file)
    
    #import abandance and phenotype data
    data_set, _, _ = read_data(phe_file,abnd_file,with_abnd)

    #create training and testing sets
    train_set = data_set.loc[[x for x in train_set_cltvrs if x in data_set.index.tolist()]]
    test_set282 = data_set.loc[[x for x in test_set282_cltvrs if x in data_set.index.tolist()]]
    #No NAM test set exists for this data
    test_setNAM = pd.DataFrame([])
    testNAM_uni_tokens = []
    print(len(train_set), len(test_set282))

    #sainity check
    if len([x for x in train_set.index if x in test_set282.index.tolist()]) != 0:
        print("TRAIN and TEST set CONTAMINATION!!!!!")
        sys.exit(1)
    return train_set, test_set282, test_setNAM


# In[ ]:


def import_data_by_fold_associated_kmers(fold, fold_sets, phe_file, abnd_train, abnd_test, with_abnd=False):
    #import set data by fold
    sets_in_fold = fold_sets[fold_sets["Fold"]==fold].copy()
    train_set_cltvrs = sets_in_fold[sets_in_fold["Set"]=="282train"]["Taxa"].tolist()
    test_set282_cltvrs = sets_in_fold[sets_in_fold["Set"]=="282test"]["Taxa"].tolist()
    
    print(abnd_train)
    print(abnd_test)
    #import abandance and phenotype data by fold
    train_set, train_set_uni_tokens, _ = read_data(phe_file, abnd_file=abnd_train, with_abnd=with_abnd)
    train_set = train_set.loc[[x for x in train_set_cltvrs if x in train_set.index.tolist()]]
    
    test_set282, test282_uni_tokens, _ = read_data(phe_file, abnd_file=abnd_test, with_abnd=with_abnd)
    test_set282 = test_set282.loc[[x for x in test_set282_cltvrs if x in test_set282.index.tolist()]]
    
    test_setNAM = pd.DataFrame([])
    testNAM_uni_tokens = []
    
    #sanity check!
    if len([x for x in train_set.index.tolist() if x not in train_set_cltvrs]) != 0:
        print("STOP!!! Train set does not match cultivar list.")
        sys.exit(1)
    if len([x for x in test_set282.index.tolist() if x not in test_set282_cltvrs]) != 0:
        print("STOP!!! Train set does not match cultivar list.")
        sys.exit(1)
    if len([x for x in train_set.index if x in test_set282.index.tolist()]) != 0:
        print("TRAIN and TEST set CONTAMINATION!!!!!")
        sys.exit(1)
    
    print (len(train_set_cltvrs), len(test_set282_cltvrs))
    print(len(train_set), len(test_set282))
    return train_set, test_set282, test_setNAM


# In[ ]:


#bring in fold splits
if len(run_params["splitting"].unique()) ==1:
    if run_params["splitting"].unique()[0] == 'Random_Kfolds_11k':
        fold_sets = pd.read_csv("../data/random_11k_folds_sets_12Dec2020.csv")
    elif run_params["splitting"].unique()[0] == 'Cluster_282_11k':
        fold_sets = pd.read_csv("../data/kmeans_11k_folds_sets_12Dec2020.csv")
        
include_folds =(0,len(fold_sets["Fold"].unique()))
for fold in fold_sets["Fold"].unique()[include_folds[0]:include_folds[1]]:
    print(fold, len(fold_sets[(fold_sets["Fold"]==fold) & (fold_sets["Set"]=="282train")]),
          len(fold_sets[(fold_sets["Fold"]==fold) & (fold_sets["Set"]=="282test")]),
          len(fold_sets[(fold_sets["Fold"]==fold) & (fold_sets["Set"]=="NAMtest")]))


# In[ ]:


print(run_params)


# In[ ]:


#assigne details for all k-folds vectorizors and sets
filt=None
fold_vec_runs = []
for fold in fold_sets["Fold"].unique()[include_folds[0]:include_folds[1]]:
    for vect in run_params["vectorizor"].unique():
        for aband in run_params["include_abandance"].unique():
            #print(fold)
            if run_params["dataSet"][0]=="1M_associated_kmeans_102421_FT":
                abnd_train = "../data/11fold_kmean_DTS_102421/train_raw/DTS"+str(fold)+"_top1M_kmean_train_raw_KOC.txt"
                abnd_test = "../data/11fold_kmean_DTS_102421/test_raw/DTS"+str(fold)+"_top1M_kmean_test_raw_KOC.txt"
            elif run_params["dataSet"][0]=="1M_associated_randomFolds_102521_FT":
                abnd_train = "../data/11_fold_random_DTS_102521/train_raw/DTS"+str(fold)+"_top1M_random_train_raw_KOC.txt"
                abnd_test = "../data/11_fold_random_DTS_102521/test_raw/DTS"+str(fold)+"_top1M_random_test_raw_KOC.txt"
            
            elif run_params["dataSet"][0]=="1M_associated_randomFolds_110121_ULA":
                abnd_train = "../data/11_fold_random_ULA_110121/train_raw/ULA"+str(fold)+"_top1M_random_train_raw_KOC.txt"
                abnd_test = "../data/11_fold_random_ULA_110121/test_raw/ULA"+str(fold)+"_top1M_random_test_raw_KOC.txt"
            elif run_params["dataSet"][0]=="1M_associated_kmeans_110121_ULA":
                abnd_train = "../data/11_fold_kmean_ULA_110121/train_raw/ULA"+str(fold)+"_top1M_kmean_train_raw_KOC.txt"
                abnd_test = "../data/11_fold_kmean_ULA_110121/test_raw/ULA"+str(fold)+"_top1M_kmean_test_raw_KOC.txt"
            
            elif run_params["dataSet"][0]=="1M_associated_kmeans_110421_Oil":
                abnd_train = "../data/11_fold_kmean_oil_110421/train_raw/oil"+str(fold)+"_top1M_kmean_train_raw_KOC.txt"
                abnd_test = "../data/11_fold_kmean_oil_110421/test_raw/oil"+str(fold)+"_top1M_kmean_test_raw_KOC.txt"
            elif run_params["dataSet"][0]=="1M_associated_randomFolds_110421_Oil":
                abnd_train = "../data/11_fold_random_oil_110421/train_raw/oil"+str(fold)+"_top1M_random_train_raw_KOC.txt"
                abnd_test = "../data/11_fold_random_oil_110421/test_raw/oil"+str(fold)+"_top1M_random_test_raw_KOC.txt"

            elif run_params["dataSet"][0]=="1M_random_121320_raw_"+phenotype:
                abnd_train=""
                abnd_test=""
                #set k-mer abandance file
                abnd_file="../data/maize282.k31.random.1M.KOC.txt"
            fold_vec_runs.append([fold, vect, run_params["splitting"].iloc[0], run_params["dataSet"].iloc[0],
                                  aband, abnd_train, abnd_test])
fold_vec_runs = pd.DataFrame(fold_vec_runs, columns = ["Fold", "vectorizor", "splitting","dataSet","include_abandance", "abnd_train", "abnd_test"])


# In[ ]:


print(fold_vec_runs)
print(fold_vec_runs["abnd_train"].str.split("/", expand=True))
print(fold_vec_runs["abnd_test"].str.split("/", expand=True))


# In[ ]:


def read_data(phe_file, abnd_file, with_abnd):
    #read in data
    if type(phe_file) == pd.core.frame.DataFrame:
        kmer_phe = phe_file.copy()
    else:
        kmer_phe = pd.read_csv(phe_file, sep="\t", index_col=[0]).dropna()
    kmer_abnd = pd.read_csv(abnd_file, sep="\t", index_col=[0])
    
    if type(filt)==int:
        pval_col = [x for x in kmer_abnd.columns if x.split("_")[-1]=="pvals"]
        if len(pval_col)==1:
            kmer_abnd = kmer_abnd.sort_values(pval_col[0]) # sort by p-value
            kmer_abnd = kmer_abnd.iloc[:filt] #take only the top "filt" number
            print("Selecting the top "+str(filt)+"k-mers based on "+pval_col[0])
            print("Max and min pvals:", kmer_abnd[pval_col[0]].max(), kmer_abnd[pval_col[0]].min())
            print("dataset size:", len(kmer_abnd))
        else:
            print("Randomly selecting "+str(filt)+" kmers.")
            rand_sample = random.sample(kmer_abnd.index.tolist(), filt)
            kmer_abnd = kmer_abnd.loc[rand_sample]
    
    #fix column name issues
    kmer_abnd.index.name = "id"
    kmer_abnd.rename(columns={"B73-1":"B73", "B97-1":"B97"}, inplace=True)
    kmer_phe.index.name = "Taxa"
    kmer_phe.columns = ["phe"]
    kmer_phe = kmer_phe.loc[[x for x in kmer_phe.index.tolist() if x in kmer_abnd.columns.tolist()]]
    
    #Create token list
    #with_abnd=True #include multiple copies of each kmer based on its abandance
    token_list=[]
    for geno in [x for x in kmer_phe.index.tolist() if x in kmer_abnd.columns.tolist()]:
        if with_abnd:
            tmp = kmer_abnd[kmer_abnd[geno]>0][geno].reset_index()
            #tmp = ((tmp["id"].str.lower()+" ") * tmp[geno]).str.strip()#.str.split(" ")
            tmp = ((tmp["id"].str.lower()+" ") * tmp[geno].round().astype(int)).str.strip()#.str.split(" ")
        else:
            tmp = kmer_abnd[kmer_abnd[geno]>0].index.str.lower().tolist()
        token_list.append(" ".join(tmp))
    kmer_phe["tokens"]=token_list
    unique_tokens = list(set(kmer_abnd.index.tolist()))
    kmer_phe.columns = ["value","sentence"]
    return kmer_phe, unique_tokens, kmer_abnd


# In[ ]:


def create_save_vec_sets(fold_vec_run, base):
    start = time.time() - base
    out_file = "@".join(fold_vec_run[["Fold","vectorizor","splitting","dataSet","include_abandance"]].astype(str).tolist())+".p"
    #print(out_file)
    if out_file in os.listdir("../data/TMP_vectorizors/"):
        print("File already exsist: "+ out_file)
    else:
        print("Loading Dataset and creating file: "+ out_file)
        if (run_params["dataSet"][0]=="1M_random_121320_raw_"+phenotype):
            train_set, test_set282, test_setNAM = import_data_by_fold(fold_vec_run["Fold"], fold_sets.copy(), with_abnd=fold_vec_run["include_abandance"],
                                                                      phe_file = all_phenos[[run_params["dataSet"][0].split("_")[-1]]].dropna(),
                                                                      abnd_file = abnd_file)
            print(abnd_file)
        else:
            train_set, test_set282, test_setNAM = import_data_by_fold_associated_kmers(fold_vec_run["Fold"], fold_sets.copy(), with_abnd=fold_vec_run["include_abandance"],
                                                                                        phe_file = all_phenos[[run_params["dataSet"][0].split("_")[-1]]].dropna(),
                                                                                        abnd_train = fold_vec_run["abnd_train"],
                                                                                        abnd_test = fold_vec_run["abnd_test"])
            print(fold_vec_run["abnd_train"])
            print(fold_vec_run["abnd_test"])
        print(len(train_set), len(test_set282), len(test_setNAM))

        print("Creating new vectorized training and testing sets. "+ out_file)
        #vectorize data sets
        X_train, y_train, X_test282, y_test282, X_testNAM, y_testNAM = vectorize(fold_vec_run, train_set, test_set282, test_setNAM)
        print(X_train.shape, len(y_train), X_test282.shape, len(y_test282), X_testNAM.shape, len(y_testNAM), out_file)

        print("saving pickle "+ out_file)

        with open("../data/TMP_vectorizors/"+out_file, "wb") as out_p:
            pickle.dump([X_train, y_train, X_test282, y_test282, X_testNAM, y_testNAM], out_p)

    stop = time.time() - base
    return start, stop, out_file


# In[ ]:


def multiprocess(func, args, workers):
    begin_time = time.time()
    with ProcessPoolExecutor(max_workers=workers) as executor:
        res = executor.map(func, args, [begin_time for i in range(len(args))])
    return list(res)

#run all in parallel and save results
args = [fold_vec_runs.loc[x] for x in fold_vec_runs.index]
results = multiprocess(create_save_vec_sets, args, n_cpus_vectorizer)


# In[ ]:


print(run_params)


# In[ ]:


#run each fold_vec_run on it own cpu
def run_save_analysis_sets(fold_vec_run, base):
    start = time.time() - base
    
    ### IMPORT SETS ###
    print("Importing new vectorized training and testing sets.")
    #print(fold_vec_run)
    in_file = "@".join(fold_vec_run[['Fold','vectorizor', 'splitting', 'dataSet', 'include_abandance']].astype(str).tolist())+".p"
    print(in_file)
    with open("../data/TMP_vectorizors/"+in_file, "rb") as f:
        X_train, y_train, X_test282, y_test282, X_testNAM, y_testNAM = pickle.load(f)
    #X_train = X_train[:,:100] #for testing
    #X_test282 = X_test282[:,:100]
    print(X_train.shape, len(y_train), X_test282.shape, len(y_test282), X_testNAM.shape, len(y_testNAM))
    
    ### Run each analysis ###
    print("Running Analysis.")
    #setup run params
    run_params_tmp = run_params[(run_params["vectorizor"]==fold_vec_run["vectorizor"]) &
                                (run_params["splitting"]==fold_vec_run["splitting"]) &
                                (run_params["dataSet"]==fold_vec_run["dataSet"]) &
                                (run_params["include_abandance"]==fold_vec_run["include_abandance"])].copy()
    for run in run_params_tmp.iterrows():
        run = run[1]
        if run["model"] == "NeuralNetwork":
            print("NeuralNetwork cannot be run in parallel. Will be run in series as end.")
            continue
        print(run["model"])
        run_results = run_analysis(run, fold_vec_run["Fold"], X_train, y_train, X_test282, y_test282, X_testNAM, y_testNAM)
        run_results["refreshed"]="vectorizor"
    
        #record and save results
        #check if results file exist
        if res_folds not in os.listdir("../results"):
            run_results.to_csv("../results/"+res_folds, mode="w")
        else:
            run_results.to_csv("../results/"+res_folds, mode="a")

    stop = time.time() - base
    return start, stop


# In[ ]:


if len([x for x in run_params["model"].unique() if x != "NeuralNetwork"]) > 1:
    print("running non-neural network models in parallel.")
    results = multiprocess(run_save_analysis_sets, args, n_cpus_models)


# In[ ]:


if "NeuralNetwork" in run_params["model"].unique():
    print("Running Neural Network models in series.")
    run_params_NN = run_params[run_params["model"]=="NeuralNetwork"].copy().reset_index(drop=True)
    print(run_params_NN)
    
    results_folds=[]
    for fold in fold_sets["Fold"].unique()[include_folds[0]:include_folds[1]]:
        print(fold)
        results_folds=[] #record which data was refreshed
        prev_run=pd.Series(index=run_params_NN.loc[0].index, dtype='str')
        for run in run_params_NN.iterrows():
            refresh=[] #record which data was refreshed
            run=run[1]
            print(run.tolist())

            #import data, and vectorizor
            #vectorize training and testing sets
            if ((run["dataSet"]== prev_run["dataSet"]) and
                (run["include_abandance"]==prev_run["include_abandance"]) and
                (run["splitting"]== prev_run["splitting"]) and 
                (run["vectorizor"]== prev_run["vectorizor"])):
                print("Using previous vectorizor.")
            else:
                print("Importing new vectorized training and testing sets.")
                in_file = "@".join([str(fold)]+run[['vectorizor', 'splitting', 'dataSet', 'include_abandance']].astype(str).tolist())+".p"
                with open("../data/TMP_vectorizors/"+in_file, "rb") as f:
                    X_train, y_train, X_test282, y_test282, X_testNAM, y_testNAM = pickle.load(f)
                print(X_train.shape, len(y_train), X_test282.shape, len(y_test282), X_testNAM.shape, len(y_testNAM))
                refresh.append("vectorizor")

            #run the analysis
            print("Running Analysis.")
            run_results = run_analysis(run, fold, X_train, y_train, X_test282, y_test282, X_testNAM, y_testNAM)
            run_results["refreshed"]=str(refresh)

            #save data
            #record and save results
            #check if results file exist
            if res_folds not in os.listdir("../results"):
                run_results.to_csv("../results/"+res_folds, mode="w")
            else:
                run_results.to_csv("../results/"+res_folds, mode="a")
            results_folds.append(run_results)

            prev_run = run.copy()
    results_folds = pd.concat(results_folds)


# In[ ]:




