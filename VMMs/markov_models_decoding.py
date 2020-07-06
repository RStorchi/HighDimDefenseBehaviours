
import vomm
import numpy as np
import csv
import itertools
import pandas as pd
import xlsxwriter

#retrieve the index of the 516 stimuli (i.e. loom, sound, flash)
indf = open('index_stimuli.csv','r')
ind_vector_raw = indf.read()
ind_vector = ind_vector_raw[:-1].split(',')
ind_vector = np.array(ind_vector).astype(int)

stimulus_0 = np.where(ind_vector == 0)
stimulus_0 = stimulus_0[0]
stimulus_1 = np.where(ind_vector == 1)
stimulus_1 = stimulus_1[0]
stimulus_2 = np.where(ind_vector == 2)
stimulus_2 = stimulus_2[0]
stimulus_idx = [list(stimulus_0), list(stimulus_1), list(stimulus_2)]
num_runs = 5 #number of runs to average classification accuracy

best_params = []
for it in range(num_runs):
    for seq_id in range(1,7):
        for nclus in range(2,11):
            # load the specific sequence generate by Matlab file
            csvf = open('training_set' + str(nclus) + '_seq' + str(seq_id) + '.csv','r')
            training_data_raw = csv.reader(csvf)
            training_data = np.array(list(training_data_raw))
            training_data_0 = training_data[stimulus_idx[0]]
            training_data_1 = training_data[stimulus_idx[1]]
            training_data_2 = training_data[stimulus_idx[2]]

            dims_0 = training_data_0.shape
            dims_1 = training_data_1.shape
            dims_2 = training_data_2.shape
            cut = 0.9 #training/test set ratio for cross-validation
            arr_0 = np.arange(dims_0[0])
            arr_1 = np.arange(dims_1[0])
            arr_2 = np.arange(dims_2[0])
            np.random.shuffle(arr_0)
            np.random.shuffle(arr_1)
            np.random.shuffle(arr_2)

            #randomly split sequences into test and train sets according to the cut ratio
            cutoff_0 = int(np.floor(cut*dims_0[0]))
            cutoff_1 = int(np.floor(cut*dims_1[0]))
            cutoff_2 = int(np.floor(cut*dims_2[0]))
            pos_train_0 = arr_0[0:cutoff_0]
            pos_train_1 = arr_1[0:cutoff_1]
            pos_train_2 = arr_2[0:cutoff_2]
            pos_test_0 = arr_0[cutoff_0:]
            pos_test_1 = arr_1[cutoff_1:]
            pos_test_2 = arr_2[cutoff_2:]

            test_set_0 = training_data_0[pos_test_0,:]
            training_set_0 = training_data_0[pos_train_0,:]
            test_set_1 = training_data_1[pos_test_1,:]
            training_set_1 = training_data_1[pos_train_1,:]
            test_set_2 = training_data_2[pos_test_2,:]
            training_set_2 = training_data_2[pos_train_2,:]

            dims_train_0 = training_set_0.shape
            dims_test_0 = test_set_0.shape
            dims_train_1 = training_set_1.shape
            dims_test_1 = test_set_1.shape
            dims_train_2 = training_set_2.shape
            dims_test_2 = test_set_2.shape

            #convert the sequences into character vectors
            train_seq_0 = []
            for row in range(dims_train_0[0]):
                train = [int(x) for x in training_set_0[row,:]]
                train.append(0)
                train_seq_0 = list(itertools.chain.from_iterable([train_seq_0, train]))

            train_seq_1 = []
            for row in range(dims_train_1[0]):
                train = [int(x) for x in training_set_1[row,:]]
                train.append(0)
                train_seq_1 = list(itertools.chain.from_iterable([train_seq_1, train]))

            train_seq_2 = []
            for row in range(dims_train_2[0]):
                train = [int(x) for x in training_set_2[row,:]]
                train.append(0)
                train_seq_2 = list(itertools.chain.from_iterable([train_seq_2, train]))

            test_set_seq = np.append(test_set_0, np.append(test_set_1, test_set_2, axis = 0), axis = 0)
            test_set_shape = test_set_seq.shape
            test_labels = np.concatenate((np.zeros(int(test_set_shape[0]/3), dtype = int), np.ones(int(test_set_shape[0]/3), dtype = int),\
                                          2*np.ones(int(test_set_shape[0]/3), dtype = int)))


            for mem in range(dims_train_0[1]+1):
                #create the VMM models for each stimulus type
                model_0 = vomm.ppm()
                model_1 = vomm.ppm()
                model_2 = vomm.ppm()

                #train the VMM models for each stimulus type
                model_0.fit(train_seq_0, d = mem)
                model_1.fit(train_seq_1, d = mem)
                model_2.fit(train_seq_2, d = mem)

                #estimate the classification accuracy by maximum likelihood maximization
                l_0 = []
                l_1 = []
                l_2 = []
                accuracy = 0
                for row in range(test_set_shape[0]):
                    test = [int(x) for x in test_set_seq[row,:]]
                    l_0 = model_0.logpdf(test)
                    l_1 = model_1.logpdf(test)
                    l_2 = model_2.logpdf(test)
                    best = np.where([l_0, l_1, l_2] == np.amax([l_0, l_1, l_2]))
                    if best[0][0] == test_labels[row]:
                        accuracy = accuracy + 1
                accuracy = accuracy/float(test_set_shape[0])

                tmp = (it, seq_id, nclus, mem, accuracy, model_0, model_1, model_2, training_data)
                best_params.append(tmp)

#inspect results and collect parameters
acc = []
itt = []
seqid = []
mems = []
ncluseq = []
for (it, seq_id, nclus, mem, accuracy, model_0, model_1, model_2, training_data) in best_params:
    acc.append(accuracy)
    itt.append(it)
    seqid.append(seq_id)
    mems.append(mem)
    ncluseq.append(nclus)

cross_accuracy = np.zeros(int(len(acc)/num_runs))
for i in range(num_runs):
    idx = np.where(np.array(itt) == i)
    for ii in range(len(idx[0])):
        cross_accuracy[ii] = cross_accuracy[ii] + acc[(i*len(idx[0]))+ii]
cross_accuracy = cross_accuracy/num_runs

#identify the best top 10 accuracy-performing models
nseqid = np.unique(seqid)
acc_seqid = dict()
ind_seqid = dict()
sel_params = dict()
for i in range(1,len(nseqid)+1):
    pos = np.where(np.array(seqid) == i)
    pos = pos[0]
    tmp = []
    for j in pos:
        tmp.append(acc[j])
    acc_seqid[i] = tmp
    ind_seqid[i] = np.argpartition(tmp, -10)[-10:]
    sel_param = []
    for j in ind_seqid[i]:
        sel_param.append(best_params[pos[j]])
    sel_params[i] = sel_param


mean_acc = []
std_acc = []
best_acc = []
for key in sel_params.keys():
    tmp = sel_params[key]
    tmp_acc = []
    tmp_mem = []
    for t in range(len(tmp)):
        tmp_acc.append(tmp[t][4])
        tmp_mem.append(tmp[t][3])
    mean_acc.append(np.mean(tmp_acc))
    std_acc.append(np.std(tmp_acc))
    pos = np.where(tmp_mem == np.amax(tmp_mem))
    pos = pos[0]
    best_acc.append(tmp[pos[0]])

#export results
accuracy_sheet = []
nclus_sheet = []
memory_sheet = []
seqid_sheet = []
runid_sheet = []
for key in sel_params.keys():
    tmp = sel_params[key]
    tmp_acc = []
    tmp_mem = []
    tmp_nclus = []
    tmp_seqid = []
    tmp_runid = []
    for t in range(len(tmp)):
        tmp_acc.append(tmp[t][4])
        tmp_mem.append(tmp[t][3])
        tmp_nclus.append(tmp[t][2])
        tmp_seqid.append(tmp[t][1])
        tmp_runid.append(tmp[t][0])
    accuracy_sheet.append(tmp_acc)
    memory_sheet.append(tmp_mem)
    nclus_sheet.append(tmp_nclus)
    seqid_sheet.append(tmp_seqid)
    runid_sheet.append(tmp_runid)

flat_accuracy = [item for sublist in accuracy_sheet for item in sublist]
flat_nclus = [item for sublist in nclus_sheet for item in sublist]
flat_memory = [item for sublist in memory_sheet for item in sublist]
flat_seqid = [item for sublist in seqid_sheet for item in sublist]
flat_runid = [item for sublist in runid_sheet for item in sublist]
length_seqid = {'1':15, '2':10, '3':6, '4':5, '5':3, '6':2}
flat_length = [length_seqid[str(item)] for item in flat_seqid]

sss = {'Length': flat_length,
       'Runs': flat_runid,
       'Accuracy': flat_accuracy,
       'Memory': flat_memory,
       'ClusterID': flat_nclus
       }

df = pd.DataFrame(sss)
df.to_excel(excel_writer = "export_vomm_reuslts.xlsx")

best_1 = best_acc[1]
best_2 = best_acc[2]

best_seq_id_1 = best_1[1]
best_seq_id_2 = best_2[1]
best_nclus_1 = best_1[2]
best_nclus_2 = best_2[2]
best_mem_1 = best_1[3]
best_mem_2 = best_2[3]
train_seq_1 = best_1[8]
train_seq_2 = best_2[8]

model0_1 = best_1[5]
model0_2 = best_1[6]
model0_3 = best_1[7]
model1_1 = best_2[5]
model1_2 = best_2[6]
model1_3 = best_2[7]

#extract the two most probable sequence from each models
gen1_1 = model0_1.generate_data()
gen1_2 = model0_2.generate_data()
gen1_3 = model0_3.generate_data()

gen2_1 = model1_1.generate_data()
gen2_2 = model1_2.generate_data()
gen2_3 = model1_3.generate_data()

gen1_1 = [i for i in gen1_1 if i != 0]
gen1_2 = [i for i in gen1_2 if i != 0]
gen1_3 = [i for i in gen1_3 if i != 0]

gen2_1 = [i for i in gen2_1 if i != 0]
gen2_2 = [i for i in gen2_2 if i != 0]
gen2_3 = [i for i in gen2_3 if i != 0]

print('best1 for stim type 0 ' + str(gen1_1[:model0_1.d]))
print('best1 for stim type 1 ' + str(gen1_2[:model0_2.d]))
print('best1 for stim type 2 ' + str(gen1_3[:model0_3.d]))
print('best2 for stim type 0 ' + str(gen2_1[:model1_1.d]))
print('best2 for stim type 1 ' + str(gen2_2[:model1_2.d]))
print('best2 for stim type 2 ' + str(gen2_3[:model1_3.d]))
