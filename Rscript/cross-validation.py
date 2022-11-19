import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split, StratifiedKFold
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import accuracy_score, recall_score,roc_auc_score, f1_score, cohen_kappa_score, precision_score, roc_curve
from sklearn.ensemble import RandomForestClassifier

# 5-fold Cross-validation on ITS1, ITS2 dataset separately
data=pd.read_csv("../ITS12_genus_data.csv",index_col=0)
ITS12_inf = ps.read_csv("../ITS12_inf.csv",index_col=0)
data['dataset'] = list(ITS12_inf.loc[data.index,'Assay.Type'])
data['enterotype'] = list(ITS12_inf.loc[data.index,'enterotype_cluster'])

driver = ["k__Fungi.p__Ascomycota.c__Eurotiomycetes.o__Eurotiales.f__Aspergillaceae.g__Aspergillus",
          "k__Fungi.p__Ascomycota.c__Saccharomycetes.o__Saccharomycetales.f__Incertae_sedis.g__Candida",
          "k__Fungi.p__Ascomycota.c__Saccharomycetes.o__Saccharomycetales.f__Saccharomycetaceae.g__Saccharomyces",
          "k__Fungi.p__Ascomycota.c__Saccharomycetes.o__Saccharomycetales.f__unidentified.g__unidentified",
          "k__Fungi.p__Ascomycota.c__unidentified.o__unidentified.f__unidentified.g__unidentified"]
train_data = data.loc[(data['Dataset']=='ITS2'),].iloc[:,:(data.shape[1]-2)]
train_label = np.array(data.loc[data['Dataset']=='ITS2','enterotype'])
selected_taxa = set(train_data.columns).difference(driver)
train_data = train_data.loc[:,selected_taxa]

train_data = data.iloc[:,:(data.shape[1]-2)]
train_label = np.array(data['enterotype'])

# Without driver
selected_taxa = set(train_data.columns).difference(driver)
train_data = train_data.loc[:,selected_taxa]


skf = StratifiedKFold(n_splits=5)
predict_pro = np.zeros([train_data.shape[0],4])
for i in range(10):
    possibility = np.zeros([train_data.shape[0], 4])
    for train_index, test_index in skf.split(train_data, train_label):
        cv_train_data = train_data.iloc[train_index,]
        cv_train_label = train_label[train_index]
        cv_test_data = train_data.iloc[test_index,]
        cv_test_label = train_label[test_index]
        e1 = cv_train_data.loc[cv_train_label=="fun_S_E",]
        e2 = cv_train_data.loc[cv_train_label=="fun_C_E",]
        e3 = cv_train_data.loc[cv_train_label=="fun_A_E",]
        e4 = pd.concat([cv_train_data.loc[cv_train_label=="fun_AS_E",],cv_train_data.loc[cv_train_label=="Ap_E",]])
        size = min(e1.shape[0],e2.shape[0],e3.shape[0],e4.shape[0])
        e1=e1.iloc[random.sample(range(e1.shape[0]),size),]
        e2=e2.iloc[random.sample(range(e2.shape[0]),size),]
        e3=e3.iloc[random.sample(range(e3.shape[0]),size),]
        e4=e4.iloc[random.sample(range(e4.shape[0]),size),]
        cv_train_balance_data = pd.concat([e1,e2,e3,e4])
        cv_train_label = np.array(data.loc[cv_train_balance_data.index,'enterotype'])
        cv_train_label[cv_train_label=='fun_S_E']  = 0
        cv_train_label[cv_train_label=='fun_C_E'] = 1
        cv_train_label[cv_train_label=='fun_A_E'] = 2
        cv_train_label[cv_train_label=='fun_AS_E'] = 3
        cv_test_label[cv_test_label=='fun_S_E']  = 0
        cv_test_label[cv_test_label=='fun_C_E'] = 1
        cv_test_label[cv_test_label=='fun_A_E'] = 2
        cv_test_label[cv_test_label=='fun_AS_E'] = 3
        clf = LogisticRegression(penalty='l1', solver='liblinear')
        clf.fit(np.array(cv_train_balance_data),list(cv_train_label))
        acc = clf.score(cv_test_data,list(cv_test_label))
        pre = clf.predict(cv_test_data)
        score = clf.predict_proba(cv_test_data)
        possibility[test_index] = score
    labels = np.zeros(len(train_label))
    labels[train_label=='fun_S_E']=1
    print("fun_S_E AUC:")
    roc_auc_score(labels,possibility[:,0],average="weighted")
    labels = np.zeros(len(train_label))
    labels[train_label=='fun_C_E']=1
    print("fun_C_E AUC:")
    roc_auc_score(labels,possibility[:,1],average="weighted")
    labels = np.zeros(len(train_label))
    labels[train_label=='fun_A_E']=1
    print("fun_A_E AUC:")
    roc_auc_score(labels,possibility[:,2],average="weighted")
    labels = np.zeros(len(train_label))
    labels[train_label=='fun_AS_E']=1
    print("fun_AS_E AUC:")
    roc_auc_score(labels,possibility[:,3],average="weighted")
    predict_pro += possibility

predict_pro = predict_pro/10
predict_pro = pd.DataFrame(predict_pro)
predict_pro['label'] = train_label
predict_pro.columns = ['fun_S_E','fun_C_E','fun_A_E','fun_AS_E','enterotype']

#predict_pro.to_csv("../ITS1_5CV_without_driver.csv")
#predict_pro.to_csv("../ITS1_5CV_with_driver.csv")
#predict_pro.to_csv("../ITS2_5CV_without_driver.csv")
#predict_pro.to_csv("../ITS2_5CV_with_driver.csv")
#predict_pro.to_csv("../ITS12_5CV_with_driver.csv")
#predict_pro.to_csv("../ITS12_5CV_without_driver.csv")

from scipy import interp
score_data = pd.read_csv("../ITS12_5CV_with_driver.csv",index_col=0)
score_data.loc[score_data['enterotype']=='fun_S_E','label']  = 0
score_data.loc[score_data['enterotype']=='fun_C_E','label']  = 1
score_data.loc[score_data['enterotype']=='fun_A_E','label']  = 2
score_data.loc[score_data['enterotype']=='fun_AS_E','label']  = 3


import numpy as np
import matplotlib.pyplot as plt
from itertools import cycle
from sklearn.metrics import roc_curve, auc

y_score = np.array(score_data.iloc[:,:4])
y_label = np.array(score_data['label'])
y_label = np.zeros([score_data.shape[0],4])
y_label[score_data['enterotype']=='fun_S_E',0] = 1
y_label[score_data['enterotype']=='fun_C_E',1] = 1
y_label[score_data['enterotype']=='fun_A_E',2] = 1
y_label[score_data['enterotype']=='fun_AS_E',3] = 1



fpr_with = dict()
tpr_with = dict()
roc_auc_with = dict()
for i in range(4):
    fpr_with[i], tpr_with[i], _ = roc_curve(y_label[:, i], y_score[:, i])
    roc_auc_with[i] = auc(fpr_with[i], tpr_with[i])

fpr_with["micro"], tpr_with["micro"], _ = roc_curve(y_label.ravel(), y_score.ravel())
roc_auc_with["micro"] = auc(fpr_with["micro"], tpr_with["micro"])


score_data = pd.read_csv("../ITS12_5CV_without_driver.csv",index_col=0)
score_data.loc[score_data['enterotype']=='fun_S_E','label']  = 0
score_data.loc[score_data['enterotype']=='fun_C_E','label']  = 1
score_data.loc[score_data['enterotype']=='fun_A_E','label']  = 2
score_data.loc[score_data['enterotype']=='fun_AS_E','label']  = 3

import numpy as np
import matplotlib.pyplot as plt
from itertools import cycle
from sklearn.metrics import roc_curve, auc

y_score = np.array(score_data.iloc[:,:4])
y_label = np.array(score_data['label'])
y_label = np.zeros([score_data.shape[0],4])
y_label[score_data['enterotype']=='fun_S_E',0] = 1
y_label[score_data['enterotype']=='fun_C_E',1] = 1
y_label[score_data['enterotype']=='fun_A_E',2] = 1
y_label[score_data['enterotype']=='fun_AS_E',3] = 1



fpr_without = dict()
tpr_without = dict()
roc_auc_without = dict()
for i in range(4):
    fpr_without[i], tpr_without[i], _ = roc_curve(y_label[:, i], y_score[:, i])
    roc_auc_without[i] = auc(fpr_without[i], tpr_without[i])

fpr_without["micro"], tpr_without["micro"], _ = roc_curve(y_label.ravel(), y_score.ravel())
roc_auc_without["micro"] = auc(fpr_without["micro"], tpr_without["micro"])


lw = 2
plt.figure(figsize=(8,5))
plt.plot(fpr_with[0],tpr_with[0],color='#C1CD24',lw=lw,label='fun_S_E (with driver: %0.2f)' % roc_auc_with[0])
plt.plot(fpr_with[1],tpr_with[1],color='#D23837',lw=lw,label='fun_C_E (with driver: %0.2f)' % roc_auc_with[1])
plt.plot(fpr_with[2],tpr_with[2],color='#439FC2',lw=lw,label='fun_A_E (with driver: %0.2f)' % roc_auc_with[2])
plt.plot(fpr_with[3],tpr_with[3],color='#5E52A0',lw=lw,label='fun_AS_E (with driver: %0.2f)' % roc_auc_with[3])
plt.plot(fpr_with["micro"],tpr_with["micro"],color='#F6C564',lw=lw,label='Average (with driver: %0.2f)' % roc_auc_with["micro"])

plt.plot(fpr_without[0],tpr_without[0],color='#C1CD24',linestyle=':',lw=lw,label='fun_S_E (without driver: %0.2f)' % roc_auc_without[0])
plt.plot(fpr_without[1],tpr_without[1],color='#D23837',linestyle=':',lw=lw,label='fun_C_E (without driver: %0.2f)' % roc_auc_without[1])
plt.plot(fpr_without[2],tpr_without[2],color='#439FC2',linestyle=':',lw=lw,label='fun_A_E (without driver: %0.2f)' % roc_auc_without[2])
plt.plot(fpr_without[3],tpr_without[3],color='#5E52A0',linestyle=':',lw=lw,label='fun_AS_E (without driver: %0.2f)' % roc_auc_without[3])
plt.plot(fpr_without["micro"],tpr_without["micro"],color='#F6C564',linestyle=':',lw=lw,label='Average (without driver: %0.2f)' % roc_auc_without["micro"])


plt.xlim([0,1])
plt.ylim([0,1])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('ROC curve(ITS1 & ITS2)')
plt.legend(loc='lower right')
#plt.show()
plt.savefig('../ITS12_5CV.pdf')
plt.close()




# Cross-dataset comparison
data=pd.read_csv("../ITS12_genus_data.csv",index_col=0)
ITS12_inf = ps.read_csv("../ITS12_inf.csv",index_col=0)
data['dataset'] = list(ITS12_inf.loc[data.index,'Assay.Type'])
data['enterotype'] = list(ITS12_inf.loc[data.index,'enterotype_cluster'])

driver = ["k__Fungi.p__Ascomycota.c__Eurotiomycetes.o__Eurotiales.f__Aspergillaceae.g__Aspergillus",
          "k__Fungi.p__Ascomycota.c__Saccharomycetes.o__Saccharomycetales.f__Incertae_sedis.g__Candida",
          "k__Fungi.p__Ascomycota.c__Saccharomycetes.o__Saccharomycetales.f__Saccharomycetaceae.g__Saccharomyces",
          "k__Fungi.p__Ascomycota.c__Saccharomycetes.o__Saccharomycetales.f__unidentified.g__unidentified",
          "k__Fungi.p__Ascomycota.c__unidentified.o__unidentified.f__unidentified.g__unidentified"]
train_data = data.loc[(data['Dataset']=='ITS1'),].iloc[:,:(data.shape[1]-2)]
test_data = data.loc[(data['Dataset']=='ITS2'),].iloc[:,:(data.shape[1]-2)]
test_data = test_data.loc[:,test_data.sum()>0]
train_data = train_data.loc[:,train_data.sum()>0]
selected_taxa = set(train_data.columns).difference(driver)
train_data = train_data.loc[:,selected_taxa]
test_data = test_data.loc[:,selected_taxa]
train_label = np.array(data.loc[data['Dataset']=='ITS1','enterotype'])
test_label = np.array(data.loc[data['Dataset']=='ITS2','enterotype'])

import random
score_result = np.zeros([test_data.shape[0],4])
for i in range(10):
    e1 = train_data.loc[train_label=="fun_S_E",]
    e2 = train_data.loc[train_label=="fun_C_E",]
    e3 = train_data.loc[train_label=="fun_A_E",]
    e4 = train_data.loc[train_label=="fun_AS_E",]
    size = min(e1.shape[0],e2.shape[0],e3.shape[0],e4.shape[0])
    e1=e1.iloc[random.sample(range(e1.shape[0]),size),]
    e2=e2.iloc[random.sample(range(e2.shape[0]),size),]
    e3=e3.iloc[random.sample(range(e3.shape[0]),size),]
    e4=e4.iloc[random.sample(range(e4.shape[0]),size),]
    train_balance_data = pd.concat([e1,e2,e3,e4])
    label = np.array(data.loc[train_balance_data.index,'enterotype'])
    label[np.array(data.loc[train_balance_data.index,'enterotype'])=='fun_S_E']=0
    label[np.array(data.loc[train_balance_data.index,'enterotype'])=='fun_C_E']=1
    label[np.array(data.loc[train_balance_data.index,'enterotype'])=='fun_A_E']=2
    label[np.array(data.loc[train_balance_data.index,'enterotype'])=='fun_AS_E']=3
    clf = LogisticRegression(penalty='l1', solver='liblinear')
    clf.fit(np.array(train_balance_data),list(label))
    test_label = np.array(data.loc[test_data.index,'enterotype'])
    test_label[np.array(data.loc[test_data.index,'enterotype'])=='fun_S_E'] = 0
    test_label[np.array(data.loc[test_data.index,'enterotype'])=='fun_C_E'] = 1
    test_label[np.array(data.loc[test_data.index,'enterotype'])=='fun_A_E'] = 2
    test_label[np.array(data.loc[test_data.index,'enterotype'])=='fun_AS_E'] = 3
    acc = clf.score(test_data,list(test_label))
    pre = clf.predict(test_data)
    score = clf.predict_proba(test_data)
    score_result += score
    labels = np.zeros(len(test_label))
    labels[test_label==0]=1
    print("S_E AUC:")
    roc_auc_score(labels,score[:,0],average="weighted")
    labels = np.zeros(len(test_label))
    labels[test_label==1]=1
    print("C_E AUC:")
    roc_auc_score(labels,score[:,1],average="weighted")
    labels = np.zeros(len(test_label))
    labels[test_label==2]=1
    print("A_E AUC:")
    roc_auc_score(labels,score[:,2],average="weighted")
    labels = np.zeros(len(test_label))
    labels[test_label==3]=1
    print("Ap_E AUC:")
    roc_auc_score(labels,score[:,3],average="weighted")


score_data = pd.DataFrame(score_result)/10
score_data.index = test_data.index
score_data['label'] = test_label
score_data['enterotype'] = data.loc[score_data.index,['enterotype']]
score_data.columns = ['fun_S_E','fun_C_E','fun_A_E','fun_AS_E','label','enterotype']
score_data.to_csv("../ITS1_transfer_ITS2_ROC_curve_without_driver.csv")


data=pd.read_csv("../ITS12_genus_data.csv",index_col=0)
ITS12_inf = ps.read_csv("../ITS12_inf.csv",index_col=0)
data['dataset'] = list(ITS12_inf.loc[data.index,'Assay.Type'])
data['enterotype'] = list(ITS12_inf.loc[data.index,'enterotype_cluster'])
driver = ["k__Fungi.p__Ascomycota.c__Eurotiomycetes.o__Eurotiales.f__Aspergillaceae.g__Aspergillus",
          "k__Fungi.p__Ascomycota.c__Saccharomycetes.o__Saccharomycetales.f__Incertae_sedis.g__Candida",
          "k__Fungi.p__Ascomycota.c__Saccharomycetes.o__Saccharomycetales.f__Saccharomycetaceae.g__Saccharomyces",
          "k__Fungi.p__Ascomycota.c__Saccharomycetes.o__Saccharomycetales.f__unidentified.g__unidentified",
          "k__Fungi.p__Ascomycota.c__unidentified.o__unidentified.f__unidentified.g__unidentified"]
train_data = data.loc[(data['Dataset']=='ITS2'),].iloc[:,:(data.shape[1]-2)]
test_data = data.loc[(data['Dataset']=='ITS1'),].iloc[:,:(data.shape[1]-2)]
test_data = test_data.loc[:,test_data.sum()>0]
train_data = train_data.loc[:,train_data.sum()>0]
selected_taxa = set(train_data.columns).difference(driver)
train_data = train_data.loc[:,selected_taxa]
test_data = test_data.loc[:,selected_taxa]
train_label = np.array(data.loc[data['Dataset']=='ITS2','enterotype'])
test_label = np.array(data.loc[data['Dataset']=='ITS1','enterotype'])

import random
score_result = np.zeros([test_data.shape[0],4])
for i in range(10):
    e1 = train_data.loc[train_label=="fun_S_E",]
    e2 = train_data.loc[train_label=="fun_C_E",]
    e3 = train_data.loc[train_label=="fun_A_E",]
    e4 = train_data.loc[train_label=="fun_AS_E",]
    size = min(e1.shape[0],e2.shape[0],e3.shape[0],e4.shape[0])
    e1=e1.iloc[random.sample(range(e1.shape[0]),size),]
    e2=e2.iloc[random.sample(range(e2.shape[0]),size),]
    e3=e3.iloc[random.sample(range(e3.shape[0]),size),]
    e4=e4.iloc[random.sample(range(e4.shape[0]),size),]
    train_balance_data = pd.concat([e1,e2,e3,e4])
    label = np.array(data.loc[train_balance_data.index,'enterotype'])
    label[np.array(data.loc[train_balance_data.index,'enterotype'])=='fun_S_E']=0
    label[np.array(data.loc[train_balance_data.index,'enterotype'])=='fun_C_E']=1
    label[np.array(data.loc[train_balance_data.index,'enterotype'])=='fun_A_E']=2
    label[np.array(data.loc[train_balance_data.index,'enterotype'])=='fun_AS_E']=3
    clf = LogisticRegression(penalty='l1', solver='liblinear')
    clf.fit(np.array(train_balance_data),list(label))
    test_label = np.array(data.loc[test_data.index,'enterotype'])
    test_label[np.array(data.loc[test_data.index,'enterotype'])=='fun_S_E'] = 0
    test_label[np.array(data.loc[test_data.index,'enterotype'])=='fun_C_E'] = 1
    test_label[np.array(data.loc[test_data.index,'enterotype'])=='fun_A_E'] = 2
    test_label[np.array(data.loc[test_data.index,'enterotype'])=='fun_AS_E'] = 3
    acc = clf.score(test_data,list(test_label))
    pre = clf.predict(test_data)
    score = clf.predict_proba(test_data)
    score_result += score
    labels = np.zeros(len(test_label))
    labels[test_label==0]=1
    print("fun_S_E AUC:")
    roc_auc_score(labels,score[:,0],average="weighted")
    labels = np.zeros(len(test_label))
    labels[test_label==1]=1
    print("fun_C_E AUC:")
    roc_auc_score(labels,score[:,1],average="weighted")
    labels = np.zeros(len(test_label))
    labels[test_label==2]=1
    print("fun_A_E AUC:")
    roc_auc_score(labels,score[:,2],average="weighted")
    labels[test_label==3]=1
    print("fun_AS_E AUC:")
    roc_auc_score(labels,score[:,3],average="weighted")


score_data = pd.DataFrame(score_result)/10
score_data.index = test_data.index
score_data['label'] = test_label
score_data['enterotype'] = data.loc[score_data.index,['enterotype']]
score_data.columns = ['fun_S_E','fun_C_E','fun_A_E','fun_AS_E','label','enterotype']
score_data.to_csv("../ITS2_transfer_ITS1_ROC_curve_without_driver.csv")



from scipy import interp
score_data = pd.read_csv("../ITS1_transfer_ITS2_ROC_curve_with_driver.csv",index_col=0)
score_data.loc[score_data['enterotype']=='fun_S_E','label']  = 0
score_data.loc[score_data['enterotype']=='fun_C_E','label']  = 1
score_data.loc[score_data['enterotype']=='fun_A_E','label']  = 2
score_data.loc[score_data['enterotype']=='fun_AS_E','label']  = 3

import numpy as np
import matplotlib.pyplot as plt
from itertools import cycle
from sklearn.metrics import roc_curve, auc

y_score = np.array(score_data.iloc[:,:4])
y_label = np.array(score_data['label'])
y_label = np.zeros([score_data.shape[0],4])
y_label[score_data['enterotype']=='fun_S_E',0] = 1
y_label[score_data['enterotype']=='fun_C_E',1] = 1
y_label[score_data['enterotype']=='fun_A_E',2] = 1
y_label[score_data['enterotype']=='fun_AS_E',3] = 1


fpr_with = dict()
tpr_with = dict()
roc_auc_with = dict()
for i in range(4):
    fpr_with[i], tpr_with[i], _ = roc_curve(y_label[:, i], y_score[:, i])
    roc_auc_with[i] = auc(fpr_with[i], tpr_with[i])

fpr_with["micro"], tpr_with["micro"], _ = roc_curve(y_label.ravel(), y_score.ravel())
roc_auc_with["micro"] = auc(fpr_with["micro"], tpr_with["micro"])


score_data = pd.read_csv("../ITS1_transfer_ITS2_ROC_curve_without_driver.csv",index_col=0)
score_data.loc[score_data['enterotype']=='fun_S_E','label']  = 0
score_data.loc[score_data['enterotype']=='fun_C_E','label']  = 1
score_data.loc[score_data['enterotype']=='fun_A_E','label']  = 2
score_data.loc[score_data['enterotype']=='fun_AS_E','label']  = 3

import numpy as np
import matplotlib.pyplot as plt
from itertools import cycle
from sklearn.metrics import roc_curve, auc

y_score = np.array(score_data.iloc[:,:4])
y_label = np.array(score_data['label'])
y_label = np.zeros([score_data.shape[0],4])
y_label[score_data['enterotype']=='fun_S_E',0] = 1
y_label[score_data['enterotype']=='fun_C_E',1] = 1
y_label[score_data['enterotype']=='fun_A_E',2] = 1
y_label[score_data['enterotype']=='fun_AS_E',3] = 1


fpr_without = dict()
tpr_without = dict()
roc_auc_without = dict()
for i in range(4):
    fpr_without[i], tpr_without[i], _ = roc_curve(y_label[:, i], y_score[:, i])
    roc_auc_without[i] = auc(fpr_without[i], tpr_without[i])

fpr_without["micro"], tpr_without["micro"], _ = roc_curve(y_label.ravel(), y_score.ravel())
roc_auc_without["micro"] = auc(fpr_without["micro"], tpr_without["micro"])


lw = 2
plt.figure(figsize=(8,5))
plt.plot(fpr_with[0],tpr_with[0],color='#C1CD24',lw=lw,label='fun_S_E (with driver: %0.2f)' % roc_auc_with[0])
plt.plot(fpr_with[1],tpr_with[1],color='#D23837',lw=lw,label='fun_C_E (with driver: %0.2f)' % roc_auc_with[1])
plt.plot(fpr_with[2],tpr_with[2],color='#439FC2',lw=lw,label='fun_A_E (with driver: %0.2f)' % roc_auc_with[2])
plt.plot(fpr_with[3],tpr_with[3],color='#5E52A0',lw=lw,label='fun_AS_E (with driver: %0.2f)' % roc_auc_with[3])
plt.plot(fpr_with["micro"],tpr_with["micro"],color='#F6C564',lw=lw,label='Average (with driver: %0.2f)' % roc_auc_with["micro"])
plt.plot(fpr_without[0],tpr_without[0],color='#C1CD24',linestyle=':',lw=lw,label='fun_S_E (without driver: %0.2f)' % roc_auc_without[0])
plt.plot(fpr_without[1],tpr_without[1],color='#D23837',linestyle=':',lw=lw,label='fun_C_E (without driver: %0.2f)' % roc_auc_without[1])
plt.plot(fpr_without[2],tpr_without[2],color='#439FC2',linestyle=':',lw=lw,label='fun_A_E (without driver: %0.2f)' % roc_auc_without[2])
plt.plot(fpr_without[3],tpr_without[3],color='#5E52A0',linestyle=':',lw=lw,label='fun_AS_E (without driver: %0.2f)' % roc_auc_without[3])
plt.plot(fpr_without["micro"],tpr_without["micro"],color='#F6C564',linestyle=':',lw=lw,label='Average (without driver: %0.2f)' % roc_auc_without["micro"])


plt.xlim([0,1])
plt.ylim([0,1])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('ROC curve(ITS1 transfer ITS2)')
plt.legend(loc='lower right')
#plt.show()
plt.savefig('../ITS1-transfer-ITS2.pdf')
plt.close()
