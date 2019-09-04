import pandas as pd
from sklearn import preprocessing, model_selection, linear_model
from matplotlib import pyplot as plt

file = 'outData/CSV/attrNets.csv'
df = pd.read_csv(file)
df = df.drop(columns=['atlas','modality'])

df[['Study.ID', 'Group']] = df[['Study.ID', 'Group']].astype('category')

selT = 0.30
selFeat = ['assortativity', 'Cp', 'density', 'diameter', 'diameter.wt', 
'E.global', 'E.global.wt', 'E.local', 'E.local.wt', 'Lp', 
'mod', 'mod.wt', 'transitivity', 'strength', 'vulnerability']

X = df.loc[df['threshold'] == selT, selFeat]
y = df.loc[df['threshold'] == selT, ['Group']]

## Preprocessing 

# Standardization

# Mean based 
X = preprocessing.scale(X)

# Median based
X = X.apply(lambda x : preprocessing.robust_scale(x))

# Train/test data
X_train, X_test, y_train, y_test = model_selection.train_test_split(
    X, y, test_size = 0.35, random_state = 666)

## Model 
clsf = linear_model.SGDClassifier(random_state = 666)

# GridSearchCV
parameters = {'alpha' : [0.00001, 0.00005, 0.0001, 0.0005, 0.001, 0.005, 0.01, 
0.05, 0.1, 0.5, 1],
'loss' : ['hinge', 'log'], 
'penalty' : ['l1', 'l2']}
searcher = model_selection.GridSearchCV(clsf, parameters, cv = 20)
searcher.fit(X_train, y_train.values.ravel())

# Results
print("Best CV params", searcher.best_params_)
print("Best CV accuracy", searcher.best_score_)
print("Test accuracy of best grid search hypers:", searcher.score(X_test, y_test))
