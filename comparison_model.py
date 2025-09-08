import pandas as pd 
import seaborn as sns
from sklearn.model_selection import train_test_split
from lazypredict.Supervised import LazyRegressor
import matplotlib.pyplot as plt

df = pd.read_csv("csv-files/papillomavirus_bioactivity_2class_descriptors_pIC50.csv")

#PARTE UGUALE IN TUTTI I MODELLI DI MACHINE LEARNING 
#In order to create the matrix of single molecule
X = df.drop("pIC50", axis=1)

Y = df.pIC50

#remove low variance features, pre-processing step to eliminate feature that don't contribute much to the model
from sklearn.feature_selection import VarianceThreshold

#initialize the VarianceThreshold filter 
#by default it removes features with variance = 0 (the same value in all rows)

selection = VarianceThreshold(threshold=(.8 * (1 - .8)))
#we choose (.8*(1-.8)) because It rapresents the variance of Bernoulli distribution, used with two element (binary)
X = selection.fit_transform(X)

#data split bewteen the percentage for training and test, in this case 80/20 ratio 80% training and 20% test
X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.2)

#PARTE DEDICATA AL CONFRONTO CON 20 MODELLI DI MACHINE LEARNING 
clf = LazyRegressor(verbose=0, ignore_warnings=True, custom_metric=None)
train,test = clf.fit(X_train, X_test, Y_train, Y_test)

# test: print on the interactive window the table for the comparison between models in the test phase
# train: print on the interactive window the table for the comparison between models in the training phase

#compared the models in based of R2 in the test part
plt.figure(figsize=(5,10))
sns.set_theme(style="whitegrid")
ax = sns.barplot(y=test.index, x="R-Squared", data=test)

ax.set(xlim=(0, 1))

#compared the models in based of R2 in the training part
plt.figure(figsize=(5,10))
sns.set_theme(style="whitegrid")
ax = sns.barplot(y=train.index, x="R-Squared", data=train)

ax.set(xlim=(0, 1))



