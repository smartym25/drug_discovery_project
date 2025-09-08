import pandas as pd 
import lazypredict
import numpy as np
import seaborn as sns
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
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

#PARTE DEDICATA AL MODELLO DEL RANDOM-FOREST-REGRESSION
np.random.seed(100)
model = RandomForestRegressor(n_estimators=100)
#n_estimators: costruisce 100 alberi decisionali indipendenti di cui viene poi svolta la media tra i risulati ottenuti 
#random_state: utilizza una fonte di casualità fissa così quando rilancio il modello i valori sono uguali 

model.fit(X_train, Y_train)
#training the model

r2 = model.score(X_test, Y_test)
#r2 percentuale di varianza del target, più il valore è alto e più la previsione è vicina alla realtà

Y_pred = model.predict(X_test)

sns.set(color_codes=True)
sns.set_style("whitegrid")

ax = sns.regplot(x=Y_test, y=Y_pred, scatter_kws={"alpha":0.4})

ax.set_xlabel("experimental pIC50", fontsize="large", fontweight="bold")
ax.set_ylabel("predicted pIC50", fontsize="large", fontweight="bold")

ax.figure.set_size_inches(5,5)

#papillomavirus 0.4 100 100
#alzheimer 0.35 300 300