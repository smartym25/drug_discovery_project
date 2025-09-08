#import libraries
import math
import numpy as np 
import pandas as pd
import seaborn as sns
from rdkit import Chem
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
from rdkit.Chem import Descriptors, Lipinski

#read the csv file that we prepared 
df = pd.read_csv("bioactivity_data_preprocessed.csv")

def lipinski(smiles):
    mol_data = []
    
    #add to the list mol_data the different mol translate from smiles
    for element in smiles:
        mol = Chem.MolFromSmiles(element)
        mol_data.append(mol)
        
    data = []    
          
    for mol in mol_data:
        desc_MolWt = Descriptors.MolWt(mol)
        desc_MolLogP = Descriptors.MolLogP(mol)
        desc_NumHDonors = Lipinski.NumHDonors(mol)
        desc_NumHAcceptors = Lipinski.NumHAcceptors(mol)
        
        rows = np.array([desc_MolWt, 
                        desc_MolLogP,
                        desc_NumHDonors,
                        desc_NumHAcceptors])
        
        data.append(rows)
        
    column_names = ["MolWT", "MolLogP", "num_Hdonors", "num_Hacceptors"]
    descriptors = pd.DataFrame(data, columns=column_names)
    
    return descriptors

df_lipinski = lipinski(df.canonical_smiles)

#combined the two collection-data and axis 1 allow to concat in the same rows
df_combined = pd.concat([df, df_lipinski], axis=1) 

df_combined = df_combined.rename(columns={"standard_value": "pIC50"})

df_combined.to_csv("bioactivity_data_descriptors.csv", index=False)

#convert IC50 in nanoMolar (standard_value) to plC50 in Molar
#valor enon deve essere superiore a 10^8 perchè se no darebbe come risultato un logaritmo negativo 
def pic50(dataframe):
    ic50 = dataframe.pIC50
    
    new_ic50 = []
    for value in ic50:
        nM_to_M = value * (10**-9)
        new_ic50.append(nM_to_M)
        
    pic50 = []
    for ic50 in new_ic50:
        value = -math.log10(ic50)
        pic50.append(value)

    return pic50 

df_combined["pIC50"] = pic50(df_combined)

# df_combined.pIC50.describe() useful to determinate the min and the max value in a column

#eliminate the rows with intermidiate activity status
df_2class = df_combined[df_combined.bioactivity_status != "intermidiate"]

#create the plot of bioactivity class with istograms
plt.figure(figsize=(5.5, 5.5))

sns.countplot(x="bioactivity_status", data=df_2class, edgecolor="black", hue="bioactivity_status")

plt.xlabel("Bioactivity class", fontsize=10, fontweight="bold")
plt.ylabel("Frequency", fontsize=10, fontweight="bold")

#create a scatter plot of MW and LogP
plt.figure(figsize=(5.5, 5.5))

sns.scatterplot(x="MolWT", y="MolLogP", data=df_2class, size="pIC50", 
                hue="bioactivity_status", edgecolor="black", alpha=0.8)

plt.legend(bbox_to_anchor=(1.05,1), loc=2, borderaxespad=0)

plt.xlabel("MW", fontsize=10, fontweight="bold")
plt.ylabel("LogP", fontsize=10, fontweight="bold")

#create the boxplot for each descriptor
def boxplot_mann_whitney(descriptor):
    plt.figure(figsize=(5.5, 5.5))
    
    sns.boxplot(x="bioactivity_status", y=descriptor, data=df_2class, hue="bioactivity_status")
    
    plt.xlabel("Bioactivity status", fontsize=10, fontweight="bold")
    plt.ylabel(descriptor, fontsize=10, fontweight="bold")

    active = df_2class[df_2class["bioactivity_status"]=="active"][descriptor]
    inactive = df_2class[df_2class["bioactivity_status"]=="inactive"][descriptor]
    
    u, p = mannwhitneyu(active, inactive, alternative="two-sided")

    print(p)

#TODO create a table with p value and the consequence of the value
#pIC50
boxplot_mann_whitney("pIC50")
#MolWT
boxplot_mann_whitney("MolWT")
#MolLogP
boxplot_mann_whitney("MolLogP")
#num_Hdonors
boxplot_mann_whitney("num_Hdonors")
#num_Hacceptors
boxplot_mann_whitney("num_Hacceptors")

df_2class.to_csv("molecule.smi", index=False)

df_2class