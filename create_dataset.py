#import libraries 
import pandas as pd 
from chembl_webresource_client.new_client import new_client

#target search for coronavius
target = new_client.target
target_query = target.search("alzheimer")
targets = pd.DataFrame.from_dict(target_query)

#select and retrieve the bioactivity data for SPECIFIC CORONAVIRUS
#The single protein corresponds to the fifth entry
selected_target = targets.target_chembl_id[0]

#retrieve the bioactivity data only for selected_target that are reported in IC50 
#value that indicate how a substance is effective to stop a biological process 

activity = new_client.activity
res = activity.filter(target_chembl_id=selected_target).filter(standard_type="IC50")

df = pd.DataFrame.from_dict(res)

#for missing data in other dataset use this code
df2 = df[df.standard_value.notna()]

#create the csv file
df2.to_csv("bioactivity_data.csv", index=False)

#data pre-processing of the bioactiivity data 
bioactivity_class = []

for i in df2.standard_value:
    if float(i) >= 10000:
        bioactivity_class.append("inactive")
    elif float(i) <= 1000:
        bioactivity_class.append("active")
    else:
        bioactivity_class.append("intermidiate")
        
df3 = pd.DataFrame(columns=["bioactivity_status", "molecule_chembl_id", "canonical_smiles", "standard_value"])  

chembl_id = []
for i in df2.molecule_chembl_id:
    chembl_id.append(i)
    
df3["molecule_chembl_id"] = chembl_id
    
canonical_smiles = []
for i in df2.canonical_smiles:
    canonical_smiles.append(i)
    
df3["canonical_smiles"] = canonical_smiles
    
standard_value = []
for i in df2.standard_value:
    standard_value.append(i)
    
df3["standard_value"] = standard_value

df3["bioactivity_status"] = bioactivity_class

df3.to_csv("bioactivity_data_preprocessed.csv", index=False)

df3