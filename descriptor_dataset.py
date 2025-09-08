import pandas as pd 
from padelpy import padeldescriptor

df = pd.read_csv("bioactivity_data_descriptors.csv")

selection = ["molecule_chembl_id", "canonical_smiles"]
df_selection = df[selection]

#convert the csv file into smi file because PaDEL can read only the SMI format
with open("molecule.smi", "w") as f: 

    for _, row in df.iterrows(): #scorre riga per riga

        f.write(f"{row['canonical_smiles']} {row['molecule_chembl_id']}\n")
        #prende nella stessa riga l'elemento della colonna canonical_smiles e molecule_id 
        #li scrive assieme nello stesso rigo seprati da uno spazio 

fingerprint = "descriptors_output"

fingerprint_output_file = ''.join([fingerprint,'.csv'])

fingerprint_directory = "PaDEL-Descriptor/descriptors.xml"

padeldescriptor(
    mol_dir="molecule.smi",                  
    d_file=fingerprint_output_file,          #output
    descriptortypes=fingerprint_directory,   #directory
    detectaromaticity=True,                  #attiva l'aromaticità dei composti
    standardizenitro=True,                   #standardizza i gruppi nitro in modo che vengano tratti in modo uguale
    standardizetautomers=True,               #standardizza i tautomeri (forme diverse della stessa molecola)
    removesalt=True,                         #rimuove i sali 
    log=False,                               #tempo impiegato per fare i singoli calcoli di ogni composto (eliminato)        
    threads=2,                               #aumenta la velocità di calcolo influenzando la CPU
    fingerprints=True                        #calcola la rappresentazione binaria della struttura chimica
)

descriptors = pd.read_csv(fingerprint_output_file)

df3_X = descriptors.drop(columns=["Name"])

df3_Y = df["pIC50"]

df3 = pd.concat([df3_X, df3_Y], axis=1)

df3.to_csv("alzheimer_bioactivity_2class_descriptors_pIC50.csv", index=False)