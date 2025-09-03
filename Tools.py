import numpy as np
import tmap as tm
import pandas as pd
from mhfp.encoder import MHFPEncoder
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem

def smi_descriptors(mol):

    descriptors = {
        'MolWt': Descriptors.MolWt,
        'TPSA': Descriptors.TPSA,
        'LogP': Descriptors.MolLogP,
        'NumHAcceptors': Descriptors.NumHAcceptors,
        'NumHDonors': Descriptors.NumHDonors
    }

    results = {desc: func(mol) for desc, func in descriptors.items()}

    return results.values()

def calculate_descriptors(df, smi_col):
    
    enc = MHFPEncoder(1024)

    fps = []
    molwt_list, tpsa_list, logp_list, ha_list, hd_list = [], [], [], [], []
    for i, row in df.iterrows():
        try:
            if i != 0 and i % 1000 == 0:
                print(f'{round(100 * i / len(df), 2)}% of the dataset processed to fps')
            
            mol = AllChem.MolFromSmiles(row[smi_col])
            fps.append(tm.VectorUint(enc.encode_mol(mol)))
            
            molwt, tpsa, logp, ha, hd = smi_descriptors(mol)
            molwt_list.append(molwt)
            tpsa_list.append(tpsa)
            logp_list.append(logp)
            ha_list.append(ha)
            hd_list.append(hd)
            
        except:
            print(f'{row[smi_col]} arised an issue and could not be included in the analysis')
    
    return fps, molwt_list, tpsa_list, logp_list, ha_list, hd_list

def fps_mhfp_gen(smi):
    
    enc = MHFPEncoder(1024)
    try:
        mol = Chem.MolFromSmiles(smi)
        fp = tm.VectorUint(enc.encode_mol(mol))
        return fp
    except:
        return np.NaN
        
