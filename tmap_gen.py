# -*- coding: utf-8 -*-
"""
Created on 2025

@author: andres.sanchez
"""

import pickle
import tmap as tm
import pandas as pd
from rdkit.Chem import AllChem
from mhfp.encoder import MHFPEncoder
from faerun import Faerun
from matplotlib import pyplot as plt
import os
import argparse
from collections import Counter
    
def calculate_fps(df, smi_col):
    enc = MHFPEncoder(1024)
    fps = []
    for i, row in df.iterrows():
        try:
            if i != 0 and i % 100 == 0:
                print(f'{round(100 * i / len(df), 2)}% of the dataset processed to fps')
            
            mol = AllChem.MolFromSmiles(row[smi_col])
            fps.append(tm.VectorUint(enc.encode_mol(mol)))
            
        except:
            print(f'{row[smi_col]} arised an issue and could not be included in the analysis')
    return fps, []

def collapse_to_top_k(labels, data):
    counts = Counter(data)
    top_k = [cat for cat, _ in counts.most_common(max_cat)] # modify in flag --max_cat

    # Create mapping: top k keep original, others -> "Other"
    new_labels = []
    mapping = {}
    new_idx = 0

    for cat, name in labels:
        if cat in top_k:
            mapping[cat] = new_idx
            new_labels.append((new_idx, name))
            new_idx += 1

    # Add "Other" as last category
    other_idx = new_idx
    new_labels.append((other_idx, "Other"))

    new_data = [mapping.get(cat, other_idx) for cat in data]
    return new_data, new_labels

def main():
    df = pd.read_csv(data_path, sep=None, engine="python")

    lf = tm.LSHForest(1024, 64)
    
    print('Let\'s start with fingerprint generation.')
    print('This might take a while...')
    if not cont_col:
        descriptors = calculate_fps(df, smi_colum_name)
        fps = descriptors[0]
        
    elif cont_col:
        fps = calculate_fps(df, smi_colum_name)[0]
        descriptors = [fps] + [df[col_name].tolist() for col_name in cont_col]
    print('Done!')
    
    lf.batch_add(fps)
    lf.index()

    lf.store(f"{analysis_folder}lf.dat")
    with open(f"{analysis_folder}props.pickle", "wb+") as f:
        pickle.dump(
            tuple(descriptors[1:]),
            f,
            protocol=pickle.HIGHEST_PROTOCOL,
        )
    
    cat_labels, tmap_labels, tmap_cat, tmap_color, descriptors_column_names = [], [], [], [], []
    if cat_col:
        for cc in cat_col:
            tl, td = Faerun.create_categories(df[cc])
            if len(tl) > 9: # only show in legend the top 9 classes + "Other" Class
                td, tl = collapse_to_top_k(tl,td)
            cat_labels.append(tl)
            tmap_labels.append(td)
            tmap_cat.append(True)
            tmap_color.append("tab10")
            descriptors_column_names.append(cc)
      
    if cont_col:
        for cc in cont_col:
            tmap_labels.append(df[cc].tolist())
            tmap_cat.append(False)
            tmap_color.append("rainbow")
            descriptors_column_names.append(cc)
    
    if not cat_col and not cont_col:
        raise ValueError(
                "Please specify at least one descriptor to be plotted. "
                "Use either --cont_col and/or --cat_col flags."
                        )
        
    cfg = tm.LayoutConfiguration()
    cfg.node_size = 1 / 26
    cfg.mmm_repeats = 2
    cfg.sl_extra_scaling_steps = 5
    cfg.k = 20
    cfg.sl_scaling_type = tm.RelativeToAvgLength
    
    print('We are almost there...')
    x, y, s, t, _ = tm.layout_from_lsh_forest(lf, cfg)

    # tab_10 = plt.cm.get_cmap("tab10")
    # colors = [i for i in tab_10.colors]
    # colors[7] = (0.17, 0.24, 0.31)
    # tab_10.colors = tuple(colors)
    
    df[smi_colum_name] = (
                        df[smi_colum_name]
                        + '__'
                        + df[molecule_id_col]   # this is the visible text inside the link
                        # + '__'                  # separator for a new field
                        # + df[molecule_id_col]   # plain FoodB ID as an extra field
                    )
    
    all_desc = cat_col+cont_col
    if len(all_desc) > 0:    
        for col in all_desc:
            df[smi_colum_name] = df[smi_colum_name] + f' | {col}: ' + df[col].astype(str) 
        
    f = Faerun(view="front", coords=False)
    f.add_scatter(
        "np_atlas",
        {
            "x": x,
            "y": y,
            "c": tmap_labels, # <- here the df columns as a list of lists containing continuous and categorical variables)
            "labels": df[smi_colum_name],
        },
        shader="smoothCircle",
        point_scale=2.0,
        max_point_size=20,
        legend_labels=cat_labels, # list of the levels of the categorical variables ("type_labels" like objects generated with Faerun.create_categories)
        categorical=tmap_cat, # booleans defining continous and categorical variables from the "c" key in the dictionary
        colormap=tmap_color, # colors that will correspond to the variables in the "c" key  dictionary
        series_title=descriptors_column_names, # name of the different labels in the visualization tree ("c" keys too) 
        has_legend=True,
        title_index=1, # shows in the tooltip the image of the molecule and the id (CHECK  LINES 79-88)
        legend_title = 'Compound_Set'
    )
    f.add_tree("np_atlas_tree", {"from": s, "to": t}, point_helper="np_atlas")
    f.plot(template="smiles", path = analysis_folder)
    
    print(f'Check the resulting t-maps (index.html file) saved on {analysis_folder}')
    
    
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="T-MAPS generator")
    
    parser.add_argument("--analysis_folder", type=str, required=True,
        help="Path to the folder where the t-map is to be saved.")
   
    parser.add_argument("--data_path", type=str, required=True,
        help="Path to input dataset file.")
    
    parser.add_argument("--smi_col", type=str, required=True,
        help="Name of the column containing SMILES strings.")
    
    parser.add_argument("--molecule_id_col", type=str, required=True,
        help="Name of the column for molecule id.")
    
    parser.add_argument("--cat_col", type=str,
        help="Name of the columns for categorical descriptors separated by coma")
   
    parser.add_argument("--cont_col", type=str,
        help="Name of the columns for continous descriptors separated by coma")
    
    parser.add_argument("--max_cat", type=int, default= 10,
        help="Name of the columns for continous descriptors separated by coma")
    
    args = parser.parse_args()
    
    analysis_folder = os.path.join(args.analysis_folder, "")
    if not os.path.exists(analysis_folder):
        os.mkdir(analysis_folder)
        
    data_path = args.data_path
    smi_colum_name = args.smi_col.strip()
    molecule_id_col = args.molecule_id_col
    max_cat = args.max_cat
    
    cont_col = args.cont_col
    cont_col = (
                 [i.strip() for i in cont_col.split(",")]
                 if args.cont_col else []
               )

    cat_col = args.cat_col
    cat_col = (
                [i.strip() for i in cat_col.split(",")]
                if args.cat_col else []
              )

    main()
    
