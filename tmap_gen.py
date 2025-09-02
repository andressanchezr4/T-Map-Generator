# -*- coding: utf-8 -*-
"""
Created on 2025

@author: andres.sanchez
"""

import pickle
import numpy as np
import tmap as tm
import pandas as pd
from rdkit.Chem import AllChem
from mhfp.encoder import MHFPEncoder
from faerun import Faerun
from matplotlib import pyplot as plt
from rdkit.Chem import Descriptors
from rdkit import Chem
import os
import argparse
from Tools import smi_descriptors, calculate_descriptors, fps_mhfp_gen
 
def main(descriptors_column_names):
    
    df = pd.read_csv(data_path, sep=",", low_memory= False)

    lf = tm.LSHForest(1024, 64)
    
    if not descriptors_column_names:
        descriptors = calculate_descriptors(df, smi_colum_name)
        fps = descriptors[0]
        descriptors_column_names = ['molwt', 'tpsa', 'logp', 'ha', 'hd']
        
    elif descriptors_column_names:
        print('Calculating fps... This might take a while')
        df['fps'] = df[smi_colum_name].apply(fps_mhfp_gen)
        df = df.dropna()
        fps = df['fps'].tolist()
        cols2descriptors = ['fps'] + descriptors_column_names
        descriptors = [df[col_name].tolist() for col_name in cols2descriptors]
        print('Done with fps generation!')
        
    lf.batch_add(fps)
    lf.index()

    lf.store(f"{analysis_folder}lf.dat")
    with open(f"{analysis_folder}props.pickle", "wb+") as f:
        pickle.dump(
            tuple(descriptors[1:]),
            f,
            protocol=pickle.HIGHEST_PROTOCOL,
        )
        
    if molecule_class_col:
        type_labels, type_data = Faerun.create_categories(df[molecule_class_col]) # here to add more categorical variables to the tmap
        tmap_labels = [type_data] + list(descriptors[1:]) # add more "type_data" objects to the first list 
        tmap_cat = [True] + [False] * len(descriptors[1:])
        tmap_color = ["tab10"] + ['rainbow'] * len(descriptors[1:])
        descriptors_column_names = [molecule_class_col] + descriptors_column_names
    else:
        type_labels = []
        tmap_labels = list(descriptors[1:])
        tmap_cat = [False] * len(descriptors[1:])
        tmap_color = ['rainbow'] * len(descriptors[1:])
        
    cfg = tm.LayoutConfiguration()
    cfg.node_size = 1 / 26
    cfg.mmm_repeats = 2
    cfg.sl_extra_scaling_steps = 5
    cfg.k = 20
    cfg.sl_scaling_type = tm.RelativeToAvgLength
    x, y, s, t, _ = tm.layout_from_lsh_forest(lf, cfg)

    tab_10 = plt.cm.get_cmap("tab10")
    colors = [i for i in tab_10.colors]
    colors[7] = (0.17, 0.24, 0.31)
    tab_10.colors = tuple(colors)
    
    f = Faerun(view="front", coords=False)
    f.add_scatter(
        "np_atlas",
        {
            "x": x,
            "y": y,
            "c": tmap_labels, # <- here the categorical and continuous variables
            "labels": df[smi_colum_name],
        },
        shader="smoothCircle",
        point_scale=2.0,
        max_point_size=20,
        legend_labels=[type_labels], # if you want to include more categorical variables, here you have to add more "type_labels" like objects to the list (generated with Faerun.create_categories)
        categorical=tmap_cat, # booleans defining continous and categorical variables from the "c" key in the dictionary
        colormap=tmap_color, # colors that will correspond to the "c" key in the dictionary
        series_title=descriptors_column_names, # name of the different labels in the visualization tree ("c" keys too) 
        has_legend=True,
    )
    f.add_tree("np_atlas_tree", {"from": s, "to": t}, point_helper="np_atlas")
    f.plot(template="smiles", path = analysis_folder)
    
    print(f'Open the resulting t-maps (index.html file) saved on {analysis_folder}')
    
    
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="T-MAPS generator")
    
    parser.add_argument("--analysis_folder", type=str, required = True,
        help="Path to the folder where the t-map is to be saved.")
   
    parser.add_argument("--data_path", type=str, required = True,
        help="Path to input dataset file (csv format separated by , ).")
    
    parser.add_argument("--smi_column_name", type=str, required = True,
        help="Name of the column containing SMILES strings.")

    parser.add_argument("--molecule_class_col", type=str,
        help="Column name for molecule class if any.")
    
    parser.add_argument("--descriptors_column_names", default = None,
        help="List of calculated descriptors (if any) separated by coma (,). CATEGORICAL DESCRIPTORS NOT HANDLED.")
    
    args = parser.parse_args()

    analysis_folder = os.path.join(args.analysis_folder, "")
    if not os.path.exists(analysis_folder):
        os.mkdir(analysis_folder)
    
    data_path = args.data_path
    smi_colum_name = args.smi_column_name.strip()
    molecule_class_col = args.molecule_class_col
    descriptors_column_names = args.descriptors_column_names
    descriptors_column_names = (
            [i.strip() for i in descriptors_column_names.split(",")]
            if args.descriptors_column_names else None
        )

    main(descriptors_column_names)
