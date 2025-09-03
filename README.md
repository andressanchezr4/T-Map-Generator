# T-Map-Generator
This code is a generalization for Chemical Space Visualization with tmap following the example "Natural Product Atlas" found [here](https://tmap.gdb.tools/#ex-npa).

The script takes a dataset with a SMILES column, calculates their descriptors if necessary, and plots all the molecules in a T-Map. 
The final figure is a index.html where the T-Map plot can be explored. 

## Requirements 
- numpy
- pandas
- matplotlib
- rdkit
- faerum
- tmap
- mhfp

## Example of usage 

```bash
python path/to/tmap_gen.py --analysis_folder /path/to/results/folder --data_path path/to/csv --smi_column_name smiles_column --molecule_class_col class_column --molecule_id_col id_column
```

## Disclaimer
When descriptor columns are not specified, descriptors are calculated for MolWt, TPSA, LogP, NumHAcceptors and NumHDonors. 

Even though more than one categorical descriptor can be included in the T-Map, it has not been implemented. Follow the comments within the script in order to know where those can be added. 

Some fingerprint generation arises issues if you are using the latest version of rdkit. Even though therese errors are handled, those molecules will not be displayed. If you want to make sure every molecule is included you should have rdkit=2022.09.1 installed.
