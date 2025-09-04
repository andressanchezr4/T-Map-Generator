# T-Map-Generator
This code is a generalization for Chemical Space Visualization with tmap following the example "Natural Product Atlas" found [here](https://tmap.gdb.tools/#ex-npa).

The script takes a dataset with a SMILES column and plots all the molecules in a T-Map.

Molecule descriptors (molwt, tpsa, logp, class) can also be represented.

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
python path/to/tmap_gen.py --analysis_folder /path/to/results/folder --data_path path/to/csv --smi_col smiles --molecule_id_col mol_id --cont_col descriptor_col1,descriptor_col1 --cat_col descriptor_col3,descriptor_col4
```

## Disclaimer
Some fingerprint generation arises issues if you are using the latest version of rdkit. Even though therese errors are handled, those molecules will not be displayed. If you want to make sure every molecule is included you should have rdkit=2022.09.1 installed.
