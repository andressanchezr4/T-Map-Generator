# T-Map-Generator
This code is a generalization of the tmap for the example "Natural Product Atlas" found [here](https://tmap.gdb.tools/#ex-npa).

The script takes a dataset with a SMILES column, calculates their descriptors if necessary, and plots all the molecules in a T-Map. 
The final figure is a index.html where the T-Map plot can be explored. 

# Requirements 
- numpy
- pandas
- matplotlib
- rdkit
- faerum
- tmap
- mhfp

# Disclaimer
When descriptor columns are not specified, descriptors are calculated for MolWt, TPSA, LogP, NumHAcceptors and NumHDonors. 

Even though more than one categorical descriptor can be included in the T-Map, it has not been implemented. Follow the comments within the script in order to know where those can be added. 
