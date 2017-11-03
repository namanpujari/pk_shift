# pk_shift
---
Python based script to visualize the surface energy of a protein. 
Can be used to analyze the nature of the change in energy. pH-Energy
dependency rules have been used to produce a matplotlib figure. 

Dependencies: numpy, matplotlib, pandas. 
Simply run `python main.py` to run the program. YOU MUST HAVE A VALID pK.out
file from which to read the data. The file source is hard-coded so pK.out must
be in the folder where main.py is located, and must not be renamed into
anything else. 

### Results

Produces matplotlib figures, one with 4 subplots indicating individual energy
changes by pH for each residue, and one indicating the aggregate for
negative/positive residues. Second figure is the aggregate of the energy
changes for the protein as a whole.
