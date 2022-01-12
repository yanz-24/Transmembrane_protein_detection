# TRANSMEMBRANE PROTEIN DETECTION

This python program identifies transmembrane regions of proteins, some of its features include:
- Detecting if a protein is globular or transmembrane.
- Finding the orientation of the protein in a membrane. 
- Indicating the position of the lipid bi-layer.


## Development environment
This program was developped by python 3.9.5. 

Conda is needed in order to set up the environment. For conda installation please check the user manual 
[https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)  
Run the following commands to create environment: 
```
$ conda env create -f TMDET.yml 
$ conda activate TMDET
```
To calculate secondary structure and accessibility, you need to download DSSP (see https://swift.cmbi.umcn.nl/gv/dssp/)  
To install DSSP on Linux run the following command:
```
$ conda activate TMDET
$ sudo apt-get install dssp
```

Usage
-----------------------------------------------------------------------------
- required arguments:  
  -i I        input file (pdb)  
  -id ID      pdb id
- optional arguments:  
  -h, --help  show this help message and exit  
  -pt PT      number of points in fibonnaci sphere (default=100)
  -o O        output pdb file (protein + bi-layer position)
run in terminal (example):
```
$ python3 main.py -id 1XYZ 
$ python3 main.py -i 1XYZ.pdb -o output.pdb
```
- Example files for 1uaz proteins are given as part of the package.
- The output pdb files  that include bi layer position can be visualized with appropriate programs such as pymol.

Authors
-----------------------------------------------------------------------------
- Sabrina SAFAR-REMALI  
- Yanyuan ZHANG

References
-----------------------------------------------------------------------------
- This work was inspired by the algorithm described in the following  article:
- Tusnády GE, Dosztányi Z, Simon I. Transmembrane proteins in the Protein Data Bank: identification and classification. Bioinformatics. 2004 Nov 22;20(17):2964-72. doi: 10.1093/bioinformatics/bth340. Epub 2004 Jun 4. PMID: 15180935.

GitLab repository
-----------------------------------------------------------------------------
You can find our source code repository and demo files on GitLab 
at [https://gitlab.com/sab.s/python2_project](https://gitlab.com/sab.s/python2_project). 


