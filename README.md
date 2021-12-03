TRANSMEMBRANE PROTEIN DETECTION
=============================================================================

This python program identifies transmembrane regions of proteins, Some of its features include:

-

Development environment                                      
-----------------------------------------------------------------------------
This program was developped by python 3.8.8 using pandas, biopython and numpy
package as well as argparse module. 

Conda is needed in order to set up the environment
for conda installation please check the user manual [https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)
Then run the following commands : \
$ **conda env create -f TMDET.yml** \
$ **conda activate TMDET** \

Usage
-----------------------------------------------------------------------------
- required arguments:  
  -i I        input file (pdb)  
  -id ID      pdb id 

- optional arguments:  
  -h, --help  show this help message and exit  
  -pt PT      number of points in fibonnaci sphere (default=1000)

run in terminal (example):
```
python3 main.py -i file.pdb -id 1XYZ
```

Authors
-----------------------------------------------------------------------------
Sabrina SAFAR-REMALI
Yanyuan ZHANG


GitLab
-----------------------------------------------------------------------------
You can find our source code repository and demo files on GitLab 
at [https://gitlab.com/sab.s/python2_project](https://gitlab.com/sab.s/python2_project). 
