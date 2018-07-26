# LigParGen_2.3


Python package to convert BOSS generated OPLS-AA/CM1A(-LBCC) & OPLS-AA/CM5 parameters to:

- OpenMM PDB AND XML FILES,
- CHARMM RTF AND PRM FILES,
- Gromacs/NAMD ITP AND GRO FILES,
- PDB2PQR PQR FILES,  
- MCPRO & BOSS ZMATRIX
- LAMMPS lmp and data FILES
- tinker xyz and key files 
- .Q.param and  lib files for package Q

### Whats New ? ###

** OPLS-AA/CM5 FF parameters  with ORCA log file as input ** 

### authors: ###

* [Leela S. Dodda](https://github.com/leelasd) - `<leela.dodda@yale.edu>`
* [Israel Cabeza de Vaca]() - `<israel.cabezadevaca@yale.edu>`
* [Ayan Bhattacharjee]() - `<abhattacharjee.me@gmail.com>`
* [Matt Robinson](https://github.com/mc-robinson) 

### REQUIREMENTS: ###
- BOSS (need to set BOSSdir in bashrc and cshrc)
- Preferably Anaconda python with following modules
- pandas 
- argparse
- numpy
- openbabel (for 2D to 3D conversion)
- RDKit when CM5 charges are requested

### Installation instructions ###

`python setup.py install`

### Usage 

if using BOSS Zmat:
`LigParGen -z phenol.z -r PHN -c 0 -o 0 -l` 

if using MOL file:
`LigParGen -m phenol.mol -r PHN -c 0 -o 0 -l`

if using PDB file:
`LigParGen -p phenol.pdb -r PHN -c 0-o  0 -l`

if using BOSS SMILES CODE: 
`LigParGen -s 'c1ccc(cc1)O' -r PHN -c 0 -o 0 -l` 



### Who do I talk to? ###

* Leela S. Dodda leela.dodda@yale.edu 
* Israel Cabeza de Vaca 
* Matthew Robinson 

