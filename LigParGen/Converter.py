from LigParGen.BOSSReader import BOSSReader, CheckForHs
from LigParGen.BOSS2OPENMM import mainBOSS2OPM
from LigParGen.BOSS2CHARMM import mainBOSS2CHARMM
from LigParGen.BOSS2GMX import mainBOSS2GMX
from LigParGen.BOSS2XPLOR import mainBOSS2XPLOR
from LigParGen.BOSS2Q import mainBOSS2Q
from LigParGen.BOSS2LAMMPS import mainBOSS2LAMMPS
from LigParGen.BOSS2DESMOND import mainBOSS2DESMOND 
from LigParGen.BOSS2TINKER import mainBOSS2TINKER 
from LigParGen.CreatZmat import GenMolRep
from LigParGen.Orca2CM5charges import LoadModel, GetLogFile, HirshfeldToCM5,AddCM5Charges
from LigParGen.fepzmat import CM5_file2zmat
import argparse
import pickle
import os

def main():

    parser = argparse.ArgumentParser(
        prog='LigParGen',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""
    Ligand Parameter Generator Based on 
    Jorgensen group's OPLS-AA/CM1A(-LBCC) FF
    Created on Mon Feb 15 15:40:05 2016
    @author: Leela S. Dodda leela.dodda@yale.edu
    @author: William L. Jorgensen Lab 

    FF formats provided : 
    --------------------
    OpenMM       .xml  
    CHARMM/NAMD  .prm & .rtf  
    GROMACS      .itp & .gro 
    CNS/X-PLOR   .param & .top
    Q            .Q.prm & .Q.lib
    DESMOND      .cms
    BOSS/MCPRO   .z
    PDB2PQR      .pqr

    Input Files supported : 
    --------------------
    SMILES code
    PDB
    MDL MOL Format

    ################################################ 
    if using MOL file 
    Usage: -m phenol.mol    -r PHN -c 0 -o 0

    if using PDB file 
    Usage: -p phenol.pdb    -r PHN -c 0 -o 0
    
    if using BOSS SMILES CODE 
    Usage: -s 'c1ccc(cc1)O' -r PHN -c 0 -o 0  
    
    REQUIREMENTS:
    BOSS (need to set BOSSdir in bashrc and cshrc)
    Preferably Anaconda python with following modules
    pandas 
    argparse
    numpy
    openbabel

    Please cite following references: 
    1. web server: an automatic OPLS-AA parameter generator for organic ligands  
       Leela S. Dodda  Israel Cabeza de Vaca  Julian Tirado-Rives William L. Jorgensen 
       Nucleic Acids Research, Volume 45, Issue W1, 3 July 2017, Pages W331-W336
    2. 1.14*CM1A-LBCC: Localized Bond-Charge Corrected CM1A Charges for Condensed-Phase Simulations
       Leela S. Dodda, Jonah Z. Vilseck, Julian Tirado-Rives , and William L. Jorgensen 
       Department of Chemistry, Yale University, New Haven, Connecticut 06520-8107, United States
       J. Phys. Chem. B, 2017, 121 (15), pp 3864-3870
    3. Accuracy of free energies of hydration using CM1 and CM3 atomic charges.
       Udier-Blagovic, M., Morales De Tirado, P., Pearlman, S. A. and Jorgensen, W. L. 
       J. Comput. Chem., 2004, 25,1322-1332. doi:10.1002/jcc.20059
    """
    )
    parser.add_argument(
        "-s", "--smiles", help="Paste SMILES code from CHEMSPIDER or PubChem", type=str)
    parser.add_argument(
        "-z", "--zmat", help="Submit Z-matrix in BOSS Zmat format", type=str)
    parser.add_argument(
        "-m", "--mol", help="Submit MOL file MDL/MOL Format with hydrogens added", type=str)
    parser.add_argument(
        "-p", "--pdb", help="Submit PDB file with Hydrogens Added", type=str)
    parser.add_argument(
        "-q", "--qorca", help="ORCA LOG FILE", type=str)
    parser.add_argument(
        "-r", "--resname", help="Residue name (Should be a 3 LETTER WORD)", type=str)
    parser.add_argument(
        "-o", "--opt", help="Optimization or Single Point Calculation", type=int, choices=[0, 1, 2, 3])
    parser.add_argument("-c", "--charge", type=int,
                        choices=[0, -1, 1, -2, 2], help="0: Neutral <0: Anion >0: Cation ")
    parser.add_argument(
        "-l", "--lbcc", help="Use 1.14*CM1A-LBCC charges instead of 1.14*CM1A", action="store_true")

    args = parser.parse_args()

    convert(**vars(args))

def convert(**kwargs):

    # set the default values
    options = {
            'opt' : 0,
            'smiles' : None,
            'zmat' : None, 
            'charge' : 0,
            'lbcc' : False,
            'mol' : None,
            'resname' : 'UNK',
            'pdb' : None,
            'qorca': None }

    # update the default values based on the arguments
    options.update(kwargs)

    # set the arguments that you would used to get from argparse
    opt = options['opt']
    smiles = options['smiles']
    zmat = options['zmat']
    charge = options['charge']
    lbcc = options['lbcc']
    resname = options['resname']
    mol = options['mol']
    pdb = options['pdb']
    qorca = options['qorca']
    if opt != None:
        optim = opt
    else:
        optim = 0

    clu = False

    # assert (which('obabel')
            # is not None), "OpenBabel is Not installed or \n the executable location is not accessable"
    if os.path.exists('/tmp/' + resname + '.xml'):
        os.system('/bin/rm /tmp/' + resname + '.*')
    if lbcc:
        if charge == 0:
            lbcc = True
            print('LBCC converter is activated')
        else:
            lbcc = False
            print(
                '1.14*CM1A-LBCC is only available for neutral molecules\n Assigning unscaled CM1A charges')
    if qorca != None:
        if optim > 0: print('CANNOT OPTIMIZE with CM5 CHARGES')
        optim = 0
        lbcc = False
        a0,rd,pt = LoadModel()
        os.system('cp  %s /tmp/'%qorca)
        os.chdir('/tmp/')
        data_cm5 = GetLogFile(qorca,pt,rd)
        qcm5 = HirshfeldToCM5(data_cm5,a0,netcharge=charge)
        os.system('cp inp_orca.pdb /tmp/%s.pdb' %resname)
        pdb = '%s.pdb'%resname

    if smiles != None:
        os.chdir('/tmp/')
        smifile = open('%s.smi' % resname, 'w+')
        smifile.write('%s' % smiles)
        smifile.close()
        GenMolRep('%s.smi' % resname, optim, resname, charge)
        mol = BOSSReader('%s.z' % resname, optim, charge, lbcc)
    elif mol != None:
        os.system('cp %s /tmp/' % mol)
        os.chdir('/tmp/')
        GenMolRep(mol, optim, resname, charge)
        mol = BOSSReader('%s.z' % resname, optim, charge, lbcc)
    elif pdb != None:
        os.system('cp %s /tmp/' % pdb)
        os.chdir('/tmp/')
        GenMolRep(pdb, optim, resname, charge)
        mol = BOSSReader('%s.z' % resname, optim, charge, lbcc)
        clu = True
    elif zmat != None:
        os.system('cp %s /tmp/%s.z' % (zmat,resname))
        os.chdir('/tmp/')
        print('THIS OPTION IS FOR SUPPLYING OPLS-AA Z-matrices only')
        if optim > 0: 
            print('CANNOT OPTIMIZE with Z-matrix Option')
            optim = 0
        if lbcc: print('CANNOT APPLY LBCC FOR A OPLS-AA Z-MATRIX')
        mol = BOSSReader('%s.z' % resname, optim, charge, lbcc=False)
    if qorca != None: 
        print('I am here')
        mol = AddCM5Charges(mol,qcm5)
        CM5_file2zmat('%s.z' % resname, qcm5.CM5_final,
                      oname='/tmp/%s_CM5.z' % resname)
        os.system('mv %s.z %s_CM1A.z' %
                  (resname,resname))
        os.system('mv %s_CM5.z %s.z' % (resname,resname))

    assert (mol.MolData['TotalQ']['Reference-Solute'] ==
            charge), "PROPOSED CHARGE IS NOT POSSIBLE: SOLUTE MAY BE AN OPEN SHELL"
    assert(CheckForHs(mol.MolData['ATOMS'])
           ), "Hydrogens are not added. Please add Hydrogens"

    pickle.dump(mol, open(resname + ".p", "wb"))
    mainBOSS2OPM(resname, clu)
    print('DONE WITH OPENMM')
    mainBOSS2Q(resname, clu)
    print('DONE WITH Q')
    mainBOSS2XPLOR(resname, clu)
    print('DONE WITH XPLOR')
    mainBOSS2CHARMM(resname, clu)
    print('DONE WITH CHARMM/NAMD')
    mainBOSS2GMX(resname, clu)
    print('DONE WITH GROMACS')
    mainBOSS2LAMMPS(resname, clu)
    print('DONE WITH LAMMPS')
    mainBOSS2DESMOND(resname, clu)
    print('DONE WITH DESMOND')
    mainBOSS2TINKER(resname, clu)
    print('DONE WITH TINKER')
    os.remove(resname + ".p")
    mol.cleanup()

if __name__ == "__main__":
  
    main()
