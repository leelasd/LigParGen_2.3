"""
SCRIPT TO WRITE GROMACS GRO AND ITP FILES 
FROM BOSS ZMATRIX
Created on Mon Feb 15 15:40:05 2016
@author: Leela S. Dodda leela.dodda@yale.edu
@author: William L. Jorgensen Lab 

Usage: python OPM_Routines.py -z phenol.z -r PHN
REQUIREMENTS:
BOSS (need to set BOSSdir in bashrc and cshrc)
Preferably Anaconda python with following modules
pandas 
argparse
numpy
"""

from LigParGen.BOSSReader import ucomb,pairing_func,bossPdbAtom2Element,bossElement2Mass
from LigParGen.BOSSReader import Refine_PDB_file,get_coos_from_pdb
import pickle
import pandas as pd
import numpy as np


def GMX_pairs(tors, ang, bond):
    atom_list = np.array(list(bond.cl1) + list(bond.cl2))
    atom_list = np.sort(atom_list)
    atom_list = np.unique(atom_list)
    dict_bond = {ano: list(bond[bond.cl1 == ano]['cl2']) +
                 list(bond[bond.cl2 == ano]['cl1']) for ano in atom_list}
    dict_ang = {}
    for ano in atom_list:
        tr_1 = []
        for bno in dict_bond[ano]:
            tr_1 += dict_bond[bno]
        tr_1 = list(set(tr_1))
        tr_1.remove(ano)
        dict_ang[ano] = tr_1
    NP_dat = []
    for a in dict_bond.keys():
        for b in dict_bond[a]:
            NP_dat.append([a, b, 1, pairing_func(a + 1, b + 1)
                           [0], pairing_func(a + 1, b + 1)[1]])
            for c in dict_ang[b]:
                NP_dat.append([a, c, 3, pairing_func(a + 1, c + 1)
                               [0], pairing_func(a + 1, c + 1)[1]])
    for a in dict_ang.keys():
        for b in dict_ang[a]:
            NP_dat.append([a, b, 2, pairing_func(a + 1, b + 1)
                           [0], pairing_func(a + 1, b + 1)[1]])
            for c in dict_bond[b]:
                NP_dat.append([a, c, 3, pairing_func(a + 1, c + 1)
                               [0], pairing_func(a + 1, c + 1)[1]])
    NP_df = pd.DataFrame(NP_dat, columns=['I', 'J', 'BSEP', 'UNQ', 'TY'])
    NP_df = NP_df[NP_df.I != NP_df.J]  # CASES WITH 3 MEMBERED RINGS
    NP_df = NP_df.sort_values(['UNQ'])
    NP_B = NP_df[NP_df.BSEP == 1]
    NP_A = NP_df[NP_df.BSEP == 2]
    NP_T = NP_df[NP_df.BSEP == 3]
    NP_T = NP_T.drop_duplicates(['UNQ'])
    NP_T = NP_T[~ NP_T.UNQ.isin(NP_B.UNQ)]
    NP_T = NP_T[~ NP_T.UNQ.isin(NP_A.UNQ)]
    return list(NP_T.TY)


def gmxDihed(df):
    '''
        PRINTING RYCKET-BELLMAN STYLE PARAMETERS BECAUSE FEP DOES NOT SUPPORT PROPER DIHEDRALS
    '''
    [f1, f2, f3, f4] = [df.V1, df.V2, df.V3, df.V4]
    cdat = [f2 + (f1 + f3) * 0.5, 1.5 * f3 - 0.5 * f1, 4.0 *
            f4 - f2, -2.0 * f3, -4.0 * f4, 0.0, 0.00]
    return ('%5s%5s%5s%5s        3      %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f\n' % (
        df['I'] + 1, df['J'] + 1, df['K'] + 1, df['L'] + 1, cdat[0], cdat[1], cdat[2], cdat[3], cdat[4], cdat[5]))


def gmxImp(df):
    odihed = []
    for pot in range(1, 5):
        if df['V' + str(pot)] != 0.00:
            odihed.append('%6d%6d%6d%6d    4     %10.3f %10.3f %5d  \n' % (
                df['I'] + 1, df['J'] + 1, df['K'] +
                1, df['L'] + 1, 180.00 * abs(pot % 2 - 1),
                float(df['V' + str(pot)]) * 0.5, pot))
    return (odihed)


def boss2opmTorsion(bnd_df, num2opls, st_no, molecule_data, itpf):
    dhd = []
    for line in molecule_data.MolData['TORSIONS']:
        dt = [float(f) for f in line]
        dhd.append(dt)
    dhd = np.array(dhd)
    dhd = dhd * 4.184  # kcal to kj conversion
    dhd = dhd  # Komm = Vopls/2
    dhd_df = pd.DataFrame(dhd, columns=['V1', 'V2', 'V3', 'V4'])
    ats = []
    for line in molecule_data.MolData['ATOMS'][3:]:
        dt = [line.split()[0], line.split()[4],
              line.split()[6], line.split()[8]]
        dt = [int(d) for d in dt]
        ats.append(dt)
    for line in molecule_data.MolData['ADD_DIHED']:
        dt = [int(l) for l in line]
        ats.append(dt)
    assert len(ats) == len(
        dhd), 'Number of Dihedral angles in Zmatrix and Out file dont match'
    ats = np.array(ats) - st_no
    for i in range(len(ats)):
        for j in range(len(ats[0])):
            if ats[i][j] < 0:
                ats[i][j] = 0
    at_df = pd.DataFrame(ats, columns=['I', 'J', 'K', 'L'])
    final_df = pd.concat([dhd_df, at_df], axis=1)
    final_df = final_df.reindex(at_df.index)
    bndlist = list(bnd_df.UR) + (list(bnd_df.UR))
    final_df['TY'] = ['Proper' if ucomb(list([final_df.I[n], final_df.J[n], final_df.K[
        n], final_df.L[n]]), bndlist) == 3 else 'Improper' for n in range(len(final_df.I))]
    final_df['SumV'] = np.abs(
        final_df.V1) + np.abs(final_df.V2) + np.abs(final_df.V3) + np.abs(final_df.V4)
    final_df['TI'] = [num2opls[j] for j in final_df.I]
    final_df['TJ'] = [num2opls[j] for j in final_df.J]
    final_df['TK'] = [num2opls[j] for j in final_df.K]
    final_df['TL'] = [num2opls[j] for j in final_df.L]
    if len(final_df.index) > 0:
        final_df['NAME'] = final_df.TI + '-' + final_df.TJ + \
            '-' + final_df.TK + '-' + final_df.TL
        final_df = final_df.sort_values(['NAME'])
        tor_bos = final_df.drop(
            ['I', 'J', 'K', 'L', 'TI', 'TJ', 'TK', 'TL'], axis=1)
        tor_bos = tor_bos.drop_duplicates()
        df = final_df.iloc[tor_bos.index]
        return final_df, df
    else:
        return final_df, final_df


def boss2gmxBond(molecule_data, st_no, itpf):
    bdat = molecule_data.MolData['BONDS']
    bdat['cl1'] = [x - st_no if not x - st_no < 0 else 0 for x in bdat['cl1']]
    bdat['cl2'] = [x - st_no if not x - st_no < 0 else 0 for x in bdat['cl2']]
    bnd_df = pd.DataFrame(bdat)
    bnd_df['KIJ'] = bnd_df['KIJ'] * 836.80
    bnd_df['RIJ'] = bnd_df['RIJ'] * 0.10
    bnd_df['UNQ'] = [pairing_func(i + 1, j + 1)[0]
                     for i, j in zip(bnd_df.cl1, bnd_df.cl2)]
    bnd_df['UNQ'] = [pairing_func(i, j)[0]
                     for i, j in zip(bnd_df.cl1, bnd_df.cl2)]
    bnd_df['UF'] = ((bnd_df.cl1 + bnd_df.cl2) *
                    (bnd_df.cl1 + bnd_df.cl2 + 1) * 0.5) + bnd_df.cl2
    bnd_df['UR'] = ((bnd_df.cl1 + bnd_df.cl2) *
                    (bnd_df.cl1 + bnd_df.cl2 + 1) * 0.5) + bnd_df.cl1
    connects = []
    for (ai, aj) in zip(bnd_df.cl1, bnd_df.cl2):
        connects.append("{:>5}{:>5}    1".format(aj + 1, ai + 1))
    full_bnd = bnd_df.copy()
    bnd_df = bnd_df.drop_duplicates(['TIJ'])
    return full_bnd, connects


def boss2gmxAngle(anglefile, num2opls, st_no, itpf):
    adat = anglefile
    adat['cl1'] = [x - st_no if not x - st_no < 0 else 0 for x in adat['cl1']]
    adat['cl2'] = [x - st_no if not x - st_no < 0 else 0 for x in adat['cl2']]
    adat['cl3'] = [x - st_no if not x - st_no < 0 else 0 for x in adat['cl3']]
    ang_df = pd.DataFrame(adat)
    ang_df = ang_df[ang_df.K > 0]
    ang_df['K'] = 8.3680 * ang_df['K']
    ang_df['R'] = ang_df['R']
    full_ang = ang_df.copy()
    ang_df['TY'] = np.array([num2opls[i] + '-' + num2opls[j] + '-' + num2opls[k]
                             for i, j, k in zip(ang_df.cl1, ang_df.cl2, ang_df.cl3)])
    ang_df = ang_df.drop_duplicates(['TY'])
    return full_ang


def bossData(molecule_data):
    ats_file = molecule_data.MolData['ATOMS']
    types = []
    for i in enumerate(ats_file):
        types.append([i[1].split()[1], 'opls_' + i[1].split()[2]])
    st_no = 3
    Qs = molecule_data.MolData['Q_LJ']
    assert len(Qs) == len(types), 'Please check the at_info and Q_LJ_dat files'
    num2opls = {}
    for i in range(0, len(types)):
        num2opls[i] = Qs[i][0]
    num2typ2symb = {i: types[i] for i in range(len(Qs))}
    for i in range(len(Qs)):
        num2typ2symb[i].append(bossPdbAtom2Element(
            num2typ2symb[i][0]) + num2typ2symb[i][1][-3:])
        num2typ2symb[i].append(bossPdbAtom2Element(num2typ2symb[i][0]))
        num2typ2symb[i].append(bossElement2Mass(num2typ2symb[i][3]))
        num2typ2symb[i].append(Qs[i][0])
    return (types, Qs, num2opls, st_no, num2typ2symb)


def pdb2gro(atoms, coos, resid):
    grof = open(resid + '.gro', 'w+')
    grof.write('LIGPARGEN GENERATED GRO FILE\n')
    grof.write('%5d\n' % len(atoms))
    num = 0
    for (i, j) in zip(atoms, coos):
        num += 1
        grof.write('%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n' %
                   (1, resid, i, num, j[0] * 0.1, j[1] * 0.1, j[2] * 0.1))
    grof.write('%10.5f%10.5f%10.5f\n' % (1.00000, 1.00000, 1.00000))
    grof.close()
    return None


def boss2gmxAtom(resid, num2typ2symb, Qs, itpf):
    itpf.write('[ atoms ]\n')
    itpf.write(
        ';   nr       type  resnr residue  atom   cgnr     charge       mass  \n')
    for i in range(len(Qs)):
        itpf.write(' %5d %10s %6d %6s %5s %6d %10s %10.4f \n' % (
            i + 1, num2typ2symb[i][1], 1, resid, num2typ2symb[i][0], ((i + 1) / 33.0) + 1, Qs[i][1], num2typ2symb[i][4]))
    return None


def boss2gmx(resid, molecule_data, pdb_file):
    types, Qs, num2opls, st_no, num2typ2symb = bossData(
        molecule_data)
    itpf = open(resid + '.itp', 'w+')
    itpf.write("""
; Writted by Leela S. Dodda leela.dodda@yale.edu
; GENERATED BY Server
; Jorgensen Lab @ Yale University 
;
""")
    itpf.write('[ atomtypes ]\n')
    atty = ['%10s %5s %10.4f     0.000    A    %10.5E   %10.5E\n' % (num2typ2symb[i][
        1], num2typ2symb[i][2], num2typ2symb[i][4],
        float(Qs[i][2]) * 0.1, float(Qs[i][3]) * 4.184) for
        i in range(len(Qs))]
    atty = list(set(atty))
    for att in atty:
        itpf.write('%s' % att)
    bnd_df, connects = boss2gmxBond(molecule_data, st_no, itpf)
    full_ang = boss2gmxAngle(molecule_data.MolData[
        'ANGLES'], num2opls, st_no, itpf)
    itpf.write('[ moleculetype ]\n')
    itpf.write('; Name               nrexcl\n')
    itpf.write('%s                   3\n' % resid)
    boss2gmxAtom(resid, num2typ2symb, Qs, itpf)
    itpf.write('[ bonds ]\n')
    for brow in bnd_df.iterrows():
        index, bnd = brow
        itpf.write('%5d %5d %5d %11.4f %10.3f\n' %
                   (bnd['cl1'] + 1, bnd['cl2'] + 1, 1, bnd['RIJ'], bnd['KIJ']))
    itpf.write('\n[ angles ]\n')
    itpf.write(
        ';  ai    aj    ak funct            c0            c1            c2            c3 \n')
    for row in full_ang.iterrows():
        index, angs = row
        itpf.write('%5d %5d %5d %5d %10.3f %10.3f\n' %
                   (angs['cl1'] + 1, angs['cl2'] + 1, angs['cl3'] + 1, 1, angs['R'], angs['K']))
    full_tor, tor_df = boss2opmTorsion(
        bnd_df, num2opls, st_no, molecule_data, itpf)
    if len(tor_df.index) != len(full_tor.index):
        itpf.write('\n[ dihedrals ]\n')
        itpf.write('; IMPROPER DIHEDRAL ANGLES \n')
        itpf.write(
            ';  ai    aj    ak    al funct            c0            c1            c2            c3            c4            c5\n')
        for row in full_tor[full_tor.TY == 'Improper'].iterrows():
            index, tors = row
            for st in gmxImp(tors):
                itpf.write('%s' % st)
        itpf.write('\n[ dihedrals ]\n')
        itpf.write('; PROPER DIHEDRAL ANGLES\n')
        itpf.write(
            ';  ai    aj    ak    al funct            c0            c1            c2            c3            c4            c5\n')
        for row in full_tor[full_tor.TY == 'Proper'].iterrows():
            index, tors = row
            for st in gmxDihed(tors):
                itpf.write('%s' % st)
    ppairs = GMX_pairs(full_tor, full_ang, bnd_df)
    itpf.write('\n[ pairs ]\n')
    for pair in ppairs:
        itpf.write('%s    1\n' % pair)
    itpf.close()
    pdblines = Refine_PDB_file(pdb_file)
    atoms, coos = get_coos_from_pdb(pdblines)
    pdb2gro(atoms, coos, resid)
    return None


def mainBOSS2GMX(resid, clu):
    mol = pickle.load(open(resid + ".p", "rb"))
    if clu:
        pdb_file = 'clu.pdb'
    else:
        pdb_file = 'plt.pdb'
    boss2gmx(resid, mol, pdb_file)
    return None
