"""
SCRIPT TO WRITE Q LIB AND Q.PRM FILES 
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

from LigParGen.BOSSReader import bossPdbAtom2Element,bossElement2Mass,ucomb
import pickle
import pandas as pd
import numpy as np


def retDihedImp(df):
    odihed = []
    if np.sum([df['V' + str(pot)] for pot in range(1, 5)]) != 0.0:
        for pot in range(1, 5):
            if (df['V' + str(pot)] != 0.0):
                odihed.append('%s %8.5f %2d %4.5f \n' % (df['NAME'].replace(
                    "-", "   "), 2.0*df['V' + str(pot)], pot, 180.00 * abs(pot % 2 - 1)))
    else:
        pot = 2
        odihed.append('%s %8.5f %2d %4.5f \n' % (df['NAME'].replace(
            "-", "   "), 2.0*df['V' + str(pot)], pot, 180.00 * abs(pot % 2 - 1)))
    return (odihed)


def retDihed(df):
    odihed = []
    for pot in range(1, 5):
        if pot <4: odihed.append('%s           %8.5f %2d %4.5f 1\n' % (df['NAME'].replace("-", "   "), 2.0*df['V' + str(pot)], int(pot*-1), 180.00 * abs(pot % 2 - 1)))
        else: odihed.append('%s           %8.5f %2d %4.5f 1\n' % (df['NAME'].replace("-", "   "), 2.0*df['V' + str(pot)], pot, 180.00 * abs(pot % 2 - 1)))
    return (odihed)


def Boss2CharmmRTF(num2typ2symb, Qs, resid, bnd_df, angs,props,imps):
    charges = [float(Qs[i][1]) for i in range(len(Qs))]
    rtf = open(resid + '.lib', 'w+')
    rtf.write('#Remarks generated for Q (by Leela Dodda)\n')
    rtf.write('\n{%s}\n'%resid)
    rtf.write('[atoms]\n')
    for i in range(len(Qs)):
        rtf.write('%8d %6s %6s %8s \n' % (i+1,num2typ2symb[i][0], bossPdbAtom2Element(
            num2typ2symb[i][0]) + num2typ2symb[i][1][-3:], Qs[i][1]))
    rtf.write('[bonds]')
    for (x, y) in zip(bnd_df.cl1, bnd_df.cl2):
        rtf.write('\n%8s %8s' % (num2typ2symb[x][0], num2typ2symb[y][0]))
    rtf.write('\n[impropers]\n')
    for i in imps:
        rtf.write('%s \n' % (i.replace("-", "    ")))
    rtf.write('[charge_groups]\n')
    for i in range(len(Qs)):rtf.write('%4s'%num2typ2symb[i][0])
    rtf.write('\n*------------------------------------------------------------------\n')
    rtf.close()
    return None


def Boss2CharmmPRM(resid, num2typ2symb, Qs, bnd_df, ang_df, tor_df):
    #### COLLECTING NONBONDING PART #######
    prm = open(resid + '.Q.prm', 'w+')
    prm.write('# generated Q-PARAM file for Aqvist group (by Leela Dodda)\n')
    prm.write('\n[options]\n')
    prm.write('\n[atom_types]\n')
    for i in range(len(Qs)): 
        eps = float(Qs[i][3])
        sig = float(Qs[i][2])
        ALJ = 2*sig**6*np.sqrt(eps)
        BLJ = 2*sig**3*np.sqrt(eps)
        prm.write('%4s %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n'%(num2typ2symb[i][2],ALJ,ALJ,BLJ,BLJ,0.000,bossElement2Mass(bossPdbAtom2Element(num2typ2symb[i][0]))))
    prm.write('\n[bonds]\n')
    for i in bnd_df.index:
        prm.write('%s %6s %8.1f %8.4f \n' % (num2typ2symb[bnd_df.cl1[i]][
            2], num2typ2symb[bnd_df.cl2[i]][2], 2.0*bnd_df.KIJ[i], bnd_df.RIJ[i]))
    prm.write('\n[angles]\n')
    for i in ang_df.index:
        prm.write('%s %5s %5s %8.1f %8.2f \n' % (num2typ2symb[ang_df.cl1[i]][2], num2typ2symb[
            ang_df.cl2[i]][2], num2typ2symb[ang_df.cl3[i]][2], 2.0*ang_df.K[i], ang_df.R[i]))
    prm.write('\n[torsions]\n')
    if len(tor_df.index) > 0:
        tor_df = tor_df.drop_duplicates(['NAME', 'TY'])
    pro_df = tor_df[tor_df.TY == 'Proper']
    for i in list(pro_df.index):
        ndf = pro_df.iloc[i]
        pro_out = retDihed(ndf.to_dict())
        for i in range(4):
            prm.write('%s' % pro_out[i])
    prm.write(
        '!   X    X    X    X    0.00000 1 0.000000 ! WILD CARD FOR MISSING TORSION PARAMETERS\n')
    prm.write('\n[impropers]\n')
    imp_df = tor_df[tor_df.TY == 'Improper']
    for i in list(imp_df.index):
        ndf = tor_df.iloc[i]
        imp_out = retDihedImp(ndf.to_dict())
        for i in range(len(imp_out)):
            prm.write('%s' % imp_out[i])
    prm.write(
        '!   X    X    X    X    0.00000 2 0.000000 ! WILD CARD FOR MISSING IMPROPER PARAMETERS \n')
    prm.close()
    return None


def Boss2CharmmTorsion(bnd_df, num2opls, st_no, molecule_data, num2typ2symb):
    dhd = []
    for line in molecule_data.MolData['TORSIONS']:
        dt = [float(l) for l in line]
        dhd.append(dt)
    dhd = np.array(dhd)
    dhd = dhd  # kcal to kj conversion
    dhd = dhd / 2.0  # Komm = Vopls/2
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
    final_df['TI'] = [num2typ2symb[j][2] for j in final_df.I]
    final_df['TJ'] = [num2typ2symb[j][2] for j in final_df.J]
    final_df['TK'] = [num2typ2symb[j][2] for j in final_df.K]
    final_df['TL'] = [num2typ2symb[j][2] for j in final_df.L]
    final_df['SYMB'] = ['   '.join([num2typ2symb[final_df.I[i]][0], num2typ2symb[final_df.J[i]][
        0], num2typ2symb[final_df.K[i]][0], num2typ2symb[final_df.L[i]][0]]) for i in final_df.index]
    if len(final_df.index) > 0:
        final_df['NAME'] = final_df.TI + '-' + final_df.TJ + \
            '-' + final_df.TK + '-' + final_df.TL
    return final_df


def boss2CharmmBond(molecule_data, st_no):
    bdat = molecule_data.MolData['BONDS']
    bdat['cl1'] = [x - st_no if not x - st_no < 0 else 0 for x in bdat['cl1']]
    bdat['cl2'] = [x - st_no if not x - st_no < 0 else 0 for x in bdat['cl2']]
    bnd_df = pd.DataFrame(bdat)
    bnd_df['UF'] = ((bnd_df.cl1 + bnd_df.cl2) *
                    (bnd_df.cl1 + bnd_df.cl2 + 1) * 0.5) + bnd_df.cl2
    bnd_df['UR'] = ((bnd_df.cl1 + bnd_df.cl2) *
                    (bnd_df.cl1 + bnd_df.cl2 + 1) * 0.5) + bnd_df.cl1
    hb_df = bnd_df.drop(['cl1', 'cl2', 'UF', 'UR'], axis=1)
    hb_df = hb_df.drop_duplicates()
    return bnd_df


def boss2CharmmAngle(anglefile, num2opls, st_no,num2typ2symb):
    adat = anglefile
    adat['cl1'] = [x - st_no if not x - st_no < 0 else 0 for x in adat['cl1']]
    adat['cl2'] = [x - st_no if not x - st_no < 0 else 0 for x in adat['cl2']]
    adat['cl3'] = [x - st_no if not x - st_no < 0 else 0 for x in adat['cl3']]
    ang_df = pd.DataFrame(adat)
    ang_df = ang_df[ang_df.K > 0]
    ang_df['TY'] = np.array([num2opls[i] + '-' + num2opls[j] + '-' + num2opls[k]
                             for i, j, k in zip(ang_df.cl1, ang_df.cl2, ang_df.cl3)])
    ang_df['TI']=[num2typ2symb[ang_df.cl1[i]][2] for i in ang_df.index]
    ang_df['TJ']=[num2typ2symb[ang_df.cl2[i]][2] for i in ang_df.index]
    ang_df['TK']=[num2typ2symb[ang_df.cl3[i]][2] for i in ang_df.index]
    ang_df['TY'] = np.array([i + '  ' + j + '  ' + k
                             for i, j, k in zip(ang_df.TI, ang_df.TJ, ang_df.TK)])
    return ang_df


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


def Boss2Charmm(resid, molecule_data):
    types, Qs, num2opls, st_no, num2typ2symb = bossData(molecule_data)
    bnd_df = boss2CharmmBond(molecule_data, st_no)
    ang_df = boss2CharmmAngle(molecule_data.MolData['ANGLES'], num2opls, st_no,num2typ2symb)
    tor_df = Boss2CharmmTorsion(bnd_df, num2opls, st_no,
                                molecule_data, num2typ2symb)
    Boss2CharmmRTF(num2typ2symb, Qs, resid, bnd_df, list(ang_df['TY']), list(
        tor_df[tor_df.TY == 'Proper']['SYMB']),list(tor_df[tor_df.TY == 'Improper']['SYMB']))
    Boss2CharmmPRM(resid, num2typ2symb, Qs, bnd_df, ang_df, tor_df)
    return None


def mainBOSS2Q(resid, clu=False):
    mol = pickle.load(open(resid + ".p", "rb"))
    Boss2Charmm(resid, mol)
    return None
