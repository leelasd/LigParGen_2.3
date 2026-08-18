"""
SCRIPT TO CONVERT WRITE CNS/X-PLOR PARAM AND TOP FILES 
FROM BOSS ZMATRIX
Created on Mon Feb 15 15:40:05 2016
@author: Leela S. Dodda leela.dodda@yale.edu
@author: William L. Jorgensen Lab 

REQUIREMENTS:
BOSS (need to set BOSSdir in bashrc and cshrc)
Preferably Anaconda python with following modules
pandas 
argparse
numpy
"""

from LigParGen.BOSSReader import bossPdbAtom2Element,bossElement2Mass,ucomb,tor_cent
import pickle
import pandas as pd
import numpy as np


def retDihedImp(df):
    odihed = []
    if np.sum([df['V' + str(pot)] for pot in range(1, 5)]) != 0.0:
        for pot in range(1, 5):
            if (df['V' + str(pot)] != 0.0):
                odihed.append('IMPRoper %s %4.5f %d %4.5f \n' % (df['NAME'].replace(
                    "-", "   "), df['V' + str(pot)], pot, 180.00 * abs(pot % 2 - 1)))
    else:
        pot = 2
        odihed.append('IMPRoper %s %4.5f %d %4.5f \n' % (df['NAME'].replace(
            "-", "   "), df['V' + str(pot)], pot, 180.00 * abs(pot % 2 - 1)))
    return (odihed)


def retDihed(df):
    odihed = []
    odihed.append('DIHEdral %s  MULT 4   %4.3f     %d   %4.2f \n' % (df['NAME'].replace("-", "   "), df['V' + str(1)], 1, 180.00 * abs(1 % 2 - 1)))
    for pot in range(2, 5):
        odihed.append('                                             %4.3f     %d   %4.2f \n' % (df['V' + str(pot)], pot, 180.00 * abs(pot % 2 - 1)))
    return (odihed)


def Boss2CharmmRTF(num2typ2symb, Qs, resid, bnd_df, angs,props,imps):
    charges = [float(Qs[i][1]) for i in range(len(Qs))]
    rtf = open(resid + '.top', 'w+')
    rtf.write('Remarks generated XPLOR-TOP file for Bonvin group (by Leela Dodda)\n')
    rtf.write('\nset echo=false end\n')
    rtf.write('\nautogenerate angles=True dihedrals=True end\n')
    rtf.write('{ atomType  mass }\n')
    Mass = ['MASS %s %3.4f \n' % (num2typ2symb[i][2], bossElement2Mass(bossPdbAtom2Element(num2typ2symb[i][0]))) for i in range(len(Qs))]
    for i in range(len(Mass)):
        rtf.write('%s' % Mass[i])
    rtf.write('\nRESIdue %5s\n' % (resid))
    rtf.write('\nGROUP\n')
    rtf.write('\n{ atomName  atomType  Charge } \n')
    for i in range(len(Qs)):
        rtf.write('ATOM %6s TYPE= %6s CHARGE= %8s END\n' % (num2typ2symb[i][0], num2typ2symb[i][2], Qs[i][1]))
    rtf.write('\n{ Bonds: atomName1  atomName2 } \n')
    for (x, y) in zip(bnd_df.cl1, bnd_df.cl2):
        rtf.write('BOND %s %s \n' % (num2typ2symb[x][0], num2typ2symb[y][0]))
    bndlist = list(bnd_df.UR) + list(bnd_df.UR)
    imp_list = []
    for i,dat in imps.iterrows():
        ndata =tor_cent([dat.I, dat.J, dat.K, dat.L], bndlist) 
        sdata = [num2typ2symb[j][0] for j in ndata]
        imp_list.append('-'.join(sdata))
    rtf.write('\n{ Improper Dihedrals: aName1 aName2 aName3 aName4 }\n')
    for i in imp_list:
        rtf.write('IMPRoper %s \n' % (i.replace("-", " ")))
    rtf.write('\nEND {RESIdue UNK}\n')
    rtf.write('\nset echo=true end\n')
    rtf.close()
    return None


def Boss2CharmmPRM(resid, num2typ2symb, Qs, bnd_df, ang_df, tor_df):
    #### COLLECTING NONBONDING PART #######
    prm = open(resid + '.param', 'w+')
    prm.write('Remarks generated XPLOR-TOP file for Bonvin group (by Leela Dodda)\n')
    prm.write('\nset echo=false end \n')
    prm.write('\n{ Bonds: atomType1 atomType2 kb r0 } \n')
    for i in bnd_df.index:
        prm.write('BOND %6s %6s %8.1f %8.4f \n' % (num2typ2symb[bnd_df.cl1[i]][
            2], num2typ2symb[bnd_df.cl2[i]][2], bnd_df.KIJ[i], bnd_df.RIJ[i]))
    prm.write('\n{ Angles: aType1 aType2 aType3 kt t0 }\n')
    for i in ang_df.index:
        prm.write('ANGLe %5s %5s %5s %8.1f %8.2f \n' % (num2typ2symb[ang_df.cl1[i]][2], num2typ2symb[
            ang_df.cl2[i]][2], num2typ2symb[ang_df.cl3[i]][2], ang_df.K[i], ang_df.R[i]))
    prm.write('\n{ Proper Dihedrals: aType1 aType2 aType3 aType4 kt period phase } \n')
    if len(tor_df.index) > 0:
        tor_df = tor_df.drop_duplicates(['NAME', 'TY'])
    pro_df = tor_df[tor_df.TY == 'Proper']
    for i in list(pro_df.index):
        ndf = pro_df.iloc[i]
        pro_out = retDihed(ndf.to_dict())
        for i in range(4):
            prm.write('%s' % pro_out[i])
    prm.write(
        'DIHEdral    X      X      X      X  MULT 1   0.000     1     0.00 ! WILD CARD FOR MISSING TORSION PARAMETERS\n')
    prm.write('\n{ Improper Dihedrals: aType1 aType2 aType3 aType4 kt period phase }\n')
    imp_df = tor_df[tor_df.TY == 'Improper']
    for i in list(imp_df.index):
        ndf = tor_df.iloc[i]
        imp_out = retDihedImp(ndf.to_dict())
        for i in range(len(imp_out)):
            prm.write('%s' % imp_out[i])
    prm.write(
        'IMPRoper    X      X      X      X  0.000     1     0.00  ! WILD CARD FOR MISSING IMPROPER PARAMETERS \n')
    prm.write(
        '\n{ Nonbonded: Type Emin sigma; (1-4): Emin/2 sigma }\n')
    Qlines = ['NONBonded %5s %11.6f %11.6f %11.6f %11.6f \n' %
              (num2typ2symb[i][2], float(Qs[i][3]), float(Qs[i][2]), float(Qs[i][3]) * 0.50,float(Qs[i][2])) for i in range(len(Qs))]
    for i in range(len(Qlines)):
        prm.write('%s' % Qlines[i])
    prm.write('\nset echo=true end \n')
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
#    bnd_df.to_csv('bos_bonds.csv', index=False)
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
#    ang_df.to_csv('bos_angles.csv', index=False)
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
        trm = bossPdbAtom2Element(num2typ2symb[i][0]) + num2typ2symb[i][1][-3:]
        num2typ2symb[i].append(trm[0:4])
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
        tor_df[tor_df.TY == 'Proper']['SYMB']),tor_df[tor_df.TY == 'Improper'])
    Boss2CharmmPRM(resid, num2typ2symb, Qs, bnd_df, ang_df, tor_df)
    return None


def mainBOSS2XPLOR(resid, clu=False):
    mol = pickle.load(open(resid + ".p", "rb"))
    Boss2Charmm(resid, mol)
    return None
