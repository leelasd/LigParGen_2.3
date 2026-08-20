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

from LigParGen.BOSSReader import ucomb
from LigParGen.boss_common import bossData, pair_declared_torsions, translate_zmat_indices
import pickle
import pandas as pd
import numpy as np


def retDihedImp(df):
    # Q's own periodic-improper energy term (bondene.f90's improper2)
    # hardcodes the multiplicity to 2 (arg = 2*phi - imp0), and its
    # [impropers] parameter line has exactly TWO numeric fields -- force
    # constant and phase -- with no periodicity column (confirmed directly
    # against Q6's own prep.f90 reader: `read(line, *, ...) taci, tacj,
    # tack, tacl, imp_prm(i)%prm`, where imp_prm(i)%prm's type -- fk, imp0
    # -- has only those two real components). Writing a third (periodicity)
    # field here, as the CHARMM/other-format writers do, silently misaligns
    # columns: Q's list-directed read binds our periodicity placeholder to
    # imp0 (the phase) and never reads our real phase value at all --
    # confirmed directly (produced ~60 kcal/mol of spurious improper energy
    # for benzene's perfectly planar ring, where BOSS's own energy is
    # exactly 0). This also means Q can only represent OPLS-AA's n=2
    # improper term (V2) -- the only term BOSS's own impropers ever
    # populate in practice.
    return ['%s %8.5f %4.5f \n' % (df['NAME'].replace("-", "   "), 2.0*df['V2'], 180.00)]


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
        rtf.write('%8d %6s %6s %8s \n' % (
            i+1, num2typ2symb[i][0], num2typ2symb[i][2], Qs[i][1]))
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
    # Q's own parameter reader (Qprep6's readprm) hard-requires vdw_rule in
    # [options] and rejects the WHOLE file without it ("vdw_rule in options
    # section not found") -- confirmed directly, an empty [options] section
    # silently drops every bond/angle/torsion/atom_type below it too, not
    # just vdW. scale_14/switch_atoms/improper_potential/
    # improper_definition match the real Qoplsaa.prm reference file
    # bundled with Q6's own test suite. improper_definition explicit is
    # required for OPLS-AA specifically: without it Qprep6 auto-generates
    # its own GROMOS-style improper for every sp2/3-connected ring atom,
    # double-counting the planarity restraint OPLS-AA already bakes into
    # the *proper* torsion Fourier series for those same ring atoms.
    prm.write('\n[options]\n')
    prm.write('name Q-OPLSAA\n')
    prm.write('type AMBER\n')
    prm.write('vdw_rule geometric\n')
    prm.write('scale_14 0.5\n')
    prm.write('switch_atoms on\n')
    prm.write('improper_potential periodic\n')
    prm.write('improper_definition explicit\n')
    prm.write('\n[atom_types]\n')
    # Q's [atom_types] reader (Qprep6) rejects a repeated type NAME
    # outright ("Could not enumerate atom type ... Duplicate name?"), and
    # that rejection corrupts every atom's type-index assignment for the
    # rest of the topology build ("Inconsistent molecule/residue start
    # atoms", every bond/angle/torsion count coming out 0) -- confirmed
    # directly. Unlike the CHARMM/TINKER/LAMMPS writers, which give every
    # atom its own row (fine there -- those readers don't reject
    # duplicates), Q needs exactly one declaration per unique OPLS type.
    seen_types = set()
    for i in range(len(Qs)):
        typename = num2typ2symb[i][2]
        if typename in seen_types:
            continue
        seen_types.add(typename)
        eps = float(Qs[i][3])
        sig = float(Qs[i][2])
        ALJ = 2*sig**6*np.sqrt(eps)
        BLJ = 2*sig**3*np.sqrt(eps)
        # Q's own [atom_types] row format (confirmed against the real
        # Qoplsaa.prm reference file and Q6's own prep.f90 parser) is
        # Avdw1 Avdw2 Bvdw1 Avdw3 Bvdw2&3 mass -- Avdw3/Bvdw2&3 are the
        # SEPARATE 1-4-scaled LJ A/B values Q uses for 1-4 pairs (its own
        # `precompute_set_values_pp` combines them by direct multiplication,
        # not sqrt, so storing them pre-scaled by sqrt(0.5) here makes the
        # combined pairwise A_14/B_14 come out scaled by the correct 0.5 --
        # matching OPLS-AA's 1-4 LJ scaling, the same 0.5 already used for
        # electrostatics via scale_14 above). Writing BLJ/0.000 here
        # (repeating the normal B value into the A_14 slot, zeroing B_14)
        # made every 1-4 LJ interaction wrong -- confirmed directly: total
        # vdW energy came out negative instead of matching BOSS.
        half_sqrt = 0.7071067811865476  # sqrt(0.5)
        prm.write('%4s %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n'%(typename,ALJ,ALJ,BLJ,half_sqrt*ALJ,half_sqrt*BLJ,num2typ2symb[i][4]))
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
    # index labels here are whatever survived drop_duplicates()/the TY
    # filter, not a contiguous 0..len-1 range, so looking a row up by
    # .iloc[i] (positional) instead of .loc[i] (by that label) goes out of
    # bounds as soon as any earlier row was filtered out (pre-existing bug,
    # previously masked by bossData() crashing before execution reached
    # here).
    for i in list(pro_df.index):
        ndf = pro_df.loc[i]
        pro_out = retDihed(ndf.to_dict())
        for i in range(4):
            prm.write('%s' % pro_out[i])
    prm.write(
        '!   X    X    X    X    0.00000 1 0.000000 ! WILD CARD FOR MISSING TORSION PARAMETERS\n')
    prm.write('\n[impropers]\n')
    imp_df = tor_df[tor_df.TY == 'Improper']
    for i in list(imp_df.index):
        ndf = tor_df.loc[i]
        imp_out = retDihedImp(ndf.to_dict())
        for i in range(len(imp_out)):
            prm.write('%s' % imp_out[i])
    prm.write(
        '!   X    X    X    X    0.00000 2 0.000000 ! WILD CARD FOR MISSING IMPROPER PARAMETERS \n')
    prm.close()
    return None


def Boss2CharmmTorsion(bnd_df, num2opls, zmat_idx_map, molecule_data, num2typ2symb):
    ats = []
    for line in molecule_data.MolData['ATOMS'][3:]:
        dt = [line.split()[0], line.split()[4],
              line.split()[6], line.split()[8]]
        dt = [int(d) for d in dt]
        ats.append(dt)
    for line in molecule_data.MolData['ADD_DIHED']:
        dt = [int(l) for l in line]
        ats.append(dt)

    paired_ats, paired_dhd = pair_declared_torsions(molecule_data, ats)

    if len(paired_dhd) == 0:
        # A molecule with no torsions at all (e.g. water, or a bare
        # monatomic ion) makes paired_dhd/paired_ats empty lists --
        # np.array([]) has shape (0,), which pd.DataFrame(..., columns=[4
        # names]) can't reshape into, so build the (correctly empty)
        # DataFrames directly instead of through the array conversion.
        dhd_df = pd.DataFrame(columns=['V1', 'V2', 'V3', 'V4'])
        at_df = pd.DataFrame(columns=['I', 'J', 'K', 'L'])
    else:
        dhd = np.array(paired_dhd)
        dhd = dhd  # kcal to kj conversion
        dhd = dhd / 2.0  # Komm = Vopls/2
        dhd_df = pd.DataFrame(dhd, columns=['V1', 'V2', 'V3', 'V4'])
        ats = np.array([translate_zmat_indices(row, zmat_idx_map) for row in paired_ats])
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


def boss2CharmmBond(molecule_data, zmat_idx_map):
    bdat = molecule_data.MolData['BONDS']
    bdat['cl1'] = translate_zmat_indices(bdat['cl1'], zmat_idx_map)
    bdat['cl2'] = translate_zmat_indices(bdat['cl2'], zmat_idx_map)
    bnd_df = pd.DataFrame(bdat)
    bnd_df['UF'] = ((bnd_df.cl1 + bnd_df.cl2) *
                    (bnd_df.cl1 + bnd_df.cl2 + 1) * 0.5) + bnd_df.cl2
    bnd_df['UR'] = ((bnd_df.cl1 + bnd_df.cl2) *
                    (bnd_df.cl1 + bnd_df.cl2 + 1) * 0.5) + bnd_df.cl1
    hb_df = bnd_df.drop(['cl1', 'cl2', 'UF', 'UR'], axis=1)
    hb_df = hb_df.drop_duplicates()
    return bnd_df


def boss2CharmmAngle(anglefile, num2opls, zmat_idx_map,num2typ2symb):
    adat = anglefile
    adat['cl1'] = translate_zmat_indices(adat['cl1'], zmat_idx_map)
    adat['cl2'] = translate_zmat_indices(adat['cl2'], zmat_idx_map)
    adat['cl3'] = translate_zmat_indices(adat['cl3'], zmat_idx_map)
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


def Boss2Charmm(resid, molecule_data):
    types, Qs, num2opls, zmat_idx_map, num2typ2symb, num2pqrtype = bossData(molecule_data)
    bnd_df = boss2CharmmBond(molecule_data, zmat_idx_map)
    ang_df = boss2CharmmAngle(molecule_data.MolData['ANGLES'], num2opls, zmat_idx_map,num2typ2symb)
    tor_df = Boss2CharmmTorsion(bnd_df, num2opls, zmat_idx_map,
                                molecule_data, num2typ2symb)
    Boss2CharmmRTF(num2typ2symb, Qs, resid, bnd_df, list(ang_df['TY']), list(
        tor_df[tor_df.TY == 'Proper']['SYMB']),list(tor_df[tor_df.TY == 'Improper']['SYMB']))
    Boss2CharmmPRM(resid, num2typ2symb, Qs, bnd_df, ang_df, tor_df)
    return None


def mainBOSS2Q(resid, clu=False):
    mol = pickle.load(open(resid + ".pkl", "rb"))
    Boss2Charmm(resid, mol)
    return None
