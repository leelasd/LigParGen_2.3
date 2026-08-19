"""
SCRIPT TO WRITE OPENMM PDB AND XML FILES 
FROM BOSS ZMATRIX
Created on Mon Feb 15 15:40:05 2016
@author: Lela S. Dodda leela.dodda@yale.edu
@author: William L. Jorgensen Lab 

Usage: python OPM_Routines.py -z phenol.z -r PHN
REQUIREMENTS:
BOSS (need to set BOSSdir in bashrc and cshrc)
Preferably Anaconda python with following modules
pandas 
argparse
numpy
"""

from collections import OrderedDict
from LigParGen.BOSSReader import Refine_PDB_file,get_coos_from_pdb
from LigParGen.BOSSReader import ucomb,bossPdbAtom2Element,bossElement2Mass,tor_cent
from rdkit import Chem
import pickle
import pandas as pd
import numpy as np

_PERIODIC_TABLE = Chem.GetPeriodicTable()


def printDihed(tdat):
    return (
        '<%s class1=\"%s\" class2=\"%s\" class3=\"%s\" class4=\"%s\" k1=\"%6.6f\" k2=\"%6.6f\" k3=\"%6.6f\" k4=\"%6.6f\" periodicity1=\"1\" periodicity2=\"2\" periodicity3=\"3\" periodicity4=\"4\" phase1=\"0.00\" phase2=\"3.141592653589793\" phase3=\"0.00\" phase4=\"3.141592653589793\"/>' % (
            tdat[-1], tdat[0], tdat[1], tdat[2], tdat[3], tdat[4], tdat[5], tdat[6], tdat[7]))


def boss2opmAtom(num2typ2symb, xmlf):
    xmlf.write('<AtomTypes>\n')
    pr_opmAtom = ['<Type name=\"%s\" class=\"%s\" element=\"%s\" mass=\"%6.6f\" />\n' %
                  (num2typ2symb[n][1], num2typ2symb[n][2], num2typ2symb[n][3], num2typ2symb[n][4]) for n in
                  num2typ2symb.keys()]
    for n in list(set(pr_opmAtom)):
        xmlf.write('%s' % n)
    xmlf.write('</AtomTypes>\n')
    return None


def boss2opmTorsion(bnd_df, num2opls, st_no, molecule_data, xmlf):
    # BOSS's own "Angle" column in the Fourier Coefficients table is a
    # 1-based index into the full declared-dihedral sequence (Variable
    # Dihedrals follow, then Additional Dihedrals follow, in that order --
    # see BOSSReader.get_ImpDat's TORinit/TORfinal slice). That table is
    # NOT guaranteed one row per declared quadruple: BOSS can skip a
    # quadruple it can't match to a known torsion type (confirmed: this
    # happens for some of BOSS's own re-derived Additional Dihedrals once
    # a molecule's ring bonds/angles are completed). Pairing by this
    # declared-order index -- MolData['TORSIONS_BY_DECL_IDX'] -- instead
    # of assuming row N always means declared quadruple N is what keeps
    # this correct when rows are skipped, rather than either mispairing
    # silently or asserting the (previously exact-count-required) match.
    ats = []
    for line in molecule_data.MolData['ATOMS'][3:]:
        dt = [line.split()[0], line.split()[4],
              line.split()[6], line.split()[8]]
        dt = [int(d) for d in dt]
        ats.append(dt)
    for line in molecule_data.MolData['ADD_DIHED']:
        dt = [int(l) for l in line]
        ats.append(dt)

    tors_by_idx = molecule_data.MolData['TORSIONS_BY_DECL_IDX']
    paired_ats, paired_dhd = [], []
    for decl_idx, quad in enumerate(ats, start=1):
        coeffs = tors_by_idx.get(decl_idx)
        if coeffs is None:
            continue
        paired_ats.append(quad)
        paired_dhd.append([float(v) for v in coeffs])

    dhd = np.array(paired_dhd)
    dhd = dhd * 4.184  # kcal to kj conversion
    dhd = dhd / 2.0  # Komm = Vopls/2
    dhd_df = pd.DataFrame(dhd, columns=['V1', 'V2', 'V3', 'V4'])
    ats = np.array(paired_ats) - st_no
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
    if len(final_df.index) >= 1:
        final_df['TI'] = [num2opls[j] for j in final_df.I]
        final_df['TJ'] = [num2opls[j] for j in final_df.J]
        final_df['TK'] = [num2opls[j] for j in final_df.K]
        final_df['TL'] = [num2opls[j] for j in final_df.L]
        final_df['NAME'] = final_df.TI + '-' + final_df.TJ + \
            '-' + final_df.TK + '-' + final_df.TL
        final_df = final_df.sort_values(['NAME'])
        tor_bos = final_df.drop(
            ['I', 'J', 'K', 'L', 'TI', 'TJ', 'TK', 'TL'], axis=1)
        tor_bos = tor_bos.drop_duplicates()
        df = final_df.iloc[tor_bos.index][['TI', 'TJ', 'TK', 'TL',
                                         'V1', 'V2', 'V3', 'V4', 'TY', 'I', 'J', 'K', 'L']]
        torlist = []
        for row in df[df.TY == 'Proper'].iterrows():
            index, dat = row
            data = dat.tolist()
            torlist.append(printDihed(data[0:-4]))
        for row in df[df.TY == 'Improper'].iterrows():
            index, dat = row
            ndata = tor_cent([dat.I, dat.J, dat.K, dat.L], bndlist)
            sdata = [num2opls[j] for j in ndata]
            vdata = [dat.V1, dat.V2, dat.V3, dat.V4, dat.TY]
            if len(set([dat.TI, dat.TJ, dat.TK, dat.TL])) > 3:
                torlist.append(printDihed(sdata + vdata))
        xmlf.write('<PeriodicTorsionForce>\n')
        for npt in torlist:
            xmlf.write("%s\n" % npt)
        xmlf.write('</PeriodicTorsionForce>\n')
    return None


def boss2opmBond(num2opls, molecule_data, st_no, xmlf):
    bdat = molecule_data.MolData['BONDS']
    bdat['cl1'] = [x - st_no if not x - st_no < 0 else 0 for x in bdat['cl1']]
    bdat['cl2'] = [x - st_no if not x - st_no < 0 else 0 for x in bdat['cl2']]
    bnd_df = pd.DataFrame(bdat)
    bnd_df['T1'] = [num2opls[x] for x in bnd_df.cl1]
    bnd_df['T2'] = [num2opls[x] for x in bnd_df.cl2]
    bnd_df['KIJ'] = bnd_df['KIJ'] * 836.80
    bnd_df['RIJ'] = bnd_df['RIJ'] * 0.10
    bnd_df['UF'] = ((bnd_df.cl1 + bnd_df.cl2) *
                    (bnd_df.cl1 + bnd_df.cl2 + 1) * 0.5) + bnd_df.cl2
    bnd_df['UR'] = ((bnd_df.cl1 + bnd_df.cl2) *
                    (bnd_df.cl1 + bnd_df.cl2 + 1) * 0.5) + bnd_df.cl1
    connects = []
    for (ai, aj) in zip(bnd_df.cl1, bnd_df.cl2):
        xmlf.write('<Bond from=\"%d\" to=\"%d\"/>\n' % (aj, ai))
        connects.append("CONECT{:>5}{:>5}".format(aj + 1, ai + 1))
    xmlf.write('</Residue>\n')
    xmlf.write('</Residues>\n')
    full_bnd = bnd_df.copy()
    xmlf.write('<HarmonicBondForce>\n')
    for i in bnd_df.index:
        xmlf.write('<Bond class1=\"%s\" class2=\"%s\" length=\"%6.6f\" k=\"%6.6f\"/>\n' %
                   (bnd_df.iloc[i]['T1'], bnd_df.iloc[i]['T2'], bnd_df.iloc[i]['RIJ'], bnd_df.iloc[i]['KIJ']))
    xmlf.write('</HarmonicBondForce>\n')
    return full_bnd, connects


def boss2opmAngle(anglefile, num2opls, st_no, xmlf):
    adat = anglefile
    adat['cl1'] = [x - st_no if not x - st_no < 0 else 0 for x in adat['cl1']]
    adat['cl2'] = [x - st_no if not x - st_no < 0 else 0 for x in adat['cl2']]
    adat['cl3'] = [x - st_no if not x - st_no < 0 else 0 for x in adat['cl3']]
    ang_df = pd.DataFrame(adat)
    ang_df = ang_df[ang_df.K > 0]
    ang_df['K'] = 8.3680 * ang_df['K']
    ang_df['R'] = (np.pi / 180.0) * ang_df['R']
    ang_df['TY'] = np.array([num2opls[i] + '-' + num2opls[j] + '-' + num2opls[k]
                             for i, j, k in zip(ang_df.cl1, ang_df.cl2, ang_df.cl3)])
    ang_df = ang_df.drop_duplicates(['TY'])
    xmlf.write('<HarmonicAngleForce>\n')
    for i in range(0, len(ang_df.index)):
        xmlf.write('<Angle class1=\"%s\" class2=\"%s\" class3=\"%s\" angle=\"%6.6f\" k=\"%6.6f\"/>\n' % (
            list(ang_df.TY)[i].split('-')[
                0], list(ang_df.TY)[i].split('-')[1], list(ang_df['TY'])[i].split('-')[2], list(ang_df['R'])[i],
            list(ang_df['K'])[i]))
    xmlf.write('</HarmonicAngleForce>\n')
    return None


def bossData(molecule_data):
    ats_file = molecule_data.MolData['ATOMS']
    # Elements come from MolData['XYZ']'s own atomic-number column, not
    # bossPdbAtom2Element's name-based guess (strip trailing char, drop
    # digits): that heuristic assumes LigParGen's own "C00"/"H0A"-style
    # generated atom names and mistypes real atom names this package
    # doesn't generate itself, e.g. an amino-acid-style "CH1" comes out as
    # element "CH" and crashes bossElement2Mass. The atomic number BOSS
    # itself reports is unambiguous regardless of naming convention.
    xyz = molecule_data.MolData['XYZ']
    types = []
    for i in enumerate(ats_file):
        types.append([i[1].split()[1], 'opls_' + i[1].split()[2]])
    # st_no offsets a bonded-pair's raw Zmat atom index down to a 0-based
    # index into types/Qs/num2opls (which start at the first REAL atom).
    # LigParGen's own auto-generated Zmats always place exactly 2 leading
    # dummy atoms, so the first real atom's raw index is always 3 -- but
    # BOSS's own reference Zmat library is inconsistent about this (some
    # files use 2 leading dummies, some 3, some place the dummies after
    # the first real atom instead of before), so hardcoding 3 silently
    # mis-indexes bonded pairs for anything that isn't the 2-leading-dummy
    # case, and can even index out of range (reproduced on his.z).
    # Reading it from the first real atom's own raw index is correct for
    # both conventions.
    st_no = int(ats_file[0].split()[0])
    Qs = molecule_data.MolData['Q_LJ']
    assert len(Qs) == len(types), 'Please check the at_info and Q_LJ_dat files'
    assert len(xyz) == len(types), 'Please check the at_info and XYZ data'
    num2typ2symb = {i: types[i] for i in range(len(Qs))}
    for i in range(len(Qs)):
        elem = _PERIODIC_TABLE.GetElementSymbol(int(xyz['at_num'][i]))
        num2typ2symb[i].append(elem + num2typ2symb[i][1][-3:])
        num2typ2symb[i].append(elem)
        num2typ2symb[i].append(bossElement2Mass(elem))
        num2typ2symb[i].append(Qs[i][0])
    num2opls = {}
    for i in num2typ2symb.keys():
        num2opls[i] = num2typ2symb[i][2]
    num2pqrtype = OrderedDict(num2typ2symb)
    for i in range(len(Qs)):
        num2pqrtype[i].append(Qs[i][1])
        num2pqrtype[i].append(Qs[i][2])
    return (types, Qs, num2opls, st_no, num2typ2symb, num2pqrtype)


def pdb_prep(atoms, coos, resid, connects, elems=None):
    # Columns 77-78 (element symbol) are optional per the PDB spec but are
    # what OpenBabel/3Dmol.js/RDKit rely on to type each atom when re-reading
    # this file back -- without them, downstream readers fall back to
    # guessing from the atom name and can mistype atoms entirely (confirmed
    # live: OpenBabel warns on every atom when BOSS2TINKER re-reads this
    # file, and the Space's 3D preview render was degraded the same way).
    # elems, when given, is the atomic-number-derived element per atom
    # (aligned with atoms/coos) computed in bossData() -- passed through
    # here rather than re-derived from the atom name via
    # bossPdbAtom2Element, which mistypes non-LigParGen-style atom names
    # (e.g. an amino-acid-style "CH1").
    opdb = open(resid + '.pdb', 'w+')
    opdb.write('REMARK LIGPARGEN GENERATED PDB FILE\n')
    num = 0
    for idx, (i, j) in enumerate(zip(atoms, coos)):
        num += 1
        elem = elems[idx] if elems is not None else bossPdbAtom2Element(i)
        opdb.write('%-6s%5d %4s %3s  %4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s\n' %
                   ('ATOM', num, i, resid, 1, j[0], j[1], j[2], 1.00, 0.00, elem))
    opdb.write('TER \n')
    for i in range(len(connects)):
        opdb.write('%s \n' % connects[i])
    opdb.write('END')
    opdb.close()
    return None


def pqr_prep(atoms, coos, resid, connects, num2pqrtype):
    at2pqrtype = {num2pqrtype[i][0]: num2pqrtype[i][-2:]
                  for i in num2pqrtype.keys()}
    opdb = open(resid + '.pqr', 'w+')
    opdb.write('REMARK LIGPARGEN GENERATED PQR FILE\n')
    num = 0
    for (i, j) in zip(atoms, coos):
        num += 1
        opdb.write('%-6s%5d %4s %3s  %4d    %8.3f%8.3f%8.3f%8.4f%7.4f\n' %
                   ('ATOM', num, i, resid, 1, j[0], j[1], j[2], float(at2pqrtype[i][0]), float(at2pqrtype[i][1]) * 0.561231))
    opdb.write('TER \n')
    for i in range(len(connects)):
        opdb.write('%s \n' % connects[i])
    opdb.write('END')
    opdb.close()
    return None


def boss2opm(resid, molecule_data, pdb_file):
    xmlf = open(resid + '.xml', 'w+')
    types, Qs, num2opls, st_no, num2typ2symb, num2pqrtype = bossData(
        molecule_data)
    #### COLLECTING NONBONDING PART #######
    nb_part = []
    res_at = []
    for i in range(len(types)):
        nb_part.append('<Atom type=\"%s\" charge=\"%3.6f\" sigma=\"%3.6f\" epsilon=\"%3.6f\" />\n' %
                       (types[i][1], float(Qs[i][1]), float(Qs[i][2]) / 10.0, float(Qs[i][3]) * 4.184))
        res_at.append('<Atom name=\"%s\" type=\"%s\" />\n' %
                      (types[i][0], types[i][1]))
    xmlf.write('<ForceField>\n')
    boss2opmAtom(num2typ2symb, xmlf)
    #### PRINTING ATOM PART #######
    xmlf.write('<Residues>\n')
    xmlf.write('<Residue name=\"%s\">\n' % resid)
    for at in res_at:
        xmlf.write("%s" % at)
    bnd_df, connects = boss2opmBond(num2opls, molecule_data, st_no, xmlf)
    # PRINTING ANGLES AND TORSIONS
    boss2opmAngle(molecule_data.MolData['ANGLES'], num2opls, st_no, xmlf)
    boss2opmTorsion(bnd_df, num2opls, st_no, molecule_data, xmlf)
    #### PRINTING NONBONDING PART #######
    nnb_part = list(set(nb_part))
    xmlf.write('<NonbondedForce coulomb14scale="0.5" lj14scale="0.5">\n')
    for nb in nnb_part:
        xmlf.write("%s" % nb)
    xmlf.write('</NonbondedForce>\n')
    xmlf.write('</ForceField>\n')
    xmlf.close()
    pdblines = Refine_PDB_file(pdb_file)
    atoms, coos = get_coos_from_pdb(pdblines)
    elems = [num2typ2symb[i][3] for i in range(len(atoms))] if len(atoms) == len(num2typ2symb) else None
    pdb_prep(atoms, coos, resid, connects, elems)
    pqr_prep(atoms, coos, resid, connects, num2pqrtype)
    return None


def mainBOSS2OPM(resid, clu):
    mol = pickle.load(open(resid + ".pkl", "rb"))
    if clu:
        pdb_file = '/tmp/clu.pdb'
    else:
        pdb_file = '/tmp/plt.pdb'
    boss2opm(resid, mol, pdb_file)
    return None
