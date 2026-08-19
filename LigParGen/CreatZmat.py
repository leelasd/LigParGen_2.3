#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
AutoZmat_VersionLSD:
A python program to create BOSS zmatrix from any molecular input format.
Need BOSS and OpenBabel executable to work
Python Modeules Needed - networkx, numpy, pandas
Created on Wed Jun 14 2017

@author: Leela Sriram Dodda
@email:  leela.dodda@yale.edu
"""
import LigParGen
import subprocess
import os
import shutil
import numpy as np
from LigParGen.Vector_algebra import pairing_func, angle, dihedral, tor_id, ang_id,bossElement2Num, Distance
import itertools
import collections
import networkx as nx
from openbabel import openbabel as ob
from openbabel import pybel


def _babel_gen3d(ifile, iform):
    # Equivalent to `babel -i<fmt> <ifile> -omol <out>.mol --gen3D`: OBBuilder +
    # MMFF94 minimization + weighted-rotor conformer search, same OBOp the CLI runs.
    # pybel.Molecule.make3D() is NOT equivalent (skips the conformer search).
    mol = pybel.readstring(iform[1], open(ifile).read())
    gen3d = ob.OBOp.FindType("Gen3D")
    gen3d.Do(mol.OBMol, "")
    mol.write("mol", "%s.mol" % iform[0], overwrite=True)


def AsitIsZmat(ifile,optim,resid):
    iform = ifile.split('.')
    # CREATE A MOL FILE FROM ANY FILE
    if iform[1] == 'smi':
        _babel_gen3d(ifile, iform)
    else:
        # Equivalent to `babel -i<fmt> <ifile> -omol <out>.mol --errorlevel 1 -b`:
        # critical-errors-only logging + dative-bond normalization.
        ob.obErrorLog.SetOutputLevel(ob.obError)
        conv = ob.OBConversion()
        conv.SetInAndOutFormats(iform[1], "mol")
        obmol = ob.OBMol()
        conv.ReadFile(obmol, ifile)
        dative = ob.OBOp.FindType("b")
        if dative is not None:
            dative.Do(obmol, "")
        conv.WriteFile(obmol, "%s.mol" % iform[0])
    mollines = open(iform[0] + '.mol', 'r').readlines()
    COOS, ATYPES, MolBonds = ReadMolFile(mollines)
    G_mol, mol_icords = make_graphs(ATYPES, COOS, MolBonds)
    print_ZMAT(ATYPES, G_mol, mol_icords, COOS, '%s.z' % resid, resid)
    return None

def CanonicaliedZmat(ifile,optim,resid):
    iform = ifile.split('.')
    # CREATE A MOL FILE FROM ANY FILE
    if iform[1] == 'smi':
        _babel_gen3d(ifile, iform)
    else:
        # Equivalent to `babel -i<fmt> <ifile> -omol --canonical <out>.mol`.
        conv = ob.OBConversion()
        conv.SetInAndOutFormats(iform[1], "mol")
        obmol = ob.OBMol()
        conv.ReadFile(obmol, ifile)
        canon = ob.OBOp.FindType("canonical")
        canon.Do(obmol, "")
        conv.WriteFile(obmol, "%s.mol" % iform[0])
    mollines = open(iform[0] + '.mol', 'r').readlines()
    COOS, ATYPES, MolBonds = ReadMolFile(mollines)
    G_mol, mol_icords = make_graphs(ATYPES, COOS, MolBonds)
    print_ZMAT(ATYPES, G_mol, mol_icords, COOS, '%s.z' % resid, resid)
    return None

def GenMolRep(ifile,optim,resid,charge):
    iform = ifile.split('.')
    try: 
        AsitIsZmat(ifile,optim,resid) 
    except (ZeroDivisionError,IndexError):
        print('Warning!!\n 1.Cannonicalising Input MOL/PDB file\n 2.Atom ordering may change \n 3.But the Coordinates remain the same')
        CanonicaliedZmat(ifile,optim,resid)
    Get_OPT('%s.z' % resid, optim, charge)
    if os.path.exists('/tmp/clu.pdb'): os.remove('/tmp/clu.pdb')
    if iform[1] == 'pdb':
        if os.environ.get('MCPROdir') is not None:
            os.system('$MCPROdir/miscexec/clu -t:f=pdb %s.pdb -r %s.z -n:f=p /tmp/clu.pdb -m ma' % (iform[0], resid))
        else:
            execfile = os.environ['BOSSdir'] + '/scripts/xSPM > /tmp/olog'
            coma = execfile + ' ' + resid
            os.system(coma)
            shutil.copyfile('/tmp/plt.pdb', '/tmp/clu.pdb')
    return(True)

def Get_OPT(zmat, optim, charge):
    assert os.path.isfile(
        zmat), 'File named %10s does not exist' % zmat
    assert 'BOSSdir' in os.environ, 'Please Make sure $BOSSdir is defined \n xZCM1A and related files are in scripts directory of BOSS'
    execs = {
        2: os.environ['BOSSdir'] + '/scripts/xZCM1A+2 > /tmp/olog',
        -2: os.environ['BOSSdir'] + '/scripts/xZCM1A-2 > /tmp/olog',
        0: os.environ['BOSSdir'] + '/scripts/xZCM1A > /tmp/olog',
        1: os.environ['BOSSdir'] + '/scripts/xZCM1A+  > /tmp/olog',
        -1: os.environ['BOSSdir'] + '/scripts/xZCM1A-  > /tmp/olog',
    }
    print('MOLECULE HAS A CHARGE of %d' % charge)
    execfile = execs[charge]
    coma = execfile + ' ' + zmat[:-2]
    os.system(coma)
    shutil.copyfile('sum', zmat)
    execfile = os.environ['BOSSdir'] + '/scripts/xSPM > /tmp/olog'
    coma = execfile + ' ' + zmat[:-2]
    os.system(coma)
    shutil.copyfile('sum', zmat)
    return (None)

def ReadMolFile(mollines):
    [nats, nbonds] = map(int, (mollines[3][0:3],mollines[3][3:6]))
    cooslines = mollines[4:4 + nats]
    coos = {}
    atypes = {}
    for i in range(nats):
        els = cooslines[i].split()
        coos[i + 1] = [float(e) for e in els[0:3]]
        atypes[i + 1] = els[3]
    bondlines = mollines[4 + nats:4 + nats + nbonds]
    bonds = {'BI': [], 'BJ': [], 'RIJ': [], 'UID': []}
    for line in bondlines:
        [bi, bj] = map(int, [line[0:3],line[3:6]])
        bonds['BI'].append(bi)
        bonds['BJ'].append(bj)
        bonds['RIJ'].append(Distance(coos[bi], coos[bj]))
        bonds['UID'].append(pairing_func(bi, bj))
    return (coos, atypes, bonds)


def make_graphs(atoms, coos, bonds):
    G = nx.DiGraph()
    # ADD NODES USING ATOM TYPES AND COORDINATES
    for i in coos.keys():
        G.add_node(i, XYZ=coos[i], elem=atoms[i],
                   atno=bossElement2Num(atoms[i]))
    for (i, j, rij) in zip(bonds['BI'], bonds['BJ'], bonds['RIJ']):
        G.add_edge(i, j, distance=rij)
        G.add_edge(j, i, distance=rij)
    # Enumerate bonds/angles/torsions by walking direct neighbors of each
    # node/edge instead of the previous all_pairs_shortest_path_length +
    # all_simple_paths(cutoff=1/2/3) approach. That O(N^2) pair-enumeration
    # also had a latent correctness gap on small rings (3-/4-/5-membered):
    # classifying a triple/quadruple by the graph-wide shortest-path
    # distance between its endpoints silently dropped angles/torsions
    # whenever a shortcut existed around the other side of the ring (e.g.
    # a bare 3-membered ring produced zero angles, since every pair of
    # atoms in a triangle is mutually 1 bond apart). Walking directly from
    # each node's (for angles) or each edge's (for torsions) neighbors,
    # exactly as expand_zmat.py's validated _complete_internals() does,
    # sidesteps that: it enumerates every angle/torsion that is actually
    # implied by the local connectivity, independent of shortcuts
    # elsewhere in the graph.
    all_bonds = [list(e) for e in G.edges()]
    new_angs = []
    for j in G.nodes():
        neigh = list(G.neighbors(j))
        for a in range(len(neigh)):
            for b in range(len(neigh)):
                if a == b:
                    continue
                new_angs.append([neigh[a], j, neigh[b]])
    new_tors = []
    for j, k in G.edges():
        for i in G.neighbors(j):
            if i == k:
                continue
            for l in G.neighbors(k):
                if l == j or l == i:
                    continue
                new_tors.append([i, j, k, l])
    dict_new_tors = {tor_id(t): t for t in new_tors}
    dict_new_angs = {ang_id(t): t for t in new_angs}
    imp_keys = [n for n in G.nodes() if G.degree(n) / 2 == 3]
    all_imps = {}
    for i in imp_keys:
        nei = list(G.neighbors(i))
        if G.nodes[i]['atno'] == 6:
            all_imps[i] = [nei[0], i, nei[1], nei[2]]
    MOL_ICOORDS = {'BONDS': all_bonds,
                   'ANGLES': dict_new_angs, 'TORSIONS': dict_new_tors, 'IMPROPERS': all_imps}
    return(G, MOL_ICOORDS)


def Get_Add_Int(mol_icords, Z_BONDS, Z_ANGLES, Z_TORSIONS):
    all_bonds_mol, all_angles_mol, all_torsions_mol = mol_icords[
        'BONDS'], mol_icords['ANGLES'], mol_icords['TORSIONS']
    Z_B = {pairing_func(i[0] - 2, i[1] - 2): [i[0] - 2, i[1] - 2]
           for i in Z_BONDS.values()}
    Z_A = {ang_id([i[0] - 2, i[1] - 2, i[2] - 2]): [i[0] - 2,
                                                    i[1] - 2, i[2] - 2] for i in Z_ANGLES.values()}
    Z_T = {tor_id([i[0] - 2, i[1] - 2, i[2] - 2, i[3] - 2]): [i[0] - 2,
                                                              i[1] - 2, i[2] - 2, i[3] - 2] for i in Z_TORSIONS.values()}
    Z_Ad_B, Z_Ad_A, Z_Ad_T = collections.OrderedDict(), collections.OrderedDict(), collections.OrderedDict()
    for b_ij in all_bonds_mol:
        uid_b_ij = pairing_func(b_ij[0], b_ij[1])
        if uid_b_ij not in list(Z_B.keys()):
            Z_Ad_B[uid_b_ij] = [b_ij[0] + 2, b_ij[1] + 2]
    for a_ij in all_angles_mol.keys():
        if a_ij not in list(Z_A.keys()):
            Z_Ad_A[a_ij] = [i + 2 for i in all_angles_mol[a_ij]]
    for t_ij in all_torsions_mol.keys():
        if t_ij not in list(Z_T.keys()):
            Z_Ad_T[t_ij] = [i + 2 for i in all_torsions_mol[t_ij]]
    for c in mol_icords['IMPROPERS'].values():
        Z_Ad_T["-".join(list(map(str, c)))] = [i + 2 for i in c]
    return(Z_Ad_B, Z_Ad_A, Z_Ad_T)


def print_ZMAT(atoms, G_mol, mol_icords, coos, zmat_name, resid):
    if not zmat_name:
        zmat_name = resid
    Z_ATOMS = {1: 'X', 2: 'X'}
    Z_NO = {1: -1, 2: -1}
    Z_BONDS = {1: (1, 0, 0.000), 2: (2, 1, 1.00), 3: (3, 2, 1.00)}
    Z_ANGLES = {1: (1, 0, 0, 0.000), 2: (2, 1, 0, 0.000),
                3: (3, 2, 1, 90.00), 4: (4, 3, 2, 90.0)}
    Z_TORSIONS = {1: (1, 0, 0, 0, 0.00), 2: (2, 1, 0, 0, 0.00), 3: (
        3, 2, 1, 0, 0.00), 4: (4, 3, 2, 1, 0.00), 5: (5, 4, 3, 2, 90.0)}
    for i in range(1, len(atoms) + 1):
        Z_ATOMS[i + 2] = atoms[i]
    for i in range(1, len(atoms) + 1):
        Z_NO[i + 2] = G_mol.nodes[i]['atno']
    n_ats = 0
    B_LINK = {}
    for i in G_mol.nodes():
        if n_ats > 0:
            neigs = np.sort(list(G_mol.neighbors(i)))
            B_LINK[i] = neigs[0]
            Z_BONDS[i + 2] = (i + 2, neigs[0] + 2, G_mol[i]
                              [neigs[0]]['distance'])
        n_ats += 1
    n_ats = 0
    A_LINK = {}
    for i in G_mol.nodes():
        if n_ats > 1:
            neigs = np.sort(list(G_mol.neighbors(B_LINK[i])))
            A_LINK[i] = neigs[0]
            ang = angle(coos[i], coos[B_LINK[i]], coos[neigs[0]])
            Z_ANGLES[i + 2] = (i + 2, B_LINK[i] + 2, neigs[0] + 2, ang)
        n_ats += 1
    n_ats = 0
    for i in G_mol.nodes():
        if n_ats > 2:
            neigs =list(G_mol.neighbors(A_LINK[i]))
            neigs = np.array([j for j in neigs if j not in [i, B_LINK[i], A_LINK[i]]])
            neigs = np.sort(neigs)
            neigs = neigs[neigs<i]
            if len(neigs)<1:
               neigs = [j for j in list(G_mol.neighbors(B_LINK[i])) if j not in [i,A_LINK[i]]]
               if (B_LINK[i] in list(mol_icords['IMPROPERS'].keys())): del mol_icords['IMPROPERS'][B_LINK[i]]
            [ti, tj, tk, tl] = [i, B_LINK[i], A_LINK[i], neigs[0]]
            dihed = dihedral(coos[ti], coos[tj], coos[tk], coos[tl])
            Z_TORSIONS[i + 2] = (ti + 2, tj + 2, tk + 2, tl + 2, dihed)
        n_ats += 1
    Z_Ad_B, Z_Ad_A, Z_Ad_T = Get_Add_Int(
        mol_icords, Z_BONDS, Z_ANGLES, Z_TORSIONS)
# PRINTING ACTUAL Z-MATRIX
    ofile = open(zmat_name, 'w+')
    ofile.write('BOSS Z-Matrix with LSDautozmat (written by Leela S. Dodda)\n')
    for i in range(1, len(atoms) + 3):
        ofile.write('%4d %-3s%5d%5d%5d%12.6f%4d%12.6f%4d%12.6f%4s%5d\n'
                    % (i, Z_ATOMS[i], Z_NO[i], Z_NO[i], Z_BONDS[i][1], Z_BONDS[i][-1], Z_ANGLES[i][-2], Z_ANGLES[i][-1], Z_TORSIONS[i][-2], Z_TORSIONS[i][-1], resid[0:3], 1)
                    )
    ofile.write(
        '''                    Geometry Variations follow    (2I4,F12.6)
                    Variable Bonds follow         (I4)\n'''
                )
    for i in range(4, len(atoms) + 3):
        ofile.write('%4d\n' % i)
    ofile.write('                    Additional Bonds follow       (2I4)\n')
    if len(Z_Ad_B) > 0:
        for i in Z_Ad_B.values():
            ofile.write('%4d%4d\n' % (i[0], i[1]))
    # CREATE A FUNCTION TO DEFINE ADDITIONAL BONDS IN CASE OF RINGS
    ofile.write('''                    Harmonic Constraints follow   (2I4,4F10.4)
                    Variable Bond Angles follow   (I4)\n''')
    for i in range(5, len(atoms) + 3):
        ofile.write('%4d\n' % i)
    ofile.write('                    Additional Bond Angles follow (3I4)\n')
    if len(Z_Ad_A) > 0:
        for i in Z_Ad_A.values():
            ofile.write('%4d%4d%4d\n' % (i[0], i[1], i[2]))
    # CREATE A FUNCTION TO DEFINE ADDITIONAL BONDS IN CASE OF RINGS
    ofile.write(
        '                    Variable Dihedrals follow     (3I4,F12.6)\n')
    for i in range(6, len(atoms) + 3):
        ofile.write('%4d%4d%4d%12.6f\n' % (i, -1, -1, 0.000))
    ofile.write('                    Additional Dihedrals follow   (6I4)\n')
    if len(Z_Ad_T) > 0:
        for k in Z_Ad_T.keys():
            torsion = Z_Ad_T[k]
            ofile.write('%4d%4d%4d%4d%4d%4d\n' %
                        (torsion[0], torsion[1], torsion[2], torsion[3], -1, -1))
    ofile.write(
        '''                    Domain Definitions follow     (4I4)
                    Conformational Search (2I4,2F12.6)
                    Local Heating Residues follow (I4 or I4-I4)
                    Final blank line
''')
    ofile.close()
    return None
