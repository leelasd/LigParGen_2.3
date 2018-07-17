"""
SCRIPT TO CONVERT WRITE CHARMM RTF AND PRM FILES 
FROM BOSS ZMATRIX
Created on Mon Feb 15 15:40:05 2016
@author: Matthew Robinson matthew.robinson@yale.edu
@author: William L. Jorgensen Lab 

Usage: python OPM_Routines.py -z phenol.z -r PHN
REQUIREMENTS:
BOSS (need to set BOSSdir in bashrc and cshrc)
Preferably Anaconda python with following modules
pandas 
argparse
numpy
"""

from LigParGen.BOSSReader import bossPdbAtom2Element,bossElement2Mass,ucomb,tor_cent
import pickle
import os
import pandas as pd
import numpy as np

ATOM_NUMBER_DICT = {'H': 1, 'He': 2, 'Li': 3, 'Be': 4,
                     'B': 5, 'C': 6, 'N': 7, 'O': 8,
                      'F': 9, 'Ne': 10, 'Na': 11, 'Mg': 12,
                       'Al': 13, 'Si': 14, 'P': 15, 'S': 16,
                        'Cl': 17, 'Ar': 18, 'k': 19, 'Ca': 20,
                         'Sc': 21, 'Ti': 22, 'v': 23, 'Cr': 24,
                          'Mn': 25, 'Fe': 26, 'Co': 27, 'Ni': 28,
                           'Cu': 29, 'Zn': 30, 'Ga': 31, 'Ge': 32,
                            'As': 33, 'Se': 34, 'Br': 35, 'Kr': 36,
                             'Rb': 37, 'Sr': 38, 'Y': 39, 'Zr': 40,
                              'Nb': 41, 'Mo': 42, 'Tc': 43, 'Ru': 45,
                               'Pd': 46, 'Ag': 47, 'Cd': 48, 'In': 49,
                                'Sn': 50, 'Sb': 51, 'Te': 52, 'I': 53,
                                 'Xe': 54, 'Cs': 55, 'Ba': 56, 'La': 57,
                                  'Hf': 72, 'Ta': 73, 'W': 74, 'Re': 75,
                                   'Os': 76, 'Ir': 77, 'Pt': 78, 'Au': 79,
                                    'Hg': 80, 'Tl': 81, 'Pb': 82, 'Bi': 83,
                                     'Po': 84, 'At': 85, 'Rn': 86, 'Fr': 87, 
                                     'Ra': 88, 'Ac': 89}

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

def Boss2Tinker(resid, molecule_data, xyz_dict):
    types, Qs, num2opls, st_no, num2typ2symb = bossData(molecule_data)
    bnd_df = boss2CharmmBond(molecule_data, st_no)
    bndlist = list(bnd_df.UR) + (list(bnd_df.UR))
    ang_df = boss2CharmmAngle(molecule_data.MolData['ANGLES'], num2opls, st_no,num2typ2symb)
    tor_df = Boss2CharmmTorsion(bnd_df, num2opls, st_no,
                                molecule_data, num2typ2symb)

    prm = open('/tmp/'+ resid + '.key', 'w+')
    prm.write(
'''

      ##############################
      ##                          ##
      ##  Force Field Definition  ##
      ##                          ##
      ##############################


forcefield              OPLS-AA

vdwindex                TYPE
vdwtype                 LENNARD-JONES
radiusrule              GEOMETRIC
radiustype              SIGMA
radiussize              DIAMETER
epsilonrule             GEOMETRIC
torsionunit             1.0
imptorunit              1.0
vdw-14-scale            2.0
chg-14-scale            2.0
electric                332.06
dielectric              1.0


      #############################
      ##                         ##
      ##  Atom Type Definitions  ##
      ##                         ##
      #############################


''')

    dict_counter = 1
    for type_list in types:
        type_num = int(type_list[1].strip('_opls'))
        prm.write('atom %10d %4d %5s %8s %10d %10.3f %5d \n' %
            (type_num, type_num, type_list[-1], '"' + type_list[0] + '"', 
                ATOM_NUMBER_DICT[type_list[-3]], type_list[-2], xyz_dict[dict_counter][2]))
        dict_counter += 1
    prm.write(
'''


      ################################
      ##                            ##
      ##  Van der Waals Parameters  ##
      ##                            ##
      ################################


''')
    types_idx = 0
    for vdw_list in Qs:
        type_num = int(types[types_idx][1].strip('_opls'))
        sigma = float(vdw_list[2])
        e_min = float(vdw_list[3])
        prm.write('vdw %11d %16.4f %8.4f \n' % 
                (type_num, sigma, e_min))
        types_idx += 1

    prm.write(
'''


      ##################################
      ##                              ##
      ##  Bond Stretching Parameters  ##
      ##                              ##
      ##################################


''' 
)
    # ask about this one
    for index, row in bnd_df.iterrows():
        atom1_type = int(types[int(row['cl1'])][1].strip('_opls'))
        atom2_type = int(types[int(row['cl2'])][1].strip('_opls'))
        R = float(row['RIJ'])
        K = float(row['KIJ'])

        prm.write('bond %10d %4d %16.2f %8.4f \n' % 
                    (atom1_type, atom2_type, K, R))

    prm.write(
'''


      ################################
      ##                            ##
      ##  Angle Bending Parameters  ##
      ##                            ##
      ################################


''')
    for index, row in ang_df.iterrows():
        atom1_type = int(types[int(row['cl1'])][1].strip('_opls'))
        atom2_type = int(types[int(row['cl2'])][1].strip('_opls'))
        atom3_type = int(types[int(row['cl3'])][1].strip('_opls'))
        R = float(row['R'])
        K = float(row['K'])

        prm.write('angle %9d %4d %4d %8.2f %8.2f \n' % 
                    (atom1_type, atom2_type, atom3_type, K, R))

    prm.write(
'''


      ################################
      ##                            ##
      ##   Urey-Bradley Parameters  ##
      ##                            ##
      ################################


ureybrad     35   34   35      38.25     1.5537


      #####################################
      ##                                 ##
      ##  Improper Torsional Parameters  ##
      ##                                 ##
      #####################################



''')
    ### Impropers ###
    for index, row in tor_df.iterrows():
        if row['TY'] == 'Improper':
            cen_nums = tor_cent([row.I,row.J,row.K,row.L],bndlist)
            atom1_type = int(num2typ2symb[cen_nums[1]][1][5:])#int[int(row['I'])][1].strip('_opls'))
            atom2_type = int(num2typ2symb[cen_nums[2]][1][5:])#int[int(row['I'])][1].strip('_opls'))
            atom3_central_type = int(num2typ2symb[cen_nums[0]][1][5:]) #int(types[int(row['J'])][1].strip('_opls'))
            atom4_type = int(num2typ2symb[cen_nums[3]][1][5:]) 

            V2 = float(row['V2'])
            gamma = 180.0
            n = 2

            # ordering for this is weird
            # see https://ryanmrichard.github.io/ForceManII/tinkerformat.html
            prm.write('imptors %7d %4d %4d %4d %12.3f %4.1f %2d \n' % 
                    (atom1_type, atom2_type, atom3_central_type, atom4_type, V2, gamma, n))

    prm.write(
'''


      ############################
      ##                        ##
      ##  Torsional Parameters  ##
      ##                        ##
      ############################


''')
    for index, row in tor_df.iterrows():
        if row['TY'] == 'Proper':
            atom1_type = int(types[int(row['I'])][1].strip('_opls'))
            atom2_type = int(types[int(row['J'])][1].strip('_opls'))
            atom3_type = int(types[int(row['K'])][1].strip('_opls'))
            atom4_type = int(types[int(row['L'])][1].strip('_opls'))

            V1 = float(row['V1'])
            gamma1 = 0.0
            n1 = 1

            V2 = float(row['V2'])
            gamma2 = 180.0
            n2 = 2

            V3 = float(row['V3'])
            gamma3 = 0.0
            n3 = 3

            prm.write('torsion %7d %4d %4d %4d %12.3f %4.1f %2d %6.3f %4.1f %2d %6.3f %4.1f %2d \n' % 
                    (atom1_type, atom2_type, atom3_type, atom4_type, V1, gamma1, n1, V2, gamma2, n2, V3, gamma3, n3))
    prm.write(
'''
torsion       0    0    0    0        0.000  0.0  1  0.000 180.0  2  0.000  0.0  3

      ########################################
      ##                                    ##
      ##  Atomic Partial Charge Parameters  ##
      ##                                    ##
      ########################################


''')
    types_idx = 0
    for vdw_list in Qs:
        type_num = int(types[types_idx][1].strip('_opls'))
        charge = float(vdw_list[1])
        prm.write('charge %11d %16.4f \n' % 
                (type_num, charge))
        types_idx += 1

    prm.close()

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
    hb_df = bnd_df.drop(['cl1', 'cl2', 'UF', 'UR'], 1)
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

def Boss2CharmmTorsion(bnd_df, num2opls, st_no, molecule_data, num2typ2symb):
    #    print num2opls
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
    final_df = pd.concat([dhd_df, at_df], axis=1, join_axes=[at_df.index])
    bndlist = list(bnd_df.UR) + (list(bnd_df.UR))
    final_df['TY'] = ['Proper' if ucomb(list([final_df.I[n], final_df.J[n], final_df.K[
        n], final_df.L[n]]), bndlist) == 3 else 'Improper' for n in range(len(final_df.I))]
#    final_df['SumV'] = np.abs(
#        final_df.V1) + np.abs(final_df.V2) + np.abs(final_df.V3) + np.abs(final_df.V4)
    #    final_df = final_df[final_df['SumV'] != 0.00]
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



def create_xyz_file(residue_name,mol):
    boss_xyz = mol.MolData['XYZ']
    # convert .pdb to Tinker style .xyz file
    os.system('obabel -ipdb /tmp/%s.pdb -otxyz -O/tmp/%s.xyz >& LLN' % (residue_name,residue_name))

    # Read in the file
    with open('/tmp/%s.xyz' % residue_name, 'r') as xyz_file:
        xyz_data = xyz_file.readlines()
    num_atoms = int(xyz_data[0][0:6])
    xyz_dict = {}
    line_counter = 1
    for line in xyz_data[1:]:
        row=(boss_xyz.iloc[line_counter-1])
        atom_number = int(line[0:6])
        element = str(line[7:9].strip())
        atom_type = int(line[49:53])
        num_bonds = len(line[55:].split())

        xyz_dict[atom_number] = [element, atom_type, num_bonds]

        # change atom type
        new_atom_type_str = ('     ' + str(799+atom_number))[-4:]
#        xyz_data[line_counter] = line[:49] + new_atom_type_str + line[53:]
        xyz_data[line_counter] = line[:12] + '%11.6f %11.6f %11.6f  '%(row.X,row.Y,row.Z) + new_atom_type_str + line[53:]
        line_counter += 1
    xyz_data[0] = '%6d %s LigParGen generated OPLS-AA/CM1A Parameters\n'%(num_atoms,residue_name)
    with open('/tmp/%s.new.xyz' % residue_name, 'w') as new_xyz_file:
        new_xyz_file.writelines(xyz_data)

    # print(xyz_dict)

    return xyz_dict



def mainBOSS2TINKER(resid, clu=False):
    mol = pickle.load(open(resid + ".p", "rb"))
    # if clu:
    #     pdb_file = '/tmp/clu.pdb'
    # else:
    #     pdb_file = '/tmp/plt.pdb'
    xyz_dict = create_xyz_file(resid,mol)
    Boss2Tinker(resid, mol, xyz_dict)
    return None


