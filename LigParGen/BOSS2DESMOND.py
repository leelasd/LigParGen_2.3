"""
SCRIPT TO WRITE DESMOND CMS FILES 
FROM BOSS ZMATRIX
Created on Mon Nov 26 15:40:05 2017
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


def boss2opmTorsion(bnd_df, num2opls, st_no, molecule_data, itpf):
    dhd = []
    for line in molecule_data.MolData['TORSIONS']:
        dt = [float(f) for f in line]
        dhd.append(dt)
    dhd = np.array(dhd)
    dhd = dhd * 1.  # kcal to kj conversion
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
    bnd_df['KIJ'] = bnd_df['KIJ'] 
    bnd_df['RIJ'] = bnd_df['RIJ'] 
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


def boss2desmond(resid, molecule_data, pdb_file):
    types, Qs, num2opls, st_no, num2typ2symb = bossData(
        molecule_data)
    itpf = open(resid + '.cms', 'w+')
    itpf.write("""{
  s_m_m2io_version
  :::
  2.0.0
}
f_m_ct {
  s_m_title
  r_chorus_box_ax
  r_chorus_box_ay
  r_chorus_box_az
  r_chorus_box_bx
  r_chorus_box_by
  r_chorus_box_bz
  r_chorus_box_cx
  r_chorus_box_cy
  r_chorus_box_cz
  s_ffio_ct_type
  :::
  "Desmond file created by (written by Leela S Dodda)"
                  70.0
                   0.0
                   0.0
                   0.0
                  70.0
                   0.0
                   0.0
                   0.0
                  70.0
  full_system
  m_atom[%d] {
    # First column is atom index #
    i_m_mmod_type
    r_m_x_coord
    r_m_y_coord
    r_m_z_coord
    i_m_residue_number
    s_m_pdb_residue_name
    i_m_atomic_number
    s_m_atom_name
    r_ffio_x_vel
    r_ffio_y_vel
    r_ffio_z_vel
    :::
"""%(len(Qs)))
    #print(molecule_data.MolData['XYZ'])
    xyz = molecule_data.MolData['XYZ']
    for i,r in xyz.iterrows(): 
        itpf.write('    %d        %d %10.8f %10.8f %10.8f      1 \"MOL\"     %d  \"%s\" 0.00000000 0.00000000 0.00000000\n'%(i+1,1,r.X,r.Y,r.Z,r.at_num,r.at_symb))
    bnd_df, connects = boss2gmxBond(molecule_data, st_no, itpf)
    itpf.write('''    :::
  }
  m_bond[%d] {
    i_m_from
    i_m_to
    i_m_order
    i_m_from_rep
    i_m_to_rep
    :::
'''%(len(bnd_df.index)))
    for rb in bnd_df.iterrows():
        i,r = rb
        itpf.write('    %d %d %d %d %d %d\n'%(i+1,r.cl2+1,r.cl1+1,1,1,1))
    itpf.write('''    :::
  }
}
f_m_ct {
  s_m_title
  r_chorus_box_ax
  r_chorus_box_ay
  r_chorus_box_az
  r_chorus_box_bx
  r_chorus_box_by
  r_chorus_box_bz
  r_chorus_box_cx
  r_chorus_box_cy
  r_chorus_box_cz
  s_ffio_ct_type
  :::
  "MOL"
                  70.0
                   0.0
                   0.0
                   0.0
                  70.0
                   0.0
                   0.0
                   0.0
                  70.0
  solute
  m_atom[%d] {
    # First column is atom index #
    i_m_mmod_type
    r_m_x_coord
    r_m_y_coord
    r_m_z_coord
    i_m_residue_number
    s_m_pdb_residue_name
    i_m_atomic_number
    s_m_atom_name
    r_ffio_x_vel
    r_ffio_y_vel
    r_ffio_z_vel
    :::
'''%(len(Qs)))
    xyz = molecule_data.MolData['XYZ']
    for i,r in xyz.iterrows():
        itpf.write('    %d        %d %10.8f %10.8f %10.8f      1 \"MOL\"     %d  \"%s\" 0.00000000 0.00000000 0.00000000\n'%(i+1,1,r.X,r.Y,r.Z,r.at_num,r.at_symb))
#    bnd_df, connects = boss2gmxBond(molecule_data, st_no, itpf)
    itpf.write('''    :::
  }
  m_bond[%d] {
    i_m_from
    i_m_to
    i_m_order
    i_m_from_rep
    i_m_to_rep
    :::
'''%(len(bnd_df.index)))
    for rb in bnd_df.iterrows():
        i,r = rb
        itpf.write('    %d %d %d %d %d %d\n'%(i+1,r.cl2+1,r.cl1+1,1,1,1))
#    for i,r in bnd_df.iterrows():itpf.write('    %d %d %d %d %d %d\n'%(i+1,r.cl2+1,r.cl1+1,1,1,1))
    itpf.write('''    :::
  }
  ffio_ff {
    s_ffio_name
    s_ffio_comb_rule
    i_ffio_version
    :::
    "MOL"
    GEOMETRIC
    1.0.0
    ffio_vdwtypes[%d] {
      s_ffio_name
      s_ffio_funct
      r_ffio_c1
      r_ffio_c2
      :::
'''%(len(Qs)))
    for i in range(len(Qs)): itpf.write('      %d %s LJ12_6_sig_epsilon %10.8f %10.8f\n'%(i+1,num2typ2symb[i][1],float(Qs[i][2]), float(Qs[i][3]))) 
    itpf.write('''      :::
    }
    ffio_sites[%d] {
      s_ffio_type
      r_ffio_charge
      r_ffio_mass
      s_ffio_vdwtype
      i_ffio_resnr
      s_ffio_residue
      :::
'''%(len(Qs)))
    for i in range(len(Qs)): itpf.write('     %d  atom %10.8f %10.8f %s 1  MOL\n'%(i+1,float(Qs[i][1]),num2typ2symb[i][4],num2typ2symb[i][1]))
    itpf.write('''      :::
    }''')
    itpf.write('''    ffio_bonds[%d] {
      i_ffio_ai
      i_ffio_aj
      s_ffio_funct
      r_ffio_c1
      r_ffio_c2
      :::
'''%(len(bnd_df.index)))
    EXCLUSIONS_12_13 = {'BI':[],'BJ':[],'UID':[]} 
    for rb in bnd_df.iterrows():
        i,r = rb
        itpf.write('      %d %d %d HARM      %10.8f    %10.8f\n'%(i+1,r.cl2+1,r.cl1+1,r.RIJ,r.KIJ))
        EXCLUSIONS_12_13['BI'].append(int(r.cl1+1));
        EXCLUSIONS_12_13['BJ'].append(int(r.cl2+1)); 
        EXCLUSIONS_12_13['UID'].append(pairing_func(r.cl2+1,r.cl1+1)[0])
    full_ang = boss2gmxAngle(molecule_data.MolData['ANGLES'], num2opls, st_no, itpf)
    itpf.write('''      :::
    }
    ffio_angles[%d] {
      i_ffio_ai
      i_ffio_aj
      i_ffio_ak
      s_ffio_funct
      r_ffio_c1
      r_ffio_c2
      ::: 
'''%(len(full_ang.index)))
    for rb in full_ang.iterrows():
        i,r = rb
        itpf.write('      %d %d %d %d HARM      %10.8f    %10.8f\n'%(i+1,r.cl1+1,r.cl2+1,r.cl3+1,r.R,r.K))
        EXCLUSIONS_12_13['BI'].append(int(r.cl1+1));
        EXCLUSIONS_12_13['BJ'].append(int(r.cl3+1)); 
        EXCLUSIONS_12_13['UID'].append(pairing_func(r.cl3+1,r.cl1+1)[0])
    full_tor, tor_df = boss2opmTorsion(bnd_df, num2opls, st_no, molecule_data, itpf)
    itpf.write('''      :::
    }
    ffio_dihedrals[%d] {
      i_ffio_ai
      i_ffio_aj
      i_ffio_ak
      i_ffio_al
      s_ffio_funct
      r_ffio_c0
      r_ffio_c1
      r_ffio_c2
      r_ffio_c3
      r_ffio_c4
      :::
'''%(len(full_tor.index)))
    full_tor.index = range(0,len(full_tor.index))
    EXCEPTION_14 = {'BI':[],'BJ':[],'UID':[]}
    for rb in full_tor.iterrows():
        i,r = rb
        if r.TY == 'Proper': 
            itpf.write('      %d %d %d %d %d Opls_proper  0   %10.8f    %10.8f       %10.8f    %10.8f\n'%(i+1,r.I+1,r.J+1,r.K+1,r.L+1,r.V1,r.V2,r.V3,r.V4))
            EXCLUSIONS_12_13['BI'].append(int(r.I+1));
            EXCLUSIONS_12_13['BJ'].append(int(r.L+1)); 
            EXCLUSIONS_12_13['UID'].append(pairing_func(r.I+1,r.L+1)[0])
            EXCEPTION_14['BI'].append(int(r.I+1));
            EXCEPTION_14['BJ'].append(int(r.L+1)); 
            EXCEPTION_14['UID'].append(pairing_func(r.I+1,r.L+1)[0])
        elif r.TY == 'Improper': itpf.write('      %d %d %d %d %d Opls_improper  0    %10.8f    %10.8f       %10.8f    %10.8f\n'%(i+1,r.I+1,r.J+1,r.K+1,r.L+1,r.V1,r.V2,r.V3,r.V4))
    excl = (pd.DataFrame(EXCLUSIONS_12_13))
    excl = excl.drop_duplicates(['UID'])
    excp = (pd.DataFrame(EXCEPTION_14))
    excp = excp.drop_duplicates(['UID'])
    itpf.write('''      :::
    }
    ffio_torsion_torsion[0] {
      :::
      :::
    }
''')
    itpf.write('''    ffio_exclusions[%d] {
      i_ffio_ai
      i_ffio_aj
      :::
'''%(len(excl.index)))
    excl.index = range(len(excl.index))
    for rb in excl.iterrows():
        i,r =rb
        itpf.write('      %d %d %d\n'%(i+1,r.BI,r.BJ))
    itpf.write('''      :::
    }
    ffio_pairs[%d] {
      i_ffio_ai
      i_ffio_aj
      s_ffio_funct
      r_ffio_c1
      r_ffio_c2
      :::
'''%(len(excp.index)*2))
#    prop_tor = full_tor[full_tor.TY=='Proper']
    excp.index = range(len(excp.index))
    for rb in excp.iterrows():
        i,r = rb
        itpf.write('      %d %d %d Coulomb 0.5 <>\n'%(2*i+1,r.BI,r.BJ))
        itpf.write('      %d %d %d LJ 0.5 <>\n'%(2*(i+1),r.BI,r.BJ))
    itpf.write('''      :::
    }
    ffio_constraints[0] {
      :::
      :::
    }
  }
}
''')
    itpf.close()
    return None


def mainBOSS2DESMOND(resid, clu):
    mol = pickle.load(open(resid + ".p", "rb"))
    if clu:
        pdb_file = 'clu.pdb'
    else:
        pdb_file = 'plt.pdb'
    boss2desmond(resid, mol, pdb_file)
    return None
