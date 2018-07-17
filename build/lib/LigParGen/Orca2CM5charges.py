'''
Orca2CM5charges.py 
Script to Get CM5 charges from Hirshfeld Charges 
Input: Log File from Orca Calculation
Output: CSV file with coordinates, raw CM5 charges, rdkit atomic FP averaged CM5 charges
Written by: Leela S. Dodda (Jorgensen Group@Yale)
Email Id: leela.dodda@yale.edu
'''
import pandas as pd
import numpy as np
import os
import argparse
from rdkit import Chem
from rdkit.Chem import AllChem

def AddCM5Charges(mol,qcm5): 
    qcm5.to_csv('CM5_charges.csv',index=False,float_format='%6.4f')
    cm1a = mol.MolData['Q_LJ'] 
    cm5_charges = list(qcm5.CM5_final)
    assert len(cm1a) == len(qcm5), 'Please check the at_info and Q_LJ_dat files'
    for c in range(len(cm1a)): 
        cm1a[c][1] = '%9.6f'%cm5_charges[c]
    mol.MolData['Q_LJ'] = cm1a
    return mol 

def AtomFPProgram(hmol,atomNum,radii=2):
    env = Chem.FindAtomEnvironmentOfRadiusN(hmol,radii,atomNum,useHs=True)
    amap = {}
    submol = Chem.PathToSubmol(hmol,env,atomMap=amap)
    atom_smi = Chem.MolToSmiles(submol)
    return('%d'%(hmol.GetAtomWithIdx(atomNum).GetAtomicNum())+atom_smi)

def GetAvals(df):
    num2a = (df.set_index(['A0_NO'])['VALUE'].to_dict())
    list_keys = list(num2a.keys())
    dvals = (np.empty([np.max(list_keys),np.max(list_keys)]))
    for i in range(dvals.shape[0]):
        for j in range(dvals.shape[1]):
            if (i != j) : dvals[i,j] = num2a[i+1]-num2a[j+1]
    dvals[ 0, 5]= 0.0502
    dvals[ 0, 6]= 0.1747
    dvals[ 0, 7]= 0.1671
    dvals[ 5, 6]= 0.0556
    dvals[ 5, 7]= 0.0234
    dvals[ 6, 7]=-0.0346
    # Repitition of the above cooefficients with a negative sign
    dvals[ 5, 0]=-0.0502
    dvals[ 6, 0]=-0.1747
    dvals[ 7, 0]=-0.1671
    dvals[ 6, 5]=-0.0556
    dvals[ 7, 5]=-0.0234
    dvals[ 7, 6]= 0.0346
    return dvals


def GetLogFile(fname,pt_df,rad_df):
    pt_df["symbol"] = pt_df["symbol"].map(str.strip)
    sym2num = pt_df.set_index(['symbol'])['atomicNumber'].to_dict() 
    num2rad = rad_df.set_index(['RAD_NO'])['VALUE'].to_dict()
    xyz_data = []
    charge_data = []
    data = open(fname).readlines()
    id_charges = False
    id_coos = False
    for line in data: 
        if 'CARTESIAN COORDINATES (ANGSTROEM)' in line: 
            id_coos = True 
        elif 'CARTESIAN COORDINATES (A.U.)' in line:
            id_coos = False 
        if 'HIRSHFELD ANALYSIS' in line: 
            id_charges=True
        elif 'TIMINGS' in line: 
            id_charges = False
        if id_charges:charge_data.append(line.strip().split()) 
        if id_coos:xyz_data.append(line.strip().split()) 
    hirCharges = pd.DataFrame(charge_data[7:-4],columns=['N','ATOM','QHir','Spin'])
    hirCharges[['N','QHir','Spin']] = hirCharges[['N','QHir','Spin']].apply(pd.to_numeric)
    hirCharges = hirCharges[['N','QHir']]
    xyzcoos = pd.DataFrame(xyz_data[2:-2],columns=['ATOM','X','Y','Z'])
    xyzcoos[['X','Y','Z']] = xyzcoos[['X','Y','Z']].apply(pd.to_numeric)
    final_data = (pd.concat([xyzcoos,hirCharges],axis=1))
    final_data['AtNum'] = [sym2num[s] for s in final_data.ATOM] 
    final_data['RAD'] = [num2rad[s] for s in final_data.AtNum] 
    return(final_data)

def Distance(a,b):
    return(np.sqrt((a[0]-b[0])**2+(a[1]-b[1])**2+(a[2]-b[2])**2))

def HirshfeldToCM5(df,a0,netcharge): 
    DVALS=GetAvals(a0)
    cm5_charges = []
    alpha = 2.474
    for i,r in df.iterrows(): 
        qcm5 = r.QHir
        for j,p in df.iterrows(): 
            if (r.AtNum != p.AtNum): 
                dist = Distance([r.X,r.Y,r.Z],[p.X,p.Y,p.Z])
                factor = np.exp(-1.0*alpha*(dist-r.RAD-p.RAD))
                qcm5=qcm5+factor*DVALS[r.AtNum-1,p.AtNum-1]
        cm5_charges.append(qcm5)
    df['QCM5']     = np.array(cm5_charges) 
    #df['QCM5']     = df.QCM5 - (netcharge-df.QCM5.sum())/len(df.QCM5) 
    mol = xyz_prep(df)
    df['FPS'] = [AtomFPProgram(mol,atomNum,radii=2) for atomNum in df.index]
    uniq_fps = list(set(df.FPS))
    df['QCM5_AVG'] = [df[df.FPS==i].QCM5.mean() for i in df.FPS]
    df['DEV_CM5'] = df.QCM5_AVG - df.QCM5
    if netcharge == 0:
        df['CM5_final'] = df.QCM5_AVG*1.20
    else: 
        df['CM5_final'] = df.QCM5_AVG*1.00
    #print(df.sum())
    return(df)

def xyz_prep(df):
    opdb = open('inp_orca.xyz', 'w+')
    opdb.write('%3d\n'%(len(df.QCM5)))
    opdb.write('REMARK LIGPARGEN GENERATED XYZ FILE\n')
    num = 0
    for (i, r) in df.iterrows(): 
        opdb.write('%-6s    %8.3f%8.3f%8.3f\n' %
                   (r.ATOM, r.X, r.Y, r.Z))
    opdb.close()
    os.system('babel -ixyz inp_orca.xyz -omol UNK.mol')
    os.system('babel -ixyz inp_orca.xyz -opdb inp_orca.pdb')
    hmol = Chem.MolFromMolFile('UNK.mol',removeHs=False,sanitize=False)
    return hmol 

def LoadModel(): 
    import json
    with open("cm5pars.json") as tweetfile:
        cm5_model = json.loads(tweetfile.read())
    a0_df = pd.DataFrame.from_dict(cm5_model['A0'])
    rd_df = pd.DataFrame.from_dict(cm5_model['radii'])
    pt_df = pd.DataFrame.from_dict(cm5_model['PeriodicTable'])
    return (a0_df,rd_df,pt_df)

