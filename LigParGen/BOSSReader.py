from __future__ import print_function
import os
import shutil
import numpy as np
from LigParGen.mol_boss import new_mol_info
import pandas as pd
from collections import OrderedDict

from LigParGen.fepzmat import BCC_file2zmat
from LigParGen.CreatZmat import GenMolRep 

def VerifyMolandSave(mol,charge,resname):
    if mol is not None: 
        import pickle
        assert (mol.MolData['TotalQ']['Reference-Solute'] ==
                charge), "PROPOSED CHARGE IS NOT POSSIBLE: SOLUTE MAY BE AN OPEN SHELL"
        pickle.dump(mol, open(resname + ".p", "wb"))
    else: 
        print('Problem Detected Molecule Object Not created')
    return None

def LinCheck(fname):
    imp_dat = 0
    zlines  = open(fname,'r').readlines()
    for l in range(len(zlines)):
        if 'Geometry Variations follow ' in zlines[l]: imp_dat = l
    Atypes = []
    for l in zlines[1:imp_dat]:Atypes.append(l.split()[2])
    Atypes = np.array(Atypes,dtype=int)
    Atypes = Atypes[Atypes<0]
    Check =False
    if len(Atypes)>2: Check = True
    return Check

def mod_add_diheds(line):
    adihed = [int(l) for l in line.split()[0:4]]+[-1,-1]
    return(adihed)

def fix_add_dihed(zmat_name):
    flines = open('%s.z'%zmat_name,'r').readlines()
    imp_lines=[]
    for l in range(len(flines)):
        if 'Additional Dihedrals follow' in flines[l]:
            imp_lines.append(l+1)
        elif 'Domain Definitions follow' in flines[l]:
            imp_lines.append(l)
    ofile = open('%s_fixed.z'%zmat_name,'w+')
    for line in flines[0:imp_lines[0]]: ofile.write('%s\n'%(line.rstrip()))
    for line in flines[imp_lines[0]:imp_lines[1]]:
        m_ad = mod_add_diheds(line)
        ofile.write('%4d%4d%4d%4d%4d%4d\n'%(m_ad[0],m_ad[1],m_ad[2],m_ad[3],m_ad[4],m_ad[5]))
    for line in flines[imp_lines[1]:]: ofile.write('%s\n'%line.rstrip())
    ofile.close()
    return(None)

def CheckForHs(atoms):
    atype = [line.split()[1][0] for line in atoms]
    ans = False
    if ('H' in atype): ans = True 
    return ans 

def bcc_db():
    '''
    19 LBCCs from 1.14*CM1A-LBCC paper
    '''
    lbcc = {'C#-C=': 0.0,
  'C-N': 0.0,
  'C-O': 0.05,
  'C-OE': 0.0,
  'C-OH': 0.0,
  'C-OS': 0.0,
  'CA-Br': 0.19,
  'CA-C': 0.0,
  'CA-C!': -0.0,
  'CA-C=': 0.0,
  'CA-CB': -0.0,
  'CA-CE': 0.0,
  'CA-CF': 0.0,
  'CA-CK': -0.0,
  'CA-CT': 0.0,
  'CA-CZ': 0.0,
  'CA-CZA': 0.0,
  'CA-Cl': 0.0,
  'CA-F': 0.13,
  'CA-I': 0.0,
  'CA-N3': 0.0,
  'CA-NC': 0.07,
  'CA-NO': -0.08,
  'CA-NP': 0.06,
  'CA-NS': 0.0,
  'CA-OH': 0.22,
  'CA-OS': -0.0,
  'CA-S': -0.0,
  'CA-SH': -0.0,
  'CAM-CA': 0.0,
  'CAM-CT': 0.0,
  'CAM-N': 0.0,
  'CAM-O': 0.0,
  'CB-C=': -0.0,
  'CB-NC': -0.0,
  'CE-O': -0.0,
  'CE-OE': 0.0,
  'CE-OS': 0.0,
  'CF-F': -0.0,
  'CF-OS': -0.0,
  'CK-O': -0.0,
  'CM-C': 0.0,
  'CM-C=': -0.0,
  'CM-CT': -0.0,
  'CM-Cl': -0.0,
  'CP-CS': 0.0,
  'CP-SA': -0.0,
  'CT-Br': 0.08,
  'CT-C': -0.0,
  'CT-C=': 0.0,
  'CT-CE': -0.0,
  'CT-CF': 0.0,
  'CT-CK': -0.0,
  'CT-CP': 0.0,
  'CT-CZ': -0.0,
  'CT-CZT': -0.0,
  'CT-Cl': 0.1,
  'CT-F': -0.0,
  'CT-I': -0.0,
  'CT-N': -0.0,
  'CT-N3': -0.0,
  'CT-NO': 0.0,
  'CT-NP': 0.04,
  'CT-NS': -0.0,
  'CT-NT': -0.0,
  'CT-OE': -0.0,
  'CT-OH': 0.1,
  'CT-OS': -0.0,
  'CT-S': 0.08,
  'CT-SH': 0.175,
  'CT-SZ': 0.0,
  'CY-C': 0.0,
  'CY-CE': 0.0,
  'CZ-NZ': -0.0,
  'CZA-NZ': 0.09,
  'CZT-NZ': 0.03,
  'H-N': -0.0,
  'H-N3': -0.0,
  'H-NP': -0.05,
  'H-NS': -0.0,
  'H-NT': -0.0,
  'HA-CA': -0.01,
  'HA-CM': 0.0,
  'HA-CP': -0.0,
  'HA-CS': -0.0,
  'HC-C': 0.0,
  'HC-C#': -0.0,
  'HC-C=': -0.0,
  'HC-CAM': 0.0,
  'HC-CE': 0.0,
  'HC-CF': -0.0,
  'HC-CM': -0.0,
  'HC-CT': 0.0,
  'HC-CY': 0.0,
  'HC-CZ': -0.0,
  'HO-OH': 0.0,
  'HS-SH': 0.0,
  'NO-ON': -0.18,
  'O-P': 0.0,
  'OS-P': 0.0,
  'OY-SZ': 0.06,
  'U-U': 0.0,
  'X-X': 0.0}
    db = OrderedDict(lbcc)
    return db


def Refine_PDB_file(fname):
    flines = open(fname, 'r+').readlines()
    pdb_lines = []
    for line in flines:
        if ('ATOM' in line) or ('HETATM' in line):
            line = line.rstrip()
            line = line.lstrip()
            if not 'DUM' in line:
                pdb_lines.append(line)
    return pdb_lines


def get_coos_from_pdb(pdb_dat):
    atoms = []
    coos = []
    for line in pdb_dat:
        atom = line.split()[2]
        x, y, z = line[28:56].split()
        atoms.append(atom)
        coos.append([float(x), float(y), float(z)])
    return (atoms, coos)

def pairing_func(a, b):
    ans = (a + b) * (a + b + 1) * 0.5
    if a > b:
        ans = ans + a
        pans = '%6d%6d' % (b, a)
    else:
        ans = ans + b
        pans = '%6d%6d' % (a, b)
    return (int(ans), pans)


def ucomb(vec, blist):
    res = 0
    for a in vec:
        vec.remove(a)
        for b in vec:
            ans = (a + b) * (a + b + 1) * 0.5
            if (ans + a in blist) or (ans + b in blist):
                res = res + 1
    return res


def tor_cent(vec, blist):
    db = {}
    for a in vec:
        na = 0
        for b in vec:
            ans = (a + b) * (a + b + 1) * 0.5
            if (ans + a in blist) or (ans + b in blist):
                na += 1
        db[a] = na
    new_vec = list(sorted(db, key=db.__getitem__, reverse=True))
    return (new_vec)


def bossPdbAtom2Element(attype):
    elem = ''.join([i for i in attype[:-1] if not i.isdigit()])
    return elem


def bossElement2Mass(elem):
    symb2mass = {
        'H': 1.008,
        'F': 18.998403163,
        'Cl': 35.45,
        'Br': 79.904,
        'I': 126.90447,
        'O': 15.999,
        'S': 32.06,
        'N': 14.007,
        'P': 30.973761998,
        'C': 12.011,
        'Si': 28.085,
        'Na': 22.98976928,
        'SOD': 22.98976928,
        'K': 39.0983,
        'Mg': 24.305,
        'Ca': 40.078,
        'Mn': 54.938044,
        'Fe': 55.845,
        'Co': 58.933194,
        'Ni': 58.6934,
        'Cu': 63.546,
        'Zn': 65.38, }
    try:
        res = symb2mass[elem]
    except NameError:
        print("Mass for atom %s is not available \n add it to symb2mass dictionary")
    return res


def Refine_file(fname):
    flines = open(fname, 'r+')
    lines = []
    for line in flines:
        if line.rstrip():
            line = line.rstrip()
            line = line.lstrip()
            lines.append(line)
    flines.close()
    return lines


# banner text -> impDat key(s) it is expected to set. Module-level (not
# inline in get_ImpDat) so tests can locate section boundaries in captured
# BOSS out/sum text without needing a real BOSS install -- see
# tests/conftest.py, which calls find_boss_sections() directly rather than
# keeping its own copy of this banner list (a prior hand-copied duplicate
# risked silently desyncing from this one).
ODAT_BANNERS = [
    ('Z-Matrix for Reference Solutes', ['ATMinit']),
    ('Net Charge', ['TotalQ']),
    ('OPLS Force Field Parameters', ['ATMfinal', 'NBDinit']),
    ('Fourier Coefficients', ['TORinit', 'NBDfinal']),
    ('Bond Stretching Parameters', ['TORfinal', 'BNDinit']),
    ('Angle Bending Parameters', ['BNDfinal', 'ANGinit']),
    ('Non-bonded Pairs List', ['ANGfinal', 'PAIRinit']),
    ('Solute 0:   X          Y          Z', ['XYZinit']),
    ('Atom I      Atom J      RIJ', ['XYZfinal']),
    ('Checking', ['PAIRfinal']),
]
SDAT_BANNERS = [
    ('Additional Dihedrals follow', ['ADDinit']),
    ('Domain Definitions follow', ['ADDfinal']),
]


def find_boss_sections(odat, sdat, zmat_name='<zmat>'):
    """Locate every BOSS out/sum section boundary get_ImpDat() needs, by
    banner-text match, and fail loudly (naming the missing banner, or an
    empty/truncated section) instead of a bare KeyError or bad slice far
    away from the actual cause."""
    impDat = {}
    for nl in range(len(odat)):
        if 'Z-Matrix for Reference Solutes' in odat[nl]:
            impDat['ATMinit'] = nl
        elif 'Net Charge' in odat[nl]:
            impDat['TotalQ'] = nl
        elif 'OPLS Force Field Parameters' in odat[nl]:
            impDat['ATMfinal'] = nl
            impDat['NBDinit'] = nl
        elif 'Fourier Coefficients' in odat[nl]:
            impDat['TORinit'] = nl
            impDat['NBDfinal'] = nl
        elif 'Bond Stretching Parameters' in odat[nl]:
            impDat['TORfinal'] = nl
            impDat['BNDinit'] = nl
        elif 'Angle Bending Parameters' in odat[nl]:
            impDat['BNDfinal'] = nl
            impDat['ANGinit'] = nl
        elif 'Non-bonded Pairs List' in odat[nl]:
            impDat['ANGfinal'] = nl
            impDat['PAIRinit'] = nl
        elif 'Solute 0:   X          Y          Z' in odat[nl]:
            impDat['XYZinit'] = nl
        elif 'Atom I      Atom J      RIJ' in odat[nl]:
            impDat['XYZfinal'] = nl
        elif 'Checking' in odat[nl]:
            impDat['PAIRfinal'] = nl
#### THIS PART IS READ FROM SUM FILE ###
    for ml in range(len(sdat)):
        # LigParGen's own auto-generated Zmats always print this banner as
        # "Additional Dihedrals follow (6I4)". BOSS's own reference Zmat
        # library (molecules/small/*.z etc.) instead prints "Additional
        # Dihedrals (6I4) - Zero types not shown" when run through xSPM
        # alone (confirmed directly against real BOSS output for
        # molecules/small/acetam.z) -- same section, same column format,
        # different banner text depending on which BOSS code path wrote
        # it. Match either.
        if 'Additional Dihedrals follow' in sdat[ml] or 'Additional Dihedrals (' in sdat[ml]:
            impDat['ADDinit'] = ml
        elif 'Domain Definitions follow' in sdat[ml]:
            impDat['ADDfinal'] = ml
#### THIS PART IS READ FROM SUM FILE ###

    # A molecule too small to have some interaction type (water has no
    # torsions; a bare monatomic ion has no bonds/angles/torsions/pairs at
    # all) makes BOSS omit that section's banner entirely from /tmp/out --
    # not print it empty, just never print it -- confirmed directly by
    # inspecting real BOSS output for 'O' (water) and '[Cl-]'. Each of
    # these four banners marks BOTH the end of the section before it and
    # the start of the section after, so when one is missing, both of its
    # keys legitimately collapse to the nearest following boundary that
    # *was* found (walking backwards from 'Checking', the fixed anchor
    # that always terminates this chain). The resulting slice then simply
    # picks up a little extra neighboring text (e.g. the 'Net Charge'
    # block), which is harmless: every get_* consumer of these slices
    # (get_QLJ/get_bonds/get_angs/get_tors) only keeps lines matching a
    # specific data-row shape and ignores everything else. This only
    # fires once 'OPLS Force Field Parameters' and 'Checking' -- always
    # present for a run that actually completed -- are both found; if
    # either is missing too, that's a real failure and falls through to
    # the strict check below unchanged.
    optional_chain = [
        ('Fourier Coefficients', ('NBDfinal', 'TORinit')),
        ('Bond Stretching Parameters', ('TORfinal', 'BNDinit')),
        ('Angle Bending Parameters', ('BNDfinal', 'ANGinit')),
        ('Non-bonded Pairs List', ('ANGfinal', 'PAIRinit')),
    ]
    collapsed_banners = set()
    if 'NBDinit' in impDat and 'PAIRfinal' in impDat:
        fallback = impDat['PAIRfinal']
        for banner, keys in reversed(optional_chain):
            if keys[0] in impDat:
                fallback = impDat[keys[0]]
            else:
                impDat[keys[0]] = fallback
                impDat[keys[1]] = fallback
                collapsed_banners.add(banner)

    missing_banners = [
        banner for banner, keys in ODAT_BANNERS + SDAT_BANNERS
        if any(key not in impDat for key in keys)
    ]
    if missing_banners:
        raise ValueError(
            "BOSSReader (%s): could not locate the following expected "
            "section banner(s) in the BOSS output (/tmp/out or "
            "/tmp/sum): %s. The BOSS run may have failed, or its output "
            "format has changed." % (
                zmat_name, ', '.join(repr(b) for b in missing_banners)))

    # Sections that must contain at least one line of content once the
    # banner-derived start/end indices are used to slice odat below.
    # ('Additional Dihedrals' is intentionally excluded: an empty block
    # there just means the solute has no additional dihedrals, which is
    # normal -- and BONDS/ANGLES/TORSIONS/PAIRS are excluded too exactly
    # when their own governing banner above was legitimately absent and
    # collapsed to a zero-width slice: that's the same "normal, not an
    # error" case, just one banner earlier in the chain.)
    non_empty_sections = [
        ('ATOMS', 'ATMinit', 'ATMfinal', None),
        ('Q_LJ (non-bonded)', 'NBDinit', 'NBDfinal', None),
        ('BONDS', 'BNDinit', 'BNDfinal', 'Bond Stretching Parameters'),
        ('ANGLES', 'ANGinit', 'ANGfinal', 'Angle Bending Parameters'),
        ('TORSIONS', 'TORinit', 'TORfinal', 'Fourier Coefficients'),
        ('XYZ', 'XYZinit', 'XYZfinal', None),
        ('PAIRS', 'PAIRinit', 'PAIRfinal', 'Non-bonded Pairs List'),
    ]
    for name, start_key, end_key, governing_banner in non_empty_sections:
        if governing_banner is not None and governing_banner in collapsed_banners:
            continue
        start, end = impDat[start_key], impDat[end_key]
        if end <= start:
            raise ValueError(
                "BOSSReader (%s): the '%s' section of the BOSS output "
                "is empty or out of order (lines %d:%d) -- the output "
                "may be malformed." % (zmat_name, name, start, end))

    if impDat['TotalQ'] + 4 > len(odat):
        raise ValueError(
            "BOSSReader (%s): the 'Net Charge' section is truncated -- "
            "expected 4 lines starting at line %d but /tmp/out only has "
            "%d lines." % (zmat_name, impDat['TotalQ'], len(odat)))

    return impDat


class BOSSReader(object):

    def __init__(self, zmatrix, optim, charge=0, lbcc=False):
        self.zmat = zmatrix
        self.impDat = {}
        self.MolData = {}
        self.refine_data(optim, charge, lbcc)

    def Get_OPT(self, optim, charge):
        assert os.path.isfile(
            self.zmat), 'File named %10s does not exist' % self.zmat
        assert ('BOSSdir' in os.environ) and os.path.isfile((os.environ[
            'BOSSdir'] + '/scripts/xZCM1A')), 'Please Make sure $BOSSdir is defined \n xZCM1A and related files are in scripts directory of BOSS'
        execs = {
             2: os.environ['BOSSdir'] + '/scripts/xZCM1A+2 > /tmp/olog',
             1: os.environ['BOSSdir'] + '/scripts/xZCM1A+  > /tmp/olog',
             0: os.environ['BOSSdir'] + '/scripts/xZCM1A > /tmp/olog',
            -1: os.environ['BOSSdir'] + '/scripts/xZCM1A-  > /tmp/olog',
            -2: os.environ['BOSSdir'] + '/scripts/xZCM1A-2 > /tmp/olog',
         'OPT': os.environ['BOSSdir'] + '/scripts/xOPT >/tmp/olog'
        }
        #print('MOLECULE HAS A CHARGE of %d' % charge)
        if optim > 0:
            print('Optimization level requested %d' % optim)
            for opt_lev in range(optim):
                print('Performing Stage %d of Charge Generation'%(opt_lev+1)) 
                execfile = execs[charge]
                coma = execfile + ' ' + self.zmat[:-2]
                os.system(coma)
                shutil.copyfile('sum', self.zmat)
                execfile = execs['OPT']
                coma = execfile + ' ' + self.zmat[:-2]
                os.system(coma)
                shutil.copyfile('sum', self.zmat)
                os.system('head -1 %s'% (self.zmat))
                #os.system('cd /tmp;/bin/cp sum %s' % (self.zmat))
        execfile = os.environ['BOSSdir'] + '/scripts/xSPM > /tmp/olog'
        coma = execfile + ' ' + self.zmat[:-2]
        os.system(coma)
        shutil.copyfile('/tmp/sum', '/tmp/' + self.zmat)
        return (None)

    def get_addihed(self, data):
        add = []
        nadd = 0
        for line in data:
            if line[0].isdigit():
                add.append(line.split()[0:4])
                nadd = nadd + 1
        return (add)

    def get_atinfo(self, data):
        ats = []
        nat = 0
        for line in data:
            if line[0].isdigit() and float(line.split()[2]) > 1:
                ats.append(line)
                nat += 1
        return (ats)

    def get_charge(self, data):
        TotQ = {}
        # get_ImpDat() hands this a fixed 4-line window (banner + Reference
        # Solute + 1st/2nd Perturbed Solute). Guard against a truncated
        # window (e.g. the banner sat too close to EOF) instead of silently
        # parsing fewer charge entries than expected.
        expected_lines = 4
        if len(data) < expected_lines:
            raise ValueError(
                "get_charge(): expected a %d-line 'Net Charge' block (banner "
                "+ 3 solute charge lines) but only got %d line(s) -- BOSS "
                "output may be truncated." % (expected_lines, len(data)))
        for line in data[1:expected_lines]:
            words = line.split()
            try:
                charge_val = round(float(words[-1]), 3)
            except (IndexError, ValueError):
                raise ValueError(
                    "get_charge(): could not parse a net charge value from "
                    "line %r" % line)
            TotQ['-'.join(words[:-1])] = charge_val
        if 'Reference-Solute' not in TotQ:
            raise ValueError(
                "get_charge(): 'Reference Solute' net charge line not found "
                "in the parsed block (got: %s) -- downstream charge checks "
                "would silently fail." % (list(TotQ.keys()),))
        return TotQ

    def get_tors(self, data):
        tors = []
        ntor = 0
        for line in data:
            if 'All Solutes' in line:
                tors.append(line.split()[4:8])
                for tor in line.split()[4:8]:
                    if abs(float(tor)) > 0.0:
                        ntor = ntor + 1
        return (tors)

    def get_tors_by_decl_idx(self, data):
        """Same 'Fourier Coefficients' section as get_tors(), but keyed by
        each row's own 'Angle' column (data[0]) -- BOSS's 1-based index
        into the full declared-dihedral sequence (Variable Dihedrals
        follow, then Additional Dihedrals follow, in that order) -- rather
        than returned as a plain list.

        This table is not guaranteed to have exactly one row per declared
        dihedral quadruple: BOSS can skip a quadruple entirely if it can't
        match its atom-type pattern against a known torsion type (observed
        directly: declaring a ring's internal bonds/angles as Additional
        entries makes BOSS re-derive more Additional Dihedrals than it
        ends up tabulating rows for). Assuming row N of this table always
        corresponds to the Nth declared quadruple -- what the plain list
        from get_tors() invites -- silently mispairs coefficients to atoms
        as soon as any row is skipped. The 'Angle' column sidesteps that:
        it's BOSS's own declaration-order index, so row N's coefficients
        can be matched to declared quadruple N directly, keyed lookup
        instead of positional zip, correct even with skipped rows.
        """
        tors = {}
        for line in data:
            if 'All Solutes' in line:
                fields = line.split()
                decl_idx = int(fields[0])
                tors[decl_idx] = fields[4:8]
        return tors

    def get_QLJ(self, data):
        qlj = []
        nqlj = 0
        for line in data:
            if 'All Solutes' in line and line[0].isalpha():
                qlj.append([line.split()[0], line.split()[2],
                            line.split()[3], line.split()[4]])
                nqlj += 1
        return (qlj)

    def get_angs(self, data):
        angs = {'cl1': [], 'cl2': [], 'cl3': [], 'R': [], 'K': []}
        nang = 0
        for line in data:
            word = line.split()
            # A real angle data row looks like:
            #   Atom1 Atom2 Atom3  A0  K0  A1  K1  A2  K2  Delta  AtomTypes
            # Header/banner lines don't have this shape (too few fields, or
            # the leading fields aren't atom-index integers). Detecting rows
            # this way -- instead of thresholding on K0's value -- means a
            # genuine (if unusual) zero force constant is still kept rather
            # than being silently treated as "not a data row".
            if len(word) < 5:
                continue
            try:
                cl1 = int(word[0])
                cl2 = int(word[1])
                cl3 = int(word[2])
                r = float(word[3])
                k = float(word[4])
            except ValueError:
                continue
            angs['cl1'].append(cl1)
            angs['cl2'].append(cl2)
            angs['cl3'].append(cl3)
            angs['R'].append(r)
            angs['K'].append(k)
            nang = nang + 1
            #        print 'Total No of Non-zero Angles in BOSS is %d' % (nang)
        return (angs)

    def get_XYZ(self, data):
        XYZ = {'at_num': [], 'X': [], 'Y': [], 'Z': [], 'at_symb': []}
        for line in data:
            if line[0].isdigit() and len(line.split()) == 5:
                word = line.split()
                if int(word[0]) > 0:
                    XYZ['at_num'].append(int(word[0]))
                    XYZ['X'].append(float(word[1]))
                    XYZ['Y'].append(float(word[2]))
                    XYZ['Z'].append(float(word[3]))
                    XYZ['at_symb'].append(word[4])
        XYZ = pd.DataFrame(XYZ)
        return XYZ

    def get_pairs(self, data):
        data = data[1:]
        plnos = []
        for i in range(0, len(data)):
            if 'Atom' in data[i]:
                plnos.append(i)
        if not plnos:
            raise ValueError(
                "get_pairs(): could not find any 'Atom N:' markers in the "
                "Non-bonded Pairs List section -- BOSS output may be "
                "malformed or its banner text may have changed.")
        natoms = len(plnos)
        plnos.append(len(data))
        pair_dat = {i: ' '.join(data[plnos[i]:plnos[i + 1]])
                    for i in range(len(plnos) - 1)}
        for nu in range(natoms):
            marker_line = data[plnos[nu]]
            if ':' not in marker_line:
                raise ValueError(
                    "get_pairs(): expected an 'Atom N:' marker line but got "
                    "%r" % marker_line)
            colon_idx = marker_line.index(':')
            # Strip exactly the 'Atom N:' prefix found on this line (whatever
            # its width, e.g. 'Atom    5:' vs 'Atom   40:') instead of
            # assuming a fixed 10-column offset. Cross-check the atom number
            # in it against the sequential 1..natoms position this parser
            # assumes -- if the markers aren't sequential, `nu` below would
            # get silently paired with the wrong atom's exclusions.
            try:
                atom_label = int(
                    marker_line[:colon_idx].replace('Atom', '').strip())
            except ValueError:
                raise ValueError(
                    "get_pairs(): could not parse an atom number out of "
                    "marker line %r" % marker_line)
            if atom_label != nu + 1:
                raise ValueError(
                    "get_pairs(): 'Atom N:' markers in the Non-bonded Pairs "
                    "List are not sequential (expected Atom %d, found Atom "
                    "%d)" % (nu + 1, atom_label))
            pair_dat[nu] = list(pair_dat[nu][colon_idx + 1:].split())
            pair_dat[nu] = np.array([int(a) - 2 for a in pair_dat[nu]])
        pairs = []
        for k in pair_dat.keys():
            for j in pair_dat[k]:
                # Every parsed partner index should resolve to one of the
                # natoms atoms this section describes -- a value outside
                # that range means the (fixed-offset) column parsing above
                # is misaligned, and should fail loudly here rather than
                # produce a bogus pair entry silently.
                if not (0 <= j < natoms):
                    raise ValueError(
                        "get_pairs(): parsed pair partner index %d for atom "
                        "%d is out of range for a %d-atom Non-bonded Pairs "
                        "List -- parsing may be misaligned." %
                        (j, k, natoms))
                pairs.append('%6d%6d%6d\n' % (k - 1, j, 1))
        return pairs

    def get_bonds(self, data):
        bnds = {'cl1': [], 'cl2': [], 'RIJ': [], 'KIJ': [], 'TIJ': []}
        nbnd = 0
        for line in data:
            word = line.split()
            # A real bond data row looks like:
            #   Atom1 Atom2  R0  K0  R1  K1  R2  K2  Delta  AtomTypes
            # (note: for single-character atom types, e.g. 'H -NT', the
            # trailing AtomTypes field itself splits into two extra tokens,
            # so column *count* isn't checked exactly -- only that there are
            # enough fields, and that the leading atom indices/R0/K0 fields
            # actually parse.) Header/banner lines fail this shape check.
            # Using shape+parseability instead of thresholding on K0's value
            # means a genuine (if unusual) zero force constant is kept
            # rather than being silently treated as "not a data row".
            if len(word) < 4:
                continue
            try:
                cl1 = int(word[0])
                cl2 = int(word[1])
                rij = float(word[2])
                kij = float(word[3])
            except ValueError:
                continue
            bnds['cl1'].append(cl1)
            bnds['cl2'].append(cl2)
            bnds['RIJ'].append(rij)
            bnds['KIJ'].append(kij)
            bnds['TIJ'].append(line[-5:])
            nbnd += 1
        return (bnds)

    def prep_lbcc(self, bond_data, qdata):
        db = bcc_db()
        bnd_df = pd.DataFrame(bond_data)
        bnd_df = bnd_df[['cl1', 'cl2']]
        bnd_df.columns = ['I', 'J']
        q_df = pd.DataFrame(columns=['TY', 'Q'])
        q_df.loc[0] = ['1', 0.000]
        q_df.loc[1] = ['2', 0.000]
        for i in range(len(qdata)):
            q_df.loc[i + 2] = [qdata[i][0], float(qdata[i][1])]
        bond, cha, QBC1 = new_mol_info(db, q_df, bnd_df)
        lbcc_qdat = []
        for i in range(len(qdata)):
            lbcc_qdat.append(
                [qdata[i][0], str(cha.QBCC.values[i]), qdata[i][2], qdata[i][3]])
        bond.to_csv('LBCC_BONDS.csv', index=False)
        cha.to_csv('LBCC_CHARGES.csv', index=False)
        return np.array(cha.QBCC), lbcc_qdat

    def cleanup(self):
        for fname in ('sum', 'log', 'olog', 'out', 'plt.pdb'):
            fpath = os.path.join('/tmp', fname)
            if os.path.exists(fpath):
                os.remove(fpath)

    def get_ImpDat(self, optim, charge):
        self.Get_OPT(optim, charge)
        odat = Refine_file('/tmp/out')
        sdat = Refine_file('/tmp/sum')
        MolData = {}
        MolData['PDB'] = Refine_file('/tmp/plt.pdb')

        impDat = find_boss_sections(odat, sdat, zmat_name=self.zmat)

        MolData['ATOMS'] = self.get_atinfo(
            odat[impDat['ATMinit']:impDat['ATMfinal']])
        MolData['Q_LJ'] = self.get_QLJ(
            odat[impDat['NBDinit']:impDat['NBDfinal']])
        MolData['BONDS'] = self.get_bonds(
            odat[impDat['BNDinit']:impDat['BNDfinal']])
        MolData['ANGLES'] = self.get_angs(
            odat[impDat['ANGinit']:impDat['ANGfinal']])
        MolData['TORSIONS'] = self.get_tors(
            odat[impDat['TORinit']:impDat['TORfinal']])
        MolData['TORSIONS_BY_DECL_IDX'] = self.get_tors_by_decl_idx(
            odat[impDat['TORinit']:impDat['TORfinal']])
        MolData['ADD_DIHED'] = self.get_addihed(
            sdat[impDat['ADDinit']:impDat['ADDfinal']])
        MolData['XYZ'] = self.get_XYZ(
            odat[impDat['XYZinit']:impDat['XYZfinal']])
        # get_pairs() requires at least one 'Atom N:' marker in its slice --
        # correct for a real (even if all-empty) Non-bonded Pairs List, but
        # a molecule too small to have one at all (find_boss_sections()
        # collapses it to a zero-width slice in that case) has no markers
        # to find, so skip straight to an empty result instead of tripping
        # get_pairs()'s own "output may be malformed" check on a slice that
        # was never expected to contain anything.
        if impDat['PAIRinit'] == impDat['PAIRfinal']:
            MolData['PAIRS'] = []
        else:
            MolData['PAIRS'] = self.get_pairs(
                odat[impDat['PAIRinit']:impDat['PAIRfinal']])
        MolData['TotalQ'] = self.get_charge(
            odat[impDat['TotalQ']:impDat['TotalQ'] + 4])
        return MolData

    def refine_data(self, optim, charge, lbcc):
        if lbcc and (charge == 0):
            lbcc_MD = self.get_ImpDat(optim, charge)
            QLBCC, DATA_Q_LJ = self.prep_lbcc(
                lbcc_MD['BONDS'], lbcc_MD['Q_LJ'])
            lbcc_MD['Q_LJ'] = DATA_Q_LJ
            BCC_file2zmat(self.zmat, QLBCC,
                          oname='/tmp/%s_BCC.z' % self.zmat[:-2])
            os.replace('%s.z' % self.zmat[:-2], '%s_NO_LBCC.z' % self.zmat[:-2])
            os.replace('%s_BCC.z' % self.zmat[:-2], '%s.z' % self.zmat[:-2])
            self.MolData = lbcc_MD
        elif lbcc and (charge != 0):
            print('LBCC IS SUPPORTED ONLY FOR NEUTRAL MOLECULES')
        else:
            self.MolData = self.get_ImpDat(optim, charge)
        return None
