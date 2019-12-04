#!/usr/bin/env python
################################################################################
################################################################################
##  I have made it for Gromacs Amber03.ff forcefied to generate CG2AA index  ###
##  CG method implemented are martini, cafemol, unres, 03/23/2016, jhpeng    ###
##  Mainchain conformations are always remained, side-chains not             ###
################################################################################
################################################################################

version="1.0"
description="Only residue based CG models are covered"
authors=["Junhui Peng"]

# Parameters are defined for the following (protein) forcefields:
forcefields = ['martini','cafca','cafcom','unres','cafsidecom']

usage = """
python cgmap.py -f <input, AA PDB>
                -x <output, optional, CG PDB>
                -ff martini|cafca|cafcom|unres|cafsidecom
                -nmap <output, AA2CG map> 
                -s <input, target CG PDB, To get initial RMSD!!>
                -o <output, the .tmi for TMD>
                -p <input, fraction of not deviated CG sites, default is 1.0>
                -o <output, the .tmi for TMD>
                -v <bool, get version>
                -h <boll, get this help>
#######################################################################
"""
import sys
    
def getoption(args):
    # options, it can work... #
    options = {
                  '-f':None,    #!input!
                  '-x':None,    #!input!
                  '-ff':None,   #!input!
                  '-nmap':None, #!output!
                  '-s':None,    #optional
                  '-o':None,    #optional
                  '-p':1.0,     #optional
                  '-v':0,       #optional
                  '-h':0        #optional
              }
    if '-h' in args:
        sys.stdout.write(usage)
        sys.exit()
    if '-v' in args:
        sys.stdout.write("cgmap_infor>>> Version:%s\n"%version)
        sys.stdout.write("cgmap_infor>>> %s\n"%description)
    if '-f' not in args and '-nmap' not in args:
        sys.stdout.write(usage)
        sys.exit()
    for i in range(len(args)):
        if args[i].startswith("-"):
            if args[i] not in options.keys():
                sys.stdout.write(usage)
                sys.exit("cgmap error>> Unkown option %s"%args[i])
            else:
                try:
                    if args[i+1].startswith("-"):
                        sys.exit("cgmap_error>> Something wrong with the option %s"%args[i])
                except:
                    IndexError
        try:
            if args[i] == '-f':
                options['-f'] = args[i+1]
            if args[i] == '-x':
                options['-x'] = args[i+1]
            if args[i] == '-nmap':
                options['-nmap'] = args[i+1]
            if args[i] == '-ff':
                options['-ff'] = args[i+1]
            if args[i] == '-s':
                options['-s'] = args[i+1]
            if args[i] == '-p':
                options['-p'] = args[i+1]
            if args[i] == '-o':
                options['-o'] = args[i+1]
        except:
            IndexError
            sys.exit("cgmap_error>> Something wring with the options")
    return options

def get_cgmapping(method):
    if method == 'martini':
        CGmapping = CGmartini
    elif method == 'cafca':
        CGmapping = CGcafca
    elif method == 'cafcom':
        CGmapping = CGcafcom
    elif method == 'unres':
        CGmapping = CGunres
    elif method == 'dppc':
        CGmapping = CGDPPC
    else:
        sys.exit("cgmap error>> CG method error. Please select a supported CG method!")
    return CGmapping

#################################################
## 3 # HELPER FUNCTIONS, CLASSES AND SHORTCUTS ##  -> @FUNC <-
#################################################

import math
import numpy as np

# Split a string                                                              
def spl(x):                                                                   
    return x.split()                                                          

# Split each argument in a list                                               
def nsplit(*x):                                                               
    return [i.split() for i in x]                                             

def joinall(*x):
    return [" ".join(x)]

# Make a dictionary from two lists                                            
def hash(x,y):                                                                
    return dict(zip(x,y))                                                     

def get_inirms(coorsini, coorstar, mass, ncg):
    if ncg != len(mass) or  ncg != len(coorsini) or  ncg != len(coorstar):
        sys.exit("cgmap error>> numbers of ini and tar CG beads are different!")
    com1 = sum(np.array([coorsini[s]*mass[s] for s in range(ncg)])) / sum(mass)
    com2 = sum(np.array([coorstar[s]*mass[s] for s in range(ncg)])) / sum(mass)
    xini = coorsini - com1
    xtar = coorstar - com2
    # correlation matrix
    a = np.dot(np.transpose(xini), xtar)
    u, d, vt = np.linalg.svd(a)
    rot = np.transpose(np.dot(np.transpose(vt), np.transpose(u)))
    # is_reflection
    if np.linalg.det(rot) < 0:
        vt[2] = -vt[2]
        rot = np.transpose(np.dot(np.transpose(vt), np.transpose(u)))
    #tran = com2 - np.dot(com1, rot)
    xini = np.dot(xini, rot)# + tran
    vec = [(xini[s]-xtar[s])**2 for s in range(ncg)]
    rms = [sum(vec[s])*mass[s] for s in xrange(ncg)]
 
    # only RMSD is desired #
    return np.sqrt(sum(rms)/sum(mass))

##########################
## 4 # FG -> CG MAPPING ##  -> @MAP <-
##########################
mass = {'H': 1.008,
        'C': 12.01,
        'N': 14.01,
        'O': 16.00,
        'S': 32.06,
        'P': 30.97,
        'M': 0.000,
        'ZN':65.0 }
elements = {
         '1H2\'' :'H',
         '1H5\'' :'H',
         '2H2\'' :'H',
         '2H5\'' :'H',
         '1HD1':'H',
         '1HD2':'H',
         '2HD1':'H',
         '2HD2':'H',
         '3HD1':'H',
         '3HD2':'H',
         '1HE2':'H',
         '2HE2':'H',
         '1HG1':'H',
         '1HG2':'H',
         '1HH1':'H',
         '1HH2':'H',
         '2HG1':'H',
         '2HG2':'H',
         '3HG1':'H',
         '3HG2':'H',
         '2HH1':'H',
         '2HH2':'H',
         '3HG2':'H',
         '0C21':'C',
         '1C21':'C',
         '2C21':'C',
         '3C21':'C',
         '4C21':'C',
         '5C21':'C',
         '6C21':'C',
         '7C21':'C',
         '8C21':'C',
         '0C31':'C',
         '1C31':'C',
         '2C31':'C',
         '3C31':'C',
         '4C31':'C',
         '5C31':'C',
         '6C31':'C',
         'ZN':'ZN'
           }
def get_elem(atname):
    if atname[0] in mass.keys():
        return atname[0]
    else:
        return elements[atname]

Protein = spl("TRP TYR PHE HIS HIH ARG LYS CYS ASP GLU ILE LEU MET ASN PRO HYP GLN SER THR VAL ALA GLY")
Lipid = spl("DPPC")
DNA = spl("DAD DCY DGU DTH ADE CYT GUA THY URA DA DC DG DT")
DNAT = spl("DA5 DC5 DG5 DT5") ### only for cafemol force field ###

class CGDPPC:
    phosphate_NC3 = nsplit("N C12 C13 C14 C15","C11 P1 O1 O2 O3 O4","C1 C2 O31 C3 HS","C21 O21 C22 O22")
    dppc1 = nsplit("C23 C24 C25 C26","C27 C28 C29 0C21","1C21 2C21 3C21 4C24","5C21 6C21 7C21 8C21")
    dppc2 = nsplit("C31 C32 C33 C34","C35 C36 C37 C38","C391 0C31 1C31 2C34","3C31 4C31 5C31 6C31")
    #phosphatidylserine       =
    mapping = {
        "DPPC": phosphate_NC3 + dppc1 + dppc2,
        }

    ### The definition may need to be changed!
    cgnames = ['NC3', 'PO4', 'GL1', 'GL2', 'C1A', 'C2A', 'C3A', 'C4A', 'C1B', 'C2B', 'C3B', 'C4B']

class CGmartini:
    # Protein backbone as Amber03.ff definition,  
    # O1 and O2 changed to OC1 and OC2
    bb        = "N CA C O H H1 H2 H3 OC1 OC2"                                                                  #@#  
    # Lipid tails
    palmitoyl1    = nsplit("C1B C1C C1D C1E","C1F C1G C1H C1I","C1J C1K C1L C1M","C1N C1O C1P")                #@#
    palmitoyl2    = nsplit("C2B C2C C2D C2E","C2F C2G C2H C2I","C2J C2K C2L C2M","C2N C2O C2P")                #@#
    oleyl1        = nsplit("C1B C1C C1D C1E","C1F C1G C1H","C1I C1J","C1K C1L C1M C1N","C1O C1P C1Q C1R")      #@#
    oleyl2        = nsplit("C2B C2C C2D C2E","C2F C2G C2H","C2I C2J","C2K C2L C2M C2N","C2O C2P C2Q C2R")      #@#
    #phoshpatidylcholine      = 
    phosphatydilethanolamine = nsplit("N H1 H2 H3 CA","CB P OA OB OC OD","CC CD OG C2A OH","CE OE C1A OF")     #@#
    phosphatidylglycerol     = nsplit("H1 O1 CA H2 O2 CB","CC P OA OB OC OD","CD CE OG C2A OH","CF OE C1A OF") #@#
    #phosphatidylserine       =

    dna_bb = "P OP1 OP2 O5' O3'","C5' O4' C4'","C3' O3' C2' C1'"

    # This is the mapping dictionary according to CG Martini
    mapping = {
        "ALA":  nsplit(bb + " CB"),
        "CYS":  nsplit(bb,"CB SG"),
        "ASP":  nsplit(bb,"CB CG OD1 OD2"),
        "GLU":  nsplit(bb,"CB CG CD OE1 OE2"),
        "PHE":  nsplit(bb,"CB CG CD1 HD1","CD2 HD2 CE2 HE2","CE1 HE1 CZ HZ"),
        "GLY":  nsplit(bb),
        "HIS":  nsplit(bb,"CB CG","CD2 HD2 NE2 HE2","ND1 HD1 CE1 HE1"),
        "HIH":  nsplit(bb,"CB CG","CD2 HD2 NE2 HE2","ND1 HD1 CE1 HE1"),     # Charged Histidine.
        "ILE":  nsplit(bb,"CB CG1 CG2 CD CD1"),
        "LYS":  nsplit(bb,"CB CG CD","CE NZ HZ1 HZ2 HZ3"),
        "LEU":  nsplit(bb,"CB CG CD1 CD2"),
        "MET":  nsplit(bb,"CB CG SD CE"),
        "ASN":  nsplit(bb,"CB CG ND1 ND2 OD1 OD2 HD11 HD12 HD21 HD22"),
        "PRO":  nsplit(bb,"CB CG CD"),
        "HYP":  nsplit(bb,"CB CG CD OD"),
        "GLN":  nsplit(bb,"CB CG CD OE1 OE2 NE1 NE2 HE11 HE12 HE21 HE22"),
        "ARG":  nsplit(bb,"CB CG CD","NE HE CZ NH1 NH2 HH11 HH12 HH21 HH22"),    
        "SER":  nsplit(bb,"CB OG HG"),
        "THR":  nsplit(bb,"CB OG1 HG1 CG2"),
        "VAL":  nsplit(bb,"CB CG1 CG2"),
        "TRP":  nsplit(bb,"CB CG CD2","CD1 HD1 NE1 HE1 CE2","CE3 HE3 CZ3 HZ3","CZ2 HZ2 CH2 HH2"),
        "TYR":  nsplit(bb,"CB CG CD1 HD1","CD2 HD2 CE2 HE2","CE1 HE1 CZ OH HH"),
        "POPE": phosphatydilethanolamine + palmitoyl1 + oleyl2,
        "DOPE": phosphatydilethanolamine + oleyl1     + oleyl2,
        "DPPE": phosphatydilethanolamine + palmitoyl1 + palmitoyl2,
        "POPG": phosphatidylglycerol     + palmitoyl1 + oleyl2,
        "DOPG": phosphatidylglycerol     + oleyl1     + oleyl2,
        "DPPG": phosphatidylglycerol     + palmitoyl1 + palmitoyl2,
        "DA": nsplit("P OP1 OP2 O5' O3' O1P O2P","C5' O4' C4'","C3' C2' C1'","N9 C4","C8 N7 C5","C6 N6 N1","C2 N3"),
        "DG": nsplit("P OP1 OP2 O5' O3' O1P O2P","C5' O4' C4'","C3' C2' C1'","N9 C4","C8 N7 C5","C6 O6 N1","C2 N2 N3"),
        "DC": nsplit("P OP1 OP2 O5' O3' O1P O2P","C5' O4' C4'","C3' C2' C1'","N1 C6","C5 C4 N4","N3 C2 O2"),
        "DT": nsplit("P OP1 OP2 O5' O3' O1P O2P","C5' O4' C4'","C3' C2' C1'","N1 C6","C5 C4 O4 C7 C5M","N3 C2 O2"),
        "DA5": nsplit("C5' O4' C4'","C3' C2' C1'","N9 C4","C8 N7 C5","C6 N6 N1","C2 N3"),
        "DG5": nsplit("C5' O4' C4'","C3' C2' C1'","N9 C4","C8 N7 C5","C6 O6 N1","C2 N2 N3"),
        "DC5": nsplit("C5' O4' C4'","C3' C2' C1'","N1 C6","C5 C4 N4","N3 C2 O2"),
        "DT5": nsplit("C5' O4' C4'","C3' C2' C1'","N1 C6","C5 C4 O4 C7 C5M","N3 C2 O2"),
        }
    # Generic names for side chain beads
    cgnames = spl("BB SC1 SC2 SC3 SC4")
    # Generic names for DNA beads
    cgnames_dna = spl("BB1 BB2 BB3 SC1 SC2 SC3 SC4")
    cgnames_dnaT = spl("BB2 BB3 SC1 SC2 SC3 SC4")

class CGcafca:
    # cafemol, CA, no lipids
    mapping = {
        "ALA":  [['CA']],
        "CYS":  [['CA']],
        "ASP":  [['CA']],
        "GLU":  [['CA']],
        "PHE":  [['CA']],
        "GLY":  [['CA']],
        "HIS":  [['CA']],
        "HIH":  [['CA']],
        "ILE":  [['CA']],
        "LYS":  [['CA']],
        "LEU":  [['CA']],
        "MET":  [['CA']],
        "ASN":  [['CA']],
        "PRO":  [['CA']],
        "HYP":  [['CA']],
        "GLN":  [['CA']],
        "ARG":  [['CA']],
        "SER":  [['CA']],
        "THR":  [['CA']],
        "VAL":  [['CA']],
        "TRP":  [['CA']],
        "TYR":  [['CA']],
        #"POPE": phosphatydilethanolamine + palmitoyl1 + oleyl2,     ####
        #"DOPE": phosphatydilethanolamine + oleyl1     + oleyl2,     ####
        #"DPPE": phosphatydilethanolamine + palmitoyl1 + palmitoyl2, #### Lipid are not supported in cafemol
        #"POPG": phosphatidylglycerol     + palmitoyl1 + oleyl2,     ####
        #"DOPG": phosphatidylglycerol     + oleyl1     + oleyl2,     ####
        #"DPPG": phosphatidylglycerol     + palmitoyl1 + palmitoyl2, ####
        "DA": nsplit("P OP1 OP2 O5' O3' O1P O2P","C5' O4' C4' C3' C2' C1'","N9 C4 C8 N7 C5 C6 N6 N1 C2 N3"),
        "DG": nsplit("P OP1 OP2 O5' O3' O1P O2P","C5' O4' C4' C3' C2' C1'","N9 C4 C8 N7 C5 C6 O6 N1 C2 N2 N3"),
        "DC": nsplit("P OP1 OP2 O5' O3' O1P O2P","C5' O4' C4' C3' C2' C1'","N1 C6 C5 C4 N4 N3 C2 O2"),
        "DT": nsplit("P OP1 OP2 O5' O3' O1P O2P","C5' O4' C4' C3' C2' C1'","N1 C6 C5 C4 O4 C7 C5M N3 C2 O2"),
        "DA5": nsplit("C5' O4' C4' C3' C2' C1'","N9 C4 C8 N7 C5 C6 N6 N1 C2 N3"),
        "DG5": nsplit("C5' O4' C4' C3' C2' C1'","N9 C4 C8 N7 C5 C6 O6 N1 C2 N2 N3"),
        "DC5": nsplit("C5' O4' C4' C3' C2' C1'","N1 C6 C5 C4 N4 N3 C2 O2"),
        "DT5": nsplit("C5' O4' C4' C3' C2' C1'","N1 C6 C5 C4 O4 C7 C5M N3 C2 O2"),
        }
    # Generic names for side chain beads
    cgnames = spl("CA")
    # Generic names for DNA beads
    cgnames_dna = spl("O S N")
    cgnames_dnaT = spl("S N")

class CGcafcom:
    #cafemol, com of residue, no lipids
    bb        = "N CA C O H H1 H2 H3 OC1 OC2"     # jhpeng                                                              #@#  
    mapping = {
        "ALA":  joinall(bb +" CB"),
        "CYS":  joinall(bb,"CB SG"),
        "ASP":  joinall(bb,"CB CG OD1 OD2"),
        "GLU":  joinall(bb,"CB CG CD OE1 OE2"),
        "PHE":  joinall(bb,"CB CG CD1 HD1","CD2 HD2 CE2 HE2","CE1 HE1 CZ HZ"),
        "GLY":  joinall(bb),
        "HIS":  joinall(bb,"CB CG","CD2 HD2 NE2 HE2","ND1 HD1 CE1 HE1"),
        "HIH":  joinall(bb,"CB CG","CD2 HD2 NE2 HE2","ND1 HD1 CE1 HE1"),     # Charged Histidine.
        "ILE":  joinall(bb,"CB CG1 CG2 CD CD1"),
        "LYS":  joinall(bb,"CB CG CD","CE NZ HZ1 HZ2 HZ3"),
        "LEU":  joinall(bb,"CB CG CD1 CD2"),
        "MET":  joinall(bb,"CB CG SD CE"),
        "ASN":  joinall(bb,"CB CG ND1 ND2 OD1 OD2 HD11 HD12 HD21 HD22"),
        "PRO":  joinall(bb,"CB CG CD"),
        "HYP":  joinall(bb,"CB CG CD OD"),
        "GLN":  joinall(bb,"CB CG CD OE1 OE2 NE1 NE2 HE11 HE12 HE21 HE22"),
        "ARG":  joinall(bb,"CB CG CD","NE HE CZ NH1 NH2 HH11 HH12 HH21 HH22"),    
        "SER":  joinall(bb,"CB OG HG"),
        "THR":  joinall(bb,"CB OG1 HG1 CG2"),
        "VAL":  joinall(bb,"CB CG1 CG2"),
        "TRP":  joinall(bb,"CB CG CD2","CD1 HD1 NE1 HE1 CE2","CE3 HE3 CZ3 HZ3","CZ2 HZ2 CH2 HH2"),
        "TYR":  joinall(bb,"CB CG CD1 HD1","CD2 HD2 CE2 HE2","CE1 HE1 CZ OH HH"),
        "DA": nsplit("P OP1 OP2 O5' O3' O1P O2P","C5' O4' C4' C3' C2' C1'","N9 C4 C8 N7 C5 C6 N6 N1 C2 N3"),
        "DG": nsplit("P OP1 OP2 O5' O3' O1P O2P","C5' O4' C4' C3' C2' C1'","N9 C4 C8 N7 C5 C6 O6 N1 C2 N2 N3"),
        "DC": nsplit("P OP1 OP2 O5' O3' O1P O2P","C5' O4' C4' C3' C2' C1'","N1 C6 C5 C4 N4 N3 C2 O2"),
        "DT": nsplit("P OP1 OP2 O5' O3' O1P O2P","C5' O4' C4' C3' C2' C1'","N1 C6 C5 C4 O4 C7 C5M N3 C2 O2"),
        "DA5": nsplit("C5' O4' C4' C3' C2' C1'","N9 C4 C8 N7 C5 C6 N6 N1 C2 N3"),
        "DG5": nsplit("C5' O4' C4' C3' C2' C1'","N9 C4 C8 N7 C5 C6 O6 N1 C2 N2 N3"),
        "DC5": nsplit("C5' O4' C4' C3' C2' C1'","N1 C6 C5 C4 N4 N3 C2 O2"),
        "DT5": nsplit("C5' O4' C4' C3' C2' C1'","N1 C6 C5 C4 O4 C7 C5M N3 C2 O2"),
        }
    # Generic names for side chain beads
    cgnames = spl("CA")
    # Generic names for DNA beads
    cgnames_dna = spl("O S N")
    cgnames_dnaT = spl("S N")

class CGunres:
    # unres have no lipid and nucleic definition
    bb        = "N CA C O H H1 H2 H3 OC1 OC2"     # jhpeng                                                              #@#  
    mapping = {
        "ALA":  nsplit(bb,"CB"),
        "CYS":  nsplit(bb,"CB SG"),
        "ASP":  nsplit(bb,"CB CG OD1 OD2"),
        "GLU":  nsplit(bb,"CB CG CD OE1 OE2"),
        "PHE":  nsplit(bb,"CB CG CD1 HD1 CD2 HD2 CE2 HE2 CE1 HE1 CZ HZ"),
        "GLY":  nsplit(bb),
        "HIS":  nsplit(bb,"CB CG CD2 HD2 NE2 HE2 ND1 HD1 CE1 HE1"),
        "HIH":  nsplit(bb,"CB CG CD2 HD2 NE2 HE2 ND1 HD1 CE1 HE1"),     # Charged Histidine.
        "ILE":  nsplit(bb,"CB CG1 CG2 CD CD1"),
        "LYS":  nsplit(bb,"CB CG CD CE NZ HZ1 HZ2 HZ3"),
        "LEU":  nsplit(bb,"CB CG CD1 CD2"),
        "MET":  nsplit(bb,"CB CG SD CE"),
        "ASN":  nsplit(bb,"CB CG ND1 ND2 OD1 OD2 HD11 HD12 HD21 HD22"),
        "PRO":  nsplit(bb,"CB CG CD"),
        "HYP":  nsplit(bb,"CB CG CD OD"),
        "GLN":  nsplit(bb,"CB CG CD OE1 OE2 NE1 NE2 HE11 HE12 HE21 HE22"),
        "ARG":  nsplit(bb,"CB CG CD NE HE CZ NH1 NH2 HH11 HH12 HH21 HH22"),    
        "SER":  nsplit(bb,"CB OG HG"),
        "THR":  nsplit(bb,"CB OG1 HG1 CG2"),
        "VAL":  nsplit(bb,"CB CG1 CG2"),
        "TRP":  nsplit(bb,"CB CG CD2 CD1 HD1 NE1 HE1 CE2 CE3 HE3 CZ3 HZ3 CZ2 HZ2 CH2 HH2"),
        "TYR":  nsplit(bb,"CB CG CD1 HD1 CD2 HD2 CE2 HE2 CE1 HE1 CZ OH HH"),
        }
    # Generic names for side chain beads
    residue_bead_names = spl("CA CB")
    cgnames = residue_bead_names

class atom:
    def __init__(self,atid,atname,resi,coori):
        self.atid = atid
        self.atname = atname
        self.resi = resi
        self.coori = coori

from collections import OrderedDict
class atomslist:
    def __init__(self,pdbatoms):
        self.atoms = pdbatoms
        residues = []
        for atom in pdbatoms:
            residues += [atom.resi]
        resnames = list(OrderedDict.fromkeys(residues))
        resnums = len(resnames)
        atids = [[] for s in xrange(resnums)]
        atnames = [[] for s in xrange(resnums)]
        atcoors = [[] for s in xrange(resnums)]
        names = []
        for atom in pdbatoms:
            residuesj = atom.resi
            j = resnames.index(residuesj)
            atids[j] += [atom.atid]
            atnames[j] += [atom.atname]
            names += [atom.atname]
            atcoors[j] += [atom.coori]
        self.resnames = resnames 
        self.residues = residues 
        self.atids = atids 
        self.atnames = atnames 
        self.atcoors = atcoors 
        self.names = names

def readpdb(inpdb):
    ##01234567890123456789012345678901234567890123456789012345678901234567890123456789
    ##ATOM   2155 HH11 ARG C 203     116.140  48.800   6.280  1.00  0.00
    try:
        fpdb = open(inpdb,'r')
    except:
        IOError
        sys.exit("cgmap error>> input PDB file not exist!")
    lines = fpdb.readlines()
    fpdb.close()
    pdbatoms = []
    nat = 0
    coords = []
    for i in range(len(lines)):
        if lines[i].startswith('ATOM '):
            atid = int(lines[i][6:11])
            atname = lines[i][12:16].strip()
            resname = lines[i][17:21].strip()
            chain = lines[i][21]
            resi = "/".join([resname,lines[i][22:26].strip(),chain])
            x = lines[i][30:38]
            y = lines[i][38:46]
            z = lines[i][46:54]
            coori = [float(s) for s in [x,y,z]]
            atomi = atom(atid,atname,resi,coori)
            coords += [coori]
            pdbatoms += [atomi]
            nat += 1
    if nat == 0:
        sys.exit("cgmap error>> input PDB file not correct!")
    return atomslist(pdbatoms),nat

def readcoords(inpdb):
    ##01234567890123456789012345678901234567890123456789012345678901234567890123456789
    ##ATOM   2155 HH11 ARG C 203     116.140  48.800   6.280  1.00  0.00
    try:
        fpdb = open(inpdb,'r')
    except:
        IOError
        sys.exit("cgmap error>> input PDB file not exist!")
    lines = fpdb.readlines()
    fpdb.close()
    coords = []
    n = 0
    for i in range(len(lines)):
        if lines[i].startswith('ATOM '):
            x = lines[i][30:38]
            y = lines[i][38:46]
            z = lines[i][46:54]
            coori = [float(s) for s in [x,y,z]]
            coords += [coori]
            n += 1
    if n == 0:
        sys.exit("cgmap error>> input PDB file not correct!")
    return coords

def mapcg(resname,atoms,CGmapping):
    mapping = CGmapping.mapping
    resnamek = resname.split('/')[0]
    if resnamek in Protein+Lipid:
        cgnames = CGmapping.cgnames
    elif resnamek in DNA:
        cgnames = CGmapping.cgnames_dna
    elif resnamek in DNAT:
        cgnames = CGmapping.cgnames_dnaT
    else:
        sys.exit("cgmap_error>> checking the residue: %s"%resname)
    resi = atoms.resnames.index(resname)
    #print mapping[resnamek]
    beadnum = len(mapping[resnamek])
    names = mapping[resnamek]

    beadsi = [[] for s in xrange(beadnum)]
    beadsielem = [[] for s in xrange(beadnum)]
    coors = [[] for s in xrange(beadnum)]
    beadsicoors = [[] for s in xrange(beadnum)]
    beadsimass = [0 for s in xrange(beadnum)]
    beadsinames = cgnames[0:beadnum]
    beadsresi = [resname for s in xrange(beadnum)]
    atnames = atoms.atnames[resi]
    atids = atoms.atids[resi]
    atcoors = atoms.atcoors[resi]
    dict_atn_id = hash(atnames,atids)
    dict_atn_coor = hash(atnames,atcoors)
    for i in xrange(beadnum):
        for atname in atnames:
            if atname in names[i]:
                massati = mass[get_elem(atname)]
                beadsi[i] += [dict_atn_id[atname]]
                coors[i] += [np.array(dict_atn_coor[atname])*massati]
                beadsimass[i] += massati
                beadsielem[i] += [atname]
    for i in xrange(beadnum):
        beadsicoors[i] = sum(coors[i])/beadsimass[i]
    return beadsi,beadsielem,beadsinames,beadsicoors,beadsresi,beadsimass

if __name__ == '__main__':

    options = getoption(sys.argv[1:])        

    sys.stdout.write("cgmap_infor>> CGMAP for BayesTMD, Version %s\n"%version)
    #######################
    ### Parsing options ###
    #######################
    fpdbi = options['-f']
    fpdbo = options['-x']
    fnmap = options['-nmap']
    method = options['-ff']
    ftarx = options['-s']
    pdev = options['-p']
    ftmi = options['-o']

    atoms,nat = readpdb(fpdbi)
    sys.stdout.write("cgmap_infor>> Read in %5d atoms from %s\n"%(nat,fpdbi))
    beads = []
    beadselem = []
    beadsname = []
    beadscoor = []
    beadsres = []
    beadsmass = []

    ############################
    ### Main Mapping process ###
    ############################
    CGmapping = get_cgmapping(method)
    ##01234567890123456789012345678901234567890123456789012345678901234567890123456789
    ##ATOM   2155 HH11 ARG C 203     116.140  48.800   6.280  1.00  0.00
    pdbAtomLine = "ATOM  %5d %4s%4s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f\n"        
    for resi in atoms.resnames:
        beadsi, beadsielem,beadsinames,beadsicoor,beadsresi,beadsimass = mapcg(resi,atoms,CGmapping)
        beads += beadsi
        beadselem += beadsielem
        beadsname += beadsinames
        beadscoor += beadsicoor
        beadsres += beadsresi
        beadsmass += beadsimass
    ncg = len(beadsres)
    sys.stdout.write("cgmap_infor>> Generated %5d CG beads for %s\n"%(ncg,fpdbo))

    #########################
    ### write CG pdb file ###
    #########################
    pdbo = open(fpdbo,'w')
    for i in range(ncg):
        resname,resid,chid = beadsres[i].split('/')
        resid = int(resid)
        name = beadsname[i]
        x,y,z = beadscoor[i]
        alt = ' '
        oc = 1.0
        bf = 0.0
        pdbo.write(pdbAtomLine%(i+1,name,resname,chid,resid,alt,x,y,z,oc,bf))
    pdbo.close()

    #######################
    ### write map index ###
    #######################
### we assume that CG backbone are more accurate than side-chains ###
    bb = "O CA BB BB1 BB2 BB3"

    pnmap = open(fnmap,'w')
    aa2cg = []
    for i in range(ncg):
        aa2cgi = []
        for j in range(len(beads[i])):
            aa2cgi += [[beads[i][j],i]]
        aa2cg += aa2cgi
    aa2cg = dict(sorted(aa2cg,key = lambda x : x[0]))
    i = 0

    aa2cgf = []
    while i < nat:
        if i+1 in aa2cg.keys():
            pnmap.write("%5d %5d\n"%(i+1,aa2cg[i+1]+1))
            aa2cgf += [[i+1,aa2cg[i+1]+1]]
        else:
            pnmap.write("%5d %5d\n"%(i+1,0))
            aa2cgf += [[i+1,0]]
        i += 1
    pnmap.close()
    sys.stdout.write("cgmap_infor>> finished with %s file\n"%fnmap)

    ##########################################################
    ### read Target CG PDB and write .tmi file for BayesTMD###
    ##########################################################
    if options['-s'] != None:
        coorstar = np.array(readcoords(ftarx))
        sys.stdout.write("cgmap_infor>> Target CG PDB detected\n")
        inirms = get_inirms(beadscoor,coorstar,beadsmass,ncg)/10
        sys.stdout.write("cgmap_infor>> RMSD between the Mapped CG PDB and the target is %5.3f nm\n"%(inirms))
        if options['-o'] == None:
            sys.exit("cgmap_error>> please specify an output for .tmi file")

        ####### now write .tmi file for BayesTMD ####
        tmptmi = [
        '#Total Nr. of atoms, CG sites and distances for TMD\n',
        '#Force constant, initial rms, and end rms\n',
        '#Freq. to do pbc\n',
        '#Residue index, atom index and mass\n',
        '#Target CG sites\n',
        '# Initial AA2CG mapping, backbone or not?\n',
        '# reference x in order to remove pbc in Rij calculation\n']

        ks = 1.0
        endrms = 0.0
        tmdfrq = 1
        ptmi = open(ftmi,'w')
        ####### first nat, ncg and ncg*(ncg-1)/2 ####
        ptmi.write(tmptmi[0])
        if options['-p'] != None:
            ptmi.write("%5d %5d %5.3f\n"%(nat,ncg,float(pdev)))
        else: 
            ptmi.write("%5d %5d %5.3f\n"%(nat,ncg,1.0))
        ptmi.write(tmptmi[1])
        #ptmi.write("%5.3f %5.3f %5.3f\n"%(ks,inirms,0.0))
        ptmi.write("tmdk inirms\n")
        ptmi.write(tmptmi[2])
        ptmi.write("%5d\n"%tmdfrq)
        ptmi.write(tmptmi[3])
        for i in range(nat):
            resid = int(atoms.residues[i].split('/')[1])
            massi = mass[get_elem(atoms.names[i])]
            ptmi.write("%5d %5d  %10.3f\n"%(resid,i+1,massi))
        ptmi.write(tmptmi[4])
        for i in range(ncg):
            x,y,z = coorstar[i]/10
            if beadsname[i] in bb:
                ptmi.write(" %8.3f %8.3f %8.3f %8d\n"%(x,y,z,1))
            else:
                ptmi.write(" %8.3f %8.3f %8.3f %8d\n"%(x,y,z,0))
        ptmi.write(tmptmi[5])
        for i in range(nat):
            ptmi.write("%5d %5d\n"%(aa2cgf[i][0],aa2cgf[i][1]))
        ptmi.write(tmptmi[6])
        xini = np.array(readcoords(fpdbi))
        for i in range(nat):
            x,y,z = xini[i]/10
            ptmi.write(" %8.3f %8.3f %8.3f\n"%(x,y,z))
        ptmi.close()
        sys.stdout.write("cgmap_infor>> finished with %s file\n"%ftmi)
