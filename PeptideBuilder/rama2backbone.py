#!/usr/bin/env python
docstring='''
rama2backbone.py seq.fasta rama.txt out.pdb
    convert amino acid sequence "seq.fasta" and ramachandran torsion
    angle table "rama.txt" to backbone-only peptide chain "out.pdb"
    rama.txt is a three column table, with the three columns 
    corresponding to omega, phi, psi respectively.

rama2backbone.py seq.fasta out.pdb
    convert amino acid sequence "seq.fasta" to extended backbone-only
    peptide chain "out.pdb"
'''

from Bio.PDB import is_aa,PDBIO
from Bio.PDB.Atom import Atom,Vector
from Bio.PDB.Residue import Residue
from Bio.PDB.Chain import Chain
from Bio.PDB.Model import Model
from Bio.PDB.Structure import Structure
from Bio.PDB.Vector import calc_dihedral,rotaxis

import math, warnings
import numpy as np
import sys

aa1to3 = { # 3 letter to 1 letter amino acid code conversion
    'A':'ALA', 'V':'VAL', 'F':'PHE', 'P':'PRO', 'M':'MET',
    'I':'ILE', 'L':'LEU', 'D':'ASP', 'E':'GLU', 'K':'LYS',
    'R':'ARG', 'S':'SER', 'T':'THR', 'Y':'TYR', 'H':'HIS',
    'C':'CYS', 'N':'ASN', 'Q':'GLN', 'W':'TRP', 'G':'GLY',
    'B':'ASX', 'Z':'GLX', 'U':'SEC', 'O':'PYL', 'X':'UNK', 'J':'UNK',
    }

class Geo():
    '''Geometry base class'''
    def __init__(self,residue_name='G'):
        '''Geometry of Glycine'''
        self.CA_N_length=1.46
        self.CA_C_length=1.52
        self.N_CA_C_angle=111.

        self.C_O_length=1.23
        self.CA_C_O_angle=120.5

        ## This dihedral angle is residue specific
        self.N_CA_C_O_diangle= 120.0 # LTRKDEQMHFYW
        if residue_name in "G":
            self.N_CA_C_O_diangle= 180.0
        elif residue_name in "A":
            self.N_CA_C_O_diangle=-60.5
        elif residue_name in "SCUVIN":
            self.N_CA_C_O_diangle= -60.0
        elif residue_name in "P":
            self.N_CA_C_O_diangle=-45.0

        self.phi=-120
        self.psi_im1=140
        self.omega=180.0
        self.peptide_bond=1.33
        self.CA_C_N_angle =116.642992978143
        self.C_N_CA_angle= 121.382215820277

        self.residue_name= residue_name

    def __repr__(self):
        return self.residue_name

def calculateCoordinates(refA, refB, refC, L, ang, di):
    AV=refA.get_vector()
    BV=refB.get_vector()
    CV=refC.get_vector()
    
    CA=AV-CV
    CB=BV-CV

    ##CA vector
    AX=CA[0]
    AY=CA[1]
    AZ=CA[2]

    ##CB vector
    BX=CB[0]
    BY=CB[1]
    BZ=CB[2]

    ##Plane Parameters
    A=(AY*BZ)-(AZ*BY)
    B=(AZ*BX)-(AX*BZ)
    G=(AX*BY)-(AY*BX)

    ##Dot Product Constant
    F= math.sqrt(BX*BX + BY*BY + BZ*BZ) * L * math.cos(ang*(math.pi/180.0))

    ##Constants
    const=math.sqrt( math.pow((B*BZ-BY*G),2) *(-(F*F)*(A*A+B*B+G*G)+(B*B*(BX*BX+BZ*BZ) + A*A*(BY*BY+BZ*BZ)- (2*A*BX*BZ*G) + (BX*BX+ BY*BY)*G*G - (2*B*BY)*(A*BX+BZ*G))*L*L))
    denom= (B*B)*(BX*BX+BZ*BZ)+ (A*A)*(BY*BY+BZ*BZ) - (2*A*BX*BZ*G) + (BX*BX+BY*BY)*(G*G) - (2*B*BY)*(A*BX+BZ*G)

    X= ((B*B*BX*F)-(A*B*BY*F)+(F*G)*(-A*BZ+BX*G)+const)/denom

    if((B==0 or BZ==0) and (BY==0 or G==0)):
        const1=math.sqrt( G*G*(-A*A*X*X+(B*B+G*G)*(L-X)*(L+X)))
        Y= ((-A*B*X)+const1)/(B*B+G*G)
        Z= -(A*G*G*X+B*const1)/(G*(B*B+G*G))
    else:
        Y= ((A*A*BY*F)*(B*BZ-BY*G)+ G*( -F*math.pow(B*BZ-BY*G,2) + BX*const) - A*( B*B*BX*BZ*F- B*BX*BY*F*G + BZ*const)) / ((B*BZ-BY*G)*denom)
        Z= ((A*A*BZ*F)*(B*BZ-BY*G) + (B*F)*math.pow(B*BZ-BY*G,2) + (A*BX*F*G)*(-B*BZ+BY*G) - B*BX*const + A*BY*const) / ((B*BZ-BY*G)*denom)

    
    #GET THE NEW VECTOR from the orgin
    D=Vector(X, Y, Z) + CV
    with warnings.catch_warnings():
        # ignore inconsequential warning
        warnings.simplefilter("ignore")
        temp=calc_dihedral(AV, BV, CV, D)*(180.0/math.pi)
    
  
    di=di-temp
    rot= rotaxis(math.pi*(di/180.0), CV-BV)
    D=(D-BV).left_multiply(rot)+BV
    
    return D.get_array()

def makeRes(segID, N, CA, C, O, geo):
    '''Creates a Glycine residue'''
    ##Create Residue Data Structure
    res= Residue((' ', segID, ' '), aa1to3[geo.residue_name], '    ')

    res.add(N)
    res.add(CA)
    res.add(C)
    res.add(O)

    return res

def initialize_res(residue):
    '''Creates a new structure containing a single amino acid. The type and
    geometry of the amino acid are determined by the argument, which has to be
    either a geometry object or a single-letter amino acid code.
    The amino acid will be placed into chain A of model 0.'''
    
    if isinstance( residue, Geo ):
        geo = residue
    else:
        geo= Geo(residue) 
    
    segID=1
    AA= geo.residue_name
    CA_N_length=geo.CA_N_length
    CA_C_length=geo.CA_C_length
    N_CA_C_angle=geo.N_CA_C_angle
    
    CA_coord= np.array([0.,0.,0.])
    C_coord= np.array([CA_C_length,0,0])
    N_coord = np.array([CA_N_length*math.cos(N_CA_C_angle*(math.pi/180.0)),CA_N_length*math.sin(N_CA_C_angle*(math.pi/180.0)),0])

    N= Atom("N", N_coord, 0.0 , 1.0, " "," N", 0, "N")
    CA=Atom("CA", CA_coord, 0.0 , 1.0, " "," CA", 0,"C")
    C= Atom("C", C_coord, 0.0, 1.0, " ", " C",0,"C")

    ##Create Carbonyl atom (to be moved later)
    C_O_length=geo.C_O_length
    CA_C_O_angle=geo.CA_C_O_angle
    N_CA_C_O_diangle=geo.N_CA_C_O_diangle
    
    carbonyl=calculateCoordinates(N, CA, C, C_O_length, CA_C_O_angle, N_CA_C_O_diangle)
    O= Atom("O",carbonyl , 0.0 , 1.0, " "," O", 0, "O")

    res=makeRes(segID, N, CA, C, O, geo)

    cha= Chain('A')
    cha.add(res)
    
    mod= Model(0)
    mod.add(cha)

    struc= Structure('X')
    struc.add(mod)
    return struc

def getReferenceResidue(structure):
    '''Returns the last residue of chain A model 0 of the given structure.
    
    This function is a helper function that should not normally be called
    directly.'''

    # If the following line doesn't work we're in trouble.
    # Likely initialize_res() wasn't called.
    resRef = structure[0]['A'].child_list[-1]
    
    # If the residue is not an amino acid we're in trouble.
    # Likely somebody is trying to append residues to an existing
    # structure that has non-amino-acid molecules in the chain.
    assert is_aa(resRef)
        
    return resRef

def add_residue_from_geo(structure, geo):
    '''Adds a residue to chain A model 0 of the given structure, and
    returns the new structure. The residue to be added is determined by
    the geometry object given as second argument.
    
    This function is a helper function and should not normally be called
    directly. Call add_residue() instead.'''
    resRef= getReferenceResidue(structure)
    AA=geo.residue_name
    segID= resRef.get_id()[1]
    segID+=1

    ##geometry to bring together residue
    peptide_bond=geo.peptide_bond
    CA_C_N_angle=geo.CA_C_N_angle
    C_N_CA_angle=geo.C_N_CA_angle

    ##Backbone Coordinages
    N_CA_C_angle=geo.N_CA_C_angle
    CA_N_length=geo.CA_N_length
    CA_C_length=geo.CA_C_length
    phi= geo.phi
    psi_im1=geo.psi_im1
    omega=geo.omega

    N_coord=calculateCoordinates(resRef['N'], resRef['CA'], resRef['C'], peptide_bond, CA_C_N_angle, psi_im1)
    N= Atom("N", N_coord, 0.0 , 1.0, " "," N", 0, "N")

    CA_coord=calculateCoordinates(resRef['CA'], resRef['C'], N, CA_N_length, C_N_CA_angle, omega)
    CA=Atom("CA", CA_coord, 0.0 , 1.0, " "," CA", 0,"C")

    C_coord=calculateCoordinates(resRef['C'], N, CA, CA_C_length, N_CA_C_angle, phi)
    C= Atom("C", C_coord, 0.0, 1.0, " ", " C",0,"C")

    ##Create Carbonyl atom (to be moved later)
    C_O_length=geo.C_O_length
    CA_C_O_angle=geo.CA_C_O_angle
    N_CA_C_O_diangle=geo.N_CA_C_O_diangle

    carbonyl=calculateCoordinates(N, CA, C, C_O_length, CA_C_O_angle, N_CA_C_O_diangle)
    O= Atom("O",carbonyl , 0.0 , 1.0, " "," O", 0, "O")
    
    res=makeRes(segID, N, CA, C, O, geo)
        
    resRef['O'].set_coord(calculateCoordinates(res['N'], resRef['CA'], resRef['C'], C_O_length, CA_C_O_angle, 180.0))

    ghost= Atom("N", calculateCoordinates(res['N'], res['CA'], res['C'], peptide_bond, CA_C_N_angle, psi_im1), 0.0 , 0.0, " ","N", 0, "N")
    res['O'].set_coord(calculateCoordinates( res['N'], res['CA'], res['C'], C_O_length, CA_C_O_angle, 180.0))

    structure[0]['A'].add(res)
    return structure


def make_extended_structure(AA_chain):
    '''Place a sequence of amino acids into a peptide in the extended
    conformation. The argument AA_chain holds the sequence of amino
    acids to be used.'''
    geo = Geo(AA_chain[0])
    struc=initialize_res(geo)
    
    for i in range(1,len(AA_chain)): 
        AA = AA_chain[i]
        geo = Geo(AA)
        add_residue(struc, geo)

    return struc
    
def add_residue(structure, residue, phi=-120, psi_im1=140, omega=-370):
    '''Adds a residue to chain A model 0 of the given structure, and
    returns the new structure. The residue to be added can be specified
    in two ways: either as a geometry object (in which case
    the remaining arguments phi, psi_im1, and omega are ignored) or as a
    single-letter amino-acid code. In the latter case, the optional
    arguments phi, psi_im1, and omega specify the corresponding backbone
    angles.
    
    When omega is specified, it needs to be a value greater than or equal
    to -360. Values below -360 are ignored.''' 
    
    if isinstance( residue, Geo ):
        geo = residue
    else:
        geo=Geo(residue) 
        geo.phi=phi
        geo.psi_im1=psi_im1
        if omega>-361:
            geo.omega=omega
    
    add_residue_from_geo(structure, geo)
    return structure
    
def make_structure(AA_chain,phi,psi_im1,omega=[]):
    '''Place a sequence of amino acids into a peptide with specified
    backbone dihedral angles. The argument AA_chain holds the
    sequence of amino acids to be used. The arguments phi and psi_im1 hold
    lists of backbone angles, one for each amino acid, *starting from
    the second amino acid in the chain*. The argument 
    omega (optional) holds a list of omega angles, also starting from
    the second amino acid in the chain.'''
    geo = Geo(AA_chain[0])
    struc=initialize_res(geo)

    if len(omega)==0:
        for i in range(1,len(AA_chain)): 
            AA = AA_chain[i]
            add_residue(struc, AA, phi[i-1], psi_im1[i-1])
    else:
        for i in range(1,len(AA_chain)): 
            AA = AA_chain[i]
            add_residue(struc, AA, phi[i-1], psi_im1[i-1], omega[i-1])

    return struc
    
def make_structure_from_geos(geos):
    '''Creates a structure out of a list of geometry objects.'''
    model_structure=initialize_res(geos[0])
    for i in range(1,len(geos)):
        model_structure=add_residue(model_structure, geos[i])

    return model_structure

if __name__=="__main__":
    if not len(sys.argv) in [3,4]:
        sys.stderr.write(docstring)
        exit()
    
    AA_chain=''
    fp=open(sys.argv[1],'rU')
    for line in fp.read().splitlines():
        if not line.startswith('>'):
            AA_chain+=line.strip()
    fp.close()
    if not len(AA_chain):
        sys.stderr.write("ERROR! Empty sequence\n")
        exit()

    if len(sys.argv)==3:
        struc=make_extended_structure(AA_chain)
    else:
        rama_table=np.loadtxt(sys.argv[2])
        if len(AA_chain)!=len(rama_table):
            sys.stderr.write("ERROR! Length of sequence length %d"%len(AA_chain
                )+" is not the same as torsion angle table length %d\n"%len(
                rama_table))
            exit()

        phi=rama_table[1:,0].tolist()
        psi_im1=rama_table[:,1].tolist()
        omega=[]
        if rama_table.shape[1]>2:
            omega=rama_table[:,0].tolist()
            phi=rama_table[1:,1].tolist()
            psi_im1=rama_table[:,2].tolist()
        struc=make_structure(AA_chain,phi=phi,psi_im1=psi_im1,omega=omega)
    
    out = PDBIO()
    out.set_structure(struc)
    out.save( sys.argv[-1] )
