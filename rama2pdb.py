#!/usr/bin/python
docstring='''rama2pdb.py seq.fasta rama.txt out.pdb
    convert amino acid sequence "seq.fasta" and ramachandran torsion
    angle table "rama.txt" to poly-glycine peptide chain "out.pdb"
'''

from PeptideBuilder import Geometry
import PeptideBuilder
import Bio.PDB
import numpy as np
import sys

def readSingleSequenceFasta(infile="seq.fasta"):
    fp=open(infile,'rU')
    txt=fp.read()
    fp.close()
    sequence=''
    for line in txt.splitlines():
        if line.startswith('>') or not line.strip():
            continue
        sequence+=line.strip()
    return sequence

if __name__=="__main__":
    if len(sys.argv)!=4:
        sys.stderr.write(docstring)
        exit()
    
    sequence=readSingleSequenceFasta(sys.argv[1])
    rama_table=np.loadtxt(sys.argv[2])

    if len(sequence)!=len(rama_table):
        sys.stderr.write("ERROR! Length of sequence %d is not the same as length of ramachandran angle table %d\n"%(len(sequence),len(rama_table)))
        exit()
    if not len(sequence):
        sys.stderr.write("ERROR! Empty sequence\n")
        exit()

    geo = Geometry.geometry(sequence[0])
    geo.phi=rama_table[0][0]
    geo.psi_im1=rama_table[0][1]
    structure = PeptideBuilder.initialize_res(geo)

    for i in range(1,len(sequence)):
        geo = Geometry.geometry(sequence[i])
        geo.phi=rama_table[i][0]
        geo.psi_im1=rama_table[i-1][1]
        structure = PeptideBuilder.add_residue(structure, geo)
    
    out = Bio.PDB.PDBIO()
    out.set_structure(structure)
    out.save( sys.argv[3] )
