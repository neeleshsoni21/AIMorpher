"""
Copyright (C) 2019 Neelesh Soni <neelesh.soni@alumni.iiserpune.ac.in>, <neeleshsoni03@gmail.com>
This file is part of AIMorpher.

AIMorpher is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

AIMorpher is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with AIMorpher.  If not, see <http://www.gnu.org/licenses/>.
"""

from addoptions import addoptions
from pdb_class import PDB
from alignment import align_structures
from transformation import morph_structure, morph_structure2

#Parsing the inputs from the command line#
arguments=addoptions()

fname1 = arguments[0]
#Load Initial state
print ("\nReading Protein 1 PDB file...")
mdlI = PDB(fname1)

#Load final state
fname2 = arguments[1]
print ("Reading Protein 2 PDB file...")
mdlF = PDB(fname2)

#Output file name
outfilename = arguments[2]

#Align the two structures
#Determine the common residues/atoms and get temperory structures
seq_I_itermap, seq_F_itermap = align_structures(mdlI,mdlF,outfilename)

#read the morphing parameters
num_intm = arguments[3]

print("\nMorphing the Aligned Residues")

#apply recursive transformations
#Temperorly deactivating similar sequence/structure morphing. 
#Morphing on non-identical sequences giving non-real morphed structures due to unusal bond lengths
#morph_structure(mdlI,seq_I_itermap,mdlF,seq_F_itermap,num_intm,outfilename)

morph_structure2(mdlI,seq_I_itermap,mdlF,seq_F_itermap,num_intm,outfilename)