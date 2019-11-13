"""
Copyright (C) 2019 Neelesh Soni <neelesh.soni@alumni.iiserpune.ac.in>, <neeleshsoni03@gmail.com>
This file is part of AIMorpher.

AIMorpher is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

AIMorpher is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with AIMorpher.  If not, see <http://www.gnu.org/licenses/>.
"""

from NMalignment import global_align

#matrix_file = './input/BLOSUM62.txt'
#AIMORPH matrix is used to get the exact alignment
matrix_file = './input/AIMORPH.txt'

def write_alignment(alignment,outfilename):
	outf = open(outfilename+'_output_alignment.txt','w')
	outf.write('\n>Sequence 1\n')
	outf.write(alignment[0])
	outf.write('\n')
	outf.write('\n>Sequence 2\n')
	outf.write(alignment[1])
	outf.close()
	return

def align_structures(mdlI,mdlF,outfilename):
	#seq_X gives sequence of mdlX and seq_X_map gives iter_length of the mdlX
	seq_I,seq_I_map = mdlI.sequence_map()
	seq_F,seq_F_map = mdlF.sequence_map()
	
	print ("\nResidues in Protein 1: ",len(seq_I))
	print ("Residues in Protein 2: ",len(seq_F))

	print("\nAligning Protein Sequences...")
	alignment = global_align(seq_I, seq_F, matrix=matrix_file)
	
	write_alignment(alignment,outfilename)

	#common sequence iter_length
	seq_I_itermap=[]
	seq_F_itermap=[]
	seqI_counter=-1
	seqF_counter=-1

	#print(alignment[0])
	#print(alignment[1])

	#find out common amino acid sequence
	for i in range(0,len(alignment[0])):
		aa1 = alignment[0][i];
		aa2 = alignment[1][i];
		
		#increase counter if not '-'
		if aa1!='-':
			seqI_counter+=1
		
		#increase counter if not '-'
		if aa2!='-':
			seqF_counter+=1

		#append iter_length variable value in map if amino acids are aligned based on substitution matrix
		if ((aa1!='-') and (aa2!='-')):
			
			#Update in sequence 1. seq_I_map[seqI_counter][0] provides the start of the residue [1] provides end
			iterl_I = seq_I_map[seqI_counter]
			seq_I_itermap.append(iterl_I);

			#update in sequence 2. [0] provides the start of the residue [1] provides end
			iterl_F = seq_F_map[seqF_counter]
			seq_F_itermap.append(iterl_F);

	print ("Residues aligned between Protein 1 and Protein 2: ",len(seq_I_itermap),len(seq_F_itermap))

	return seq_I_itermap,seq_F_itermap