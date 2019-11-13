
"""
Copyright (C) 2019 Neelesh Soni <neelesh.soni@alumni.iiserpune.ac.in>, <neeleshsoni03@gmail.com>
This file is part of AIMorpher.

AIMorpher is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

AIMorpher is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with AIMorpher.  If not, see <http://www.gnu.org/licenses/>.
"""

from copy import deepcopy

'''
This module does the linear transformation on all atoms independently
'''
def morph_structure2(mdlI,seq_I_itermap,mdlF,seq_F_itermap,num_intm,outfilename):

	#Open filename for writing
	fout = open(outfilename,'w')
	fout.close()

	#Copy the class object
	mdlM = deepcopy(mdlI)

	#Iterate over number of iterations
	for iteration in range(0,num_intm+1):

		#dist = 1.0 / num_intm;
		distfrac = float(iteration) / num_intm;

		output_selection = []

		#Iterate over I itermap. seq_I_itermap[i] contains iter_idx of residue i
		for i in range(0,len(seq_I_itermap)):
			itrI,itrF = seq_I_itermap[i][0],seq_F_itermap[i][0]
			
			curr_resno = mdlI.resSeq()[itrI]
			new_resno = curr_resno
			#residue atoms counter
			j1 = itrI
			j2 = itrF

			#make sure all atoms of residue moves together
			while curr_resno==new_resno:
				
				x1,y1,z1 = mdlI.x()[j1],mdlI.y()[j1],mdlI.z()[j1]
				x2,y2,z2 = mdlF.x()[j2],mdlF.y()[j2],mdlF.z()[j2]

				mdlM.x()[j1] = mdlI.x()[j1] + (x2-x1)*distfrac 
				mdlM.y()[j1] = mdlI.y()[j1] + (y2-y1)*distfrac 
				mdlM.z()[j1] = mdlI.z()[j1] + (z2-z1)*distfrac 

				output_selection.append(j1)
				j1 = j1 + 1
				j2 = j2 + 1

				#print(len(mdlI.resSeq()),j1)
				#new_resno = mdlI.resSeq()[j1]
				#for last residue j1 can exceed mdlI.resSeq() length. thus try statement
				try:
					new_resno = mdlI.resSeq()[j1]
				except:
					new_resno = -100000
			
		#Write the models
		mdlM.write_selection_multimdl(outfilename,output_selection,str(iteration))

	return

'''
This module does the linear transformation with each residue as block
'''
def morph_structure(mdlI,seq_I_itermap,mdlF,seq_F_itermap,num_intm,outfilename):

	#Open filename for writing
	fout = open(outfilename,'w')
	fout.close()

	#Copy the class object
	mdlM = deepcopy(mdlI)

	#Iterate over number of iterations
	for iteration in range(0,num_intm+1):

		#dist = 1.0 / num_intm;
		distfrac = float(iteration) / num_intm;

		output_selection = []

		#Iterate over I itermap. seq_I_itermap[i] contains iter_idx of residue i
		for i in range(0,len(seq_I_itermap)):
			
			#residue iter counters. j1s and j1e provides residue start and end iter index.
			j1s,j1e,j2s,j2e = seq_I_itermap[i][0],seq_I_itermap[i][1],seq_F_itermap[i][0],seq_F_itermap[i][1]

			#print(i,j1s,j1e,j2s,j2e)

			#Get the average coordinate of ith residues in the Initial structure
			xcoord1=[];ycoord1=[];zcoord1=[]
			for j in range(j1s,j1e+1):
				
				x1,y1,z1 = mdlI.x()[j],mdlI.y()[j],mdlI.z()[j]
				
				xcoord1.append(x1);ycoord1.append(y1);zcoord1.append(z1);
			
			#determine the average	
			avg_xcoord1 = sum(xcoord1)/float(len(xcoord1))
			avg_ycoord1 = sum(ycoord1)/float(len(ycoord1))
			avg_zcoord1 = sum(zcoord1)/float(len(zcoord1))

			#Get the average coordinate of ith residues in the Final structure
			xcoord2=[];ycoord2=[];zcoord2=[]
			for j in range(j2s,j2e+1):
				
				x2,y2,z2 = mdlF.x()[j],mdlF.y()[j],mdlF.z()[j]
				
				xcoord2.append(x2);ycoord2.append(y2);zcoord2.append(z2);
			
			#determine the average
			avg_xcoord2 = sum(xcoord2)/float(len(xcoord2))
			avg_ycoord2 = sum(ycoord2)/float(len(ycoord2))
			avg_zcoord2 = sum(zcoord2)/float(len(zcoord2))

			#Get the distance fraction in each direction (x,y,z) to move for this iteration
			xtrans = (avg_xcoord2 - avg_xcoord1)*distfrac
			ytrans = (avg_ycoord2 - avg_ycoord1)*distfrac
			ztrans = (avg_zcoord2 - avg_zcoord1)*distfrac

			#Move the coordinates in each directions
			for j in range(j1s,j1e+1):
				x1,y1,z1 = mdlI.x()[j],mdlI.y()[j],mdlI.z()[j]
				mdlM.x()[j] = x1 + xtrans;
				mdlM.y()[j] = y1 + ytrans;
				mdlM.z()[j] = z1 + ztrans;

				#Record the iter index within model M for writing later
				output_selection.append(j)
			
		#Write the models
		mdlM.write_selection_multimdl(outfilename,output_selection,str(iteration))
		
	return



