"""
Copyright (C) 2019 Neelesh Soni <neelesh.soni@alumni.iiserpune.ac.in>, <neeleshsoni03@gmail.com>
This file is part of AIMorpher.

AIMorpher is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

AIMorpher is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with AIMorpher.  If not, see <http://www.gnu.org/licenses/>.
"""

from os import system

#dictionary for converting three letter code to one letter code
THREEONE={'LEU':'L','PHE':'F','TRP':'W','TYR':'Y','ASN':'N','ASP':'D','HIS':'H','VAL':'V','ILE':'I',
'PRO':'P','THR':'T','GLY':'G','ALA':'A','SER':'S','CYS':'C','GLN':'Q','GLU':'E','MET':'M','LYS':'K','ARG':'R'}

#Class definition
class PDB:

	#init function for parsing the pdb format file
	def __init__(self,filename):
		
		#open the file for reading
		pdbfile = open(filename,'r')
		pdblines = pdbfile.readlines()
		pdbfile.close();

		#declaring variables as list for storing values of each column in the pdb file
		self.__serial = []
		self.__name = []
		self.__altLoc = []
		self.__resName = []
		self.__chainID = []
		self.__resSeq  = []
		self.__iCode = []
		self.__x = []
		self.__y = []
		self.__z = []

		#Ierating over lines in the pdb file
		for line in pdblines:				
			atom = line[0:6]

			#Skip all HETATM/ANISOU etc
			if atom.strip() != 'ATOM':
				continue
			
			#Skip all Hydrogens
			element = line[76:78]
			if element.strip() == "H":
				continue

			#Skip all iCodes
			iCode = line[26:27]
			if iCode.strip() != "":
				continue

			#Skip all alternate locations
			altLoc = line[16:17]
			if (altLoc.strip() != "") and (altLoc != "A") and (altLoc != 1):
				continue

			#Initialize class variables with values from the pdbline
			self.__serial.append(int(float(line[6:11])))
			self.__name.append(line[12:16].strip())
			self.__altLoc.append(altLoc.strip())
			self.__resName.append(line[17:20].strip())
			self.__chainID.append(line[21:22].strip())
			self.__resSeq.append(int(line[22:26]))
			self.__iCode.append(line[26:27].strip())
			self.__x.append(float(line[30:38]))
			self.__y.append(float(line[38:46]))
			self.__z.append(float(line[46:54]))
		
		#get the total number of amino acids
		self.__aa_length = len(set(self.__resSeq))
		#Get the total length of the pdbfile
		self.__iter_length = len(self.__x)
		return
	
	#Return all atoms or iteration length
	def iter_length(self):
		return self.__iter_length
	
	#return amino acid length
	def aa_length(self):
		return self.__aa_length
	
	#return serial if column if the pdb file
	def serial(self):
		return self.__serial
	
	#return the name of the atom
	def name(self):
		return self.__name

	#return altLoc
	def altLoc(self):
		return self.__altLoc
	
	#return the residue name
	def resName(self):
		return self.__resName

	#return the sequence extracted from the structure
	def sequence(self):
		seq = ""
		#iterate over the total iteration length
		for i in range(self.__iter_length):
			#check if the atom is CA, then add one leter code in the seq variable
			if self.__name[i] == 'CA':
				seq = seq + THREEONE[self.__resName[i]]
		return seq

	#return the sequence map that consists of start and end of the iter variable for every residue
	def sequence_map(self):
		seq = ''
		seq_map = []
		
		#Get the iter index for the first residue
		cres=0;
		if self.__name[cres] == 'N':
			seq = seq + THREEONE[self.__resName[cres]]
			seq_map.append([cres]);

		#get the iter index for following residues
		for i in range(1,self.__iter_length):
			if self.__name[i] == 'N':
				seq = seq + THREEONE[self.__resName[i]]
				#append the previous residue iter in the previous residue map
				seq_map[-1].append(i-1);
				#append the current residue
				seq_map.append([i]); 

		seq_map[-1].append(self.__iter_length-1);
			
		return seq, seq_map
	
	#return the chain ID
	def chainID(self):
		return self.__chainID

	#return residue number
	def resSeq(self):
		return self.__resSeq

	#return x coordinate
	def x(self):
		return self.__x
	
	#return y coordinate
	def y(self):
		return self.__y
	
	#return z coordinate
	def z(self):
		return self.__z

	#write pdb file from class object
	def write(self, outfilename):
		fout = open(outfilename,'w')
		for i in range(self.iter_length()):
			serial = str(self.__serial[i])
			name = self.__name[i]
			resName = self.__resName[i]
			chainID = self.__chainID[i]
			resSeq  = str(self.__resSeq[i])
			x = str("%.3f" % self.__x[i])
			y = str("%.3f" % self.__y[i])
			z = str("%.3f" % self.__z[i])
			fout.writelines(pdb_line(serial, name, resName, chainID, resSeq, x, y, z))
		fout.close()
		return

	#write pdb file from class object for selected atoms
	def write_selection(self, outfilename, selection):
		fout = open(outfilename,'w')
		for i in selection:
			serial = str(self.__serial[i])
			name = self.__name[i]
			resName = self.__resName[i]
			chainID = self.__chainID[i]
			resSeq  = str(self.__resSeq[i])
			x = str("%.3f" % self.__x[i])
			y = str("%.3f" % self.__y[i])
			z = str("%.3f" % self.__z[i])
			fout.writelines(pdb_line(serial, name, resName, chainID, resSeq, x, y, z))
		fout.close()
		return

	#write pdb file from class object in multiple models
	def write_selection_multimdl(self, outfilename, selection, modelid):
		fout = open(outfilename,'a')
		fout.writelines('MODEL        '+modelid+'\n')
		for i in selection:
			serial = str(self.__serial[i])
			name = self.__name[i]
			resName = self.__resName[i]
			chainID = self.__chainID[i]
			resSeq  = str(self.__resSeq[i])
			x = str("%.3f" % self.__x[i])
			y = str("%.3f" % self.__y[i])
			z = str("%.3f" % self.__z[i])
			fout.writelines(pdb_line(serial, name, resName, chainID, resSeq, x, y, z))
		fout.writelines('ENDMDL'+'\n')
		fout.close()
		return
	
	#write pdb file from class object for different chains
	def write_all_chains(self, outfilename):
		fout = open(outfilename,'a')
		for i in range(self.iter_length()):
			serial = str(self.__serial[i])
			name = self.__name[i]
			resName = self.__resName[i]
			chainID = self.__chainID[i]
			resSeq  = str(self.__resSeq[i])
			x = str("%.3f" % self.__x[i])
			y = str("%.3f" % self.__y[i])
			z = str("%.3f" % self.__z[i])
			fout.writelines(pdb_line(serial, name, resName, chainID, resSeq, x, y, z))
		fout.close()
		return

#write pdb line in a file according to pdb format
def pdb_line(serial = "", name = "", resName = "", chainID = "", resSeq  = "", x = "", y = "", z = "", occupancy = "", T = "", element = "", charge = "", altLoc = "", iCode = ""):
	'''serial, name, resName, resSeq , x, y, z, occupancy, T, element, charge, altLoc, chainID, iCode'''

	serial = str(serial)
	resSeq = str(resSeq)
	x = str(x)
	y = str(y)
	z = str(z)
	occupancy = str(occupancy)
	T = str(T)
	charge = str(charge)

	atom = 'ATOM  '
	serial = serial.rjust(len(range( 6,11)))
	name = name.ljust(len(range(12,15)))
	altLoc = altLoc.ljust(len(range(16,17)))
	resName = resName.ljust(len(range(17,20)))
	chainID = chainID.ljust(len(range(21,22)))
	resSeq = resSeq .rjust(len(range(22,26)))
	iCode = iCode.ljust(len(range(26,27)))
	x = x.rjust(len(range(30,38)))
	y = y.rjust(len(range(38,46)))
	z = z.rjust(len(range(46,54)))
	occupancy = occupancy.rjust(len(range(54,60)))
	T = T.rjust(len(range(60,66)))
	element = element.ljust(len(range(76,78)))
	charge = charge.rjust(len(range(78,80)))

	line = atom + serial + '  ' + name + altLoc + resName + ' ' + chainID + resSeq + iCode + ' '*3 + x + y + z + occupancy + T + ' '*11 + element + charge
	line = line[:80].strip() + '\n'

	return line
