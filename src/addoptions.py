"""
Copyright (C) 2019 Neelesh Soni <neelesh.soni@alumni.iiserpune.ac.in>, <neeleshsoni03@gmail.com>
This file is part of AIMorpher.

AIMorpher is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

AIMorpher is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with AIMorpher.  If not, see <http://www.gnu.org/licenses/>.
"""

'''
This function add the options to the interpreter from the command line

'''
def addoptions():
	
	import optparse
	
	parser = optparse.OptionParser()

	exm_i = "./example/F_Protein_PRE_STATE_TRIMER.pdb"
	exm_f = "./example/F_Protein_POST_STATE_TRIMER.pdb"
	exm_o = "./output/Morph_PRE_POST_Transition.pdb"
	exm_n = 5
	
	#Adding the options
	parser.add_option('-i', default=exm_i, type="string", help='Initial state input pdb file. For Example: -i ./input/state_initial.pdb' )
	parser.add_option('-f', default=exm_f, type="string", help='Final state input pdb file. For Example: -f ./input/state_final.pdb' )
	parser.add_option('-o', default=exm_o, type="string", help='Morphed output pdb file name For example: -o ./output/Morph_Initial_Final.pdb' )
	parser.add_option('-n', default=exm_n, type="int", help='Total number of intermediate structures to generate');
	
	args, remainder = parser.parse_args()
	
	#Check if number of arguments are correct or not.
	if (args.f is None) or (args.o is None):
		parser.error("Not enough number of arguments\nUse -h for help")
	
	#Assigning the command line options to variables
	fname1 = args.i;
	fname2 = args.f;
	oname= args.o;
	num_intm= int(args.n)
		
	#return the input and output files
	return fname1,fname2,oname,num_intm

