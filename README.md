# AIMorpher
Rigid body transformation and Morphing tool


# Help on AIMorpher
python ./src/main.py -h

Usage: main.py [options]

Options:
  -h, --help  show this help message and exit
  -i I        Initial state input pdb file. For Example: -i
              ./input/state_initial.pdb
  -f F        Final state input pdb file. For Example: -f
              ./input/state_final.pdb
  -o O        Morphed output pdb file name For example: -o
              ./output/Morph_Initial_Final.pdb
  -n N        Total number of intermediate structures to generate


# Test Usage with Default Parameters and files provided in example folder
python ./src/main.py

Above will run with the following default parameters.
Input files: 'F_Protein_PRE_STATE_TRIMER.pdb' and 'F_Protein_POST_STATE_TRIMER.pdb' 

Output is generated in ./output folder. 
./output/Morph_PRE_POST_Transition.pdb
./output/Morph_PRE_POST_Transition.pdb_output_alignment.txt

Number of Intermediate States.
5