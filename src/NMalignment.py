"""
Copyright (C) 2019 Neelesh Soni <neelesh.soni@alumni.iiserpune.ac.in>, <neeleshsoni03@gmail.com>
This file is part of AIMorpher.

AIMorpher is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

AIMorpher is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with AIMorpher.  If not, see <http://www.gnu.org/licenses/>.
"""


'''
This file codes for global alignment of two sequences. 
1. The code first reads a substitution matrix file for scoring a substitution.
2. The matrix could be provided by the user or as a default option of the software.
3. The alignment scoring matrix was filled for global alignment using the algoirithm.
4. The tracback of the best scoring aligned sequences will be doine using recursive traceback function.
5. Finally, the aligned sequences were returned by the module.
'''

from sys import setrecursionlimit
'''
This module read the substitution matrix file.
'''
def read_matrix_file(fname):

	#Open the file in reading mode
	inf = open(fname,'r')
	lines = inf.readlines();
	inf.close();

	#temperory list for storing
	tempmat=[];

	#iterate over all lines
	for l in lines:
		#skip if first character in the file is '#'
		if l[0]=='#':
			continue;

		#first remove the trailing '\n' and then split the lines 
		#for removing multiple spaces as delimeter
		toks = l.strip().split();

		#append the lins converted to a list
		tempmat.append(toks)

	#Score matrix variable for storing substitution scores
	scoremat={}
	#Get all rows except first that contains matrix headers
	submat = tempmat[1:]

	#Iterate over all rows
	for i in range(0,len(submat)):
		#Iterate over columns
		for j in range(0,len(tempmat[0])):

			#get column res
			colid = tempmat[0][j]
			#get row res
			rowid = submat[i][0]
			#get the value of res-res pair from substitution matrix
			value = submat[i][j+1]
			
			#Append the score matrix variable with key-value pairs
			scoremat[(rowid,colid)]=float(value);

	#return the scoring matrix variable dictionary
	return scoremat

'''
This module traceback the best scoring sequence using recursion on the alignment scoring matrix.
'''
def recursive_traceback(i,j,aligned_score,score_matrix,seq1,seq2,ali1,ali2):

	#Stop the recursion of the index becomes 0 or less. 
	#This is to ignore the first row and first column.
	#Return the alignment sequences as it is.
	if ((i <= 0) or (j <= 0)):
		return ali1,ali2

	#Get the scores of current cell of the aligned scoring matrix.
	#Get corresponding diagonal, left and top scores
	score_current = aligned_score[i][j]
	score_diagonal = aligned_score[i-1][j-1]
	score_left = aligned_score[i][j-1]
	score_up = aligned_score[i-1][j]

	#Check if current score is obatined from the diagonal element
	if  score_current == score_diagonal + score_matrix[(seq1[j-1], seq2[i-1])]:
		#If true then append the updated alignment
		ali1 += seq1[j-1]; ali2 += seq2[i-1]
		#procede recursion towards diagonal for next iteration
		ali1,ali2=recursive_traceback(i-1,j-1,aligned_score,score_matrix,seq1,seq2,ali1,ali2)
		return ali1,ali2
	
	#Check if current score is obatined from the left element
	elif  score_current == score_left + score_matrix[(seq1[j-1], '*')]:
		#If true then append the updated alignment
		ali1 += seq1[j-1]; ali2 += '-'
		#procede recursion towards left for next iteration
		ali1,ali2 = recursive_traceback(i,j-1,aligned_score,score_matrix,seq1,seq2,ali1,ali2)
		return ali1,ali2
	
	#Check if current score is obatined from the top element
	elif  score_current == score_up + score_matrix[('*', seq2[i-1])]:
		#If true then append the updated alignment
		ali1 += '-'; ali2 += seq2[i-1]
		#procede recursion towards top for next iteration
		ali1,ali2 = recursive_traceback(i-1,j,aligned_score,score_matrix,seq1,seq2,ali1,ali2)
		return ali1,ali2

	return ali1,ali2

'''
This is the main function of the module
This module provide module for global alignment of two sequences.
'''
def global_align(seq1, seq2, matrix):

	#get the length of the sequences in two variables, n and m
	n = len(seq1)
	m = len(seq2)
	
	#set the recursion limit according to the longest sequence + offset (normally 100)
	setrecursionlimit(max(n,m)+100)

	#Get the substitution matrix
	score_matrix = read_matrix_file(matrix)

	#Initialize the alignment scoring matrix with all zeros. 
	#Its a rowsXcolumn => (m+1)X(n+1) matrix
	aligned_score = [[0 for x in range(0,n+1)] for y in range(0,m+1)] 

	#Iterate over first column
	for j in range(0, n + 1):
		#substitution matrix has '*' for gaps
		gap_score = score_matrix[(seq1[j-1], '*')]
		#Linearly increment the gap score
		aligned_score[0][j] = gap_score * j

	#Iterate over first row
	for i in range(0, m + 1):
		#substitution matrix has '*' for gaps
		gap_score = score_matrix[('*', seq2[i-1])]
		#Linearly increment the gap score
		aligned_score[i][0] = gap_score * i
	
	#Iterate overs all rows
	for i in range(1, m + 1):
		#Iterate overs all columns
		for j in range(1, n + 1):
			
			# Calculating the score of top (seq2 aligned with gap), left(seq1 aligned with gap), and diagonal (match) cells
			# match, delete and insert is wrt seq 2
			
			#seq1 aligned with seq2
			match_score = aligned_score[i - 1][j - 1] + score_matrix[(seq1[j-1], seq2[i-1])]
			
			#seq2 aligned with gap
			delete_score = aligned_score[i - 1][j] + score_matrix[('*', seq2[i-1])]
			
			#seq1 aligned with gap
			insert_score = aligned_score[i][j - 1] + score_matrix[(seq1[j-1], '*')]
			
			#Get the maximum score among diagonal, left and top
			aligned_score[i][j] = max(match_score, delete_score, insert_score)

	#Initialize aligned sequences for seq1 and seq2
	ali1 = ""
	ali2 = ""

	#Call the recursive functions to traceback the aligned scores
	ali1,ali2 = recursive_traceback(m,n,aligned_score,score_matrix,seq1,seq2,ali1,ali2)
	
	#De-reversing the sequence due to traceback
	ali1 = ali1[::-1]
	ali2 = ali2[::-1]
	
	return(ali1, ali2)