#!/usr/bin/python
#filename: pdbminer.py

#####################################################################
#		Search the protein sequence fasta file
#	    This code requires working biopython installed			
#			written by Bilal Nizami			
#
#			    UKZN, Durban
#				2016
#
#
#####################################################################
#
# import functions
from __future__ import print_function
from optparse import OptionParser
from Bio import SeqIO
import sys
import re
from Bio.PDB import *

## Welcome message
print ('\nPDBMINER\n', 'written by Bilal Nizami\n', 'UKZN, Durban 2016\n\n\n', end="")

#==============================================================================
# Setting the options
#==============================================================================
parser = OptionParser("Usage: seqsearch.py -i --fasta fasta file name ")

parser.add_option("-i", "--fasta", type='string', dest="fasta",
                  help="A fasta file having amino acid sequence")

parser.add_option("-o", "--out", type='string', dest="output",
                  help="output text file")
                
parser.add_option("-p", "--path", type='string', dest="path",
                  help="string for path to the pdb files, it will be used for accesing the pdb file of parsed \
                  pdbid")                                                  

(options, args) = parser.parse_args()
                     

#====================================================================
# if no arguments are passed
if options.fasta is None: # where fasta is obviously your required option
    parser.print_help()
    sys.exit(1)
    
#===========================================================================
#
#	function to get the pdbid from fasta file
#===========================================================================
def fasta_get_pdbid(filename):
	"fasta_get_pdbid() parse the pdbid from single fasta file having multiple sequences"
	fastaf = []
	pdbid = []
	fp = open(filename, 'r')
	fastaf = fp.read()
	fastaf = fastaf.split('>')
	fastaf = fastaf[1:len(fastaf)] ## knock off 0th index from list, this is done becuase preveious line of split insert a null at index 0
	for line in fastaf:
		line = line[0:4]
		pdbid.append(line)
	fp.close() 
	#print (pdbid)
	return pdbid;
#===========================================================================

# read the single fasta file containing the multiple sequences. #

pdb = []
fastafile = ''
seq = ''
aacid = []

A = R = N = D = B = C = E = Q = Z = G = H = I = L = K = M = F = P = S = T = W = Y = V = 0
#aminoacid = ['A', 'R', 'N', 'D', 'B', 'C', 'E', 'Q', 'Z', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']


pdb_id = []
pdb_id_chain = []
pdbid = []
pdbchain = []
count = []
counter=0
sf = open(options.fasta, "rU")
pdbid = fasta_get_pdbid(options.fasta)
#print (pdbid)
#for pdb in pdbid:
	
for p in SeqIO.parse(sf, "fasta"): ## takes a file handle and format name, and returns a SeqRecord iterator using biopython SeqIO
	#print (p.id)
	#pdb_id_chain = re.findall("\d.{3}:[A-Z]", p.id)
	#print (str(pdb_id_chain))
	#pdb_id = re.findall("\d.{3}", str(pdb_id_chain))
	#print (re.findall("\d.{3}", pdb_id[0]), end='')
	#print ('PDBID"\t', pdb_id[0], end='')
	#print ('PDBID and chain:\t', pdb_id_chain)
	
	#pdbchain.append(pdb_id_chain)
	#pdbid.append(pdb_id[0])
	aacid = list(p.seq)
 	for aa in aacid:
		if aa == 'A':
			A=A+1
		elif aa == 'R':
			R=R+1
		elif aa == 'N':
			N=N+1
		elif aa == 'D':
			D=D+1
		elif aa == 'B':
			B=B+1
		elif aa == 'C':
			C=C+1
		elif aa == 'E':
			E=E+1
		elif aa == 'Q':
			Q=Q+1
		elif aa == 'Z':
			Z=Z+1
		elif aa == 'G':
			G=G+1
		elif aa == 'H':
			H=H+1
		elif aa == 'I':
			I=I+1
		elif aa == 'L':
			L=L+1
		elif aa == 'K':
			K=K+1
		elif aa == 'M':
			M=M+1
		elif aa == 'F':
			F=F+1
		elif aa == 'P':
			P=P+1
		elif aa == 'S':
			S=S+1
		elif aa == 'T':
			T=T+1
		elif aa == 'W':
			W=W+1
		elif aa == 'V':
			V=V+1
		elif aa == 'Y':
			Y=Y+1
	counter=counter+1
	#print ('. ', end="")
print ('total', counter, 'processed') 
#print ('PDBID "\t', pdbid, pdbchain, end='')	
sf.close()
#=============================================================================
# open Summary file and print into it

sf = open(options.output, 'w')
print('Amino Acid,', 'Count,', file=sf)


print('A\t', A, '\nR\t', R, '\nN\t', N, '\nD\t', D, '\nB\t', B, '\nC\t', C, '\nE\t', E, '\nQ\t', Q, '\nZ\t', Z, '\nG\t', G, '\nH\t', H, '\nI\t', I, '\nL\t', L, '\nK\t', K, '\nM\t', M, '\nF\t', F, '\nP\t', P, '\nS\t', S, '\nT\t', T, '\nW\t', W, '\nY\t', Y, '\nV\t', V, file=sf)

sf.close()

#===========================================================================================================
#
# sort and store only unique pdb ids, there might be redundent parsed pdbids due to seperate entry in input 
# fasta file for multiple chains of a pdbid. it uses built in set() method
#

pdbid = list(set(pdbid)) # remove duplicates pdbids 
#print (pdbid)

## opening file for dumping resolutions
sf = open('res.txt', 'w')
if options.path != None:		
	print ('Resolution', file=sf)
	resolution = ''
	pdb_path = []
	parser = PDBParser()
	for pdb in pdbid:
		try:
			pdb_dir = pdb[1:3]
			pdb_dir = pdb_dir.lower()
			pdb = pdb.lower() ## lower case, becuase last character fo pdbid from fasta file is in uppercase
			#print (pdb)
			pdb_path1 = str(options.path) + pdb_dir + '/pdb' + str(pdb) + '.ent' ## full path to the pdb.ent file
			pdb_path.append(pdb_path1)
			structure = parser.get_structure('PHA-L', pdb_path1)
			resolution = str(structure.header['resolution'])
			print ('for ', pdb)
			print (resolution, file=sf)
			#print (pdb_path)
		except IOError:
			print ('Can not find', pdb_path1, ' file')	
	#for pdb in pdb_path:
	#	try:
	#		#print (pdb)
	#		structure = parser.get_structure('PHA-L', pdb)	
	#		resolution = str(structure.header['resolution'])
	#		print ('for ', pdb)
	#		print (resolution, file=sf)
	#		#print (resolution)
	#	except IOError:
	#		print ('Can not find', pdb ' file')		
			

sf.close()
#	open file with the path to pdbs
#if options.pdblist != None:
#	pdbfl=[]
#	fp = open(options.pdblist, 'r')
#	pdbfl =	fp.read()
#	pdbfl = pdbfl.split('\n')
#	fp.close()

### iterate over pdbfl
#test = ''
#test1 = ''
#j = 0
#for i in pdbfl:
#	test1 = i
#	test = str(pdbid[j])
	#print (test, test1)
#	j=j+1
#	if test.find(test1):
#	 print ('found', test)
