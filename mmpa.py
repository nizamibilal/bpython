#!/usr/bin/python
#filename: mmpa.py

#####################################################################
#		Matched molecular Pair Analysis
#	    This code requires working biopython installed			
#			written by Bilal Nizami			
#
#			    UKZN, Durban
#				2016
#
#
#####################################################################
#
#=======================================================================
# import functions

from __future__ import print_function
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
#========================================================================

## Welcome message
print ('\n', '*'*60, '*\n\n\t\tMatched Molecular Pair Analysis\t*\n', '*\t\tBy: Bilal Nizami*\n', '*\t\tnizamibilal1064@gmail.com*\n','*\t\tUKZN, Durban 2016*\n\n\n','*'*60, '\n', end="")

#===================
# output file
outputfile = 'out.txt'
of = open(outputfile, 'w')
print('SMILE', ',','Atom,', 'Index,','Neighbor,', 'Bond type,', 'if in ring,', file=of)

fp = open('single_smiles.txt', 'r')
smfile = fp.read()
smfile = smfile.split('\n')
fp.close
#print (smfile)
for f in smfile:
	#print (f)
	mols = Chem.MolFromSmiles(f)
	atom = mols.GetAtomWithIdx(0)
	atoms = mols.GetAtoms()
	#print ('\n', f)
	for at in atoms:
		atsym = at.GetSymbol()
		print (atsym)
		#print (at.GetBonds()[0].GetBondType())
		#print ('neighbour', at.GetNeighbors())
		##atom = m.GetAtomWithIdx(0)
		#print ([x.GetSymbol() for x in at.GetNeighbors()])
		for x in at.GetNeighbors():
			neighatm = x.GetSymbol()
			print ('neighbor atom', neighatm)
			
			## return the bond object between two atom index, first atom is neighbor of at in atoms
			bondbtw = mols.GetBondBetweenAtoms(x.GetIdx(), at.GetIdx()) 
			bondtp = bondbtw.GetBondType()
			print (bondtp)
			
			## check if atom is in ring formation with its neighboring atoms
			if x.IsInRing() and at.IsInRing():
				print (x.GetSymbol(), 'and', at.GetSymbol(), 'are in ring')
				rstat = 'in ring'
			else:
				rstat = 'not in ring'
			print (f, ',',atsym, ',', at.GetIdx(),',',neighatm, ',',bondtp, ',',rstat, file=of)
		#print (len(at.GetNeighbors()[0].GetBonds()))
		print ('\n')

		#if mols.GetAtomWithIdx(0).IsInRing():
		#	print ('in ring')
		#print (mols.GetBonds()[0].GetEndAtom())
#mols = Chem.MolFromMolFile('test_smiles.mol')
#atoms = mols.getAtoms()
#mols = [x for x in suppl]
#index = 0
#print ('molecule ',mols)
#for mol in suppl:
#for a in atoms:
#	atnum = a.getAtomiNum()
#	print ('atomic num', atnum)
	#print(m.GetNumAtoms())
	#atoms = m.GetAtoms()
	#atoms = atoms.GetSymbol()
	#atom = m.GetAtomWithIdx(index).GetSymbol()
	#print ('', atom, atoms)
	#index = index+1
of.close()
	
