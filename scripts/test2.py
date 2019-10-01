#!/usr/bin/env python

"""
Ce script cherche a tester les nouvelles fn de la classe GeneCluster

Usage : ./test2.py

"""

__auteur__ = "Jean-Simon Brouard"
__date__ = "2019-09-26"

# Importation de modules standards
import sys
import os
from pathlib import Path
#from Bio.Graphics import GenomeDiagram
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from io import StringIO
from collections import OrderedDict
#import pandas as pd
import numpy
import re
import csv
import glob

# Importation des modules locaux
from foreign_code import prepare_color_dict
from foreign_code import digest_features
from foreign_code import get_coord_gene
from geneCluster import geneCluster
from dock import dock

# Programme principal

if __name__ == "__main__":




	# Initialisons quelques gene clusters

	

	# a0
	a0 = geneCluster("a0")
	geneList = ['A', 'B', 'C']
	strandList = ['+','+','+']
	a0.load(geneList, strandList)
	print(a0.name + '\t' + a0.read())


	# a1
	a1 = geneCluster("a1")
	geneList = ['A', 'B', 'C']
	strandList = ['+','+','+']
	a1.load(geneList, strandList)
	print(a1.name + '\t' + a1.read())

	# r1
	r1 = geneCluster("r1")
	geneList = ['C', 'B', 'A']
	strandList = ['-','-','-']
	r1.load(geneList, strandList)
	print(a1.name + '\t' + a1.read())





	# b1
	b1 = geneCluster("b1")
	geneList = ['A', 'B', 'C', 'D']
	strandList = ['+','+','+','+']
	b1.load(geneList, strandList)

	# c1
	c1 = geneCluster("c1")
	geneList = ['A', 'B', 'C', 'D']
	strandList = ['+','+','+','+']
	c1.load(geneList, strandList)
	c1.add('E', '+')
	print(c1.name + '\t' + c1.read() + ' after the addition')
	

	# d1
	d1 = geneCluster("d1")
	geneList = ['A', 'C', 'D']
	strandList = ['+','+','+']
	d1.load(geneList, strandList)
	print(d1.name + '\t' + d1.read())


	print ('--Testing operator overloading--')

	if a1 == b1:
		print('c1 equal a1')
	else:
		print('c1 not equal a1')

	if a0 != 1:
		print('a0 not equal a1')
	else:
		print('a0 equal a1')



	print ('--Testing isin fn--')


	# a1 vs r1
	if a1.isin(r1):
		print ('yes a1 is in r1')
	else:
		print ('nope a1 is not in r1')
	

	# d1 vs a1
	if d1.isin(a1):
		print ('yes d1 is in a1')
	else:
		print ('nope d1 is not in a1')



	# Tentons de travailler avec notre nouvelle classe : dock

	my_dock = dock()


	my_dock.add(a0)
	my_dock.add(a1)
	my_dock.add(c1)
	my_dock.add(r1)

	print('dock size\t' + str(my_dock.size))
	print(my_dock.elems[3].cluster)		# marche!


	# pour avoir le header
	print(my_dock.head.cluster)		# marche!
















