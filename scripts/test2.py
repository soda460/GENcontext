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


# Programme principal

if __name__ == "__main__":




	# Initialisons quelques gene clusters

	

	b1 = geneCluster("medium")
	b2 = geneCluster("mediumCopy")

	c1 = geneCluster("large")
	d1 = geneCluster("test")

	r1 = geneCluster("envers")

	# a1
	a1 = geneCluster("small")
	geneList = ['A', 'B', 'C', 'D', 'E', 'F']
	strandList = ['+','+','+','+','+','+']
	a1.load(geneList, strandList)

	# b1
	b1 = geneCluster("small like a1")
	geneList = ['A', 'B', 'C1', 'D', 'E', 'F']
	strandList = ['+','+','+','+','+','+']
	b1.load(geneList, strandList)


	if b1 != a1:
		print('not equal')



	print(a1.read())


























