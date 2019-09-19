#!/usr/bin/env python

"""
Ce script cherche a identifier les genes a proximite d'un gene de resistance passe
en argument.

Usage : ./context.py amr_gene gbk_folder_path number_of_genes_in_context

"""

__auteur__ = "Jean-Simon Brouard"
__date__ = "2019-05-08"

# Importation de modules standards
import sys
import os
from Bio.Graphics import GenomeDiagram
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from io import StringIO
from collections import OrderedDict
import pandas as pd
import numpy
import re
import csv
import glob

# Importation des modules locaux
from foreign_code import prepare_color_dict
from foreign_code import digest_features
from foreign_code import get_coord_gene

from geneCluster import geneCluster


# try to use gene cluster
a1 = geneCluster("small")
b1 = geneCluster("medium")
b2 = geneCluster("mediumCopy")
c1 = geneCluster("large")
d1 = geneCluster("test")

r1 = geneCluster("envers")


		

# Disons que l'on veut faire grossir nos regroupements de genes

a1.add("rci", "+")
a1.add("CTX-M-1", "+")



# 2 identiques
b1.add("rci", "+")
b1.add("CTX-M-1", "+")
b1.add("int1", "+")

b2.add("rci", "+")
b2.add("CTX-M-1", "+")
b2.add("int1", "+")



c1.add("pilA", "-")
c1.add("int1", "-")

d1.add("rci", "+")
d1.add("CTX-M-1", "+")
d1.add("int1", "+")
d1.add("pilA", "+")
d1.add("pilV", "+")
d1.add("pilVI", "+")
d1.add("pilU", "+")




# try to add a small contig, but on the reverse strand

r1.add("CTX-M-1", "-")
r1.add("rci", "-")


def cluster_Gene_clusters(gene_cluster):
	"""This function will try to classify a cluster into an appropriate entry in a dictionnary d
	Pour que cette fn soit efficace, il faut 'envoyer' les gene clusters a regrouper en ordre
	decroissant de taille, de maniere a ce que les plus gros regroupements se forment, puis que
	ceux qui sont identiques ou plus petits mais compatibles s'y ajoutent.
	"""

	s = ""
	s = gene_cluster.read()

	# check if we have an entry for that particular gene cluster
	if s in d:
		d[s]['n'] += 1
		d[s]['geneClusterObjects'].append(gene_cluster)
		print ("--> Incrementing an existing gene cluster")
		return

	else:

		# Try if the reverse complement of the cluster match an entry in the dict
		reverse_gene_cluster = gene_cluster.reverse()
		r = ""
		r = reverse_gene_cluster.read()
		if r in d:
			print ('--> Adding a reverse gene cluster to an existing category yeah!')

			# We count this cluster and add relative info in the entry of the 'unreversed' cluster
			d[r]['n'] += 1
			d[r]['geneClusterObjects'].append(gene_cluster)
			return

		# check if the cluster is included in another cluster(larger)
		else:

			# Ce bloc de code tente d'inclure des regroupements de genes
			# dans des regroupement compatibles plus grands 
			# boucler sur tous les clusters presents dans le dict d
			for key, value in d.items():
				if isinstance(value, dict):
					for sub_key, sub_value in value.items():
						if sub_key == 'geneClusterObjects':
							for clus in sub_value:
								if gene_cluster.isin(clus):
									print (" --> Including a gene cluster into an existing cluster of gene clusters")
									# la methode read donne les genes du gene cluster, c'est aussi la clef de notre
									# dict de regroupements de gene clusters
									d[clus.read()]['n'] += 1
									d[clus.read()]['geneClusterObjects'].append(gene_cluster)

									print('parent -->' + gene_cluster.parent + '<--parent')

									return
	
								# Tentative en retournant le custer
								else:
									if reverse_gene_cluster.isin(clus):
										print (" --> Including a reversed gene cluster into an existing cluster of gene clusters")
										d[clus.read()]['n'] += 1
										d[clus.read()]['geneClusterObjects'].append(gene_cluster)
										return
			
			print ('--> Adding a new gene cluster')
			d[s] = {}
			d[s]['n'] = 1
			d[s]['geneClusterObjects'] = [gene_cluster]
			return







# Programme principal

if __name__ == "__main__":

	d = {}


	print("test new function")





	print(c1.isin(d1))
	print(c1.read())
	print(d1.read())


	L=[a1,b1,b2,c1, r1]

	L.append(d1)

	# This command allow to sort a list of objects according to an attribute
	newlist = sorted(L, key=lambda x: x.nb_genes, reverse=True)

	# Verification que mes objects sont bien en ordre descendants de taille
	for i in newlist:
		print('dealing with this cluster:(', i.cluster, ') ', 'containing ',  i.nb_genes, ' genes and named: ', i.name)

		cluster_Gene_clusters(i)	


	# Passe au travers du dictio pour 	
	for key, value in d.items():
		if isinstance(value, dict):
			for sub_key, sub_value in value.items():

				if sub_key == 'n':
					print('regroupement de taille', sub_value)

				if sub_key == 'geneClusterObjects':
					for clus in sub_value:

						print(clus.read(), clus.name, clus.nb_genes)













