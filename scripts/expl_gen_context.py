#!/usr/bin/env python

"""
Ce script cherche a identifier les genes a proximite d'un gene de resistance passe
en argument.

Usage : ./context.py (i)amr_gene (ii)genbank_folder_path (iii)number_of_genes_in_context

"""

__auteur__ = "Jean-Simon Brouard"
__date__ = "2019-09-17"

# Importation de modules standards
import sys
import os
from pathlib import Path
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



def get_gbk_files(start_dir):
	""" Given a directory, make a list of genbank files
	    This function will find all files recursively.
	"""
	gbk_files = []
	for filename in Path(start_dir).glob('**/*.gbk'):
		print(filename)
		gbk_files.append(filename)
	return gbk_files








def examine_gbk_file (gene, gbk_f):
	""" Parse the content of a genbank file with BioPython SeqFeature module
		Decouple the parsing with the SeqFeature engine and the function
		that will check gene presence for a particular gene.
	""" 

	hits = []	# a list that may contain features and feature_iterators for that particular file

	# Get the gbk records with BioPython | many files have more than 1 record (ended by //) !
	for i_gbk, gbk_rec in enumerate(SeqIO.parse(gbk_f, "gb"), start=1):

		# If required print stuff to standard output
		# s = "Handling GBK record no. " + repr(i_gbk) + ":" + gbk_rec.name	# car i_gbk est un int
		# print (s)

		# Check in each record for the gene of interest
		abc = check_gene_presence(gene, gbk_rec)

		hits.extend(abc)

	return (hits)


def check_gene_presence(gene, record):
	""" Check the presence of the AMR gene given in argument of the script 
		This should be a gene referenced in th CARD database and annotated in an inference fied by PROKKA
		This function return a 4-items list : an index, the feature, the AMR gene that have been found, and the entire record
	"""

	record_out= []

	for i_feature, feature in enumerate(record.features):		
		if feature.type == 'CDS':
			if 'inference' in feature.qualifiers:
				for inf in feature.qualifiers['inference']:
					if ("protein_fasta_protein_homolog_model" or 
						"protein_fasta_protein_knockout_model" or
						"protein_fasta_protein_overexpression_model" or
						"protein_fasta_protein_variant_model") in inf:
						inf2 = inf.replace(" ", '') 	# To remove extra space converted from \n in genbank
						b = inf2.split("|")				# To get the gene name
						my_current_gene = (b[3])
						
						# Fill an array with feature and iterator if current gene is the one that we are looking for 
						if gene in my_current_gene:
							record_out.append([i_feature, feature, my_current_gene, record])
	return (record_out)
					
						
	
def get_metadata(f, out_list):
	""" Get and print important information about the strain, molecule, WGS contig, locus where
		the amr gene was found 
	"""

	# This is related to the organization of files
	# For now, this script will work with the following files organization:
	# strains folders then molecule folder then gbk files /Res13-blabla-strain1/plasmid_476/Res13-blabla_plasmid_476.gbk
	gbk_file_name = (f.split('/')[-1])
	molecule_name = (f.split('/')[-2])
	strain_name = (f.split('/')[-3])
	meta = strain_name + "\t" + molecule_name + "\t" + gbk_file_name

	# First of all, try to grab the amr gene name 
	# In item we should have [index, feature, AMR_gene, gbk_record_name (== wgs contig)]
	for item in out_list:

		feature_position_in_record = repr(item[0])	# un index

		# Here we have a seqFeature object for our desired AMR gene
		amr_found_feature = item[1]
		if 'locus_tag' in amr_found_feature.qualifiers:
			locus = amr_found_feature.qualifiers['locus_tag'][0]

		amr_found = item[2]
		gbk_record_name = item[3].name


		# Now we are ready to not print this, the infos is returned anyway
		#print (amr_found + "\t"
		#		+ meta + "\t"
		#		+ gbk_record_name + "\t"
		#		+ locus + "\t"
		#		+ feature_position_in_record + "\n")

	# Return a list, not a tuple
	return [amr_found, gbk_file_name, molecule_name, strain_name, gbk_record_name, locus]


def splice_record(index, record, features_range):
	"""This function will try to get features around a given feature"""


	# to print the entire record
	# print(record.features[index])

	nb_features = len(record.features)		# number of features in the record

	# Math pour savoir la plage de feature qui sera consideree
	upper_bound = index + features_range

	# if the index is too close from the end, the upper_bound will be nb_features minus 1 (0-indexed)
	upper_bound = min(upper_bound, (nb_features - 1))	

	lower_bound = index - features_range

	# if the index is too close from the start, the lower_bound will be 0
	if (lower_bound <= 0):
		lower_bound = 1		# !!! Because the first feature is usually source and it describes the whole molecule


	# s = "Feature range will be " + repr(lower_bound) + ":" + repr(upper_bound)
	# print(s)

	
	'''
	Voici les qualifiers disponibles
	locus_tag
	codon_start
	inference
	translation
	transl_table

	location n'en fait pas partie

	Pour acceder au coordonnes en nucletotide du
	print (record.features[index].location.start)
	
	https://biopython.org/DIST/docs/api/Bio.SeqFeature.FeatureLocation-class.html
	En somme pour acceder aux elements on a les methodes start, end, et strand
	'''

	'''Note that the start and end location numbering follow Python's scheme, thus a GenBank entry of 123..150 (one based counting)
	becomes a location of [122:150] (zero based counting).'''

	lower_pos = record.features[lower_bound].location.start
	upper_pos = record.features[upper_bound].location.end

	# s = "Original record will be spliced at positions " + repr(lower_pos) + ":" + repr(upper_pos)
	# print(s)

	sub_record = record[lower_pos:upper_pos]

	return (sub_record)


def extract_gene_cluster (record):
	"""This function will give gene cluster and orientations"""

	# Use our own class geneCluster

	a = geneCluster(record.name)	

	for feature in record.features:
		if feature.type == 'CDS':		# We use the CDS fields
			flagGeneDraw = 'FALSE'		# Initialize this 'Boolean'!

			if 'inference' in feature.qualifiers:					
				for inf in feature.qualifiers['inference']:		# In some cases, more than one inference tag are present in gbk file.
					if "protein_fasta_protein_homolog_model" in inf:
						inf2 = inf.replace(" ", '') 	# To remove extra space converted from \n in genbank
						b = inf2.split("|")				# To get the gene name, 
						my_gene = (b[3])
						a.add(my_gene, feature.strand)
						flagGeneDraw = 'TRUE'
						break
			
					if "ISfinder" in inf:		
						inf2 = inf.replace(" ", '') 	# To remove extra space converted from \n in genbank
						b = inf2.split(":")				# To get the gene name, 
						my_gene = (b[2])
						a.add(my_gene, feature.strand)
						flagGeneDraw = 'TRUE'
						break

			# Continue searching...
			if flagGeneDraw == 'FALSE':
				if 'locus_tag' in feature.qualifiers:
					if 'gene' in feature.qualifiers:
						my_gene = feature.qualifiers['gene'][0]
						a.add(my_gene, feature.strand)
					else:
						my_gene= '?'
						a.add(my_gene, feature.strand)
	return (a)





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
			d[s]['amr_found'] = gene_cluster.amr_found
			d[s]['n'] = 1
			d[s]['geneClusterObjects'] = [gene_cluster]
			return













# Programme principal

if __name__ == "__main__":

	""" The program will find gbk files and will produce one features list per gbk file""" 

	# Dict of gene clusters
	d = {}
	
	# List of geneCluster objects
	L = []

	# Parsing arguments
	amr_in = sys.argv[1]	# There is a subtle diff between amr_in (ex. 'TEM' and amr_found (TEM-1 or TEM-212)
	folder = sys.argv[2]
	nb_genes_in_context = int(sys.argv[3])

	gbk_files = get_txt_files(folder)		

	for f in gbk_files:
		hits = examine_gbk_file(amr_in, f)

		'''
			Note that the hits list is a list of lists, which contains distinct objects
			[0], is the index of the amr_gene, i.e. the feature position of the amr_gene in the record
			[1], is the feature of the amr_gene 
			[2], is the amr_gene (string)
			[3], is the complete record
 
		'''

		if hits:
			for my_hit in hits:

				# Get the metadata
				my_meta_stuff = []
				my_meta_stuff = get_metadata(f, hits)

				# Give intelligent variables names to metadata
				amr_found = my_meta_stuff.pop(0)
				gbk_file_name = my_meta_stuff.pop(0)
				molecule_name = my_meta_stuff.pop(0)
				strain_name = my_meta_stuff.pop(0)
				gbk_record_name = my_meta_stuff.pop(0)
				locus = my_meta_stuff.pop(0)

				
	
				# The index of the amr_gene, i.e. the feature position of the amr_gene in the record
				amr_index = my_hit[0]
				
				# Correspond to the complete record
				amr_record = my_hit[3]

				# Splice the record where an AMR gene has been found
				sub_rec = splice_record(amr_index, amr_record, nb_genes_in_context)

				# Obtain a GeneCluster object from a subrecord 
				my_gene_cluster = extract_gene_cluster(sub_rec)


				# Add metadata to our gene cluster object
				my_gene_cluster.locus=locus
				my_gene_cluster.amr_found=amr_found
				my_gene_cluster.molecule_name=molecule_name
				my_gene_cluster.strain_name=strain_name


				# Simply put the GeneCluster object in a list
				L.append(my_gene_cluster)


	# Sort geneCluster objects according to their size
	newlist = sorted(L, key=lambda x: x.nb_genes, reverse=True)
		
				

	for c in newlist:

		# Verification que mes objects sont bien en ordre descendants de taille
		# print(c.nb_genes, c.name, c.cluster)

		# This will modify d a dictionnary
		# Here we cluster similar or identical gene clusters
		cluster_Gene_clusters(c)
				


				


	# we have a cool dictionnay d and we want to sorted it by the reverse order of the frequency of clusters
	# We can do the following, but this will give just a list with our desired sorted keys
	# Used at line 505
	sorted_keys = sorted(d, key=lambda x: (d[x]['n']), reverse=True)


	# After some stackoverflow searches, I,ve found the magic command that allow to create an OrderedDict object
	# that will remember the insertion order of the items. We want to create the keys in reverse order of ['n']
	# which is the number of occurence of the clusters

	d2 = OrderedDict(sorted(d.items(), key=lambda x: x[1]['n'], reverse=True))
	

	# For output...
	detailed='Gene found' + '\t' + 'Number of instances' + '\t' + 'Cluster' + '\t' + 'Strain name' + '\t' + 'Molecule' + '\t' + 'Locus' + '\t' + 'Nb genes in cluster' + '\n'
	summary='Gene found' + '\t' + 'Number of instances' + '\t' + 'Cluster' + '\n'


	# Passe au travers du dictio 	
	for key, value in d2.items():
		if isinstance(value, dict):
			for sub_key, sub_value in value.items():

				if sub_key == 'amr_found':
					summary += sub_value + '\t'
					detailed += '\n' + sub_value + '\t'

				if sub_key == 'n':
					print('regroupement de taille', sub_value)
					summary += repr(sub_value) + '\t' + key + '\n'
					detailed += repr(sub_value) + '\t' + key + '\n'
		
				if sub_key == 'geneClusterObjects':
					for clus in sub_value:

						print(clus.amr_found, clus.read(), clus.name, clus.nb_genes)
						detailed += '\t\t' + clus.read() + '\t' + clus.strain_name + '\t' + clus.molecule_name + '\t' + clus.locus + '\t' + repr(clus.nb_genes) + '\n'


				


	with open('summary.txt', 'w') as file:
		file.write(summary)



	
	with open('detailed.txt', 'w') as file:
		file.write(detailed)



