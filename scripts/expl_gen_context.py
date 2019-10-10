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


def get_gbk_files(start_dir):
	""" Given a directory, make a list of genbank files
	    This function will find all files recursively.
	"""
	gbk_files = []
	for filename in Path(start_dir).glob('**/*.gbk'):
		gbk_files.append(str(filename))
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





# Main program

if __name__ == "__main__":

	""" The program will find gbk files and will produce one features list per gbk file""" 

	
	L = []				# List of geneCluster objects	
	DockList = []		# List of dock objects. A dock is a group of related gene clusters (see geneCluster.py)!

	# Parsing arguments

	amr_in = sys.argv[1]	# There is a subtle diff between amr_in (ex. 'TEM' and amr_found (TEM-1 or TEM-212)
	folder = sys.argv[2]
	nb_genes_in_context = int(sys.argv[3])

	gbk_files = get_gbk_files(folder)		

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


	# Sorting geneCluster objects according to their size
	L2 = sorted(L, key=lambda x: x.nb_genes, reverse=True)
		

	# Step 1 ; Place the first cluster in an first Dock object
	very_first_cluster = L2.pop(0)
	my_dock = dock()
	my_dock.add(very_first_cluster)	
	DockList.append(my_dock)		# Add the dock to the dock list

	# Step 2 - All other clusters need to be grouped
	for c_index, c_value in enumerate(L2):
		
		# Try to put that cluster in an existing dock object
		for d in DockList:

			if c_value.isin(d.head):	# Check if the cluster examined can be inclued in the HEAD of a dock instance
				d.add(c_value)			# If True, the gene cluster is add to the dock

				print ('-->' + c_value.read() + ' has been add to ' + repr(d.name) + ' whose HEAD is : ' + d.head.read())
				break
		
		else:
			print ('--> Create new dock items with' + c_value.read())
			new_dock = dock()
			new_dock.add(c_value)
			DockList.append(new_dock)	# Add the dock object to the list

	print ('---Final output----')
	for dock_item in DockList:
		print (dock_item.name + '\t' + 'size:' + repr(dock_item.size) + '\t' + 'HEAD: ' + dock_item.head.pretty_read())
		for c in dock_item.elems:
			print (c.read() + '\t' + c.locus + '\t' + c.strain_name + '\t' + c.name)

	print ('---Very summarized output----')
	for dock_item in DockList:
		print (dock_item.name + '\t' + 'size:' + repr(dock_item.size) + '\t' + 'HEAD: ' + dock_item.head.pretty_read())

	print ('---Just HEAD names ----')
	for dock_item in DockList:
		print (dock_item.head.pretty_read())


