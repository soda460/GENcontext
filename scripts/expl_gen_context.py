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
#import pandas as pd
import numpy
import re
import csv
import glob
import argparse

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



def examine_gbk_file (gene, gbk_f, card_gene):
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
		abc = check_gene_presence(gene, gbk_rec, card_gene)

		hits.extend(abc)

	return (hits)


def check_gene_presence(gene, record, card_gene):
	""" Check the presence of the AMR gene given in argument of the script 
		This should be a gene referenced in th CARD database and annotated in an inference fied by PROKKA
		This function return a 4-items list : an index, the feature, the AMR gene that have been found, and the entire record
	"""

	record_out= []

	if card_gene == 'IS':
		#print ('Dealing with an IS')
		# implement...
		return

	if card_gene == 'normal':
		#print ('Dealing with a normal gene')
		# implement...
		return

	if card_gene == 'card':


	

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
	"""This function get features around a given target feature index ignoring unknown genes in the same way than the
	the extract_gene_cluster function	
	"""
	nb_features = len(record.features)
	#print("inside fn \n nb feats:" + repr(nb_features) + ' ' + record.name)

	''' The following code will count 'true' genes after the targeted AMR gene (index of record is given actually)'''
	right_index = 0
	left_index = 0
	upstream_genes = 0
	downstream_genes = 0


	# Block A: Searching downstream of the targeted index
	for i_feat, feature in enumerate(record.features[index + 1:], start=index + 1):

		if downstream_genes >= features_range:		# End condition		
			right_index = i_feat -1					# Because it iterate one more time before break
			break

		if feature.type == 'CDS':					# We focus on high quality CDS
			flagGeneDraw = 'FALSE'
			#print('ifeat:' + repr(i_feat) + ' feature: ' + repr(feature))

			# Block A1 'True' genes are grabed if annotated by CARD or labeled with ISfinder
			if 'inference' in feature.qualifiers:					
				for inf in feature.qualifiers['inference']:		# In some cases, more than one inference tag are present in gbk file.
					if "protein_fasta_protein_homolog_model" in inf:
						flagGeneDraw = 'TRUE'
						downstream_genes += 1
						break
			
					if "ISfinder" in inf:		
						flagGeneDraw = 'TRUE'
						downstream_genes += 1
						break

			# Otherwise, continue searching for high-quality CDS, which have a \gene tag
			if flagGeneDraw == 'FALSE':
				if 'locus_tag' in feature.qualifiers:
					if 'gene' in feature.qualifiers:
						flagGeneDraw = 'TRUE'
						downstream_genes += 1

	# Block B: Searching upstream of the targeted index
	''' We use enumarate on the elements before the targeted genes, then the built-in reversed fn and the orginal iterators!'''
	for i_feat, feature in reversed(list(enumerate(record.features[:index -1] ))):

		if upstream_genes >= features_range:		# End condition		
			left_index = i_feat +1	# because the list is reversed and that the indices are the orginal ones 
			break
	
		if feature.type == 'CDS':
			flagGeneDraw = 'FALSE'
			#print('ifeat:' + repr(i_feat) + ' feature: ' + repr(feature))
			if 'inference' in feature.qualifiers:					
				for inf in feature.qualifiers['inference']:		# In some cases, more than one inference tag are present in gbk file.
					if "protein_fasta_protein_homolog_model" in inf:
						flagGeneDraw = 'TRUE'
						upstream_genes += 1
						break
			
					if "ISfinder" in inf:		
						flagGeneDraw = 'TRUE'
						upstream_genes += 1
						break

			# Continue searching if not found in the last block
			if flagGeneDraw == 'FALSE':
				if 'locus_tag' in feature.qualifiers:
					if 'gene' in feature.qualifiers:
						flagGeneDraw = 'TRUE'
						upstream_genes += 1

	'''
	To access to nucletotide coords of the feature; do : 
	print (record.features[index].location.start)
	https://biopython.org/DIST/docs/api/Bio.SeqFeature.FeatureLocation-class.html
	Other qualifiers for location are : start, end, strand
	Note that the start and end location numbering follow Python's scheme, thus a GenBank entry of 123..150 (one based counting)
	becomes a location of [122:150] (zero based counting).
	'''
	lower_pos = record.features[left_index].location.start
	upper_pos = record.features[right_index].location.end

	s = "Original record will be spliced at positions " + repr(lower_pos) + ":" + repr(upper_pos)
	#print(s)

	sub_record = record[lower_pos:upper_pos]
	return (sub_record)







def extract_gene_cluster (record):
	"""This function will give gene cluster and orientations
	The function will also put the record into an attribute of the cluster object
	"""

	# Use our own class geneCluster

	a = geneCluster(record.name)	
	a.record = record	# We put a seqRecord object in our cluster object
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



def check_arg_orientation(gene_cluster, amr_found):
	"""This function will check if the antibiotic resistance gene of the given gene cluster is in the 5' 3' orientation."""

	s='5\'-' + amr_found + '-3\''

	if s in gene_cluster.read():
			return	gene_cluster	# Already in the right orientation
	else:
		if s in gene_cluster.reverse_read():
			reversed_gene_cluster = gene_cluster.reverse()
			return reversed_gene_cluster
		else:
			print ('This should not happen!')


def add_blank_tracks(gd_diagram, n):
	""" This function will create additionnal blank tracks to equilibrate the diagrams.
		When few track are present, this prevent tracks to use a lot of vertical space in the page.
	"""
	for i in range (n):
		gd_track_for_features = gd_diagram.new_track(1, name='void', greytrack=False, start=0, height=1, end = 0, scale=1, hide=1)

		# Create an empty set of features. The newly created set is linked to the abovementionned track
		gd_feature_set = gd_track_for_features.new_set()
	return


def recursive_draw_with_genomeDiagram(dock, clusters_per_page):
	""" This recursive function will draw one genome diagram each until all gene clusters are draw."""
	if dock.draw >= dock.size:	# End condition
		return

	cwd = os.getcwd()
	# Create alphabet list of uppercase letters
	alphabet = []
	for letter in range(65, 91):
		alphabet.append(chr(letter))

	fn_colors = prepare_color_dict()	# a small dictionnary that associate colors to gene categories
	iterare = round(dock.draw/clusters_per_page) + 1	# to known how many times the function have been called

	''' Get the next clusters that were not already draw. They are stored (in reverse order of size) in the dock object'''
	next_clusters = []
	next_clusters = dock.elems[dock.draw:(dock.draw + clusters_per_page)]

	# Creation of the genomeDiagram object for that call
	gd_diagram = GenomeDiagram.Diagram(dock.name + '_' + repr(iterare))

	# For each geneClusters, we will create a track, add features, etc!
	for geneCluster in next_clusters:

		name = geneCluster.record.annotations["source"]		# !need to fix this to ensure that we have a good label above each track!
	
		''' It is important to get the longest gbk record to keep proportionarity between clusters present on different diagrams
		but belonging to the same dock. It is why, this value is put into a dock attribute.
		'''
		dock.max_cluster_len = max(dock.max_cluster_len, len(geneCluster.record))

		# Adding a track; each track being a different geneCluster
		gd_track_for_features = gd_diagram.new_track(1, name=geneCluster.record.name, greytrack=True, start=0, height=1, end = len(geneCluster.record), scale=1, scaleticks=0)

		# Create an empty set of features linked to the newly created track
		gd_feature_set = gd_track_for_features.new_set()

		''' This is our home-made fn to fill the features set with stuff like genes colored according to categories,
		specific label for AMR genes, etc! Inside the function digest_feature, there are many calls to gd_feature_set.add_feature
		'''
		gd_feature_set = digest_features(geneCluster.record, fn_colors, gd_feature_set)

	# Calculate the number of needed blanck tracks and add them
	n_blank_tracks = clusters_per_page - len(next_clusters)
	add_blank_tracks(gd_diagram, n_blank_tracks)

	gd_diagram.draw(format="linear", orientation="landscape", fragments=1, start=0, end=dock.max_cluster_len, track_size=0.75, tracklines=0)
	gd_diagram.write(cwd + '/output_' + target_gene + '/maps/' + dock.name  + '_page' + repr(iterare) + '_linear.pdf', "PDF")

	# Update the dock.draw counter
	dock.draw += clusters_per_page

	recursive_draw_with_genomeDiagram(dock, clusters_per_page)




if __name__ == "__main__":

	""" The program will find gbk files and will produce one features list per gbk file""" 

	# Construct the argument parser
	ap = argparse.ArgumentParser()

	# Add the arguments to the parser
	ap.add_argument("-t", "--target_gene", required=True, help="Gene targeted by the program")
	ap.add_argument("-c", "--card", required=True, help="Boolean indicating if the target gene is referenced in CARD database")
	ap.add_argument("-p", "--path", required=True, help="PATH of the folder containing the data ")
	ap.add_argument("-n", "--genes_around_target", required=True, help="Integer specifying the number of genes explored around the target")

	# Parsing arguments
	args = vars(ap.parse_args())	# args a dict
	target_gene = args['target_gene']	# There is a subtle diff between target_gene (ex. 'TEM' and target_found (TEM-1 or TEM-212)
	folder = args['path']
	card_gene = args['card']
	nb_genes_in_context = int(args['genes_around_target'])

	# Declaration of variables
	L = []				# List of geneCluster objects	
	DockList = []		# List of dock objects. A dock is a group of related gene clusters (see geneCluster.py)!

	# Data treatment 
	gbk_files = get_gbk_files(folder)		

	for f in gbk_files:
		hits = examine_gbk_file(target_gene, f, card_gene)

		'''Note that the hits list is a list of lists, which contains distinct objects
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
				target_found = my_meta_stuff.pop(0)
				gbk_file_name = my_meta_stuff.pop(0)
				molecule_name = my_meta_stuff.pop(0)
				strain_name = my_meta_stuff.pop(0)
				gbk_record_name = my_meta_stuff.pop(0)
				locus = my_meta_stuff.pop(0)

				# The index of the targeted gene, i.e. the feature position of the  targeted gene in the record
				target_gene_index = my_hit[0]
				
				# Correspond to the complete record
				target_record = my_hit[3]

				# Splice the record where an AMR gene has been found
				sub_rec = splice_record(target_gene_index, target_record, nb_genes_in_context)
				
				""" As mentionned in the Biopython tutorial : "The SeqRecord slicing step is cautious in what annotation it preserves [...]
				To keep the database cross references or the annotations dictionary, this must be done explicitly"
				"""
				sub_rec.annotations = target_record.annotations.copy()

				# Obtain a GeneCluster object from a subrecord 
				my_gene_cluster = extract_gene_cluster(sub_rec)

				# Add metadata to our gene cluster object
				my_gene_cluster.locus=locus
				my_gene_cluster.target_found=target_found
				my_gene_cluster.molecule_name=molecule_name
				my_gene_cluster.strain_name=strain_name

				# if not put my gene cluster in the proper orientation relative to target_gene
				returned_cluster = check_arg_orientation(my_gene_cluster, target_found)

				# Simply put the GeneCluster object in a list
				L.append(returned_cluster)


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

				#print ('-->' + c_value.read() + ' has been add to ' + repr(d.name) + ' whose HEAD is : ' + d.head.read())
				break
		
		else:
			#print ('--> Create new dock items with' + c_value.read())
			new_dock = dock()
			new_dock.add(c_value)
			DockList.append(new_dock)	# Add the dock object to the list

	# Make folders if folder for strain does not exists
	cwd = os.getcwd()

	# Create output directory
	if not os.path.exists(cwd + '/output_' + target_gene):
		os.mkdir(cwd + '/output_' + target_gene)


	if not os.path.exists(cwd + '/output_' + target_gene + '/maps'):
		os.mkdir(cwd + '/output_' + target_gene + '/maps')


	""" We will pass a dock object to a fn that will draw nice diagrams of the geneClusters present in that dock"""
	for dock_item in DockList:
		recursive_draw_with_genomeDiagram(dock_item, 8)		# 8 a good value for letter

	# Output section
	
	# First file -- Detailed output
	print ('---Detailled output----')
	f1 = open(cwd + '/output_' + target_gene + '/detailed_output.csv','w')
	outlines = 'amr\tdock_name\tdock_size\tStrain\tMolecule\tLocus\tReversed\tCluster\tCluster(Pretty) \tdock_HEAD(Pretty)\n'

	for dock_item in DockList:
		for c in dock_item.elems:
			for i in [
				c.target_found, dock_item.name, repr(dock_item.size),
				c.strain_name, c.molecule_name, c.locus, dock_item.reversed, c.read(), c.pretty_read(),
				dock_item.head.pretty_read()
				]:
					outlines += i + '\t'
			outlines += '\n'
	f1.write(outlines)
	f1.close()

	# Second file -- Brief summary
	print ('---Very summarized output----')
	f2 = open(cwd + '/output_' + target_gene + '/summary.csv','w')
	summary = "Dock_name\tDock_size\tdock_HEAD\n"

	for dock_item in DockList:
		print (dock_item.name + '\t' + 'size:' + repr(dock_item.size) + '\t' + 'HEAD: ' + dock_item.head.pretty_read())
		summary += dock_item.name + '\t' + repr(dock_item.size) + '\t' + dock_item.head.pretty_read() + '\n'	
	f2.write(summary)
	f2.close()













