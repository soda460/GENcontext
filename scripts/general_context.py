#!/usr/bin/env python

"""
This command list antibiotic resistance genes referenced in the CARD database
as well as IS that co-localized on the same molecule (e.g. a plasmid).
As input the user should give a 'focal_gene' that is searched in 
annotation files and from which relative distances are computed

Usage : ./general_context.py -t 'CTX-M-1' -c 'card' -p '../mini_prokka/'

"""
__auteur__ = "Jean-Simon Brouard"
__date__ = "2020-02-10"

# Importation of standard modules
import sys
import os
from pathlib import Path
from Bio.Graphics import GenomeDiagram
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from io import StringIO
from collections import OrderedDict
import numpy
import re
import csv
import glob
import argparse

# Importation of local modules
from foreign_code import prepare_color_dict
from foreign_code import digest_features
from foreign_code import get_coord_gene
from geneCluster import geneCluster
from cisGenes import cisGenes
from dock import dock

# Importation of functions from the main program
from expl_gen_context import get_gbk_files, examine_gbk_file, \
							 get_metadata_and_record, \
							 check_excluded_terms


def fill_cis_genes(cis_genes, record):
	"""This function will return all amr (card-labeled) and IS genes
	that are present on the record. The positions are given relative to the
	most proximal one of the focal gene.
	"""

	# Adding the organism info to the object
	for feature in record.features:
		if feature.type == 'source':
			if 'organism' in feature.qualifiers:				
				cis_genes.organism = feature.qualifiers['organism']
			break

	index = cis_genes.focal_gene[1]

	# Searching the gbk record *downstream* of the targeted index
	for i_feat, feature in enumerate(record.features[index + 1:], start=index + 1):
		if feature.type == 'CDS':
			# Block A1 'True' genes are grabbed if annotated by CARD or labeled with ISfinder
			if 'inference' in feature.qualifiers:					
				for inf in feature.qualifiers['inference']:		# In some cases, more than one inference tag are present in gbk file.
					if "protein_fasta_protein_homolog_model" in inf:
						inf2 = inf.replace(' ', '') 	# To remove extra space converted from \n in genbank
						b = inf2.split('|')				# To get the gene name, 
						my_gene = (b[3])
						# Distance between the end of the focal gene and this one
						bp_distance = feature.location.start - cis_genes.focal_gene[3]
						feature_distance = i_feat - index
						cis_genes.add(my_gene, 'card', '+' + repr(bp_distance), feature_distance, feature.strand)	
						break
					if "ISfinder" in inf:
						inf2 = inf.replace(' ', '')
						b = inf2.split(":")
						my_gene = (b[2])
						bp_distance = feature.location.start - cis_genes.focal_gene[3]
						feature_distance = i_feat -index
						cis_genes.add(my_gene, 'IS', '+' + repr(bp_distance), feature_distance, feature.strand)
						break

	# Searching the gbk record *upstream* of the targeted index
	''' We use enumarate on the elements before the targeted genes, then the built-in reversed fn and the orginal iterators!'''
	for i_feat, feature in reversed(list(enumerate(record.features[:index -1] ))):
		if feature.type == 'CDS':
			if 'inference' in feature.qualifiers:					
				for inf in feature.qualifiers['inference']:
					if "protein_fasta_protein_homolog_model" in inf:
						inf2 = inf.replace(' ', '') 	
						b = inf2.split('|')
						my_gene = (b[3])
						# Distance between the start of the focal gene and this one
						bp_distance = cis_genes.focal_gene[2] - feature.location.end
						feature_distance = index - i_feat
						cis_genes.add(my_gene, 'card', '-' + repr(bp_distance), -(feature_distance), feature.strand)	
						break
			
					if "ISfinder" in inf:
						inf2 = inf.replace(' ', '')
						b = inf2.split(":")
						my_gene = (b[2])
						bp_distance = cis_genes.focal_gene[2] - feature.location.end 
						feature_distance = index - i_feat
						cis_genes.add(my_gene, 'IS', '-' + repr(bp_distance), -(feature_distance), feature.strand)
						break
	return (cis_genes)

if __name__ == "__main__":

	# Construct the argument parser
	ap = argparse.ArgumentParser()

	# Add the arguments to the parser
	ap.add_argument("-t", "--target_gene", required=True, help="Gene targeted by the program")
	ap.add_argument("-c", "--card", required=True, help=" Gene category. \
						    Two choices: card, which indicate that the \
							target gene is referenced in CARD database \
							and IS, for Insertion Sequence")
	ap.add_argument("-p", "--path", required=True, help="PATH of the folder containing the data ")
	ap.add_argument("-e", "--exclude", required=False, nargs='+', help="To exclude genbank files containing a specific term (e.g. chromosome")

	# Parsing arguments
	args = vars(ap.parse_args())	# args is a dict
	target_gene = args['target_gene']
	folder = args['path']
	card_gene = args['card']
	excluded_terms = args['exclude']

	# Declaration of variables
	L = []				# List of geneContext objects	
	DockList = []		# List of dock objects. A dock is a group of related gene clusters (see geneCluster.py)!
	cwd = os.getcwd()
	header = 'Organism\tStrain\tMolecule\tLocus\tGene\tGenetype\tDistanceFromFocalGeneBp\tDistanceFromFocalGeneFeatures\tStrand\n'
	outlines = ""

	# Data treatment 
	gbk_files = get_gbk_files(folder)

	if not os.path.exists(cwd + '/general_context'):
		os.mkdir(cwd + '/general_context')
	if not os.path.exists(cwd + '/general_context/' + target_gene):
		os.mkdir(cwd + '/general_context/' + target_gene)

	# Prepare the file for output
	f1 = open(cwd + '/general_context/' + target_gene + '/detailed_output.csv','w')
	f1.write(header)

	for f in gbk_files:
		if not check_excluded_terms(excluded_terms, f):
			hits = examine_gbk_file(target_gene, f, card_gene)
			'''Note that the hits list is a list of lists, which contains distinct objects
				[0], is the index of the amr_gene, i.e. the feature position of the amr_gene in the record
				[1], is the feature of the amr_gene 
				[2], is the amr_gene (string)
				[3], is the complete record
				[4], is a FeatureLocation object containing the start, end and strand of the amr_gene
			'''
			if hits:
				for my_hit in hits:
					my_meta_stuff = []
					my_meta_stuff = get_metadata_and_record(f, my_hit)
					target_found = my_meta_stuff.pop(0)
					gbk_file_name = my_meta_stuff.pop(0)
					molecule_name = my_meta_stuff.pop(0)
					strain_name = my_meta_stuff.pop(0)		# Will be override later by  the extract_cis_genes fn
					gbk_record_name = my_meta_stuff.pop(0)
					locus = my_meta_stuff.pop(0)
					target_gene_index = my_hit[0]    # The feature position of the targeted gene in the record
					target_record = my_hit[3]        # Correspond to the complete record
					FeatLoc = my_hit[4]			     # A tuple with SeqFeature Location infos (start, end, strand)

					# Instantiate a new object belonging to the new class cisGenes, 
					my_cis_genes = cisGenes(target_found, target_gene_index, FeatLoc[0], FeatLoc[1], FeatLoc[2])	

					# Add some metadata fields to our cisGenes object.
					# The organism field will be fill later.
					my_cis_genes.locus=locus
					my_cis_genes.molecule_name=molecule_name
					my_cis_genes.strain_name=strain_name

					# Now use a function to gently parse the actual record and populate the cis_genes object
					b = fill_cis_genes(my_cis_genes, target_record)
					b.read_metadata()
					b.read()
					outlines += b.get_results()
					L.append(b)
	if not L:
		print ("The target gene was not found")
		sys.exit()
	f1.write(outlines)

					


			







