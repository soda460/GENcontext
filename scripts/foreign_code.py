#!/usr/bin/env python
# -*- coding: utf-8 -*-


import sys
from Bio import SeqIO
from reportlab.lib import colors
from reportlab.lib.units import cm
from Bio.Graphics import GenomeDiagram
#from Bio import SeqIO


# Declaration des fonctions

def prepare_color_dict():

	""" Cette fonction permet de preparer un dictionnaire permettant d'associer un nom de gene a une couleur (groupe fonctionnel)
		On a les groupes suivants :
		i) antibiotic resistance : n'apparaît pas ici, car il est hardcode ailleurs dans le script (tag megares)
		ii) conjugation | yellow
		iii) replication/transposase | blue
		iv) maintenance and stability | green
		v) pil genes | lightblue
		vi) shufflon | orange"""

	fn_colors = {}		# Dictionnaire pour associer les genes a une couleur refletant la categorie fonctionnelle; voir fn prepare_color_dict


	#			 0					 4						  9						  14					   19				        24    25
	alphabet = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z"]


	pil = alphabet[8:22] + ["VA"]	# les pil genes connus vont de pilI a pilV et on a aussi un gene "pilVA"  
	
	# lightblue genes

	for p in pil:
		fn_colors["pil" + p] = colors.lightblue


	# pink genes

	tra = alphabet[:-1]				# les genes tra vont de traA a traX

	for i in tra:
		fn_colors["tra" + i] = colors.pink

	trb = alphabet[:3]				# traA, traB et traC


	for i in trb:
		fn_colors["trb" + i] = colors.pink


	otherPinkGenes = ["sogL", "sogS", "yggA", "nikA", "nikB"]

	for i in otherPinkGenes:	
		fn_colors[i] = colors.pink


	# orange genes

	fn_colors["rci"] = colors.orange	# un seul gene dans cette categorie


	# blue genes

	blueGenes = ["repZ", "repC", "ISEcp1", "Tn1721", "IS1S", "tnpA", "intI1", "IS5075", "IS5", "Tn2", "hin", "TnAs1", "ISSbo1"]


	for i in blueGenes:
	
		fn_colors[i] = colors.blue


	# green genes

	greenGenes = ["ibfA", "ydeA", "ydfA", "ydfB", 
				  "mck", "kor",
				  "ydiA", "ydjA", "pifA",
				  "yefA","yegA", "parA", "parB", "impB", "impA", "impC",
				  "yfaA","yfaB", "yfbA", "yfbB", "yfcA", "yfcB", "yfeA","yfeB","yfeC",
				  "yfhA", "psiB", "psiA", "ardA"]


	for i in greenGenes:
		
		fn_colors[i] = colors.green

	return fn_colors



def digest_features(record, fn_colors, gd_feature_set):

	for feature in record.features:

		if feature.type == 'CDS':		# We use the CDS fields

			flagGeneDraw = 'FALSE'		# Initialize this 'Boolean'!

			#print ("\n#### We have in our hands a CDS feature", "at locus", feature.qualifiers['locus_tag'][0], "####\n")	# 

			# We removed it because PROKKA doesn't always add a /product tag!
			#if 'product' in feature.qualifiers: # Verify it codes for a real protein (not pseudogene)

			if 'inference' in feature.qualifiers:
					
				for inf in feature.qualifiers['inference']:		# In some cases, more than one inference tag are present in gbk file.
					
					# Case 1 : a resistance gene tagged in an inference field
					# We are looking for resistance genes;
					# they were annotated using megares db with PROKKA.
					# Therefore they be found in our gbk files by searching for the megares tag in the inference field
	      			#         /locus_tag="MECPGFAM_00008"
           			#         /inference="ab initio prediction:Prodigal:2.6"
           			#         /inference="similar to AA
           			#         sequence:megares.pep:Tmt|DfrA17|AB126604|98-
           			#         571|474|Trimethoprim|Dihydrofolate_reductase|DHFR|Requires
           			#         SNPConfirmation"
            		#         /codon_start=1
            		#         /transl_table=11
            		#         /product="hypothetical protein"

					if "protein_fasta_protein_homolog_model" in inf:
							
						inf2 = inf.replace(" ", '') 	# To remove extra space converted from \n in genbank

						b = inf2.split("|")				# To get the gene name, 

						#megares_gene = (b[1])			# located after the first pipe, e.g. : sequence:megares.pep:Tmt|DfrA17|AB126604|98-
						card_gene = (b[3])			# located after the 3rd pipe, e.g. : similar to AA sequence:protein_fasta_protein_homolog_model.fasta:gb
															# |CAA63262.1|ARO:3001864|CTX-M-1"
							
						gd_feature_set.add_feature(feature,
													sigil="ARROW",name = card_gene,
													color=colors.red,
													label_size=8,
													label_color=colors.red,
													label_angle=90,
													label=True)

						# print ("We have found a CARD AMR gene and will draw its CDS in RED with gene label", card_gene)

						flagGeneDraw = 'TRUE'

						break
			
					if "ISfinder" in inf:
							
						inf2 = inf.replace(" ", '') 	# To remove extra space converted from \n in genbank

						b = inf2.split(":")				# To get the gene name, 

						isFinder_gene = (b[2])			# located after the second : e.g. : similar to AA sequence:ISfinder:ISVsa5"

						gd_feature_set.add_feature(feature,
													sigil="ARROW",name = isFinder_gene,
													color=colors.purple,
													label_size=8,
													label_angle=90,
													label_color=colors.purple,
													label=True)

						#print ("We have found a ISfinder gene (ISfinder from PROKKA db) and will draw its CDS in purple with gene label", ISfinder_gene)
						
						flagGeneDraw = 'TRUE'

						break





				# To control the flow process; we want to continue searching if flagGeneDraw != TRUE

				if flagGeneDraw == 'TRUE':
					continue
			


				# Normally there is always a locus_tag in features.qualifiers
				# Case 2 : the clearest case : a CDS with a gene tag in qualifiers

				if 'locus_tag' in feature.qualifiers:

					#print('we still work with', feature.qualifiers['locus_tag'][0])

					if 'gene' in feature.qualifiers:

						# data in qualifiers are all lists, even if only 1 string, so [0] converts to string
        				# use lower() to make sure weirdly capitilized strings get selected as well
						# my_gene = feature.qualifiers['gene'][0].lower()
						# interesting, but I want my genes like rpoA

						my_gene = feature.qualifiers['gene'][0]


						# coloring the gene accroding to colors defined in a dictionnary


						# avec la methode get, va chercher la colour associe au gene dans le dict fn_colors, si absent : gris
						color = fn_colors.get(my_gene, colors.gray)


						# Adding this gene to our feature set
			

						gd_feature_set.add_feature(feature,
											sigil="ARROW",
											color = color,
											label_color = color,
											label_size=8,
											label_angle=90,
											label=True)

						#print ("We will draw a CDS with gene label", my_gene)
				
						continue

					# When CDS does not contains the gene tag, we want to put them on the map in gray
		
					else:

						#print ("We do not found a CDS with a gene label for that CDS feature")

						# Adding this gene to our feature set

						gd_feature_set.add_feature(feature, name="",
											sigil="ARROW",
											color=colors.gray,
											label_size=8,
											label_angle=90,
											label=True)

						#print ("We will draw a CDS in gray with no gene name")

	return gd_feature_set


def get_coord_gene(record, geneName):

	for feature in record.features:

		if feature.type == 'gene':

			if 'gene' in feature.qualifiers:
	
				if geneName in feature.qualifiers['gene']:

					location = feature.location.start
					break

		

	return location








