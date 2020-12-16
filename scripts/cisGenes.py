#!/usr/bin/env python
# -*- coding: utf-8 -*-

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

class cisGenes:
	"""This class define cisGenes objects characterized by :
	- the presence of specific (e.g. AMR od IS) genes localized on the same molecule
	- their distance from a focal gene (bp)	
	- the number of genes present
	- the locus
	- the record (in genbank format) that have been used to retreive this information 
	- etc """

	def __init__(self, focal_gene, focal_gene_index, start, end, strand):
		self.name = focal_gene
		# At instanciation, the object contains only the focal gene
		self.nb_genes = 1
		self.genes = [(focal_gene, 'Focal_gene', '', 0, strand)]  
		self.focal_gene = (focal_gene, focal_gene_index, start, end, strand)
		self.molecule_name = []
		self.strain_name = []
		self.locus = []
		self.organism = []

	def add(self, gene_name, gene_type, bp_distance_from_focal_gene, feat_numbers_distance_from_focal_gene, strand):
		"""Method to add one gene at a time. This method is largely used when
 			parsing an annotation file.
		"""
		self.nb_genes += 1
		self.genes.append((gene_name,
							gene_type,
							bp_distance_from_focal_gene,
							feat_numbers_distance_from_focal_gene, strand))

	def read_metadata(self):
		"""This method just give the names of the strain, the locus, etc """
		print('{}\t{}\t{}\t{}'.format(self.organism, self.strain_name, self.molecule_name, self.locus), end='\n')

	def read_focal(self):
		"""This method just throw up basic infos from the focal gene """
		print('focal gene:{} index:{} start:{} end:{} strand:{}'.format(
				self.name, self.focal_gene[1], self.focal_gene[2],
				self.focal_gene[3], self.focal_gene[4]))

	def read(self):
		"""This method allow to known the amr genes and IS that lies in the same molecule (record) than the focal gene"""
		self.genes.sort(key=lambda tup: tup[3])    # Sort with the distance in features
		for gene in self.genes:
			print('\t\t\t\t', end='')
			for v in gene:
				print('{}\t'.format(v), end='')
			print('\n', end='')

	def get_results(self):
		"""This method return a string that can be concatenated with others
		similar results to produce a csv output with metadata and the cis
		genes found on all molecules that have the target gene
		"""
		s = ""
		s += self.organism[0] + '\t' + self.strain_name + '\t' + self.molecule_name + '\t' + self.locus + '\n'
		self.genes.sort(key=lambda tup: tup[3])
		for gene in self.genes:
			s += '\t\t\t'
			for v in gene:
				s += '\t' + str(v)
			s += '\n'
		return (s)		












