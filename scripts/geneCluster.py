#!/usr/bin/env python
# -*- coding: utf-8 -*-

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

class geneCluster:
	"""This class define gene clusters objects characterized by :
	- their name
	- their genes
	- the orientation of the genes
	- the number of genes present
	- the locus
	- the record (in genbank format) that have been used to produce this gene cluster 
	- etc """

	def __init__(self, name):
		self.name = name
		self.nb_genes = 0		# At instanciation, the object contains no genes
		self.cluster = []
		self.strand = []
		self.target_found = []
		self.molecule_name = []
		self.strain_name = []
		self.locus = []
		self.record = SeqRecord(Seq(""))		# Create a minimal SeqRecord objects from scratch (see biopython doc)

	def add(self, gene, strand):
		"""Method to add one gene at a time. This method is largely used when parsing an annotation file. """

		self.nb_genes += 1
		self.cluster.append(gene)
		self.strand.append(strand)


	def read(self):
		"""This method allow to read the content of a cluster. """
		s = ""
		for i,j in zip(self.cluster, self.strand):

			if (j == 1 or j == '+'):
				s = s + '5\'-' + i + '-3\' '

			if (j == -1 or j == '-'):
				s = s + '3\'-' + i + '-5\' '

		return(s)


	def pretty_read(self):
		"""This method allow to read the gene content of a cluster ignoring unknown genes . """
		s = ""
		for i,j in zip(self.cluster, self.strand):
			if i == '?':
				continue
			

			if (j == 1 or j == '+'):
				s = s + '5\'-' + i + '-3\' '

			if (j == -1 or j == '-'):
				s = s + '3\'-' + i + '-5\' '

		return(s)


	def __eq__(self, other):
		"""Override the default Equals behavior"""
		if isinstance(other, self.__class__):
			return self.cluster == other.cluster and self.strand == other.strand
		return False


	def __ne__(self, other):
		"""Override the default Unequals behavior"""
		if isinstance(other, self.__class__):
			return self.cluster != other.cluster or self.strand != other.strand
		return False


	def isin(self, other):
		"""Method to test if a gene cluster is inclued in another gene cluster. The method will also test the reverse complement.
			Note that the presence of unknown genes (e.g. 5'-?-3' or 3'-?-5' will be ignored.)
		"""
		if isinstance(other, self.__class__):	# To be sure that the passed argument belong to class geneCluster

			# Create temporary objects to test the inclustion in clusters where Unknown genes have been replaced
			self_like = geneCluster('bid')
			other_like = geneCluster('bid')
			
			# proceed to the removal of unknown genes with the remove method
			self_like = self.remove('?')
			other_like = other.remove('?')		


			self_string = self_like.read()
			other_string = other_like.read()
			r_self_string = self_like.reverse_read()

			if self_string in other_string:
				return True
			else:		

				if r_self_string in other_string:
					return True
				else:
					return False
	

	def load(self, gene_list, strand_list):
		"""Method to load genes from a list"""
		self.cluster = gene_list
		self.strand = strand_list
		self.nb_genes = len(strand_list)


	def remove(self, gene_to_remove):
		"""Method to remove correctly sepcific genes from gene clusters, typically unknown genes"""
		new = geneCluster(self.name + '_edited')
		for i, gene in enumerate(self.cluster):
			if gene != gene_to_remove:
				new.add(gene, self.strand[i])
		return(new)









	def get_size(self):
		"""Method to get the number of genes included in the gene cluster. """
		return(self.nb_genes)


	def reverse(self):
		"""Method to reverse a gene cluster. """
		r_cluster=[]	# To contain the inverted cluster
		r_strand=[]		# To contain the new polarity
		
		# We use the buid-in reversed function to reverse the list
		r_cluster = list(reversed(self.cluster))
		r_strand = list(reversed(self.strand))

		# Polarity need to be changed (reverse complement)
		r_strand_2 = []
		for i in r_strand:
			if i == '+' or i == 1:
				r_strand_2.append(-1)

			if i == '-' or i == -1:
				r_strand_2.append(1)

		new = geneCluster(self.name + '_reversed')	# A new instance is created
		
		# Keep other atributes
		new.nb_genes = self.nb_genes
		new.target_found =self.target_found
		new.molecule_name = self.molecule_name
		new.strain_name = self.strain_name
		new.locus = self.locus
		new.load(r_cluster, r_strand_2)
		# Reverse complement the seqFeature object and keep some other annotations from the original record which otherwise are lost
		new.record = self.record.reverse_complement(id=self.record.id, name=self.record.name, description=self.record.description, annotations=self.record.annotations)
		# Note that by default there is no dict annotations in spliced seqFeature objects
		return(new)

		
	def reverse_read(self):
		"""A simple method to read the reverse complement of a gene cluster. """
		r_string=""		
		reverse_object = self.reverse()
		r_string = reverse_object.read()
		return (r_string)



	






