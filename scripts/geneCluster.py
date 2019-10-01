#!/usr/bin/env python
# -*- coding: utf-8 -*-

class geneCluster:
	"""This class define gene clusters objects characterized by :
	- their name
	- their genes
	- the orientation of the genes
	- the number of genes present
	- the locus
	- etc """

	def __init__(self, name):
		self.name = name
		self.nb_genes = 0		# At instanciation, the object contains no genes
		self.cluster = []
		self.strand = []
		self.amr_found = []
		self.molecule_name = []
		self.strain_name = []
		self.locus = []

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
		"""Method to test if a gene cluster is inclued in another gene cluster. The method will also test the reverse complement. """
		if isinstance(other, self.__class__):	# To be sure that the passed argument belong to class geneCluster
			self_string = self.read()
			other_string = other.read()
			if self_string in other_string:
				return True
			else:		
				r_self_string = self.reverse_read()
				if r_self_string in other_string:
					return True
				else:
					return False
	

	def load(self, gene_list, strand_list):
		"""Method to load genes from a list"""
		self.cluster = gene_list
		self.strand = strand_list
		self.nb_genes = len(strand_list)


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
		new.amr_found =self.amr_found
		new.molecule_name = self.molecule_name
		new.strain_name = self.strain_name
		new.locus = self.locus
		new.load(r_cluster, r_strand_2)
		return(new)

		
	def reverse_read(self):
		"""A simple method to read the reverse complement of a gene cluster. """
		r_string=""		
		reverse_object = self.reverse()
		r_string = reverse_object.read()
		return (r_string)



	






