#!/usr/bin/env python
# -*- coding: utf-8 -*-

class geneCluster:
	"""Classe définissant une regroupement de genes caractérisé par :
	- son nom
	- ses genes
	- leur orientation
	- le nombre de genes"""

	def __init__(self, name):				# Notre méthode constructeur		
		self.name = name
		self.nb_genes = 0 					# Au depart il est vide
		self.cluster = []
		self.strand = []
		self.amr_found = []
		self.molecule_name = []
		self.strain_name = []
		self.locus = []
		self.parent = ""	# Pour garder la trace des inclusions dans d'autres gene clusters


	def add(self, gene, strand):
		"""Methode pour ajouter un gene a la fois """

		self.nb_genes += 1
		self.cluster.append(gene)
		self.strand.append(strand)


	def read(self):
		"""Une methode pour lire le contenu du cluster"""
		s = ""
		#s = self.name + ": "
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
		"""Methode pour tester si un regroupement de gene est inclus dans un autre"""
		l3=[]	# a list to store list index :)

		if isinstance(other, self.__class__):
			"""Pour s'assurer que le parametre passe a la fn appartienne a la classe geneCluster"""

			for i, val in enumerate(self.cluster):
				if val in other.cluster:
			
					# Les indices des element du deuxieme cluster (other) communs a ceux du premier gene cluster sont gardes
					l3.append(other.cluster.index(val))					
				else:
					return False

			# Here we test if the indices are consecutive 
			if ( (l3[-1] - l3[0]) == (len(l3) - 1) ):
				self.parent = other.read()
				return True	
			else:
				# Try again with the reversed small cluster 
				l3=[]				
				rev_small = self.reverse()
				for i, val in enumerate(rev_small.cluster):
					if val in other.cluster:
						l3.append(other.cluster.index(val))
					else:
						return False
				# Here we test if the indices are consecutive 
				if ( (l3[-1] - l3[0]) == (len(l3) - 1) ):
					self.parent = other.read()
					return True
				else:
					return False

	
	def load(self, gene_list, strand_list):
		"""Methode pour loader des genes from scratch a partir d'une liste"""
		self.cluster = gene_list
		self.strand = strand_list

	def get_size(self):
		"""Methode pour avoir le nombre de genes"""
		return(self.nb_genes)


	def reverse(self):
		"""Methode pour inverser un regroupement de genes"""
		r_cluster=[]	# to hold the inverted cluster
		r_strand=[]		# to hold the new polarity
		
		# use the reversed function to reverse the list
		r_cluster = list(reversed(self.cluster))
		r_strand = list(reversed(self.strand))

		# polariyy need to be changed (reverse complement)
		r_strand_2 = []
		for i in r_strand:
			if i == '+' or i == 1:
				r_strand_2.append(-1)
			if i == '-' or i == -1:
				r_strand_2.append(1)

		new = geneCluster(self.name + '_reversed')
		

		# keep other atributes
		new.nb_genes = self.nb_genes
		new.amr_found =self.amr_found
		new.molecule_name = self.molecule_name
		new.strain_name = self.strain_name
		new.locus = self.locus
		new.parent = self.parent	# Pour garder la trace des inclusions dans d'autres gene clusters
	
		new.load(r_cluster, r_strand_2)

		return(new)

		
	def reverse_read(self):
		"""Une methode pour lire le reverse complement du cluster"""
		r_string=""		
		reverse_object = self.reverse()
		r_string = reverse_object.read()
		return (r_string)



