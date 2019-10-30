#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Jean-Simon Brouard"
__date__ = "2019-10-01"

class dock:
	"""This simple class define dock objects which can be used to regroup similar gene clusters objects.
	This class is defined by
	- its name
	- an HEAD, which is simply the geneCluster with the greatest number of genes"""

	dock_id = 0		# This is a class attribute which is incremented when an instance of a dock object is created

	def __init__(self):
		dock.dock_id += 1
		self.name = 'dock_' + repr(dock.dock_id)
		self.head = []
		self.elems = []		# All related gene clusters objects will 
		self.size = 0
		self.reversed = 'no'	# by default the gene cluster is not reversed!
		self.draw = 0			# just to count the nunmber of draw items
		self.max_cluster_len = 0	# useful to draw clusters proportionately
		
	def add(self, geneCluster):

		"""This method allow to add a gene cluster to the dock. After the addition of a gene cluster, the list which contains the cluster
			objects is reversely sorted according to the number of their genes. The cluster at the begining of the list become the new HEAD.
			This is extremely useful because to regroup gene clusters by relatedness, we just have to test the inclusion of a gene cluster
			in the HEAD. """

		L = []
		self.size += 1
		self.elems.append(geneCluster)
		L = sorted(self.elems, key=lambda x: x.nb_genes, reverse=True)		
		self.elems = L		# Replace the list in the object by the reversed sorted list
		self.head = L[0]	# Place the HEAD
		return


	# We could add method get_name, get_size and pretty_read in the future.
