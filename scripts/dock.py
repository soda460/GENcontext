#!/usr/bin/env python
# -*- coding: utf-8 -*-

class dock:
	"""Classe définissant un regroupement caractérisé par :
	- son nom
	- un HEAD, i.e. l'element dont un des attributs est le plus grand dans la liste (.size)"""


	"""Cette classe possède un attribut de classe qui s'incrémente à chaque fois que l'on crée un objet de ce type"""

	dock_id = 0

	



	def __init__(self):				# Notre méthode constructeur

		dock.dock_id += 1

		self.name = 'cluster_' + repr(dock.dock_id)
		self.head = []
		self.elems = []
		self.size = 0
		
	def add(self, geneCluster):
		"""Methode pour ajouter un regroupement de genes a ce dock. Un tri se fait apres l'ajout
		pour que les plus gros regroupement soient au debut (head). Cela empeche d'avoir a tester
		tous les regroupements de gene. Si un regroupement est inclus dans le head (ou son inverse)
		il le sera aussi dans tous les representants du dock """

		L = []	# We will use this tp put the reverse sort result 
			
		self.size += 1
		self.elems.append(geneCluster)

		L = sorted(self.elems, key=lambda x: x.nb_genes, reverse=True)		
		self.elems = L
		self.head = L[0]
		return
