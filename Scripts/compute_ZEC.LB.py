#!/usr/bin/env python
# -*- coding: utf-8 -*-
__author__ = 'bastianelli'
import sys
import re
import os



# ========================================================================================================================  
# INDEXATION(GTF) :
# fonction d'indexation d'un fichier gtf en hashtables :
# argument "gtf" (str, nom de fichier)
# ========================================================================================================================  
# STRUCTURE D'UN FICHIER GTF :
# 0:chromosome    1:source    2:feature   3:start  4:end    5:score     6:strand    7:frame  8:attributes
# ---------------------------------------------------------------------------------------------------------------------------------------------------
def indexation_ex(gtf) :
   
	ref=open(gtf, 'r')    	# on ouvre le fichier gtf, en mode r=reading, dans la variable ref
	index={}          		# on crée une hashtable vide

	print "Indexation du GTF"
                 
	for line in ref :       # pour chaque élément du fichier lu, qui correspondent en fait à une ligne chacun
		l=re.split('\s+',line)     # on split la ligne au niveau des espaces
		last_elt = ("").join(l[8:])
	
		if l[2]=="exon" :

			attributes = re.split('\s*;\s*',last_elt)

			for attribute in attributes :
				gene = re.search('^gene_id"(\w+\.*\w*)"',attribute)
				if gene:
					gene_id = gene.groups()[0]


		# l = [l[0], l[3], l[4], l[6], l[1], l[2], l[8]]

		# if not index.has_key(chrom) : # si on n'a pas encore ce chromosome dans l'index
		# 	index[chrom]={}        # on l'ajoute en tant que clé, avec pour valeur une nouvelle hashtable vide
		# 	print chrom
					if not index.has_key(gene_id) :
						index[gene_id]=[]         # on stocke la ligne dans la sous-hashtable du bon chromosome, clé = transcript ID, valeur = ligne dans sa totalité
					index[gene_id].append((int(l[3]),int(l[4])))

	print "Indexation terminée, "+str(len(index))+" gènes indexés."
	#print len(index)
	#print index
	ref.close()
	return(index)

# ========================================================================================================================  
# MERGE(EXONS):
# fonction de calcul de la Zone Exonique Codante :
# union de tous les exons possibles sur un même gène, pour mesurer la longueur de la zone potentiellement codante qqsoit le transcrit
# argument "exons" : gtf indexé par la fonction "indexation_ex"
# ========================================================================================================================  
def merge(exons, outfile) :

	out = open(outfile,'w')
	for gene in exons :
		a = exons[gene]
		#print (gene, a)
		a.sort()		# sort() par défaut trie les tuples en fonction du 1er élément puis du 2nd
		#print (gene, a)

		zec = 0
		if len(a) > 1 :
			for i in xrange(0,len(a)-1) :
				if a[i+1][0] > a[i][1] :
					zec += a[i][1] - a[i][0] +1 # le stop est inclus (+1)
				else :
					if a[i+1][1] > a[i][1] :					
						zec += a[i+1][0] - a[i][0]	# le stop n'a pas besoin d'être inclus, c'est le start du suivant
					else :
						a[i+1]=a[i]
				#print a, zec
			if not (a[-1][1] < a[-2][1] and a[-1][0] > a[-2][0]) :
				zec += a[-1][1] - a[-1][0] +1 # on compte le dernier exon
												# il ne sera pas compté s'il est inclus dans l'avant dernier (> et <) ou s'il a été remplacé par l'avant dernier (=)
		else :
			zec = a[0][1]-a[0][0]
		#print zec
		out.write(str(gene)+"\t"+str(zec)+"\n")
	out.close()
	return




# ========================================================================================================================  
# ======================================================================================================================== 
# ============== MAIN ====================================================================================================
# ======================================================================================================================== 
# ========================================================================================================================
# arguments : gtf, outfile

if len(sys.argv) == 3 :
	b = indexation_ex(sys.argv[1])
	merge(b, sys.argv[2])

else :
	print "compute_ZEC.py"
	print "syntax: python compute_ZEC.py annotation.gtf outfile"
	print "		- annotation in GTF format"
	print "		- outfile = tabulated 	1: gene_id 		2: length ZEC"
