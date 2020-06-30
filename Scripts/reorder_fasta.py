from Bio import SeqIO
import sys

handle = open(sys.argv[1], "rU")
l = SeqIO.parse(handle, "fasta")
sortedList = [f for f in sorted(l, key=lambda x : x.id)]
for s in sortedList:
   print '>' + s.description
   print str(s.seq)
