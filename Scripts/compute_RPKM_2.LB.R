#!/usr/bin/env Rscript
#Author: Leila Bastianelli

###########################################################################################################
# mRNA : HTSeq-count
# ZEC: compute_ZEC.py
###########################################################################################################

args <- commandArgs(TRUE)

zec<-read.csv(args[1],h=F,sep="\t",stringsAsFactors=FALSE)
rna<-read.csv(args[2],h=F,sep="\t",stringsAsFactors=FALSE)
xp<-args[3]
#million<-as.numeric(args[4])

# zec<-read.csv("~/Desktop/test.out",h=F,sep="\t",stringsAsFactors=FALSE)
# rna<-read.csv("~/Desktop/Projets/MBMT/HTSeq-count_ensembl75/data/RNA_ReadsPer10M/RNA_MB.mean.Counts.RP10M.txt",h=F,sep="\t",stringsAsFactors=FALSE)
# xp<-"MB"

# just to check
# zec <- zec[ order(zec$V1), ]

a<-merge(zec,rna,by.x=1,by.y=1)

names<- c("GeneID","ZEC", "mean_RNA_counts")
colnames(a)<-names
a$expression <- a$mean_RNA_counts/a$ZEC

# names<- c("GeneID","ZEC", "mean_RNA_counts", "expression")
# colnames(a)<-names

#head(a$mean_RNA_counts)

reads <- sum(a$mean_RNA_counts)


a$expression <- (a$expression*1000)
a$expression <- a$expression/reads
a$expression <- a$expression*1000000


matrice<-data.frame(a$GeneID, a$expression)

# Library size
#print(reads)
write.table(matrice, paste(xp,".RPKM.txt",sep="") ,sep="\t",col.names = FALSE, row.names = FALSE, quote = FALSE)
