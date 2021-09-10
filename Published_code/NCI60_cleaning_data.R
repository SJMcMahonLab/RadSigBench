## ---------------------------
##
## Script name: NCI60_cleaning_data.R
##
## Purpose of script: To clean NCI60 microarray expression data and match with 
## radiosensitvity metrics
##
## Author: Dr. John O'Connor
##
## Date Created: 03-09-2021
##
## Email: john.oconnor@qub.ac.uk
##
## ---------------------------
##
## Notes: 
##
##
##
## ---------------------------

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("GEOquery")
#BiocManager::install("affy")
#BiocManager::install("hgu133plus2.db")
#BiocManager::install("annotate")
#BiocManager::install("limma")

library(GEOquery)
library(affy)
library(hgu133plus2.db)
library(annotate)
library(limma)
library(stringr)

homedir="C:/Users/HOMEDIRECTORY/"
dir.create(paste(homedir,"Cleaned_data/", sep=""))
setwd(paste(homedir, "Raw_data/", sep=""))
getGEOSuppFiles("GSE32474")
setwd("GSE32474/")
untar("GSE32474_RAW.tar")

##GEO dataset has 3 repeats, use only first set. 
fileids=sprintf("GSM8036%s",seq(15,73))
masterlist=list.files(pattern="GSM")
shortlist=substr(masterlist, 1,9)
m=match(shortlist, fileids) 
unlink(masterlist[is.na(m)])

##RMA normalization
dat <- ReadAffy()
eset=rma(dat)
data=data.frame(exprs(eset))

##read labels and SF2
a=read.table(paste(homedir,"Raw_data/SF2_and_labels.csv", sep=""),sep=",", header=T)
colnames(data)=a$cellline
genenames=read.table(paste(homedir,"Raw_data/gene_symbols_U133.txt", sep=""))

##Read RSI and hall probes
rsiprobes=read.table(paste(homedir,"Signatures/Probes_Eschrich.csv", sep=""), header=T, sep=",")
hallprobes=read.table(paste(homedir,"Signatures/Probes_Hall_HNSCC.csv", sep=""), header=T,sep=",")

g=match(c(rsiprobes[,1],hallprobes[,1]), rownames(data))
data1=data.frame(t(data[g,]))
data1$SF2=a$SF2
data1$site=a$tissue

##Full dataset with gene symbols
data$names=genenames$x2.Description

##Average rows with same gene name
data2=data.frame(avereps(data[,1:ncol(data)-1], ID=data$names))

##Remove empty gene names
data3=data.frame(t(data2[-which(row.names(data2)==""),]))
data3$SF2=a$SF2
data3$site=a$tissue

##Split Eschrich and Hall data

Eschdata=data1[,-(11:107)]
Halldata=data1[,-(1:10)]

write.table(Eschdata, paste(homedir, "Cleaned_data/Eschdata_NCI60.txt", sep=""))
write.table(Halldata, paste(homedir, "Cleaned_data/Halldata_NCI60.txt", sep=""))

write.table(data3, paste(homedir, "Cleaned_data/NCI60_full_gene_averaged.txt", sep=""))

