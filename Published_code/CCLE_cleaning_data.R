## ---------------------------
##
## Script name: CCLE_cleaning_data.R
##
## Purpose of script: To clean CCLE microarray expression data and match with 
## radiosensitvity metrics
##
## Author: Dr. John O'Connor
##
## Date Created: 18-08-2021
##
## Email: john.oconnor@qub.ac.uk
##
## ---------------------------
##
## Notes: 
##
## This script requires the file "CCLE_Expression_2012-09-29.res" from 
## https://sites.broadinstitute.org/ccle/datasets and 
## Supplementary Data 1 from https://www.nature.com/articles/ncomms11428 
##
##
## ---------------------------


##Loading and installing packages (if required)
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install(version = "3.13")
#BiocManager::install("limma")
library(data.table)
library(limma)

##Data I/O

#Directory for raw data (adjust this to your directory name)
homedir="C:/Users/HOMEDIRECTORY/"

x=fread(paste(homedir,"Raw_data/CCLE_Expression_2012-09-29.res", sep=""), fill=T)
datayard=read.table(paste(homedir,"Raw_data/Supplementary_Data_1_Yard_et_al.csv", sep=""), sep=",", header=T, check.names=F)

##Cleaning

##Remove A and P calls 
g=c(1,2,3,seq(5,2076,2))

##Remove 2 blank rows
x2=data.frame(x[3:54677,..g])

##Average rows with same gene name
x3=data.frame(avereps(x2[,3:ncol(x2)], ID=x2$Description))

##Remove empty gene names
x4=x3[-which(row.names(x3)==""),]

##Removing and aggregating a few problematic columns
##There is one duplicate in the original data "NCI-H292". They will be averaged here. 
dups=which(startsWith(colnames(x4),"NCIH292_LUNG"))

newcol=(x4[,dups[1]]+x4[,dups[2]])/2
x4[,dups[1]]=newcol
x4=x4[,-dups[2]]

##Remove TT_oesophagus and TT_thyroid (they don't match RS data), cause problems with matching
tts=grep("^TT_",colnames(x4))
x4=x4[,c(-tts[1],-tts[2])]

##Duplicated or missing gene names

#Transpose
x5=data.frame(t(x4))
#colnames(x5)=x4$Description

#Shorten row names and remove leading X
p=sub("\\_.*", "", rownames(x5))
p=sub("^X", "", p)
rownames(x5)=p

##Match to supplemental data from Yard
mp=match(p,toupper(datayard$`Cell Line`))
sum(is.na(mp))

##Match AUC data
x6=x5
x6$aucp=datayard$AUC[mp]
x6$site=datayard$Site[mp]
x7=subset(x6, !is.na(aucp))


##Probesets data
##Read and reduce to the two signatures which quote ProbeSet IDs

rsiprobes=read.table(paste(homedir,"Signatures/Probes_Eschrich.csv", sep=""), header=T, sep=",")
hallprobes=read.table(paste(homedir,"Signatures/Probes_Hall_HNSCC.csv", sep=""), header=T,sep=",")

g2=match(c(rsiprobes[,1],hallprobes[,1]), x2$Accession)

y1=x2[g2,]

##There is one duplicate in the original data "NCI-H292". They will be averaged here. 
newcol=(y1[,965]+y1[,993])/2
y1[,965]=newcol
y1=y1[,-993]

##Remove TT_oesophagus and TT_thyroid (they don't match RS data), cause problems with matching

y1=y1[,c(-567,-684)]

#Transpose
y2=data.frame(t(y1[,3:ncol(y1)]))

k=gsub("_", "", y1[,2], fixed=T)
k2=gsub("-", "", k, fixed=T)
k3=gsub("/", "", k2, fixed=T)
colnames(y2)=k3

#Shorten row names
p=sub("\\_.*", "", rownames(y2))
p=sub("^X", "", p)
rownames(y2)=p

##Match to supplemental data from Yard
mp=match(p,toupper(datayard$`Cell Line`))
sum(is.na(mp))

##Match AUC data
y3=y2
y3$aucp=datayard$AUC[mp]
y3$site=datayard$Site[mp]
raw_data=subset(y3, !is.na(aucp))

##Split Eschrich and Hall data

Eschdata=raw_data[,-(11:107)]
colnames(Eschdata)[1:10]=y1$Description[1:10]

Halldata=raw_data[,-(1:10)]
colnames(Halldata)[1:97]=colnames(y2)[11:107]

write.table(Eschdata, paste(homedir,"Cleaned_data/Eschdata_CCLE.txt", sep=""))
write.table(Halldata, paste(homedir,"Cleaned_data/Halldata_CCLE.txt", sep=""))

##Write Gene averaged dataset

write.table(x7, paste(homedir,"Cleaned_data/CCLE_Expression_2012-09-29_gene_averaged", sep=""))
