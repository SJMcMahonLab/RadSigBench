## ---------------------------
##
## Script name: CCLE_random_signature_testing_V1.R
##
## Purpose of script: To generate random signatures and test their accuracy 
##
## Author: Dr. John O'Connor
##
## Date Created: 19-08-2021
##
## Email: john.oconnor@qub.ac.uk
##
## ---------------------------
##
## Notes: 
##
## This script requires the cleaned data from the CCLE_cleaning_data.R script
##
##
## ---------------------------

##Data I/O
homedir="C:/Users/HOMEDIRECTORY/"
outdir=paste(homedir,"Random_signatures/CCLE/",sep="")
raw_data=read.table(paste(homedir,"Cleaned_data/CCLE_Expression_2012-09-29_gene_averaged",sep=""))

##PCR fitting function

pca.fun<-function(train,test, varthres, sigsize) {
  
  pca_model=prcomp(train[,1:sigsize],scale =TRUE)
  
  pca.var=pca_model$sdev^2
  pervar=cumsum(round(pca.var/sum(pca.var)*100,1))
  
  if ((varthres-pervar[sum(pervar<varthres)])<((pervar[sum(pervar<varthres)+1])-varthres)){
    pccut=sum(pervar<varthres)} else {pccut=sum(pervar<varthres)+1}
  
  scores=data.frame(pca_model$x)
  scores$aucp=train$aucp
  scores$metapcna=train$metapcna
  mod1=lm(as.formula(paste("aucp", "~",paste(colnames(scores)[c(1:pccut)], collapse = "+"),sep = "")),data = scores)
  
  #Rotate, scale and centre test data
  newscores=data.frame(scale(test[,1:sigsize], pca_model$center, pca_model$scale) %*% pca_model$rotation)
  newscores$metapcna=test$metapcna
  
  pca_pred <- predict(mod1, newscores, ncomp = pccut)
  output=list(pca_pred, abs(test$aucp-as.numeric(pca_pred)))
  
  return(output)
}

##Parameters
nresamp=500
sizes=c(10, 19, 31, 49, 97, 114, 127, 129, 131, 146, 168, 182, 272, 289, 306, 546)
set.seed(1001)
dt=sample(nrow(raw_data),nrow(raw_data))
varthres=80
medsigs=data.frame(matrix(NA, nrow = dim(raw_data)[1], ncol = length(sizes)))
medsigsae=data.frame(matrix(NA, nrow = dim(raw_data)[1], ncol = length(sizes)))

for (i in 1:length(sizes)){
  
  x=data.frame(matrix(NA, nrow = nresamp, ncol = sizes[i]+1))
  randpred=data.frame(matrix(NA, nrow = dim(raw_data)[1], ncol = nresamp))
  randae=data.frame(matrix(NA, nrow = dim(raw_data)[1], ncol = nresamp))
  
  for (j in 1:nresamp){
    
    x[j,1:sizes[i]]=sample(colnames(raw_data)[1:(ncol(raw_data)-2)],sizes[i])
    m1=match(x[j,1:sizes[i]], colnames(raw_data))
    
    reddata=cbind(raw_data[,m1], raw_data[,which(colnames(raw_data)=="aucp"|colnames(raw_data)=="site")])
    
    fold1=reddata[dt[1:174],]
    fold2=reddata[dt[175:349],]
    fold3=reddata[dt[350:523],]
    
    sigsize=sizes[i]
    
    ##CV1
    test=fold3
    train=rbind(fold1, fold2)
    
    f1=pca.fun(train,test, varthres, sigsize)
    
    ##CV2
    test=fold2
    train=rbind(fold1, fold3)
    
    f2=pca.fun(train,test, varthres, sigsize)
    
    ##CV3
    test=fold1
    train=rbind(fold2, fold3)
    
    f3=pca.fun(train,test, varthres, sigsize)
    
    sig_pred=rbind(as.matrix(f1[[1]]),as.matrix(f2[[1]]),as.matrix(f3[[1]]))
    sig_ae=rbind(as.matrix(f1[[2]]),as.matrix(f2[[2]]),as.matrix(f3[[2]]))
    
    x[j,sigsize+1]=mean(sig_ae)
    randpred[1:dim(randpred)[1],j]=sig_pred
    randae[1:dim(randae)[1],j]=sig_ae
    
    print(paste("SigSize=",sizes[i]))
    print(paste("Resample#=",j))
    
  }
  
  #find signature closest to the median 
  x[,length(x)+1]=(x[,length(x)]-median(x[,length(x)]))
  medind=which(abs(x[,length(x)])==min(abs(x[,length(x)])))[1]
  
  medsigs[1:dim(randpred)[1],i]=randpred[,medind]
  medsigsae[1:dim(randae)[1],i]=randae[,medind]
  
  write.table(x,paste(outdir, "sigs", sizes[i], ".txt", sep=""))
  write.table(randpred,paste(outdir, "sigs_preds", sizes[i], ".txt", sep=""))
  
}

write.table(medsigsae, paste(outdir, "medsigsae.txt", sep=""))
write.table(medsigs, paste(outdir, "medsigs.txt", sep=""))

