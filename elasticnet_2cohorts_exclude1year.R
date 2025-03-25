#!/usr/bin/env Rscript
#use all common metabolites, without adjustments in elasticnet,use 2 cohorts. one for training and the other for validation. 

library(pROC)
library(boot)
#set fold ids for cross-validation based on STUDY and pairid
#samples from a pairid should be assigned to the same fold; samples from STUD should be evenly distributed across folds
generate_foldid=function(X,Y,pairid,nfolds=5)
{
  foldid=rep(NA,nrow(X))
  uniqstudy=unique(X$STUDY)
  #assign 1..nfold for samples from each study
  for (i in 1:length(uniqstudy))
  {
    idx=which(X$STUDY==uniqstudy[i])
    uniqpair=unique(pairid[idx])
    #fold number n
    n=1
    for (j in 1:length(uniqpair))
    {
      idx1=which(pairid==uniqpair[j])
      #make sure samples from same pair in the same fold n
      foldid[idx1]=n
      n=n+1
      if (n==nfolds+1) n=1
    }
  }
  #table(foldid,X$STUDY)
  #table(foldid,Y)
  return(foldid)
}

AUCBoot1 = function(data,indices){
  boot_data = data[indices, ]
  model1=glm("CASE~BMI+diabetes+packyear",data=boot_data,family = binomial)
  predict1 <- predict(model1,newdata=boot_data, type="response")
  auc1 = as.numeric(pROC::auc(boot_data$CASE,predict1,quiet=T))
  return(c(auc1))
}

AUCBoot2 = function(data,indices){
  boot_data = data[indices, ]
  model2=glm("CASE~predictY",data=boot_data,family = binomial)
  predict2 <- predict(model2,newdata=boot_data, type="response")
  auc2 = as.numeric(pROC::auc(boot_data$CASE,predict2,quiet=T))
  return(c(auc2))
}

AUCBoot3 = function(data,indices){
  boot_data = data[indices, ]
  model3=glm("CASE~BMI+diabetes+packyear+predictY",data=boot_data,family = binomial)
  predict3 <- predict(model3,newdata=boot_data, type="response")
  auc3 = as.numeric(pROC::auc(boot_data$CASE,predict3,quiet=T))
  return(c(auc3))
}

library(glmnet)
###need to specify your input file for dataset
alldat=read.csv("../data/Dataset for LASSO/Dataset for LASSO-7Jun23.csv")
#change race 
alldat$RACE[which(alldat$RACE!=1)]="Other"
alldat$diabetes[which(alldat$diabetes=="M")]=NA
alldat$diabetes=factor(alldat$diabetes)

#make some variables to be factors
factors=c("COHORT","STUDY","BATCH_NUM","CELL_SET","GENDER","RACE")
for (i in 1:length(factors)) alldat[,factors[i]]=as.factor(alldat[,factors[i]])
# ###need to specify your input file for FDR significant metabolites
FDRmetabolite=read.csv("../data/Dataset for LASSO/FDR significant metabolites.csv")
FDRmetabolite=FDRmetabolite$Metabolite
tmp=alldat[,FDRmetabolite]
idx=complete.cases(t(tmp))
#common 41 markers
FDRmetabolite=FDRmetabolite[idx]

#all 539 common markers
allmetabolite=alldat[,20:ncol(alldat)]
idx=complete.cases(t(allmetabolite))
allmarkers=colnames(allmetabolite)[idx]
unknownmarkers=read.csv("../data/Unknown metabolite list.csv")
table(allmarkers %in% unknownmarkers$x)
# FALSE  TRUE 
# 472    67
allmarkers=allmarkers[!allmarkers %in% unknownmarkers$x]
length(allmarkers)
#[1] 472
#2 cohorts
alldat_followup=read.csv("../data/Dataset for LASSO-22Jul24.csv")
clientid=alldat_followup$CLIENT_IDENTIFIER[which(alldat_followup$follow_time>=1)]
alldat1=alldat[alldat$CLIENT_IDENTIFIER %in% clientid,]
ATBCdat=alldat1[alldat1$COHORT=="ATBC",]
PLCOdat=alldat1[alldat1$COHORT=="PLCO",]

#compute 3 AUCs (riskfactors,markers, and combined)
get3AUC=function(testdat,predictY=predictY_min)
{
  testdat$predictY=predictY
  model1=glm("CASE~BMI+diabetes+packyear",data=testdat,family = binomial)
  predict1 <- predict(model1,newdata=testdat, type="response")
  auc1 = as.numeric(pROC::auc(testdat$CASE,predict1,quiet=T))
  boot_auc = boot(data =testdat, statistic = AUCBoot1, R = 1000)
  tmp=boot.ci(boot_auc,type="perc")
  auc1low=tmp$percent[4] #low
  auc1high=tmp$percent[5] #high
  auc1=paste0(round(auc1,2),"(",round(auc1low,2),",",round(auc1high,2),")")
  
  model2=glm("CASE~predictY",data=testdat,family = binomial)
  predict2 <- predict(model2,newdata=testdat, type="response")
  auc2 = as.numeric(pROC::auc(testdat$CASE,predict2,quiet=T))
  boot_auc = boot(data =testdat, statistic = AUCBoot2, R = 1000)
  tmp=boot.ci(boot_auc,type="perc")
  auc2low=tmp$percent[4] #low
  auc2high=tmp$percent[5] #high
  auc2=paste0(round(auc2,2),"(",round(auc2low,2),",",round(auc2high,2),")")

  model3=glm("CASE~BMI+diabetes+packyear+predictY",data=testdat,family = binomial)
  predict3 <- predict(model3,newdata=testdat, type="response")
  auc3 = as.numeric(pROC::auc(testdat$CASE,predict3,quiet=T))
  boot_auc = boot(data =testdat, statistic = AUCBoot3, R = 1000)
  tmp=boot.ci(boot_auc,type="perc")
  auc3low=tmp$percent[4] #low
  auc3high=tmp$percent[5] #high
  auc3=paste0(round(auc3,2),"(",round(auc3low,2),",",round(auc3high,2),")")
  res=list(auc1=auc1,auc2=auc2,auc3=auc3)
  return(res)
}

# #version 2: add age into models
# get3AUC=function(testdat,predictY=predictY_min)
# {
#   testdat$predictY=predictY
#   model1=glm("CASE~BMI+diabetes+packyear+AGE",data=testdat,family = binomial)
#   predict1 <- predict(model1,newdata=testdat, type="response")
#   auc1 = as.numeric(pROC::auc(testdat$CASE,predict1,quiet=T))
#   model2=glm("CASE~predictY",data=testdat,family = binomial)
#   predict2 <- predict(model2,newdata=testdat, type="response")
#   auc2 = as.numeric(pROC::auc(testdat$CASE,predict2,quiet=T))
#   model3=glm("CASE~BMI+diabetes+packyear+AGE+predictY",data=testdat,family = binomial)
#   predict3 <- predict(model3,newdata=testdat, type="response")
#   auc3 = as.numeric(pROC::auc(testdat$CASE,predict3,quiet=T))
#   res=list(auc1=auc1,auc2=auc2,auc3=auc3)
#   return(res)
# }
#select markers using training data and compute AUC on validation data
computeAUC=function(traindat=ATBCdat,testdat=PLCOdat,selmarkers=allmarkers,prefix="ATBC")
{
  if (!is.na(traindat$ID[1]))
  {
    traindat$allID=traindat$ID
  }else
  {
    traindat$allID=traindat$LIID
  }
  if (!is.na(testdat$ID[1]))
  {
    testdat$allID=testdat$ID
  }else
  {
    testdat$allID=testdat$LIID
  }
  #elastic net
  trainX=traindat[,selmarkers]
  trainY=traindat$CASE
  testX=testdat[,selmarkers]
  testY=testdat$CASE
  pairid=traindat$CELL_SET
  foldid=generate_foldid(X=traindat,Y=trainY,pairid,nfolds=10)
  set.seed(1000)
  cvfit=cv.glmnet(data.matrix(trainX),trainY,foldid = foldid, alpha=0.5,family=quasibinomial)
  coeffs_min=coef(cvfit, s=cvfit$lambda.min) 
  coeffs_min=data.frame(name=coeffs_min@Dimnames[[1]][coeffs_min@i + 1], coefficient=coeffs_min@x) 
  # coeffs_1se=coef(cvfit, s=cvfit$lambda.1se) 
  # coeffs_1se=data.frame(name=coeffs_1se@Dimnames[[1]][coeffs_1se@i + 1], coefficient=coeffs_1se@x) 
  # 
  #apply on validation
  #lasso
  predictY_min=predict(cvfit,newx=data.matrix(testX),s=cvfit$lambda.min)
  predictY_train_min=predict(cvfit,newx=data.matrix(trainX),s=cvfit$lambda.min)
  # predictY_1se=predict(cvfit,newx=data.matrix(testX),s=cvfit$lambda.1se)
  model1=glm("CASE~BMI+diabetes+packyear",data=testdat,family = binomial)
  predictrisk <- predict(model1,newdata=testdat, type="response")
  testdat$predictY=predictY_min
  model3=glm("CASE~BMI+diabetes+packyear+predictY",data=testdat,family = binomial)
  predictcombined <- predict(model3,newdata=testdat, type="response")
  #riskfactor vs combined
  auc.pvalue1=roc.test(response=testdat$CASE,predictor1 = predictrisk,predictor2 = predictcombined,method="bootstrap",alternative="less")$p.value
  print(paste0("AUC: riskfactors vs combined p:",auc.pvalue1))
  #marker vs combined
  auc.pvalue2=roc.test(response=testdat$CASE,predictor1 = testdat$predictY[,1],predictor2 = predictcombined,method="bootstrap",alternative="less")$p.value
  print(paste0("AUC: marker vs combined p:",auc.pvalue2))
  #Add comparison between random of AUC 0.5
  model4=glm("CASE~1",data=testdat,family = binomial)
  predictrand <- predict(model4,newdata=testdat, type="response")
  testdat$predictrand=predictrand
  aucrand=as.numeric(pROC::auc(testdat$CASE,predictrand,quiet=T))
  auc.pvalue3=roc.test(response=testdat$CASE,predictor1 = testdat$predictrand,predictor2 = testdat$predictY[,1],method="bootstrap",alternative="less")$p.value
  print(paste0("AUC: marker vs intercept (AUC 0.5) p:",auc.pvalue3))
  auc.pvalue4=roc.test(response=testdat$CASE,predictor1 = testdat$predictrand,predictor2 = predictrisk,method="bootstrap",alternative="less")$p.value
  print(paste0("AUC: riskfactros vs intercept (AUC 0.5) p:",auc.pvalue4))
  auc.pvalue5=roc.test(response=testdat$CASE,predictor1 = testdat$predictrand,predictor2 = predictcombined,method="bootstrap",alternative="less")$p.value
  print(paste0("AUC: combined vs intercept (AUC 0.5) p:",auc.pvalue5))
  #refitting linear models
  selmarkers1=coeffs_min$name[!grepl("Intercept",coeffs_min$name)]
  fm=c("CASE ~",paste0(selmarkers1,collapse = "+"))
  model2=glm(formula(paste(fm, collapse = " ")),data=traindat,family = binomial)
  coeffs_refit=summary(model2)$coefficients[,1]
  coeffs_refit=data.frame(name=names(coeffs_refit),coefficient=unname(coeffs_refit))
  predictY_min1=predict(model2,newdata=testdat,type="response")
  print(cor(predictY_min,predictY_min1)) #0.84
  
  # selmarkers1=coeffs_1se$name[!grepl("Intercept",coeffs_1se$name)]
  # fm=c("CASE ~",paste0(selmarkers1,collapse = "+"))
  # model2=glm(formula(paste(fm, collapse = " ")),data=traindat,family = binomial)
  # predictY_1se1=predict(model2,newdata=testdat,type="response")
  #print(cor(predictY_1se,predictY_1se1)) #0.92
  #only use elesticnet results
  auc=data.frame(matrix(0,nrow=2,ncol=3))
  rownames(auc)=c("elsticnet","refit")
  colnames(auc)=c("risk","maker","combine")
  auc[1,]=get3AUC(testdat=testdat,predictY = predictY_min)
  #auc[2,]=get3AUC(testdat=testdat,predictY = predictY_1se)
  auc[2,]=get3AUC(testdat=testdat,predictY = predictY_min1)
  #auc[4,]=get3AUC(testdat=testdat,predictY = predictY_1se1)
  trainres=data.frame(ID=traindat$allID,predictY=predictY_train_min)
  testres=data.frame(ID=testdat$allID,predictY=predictY_min)
  resall=list(auc=auc,coeffs_min=coeffs_min,coeffs_refit=coeffs_refit,
              predictY_min=predictY_min,predictY_train_min=predictY_train_min,
              trainres=trainres,testres=testres,predictrisk=predictrisk,predictcombined=predictcombined)
  #to save weights
  weightfile=paste0("../result/",prefix,"_elsticnetweights_rmunknown.csv")
  tmp=coeffs_min
  colnames(tmp)[1]="Marker"
  tmp=tmp[-1,]
  write.csv(tmp,file=weightfile)
  return(resall)
}

#Use ATBC as training
trainATBC=computeAUC(prefix = "ATBC_exclue1year")
#AUC: riskfactors vs combined p:0.00162269266296895
#AUC: marker vs combined p:0.122206044738911
#AUC: marker vs intercept (AUC 0.5) p:7.54647801909694e-10
#AUC: riskfactros vs intercept (AUC 0.5) p:2.71962590457493e-05
#AUC: combined vs intercept (AUC 0.5) p:2.06655375540113e-12
trainATBC$auc
# risk              maker            combine
# elsticnet 0.59(0.56,0.64) 0.63(0.59,0.67) 0.65(0.61,0.69)
# refit     0.59(0.56,0.64)  0.64(0.6,0.68)  0.66(0.62,0.7)
#Use PLCO as training
trainPLCO=computeAUC(traindat=PLCOdat,testdat=ATBCdat,prefix="PLCO_exclude1year")
#AUC: AUC: riskfactors vs combined p:0.0404032222064678
#AUC: marker vs combined p:0.239963320934529
#AUC: marker vs intercept (AUC 0.5) p:4.40848357369694e-06
#AUC: riskfactros vs intercept (AUC 0.5) p:0.00157548500391388
#AUC: combined vs intercept (AUC 0.5) p:1.07179924209814e-06
trainPLCO$auc
# risk              maker           combine
# elsticnet 0.57(0.52,0.61) 0.59(0.55,0.63)  0.6(0.57,0.65)
# refit     0.57(0.52,0.61) 0.58(0.54,0.62) 0.59(0.55,0.63)
#The results
#the weights from elastic net are stored in coeffs_min
save(trainATBC,trainPLCO,PLCOdat,ATBCdat,file="../result/elasticnet_exclude1year_2cohorts.RData")
head(trainATBC$coeffs_min)
#the weights from refit-logistic regression are stored in coeffs_min
head(trainATBC$coeffs_refit)
#the AUCs are stored in auc
head(trainATBC$auc)
head(trainPLCO$coeffs_min)
length(intersect(trainATBC$coeffs_min$name,trainPLCO$coeffs_min$name)) #4-1

#to draw ROC
library(ROCR) #prediction
library(pROC) #roc,auc,plot.roc
calauc=function(pfit=NULL,y=NULL)
{
  rocdata=function(pfit,y)
  {
    yy=as.numeric(as.factor(y[!is.na(pfit)]))
    pfit=pfit[!is.na(pfit)]
    xx <- cbind(pfit,yy)
    if (max(xx[,2])>=2)
      xx[,2] <- ifelse(xx[,2]>=2,1,0)
    xx <- xx[order(xx[,1]),]
    pred <- prediction(xx[,1],xx[,2]) #ROCR
    auc.tmp <- performance(pred,"auc")
    if (auc.tmp@y.values<0.5)
    {
      xx[,1]=1-xx[,1]
      pred <- prediction(xx[,1],xx[,2])
      auc.tmp <- performance(pred,"auc")
    }
    auc <- round(as.numeric(auc.tmp@y.values),2)
    result=list(pred=pred,auc=auc)
  }
  res_train=rocdata(pfit,y)
  res=as.numeric(res_train$auc)
}

#precit1: estimates based on risk factors
#predict2: estimates based on markers
#predict3: combined
#aucres=NULL: legend without CI; aucres (AUC with CI)
#cex: font size
plotroc3=function(predict1=trainATBC$predictrisk,predict2=trainATBC$predictY_min[,1],predict3=trainATBC$predictcombined,y=PLCOdat$CASE,main="",
                  prefix="trainATBC",aucres=NULL,cex=1.3)
{
  colors=c("#4DBBD5FF","#3C5488FF","#E64B35FF")
  predict1_=predict1[order(predict1)]
  TP1=FP1=rep(0,length(y))
  for (i in 1:length(y)) {
    TP1[i] <- mean(predict1[y>0]>=predict1_[i],na.rm=T) 
    FP1[i] <- mean(predict1[y==0]>=predict1_[i],na.rm=T) 
  }
  
  pdf(paste0("../result/",prefix,"_AUC.pdf"),width = 8,height = 8)
  par(mar=c(6,5,2,1))
  plot(c(FP1,0),c(TP1,0),type="l",xlim=c(0,1),ylim=c(0,1),lwd=3,xlab="1-specificity",ylab="Sensitivity",main=main,cex.lab=cex,cex.axis=cex,col=colors[1],yaxs="i")
  
  #abline(0,1,lty=2)
  auc1=calauc(predict1,y)
  #auc1= as.numeric(pROC::auc(y,predict1,quiet=T))
  print(paste0("auc1=",auc1))
  #pauc1=pauc(marker = predict1,status = y,fpr = 0.05)
  #print(paste0("pauc1=",round(pauc1,4)))
  #text(0.07,0.94,paste0("AUC=",auc1),cex=1.2)
  
  predict2_=predict2[order(predict2)]
  TP2=FP2=rep(0,length(y))
  for (i in 1:length(y)) {
    TP2[i] <- mean(predict2[y>0]>=predict2_[i]) 
    FP2[i] <- mean(predict2[y==0]>=predict2_[i]) 
  }
  lines(c(FP2,0),c(TP2,0),lwd=2,col=colors[2])
  #lines(c(0.05,0.05),c(0,1),lty=2)
  auc2=calauc(predict2,y)
  print(paste0("auc2=",auc2))
  #pauc2=pauc(marker = predict2,status = y,fpr = 0.05)
  #print(paste0("pauc2=",round(pauc2,4)))
  
  predict3_=predict3[order(predict3)]
  TP3=FP3=rep(0,length(y))
  for (i in 1:length(y)) {
    TP3[i] <- mean(predict3[y>0]>=predict3_[i],na.rm=T) 
    FP3[i] <- mean(predict3[y==0]>=predict3_[i],na.rm=T) 
  }
  lines(c(FP3,0),c(TP3,0),lwd=2,col=colors[3])
  auc3=calauc(predict3,y)
  print(paste0("auc3=",auc3))
  
  if (is.null(aucres))
  {
    legend=c(paste0("BMI+Diabetes+Packyear:AUC=",round(auc1,3)),paste0("Metabolite markers:AUC=",round(auc2,3)),paste0("Combined:AUC=",round(auc3,3)))
  }else
  {
    legend=c(paste0("BMI+Diabetes+Packyear:AUC=",aucres$auc[1,1]),paste0("Metabolite markers:AUC=",aucres$auc[1,2]),paste0("Combined:AUC=",aucres$auc[1,3]))
  }
  legend("bottomright",legend=legend,col=colors,lty=1,cex=cex)
  dev.off()
  #res=data.frame(auc1=auc1,pauc1=pauc1,auc2=auc2,pauc2=pauc2,auc3=auc3,pauc3=pauc3,auc.pvalue=auc.pvalue,pauc.pvalue=pauc.pvalue,stringsAsFactors = F)
  #return(res)
}

plotroc3()
plotroc3(predict1=trainPLCO$predictrisk,predict2=trainPLCO$predictY_min[,1],predict3=trainPLCO$predictcombined,y=ATBCdat$CASE,main="",
                  prefix="trainPLCO")
plotroc3(aucres = trainATBC,prefix="trainATBC1")
plotroc3(predict1=trainPLCO$predictrisk,predict2=trainPLCO$predictY_min[,1],predict3=trainPLCO$predictcombined,y=ATBCdat$CASE,main="",
         prefix="trainPLCO1",aucres = trainPLCO)
#include CI in AUC
plotroc3(aucres = trainATBC,prefix="trainATBC2",cex=1.6)
plotroc3(predict1=trainPLCO$predictrisk,predict2=trainPLCO$predictY_min[,1],predict3=trainPLCO$predictcombined,y=ATBCdat$CASE,main="",
         prefix="trainPLCO2",aucres = trainPLCO,cex=1.6)
#not include CI
plotroc3(aucres = NULL,prefix="trainATBC2",cex=1.6)
plotroc3(predict1=trainPLCO$predictrisk,predict2=trainPLCO$predictY_min[,1],predict3=trainPLCO$predictcombined,y=ATBCdat$CASE,main="",
         prefix="trainPLCO2",aucres = NULL,cex=1.6)
