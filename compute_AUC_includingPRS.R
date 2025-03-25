#!/usr/bin/env Rscript
library(data.table)
library(tidyverse)
plink="/usr/local/apps/plink/1.9.0-beta4.4/plink"
plink2="/usr/local/apps/plink/2.3-alpha/plink2"
IDtable=read.csv("../data/063_Rachel_09082022/IDs for Dataset for LASSO-8Nov24.csv")
table(IDtable$COHORT)
# ATBC PLCO 
# 744 1222
table(is.na(IDtable$ID),is.na(IDtable$LIID))
IDtable$allID=IDtable$ID
idx=is.na(IDtable$allID)
IDtable$allID[idx]=IDtable$LIID[idx]
table(duplicated(IDtable$allID))
# FALSE  TRUE 
# 1464   502
length(unique(IDtable$allID)) #1464
idx=duplicated(IDtable$allID)
table(IDtable$COHORT[!idx])
# ATBC PLCO 
# 744  720
IDtable=IDtable %>% mutate(atbcid=str_split(Masked.ID.associated.with.genetic.data,"_",simplify = T)[,1])
IDtable$atbcid[which(IDtable$atbcid=="")]=NA
sum(!is.na(IDtable$atbcid)) #584

pheno=read.csv("../data/Dataset for LASSO/Dataset for LASSO-7Jun23.csv")
pheno$diabetes[which(pheno$diabetes=="M")]=NA
pheno$diabetes=factor(pheno$diabetes)
table(is.na(pheno$ID),is.na(pheno$LIID))
pheno$allID=pheno$ID
idx=is.na(pheno$allID)
pheno$allID[idx]=pheno$LIID[idx]
table(IDtable$allID %in% pheno$allID )
table(pheno$COHORT)
# ATBC PLCO 
# 744  720

check_plcosamples=function()
{
  allplcodat=NULL
  allpsamfiles=list.files("../data/063_Rachel_09082022/PLCO/",".psam")
  for(i in 1:length(allpsamfiles))
  {
    plcopsam=as.data.frame(fread(paste0('../data/063_Rachel_09082022/PLCO/',allpsamfiles[i])))
    colnames(plcopsam)[1]="IID"
    plcopsam=plcopsam %>% mutate(plcoid=str_split(IID,"_",simplify = T)[,1])
    if (sum(plcopsam$plcoid %in% IDtable$plco_id)>0)
    {
      tmp=intersect(plcopsam$plcoid,IDtable$plco_id)
      idx=match(tmp,plcopsam$plcoid)
      tmp=data.frame(plcoid=tmp,IID=plcopsam$IID[idx],file=allpsamfiles[i])
      allplcodat=rbind(allplcodat,tmp)
      print(allpsamfiles[i])
      print(sum(plcopsam$plcoid %in% IDtable$plco_id))
    }
  }
  write.table(allplcodat$IID,file="../result/PGS002740_pancreatic_plcoIID.txt",row.names = F,col.names = F,quote=F)
}
PGS002740score=as.data.frame(fread("../data/063_Rachel_09082022/PGS002740_hmPOS_GRCh38.txt"))
scoredat=NULL
for(i in 1:nrow(PGS002740score))
{
  a2=unlist(strsplit(PGS002740score$hm_inferOtherAllele[i],"/"))
  tmp1=data.frame(SNP=paste0("chr",PGS002740score$hm_chr[i],":",PGS002740score$hm_pos[i],":",PGS002740score$effect_allele[i],":",a2),A1=PGS002740score$effect_allele[i],beta=PGS002740score$effect_weight[i])
  tmp2=data.frame(SNP=paste0("chr",PGS002740score$hm_chr[i],":",PGS002740score$hm_pos[i],":",a2,":",PGS002740score$effect_allele[i]),A1=PGS002740score$effect_allele[i],beta=PGS002740score$effect_weight[i])
  scoredat=rbind(scoredat,tmp1,tmp2)
}
write.table(scoredat,file="../result/PGS002740score.txt",row.names = F,sep="\t",quote=F)

compute_plcoprs=function()
{
  allplcodat$nsnps=NA
  allplcodat$PRS=NA
  uniqfiles=unique(allplcodat$file)
  for(i in 1:length(uniqfiles))
  {
    prefix=paste0("../data/063_Rachel_09082022/PLCO/",uniqfiles[i])
    prefix=gsub(".psam","",prefix,fixed = T)
    cmd=paste0(plink2," --pfile ",prefix," --keep ../result/PGS002740_pancreatic_plcoIID.txt "," --allow-no-sex --score ../result/PGS002740score.txt "," cols=+scoresums,-scoreavgs header no-mean-imputation "," --out ../result/tmp")
    system(cmd)
    tmp=as.data.frame(fread("../result/tmp.sscore"))
    system("rm ../result/tmp.sscore")
    idx=match(tmp$`#IID`,allplcodat$IID)
    if (sum(is.na(idx))>0) stop()
    allplcodat$nsnps[idx]=tmp$ALLELE_CT/2
    allplcodat$PRS[idx]=tmp$SCORE1_SUM
  }
}

check_atbcsamples=function()
{
  allatbcdat=NULL
  allpsamfiles=list.files("../data/063_Rachel_09082022/ATBC/",".psam")
  for(i in 1:length(allpsamfiles))
  {
    atbcpsam=as.data.frame(fread(paste0('../data/063_Rachel_09082022/ATBC/',allpsamfiles[i])))
    colnames(atbcpsam)[1]="IID"
    atbcpsam=atbcpsam %>% mutate(atbcid=str_split(IID,"_",simplify = T)[,1])
    if (sum(atbcpsam$atbcid %in% IDtable$atbcid)>0)
    {
      tmp=intersect(atbcpsam$atbcid,IDtable$atbcid)
      idx=match(tmp,atbcpsam$atbcid)
      tmp=data.frame(atbcid=tmp,IID=atbcpsam$IID[idx],file=allpsamfiles[i])
      allatbcdat=rbind(allatbcdat,tmp)
      print(allpsamfiles[i])
      print(sum(atbcpsam$atbcid %in% IDtable$atbcid))
    }
  }
  write.table(allatbcdat$IID,file="../result/PGS002740_pancreatic_atbcIID.txt",row.names = F,col.names = F,quote=F)
}

compute_atbcprs=function()
{
  allatbcdat$nsnps=NA
  allatbcdat$PRS=NA
  uniqfiles=unique(allatbcdat$file)
  for(i in 1:length(uniqfiles))
  {
    prefix=paste0("../data/063_Rachel_09082022/ATBC/",uniqfiles[i])
    prefix=gsub(".psam","",prefix,fixed = T)
    cmd=paste0(plink2," --pfile ",prefix," --keep ../result/PGS002740_pancreatic_atbcIID.txt "," --allow-no-sex --score ../result/PGS002740score.txt "," cols=+scoresums,-scoreavgs header no-mean-imputation "," --out ../result/tmp")
    system(cmd)
    tmp=as.data.frame(fread("../result/tmp.sscore"))
    system("rm ../result/tmp.sscore")
    idx=match(tmp$`#IID`,allatbcdat$IID)
    if (sum(is.na(idx))>0) stop()
    allatbcdat$nsnps[idx]=tmp$ALLELE_CT/2
    allatbcdat$PRS[idx]=tmp$SCORE1_SUM
  }
}
allprsdat=allplcodat
allprsdat$cohort="PLCO"
colnames(allprsdat)[1]="ID"
tmp=allatbcdat
tmp$cohort="ATBC"
colnames(tmp)[1]="ID"
allprsdat=rbind(allprsdat,tmp)
table(allprsdat$cohort)
# ATBC PLCO 
# 552  636
allprsdat$allID=NA
tmp=intersect(allprsdat$ID,IDtable$plco_id)
idx1=match(tmp,allprsdat$ID)
idx2=match(tmp,IDtable$plco_id)
allprsdat$allID[idx1]=IDtable$allID[idx2]
tmp=intersect(allprsdat$ID,IDtable$atbcid)
idx1=match(tmp,allprsdat$ID)
idx2=match(tmp,IDtable$atbcid)
allprsdat$allID[idx1]=IDtable$allID[idx2]
all(allprsdat$allID %in% pheno$allID)
idx=match(allprsdat$allID,pheno$allID)
allprsdat$CASE=pheno$CASE[idx]
table(allprsdat$cohort,allprsdat$CASE)
#        0   1
# ATBC 330 222
# PLCO 341 295

alldat=read.csv("../data/Dataset for LASSO-22Jul24.csv")
y5_1 <- alldat[alldat$SAMPLE_TYPE==1&alldat$follow_time<=5,]  #358 samples
y5_2 <- alldat[alldat$SAMPLE_TYPE==2&alldat$follow_time<=5,]  #270 samples
y5_2 <- y5_2[!y5_2$LIID %in% y5_1$LIID,]                                           #166 samples
y5 <- rbind(y5_1, y5_2)
y5$allID=y5$ID
idx=which(is.na(y5$allID))
y5$allID[idx]=y5$LIID[idx]
table(y5$allID %in% allprsdat$allID)
# FALSE  TRUE 
# 121   403
allprsdat$y5=NA
#case
allprsdat$y5[allprsdat$allID %in% y5$allID[which(y5$CASE==1)]]=1
allprsdat$y5[allprsdat$allID %in% y5$allID[which(y5$CASE==0)]]=0
table(allprsdat$cohort,allprsdat$y5)
#          0   1
# # ATBC  59   3
# # PLCO 185 156
allprsdat$y5_SAMPLE_TYPE=NA
allprsdat$y5_SAMPLE_TYPE[allprsdat$allID %in% y5$allID[which(y5$SAMPLE_TYPE==1)]]=1
allprsdat$y5_SAMPLE_TYPE[allprsdat$allID %in% y5$allID[which(y5$SAMPLE_TYPE==2)]]=2
table(allprsdat$cohort,allprsdat$y5_SAMPLE_TYPE)
#       1   2
# ATBC  62   0
# PLCO 176 165
table(y5$COHORT,y5$SAMPLE_TYPE)
#        1   2
# ATBC 138   0
# PLCO 220 166
idx=match(y5$allID,pheno$allID)
table(pheno$COHORT[idx],pheno$CASE[idx])
#        0   1
# ATBC  69  69
# PLCO 193 193
idx=grepl("Eur",allprsdat$file,ignore.case = T) | allprsdat$cohort=="ATBC"
table(idx)
# FALSE  TRUE 
# 58  1130
allprsdat$EUR[!idx]=0
write.csv(allprsdat,file="../result/PGS002740prs.csv",row.names = F,quote=F)
allprsdat=read.csv("../result/PGS002740prs.csv")

quantile(allprsdat$PRS[allprsdat$cohort=="PLCO"])
# 0%       25%       50%       75%      100% 
# -3.266160 -2.227182 -1.873900 -1.375938  0.420092 
quantile(allprsdat$PRS[allprsdat$cohort=="ATBC"])
# 0%       25%       50%       75%      100% 
# -3.170100 -2.267583 -1.997955 -1.710317 -0.771578
idx=match(allprsdat$allID,pheno$allID)
pheno$EUR=NA
pheno$EUR[idx]=allprsdat$EUR
pheno$prs=NA
pheno$prs[idx]=allprsdat$PRS
table(pheno$COHORT,is.na(pheno$prs))
boxplot(pheno$prs~pheno$CASE)
# FALSE TRUE
# ATBC   552  192
# PLCO   636   84
t.test(pheno$prs[pheno$CASE==1],pheno$prs[pheno$CASE==0])
#t = 5.2518, df = 1081.2, p-value = 1.813e-07
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   0.1096984 0.2405631
# sample estimates:
#   mean of x mean of y 
# -1.773415 -1.948545
idx=which(is.na(pheno$prs))
pheno$naprs=0
pheno$naprs[idx]=1
pheno$naprs=as.factor(pheno$naprs)
pheno$prs1=pheno$prs
pheno$prs1[idx]=mean(pheno$prs,na.rm=T)
idx=pheno$COHORT=="ATBC"
table(pheno$CASE[idx],pheno$naprs[idx])
#    0   1
# 0 330  42
# 1 222 150
table(pheno$CASE[!idx],pheno$naprs[!idx])
#    0   1
# 0 341  19
# 1 295  65

#opt: EUR or all
#opt1: compelte or all
library(pROC)
library(boot)
AUCBoot0 = function(data,indices){
  boot_data = data[indices, ]
  model1=glm("CASE~BMI+diabetes+packyear",data=boot_data,family = binomial)
  predict1 <- predict(model1,newdata=boot_data, type="response")
  auc1 = as.numeric(pROC::auc(boot_data$CASE,predict1,quiet=T))
  return(c(auc1))
}
AUCBoot1 = function(data,indices){
  boot_data = data[indices, ]
  model1=glm("CASE~metabscore",data=boot_data,family = binomial)
  predict1 <- predict(model1,newdata=boot_data, type="response")
  auc1 = as.numeric(pROC::auc(boot_data$CASE,predict1,quiet=T))
  return(c(auc1))
}
AUCBoot2 = function(data,indices){
  boot_data = data[indices, ]
  model1=glm("CASE~prs",data=boot_data,family = binomial)
  predict1 <- predict(model1,newdata=boot_data, type="response")
  auc1 = as.numeric(pROC::auc(boot_data$CASE,predict1,quiet=T))
  return(c(auc1))
}
AUCBoot3 = function(data,indices){
  boot_data = data[indices, ]
  model1=glm("CASE~metabscore+prs",data=boot_data,family = binomial)
  predict1 <- predict(model1,newdata=boot_data, type="response")
  auc1 = as.numeric(pROC::auc(boot_data$CASE,predict1,quiet=T))
  return(c(auc1))
}
AUCBoot5 = function(data,indices){
  boot_data = data[indices, ]
  model1=glm("CASE~metabscore+prs+BMI+diabetes+packyear",data=boot_data,family = binomial)
  predict1 <- predict(model1,newdata=boot_data, type="response")
  auc1 = as.numeric(pROC::auc(boot_data$CASE,predict1,quiet=T))
  return(c(auc1))
}

compute_auccohort=function(traincohort="ATBC",opt="all",opt1="all")
{
  if (traincohort=="ATBC")
  {
    trainmetabres=trainATBCres
    testmetabres=testATBCres
  }else
  {
    trainmetabres=trainPLCOres
    testmetabres=testPLCOres
  }
  if (opt=="EUR")
  {
    trainmetabres=trainmetabres[!trainmetabres$ID %in% pheno$allID[which(pheno$EUR==0)],]
    testmetabres=testmetabres[!testmetabres$ID %in% pheno$allID[which(pheno$EUR==0)],]
  }
  idx=match(trainmetabres$ID,pheno$allID)
  traindat=pheno[idx,]
  traindat$metabscore=trainmetabres[,2]
  idx=match(testmetabres$ID,pheno$allID)
  testdat=pheno[idx,]
  testdat$metabscore=testmetabres[,2]
  if (opt1=="complete")
  {
    traindat=traindat[!is.na(traindat$prs),]
    testdat=testdat[!is.na(testdat$prs),]
  }
  
  idx=is.na(traindat$prs)
  model0=glm("CASE~BMI+diabetes+packyear",data=traindat[!idx,],family = binomial)
  predict0=predict(model0,newdata=testdat,type="response")
  auc0 = round(as.numeric(pROC::auc(testdat$CASE,predict0,quiet=T)),2)
  boot_auc = boot(data =testdat, statistic = AUCBoot0, R = 1000)
  tmp=boot.ci(boot_auc,type="perc")
  auc0low=tmp$percent[4] #low
  auc0high=tmp$percent[5] #high
  auc0=paste0(round(auc0,2),"(",round(auc0low,2),",",round(auc0high,2),")")
  
  model1=glm("CASE~metabscore",data=traindat[!idx,],family = binomial)
  predict1=predict(model1,newdata=testdat,type="response")
  auc1 = as.numeric(pROC::auc(testdat$CASE,predict1,quiet=T))
  boot_auc = boot(data =testdat, statistic = AUCBoot1, R = 1000)
  tmp=boot.ci(boot_auc,type="perc")
  auc1low=tmp$percent[4] #low
  auc1high=tmp$percent[5] #high
  auc1=paste0(round(auc1,2),"(",round(auc1low,2),",",round(auc1high,2),")")
  
  model2=glm("CASE~prs",data=traindat,family = binomial)
  predict2=predict(model2,newdata=testdat,type="response")
  auc2 = as.numeric(pROC::auc(testdat$CASE,predict2,quiet=T))
  boot_auc = boot(data =testdat, statistic = AUCBoot2, R = 1000)
  tmp=boot.ci(boot_auc,type="perc")
  auc2low=tmp$percent[4] #low
  auc2high=tmp$percent[5] #high
  auc2=paste0(round(auc2,2),"(",round(auc2low,2),",",round(auc2high,2),")")
  
  model3=glm("CASE~metabscore+prs",data=traindat,family = binomial)
  predict3=predict(model3,newdata=testdat,type="response")
  idx1=testdat$naprs==1
  predict3[idx1]=predict1[idx1]
  auc3 = as.numeric(pROC::auc(testdat$CASE,predict3,quiet=T))
  boot_auc = boot(data =testdat, statistic = AUCBoot3, R = 1000)
  tmp=boot.ci(boot_auc,type="perc")
  auc3low=tmp$percent[4] #low
  auc3high=tmp$percent[5] #high
  auc3=paste0(round(auc3,2),"(",round(auc3low,2),",",round(auc3high,2),")")
  
  
  model4=glm("CASE~metabscore+BMI+diabetes+packyear",data=traindat[!idx,],family = binomial)
  predict4=predict(model4,newdata=testdat,type="response")
  model5=glm("CASE~metabscore+prs+BMI+diabetes+packyear",data=traindat,family = binomial)
  predict5=predict(model5,newdata=testdat,type="response")
  idx1=is.na(predict5)
  predict5[idx1]=predict4[idx1]
  auc5 = as.numeric(pROC::auc(testdat$CASE,predict5,quiet=T))
  boot_auc = boot(data =testdat, statistic = AUCBoot5, R = 1000)
  tmp=boot.ci(boot_auc,type="perc")
  auc5low=tmp$percent[4] #low
  auc5high=tmp$percent[5] #high
  auc5=paste0(round(auc5,2),"(",round(auc5low,2),",",round(auc5high,2),")")
  auc.pvalue0=roc.test(response=testdat$CASE,predictor1 = predict0,predictor2 = predict5,method="bootstrap",alternative="less")$p.value
  auc.pvalue1=roc.test(response=testdat$CASE,predictor1 = predict1,predictor2 = predict5,method="bootstrap",alternative="less")$p.value
  auc.pvalue2=roc.test(response=testdat$CASE,predictor1 = predict2,predictor2 = predict5,method="bootstrap",alternative="less")$p.value
  auc.pvalue3=roc.test(response=testdat$CASE,predictor1 = predict3,predictor2 = predict5,method="bootstrap",alternative="less")$p.value
  
  allauc=data.frame(riskfactors=auc0,metab=auc1,prs=auc2,metabprs=auc3,combined=auc5,
                    auc.pvalue0=auc.pvalue0,auc.pvalue1=auc.pvalue1,auc.pvalue2=auc.pvalue2,auc.pvalue3=auc.pvalue3)
  return(allauc)
}


#train using 1 cohort, validate using the other
#add metabolite score
load("../result/elasticnet_2cohorts.RData")
#744
trainATBCres=trainATBC$trainres
#720
testATBCres=trainATBC$testres
metabID=c(trainATBCres$ID,testATBCres$ID)
allavai=data.frame(ID=metabID)
idx=match(allavai$ID,pheno$allID)
allavai$COHORT=pheno$COHORT[idx]
allavai$CASE=pheno$CASE[idx]
allavai$PRS=NA
allavai$PRS[which(allavai$ID %in% allprsdat$allID)]=1
table(allavai$COHORT,allavai$CASE)
#        0   1
# ATBC 372 372
# PLCO 360 360
idx=which(allavai$PRS==1)
table(allavai$COHORT[idx],allavai$CASE[idx])
#        0   1
# ATBC 330 222
# PLCO 341 295
trainPLCOres=trainPLCO$trainres
testPLCOres=trainPLCO$testres
# auc2cohort_ATBC=compute_auccohort(traincohort="ATBC")
# # riskfactors      metab      prs  metabprs  combined
# #   0.5670799  0.6378009 0.6045827 0.6580324 0.6689728
auc2cohort_ATBC1=compute_auccohort(traincohort="ATBC",opt1="complete")
#   riskfactors     metab       prs metabprs combined
# 1   0.5611199 0.6340772 0.6045827 0.662528 0.670716
# auc2cohort_PLCO=compute_auccohort(traincohort="PLCO")
# #  riskfactors     metab       prs metabprs  combined
# #     0.557232 0.6067826 0.5650969 0.6262646 0.6236631
auc2cohort_PLCO1=compute_auccohort(traincohort="PLCO",opt1="complete")
save(auc2cohort_ATBC1,auc2cohort_PLCO1,file="../result/compute_AUC_includingPRS.RData")
