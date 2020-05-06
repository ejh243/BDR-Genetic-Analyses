## Script to perform survival analysis comparing age at death by APOE status, and AD PRS
## use cox's proportional hazards models as allows inclusion of covariates and continuous variables
## all models include sex, centre as covariates and first 8 genetics PCs
## PRS is standardized to mean = 0 and sd =1 so interpretation is in units of SD PRS.


library("survival")
library("survminer")

setwd("/mnt/data1/BDR/")
prs.imp.noapoe<-read.table("Genetic/PRS/Imputed/BDR_kunkle_Imputed_EUR_Unrelated_no_APOE.all.score",  stringsAsFactors = FALSE, header = TRUE)

nprs.imp.noapoe<-nrow(prs.imp.noapoe)
genoPCs<-read.table("Genetic/Imputed/FINAL/BDR_imputed_PCs.eigenvec", stringsAsFactors = FALSE, header = TRUE)

## pathology data
pathDat<-read.csv("pheno_files/pathology/BDR_pathology_data_for_analysis.csv", stringsAsFactors = FALSE, na.strings = c("NA", "<NA>"))
nPath<-nrow(pathDat)
#exclude indidivuals with no DNA_ID
pathDat<-pathDat[which(pathDat$DNA_ID != ""),]

## where APOE is missing replace with status derived from SNP chip data
apoeGeno<-read.csv("/mnt/data1/BDR/Genetic/BDR_all/FINAL/APOEstatusfromSNPdata.csv", stringsAsFactors = FALSE, header = TRUE, row.names = 1)
apoeGeno<-apoeGeno[match(pathDat$DNA_ID, apoeGeno$IID),]
pathDat$APOE_geno[which(is.na(pathDat$APOE_geno))]<-apoeGeno$apoeStatus3[which(is.na(pathDat$APOE_geno))]


## want to keep samples with either APOE geno or SNP data 
keepSamples<-pathDat$DNA_ID[!is.na(pathDat$APOE_geno) | pathDat$DNA_ID %in% prs.imp.noapoe$FID]  

## filter path data to analytical sample
pathDat<-pathDat[match(keepSamples, pathDat$DNA_ID),]
pathDat$APOE_geno <- as.factor(pathDat$APOE_geno)
pathDat$BDR_Centre_key <- as.factor(pathDat$BDR_Centre_key)

## count number of e2 and e4 alleles
## as e2/e2 very rare group together n=1 and n =2
n_e2<-rep(0, nrow(pathDat))
n_e2[pathDat$APOE_geno %in% c("e2/e3", "e2/e4")]<-1
n_e2[pathDat$APOE_geno %in% c("e2/e2")]<-1
n_e2[is.na(pathDat$APOE_geno)]<-NA
n_e4<-rep(0, nrow(pathDat))
n_e4[pathDat$APOE_geno %in% c("e2/e4", "e3/e4")]<-1
n_e4[pathDat$APOE_geno %in% c("e4/e4")]<-2
n_e4[is.na(pathDat$APOE_geno)]<-NA
pathDat<-cbind(pathDat, n_e2, n_e4)

## refactor APOE so e3/e3 the baseline group
pathDat$APOE_geno<-factor(pathDat$APOE_geno, levels =c("e3/e3", "e2/e2", "e2/e3", "e2/e4", "e3/e4", "e4/e4"))

## filter genetic data to analytical sample
prs.imp.noapoe<-prs.imp.noapoe[match(keepSamples, prs.imp.noapoe$FID),]
genoPCs<-genoPCs[match(keepSamples, genoPCs$FID) ,] 

## scale prs variable
pathDat$PRS<-scale(prs.imp.noapoe$X5e.08)

## summary of sample numbers
table(c("NA", "APOE")[as.factor(!is.na(pathDat$APOE_geno))], c("NA", "PRS")[as.factor(!is.na(prs.imp.noapoe$FID))])

## create variable indicated censor status all identical as none of the data is censored
status<-rep(1, nrow(pathDat))
pathDat<-cbind(pathDat, status)

## survival analysis of age at death by apoe status
apoeModel<-coxph(Surv(Age, status) ~ n_e2 + n_e4 + Sex + BDR_Centre_key + genoPCs$C1 + genoPCs$C2 + genoPCs$C3 + genoPCs$C4 + genoPCs$C5 + genoPCs$C6 + genoPCs$C7 + genoPCs$C8, data = pathDat)
prsModel<-coxph(Surv(Age, status) ~ PRS + Sex + BDR_Centre_key + genoPCs$C1 + genoPCs$C2 + genoPCs$C3 + genoPCs$C4 + genoPCs$C5 + genoPCs$C6 + genoPCs$C7 + genoPCs$C8, data = pathDat)
jointModel<-coxph(Surv(Age, status) ~ n_e2 + n_e4 + PRS + Sex + BDR_Centre_key + genoPCs$C1 + genoPCs$C2 + genoPCs$C3 + genoPCs$C4 + genoPCs$C5 + genoPCs$C6 + genoPCs$C7 + genoPCs$C8, data = pathDat)

res.tab<-matrix(data = NA, nrow = 6, ncol = 4)
colnames(res.tab)<-c("HR", "SE", "Pval", "N")
rownames(res.tab)<-rep(c("nE2", "nE4", "PRS"),2)

res.tab[1,1:3]<-summary(apoeModel)$coefficients["n_e2",c(2,3,5)]
res.tab[2,1:3]<-summary(apoeModel)$coefficients["n_e4",c(2,3,5)]
res.tab[3,1:3]<-summary(prsModel)$coefficients["PRS",c(2,3,5)]

res.tab[1,4]<-apoeModel$n
res.tab[2,4]<-apoeModel$n
res.tab[3,4]<-prsModel$n

res.tab[4,1:3]<-summary(jointModel)$coefficients["n_e2",c(2,3,5)]
res.tab[5,1:3]<-summary(jointModel)$coefficients["n_e4",c(2,3,5)]
res.tab[6,1:3]<-summary(jointModel)$coefficients["PRS",c(2,3,5)]

res.tab[4,4]<-jointModel$n
res.tab[5,4]<-jointModel$n
res.tab[6,4]<-jointModel$n

write.csv(res.tab, "/mnt/data1/Gemma/BDR/BDR_genetic/Analysis/SurvivalAnalysisAgeAtDeath.csv")

## convert to a meaningful interpretation
quantile(pathDat$PRS, seq(0.1,0.9,0.1), na.rm = TRUE)