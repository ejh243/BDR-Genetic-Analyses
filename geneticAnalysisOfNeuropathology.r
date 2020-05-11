## script to perform genetic analysis of quantitative neuropathology variables
## fits following 3 models
## model 1: neuropathology ~ apoe + covariates + genetic PCs
## model 2: neuropathology ~ prs + covariates + genetic PCs
## model 3: neuropathology ~ apoe + prs + covariates + genetic PCs
## model 4: neuropathology ~ apoe * prs + covariates + genetic PCs
## apoe is primarily modelled as number of E2 and E4 alleles, where e2/e2 is counted as n_e2 = 1
## sensitivity analyses i) modelling apoe as a factor
rsquared <- function(created_model, null) {
  dev <- created_model$deviance
  null_dev <- null$deviance
  model_n <- length(created_model$fitted.values)
  R_l <- 1 - dev / null_dev
  R_cs <- 1 - exp(-(null_dev - dev) / model_n)
  R_n <- R_cs / (1 - exp(-(null_dev / model_n)))
  return(R_n)
}

pMT<-0.05/33

library(corrplot)

quantPath<-c("Braak_tangle","Thal_amyloid","cerad_numeric","LB_stage", "TDP43")

setwd("/mnt/data1/BDR/")
prs.imp.noapoe<-read.table("Genetic/PRS/Imputed/BDR_kunkle_Imputed_EUR_Unrelated_no_APOE.all.score",  stringsAsFactors = FALSE, header = TRUE)

genoPCs<-read.table("Genetic/Imputed/FINAL/BDR_imputed_PCs.eigenvec", stringsAsFactors = FALSE, header = TRUE)

## pathology data
pathDat<-read.csv("pheno_files/pathology/BDR_pathology_data_for_analysis.csv", stringsAsFactors = FALSE, na.strings = c("NA", "<NA>"))
#exclude indidivuals with no DNA_ID
pathDat<-pathDat[which(pathDat$DNA_ID != ""),]

## where APOE is missing replace with status derived from SNP chip data
apoeGeno<-read.csv("/mnt/data1/BDR/Genetic/BDR_all/FINAL/APOEstatusfromSNPdata.csv", stringsAsFactors = FALSE, header = TRUE, row.names = 1)
apoeGeno<-apoeGeno[match(pathDat$DNA_ID, apoeGeno$IID),]
pathDat$APOE_geno[which(is.na(pathDat$APOE_geno))]<-apoeGeno$apoeStatus3[which(is.na(pathDat$APOE_geno))]

## only keep samples with SNP data as needed for genetic PCs 
keepSamples<-pathDat$DNA_ID[pathDat$DNA_ID %in% prs.imp.noapoe$FID]  

## filter path data to analytical sample
pathDat<-pathDat[match(keepSamples, pathDat$DNA_ID),]

pathDat$APOE_geno <- as.factor(pathDat$APOE_geno)
pathDat$BDR_Centre_key <- as.factor(pathDat$BDR_Centre_key)

## conver cerad density to a number
cerad_numeric<-rep(NA, nrow(pathDat))
cerad_numeric[which(pathDat$cerad_densitiy == "no")]<-0
cerad_numeric[which(pathDat$cerad_densitiy == "sparse")]<-1
cerad_numeric[which(pathDat$cerad_densitiy == "moderate")]<-2
cerad_numeric[which(pathDat$cerad_densitiy == "high")]<-3
pathDat<-cbind(pathDat, cerad_numeric)

## convert TDP-43 to factor
pathDat$TDP43<-as.factor(pathDat[,"TDP43"])

## count number of e2 and e4 alleles
## as e2/e2 very r are group together n=1 and n =2
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

## standardize PRS
PT<-3 ## column with PRS in
prs.scale<-scale(prs.imp.noapoe[,PT])

## correlation heatmap of neuropathological variables
corDat<-as.matrix(pathDat[,quantPath])
mode(corDat)<-"numeric"
corMat<-cor(corDat, use = "p")
colnames(corMat)<-c("Braak NFT stage", "Thal amyloid stage", "Braak LB stage", "CERAD stage", "TDP-43")
rownames(corMat)<-c("Braak NFT stage", "Thal amyloid stage", "Braak LB stage", "CERAD stage", "TDP-43")
pdf("/mnt/data1/Gemma/BDR/BDR_genetic/Analysis/SupplementaryFigure1.pdf")
corrplot(corMat, method = "number", type = "lower")
dev.off()

## summarise demographics and variables


sum.demo<-rbind(c(sum(pathDat[,"Sex"] == "male")/sum(!is.na(pathDat[,"Sex"]))*100,NA, NA, sum(!is.na(pathDat[,"Sex"]))),
c(NA, mean(pathDat[,"Age"]), sd(pathDat[,"Age"]), sum(!is.na(pathDat[,"Age"]))))
for(each in quantPath){
	if(class(pathDat[,each]) == "numeric"){
		sum.demo<-rbind(sum.demo, 
		c(NA, mean(pathDat[,each], na.rm = TRUE), sd(pathDat[,each], na.rm = TRUE), sum(!is.na(pathDat[,each]))))
	} else {
		sum.demo<-rbind(sum.demo, 
		c(sum(pathDat[,each] == 1, na.rm = TRUE)/sum(!is.na(pathDat[,each]))*100,NA, NA, sum(!is.na(pathDat[,each]))))
	}
}
rownames(sum.demo)<-c("Sex", "Age", quantPath)
colnames(sum.demo)<-c("%", "Mean", "SD", "N")
write.csv(signif(sum.demo, 3), "/mnt/data1/Gemma/BDR/BDR_genetic/Analysis/CohortSummary.csv")

apoeP<-matrix(NA, nrow = length(quantPath),ncol = 6)
rownames(apoeP)<-quantPath
colnames(apoeP)<-c("e2:Pvalue", "e2:Coeff","e2:%VarExp", "e4:Pvalue", "e4:Coeff","e4:%VarExp")
for(each in quantPath){
	if(class(pathDat[,each]) != "factor"){
		model<-lm(pathDat[,each] ~ pathDat$n_e2 + pathDat$n_e4 + pathDat$Age + pathDat$Sex + pathDat$BDR_Centre_key + genoPCs$C1 + genoPCs$C2 + genoPCs$C3 + genoPCs$C4 + genoPCs$C5 + genoPCs$C6 + genoPCs$C7 + genoPCs$C8)
		sumSq<-anova(model)["Sum Sq"]
		apoeP[each,3]<-sumSq[1,1]/sum(sumSq)*100
		apoeP[each,6]<-sumSq[2,1]/sum(sumSq)*100
	} else {
		model<-glm(pathDat[,each] ~ pathDat$n_e2 + pathDat$n_e4 + pathDat$Age + pathDat$Sex + pathDat$BDR_Centre_key + genoPCs$C1 + genoPCs$C2 + genoPCs$C3 + genoPCs$C4 + genoPCs$C5 + genoPCs$C6 + genoPCs$C7 + genoPCs$C8, family = binomial("logit"))
		nulle2<-glm(pathDat[,each] ~ pathDat$n_e4 + pathDat$Age + pathDat$Sex + pathDat$BDR_Centre_key + genoPCs$C1 + genoPCs$C2 + genoPCs$C3 + genoPCs$C4 + genoPCs$C5 + genoPCs$C6 + genoPCs$C7 + genoPCs$C8, family = binomial("logit"))
		nulle4<-glm(pathDat[,each] ~ pathDat$n_e2 + pathDat$Age + pathDat$Sex + pathDat$BDR_Centre_key + genoPCs$C1 + genoPCs$C2 + genoPCs$C3 + genoPCs$C4 + genoPCs$C5 + genoPCs$C6 + genoPCs$C7 + genoPCs$C8, family = binomial("logit"))
		apoeP[each,3]<-rsquared(model, nulle2)*100
		apoeP[each,6]<-rsquared(model, nulle4)*100
	}
	apoeP[each,c(2,1)]<-summary(model)$coefficients["pathDat$n_e2",c(1,4)]
	apoeP[each,c(5,4)]<-summary(model)$coefficients["pathDat$n_e4",c(1,4)]
}

write.csv(signif(apoeP,3), "/mnt/data1/Gemma/BDR/BDR_genetic/Analysis/APOEAnalysisNumeric.csv")

apoeP2<-matrix(NA, nrow = length(quantPath),ncol = 12)
rownames(apoeP2)<-quantPath
colnames(apoeP2)<-c("APOE:Pvalue", "APOE:%VarExp", "e2/e2:Pvalue","e2/e2:Coeff", "e2/e3:Pvalue","e2/e3:Coeff","e2/e4:Pvalue","e2/e4:Coeff", "e3/e4:Pvalue","e3/e4:Coeff", "e4/e4:Pvalue","e4/e4:Coeff")
for(each in quantPath){
	if(class(pathDat[,each]) != "factor"){
		model<-lm(pathDat[,each] ~ pathDat$APOE_geno + pathDat$Age + pathDat$Sex + pathDat$BDR_Centre_key + genoPCs$C1 + genoPCs$C2 + genoPCs$C3 + genoPCs$C4 + genoPCs$C5 + genoPCs$C6 + genoPCs$C7 + genoPCs$C8)
		null<-lm(pathDat[,each] ~ pathDat$Age + pathDat$Sex + pathDat$BDR_Centre_key + genoPCs$C1 + genoPCs$C2 + genoPCs$C3 + genoPCs$C4 + genoPCs$C5 + genoPCs$C6 + genoPCs$C7 + genoPCs$C8, subset = !is.na(pathDat$APOE_geno))
		sumSq<-anova(model)["Sum Sq"]
		apoeP2[each,2]<-sumSq[1,1]/sum(sumSq)*100
		apoeP2[each,1]<-anova(model, null)[2,6]
	} else {
		model<-glm(pathDat[,each] ~ pathDat$APOE_geno + pathDat$Age + pathDat$Sex + pathDat$BDR_Centre_key + genoPCs$C1 + genoPCs$C2 + genoPCs$C3 + genoPCs$C4 + genoPCs$C5 + genoPCs$C6 + genoPCs$C7 + genoPCs$C8, family = binomial("logit"))
		null<-glm(pathDat[,each] ~ pathDat$Age + pathDat$Sex + pathDat$BDR_Centre_key + genoPCs$C1 + genoPCs$C2 + genoPCs$C3 + genoPCs$C4 + genoPCs$C5 + genoPCs$C6 + genoPCs$C7 + genoPCs$C8, subset = !is.na(pathDat$APOE_geno), family = binomial("logit"))
		
		apoeP2[each,2]<-rsquared(model, null)*100
		apoeP2[each,1]<-anova(model, null, test = "Chisq")[2,5]
	}
	
	apoeP2[each,c(3,5,7,9,11)]<-summary(model)$coefficients[c("pathDat$APOE_genoe2/e2","pathDat$APOE_genoe2/e3","pathDat$APOE_genoe2/e4","pathDat$APOE_genoe3/e4","pathDat$APOE_genoe4/e4"),4]
	apoeP2[each,c(4,6,8,10,12)]<-summary(model)$coefficients[c("pathDat$APOE_genoe2/e2","pathDat$APOE_genoe2/e3","pathDat$APOE_genoe2/e4","pathDat$APOE_genoe3/e4","pathDat$APOE_genoe4/e4"),1]
 
}

write.csv(signif(apoeP2,3), "/mnt/data1/Gemma/BDR/BDR_genetic/Analysis/APOEAnalysisGenotype.csv")


prsP<-matrix(NA, nrow = length(quantPath),ncol = 3)
rownames(prsP)<-quantPath
colnames(prsP)<-c("Pvalue", "%VarExp", "Coeff")

for(each in quantPath){
	if(class(pathDat[,each]) != "factor"){
	
		model<-lm(pathDat[,each] ~ prs.scale + pathDat$Sex  + pathDat$Age + pathDat$BDR_Centre_key + genoPCs$C1 + genoPCs$C2 + genoPCs$C3 + genoPCs$C4 + genoPCs$C5 + genoPCs$C6 + genoPCs$C7 + genoPCs$C8)
		sumSq<-anova(model)["Sum Sq"]
		prsP[each,2]<-sumSq[1,1]/sum(sumSq)*100
	} else {
		model<-glm(pathDat[,each] ~ prs.scale + pathDat$Sex  + pathDat$Age + pathDat$BDR_Centre_key + genoPCs$C1 + genoPCs$C2 + genoPCs$C3 + genoPCs$C4 + genoPCs$C5 + genoPCs$C6 + genoPCs$C7 + genoPCs$C8, family = binomial("logit"))
		null<-glm(pathDat[,each] ~ pathDat$Sex  + pathDat$Age + pathDat$BDR_Centre_key + genoPCs$C1 + genoPCs$C2 + genoPCs$C3 + genoPCs$C4 + genoPCs$C5 + genoPCs$C6 + genoPCs$C7 + genoPCs$C8, family = binomial("logit"))
		prsP[each,2]<-rsquared(model,null)*100
	}
	prsP[each,c(3,1)]<-summary(model)$coefficients["prs.scale",c(1,4)]
}

write.csv(signif(prsP,3), "/mnt/data1/Gemma/BDR/BDR_genetic/Analysis/PRSAnalysis.csv")

## fit joint model with both APOE and PRS
## only test those with significant apoe and prs effect
jointPath<-quantPath[which((apoeP[,1] < pMT | apoeP[,4] < pMT) & prsP[,1] < pMT)]

jointP<-matrix(NA, nrow = length(jointPath),ncol = 9)
rownames(jointP)<-jointPath
colnames(jointP)<-c("e2:Pvalue", "e2:Coeff","e2:%VarExp", "e4:Pvalue", "e4:Coeff","e4:%VarExp","PRS:Pvalue", "PRS:Coeff", "PRS:%VarExp")
for(each in jointPath){
	model<-lm(pathDat[,each] ~ pathDat$n_e2 + pathDat$n_e4 + prs.scale + pathDat$Age + pathDat$Sex + pathDat$BDR_Centre_key + genoPCs$C1 + genoPCs$C2 + genoPCs$C3 + genoPCs$C4 + genoPCs$C5 + genoPCs$C6 + genoPCs$C7 + genoPCs$C8)
	sumSq<-anova(model)["Sum Sq"]
	jointP[each,3]<-sumSq[1,1]/sum(sumSq)*100
	jointP[each,6]<-sumSq[2,1]/sum(sumSq)*100
	jointP[each,9]<-sumSq[3,1]/sum(sumSq)*100
	jointP[each,c(8,7)]<-summary(model)$coefficients["prs.scale",c(1,4)]
	jointP[each,c(2,1)]<-summary(model)$coefficients["pathDat$n_e2",c(1,4)]
	jointP[each,c(5,4)]<-summary(model)$coefficients["pathDat$n_e4",c(1,4)]
}

write.csv(signif(jointP,3), "/mnt/data1/Gemma/BDR/BDR_genetic/Analysis/JointAnalysisAPOENumeric.csv")

jointP2<-matrix(NA, nrow = length(jointPath),ncol = 15)
rownames(jointP2)<-jointPath
colnames(jointP2)<-c("APOE:Pvalue", "APOE:%VarExp", "e2/e2:Pvalue","e2/e2:Coeff", "e2/e3:Pvalue","e2/e3:Coeff","e2/e4:Pvalue","e2/e4:Coeff", "e3/e4:Pvalue","e3/e4:Coeff", "e4/e4:Pvalue","e4/e4:Coeff","PRS:Pvalue", "PRS:Coeff", "PRS:%VarExp")
for(each in jointPath){
	model<-lm(pathDat[,each] ~ pathDat$APOE_geno  + prs.scale + pathDat$Age + pathDat$Sex + pathDat$BDR_Centre_key + genoPCs$C1 + genoPCs$C2 + genoPCs$C3 + genoPCs$C4 + genoPCs$C5 + genoPCs$C6 + genoPCs$C7 + genoPCs$C8)
	
	null<-lm(pathDat[,each] ~ prs.scale + pathDat$Age + pathDat$Sex + pathDat$BDR_Centre_key + genoPCs$C1 + genoPCs$C2 + genoPCs$C3 + genoPCs$C4 + genoPCs$C5 + genoPCs$C6 + genoPCs$C7 + genoPCs$C8, subset = !is.na(pathDat$APOE_geno))
	jointP2[each,1]<-anova(model, null)[2,6]
	sumSq<-anova(model)["Sum Sq"]
	jointP2[each,2]<-sumSq[1,1]/sum(sumSq)*100
	jointP2[each,15]<-sumSq[2,1]/sum(sumSq)*100
	
	jointP2[each,c(3,5,7,9,11)]<-summary(model)$coefficients[c("pathDat$APOE_genoe2/e2","pathDat$APOE_genoe2/e3","pathDat$APOE_genoe2/e4","pathDat$APOE_genoe3/e4","pathDat$APOE_genoe4/e4"),4]
	jointP2[each,c(4,6,8,10,12)]<-summary(model)$coefficients[c("pathDat$APOE_genoe2/e2","pathDat$APOE_genoe2/e3","pathDat$APOE_genoe2/e4","pathDat$APOE_genoe3/e4","pathDat$APOE_genoe4/e4"),1]
	jointP2[each,c(14,13)]<-summary(model)$coefficients["prs.scale",c(1,4)]

}

write.csv(signif(jointP2,3), "/mnt/data1/Gemma/BDR/BDR_genetic/Analysis/JointAnalysisAPOEGenotype.csv")

## fit multiplicative model

intP<-matrix(NA, nrow = length(quantPath),ncol = 13)
rownames(intP)<-quantPath
colnames(intP)<-c("e2:Pvalue", "e2:Coeff","e2:%VarExp", "e4:Pvalue", "e4:Coeff","e4:%VarExp","PRS:Pvalue", "PRS:Coeff", "PRS:%VarExp", "PRSxE2:P", "PRSxE2:%VarExp", "PRSxE4:P", "PRSxE4:%VarExp")

for(each in quantPath){
	if(class(pathDat[,each]) != "factor"){
		model<-lm(pathDat[,each] ~ pathDat$n_e2 + pathDat$n_e4 + prs.scale + prs.scale*pathDat$n_e2 + prs.scale*pathDat$n_e4 + pathDat$Age + pathDat$Sex + pathDat$BDR_Centre_key + genoPCs$C1 + genoPCs$C2 + genoPCs$C3 + genoPCs$C4 + genoPCs$C5 + genoPCs$C6 + genoPCs$C7 + genoPCs$C8)
		
		sumSq<-anova(model)["Sum Sq"]
		intP[each,3]<-sumSq[1,1]/sum(sumSq)*100
		intP[each,6]<-sumSq[2,1]/sum(sumSq)*100
		intP[each,9]<-sumSq[3,1]/sum(sumSq)*100
		intP[each,11]<-sumSq[7,1]/sum(sumSq)*100
		intP[each,13]<-sumSq[8,1]/sum(sumSq)*100
		null<-lm(pathDat[,each] ~ prs.scale + pathDat$Age + pathDat$Sex + pathDat$BDR_Centre_key + genoPCs$C1 + genoPCs$C2 + genoPCs$C3 + genoPCs$C4 + genoPCs$C5 + genoPCs$C6 + genoPCs$C7 + genoPCs$C8, subset = !is.na(pathDat$APOE_geno))
		intP[each,1]<-anova(model, null)[2,6]
	} else {
		model<-glm(pathDat[,each] ~ pathDat$n_e2 + pathDat$n_e4 + prs.scale + prs.scale*pathDat$n_e2 + prs.scale*pathDat$n_e4 + pathDat$Age + pathDat$Sex + pathDat$BDR_Centre_key + genoPCs$C1 + genoPCs$C2 + genoPCs$C3 + genoPCs$C4 + genoPCs$C5 + genoPCs$C6 + genoPCs$C7 + genoPCs$C8, family = binomial("logit"))
		
		null<-glm(pathDat[,each] ~ pathDat$n_e4 + prs.scale + prs.scale*pathDat$n_e2 + prs.scale*pathDat$n_e4+ pathDat$Age + pathDat$Sex + pathDat$BDR_Centre_key + genoPCs$C1 + genoPCs$C2 + genoPCs$C3 + genoPCs$C4 + genoPCs$C5 + genoPCs$C6 + genoPCs$C7 + genoPCs$C8, subset = !is.na(pathDat$APOE_geno), family = binomial("logit"))
		intP[each,3]<-rsquared(model, null)*100
		null<-glm(pathDat[,each] ~ pathDat$n_e2 + prs.scale + prs.scale*pathDat$n_e2 + prs.scale*pathDat$n_e4 + pathDat$Age + pathDat$Sex + pathDat$BDR_Centre_key + genoPCs$C1 + genoPCs$C2 + genoPCs$C3 + genoPCs$C4 + genoPCs$C5 + genoPCs$C6 + genoPCs$C7 + genoPCs$C8, subset = !is.na(pathDat$APOE_geno), family = binomial("logit"))
		intP[each,6]<-rsquared(model, null)*100
		null<-glm(pathDat[,each] ~ pathDat$n_e2 + pathDat$n_e4 + prs.scale*pathDat$n_e2 + prs.scale*pathDat$n_e4 + pathDat$Age + pathDat$Sex + pathDat$BDR_Centre_key + genoPCs$C1 + genoPCs$C2 + genoPCs$C3 + genoPCs$C4 + genoPCs$C5 + genoPCs$C6 + genoPCs$C7 + genoPCs$C8, subset = !is.na(pathDat$APOE_geno), family = binomial("logit"))	
		intP[each,9]<-rsquared(model, null)*100
		null<-glm(pathDat[,each] ~ pathDat$n_e2 + pathDat$n_e4 + prs.scale + prs.scale*pathDat$n_e4 + pathDat$Age + pathDat$Sex + pathDat$BDR_Centre_key + genoPCs$C1 + genoPCs$C2 + genoPCs$C3 + genoPCs$C4 + genoPCs$C5 + genoPCs$C6 + genoPCs$C7 + genoPCs$C8, subset = !is.na(pathDat$APOE_geno), family = binomial("logit"))
		intP[each,11]<-rsquared(model, null)*100
		null<-glm(pathDat[,each] ~ pathDat$n_e2 + pathDat$n_e4 + prs.scale + prs.scale*pathDat$n_e2 + pathDat$Age + pathDat$Sex + pathDat$BDR_Centre_key + genoPCs$C1 + genoPCs$C2 + genoPCs$C3 + genoPCs$C4 + genoPCs$C5 + genoPCs$C6 + genoPCs$C7 + genoPCs$C8, subset = !is.na(pathDat$APOE_geno), family = binomial("logit"))
		intP[each,13]<-rsquared(model, null)*100
		}
	intP[each,c(8,7)]<-summary(model)$coefficients["prs.scale",c(1,4)]
	intP[each,c(2,1)]<-summary(model)$coefficients["pathDat$n_e2",c(1,4)]
	intP[each,c(5,4)]<-summary(model)$coefficients["pathDat$n_e4",c(1,4)]
	intP[each,10]<-summary(model)$coefficients["pathDat$n_e2:prs.scale",4]
	intP[each,12]<-summary(model)$coefficients["pathDat$n_e4:prs.scale",4]
}
write.csv(signif(intP,3), "/mnt/data1/Gemma/BDR/BDR_genetic/Analysis/InteractionAnalysisAPOENumeric.csv")


int2P<-matrix(NA, nrow = length(quantPath),ncol = 7)
rownames(int2P)<-quantPath
colnames(int2P)<-c("APOE:Pvalue","APOE:%VarExp","PRS:Pvalue", "PRS:Coeff", "PRS:%VarExp", "PRSxAPOE:P", "PRSxAPOE:%VarExp")

for(each in quantPath){
	if(class(pathDat[,each]) != "factor"){
		model<-lm(pathDat[,each] ~ pathDat$APOE_geno + prs.scale + prs.scale*pathDat$APOE_geno  + pathDat$Age + pathDat$Sex + pathDat$BDR_Centre_key + genoPCs$C1 + genoPCs$C2 + genoPCs$C3 + genoPCs$C4 + genoPCs$C5 + genoPCs$C6 + genoPCs$C7 + genoPCs$C8)
		
		sumSq<-anova(model)["Sum Sq"]
		int2P[each,2]<-sumSq[1,1]/sum(sumSq)*100
		int2P[each,5]<-sumSq[2,1]/sum(sumSq)*100
		int2P[each,7]<-sumSq["pathDat$APOE_geno:prs.scale",1]/sum(sumSq)*100

		null<-lm(pathDat[,each] ~ pathDat$APOE_geno + prs.scale + pathDat$Age + pathDat$Sex + pathDat$BDR_Centre_key + genoPCs$C1 + genoPCs$C2 + genoPCs$C3 + genoPCs$C4 + genoPCs$C5 + genoPCs$C6 + genoPCs$C7 + genoPCs$C8, subset = !is.na(pathDat$APOE_geno))
		int2P[each,1]<-anova(model)["pathDat$APOE_geno",5]
		int2P[each,6]<-anova(model, null)[2,6]
	} else{
		model<-glm(pathDat[,each] ~ pathDat$APOE_geno + prs.scale + prs.scale*pathDat$APOE_geno  + pathDat$Age + pathDat$Sex + pathDat$BDR_Centre_key + genoPCs$C1 + genoPCs$C2 + genoPCs$C3 + genoPCs$C4 + genoPCs$C5 + genoPCs$C6 + genoPCs$C7 + genoPCs$C8, family = binomial("logit"))
		null<-glm(pathDat[,each] ~  prs.scale + prs.scale:pathDat$APOE_geno + pathDat$Age + pathDat$Sex + pathDat$BDR_Centre_key + genoPCs$C1 + genoPCs$C2 + genoPCs$C3 + genoPCs$C4 + genoPCs$C5 + genoPCs$C6 + genoPCs$C7 + genoPCs$C8, subset = !is.na(pathDat$APOE_geno), family = binomial("logit"))
		int2P[each,2]<-rsquared(model, null)*100
		int2P[each,1]<-anova(model, null, test = "Chisq")[2,5]
		null<-glm(pathDat[,each] ~  pathDat$APOE_geno + prs.scale*pathDat$APOE_geno + pathDat$Age + pathDat$Sex + pathDat$BDR_Centre_key + genoPCs$C1 + genoPCs$C2 + genoPCs$C3 + genoPCs$C4 + genoPCs$C5 + genoPCs$C6 + genoPCs$C7 + genoPCs$C8, subset = !is.na(pathDat$APOE_geno), family = binomial("logit"))
		int2P[each,5]<-rsquared(model, null)*100
		null<-glm(pathDat[,each] ~  pathDat$APOE_geno + prs.scale + pathDat$Age + pathDat$Sex + pathDat$BDR_Centre_key + genoPCs$C1 + genoPCs$C2 + genoPCs$C3 + genoPCs$C4 + genoPCs$C5 + genoPCs$C6 + genoPCs$C7 + genoPCs$C8, subset = !is.na(pathDat$APOE_geno), family = binomial("logit"))
		int2P[each,7]<-rsquared(model, null)*100
		int2P[each,6]<-anova(model, null, test = "Chisq")[2,5]
	}
	int2P[each,c(4,3)]<-summary(model)$coefficients["prs.scale",c(1,4)]
}

write.csv(signif(int2P,3), "/mnt/data1/Gemma/BDR/BDR_genetic/Analysis/InteractionAnalysisAPOEGenotype.csv")

## as Braak stage is most significant, correct tests of other variables for Braak stage
## also correct analysis of BRaak stage co
apoePAdj<-matrix(NA, nrow = length(quantPath)^2,ncol = 6)
rownames(apoePAdj)<-rep(quantPath, length(quantPath))
colnames(apoePAdj)<-c("e2:Pvalue", "e2:Coeff","e2:%VarExp", "e4:Pvalue", "e4:Coeff","e4:%VarExp")

prsPAdj<-matrix(NA, nrow = length(quantPath)^2,ncol = 3)
rownames(prsPAdj)<-rep(quantPath, length(quantPath))
colnames(prsPAdj)<-c("Pvalue", "%VarExp", "Coeff")

jointPAdj<-matrix(NA, nrow = length(quantPath)^2,ncol = 9)
rownames(jointPAdj)<-rep(quantPath, length(quantPath))
colnames(jointPAdj)<-c("e2:Pvalue", "e2:Coeff","e2:%VarExp", "e4:Pvalue", "e4:Coeff","e4:%VarExp","PRS:Pvalue", "PRS:Coeff", "PRS:%VarExp")
rowNum<-1
for(each in quantPath){
	for(condition in quantPath){
		if(each != condition){
			if(class(pathDat[,each]) != "factor"){
				model<-lm(pathDat[,condition] ~ pathDat$n_e2 + pathDat$n_e4 + pathDat[,each] + pathDat$Age + pathDat$Sex + pathDat$BDR_Centre_key + genoPCs$C1 + genoPCs$C2 + genoPCs$C3 + genoPCs$C4 + genoPCs$C5 + genoPCs$C6 + genoPCs$C7 + genoPCs$C8)
				sumSq<-anova(model)["Sum Sq"]
				apoePAdj[rowNum,3]<-sumSq[1,1]/sum(sumSq)*100
				apoePAdj[rowNum,6]<-sumSq[2,1]/sum(sumSq)*100
				
			} else {
				model<-glm(pathDat[,condition] ~ pathDat$n_e2 + pathDat$n_e4 + pathDat[,each] + pathDat$Age + pathDat$Sex + pathDat$BDR_Centre_key + genoPCs$C1 + genoPCs$C2 + genoPCs$C3 + genoPCs$C4 + genoPCs$C5 + genoPCs$C6 + genoPCs$C7 + genoPCs$C8,  family = binomial("logit"))
				nulle2<-glm(pathDat[,condition] ~ pathDat$n_e4 + pathDat[,each] + pathDat$Age + pathDat$Sex + pathDat$BDR_Centre_key + genoPCs$C1 + genoPCs$C2 + genoPCs$C3 + genoPCs$C4 + genoPCs$C5 + genoPCs$C6 + genoPCs$C7 + genoPCs$C8,  family = binomial("logit"))
				nulle4<-glm(pathDat[,condition] ~ pathDat$n_e2 + pathDat[,each] + pathDat$Age + pathDat$Sex + pathDat$BDR_Centre_key + genoPCs$C1 + genoPCs$C2 + genoPCs$C3 + genoPCs$C4 + genoPCs$C5 + genoPCs$C6 + genoPCs$C7 + genoPCs$C8,  family = binomial("logit"))
				apoePAdj[rowNum,3]<-rsquared(model, nulle2)*100
				apoePAdj[rowNum,6]<-rsquared(model, nulle4)*100
			}
			apoePAdj[rowNum,c(2,1)]<-summary(model)$coefficients["pathDat$n_e2",c(1,4)]
			apoePAdj[rowNum,c(5,4)]<-summary(model)$coefficients["pathDat$n_e4",c(1,4)]

			if(class(pathDat[,each]) != "factor"){
			model<-lm(pathDat[,condition] ~ prs.scale + pathDat[,each] + pathDat$Sex  + pathDat$Age + pathDat$BDR_Centre_key + genoPCs$C1 + genoPCs$C2 + genoPCs$C3 + genoPCs$C4 + genoPCs$C5 + genoPCs$C6 + genoPCs$C7 + genoPCs$C8)
			sumSq<-anova(model)["Sum Sq"]
			prsPAdj[rowNum,2]<-sumSq[1,1]/sum(sumSq)*100
			} else {
				model<-glm(pathDat[,condition] ~ prs.scale + pathDat[,each] + pathDat$Sex  + pathDat$Age + pathDat$BDR_Centre_key + genoPCs$C1 + genoPCs$C2 + genoPCs$C3 + genoPCs$C4 + genoPCs$C5 + genoPCs$C6 + genoPCs$C7 + genoPCs$C8,  family = binomial("logit"))
				null<-glm(pathDat[,condition] ~ pathDat[,each] + pathDat$Sex  + pathDat$Age + pathDat$BDR_Centre_key + genoPCs$C1 + genoPCs$C2 + genoPCs$C3 + genoPCs$C4 + genoPCs$C5 + genoPCs$C6 + genoPCs$C7 + genoPCs$C8,  family = binomial("logit"))
				prsPAdj[rowNum,2]<-rsquared(model, null)*100
			}
			prsPAdj[rowNum,c(3,1)]<-summary(model)$coefficients["prs.scale",c(1,4)]

			if(class(pathDat[,each]) != "factor"){
				model<-lm(pathDat[,condition] ~ pathDat$n_e2 + pathDat$n_e4 + prs.scale + pathDat[,each] + pathDat$Age + pathDat$Sex + pathDat$BDR_Centre_key + genoPCs$C1 + genoPCs$C2 + genoPCs$C3 + genoPCs$C4 + genoPCs$C5 + genoPCs$C6 + genoPCs$C7 + genoPCs$C8)
				sumSq<-anova(model)["Sum Sq"]
				jointPAdj[rowNum,3]<-sumSq[1,1]/sum(sumSq)*100
				jointPAdj[rowNum,6]<-sumSq[2,1]/sum(sumSq)*100
				jointPAdj[rowNum,9]<-sumSq[3,1]/sum(sumSq)*100
			} else {
				model<-glm(pathDat[,condition] ~ pathDat$n_e2 + pathDat$n_e4 + prs.scale + pathDat[,each] + pathDat$Age + pathDat$Sex + pathDat$BDR_Centre_key + genoPCs$C1 + genoPCs$C2 + genoPCs$C3 + genoPCs$C4 + genoPCs$C5 + genoPCs$C6 + genoPCs$C7 + genoPCs$C8,  family = binomial("logit"))
				null<-glm(pathDat[,condition] ~ pathDat$n_e4 + prs.scale + pathDat[,each] + pathDat$Age + pathDat$Sex + pathDat$BDR_Centre_key + genoPCs$C1 + genoPCs$C2 + genoPCs$C3 + genoPCs$C4 + genoPCs$C5 + genoPCs$C6 + genoPCs$C7 + genoPCs$C8,  family = binomial("logit"))
				jointPAdj[rowNum,3]<-rsquared(model, null)*100
				null<-glm(pathDat[,condition] ~ pathDat$n_e2 + prs.scale + pathDat[,each] + pathDat$Age + pathDat$Sex + pathDat$BDR_Centre_key + genoPCs$C1 + genoPCs$C2 + genoPCs$C3 + genoPCs$C4 + genoPCs$C5 + genoPCs$C6 + genoPCs$C7 + genoPCs$C8,  family = binomial("logit"))
				jointPAdj[rowNum,6]<-rsquared(model, null)*100
				null<-glm(pathDat[,condition] ~ pathDat$n_e2 + pathDat$n_e4 + pathDat[,each] + pathDat$Age + pathDat$Sex + pathDat$BDR_Centre_key + genoPCs$C1 + genoPCs$C2 + genoPCs$C3 + genoPCs$C4 + genoPCs$C5 + genoPCs$C6 + genoPCs$C7 + genoPCs$C8,  family = binomial("logit"))
				jointPAdj[rowNum,9]<-rsquared(model, null)*100
			}
			jointPAdj[rowNum,c(8,7)]<-summary(model)$coefficients["prs.scale",c(1,4)]
			jointPAdj[rowNum,c(2,1)]<-summary(model)$coefficients["pathDat$n_e2",c(1,4)]
			jointPAdj[rowNum,c(5,4)]<-summary(model)$coefficients["pathDat$n_e4",c(1,4)]
			
	}
	rowNum<-rowNum+1
	}
}
write.csv(signif(apoePAdj,3), "/mnt/data1/Gemma/BDR/BDR_genetic/Analysis/APOEAnalysisControlForNeuropath.csv")
write.csv(signif(prsPAdj,3), "/mnt/data1/Gemma/BDR/BDR_genetic/Analysis/PRSAnalysisControlForNeuropath.csv")
write.csv(signif(jointPAdj,3), "/mnt/data1/Gemma/BDR/BDR_genetic/Analysis/JointAnalysisControlForNeuropath.csv")
