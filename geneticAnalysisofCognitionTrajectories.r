# Longitudinal analysis of cognition
# Will consider time in three different ways: i) standardized to baseline (i.e. study entry), ii) standardized to death and iii ) age. Analyses performed using a multi-level model, adjusted for age and sex with individual as a random effect.

library(lme4)
library(lmerTest)

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

## want to keep samples with either APOE geno or SNP data 
keepSamples<-pathDat$DNA_ID[!is.na(pathDat$APOE_geno) | pathDat$DNA_ID %in% prs.imp.noapoe$FID]  


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

## CAD data
cad<-read.csv("/mnt/data1/Gemma/BDR/CAD/CADAssessments.csv", stringsAsFactors = FALSE)

cad<-cad[cad$BBNId %in% pathDat$BBNId,]
cad$VISIT_DATE<-as.Date(cad$VISIT_DATE, format = "%d-%b-%Y")
cad$MRC_DoD<-as.Date(cad$MRC_DoD, format = "%d-%b-%Y") ## nb dateofdeath variable appears to have typos
cad$TimePriorDeath<-as.numeric(cad$MRC_DoD-cad$VISIT_DATE)

## summarise number of visits per individual
nVisit<-table(cad$BBNId)
mean(as.numeric(nVisit))
sd(as.numeric(nVisit))

cad.first<-cad[which(cad$VisitNumber == 1),]

## due to prior filtering can't use VisitNumberDesc column
cad.last<-cad[match(paste(names(nVisit), nVisit), paste(cad$BBNId, cad$VisitNumber)),]

cad.first<-cad.first[match(intersect(cad.first$BBNId,cad.last$BBNId), cad.first$BBNId),]
cad.last<-cad.last[match(intersect(cad.first$BBNId,cad.last$BBNId), cad.last$BBNId),]

## summarise mean gap between visits per person
nVisit<-nVisit[cad.first$BBNId]
spanVisits<-as.numeric(cad.last$VISIT_DATE - cad.first$VISIT_DATE)
meangapVisits<-spanVisits/(nVisit-1)
spanVisits[which(spanVisits == 0)]<-NA ## exclude people who only had 1 visit
meangapVisits[is.na(spanVisits)]<-NA ## exclude people who only had 1 visit

mean(as.numeric(spanVisits), na.rm = TRUE)/365
sd(as.numeric(spanVisits), na.rm = TRUE)/365
mean(as.numeric(meangapVisits), na.rm = TRUE)/365
sd(as.numeric(meangapVisits), na.rm = TRUE)/365

cogVar<-c("CDR_Global_Score","MMSE", "MOCA")

## profile how scores change over course of study

timeP<-matrix(data = NA, ncol = 8, nrow = length(cogVar))
colnames(timeP)<-c("Baseline:Coeff", "Baseline:Pvalue", "Death:Coeff", "Death:Pvalue", "Age:Coeff", "Age:Pvalue", "NIDs", "NObs")
rownames(timeP)<-cogVar

for(each in cogVar){
	cad.sub<-cad[!is.na(cad[,each]),c("BBNId", "VISIT_DATE", "Age", "GENDER", each, "VisitNumber")]
	colnames(cad.sub)[5]<-"each"
	## identify first visit with this cognative measure for each individual
	visitNum<-aggregate(cad.sub$VisitNumber, by = list(cad.sub$BBNId), FUN = range)
	cad.sub.first<-cad[match(paste(visitNum[,1], visitNum$x[,1]), paste(cad$BBNId, cad$VisitNumber)),]
	cad.sub.last<-cad[match(paste(visitNum[,1], visitNum$x[,2]), paste(cad$BBNId, cad$VisitNumber)),]
	
	## set baseline to first visit with this cognitive measure NB may not be first visit
	timeBaseline<-as.numeric(cad.sub$VISIT_DATE - cad.sub.first$VISIT_DATE[match(cad.sub$BBNId, cad.sub.first$BBNId)])
	timeDeath<-as.numeric(cad.sub$VISIT_DATE - cad.sub.last$MRC_DoD[match(cad.sub$BBNId, cad.sub.last$BBNId)])

	index<-match(cad.sub$BBNId, pathDat$BBNId)
	model1<-lmer(each ~ timeBaseline + Age + GENDER + (1 | BBNId), data = cad.sub)
	model2<-lmer(each ~ timeDeath + Age + GENDER + (1 | BBNId), data = cad.sub)
	model3<-lmer(each ~ Age + GENDER + (1 | BBNId), data = cad.sub)

	timeP[each,1:2]<-summary(model1)$coefficients["timeBaseline", c(1,5)]
	timeP[each,3:4]<-summary(model2)$coefficients["timeDeath", c(1,5)]
	timeP[each,5:6]<-summary(model3)$coefficients["Age", c(1,5)]
	
	timeP[each,7]<-summary(model1)$ngrps[1]
	timeP[each,8]<-length(summary(model1)$residuals)
}

## APOE Analysis of Cognitive Trajectory

timeP<-matrix(data = NA, ncol = 32, nrow = length(cogVar))
colnames(timeP)<-c("Baseline:Coeff", "Baseline:Pvalue", "E2:Coeff", "E2:Pvalue", "E4:Coeff", "E4:Pvalue", "BaselineIntE2:Coeff", "BaselineIntE2:Pvalue", "BaselineIntE4:Coeff", "BaselineIntE4:Pvalue", "Death:Coeff", "Death:Pvalue", "E2:Coeff", "E2:Pvalue", "E4:Coeff", "E4:Pvalue", "DeathIntE2:Coeff", "DeathIntE2:Pvalue", "DeathIntE4:Coeff", "DeathIntE4:Pvalue", "Age:Coeff", "Age:Pvalue",  "E2:Coeff", "E2:Pvalue","E4:Coeff", "E4:Pvalue","AgeIntE2:Coeff", "AgeIntE2:Pvalue","AgeIntE4:Coeff", "AgeIntE4:Pvalue", "NIDs", "NObs")
rownames(timeP)<-cogVar

for(each in cogVar){
	cad.sub<-cad[!is.na(cad[,each]),c("BBNId", "VISIT_DATE", "Age", "GENDER", each, "VisitNumber")]
	colnames(cad.sub)[5]<-"each"
	## identify first visit with this cognative measure for each individual
	visitNum<-aggregate(cad.sub$VisitNumber, by = list(cad.sub$BBNId), FUN = range)
	cad.sub.first<-cad[match(paste(visitNum[,1], visitNum$x[,1]), paste(cad$BBNId, cad$VisitNumber)),]
	cad.sub.last<-cad[match(paste(visitNum[,1], visitNum$x[,2]), paste(cad$BBNId, cad$VisitNumber)),]
	
	## set baseline to first visit with this cognitive measure NB may not be first visit
	timeBaseline<-as.numeric(cad.sub$VISIT_DATE - cad.sub.first$VISIT_DATE[match(cad.sub$BBNId, cad.sub.first$BBNId)])
	timeDeath<-as.numeric(cad.sub$VISIT_DATE - cad.sub.last$MRC_DoD[match(cad.sub$BBNId, cad.sub.last$BBNId)])

	index<-match(cad.sub$BBNId, pathDat$BBNId)
	model1<-lmer(each ~ timeBaseline + Age + GENDER + pathDat$n_e2[index] + pathDat$n_e2[index] * timeBaseline + pathDat$n_e4[index] + pathDat$n_e4[index] * timeBaseline + (1 | BBNId), data = cad.sub)
	model2<-lmer(each ~ timeDeath + Age + GENDER + pathDat$n_e2[index] + pathDat$n_e2[index] * timeDeath + pathDat$n_e4[index] + pathDat$n_e4[index] * timeDeath+ (1 | BBNId), data = cad.sub)
	model3<-lmer(each ~ Age + GENDER + pathDat$n_e2[index] + pathDat$n_e2[index] * Age + pathDat$n_e4[index] + pathDat$n_e4[index] * Age  + (1 | BBNId), data = cad.sub)

	timeP[each,1:2]<-summary(model1)$coefficients["timeBaseline", c(1,5)]
	timeP[each,3:4]<-summary(model1)$coefficients["pathDat$n_e2[index]", c(1,5)]
	timeP[each,5:6]<-summary(model1)$coefficients["pathDat$n_e4[index]", c(1,5)]
	timeP[each,7:8]<-summary(model1)$coefficients["timeBaseline:pathDat$n_e2[index]", c(1,5)]
	timeP[each,9:10]<-summary(model1)$coefficients["timeBaseline:pathDat$n_e4[index]", c(1,5)]

	timeP[each,11:12]<-summary(model2)$coefficients["timeDeath", c(1,5)]
	timeP[each,13:14]<-summary(model2)$coefficients["pathDat$n_e2[index]", c(1,5)]
	timeP[each,15:16]<-summary(model2)$coefficients["pathDat$n_e4[index]", c(1,5)]
	timeP[each,17:18]<-summary(model2)$coefficients["timeDeath:pathDat$n_e2[index]", c(1,5)]
	timeP[each,19:20]<-summary(model2)$coefficients["timeDeath:pathDat$n_e4[index]", c(1,5)]

	timeP[each,21:22]<-summary(model3)$coefficients["Age", c(1,5)]
	timeP[each,23:24]<-summary(model3)$coefficients["pathDat$n_e2[index]", c(1,5)]
	timeP[each,25:26]<-summary(model3)$coefficients["pathDat$n_e4[index]", c(1,5)]
	timeP[each,27:28]<-summary(model3)$coefficients["Age:pathDat$n_e2[index]", c(1,5)]
	timeP[each,29:30]<-summary(model3)$coefficients["Age:pathDat$n_e4[index]", c(1,5)]
	
	timeP[each,31]<-summary(model1)$ngrps[1]
	timeP[each,32]<-length(summary(model1)$residuals)
}
write.csv(timeP, "/mnt/data1/Gemma/BDR/BDR_genetic/Analysis/CognitiveTrajectoryInteractionWithAPOE.csv")


### test for different trajectories by AD PRS.

timeP<-matrix(data = NA, ncol = 20, nrow = length(cogVar))
colnames(timeP)<-c("Baseline:Coeff", "Baseline:Pvalue", "PRS:Coeff", "PRS:Pvalue", "BaselineIntPRS:Coeff", "BaselineIntPRS:Pvalue",  "Death:Coeff", "Death:Pvalue", "PRS:Coeff", "PRS:Pvalue","DeathIntPRS:Coeff", "DeathIntPRS:Pvalue", "Age:Coeff", "Age:Pvalue",  "PRS:Coeff", "PRS:Pvalue","AgeIntPRS:Coeff", "AgeIntPRS:Pvalue", "NIDs", "NObs")
rownames(timeP)<-cogVar

for(each in cogVar){
	cad.sub<-cad[!is.na(cad[,each]),c("BBNId", "VISIT_DATE", "Age", "GENDER", each, "VisitNumber")]
	colnames(cad.sub)[5]<-"each"
	## identify first visit with this cognative measure for each individual
	visitNum<-aggregate(cad.sub$VisitNumber, by = list(cad.sub$BBNId), FUN = range)
	cad.sub.first<-cad[match(paste(visitNum[,1], visitNum$x[,1]), paste(cad$BBNId, cad$VisitNumber)),]
	cad.sub.last<-cad[match(paste(visitNum[,1], visitNum$x[,2]), paste(cad$BBNId, cad$VisitNumber)),]
	
	## set baseline to first visit with this cognitive measure NB may not be first visit
	timeBaseline<-as.numeric(cad.sub$VISIT_DATE - cad.sub.first$VISIT_DATE[match(cad.sub$BBNId, cad.sub.first$BBNId)])
	timeDeath<-as.numeric(cad.sub$VISIT_DATE - cad.sub.last$MRC_DoD[match(cad.sub$BBNId, cad.sub.last$BBNId)])

	index<-match(cad.sub$BBNId, pathDat$BBNId)
	model1<-lmer(each ~ timeBaseline + Age + GENDER + prs.scale[index] + prs.scale[index] * timeBaseline + (1 | BBNId), data = cad.sub)
	model2<-lmer(each ~ timeDeath + Age + GENDER + prs.scale[index] + prs.scale[index] * timeDeath + (1 | BBNId), data = cad.sub)
	model3<-lmer(each ~ Age + GENDER + prs.scale[index] + prs.scale[index] * Age + (1 | BBNId), data = cad.sub)

	timeP[each,1:2]<-summary(model1)$coefficients["timeBaseline", c(1,5)]
	timeP[each,3:4]<-summary(model1)$coefficients["prs.scale[index]", c(1,5)]
	timeP[each,5:6]<-summary(model1)$coefficients["timeBaseline:prs.scale[index]", c(1,5)]

	timeP[each,7:8]<-summary(model2)$coefficients["timeDeath", c(1,5)]
	timeP[each,9:10]<-summary(model2)$coefficients["prs.scale[index]", c(1,5)]
	timeP[each,11:12]<-summary(model2)$coefficients["timeDeath:prs.scale[index]", c(1,5)]

	timeP[each,13:14]<-summary(model3)$coefficients["Age", c(1,5)]
	timeP[each,15:16]<-summary(model3)$coefficients["prs.scale[index]", c(1,5)]
	timeP[each,17:18]<-summary(model3)$coefficients["Age:prs.scale[index]", c(1,5)]
	
	timeP[each,19]<-summary(model1)$ngrps[1]
	timeP[each,20]<-length(summary(model1)$residuals)
}

write.csv(timeP, "/mnt/data1/Gemma/BDR/BDR_genetic/Analysis/CognitiveTrajectoryInteractionWithPRS.csv")


## how does neuropathology influence cognitive trajectory over the clinical visits

timeP<-matrix(data = NA, ncol = 20, nrow = length(cogVar))
colnames(timeP)<-c("Baseline:Coeff", "Baseline:Pvalue", "Braak:Coeff", "Braak:Pvalue", "BaselineInt:Coeff", "BaselineInt:Pvalue", "Death:Coeff", "Death:Pvalue", "Braak:Coeff", "Braak:Pvalue", "DeathInt:Coeff", "DeathInt:Pvalue", "Age:Coeff", "Age:Pvalue",  "Braak:Coeff", "Braak:Pvalue","AgeInt:Coeff", "AgeInt:Pvalue", "NIDs", "NObs")
rownames(timeP)<-cogVar

for(each in cogVar){
	cad.sub<-cad[!is.na(cad[,each]),c("BBNId", "VISIT_DATE", "Age", "GENDER", each, "VisitNumber")]
	colnames(cad.sub)[5]<-"each"
	## identify first visit with this cognative measure for each individual
	visitNum<-aggregate(cad.sub$VisitNumber, by = list(cad.sub$BBNId), FUN = range)
	cad.sub.first<-cad[match(paste(visitNum[,1], visitNum$x[,1]), paste(cad$BBNId, cad$VisitNumber)),]
	cad.sub.last<-cad[match(paste(visitNum[,1], visitNum$x[,2]), paste(cad$BBNId, cad$VisitNumber)),]
	
	## set baseline to first visit with this cognitive measure NB may not be first visit
	timeBaseline<-as.numeric(cad.sub$VISIT_DATE - cad.sub.first$VISIT_DATE[match(cad.sub$BBNId, cad.sub.first$BBNId)])
	timeDeath<-as.numeric(cad.sub$VISIT_DATE - cad.sub.last$MRC_DoD[match(cad.sub$BBNId, cad.sub.last$BBNId)])

	index<-match(cad.sub$BBNId, pathDat$BBNId)
	model1<-lmer(each ~ timeBaseline + Age + GENDER + pathDat$Braak_tangle[index] + pathDat$Braak_tangle[index] * timeBaseline + (1 | BBNId), data = cad.sub)
	model2<-lmer(each ~ timeDeath + Age + GENDER + pathDat$Braak_tangle[index] + pathDat$Braak_tangle[index] * timeDeath + (1 | BBNId), data = cad.sub)
	model3<-lmer(each ~ Age + GENDER + pathDat$Braak_tangle[index] + pathDat$Braak_tangle[index] * Age + (1 | BBNId), data = cad.sub)

	timeP[each,1:2]<-summary(model1)$coefficients["timeBaseline", c(1,5)]
	timeP[each,3:4]<-summary(model1)$coefficients["pathDat$Braak_tangle[index]", c(1,5)]
	timeP[each,5:6]<-summary(model1)$coefficients["timeBaseline:pathDat$Braak_tangle[index]", c(1,5)]
	timeP[each,7:8]<-summary(model2)$coefficients["timeDeath", c(1,5)]
	timeP[each,9:10]<-summary(model2)$coefficients["pathDat$Braak_tangle[index]", c(1,5)]
	timeP[each,11:12]<-summary(model2)$coefficients["timeDeath:pathDat$Braak_tangle[index]", c(1,5)]
	timeP[each,13:14]<-summary(model3)$coefficients["Age", c(1,5)]
	timeP[each,15:16]<-summary(model3)$coefficients["pathDat$Braak_tangle[index]", c(1,5)]
	timeP[each,17:18]<-summary(model3)$coefficients["Age:pathDat$Braak_tangle[index]", c(1,5)]
	
	timeP[each,19]<-summary(model1)$ngrps[1]
	timeP[each,20]<-length(summary(model1)$residuals)
}

write.csv(timeP, "/mnt/data1/Gemma/BDR/BDR_genetic/Analysis/CognitiveTrajectoryInteractionWithBraak.csv")

## does neuropathology drive cognitive differences by APOE


timeP<-matrix(data = NA, ncol = 32, nrow = length(cogVar))
colnames(timeP)<-c("Baseline:Coeff", "Baseline:Pvalue", "E2:Coeff", "E2:Pvalue", "E4:Coeff", "E4:Pvalue", "BaselineIntE2:Coeff", "BaselineIntE2:Pvalue", "BaselineIntE4:Coeff", "BaselineIntE4:Pvalue", "Death:Coeff", "Death:Pvalue", "E2:Coeff", "E2:Pvalue", "E4:Coeff", "E4:Pvalue", "DeathIntE2:Coeff", "DeathIntE2:Pvalue", "DeathIntE4:Coeff", "DeathIntE4:Pvalue", "Age:Coeff", "Age:Pvalue",  "E2:Coeff", "E2:Pvalue","E4:Coeff", "E4:Pvalue","AgeIntE2:Coeff", "AgeIntE2:Pvalue","AgeIntE4:Coeff", "AgeIntE4:Pvalue", "NIDs", "NObs")
rownames(timeP)<-cogVar

for(each in cogVar){
	cad.sub<-cad[!is.na(cad[,each]),c("BBNId", "VISIT_DATE", "Age", "GENDER", each, "VisitNumber")]
	colnames(cad.sub)[5]<-"each"
	## identify first visit with this cognative measure for each individual
	visitNum<-aggregate(cad.sub$VisitNumber, by = list(cad.sub$BBNId), FUN = range)
	cad.sub.first<-cad[match(paste(visitNum[,1], visitNum$x[,1]), paste(cad$BBNId, cad$VisitNumber)),]
	cad.sub.last<-cad[match(paste(visitNum[,1], visitNum$x[,2]), paste(cad$BBNId, cad$VisitNumber)),]
	
	## set baseline to first visit with this cognitive measure NB may not be first visit
	timeBaseline<-as.numeric(cad.sub$VISIT_DATE - cad.sub.first$VISIT_DATE[match(cad.sub$BBNId, cad.sub.first$BBNId)])
	timeDeath<-as.numeric(cad.sub$VISIT_DATE - cad.sub.last$MRC_DoD[match(cad.sub$BBNId, cad.sub.last$BBNId)])

	index<-match(cad.sub$BBNId, pathDat$BBNId)
	model1<-lmer(each ~ timeBaseline + Age + GENDER + pathDat$n_e2[index] + pathDat$n_e2[index] * timeBaseline + pathDat$n_e4[index] + pathDat$n_e4[index] * timeBaseline + pathDat$Braak_tangle[index] + pathDat$Braak_tangle[index] * timeBaseline + (1 | BBNId), data = cad.sub)
	
	model3<-lmer(each ~ Age + GENDER + pathDat$n_e2[index] + pathDat$n_e2[index] * Age + pathDat$n_e4[index] + pathDat$n_e4[index] * Age  + pathDat$Braak_tangle[index] + pathDat$Braak_tangle[index] * Age +(1 | BBNId), data = cad.sub)

	timeP[each,1:2]<-summary(model1)$coefficients["timeBaseline", c(1,5)]
	timeP[each,3:4]<-summary(model1)$coefficients["pathDat$n_e2[index]", c(1,5)]
	timeP[each,5:6]<-summary(model1)$coefficients["pathDat$n_e4[index]", c(1,5)]
	timeP[each,7:8]<-summary(model1)$coefficients["timeBaseline:pathDat$n_e2[index]", c(1,5)]
	timeP[each,9:10]<-summary(model1)$coefficients["timeBaseline:pathDat$n_e4[index]", c(1,5)]

	timeP[each,21:22]<-summary(model3)$coefficients["Age", c(1,5)]
	timeP[each,23:24]<-summary(model3)$coefficients["pathDat$n_e2[index]", c(1,5)]
	timeP[each,25:26]<-summary(model3)$coefficients["pathDat$n_e4[index]", c(1,5)]
	timeP[each,27:28]<-summary(model3)$coefficients["Age:pathDat$n_e2[index]", c(1,5)]
	timeP[each,29:30]<-summary(model3)$coefficients["Age:pathDat$n_e4[index]", c(1,5)]
	
	timeP[each,31]<-summary(model1)$ngrps[1]
	timeP[each,32]<-length(summary(model1)$residuals)
}
write.csv(timeP, "/mnt/data1/Gemma/BDR/BDR_genetic/Analysis/CognitiveTrajectoryInteractionWithAPOEControlBraak.csv")
