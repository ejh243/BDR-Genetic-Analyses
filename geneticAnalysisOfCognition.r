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

## aggregate CAD data
demo<-unique(pathDat[,c("BBNId", "DNA_ID", "Age", "Sex", "APOE_geno", "n_e2", "n_e4", "BDR_Centre_key", "Braak_tangle")])
cogVar<-c("MMSE", "MOCA", "CDR_Global_Score")

tabCogGenetics.apoe<-matrix(data = NA, ncol = 9, nrow = length(cogVar))
colnames(tabCogGenetics.apoe)<-c("N", "MeanTimeSinceVisit", "SETimeSinceVisit", "APOE:E2:Coeff", "APOE:E2:P", "APOE:E2:VarExp", "APOE:E4:Coeff", "APOE:E4:P", "APOE:E4:VarExp")
rownames(tabCogGenetics.apoe)<-cogVar

tabCogGenetics.apoe2<-matrix(data = NA, ncol = 15, nrow = length(cogVar))
colnames(tabCogGenetics.apoe2)<-c("N", "MeanTimeSinceVisit", "SETimeSinceVisit", "APOE:Pvalue", "APOE:%VarExp", "e2/e2:Pvalue","e2/e2:Coeff", "e2/e3:Pvalue","e2/e3:Coeff","e2/e4:Pvalue","e2/e4:Coeff", "e3/e4:Pvalue","e3/e4:Coeff", "e4/e4:Pvalue","e4/e4:Coeff")
rownames(tabCogGenetics.apoe2)<-cogVar

tabCogGenetics.prs<-matrix(data = NA, ncol = 6, nrow = length(cogVar))
colnames(tabCogGenetics.prs)<-c("N", "MeanTimeSinceVisit", "SETimeSinceVisit", "PRS:Coeff", "PRS:P", "PRS:VarExp")
rownames(tabCogGenetics.prs)<-cogVar

tabCogGenetics.joint<-matrix(data = NA, ncol = 12, nrow = length(cogVar))
colnames(tabCogGenetics.joint)<-c("N", "MeanTimeSinceVisit", "SETimeSinceVisit", "PRS:Coeff", "PRS:P", "PRS:VarExp","APOE:E2:Coeff", "APOE:E2:P", "APOE:E2:VarExp", "APOE:E4:Coeff", "APOE:E4:P", "APOE:E4:VarExp")
rownames(tabCogGenetics.joint)<-cogVar

tabCogGenetics.int<-matrix(data = NA, ncol = 16, nrow = length(cogVar))
colnames(tabCogGenetics.int)<-c("N", "MeanTimeSinceVisit", "SETimeSinceVisit", "e2:Pvalue", "e2:Coeff","e2:%VarExp", "e4:Pvalue", "e4:Coeff","e4:%VarExp","PRS:Pvalue", "PRS:Coeff", "PRS:%VarExp", "PRSxE2:P", "PRSxE2:%VarExp", "PRSxE4:P", "PRSxE4:%VarExp")
rownames(tabCogGenetics.int)<-cogVar


sum.demo<-rbind(c(sum(pathDat[,"Sex"] == "male")/sum(!is.na(pathDat[,"Sex"]))*100,NA, NA, sum(!is.na(pathDat[,"Sex"]))),
c(NA, mean(pathDat[,"Age"]), sd(pathDat[,"Age"]), sum(!is.na(pathDat[,"Age"]))))



cad<-cad[order(cad$VISIT_DATE),]
par(mfrow = c(2,3))
for(each in cogVar){

cad.tmp<-cad[!is.na(cad[,each]),]
nIds<-length(unique(cad.tmp$BBNId))
nVisitCog<-table(cad.tmp$BBNId)
scoreLast<-aggregate(cad.tmp[,each], by = list(cad.tmp$BBNId), FUN = tail, n = 1)
timeLast<-aggregate(cad.tmp$TimePriorDeath, by = list(cad.tmp$BBNId), FUN = tail, n = 1)
indexDemo<-match(scoreLast$Group.1, demo$BBNId)
dnaIds<-demo$DNA_ID[indexDemo]
indexGenetic<-match(dnaIds, prs.imp.noapoe$FID)

hist(scoreLast$x, xlab = each, ylab = "Number of Individuals", main = "")
  mtext(paste("mean =", signif(mean(scoreLast$x, na.rm = TRUE),3), "sd =", signif(sd(scoreLast$x, na.rm = TRUE),3), "N =", sum(!is.na(scoreLast$x))), side = 3, line = 0.5)
hist(timeLast$x, xlab = each, ylab = "Number of Individuals", main = "")
  mtext(paste("mean =", signif(mean(timeLast$x, na.rm = TRUE),3), "sd =", signif(sd(timeLast$x, na.rm = TRUE),3), "N =", sum(!is.na(timeLast$x))), side = 3, line = 0.5)
  
sum.demo<-rbind(sum.demo, 
	c(NA, mean(scoreLast$x, na.rm = TRUE), sd(scoreLast$x, na.rm = TRUE), sum(!is.na(scoreLast$x))))

##test apoe as number of e2 and number of e4/e4
model<-lm(scoreLast$x ~ demo$n_e2[indexDemo] + demo$n_e4[indexDemo] + timeLast$x + demo$Age[indexDemo] + demo$Sex[indexDemo]+ demo$BDR_Centre_key[indexDemo] + genoPCs$C1[indexGenetic] + genoPCs$C2[indexGenetic] + genoPCs$C3[indexGenetic] + genoPCs$C4[indexGenetic] + genoPCs$C5[indexGenetic] + genoPCs$C6[indexGenetic] + genoPCs$C7[indexGenetic] + genoPCs$C8[indexGenetic])
tabCogGenetics.apoe[each,1]<-length(residuals(model))
tabCogGenetics.apoe[each,2]<-mean(timeLast$x[!is.na(demo$n_e2[indexDemo])], na.rm = TRUE)
tabCogGenetics.apoe[each,3]<-sd(timeLast$x[!is.na(demo$n_e2[indexDemo])], na.rm = TRUE)
tabCogGenetics.apoe[each,c(4,5)]<-summary(model)$coefficients[2,c(1,4)]
tabCogGenetics.apoe[each,c(7,8)]<-summary(model)$coefficients[3,c(1,4)]
sumSq<-anova(model)["Sum Sq"]
tabCogGenetics.apoe[each,c(6,9)]<-sumSq[1:2,1]/sum(sumSq)*100


##test apoe as factor
model<-lm(scoreLast$x ~ demo$APOE_geno[indexDemo] + timeLast$x + demo$Age[indexDemo] + demo$Sex[indexDemo] + demo$BDR_Centre_key[indexDemo] + genoPCs$C1[indexGenetic] + genoPCs$C2[indexGenetic] + genoPCs$C3[indexGenetic] + genoPCs$C4[indexGenetic] + genoPCs$C5[indexGenetic] + genoPCs$C6[indexGenetic] + genoPCs$C7[indexGenetic] + genoPCs$C8[indexGenetic])
null<-lm(scoreLast$x ~ timeLast$x + demo$Age[indexDemo] + demo$Sex[indexDemo]+ demo$BDR_Centre_key[indexDemo] + genoPCs$C1[indexGenetic] + genoPCs$C2[indexGenetic] + genoPCs$C3[indexGenetic] + genoPCs$C4[indexGenetic] + genoPCs$C5[indexGenetic] + genoPCs$C6[indexGenetic] + genoPCs$C7[indexGenetic] + genoPCs$C8[indexGenetic], subset = !is.na(demo$APOE_geno[indexDemo]))
	tabCogGenetics.apoe2[each,1]<-length(residuals(model))
	tabCogGenetics.apoe2[each,2]<-mean(timeLast$x[!is.na(demo$n_e2[indexDemo])], na.rm = TRUE)
	tabCogGenetics.apoe2[each,3]<-sd(timeLast$x[!is.na(demo$n_e2[indexDemo])], na.rm = TRUE)
	
	tabCogGenetics.apoe2[each,4]<-anova(model, null)[2,6]
	sumSq<-anova(model)["Sum Sq"]
	tabCogGenetics.apoe2[each,5]<-sumSq[1,1]/sum(sumSq)*100
	tabCogGenetics.apoe2[each,c(6,8,10,12,14)]<-summary(model)$coefficients[c("demo$APOE_geno[indexDemo]e2/e2","demo$APOE_geno[indexDemo]e2/e3","demo$APOE_geno[indexDemo]e2/e4","demo$APOE_geno[indexDemo]e3/e4","demo$APOE_geno[indexDemo]e4/e4"),4]
	tabCogGenetics.apoe2[each,c(7,9,11,13,15)]<-summary(model)$coefficients[c("demo$APOE_geno[indexDemo]e2/e2","demo$APOE_geno[indexDemo]e2/e3","demo$APOE_geno[indexDemo]e2/e4","demo$APOE_geno[indexDemo]e3/e4","demo$APOE_geno[indexDemo]e4/e4"),1]


##test prs
model<-lm(scoreLast$x ~ prs.scale[indexGenetic] + timeLast$x + demo$Age[indexDemo] + demo$Sex[indexDemo]+ demo$BDR_Centre_key[indexDemo] + genoPCs$C1[indexGenetic] + genoPCs$C2[indexGenetic] + genoPCs$C3[indexGenetic] + genoPCs$C4[indexGenetic] + genoPCs$C5[indexGenetic] + genoPCs$C6[indexGenetic] + genoPCs$C7[indexGenetic] + genoPCs$C8[indexGenetic])
sumSq<-anova(model)["Sum Sq"]
tabCogGenetics.prs[each,1]<-length(residuals(model))
tabCogGenetics.prs[each,2]<-mean(timeLast$x[!is.na(prs.scale[indexGenetic])], na.rm = TRUE)
tabCogGenetics.prs[each,3]<-sd(timeLast$x[!is.na(prs.scale[indexGenetic])], na.rm = TRUE)
tabCogGenetics.prs[each,6]<-sumSq[1,1]/sum(sumSq)*100
tabCogGenetics.prs[each,c(4,5)]<-summary(model)$coefficients[2,c(1,4)]

## test prs and apoe together
model<-lm(scoreLast$x ~ prs.scale[indexGenetic] + demo$n_e2[indexDemo] + demo$n_e4[indexDemo] + timeLast$x + demo$Age[indexDemo] + demo$Sex[indexDemo]+ demo$BDR_Centre_key[indexDemo] + genoPCs$C1[indexGenetic] + genoPCs$C2[indexGenetic] + genoPCs$C3[indexGenetic] + genoPCs$C4[indexGenetic] + genoPCs$C5[indexGenetic] + genoPCs$C6[indexGenetic] + genoPCs$C7[indexGenetic] + genoPCs$C8[indexGenetic])
tabCogGenetics.joint[each,1]<-length(residuals(model))
tabCogGenetics.joint[each,2]<-mean(timeLast$x[!is.na(prs.scale[indexGenetic]) & !is.na(demo$n_e2[indexDemo])], na.rm = TRUE)
tabCogGenetics.joint[each,3]<-sd(timeLast$x[!is.na(prs.scale[indexGenetic]) & !is.na(demo$n_e2[indexDemo])], na.rm = TRUE)
tabCogGenetics.joint[each,c(4,5)]<-summary(model)$coefficients[2,c(1,4)]
tabCogGenetics.joint[each,c(7,8)]<-summary(model)$coefficients[3,c(1,4)]
tabCogGenetics.joint[each,c(10,11)]<-summary(model)$coefficients[4,c(1,4)]
sumSq<-anova(model)["Sum Sq"]
tabCogGenetics.joint[each,c(6,9,12)]<-sumSq[1:3,1]/sum(sumSq)*100

## test interaction
model<-lm(scoreLast$x ~ demo$n_e2[indexDemo] + demo$n_e4[indexDemo] + prs.scale[indexGenetic] + prs.scale[indexGenetic] * demo$n_e2[indexDemo] +  prs.scale[indexGenetic] * demo$n_e4[indexDemo] + timeLast$x + demo$Age[indexDemo] + demo$Sex[indexDemo]+ demo$BDR_Centre_key[indexDemo] + genoPCs$C1[indexGenetic] + genoPCs$C2[indexGenetic] + genoPCs$C3[indexGenetic] + genoPCs$C4[indexGenetic] + genoPCs$C5[indexGenetic] + genoPCs$C6[indexGenetic] + genoPCs$C7[indexGenetic] + genoPCs$C8[indexGenetic])

tabCogGenetics.int[each,1]<-length(residuals(model))
tabCogGenetics.int[each,2]<-mean(timeLast$x[!is.na(prs.scale[indexGenetic]) & !is.na(demo$n_e2[indexDemo])], na.rm = TRUE)
tabCogGenetics.int[each,3]<-sd(timeLast$x[!is.na(prs.scale[indexGenetic]) & !is.na(demo$n_e2[indexDemo])], na.rm = TRUE)
tabCogGenetics.int[each,c(5,4)]<-summary(model)$coefficients[2,c(1,4)]
tabCogGenetics.int[each,c(8,7)]<-summary(model)$coefficients[3,c(1,4)]
tabCogGenetics.int[each,c(11,10)]<-summary(model)$coefficients[4,c(1,4)]
tabCogGenetics.int[each,c(13)]<-summary(model)$coefficients["demo$n_e2[indexDemo]:prs.scale[indexGenetic]",c(4)]
tabCogGenetics.int[each,c(15)]<-summary(model)$coefficients["demo$n_e4[indexDemo]:prs.scale[indexGenetic]",c(4)]
sumSq<-anova(model)["Sum Sq"]
tabCogGenetics.int[each,c(6,9,12)]<-sumSq[1:3,1]/sum(sumSq)*100
tabCogGenetics.int[each,c(14)]<-sumSq["demo$n_e2[indexDemo]:prs.scale[indexGenetic]",1]/sum(sumSq)*100
tabCogGenetics.int[each,c(16)]<-sumSq["demo$n_e4[indexDemo]:prs.scale[indexGenetic]",1]/sum(sumSq)*100

}

rownames(sum.demo)<-c("Sex", "Age", cogVar)
colnames(sum.demo)<-c("%", "Mean", "SD", "N")
write.csv(signif(sum.demo, 3), "/mnt/data1/Gemma/BDR/BDR_genetic/Analysis/CohortSummaryCognitiveStatusAtDeath.csv")


write.csv(tabCogGenetics.apoe, "/mnt/data1/Gemma/BDR/BDR_genetic/Analysis/APOEAnalysisCognitiveStatusAtDeathNumeric.csv")
write.csv(tabCogGenetics.apoe2, "/mnt/data1/Gemma/BDR/BDR_genetic/Analysis/APOEAnalysisCognitiveStatusAtDeathGenotype.csv")
write.csv(tabCogGenetics.prs, "/mnt/data1/Gemma/BDR/BDR_genetic/Analysis/PRSAnalysisCognitiveStatusAtDeath.csv")
write.csv(tabCogGenetics.joint, "/mnt/data1/Gemma/BDR/BDR_genetic/Analysis/JointAnalysisCognitiveStatusAtDeathNumeric.csv")
 
write.csv(tabCogGenetics.int, "/mnt/data1/Gemma/BDR/BDR_genetic/Analysis/IntAnalysisCognitiveStatusAtDeathNumeric.csv")

## test controlling for neuropathology
tabCogGenetics.apoeAdj<-matrix(data = NA, ncol = 9, nrow = length(cogVar))
colnames(tabCogGenetics.apoeAdj)<-c("N", "MeanTimeSinceVisit", "SETimeSinceVisit", "APOE:E2:Coeff", "APOE:E2:P", "APOE:E2:VarExp", "APOE:E4:Coeff", "APOE:E4:P", "APOE:E4:VarExp")
rownames(tabCogGenetics.apoeAdj)<-cogVar

tabBraak.apoeAdj<-matrix(data = NA, ncol = 9, nrow = length(cogVar))
colnames(tabBraak.apoeAdj)<-c("N", "MeanTimeSinceVisit", "SETimeSinceVisit", "APOE:E2:Coeff", "APOE:E2:P", "APOE:E2:VarExp", "APOE:E4:Coeff", "APOE:E4:P", "APOE:E4:VarExp")
rownames(tabBraak.apoeAdj)<-cogVar


tabCogGenetics.prsAdj<-matrix(data = NA, ncol = 6, nrow = length(cogVar))
colnames(tabCogGenetics.prsAdj)<-c("N", "MeanTimeSinceVisit", "SETimeSinceVisit", "PRS:Coeff", "PRS:P", "PRS:VarExp")
rownames(tabCogGenetics.prsAdj)<-cogVar

tabCogGenetics.jointAdj<-matrix(data = NA, ncol = 12, nrow = length(cogVar))
colnames(tabCogGenetics.jointAdj)<-c("N", "MeanTimeSinceVisit", "SETimeSinceVisit", "PRS:Coeff", "PRS:P", "PRS:VarExp","APOE:E2:Coeff", "APOE:E2:P", "APOE:E2:VarExp", "APOE:E4:Coeff", "APOE:E4:P", "APOE:E4:VarExp")
rownames(tabCogGenetics.jointAdj)<-cogVar


for(each in cogVar){


cad.tmp<-cad[!is.na(cad[,each]),]
nIds<-length(unique(cad.tmp$BBNId))
nVisitCog<-table(cad.tmp$BBNId)
scoreLast<-aggregate(cad.tmp[,each], by = list(cad.tmp$BBNId), FUN = tail, n = 1)
timeLast<-aggregate(cad.tmp$TimePriorDeath, by = list(cad.tmp$BBNId), FUN = tail, n = 1)
indexDemo<-match(scoreLast$Group.1, demo$BBNId)
dnaIds<-demo$DNA_ID[indexDemo]
indexGenetic<-match(dnaIds, prs.imp.noapoe$FID)

##test apoe as number of e2 and number of e4/e4
model<-lm(scoreLast$x ~ demo$n_e2[indexDemo] + demo$n_e4[indexDemo] + demo$Braak_tangle[indexDemo] + timeLast$x + demo$Age[indexDemo] + demo$Sex[indexDemo]+ demo$BDR_Centre_key[indexDemo] + genoPCs$C1[indexGenetic] + genoPCs$C2[indexGenetic] + genoPCs$C3[indexGenetic] + genoPCs$C4[indexGenetic] + genoPCs$C5[indexGenetic] + genoPCs$C6[indexGenetic] + genoPCs$C7[indexGenetic] + genoPCs$C8[indexGenetic])
tabCogGenetics.apoeAdj[each,1]<-length(residuals(model))
tabCogGenetics.apoeAdj[each,2]<-mean(timeLast$x[!is.na(demo$n_e2[indexDemo])], na.rm = TRUE)
tabCogGenetics.apoeAdj[each,3]<-sd(timeLast$x[!is.na(demo$n_e2[indexDemo])], na.rm = TRUE)
tabCogGenetics.apoeAdj[each,c(4,5)]<-summary(model)$coefficients[2,c(1,4)]
tabCogGenetics.apoeAdj[each,c(7,8)]<-summary(model)$coefficients[3,c(1,4)]
sumSq<-anova(model)["Sum Sq"]
tabCogGenetics.apoeAdj[each,c(6,9)]<-sumSq[1:2,1]/sum(sumSq)*100

## test Braak tangle stage with cognition as a confounder
model<-lm(demo$Braak_tangle[indexDemo] ~ demo$n_e2[indexDemo] + demo$n_e4[indexDemo] + scoreLast$x  + timeLast$x + demo$Age[indexDemo] + demo$Sex[indexDemo]+ demo$BDR_Centre_key[indexDemo] + genoPCs$C1[indexGenetic] + genoPCs$C2[indexGenetic] + genoPCs$C3[indexGenetic] + genoPCs$C4[indexGenetic] + genoPCs$C5[indexGenetic] + genoPCs$C6[indexGenetic] + genoPCs$C7[indexGenetic] + genoPCs$C8[indexGenetic])
tabBraak.apoeAdj[each,1]<-length(residuals(model))
tabBraak.apoeAdj[each,2]<-mean(timeLast$x[!is.na(demo$n_e2[indexDemo])], na.rm = TRUE)
tabBraak.apoeAdj[each,3]<-sd(timeLast$x[!is.na(demo$n_e2[indexDemo])], na.rm = TRUE)
tabBraak.apoeAdj[each,c(4,5)]<-summary(model)$coefficients[2,c(1,4)]
tabBraak.apoeAdj[each,c(7,8)]<-summary(model)$coefficients[3,c(1,4)]
sumSq<-anova(model)["Sum Sq"]
tabBraak.apoeAdj[each,c(6,9)]<-sumSq[1:2,1]/sum(sumSq)*100


##test prs
model<-lm(scoreLast$x ~ prs.scale[indexGenetic] + demo$Braak_tangle[indexDemo] + timeLast$x + demo$Age[indexDemo] + demo$Sex[indexDemo]+ demo$BDR_Centre_key[indexDemo] + genoPCs$C1[indexGenetic] + genoPCs$C2[indexGenetic] + genoPCs$C3[indexGenetic] + genoPCs$C4[indexGenetic] + genoPCs$C5[indexGenetic] + genoPCs$C6[indexGenetic] + genoPCs$C7[indexGenetic] + genoPCs$C8[indexGenetic])
sumSq<-anova(model)["Sum Sq"]
tabCogGenetics.prsAdj[each,1]<-length(residuals(model))
tabCogGenetics.prsAdj[each,2]<-mean(timeLast$x[!is.na(prs.scale[indexGenetic])], na.rm = TRUE)
tabCogGenetics.prsAdj[each,3]<-sd(timeLast$x[!is.na(prs.scale[indexGenetic])], na.rm = TRUE)
tabCogGenetics.prsAdj[each,6]<-sumSq[1,1]/sum(sumSq)*100
tabCogGenetics.prsAdj[each,c(4,5)]<-summary(model)$coefficients[2,c(1,4)]

## test prs and apoe together
model<-lm(scoreLast$x ~ prs.scale[indexGenetic] + demo$n_e2[indexDemo] + demo$n_e4[indexDemo] + demo$Braak_tangle[indexDemo] + timeLast$x + demo$Age[indexDemo] + demo$Sex[indexDemo]+ demo$BDR_Centre_key[indexDemo] + genoPCs$C1[indexGenetic] + genoPCs$C2[indexGenetic] + genoPCs$C3[indexGenetic] + genoPCs$C4[indexGenetic] + genoPCs$C5[indexGenetic] + genoPCs$C6[indexGenetic] + genoPCs$C7[indexGenetic] + genoPCs$C8[indexGenetic])
tabCogGenetics.jointAdj[each,1]<-length(residuals(model))
tabCogGenetics.jointAdj[each,2]<-mean(timeLast$x[!is.na(prs.scale[indexGenetic]) & !is.na(demo$n_e2[indexDemo])], na.rm = TRUE)
tabCogGenetics.jointAdj[each,3]<-sd(timeLast$x[!is.na(prs.scale[indexGenetic]) & !is.na(demo$n_e2[indexDemo])], na.rm = TRUE)
tabCogGenetics.jointAdj[each,c(4,5)]<-summary(model)$coefficients[2,c(1,4)]
tabCogGenetics.jointAdj[each,c(7,8)]<-summary(model)$coefficients[3,c(1,4)]
tabCogGenetics.jointAdj[each,c(10,11)]<-summary(model)$coefficients[4,c(1,4)]
sumSq<-anova(model)["Sum Sq"]
tabCogGenetics.jointAdj[each,c(6,9,12)]<-sumSq[1:3,1]/sum(sumSq)*100


}

write.csv(tabCogGenetics.apoeAdj, "/mnt/data1/Gemma/BDR/BDR_genetic/Analysis/APOEAnalysisCognitiveStatusAtDeathNumericControlBraak.csv")
write.csv(tabCogGenetics.prsAdj, "/mnt/data1/Gemma/BDR/BDR_genetic/Analysis/PRSAnalysisCognitiveStatusAtDeathControlBraak.csv")
write.csv(tabCogGenetics.jointAdj, "/mnt/data1/Gemma/BDR/BDR_genetic/Analysis/JointAnalysisCognitiveStatusAtDeathNumericControlBraak.csv")
 
write.csv(tabBraak.apoeAdj, "/mnt/data1/Gemma/BDR/BDR_genetic/Analysis/APOEAnalysisBraakNumericControlCognitiveStatusAtDeath.csv")



