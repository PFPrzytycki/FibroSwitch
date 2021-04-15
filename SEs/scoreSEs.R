superEnhancerScores = data.frame(chr=allSuperEnhancers$V2, start=allSuperEnhancers$V3, end=allSuperEnhancers$V4)
for(cType in c("Fibroblast", "Endothelial", "Myeloid")){
scoreSuperEnhancersEF_LM_pval = c()
scoreSuperEnhancersEF_LM = c()
for(i in 1:length(SuperEnhancersMM10)){
    if(! (allSuperEnhancers$V1[i] %in% normPeakFrame$Enhancer[normPeakFrame$Condition=="Sham" & normPeakFrame$Type==cType] & allSuperEnhancers$V1[i] %in% normPeakFrame$Enhancer[normPeakFrame$Condition=="TAC" & normPeakFrame$Type==cType] & allSuperEnhancers$V1[i] %in% normPeakFrame$Enhancer[normPeakFrame$Condition=="JQ1" & normPeakFrame$Type==cType] & allSuperEnhancers$V1[i] %in% normPeakFrame$Enhancer[normPeakFrame$Condition=="JQ1 withdrawn" & normPeakFrame$Type==cType])){
    		scoreSuperEnhancersEF_LM_pval = c(scoreSuperEnhancersEF_LM_pval, 1)
    		scoreSuperEnhancersEF_LM = c(scoreSuperEnhancersEF_LM, 0)
    	}
    else{
        meanScoresAll = normPeakFrame[normPeakFrame$Enhancer==allSuperEnhancers$V1[i],] %>% group_by(Condition) %>% summarise(condScore=mean(normScore),condsd=sd(normScore))
        meanScores = normPeakFrame[normPeakFrame$Enhancer==allSuperEnhancers$V1[i] & normPeakFrame$Type==cType,] %>% group_by(Condition) %>% summarise(condScore=mean(normScore))
        

        meanScoresAll = normPeakFrame[normPeakFrame$Enhancer==allSuperEnhancers$V1[i] & normPeakFrame$Type!=cType,] %>% group_by(Condition) %>% summarise(condScore=mean(normScore))
        fitFrame = data.frame(EF=EF$V2, Acc=rep(meanScoresAll$condScore, times=c(4,6,4,6)), Acc_F=rep(meanScores$condScore, times=c(4,6,4,6)))
        fit = lm(EF ~ Acc+Acc_F, data=fitFrame)
        scoreSuperEnhancersEF_LM_pval = c(scoreSuperEnhancersEF_LM_pval, summary(fit)$coefficients["Acc_F",4])
        scoreSuperEnhancersEF_LM = c(scoreSuperEnhancersEF_LM, summary(fit)$coefficients["Acc_F",3])
    }

}
superEnhancerScores = cbind(superEnhancerScores, scoreSuperEnhancersEF_LM_pval, scoreSuperEnhancersEF_LM)
}
colnames(superEnhancerScores) = c("chr", "start", "end", "Fibroblast_pval", "Fibroblast_coef", "Endothelial_pval", "Endothelial_coef", "Myeloid_pval", "Myeloid_coef")

write.table(superEnhancerScores, "superEnhancerScores.txt", quote = FALSE, row.names = FALSE, sep="\t")
