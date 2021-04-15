#Compare peaks in two samples in a fixed region
normPeakFrame = data.frame(Condition=character(),Type=character(),normScore=numeric())
for(samp in paste0("Sample",c(1:4,7,8,5,6))) {
    mat = as.matrix(get(paste0("atac_",samp)))
    for(cType in c("Fibroblast", "Endothelial", "Myeloid")){
	    fCells = get(paste0(cType,"Cells_",samp))
	    normPeaks = colSums(mat[countOverlaps(as(rownames(mat), "GRanges"),as("chr11:101805576-101876751", "GRanges"))>0,fCells]>0) / colSums(mat[,fCells]>0)
	    if(samp == "Sample1" | samp == "Sample2"){cond = "Sham"}
	    if(samp == "Sample3" | samp == "Sample4"){cond = "TAC"}
	    if(samp == "Sample7" | samp == "Sample8"){cond = "JQ1"}
	    if(samp == "Sample5" | samp == "Sample6"){cond = "JQ1 withdrawn"}
	    normPeakFrame = rbind(normPeakFrame, data.frame(Condition=cond, Type=cType,normScore=normPeaks))
	    }
}
#Genome-wide normpeak frame with pseudoRep
normPeakFrame = data.frame(Condition=character(),Type=character(),normScore=numeric())
for(cType in cTypes){
    fPeaksTest = fread(paste0("all",cType,"Peaks_intersect.bed"))
    fPeaksTest = reduce(GRanges(fPeaksTest$V1, IRanges(fPeaksTest$V2, fPeaksTest$V3)))
    for(cond in c("Sham","TAC","JQ1","JQ1 Withdrawn")){
        if(cond=="Sham"){samps = paste0("Sample",1:2)}
        if(cond=="TAC"){samps = paste0("Sample",3:4)}
        if(cond=="JQ1"){samps = paste0("Sample",7:8)}
        if(cond=="JQ1 Withdrawn"){samps = paste0("Sample",5:6)}
        for(samp in samps){
            mat = as.matrix(get(paste0("atac_",samp)))
            fCells = get(paste0(cType,"Cells_",samp))
            normPeaks = colSums(mat[countOverlaps(as(rownames(mat), "GRanges"),fPeaksTest)>0,fCells]>0) / colSums(mat[,fCells]>0)
            normPeakFrame = rbind(normPeakFrame, data.frame(Condition=cond, Type=cType,normScore=normPeaks))
        }
    }
}
ggplot(normPeakFrame) + geom_boxplot(aes(Condition,normScore, color=Type)) + facet_wrap(~Type)
#every individual peak in meox1 SE:
cType = "Fibroblast"
fPeaksTest = fread(paste0("all",cType,"Peaks.bed"))
fPeaksTest = reduce(GRanges(fPeaksTest$V1, IRanges(fPeaksTest$V2, fPeaksTest$V3)))
fPeaksTest = fPeaksTest[countOverlaps(fPeaksTest,as("chr11:101805576-101876751", "GRanges"))>0]
normPeakFrame = data.frame(Condition=character(),Type=character(),normScore=numeric(), PeakNumber=numeric())
for(peak in 15:1){
    thisPeak = fPeaksTest[peak]
    for(samp in paste0("Sample",c(1:4,7,8,5,6))) {
        mat = as.matrix(get(paste0("atac_",samp)))
        for(cType in c("Fibroblast")){
            fCells = get(paste0(cType,"Cells_",samp))
            whichPeaks = countOverlaps(as(rownames(mat), "GRanges"),thisPeak)>0
            if(sum(whichPeaks)==0){
                normPeaks = rep(0, length(fCells))
            }
            else if(sum(whichPeaks)==1){
                normPeaks = as.numeric(mat[whichPeaks,fCells]>0) / colSums(mat[,fCells]>0)
            }
            else{
                normPeaks = colSums(mat[,fCells]>0) / colSums(mat[,fCells]>0)
            }
            if(samp == "Sample1" | samp == "Sample2"){cond = "Sham"}
            if(samp == "Sample3" | samp == "Sample4"){cond = "TAC"}
            if(samp == "Sample7" | samp == "Sample8"){cond = "JQ1"}
            if(samp == "Sample5" | samp == "Sample6"){cond = "JQ1 withdrawn"}
            normPeakFrame = rbind(normPeakFrame, data.frame(Condition=cond, Type=cType,normScore=normPeaks, PeakNumber=as.character(thisPeak)))
        }
    }
}
ggplot(normPeakFrame) + geom_boxplot(aes(Condition,normScore)) + facet_wrap(~PeakNumber, ncol = 3)
ggplot(normPeakFrame) + stat_summary(aes(Condition,normScore), fun.y=mean, geom="point", shape=20, size=1, color="red", fill="red") + facet_wrap(~PeakNumber, ncol = 3, scales="free_y")
ggplot(normPeakFrame) + geom_smooth(aes(as.numeric(Condition),normScore), method="loess", span=1) + facet_wrap(~PeakNumber, ncol = 3, scales="free_y")

conds = c("Sham","TAC","JQ1","JQ1 withdrawn")
pvalFrame = data.frame()
for(peak in unique(normPeakFrame$PeakNumber)){
    for(cond1 in conds[1:3]){
        for(cond2 in conds[(which(conds==cond1)+1):4]){
            scores1 = normPeakFrame$normScore[normPeakFrame$PeakNumber==peak & normPeakFrame$Condition==cond1]
            scores2 = normPeakFrame$normScore[normPeakFrame$PeakNumber==peak & normPeakFrame$Condition==cond2]
            if(length(scores1)==0 | length(scores2)==0){next}
            pvalFrame = rbind(pvalFrame, data.frame(peak,cond1,cond2,pval=wilcox.test(scores1, scores2)$p.value))
        }
    }
}
pvalFrame$pval[is.na(pvalFrame$pval)] = 1
