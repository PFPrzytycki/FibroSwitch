##SEs per cell-types per condition:
for(cType in cTypes){
	for(cond in c("Sham","TAC")){
		if(cond=="Sham"){samps = paste0("Sample",1:2)}
		if(cond=="TAC"){samps = paste0("Sample",3:4)}
		allPeaks = GRanges()
		for(samp in samps){
		    mat = get(paste0("atac_",samp))
		    fCells = get(paste0(cType,"Cells_",samp))
		    allPeaks = c(allPeaks, rownames(mat[,fCells]))
		}
		allPeaks = reduce(allPeaks)
		allPeaks = sort(allPeaks)
		####### V1, completely new
		newSEs = GRanges()
		thisSE = allPeaks[1]
		for(i in 2:(length(allPeaks))){
			if(!is.na(distance(thisSE, allPeaks[i])) & distance(thisSE, allPeaks[i])<12500) {
				thisSE = punion(thisSE, allPeaks[i], fill.gap=TRUE)
			}
			else{
				newSEs = c(newSEs, thisSE)
				thisSE = allPeaks[i]
			}
		}
		newEnhancers = data.frame(newSEs[width(newSEs)>500])
		write.table(newEnhancers, paste0("newEnhancers_",cType,cond,".txt"), quote = FALSE, row.names = FALSE)
	}
}

#select super enhancers in each condition
for(cType in cTypes){
    for(cond in c("Sham","TAC")){
        tissue = paste0(cType,cond)
        allEnhancerAcc = read.table(paste0("acc_frame_",tissue,"_E.txt"), header=TRUE)
        SEcut = calculate_cutoff(allEnhancerAcc$acc, plotName = paste0(cType,cond,".pdf")) #ROSE
        write.table(allEnhancerAcc[allEnhancerAcc$acc>SEcut$absolute,], paste0("acc_frame_",tissue,"_SE.txt"), row.names=FALSE, quote=FALSE, sep="\t")
    }
}
