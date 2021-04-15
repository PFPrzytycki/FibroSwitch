#Cell-type specific peaks for sample 2:
#Fib Cells
assign(paste0("fibCells_",samp), clusters_split==7)
fib_sums = colSums(as.matrix(atac_Sample2[,fibCells_Sample2])>0)
nonFib_sums = colSums(as.matrix(atac_Sample2[,!fibCells_Sample2])>0)

m2 = ((as.matrix(atac_Sample2)[,fibCells_Sample2])>0)

select_nonFib_cells = sapply(fib_sums, function(x) sample(which(nonFib_sums<=x+1500 & nonFib_sums>=x-1500), 1))
m1 = ((as.matrix(atac_Sample2)[,!fibCells_Sample2][,select_nonFib_cells])>0)
peakScores1 = sapply(1:length(rownames(atac_Sample2)), function(x) wilcox.test(as.numeric(m1[x,]) , as.numeric(m2[x,]), alternative = "less")$p.value)
select_nonFib_cells = sapply(fib_sums, function(x) sample(which(nonFib_sums<=x+1500 & nonFib_sums>=x-1500), 1))
m1 = ((as.matrix(atac_Sample2)[,!fibCells_Sample2][,select_nonFib_cells])>0)
peakScores2 = sapply(1:length(rownames(atac_Sample2)), function(x) wilcox.test(as.numeric(m1[x,]) , as.numeric(m2[x,]), alternative = "less")$p.value)
select_nonFib_cells = sapply(fib_sums, function(x) sample(which(nonFib_sums<=x+1500 & nonFib_sums>=x-1500), 1))
m1 = ((as.matrix(atac_Sample2)[,!fibCells_Sample2][,select_nonFib_cells])>0)
peakScores3 = sapply(1:length(rownames(atac_Sample2)), function(x) wilcox.test(as.numeric(m1[x,]) , as.numeric(m2[x,]), alternative = "less")$p.value)

combinedPeakScores_Fib = p.adjust(sapply(1:length(peakScores1), function(x) sumlog(c(peakScores1[x],peakScores2[x],peakScores3[x]))[[3]]))
# write.table(data.frame(as(rownames(atac_Sample2), "GRanges"))[combinedPeakScores_Fib<.1,1:3], "fibro_sham.bed", quote = FALSE, col.names = FALSE, row.names = FALSE)


#Endothelial Cells
assign(paste0("endCells_",samp), clusters_split==3)
end_sums = colSums(as.matrix(atac_Sample2[,endCells_Sample2])>0)
nonEnd_sums = colSums(as.matrix(atac_Sample2[,!endCells_Sample2])>0)

m2 = ((as.matrix(atac_Sample2)[,endCells_Sample2])>0)

select_nonEnd_cells = sapply(end_sums, function(x) sample(which(nonEnd_sums<=x+1500 & nonEnd_sums>=x-1500), 1))
m1 = ((as.matrix(atac_Sample2)[,!endCells_Sample2][,select_nonEnd_cells])>0)
peakScores1 = sapply(1:length(rownames(atac_Sample2)), function(x) wilcox.test(as.numeric(m1[x,]) , as.numeric(m2[x,]), alternative = "less")$p.value)
select_nonEnd_cells = sapply(end_sums, function(x) sample(which(nonEnd_sums<=x+1500 & nonEnd_sums>=x-1500), 1))
m1 = ((as.matrix(atac_Sample2)[,!endCells_Sample2][,select_nonEnd_cells])>0)
peakScores2 = sapply(1:length(rownames(atac_Sample2)), function(x) wilcox.test(as.numeric(m1[x,]) , as.numeric(m2[x,]), alternative = "less")$p.value)
select_nonEnd_cells = sapply(end_sums, function(x) sample(which(nonEnd_sums<=x+1500 & nonEnd_sums>=x-1500), 1))
m1 = ((as.matrix(atac_Sample2)[,!endCells_Sample2][,select_nonEnd_cells])>0)
peakScores3 = sapply(1:length(rownames(atac_Sample2)), function(x) wilcox.test(as.numeric(m1[x,]) , as.numeric(m2[x,]), alternative = "less")$p.value)

combinedPeakScores_End = p.adjust(sapply(1:length(peakScores1), function(x) sumlog(c(peakScores1[x],peakScores2[x],peakScores3[x]))[[3]]))

#Myeloid Cells
assign(paste0("myeCells_",samp), clusters_split==8)
mye_sums = colSums(as.matrix(atac_Sample2[,myeCells_Sample2])>0)
nonMye_sums = colSums(as.matrix(atac_Sample2[,!myeCells_Sample2])>0)

m2 = ((as.matrix(atac_Sample2)[,myeCells_Sample2])>0)

select_nonMye_cells = sapply(mye_sums, function(x) sample(which(nonMye_sums<=x+1500 & nonMye_sums>=x-1500), 1))
m1 = ((as.matrix(atac_Sample2)[,!myeCells_Sample2][,select_nonMye_cells])>0)
peakScores1 = sapply(1:length(rownames(atac_Sample2)), function(x) wilcox.test(as.numeric(m1[x,]) , as.numeric(m2[x,]), alternative = "less")$p.value)
select_nonMye_cells = sapply(mye_sums, function(x) sample(which(nonMye_sums<=x+1500 & nonMye_sums>=x-1500), 1))
m1 = ((as.matrix(atac_Sample2)[,!myeCells_Sample2][,select_nonMye_cells])>0)
peakScores2 = sapply(1:length(rownames(atac_Sample2)), function(x) wilcox.test(as.numeric(m1[x,]) , as.numeric(m2[x,]), alternative = "less")$p.value)
select_nonMye_cells = sapply(mye_sums, function(x) sample(which(nonMye_sums<=x+1500 & nonMye_sums>=x-1500), 1))
m1 = ((as.matrix(atac_Sample2)[,!myeCells_Sample2][,select_nonMye_cells])>0)
peakScores3 = sapply(1:length(rownames(atac_Sample2)), function(x) wilcox.test(as.numeric(m1[x,]) , as.numeric(m2[x,]), alternative = "less")$p.value)

combinedPeakScores_Mye = p.adjust(sapply(1:length(peakScores1), function(x) sumlog(c(peakScores1[x],peakScores2[x],peakScores3[x]))[[3]]))

write.table(c("track name=Fibroblast color=0,225,0"), "sham.bed", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(data.frame(as(rownames(atac_Sample2), "GRanges"))[combinedPeakScores_Fib<.1,1:3], "sham.bed", quote = FALSE, col.names = FALSE, row.names = FALSE, append = TRUE)
write.table(c("track name=Myeloid color=0,0,225"), "sham.bed", quote = FALSE, col.names = FALSE, row.names = FALSE, append = TRUE)
write.table(data.frame(as(rownames(atac_Sample2), "GRanges"))[combinedPeakScores_Mye<.1,1:3], "sham.bed", quote = FALSE, col.names = FALSE, row.names = FALSE, append = TRUE)
write.table(c("track name=Endothelial color=225,0,0"), "sham.bed", quote = FALSE, col.names = FALSE, row.names = FALSE, append = TRUE)
write.table(data.frame(as(rownames(atac_Sample2), "GRanges"))[combinedPeakScores_End<.1,1:3], "sham.bed", quote = FALSE, col.names = FALSE, row.names = FALSE, append = TRUE)



#Cell-type specific peaks for sample 4:
#Fib Cells
assign(paste0("fibCells_",samp), clusters_split==4)
fib_sums = colSums(as.matrix(atac_Sample4[,fibCells_Sample4])>0)
nonFib_sums = colSums(as.matrix(atac_Sample4[,!fibCells_Sample4])>0)

m2 = ((as.matrix(atac_Sample4)[,fibCells_Sample4])>0)

select_nonFib_cells = sapply(fib_sums, function(x) sample(which(nonFib_sums<=x+1500 & nonFib_sums>=x-1500), 1))
m1 = ((as.matrix(atac_Sample4)[,!fibCells_Sample4][,select_nonFib_cells])>0)
peakScores1 = sapply(1:length(rownames(atac_Sample4)), function(x) wilcox.test(as.numeric(m1[x,]) , as.numeric(m2[x,]), alternative = "less")$p.value)
select_nonFib_cells = sapply(fib_sums, function(x) sample(which(nonFib_sums<=x+1500 & nonFib_sums>=x-1500), 1))
m1 = ((as.matrix(atac_Sample4)[,!fibCells_Sample4][,select_nonFib_cells])>0)
peakScores2 = sapply(1:length(rownames(atac_Sample4)), function(x) wilcox.test(as.numeric(m1[x,]) , as.numeric(m2[x,]), alternative = "less")$p.value)
select_nonFib_cells = sapply(fib_sums, function(x) sample(which(nonFib_sums<=x+1500 & nonFib_sums>=x-1500), 1))
m1 = ((as.matrix(atac_Sample4)[,!fibCells_Sample4][,select_nonFib_cells])>0)
peakScores3 = sapply(1:length(rownames(atac_Sample4)), function(x) wilcox.test(as.numeric(m1[x,]) , as.numeric(m2[x,]), alternative = "less")$p.value)

combinedPeakScores_Fib = p.adjust(sapply(1:length(peakScores1), function(x) sumlog(c(peakScores1[x],peakScores2[x],peakScores3[x]))[[3]]))
# write.table(data.frame(as(rownames(atac_Sample4), "GRanges"))[combinedPeakScores_Fib<.1,1:3], "fibro_sham.bed", quote = FALSE, col.names = FALSE, row.names = FALSE)


#Endothelial Cells
assign(paste0("endCells_",samp), clusters_split==5)
end_sums = colSums(as.matrix(atac_Sample4[,endCells_Sample4])>0)
nonEnd_sums = colSums(as.matrix(atac_Sample4[,!endCells_Sample4])>0)

m2 = ((as.matrix(atac_Sample4)[,endCells_Sample4])>0)

select_nonEnd_cells = sapply(end_sums, function(x) {
    matchSet = which(nonEnd_sums<=x+1000 & nonEnd_sums>=x-1000)
    if(length(matchSet)==0){matchSet = order(abs(nonEnd_sums-x))[1]}
    sample(matchSet, 1)
    })
m1 = ((as.matrix(atac_Sample4)[,!endCells_Sample4][,select_nonEnd_cells])>0)
peakScores1 = sapply(1:length(rownames(atac_Sample4)), function(x) wilcox.test(as.numeric(m1[x,]) , as.numeric(m2[x,]), alternative = "less")$p.value)
select_nonEnd_cells = sapply(end_sums, function(x) {
    matchSet = which(nonEnd_sums<=x+1000 & nonEnd_sums>=x-1000)
    if(length(matchSet)==0){matchSet = order(abs(nonEnd_sums-x))[1]}
    sample(matchSet, 1)
    })
m1 = ((as.matrix(atac_Sample4)[,!endCells_Sample4][,select_nonEnd_cells])>0)
peakScores2 = sapply(1:length(rownames(atac_Sample4)), function(x) wilcox.test(as.numeric(m1[x,]) , as.numeric(m2[x,]), alternative = "less")$p.value)
select_nonEnd_cells = sapply(end_sums, function(x) {
    matchSet = which(nonEnd_sums<=x+1000 & nonEnd_sums>=x-1000)
    if(length(matchSet)==0){matchSet = order(abs(nonEnd_sums-x))[1]}
    sample(matchSet, 1)
    })
m1 = ((as.matrix(atac_Sample4)[,!endCells_Sample4][,select_nonEnd_cells])>0)
peakScores3 = sapply(1:length(rownames(atac_Sample4)), function(x) wilcox.test(as.numeric(m1[x,]) , as.numeric(m2[x,]), alternative = "less")$p.value)

combinedPeakScores_End = p.adjust(sapply(1:length(peakScores1), function(x) sumlog(c(peakScores1[x],peakScores2[x],peakScores3[x]))[[3]]))

#Myeloid Cells
assign(paste0("myeCells_",samp), clusters_split==6)
mye_sums = colSums(as.matrix(atac_Sample4[,myeCells_Sample4])>0)
nonMye_sums = colSums(as.matrix(atac_Sample4[,!myeCells_Sample4])>0)

m2 = ((as.matrix(atac_Sample4)[,myeCells_Sample4])>0)

select_nonMye_cells = sapply(mye_sums, function(x) sample(which(nonMye_sums<=x+1500 & nonMye_sums>=x-1500), 1))
m1 = ((as.matrix(atac_Sample4)[,!myeCells_Sample4][,select_nonMye_cells])>0)
peakScores1 = sapply(1:length(rownames(atac_Sample4)), function(x) wilcox.test(as.numeric(m1[x,]) , as.numeric(m2[x,]), alternative = "less")$p.value)
select_nonMye_cells = sapply(mye_sums, function(x) sample(which(nonMye_sums<=x+1500 & nonMye_sums>=x-1500), 1))
m1 = ((as.matrix(atac_Sample4)[,!myeCells_Sample4][,select_nonMye_cells])>0)
peakScores2 = sapply(1:length(rownames(atac_Sample4)), function(x) wilcox.test(as.numeric(m1[x,]) , as.numeric(m2[x,]), alternative = "less")$p.value)
select_nonMye_cells = sapply(mye_sums, function(x) sample(which(nonMye_sums<=x+1500 & nonMye_sums>=x-1500), 1))
m1 = ((as.matrix(atac_Sample4)[,!myeCells_Sample4][,select_nonMye_cells])>0)
peakScores3 = sapply(1:length(rownames(atac_Sample4)), function(x) wilcox.test(as.numeric(m1[x,]) , as.numeric(m2[x,]), alternative = "less")$p.value)

combinedPeakScores_Mye = p.adjust(sapply(1:length(peakScores1), function(x) sumlog(c(peakScores1[x],peakScores2[x],peakScores3[x]))[[3]]))

write.table(c("track name=Fibroblast_TAC color=0,225,0"), "TAC.bed", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(data.frame(as(rownames(atac_Sample4), "GRanges"))[combinedPeakScores_Fib<.1,1:3], "TAC.bed", quote = FALSE, col.names = FALSE, row.names = FALSE, append = TRUE)
write.table(c("track name=Myeloid_TAC color=0,0,225"), "TAC.bed", quote = FALSE, col.names = FALSE, row.names = FALSE, append = TRUE)
write.table(data.frame(as(rownames(atac_Sample4), "GRanges"))[combinedPeakScores_Mye<.1,1:3], "TAC.bed", quote = FALSE, col.names = FALSE, row.names = FALSE, append = TRUE)
write.table(c("track name=Endothelial_TAC color=225,0,0"), "TAC.bed", quote = FALSE, col.names = FALSE, row.names = FALSE, append = TRUE)
write.table(data.frame(as(rownames(atac_Sample4), "GRanges"))[combinedPeakScores_End<.1,1:3], "TAC.bed", quote = FALSE, col.names = FALSE, row.names = FALSE, append = TRUE)

#Generic function for the above:
samp = "Sample1"
mat = as.matrix(get(paste0("atac_",samp)))
projTestMapSNE = get(paste0("atac_tNSE_",samp))
clusters = hclust(dist(projTestMapSNE$Y), method = "ward.D")
clusters_split = factor(cutree(clusters, 9))
qplot(projTestMapSNE$Y[,1],projTestMapSNE$Y[,2]) + geom_point(aes(projTestMapSNE$Y[mat[which(countOverlaps(as(rownames(mat), "GRanges"),mmProm[mmGenes[["ENSMUSG00000062515"]]$tx_name])>0),]>0,1],projTestMapSNE$Y[mat[which(countOverlaps(as(rownames(mat), "GRanges"),mmProm[mmGenes[["ENSMUSG00000062515"]]$tx_name])>0),]>0,2]), color="red") + geom_point(aes(projTestMapSNE$Y[colSums(mat[which(countOverlaps(as(rownames(mat), "GRanges"),mmProm[mmGenes[["ENSMUSG00000019929"]]$tx_name])>0),])>0,1],projTestMapSNE$Y[colSums(mat[which(countOverlaps(as(rownames(mat), "GRanges"),mmProm[mmGenes[["ENSMUSG00000019929"]]$tx_name])>0),])>0,2]), color="green") + geom_point(aes(projTestMapSNE$Y[mat[which(countOverlaps(as(rownames(mat), "GRanges"),mmProm[mmGenes[["ENSMUSG00000069516"]]$tx_name])>0),]>0,1],projTestMapSNE$Y[mat[which(countOverlaps(as(rownames(mat), "GRanges"),mmProm[mmGenes[["ENSMUSG00000069516"]]$tx_name])>0),]>0,2]), color="blue")
qplot(projTestMapSNE$Y[,1],projTestMapSNE$Y[,2], color=clusters_split)

for(samp in paste0("Sample", 1:8)){
    mat = as.matrix(get(paste0("atac_",samp)))
    projTestMapSNE = get(paste0("atac_tNSE_",samp))
    clusters = hclust(dist(projTestMapSNE$Y), method = "ward.D")
    clusters_split = factor(cutree(clusters, 9))
    overlaps = which(countOverlaps(as(rownames(mat), "GRanges"),mmProm[mmGenes[["ENSMUSG00000019929"]]$tx_name])>0)
    if(length(overlaps)>1){
        print(qplot(projTestMapSNE$Y[,1],projTestMapSNE$Y[,2]) + geom_point(aes(projTestMapSNE$Y[mat[which(countOverlaps(as(rownames(mat), "GRanges"),mmProm[mmGenes[["ENSMUSG00000062515"]]$tx_name])>0),]>0,1],projTestMapSNE$Y[mat[which(countOverlaps(as(rownames(mat), "GRanges"),mmProm[mmGenes[["ENSMUSG00000062515"]]$tx_name])>0),]>0,2]), color="red") + geom_point(aes(projTestMapSNE$Y[colSums(mat[which(countOverlaps(as(rownames(mat), "GRanges"),mmProm[mmGenes[["ENSMUSG00000019929"]]$tx_name])>0),])>0,1],projTestMapSNE$Y[colSums(mat[which(countOverlaps(as(rownames(mat), "GRanges"),mmProm[mmGenes[["ENSMUSG00000019929"]]$tx_name])>0),])>0,2]), color="green") + geom_point(aes(projTestMapSNE$Y[mat[which(countOverlaps(as(rownames(mat), "GRanges"),mmProm[mmGenes[["ENSMUSG00000069516"]]$tx_name])>0),]>0,1],projTestMapSNE$Y[mat[which(countOverlaps(as(rownames(mat), "GRanges"),mmProm[mmGenes[["ENSMUSG00000069516"]]$tx_name])>0),]>0,2]), color="blue"))
        print(sapply(1:9, function(x) length(which(clusters_split==x & colSums(mat[which(countOverlaps(as(rownames(mat), "GRanges"),mmProm[mmGenes[["ENSMUSG00000019929"]]$tx_name])>0),])>0)) / length(which(clusters_split==x))))
    }
    else{
        print(qplot(projTestMapSNE$Y[,1],projTestMapSNE$Y[,2]) + geom_point(aes(projTestMapSNE$Y[mat[which(countOverlaps(as(rownames(mat), "GRanges"),mmProm[mmGenes[["ENSMUSG00000062515"]]$tx_name])>0),]>0,1],projTestMapSNE$Y[mat[which(countOverlaps(as(rownames(mat), "GRanges"),mmProm[mmGenes[["ENSMUSG00000062515"]]$tx_name])>0),]>0,2]), color="red") + geom_point(aes(projTestMapSNE$Y[(mat[which(countOverlaps(as(rownames(mat), "GRanges"),mmProm[mmGenes[["ENSMUSG00000019929"]]$tx_name])>0),])>0,1],projTestMapSNE$Y[(mat[which(countOverlaps(as(rownames(mat), "GRanges"),mmProm[mmGenes[["ENSMUSG00000019929"]]$tx_name])>0),])>0,2]), color="green") + geom_point(aes(projTestMapSNE$Y[mat[which(countOverlaps(as(rownames(mat), "GRanges"),mmProm[mmGenes[["ENSMUSG00000069516"]]$tx_name])>0),]>0,1],projTestMapSNE$Y[mat[which(countOverlaps(as(rownames(mat), "GRanges"),mmProm[mmGenes[["ENSMUSG00000069516"]]$tx_name])>0),]>0,2]), color="blue"))
        print(sapply(1:9, function(x) length(which(clusters_split==x & mat[which(countOverlaps(as(rownames(mat), "GRanges"),mmProm[mmGenes[["ENSMUSG00000019929"]]$tx_name])>0),]>0)) / length(which(clusters_split==x))))
    }
    print(qplot(projTestMapSNE$Y[,1],projTestMapSNE$Y[,2], color=clusters_split))
}

for(samp in paste0("Sample", 1:8)){
	mat = as.matrix(get(paste0("atac_",samp)))
    projTestMapSNE = get(paste0("atac_tNSE_",samp))
    clusters = hclust(dist(projTestMapSNE$Y), method = "ward.D")
    clusters_split = factor(cutree(clusters, 9))
	if(samp == "Sample1"){clustIds = c(9,7,3)} #Sample1
	if(samp == "Sample2"){clustIds = c(7,3,5)} #Sample2
	if(samp == "Sample3"){clustIds = c(6,3,5)} #Sample3
	if(samp == "Sample4"){clustIds = c(7,4,9)} #Sample4
	if(samp == "Sample5"){clustIds = c(2,6,8)} #Sample5
	if(samp == "Sample6"){clustIds = c(2,6,8)} #Sample6
	if(samp == "Sample7"){clustIds = c(1,4,2)} #Sample7
	if(samp == "Sample8"){clustIds = c(list(c(5,2)),4,8)} #Sample8
	names(clustIds) = c("Fibroblast", "Endothelial", "Myeloid")

for(cType in c("Fibroblast", "Endothelial", "Myeloid")){

	assign(paste0(cType,"Cells_",samp), clusters_split %in% clustIds[[cType]])
	# assign(paste0(cType,"Cells_",samp), clusters_split==clustIds[cType])
	theseCells = get(paste0(cType,"Cells_",samp))
	cell_sums = colSums(mat[,theseCells]>0)
	nonCell_sums = colSums(mat[,!theseCells]>0)

	m2 = mat[,theseCells]>0

	select_cells = sapply(cell_sums, function(x) {
	    matchSet = which(nonCell_sums<=x+1000 & nonCell_sums>=x-1000)
	    if(length(matchSet)==0){matchSet = order(abs(nonCell_sums-x))[1]}
	    sample(matchSet, 1)
    })
	m1 = mat[,!theseCells][,select_cells]>0
	peakScores1 = sapply(1:length(rownames(mat)), function(x) wilcox.test(as.numeric(m1[x,]) , as.numeric(m2[x,]), alternative = "less")$p.value)
	select_cells = sapply(cell_sums, function(x) {
	    matchSet = which(nonCell_sums<=x+1000 & nonCell_sums>=x-1000)
	    if(length(matchSet)==0){matchSet = order(abs(nonCell_sums-x))[1]}
	    sample(matchSet, 1)
    })
	m1 = mat[,!theseCells][,select_cells]>0
	peakScores2 = sapply(1:length(rownames(mat)), function(x) wilcox.test(as.numeric(m1[x,]) , as.numeric(m2[x,]), alternative = "less")$p.value)
	select_cells = sapply(cell_sums, function(x) {
	    matchSet = which(nonCell_sums<=x+1000 & nonCell_sums>=x-1000)
	    if(length(matchSet)==0){matchSet = order(abs(nonCell_sums-x))[1]}
	    sample(matchSet, 1)
    })
	m1 = mat[,!theseCells][,select_cells]>0
	peakScores3 = sapply(1:length(rownames(mat)), function(x) wilcox.test(as.numeric(m1[x,]) , as.numeric(m2[x,]), alternative = "less")$p.value)
	
	combinedPeakScores = p.adjust(sapply(1:length(peakScores1), function(x) sumlog(c(peakScores1[x],peakScores2[x],peakScores3[x]))[[3]]))

	if(cType=="Fibroblast"){
		write.table(c(paste0("track name=",cType,"_",samp," color=0,225,0")), paste0(samp,".bed"), quote = FALSE, col.names = FALSE, row.names = FALSE)
		write.table(data.frame(as(rownames(mat), "GRanges"))[combinedPeakScores<.1,1:3], paste0(samp,".bed"), quote = FALSE, col.names = FALSE, row.names = FALSE, append = TRUE)
	}
	if(cType=="Endothelial"){
		write.table(c(paste0("track name=",cType,"_",samp," color=225,0,0")), paste0(samp,".bed"), quote = FALSE, col.names = FALSE, row.names = FALSE, append = TRUE)
		write.table(data.frame(as(rownames(mat), "GRanges"))[combinedPeakScores<.1,1:3], paste0(samp,".bed"), quote = FALSE, col.names = FALSE, row.names = FALSE, append = TRUE)
	}
	if(cType=="Myeloid"){
		write.table(c(paste0("track name=",cType,"_",samp," color=0,0,225")), paste0(samp,".bed"), quote = FALSE, col.names = FALSE, row.names = FALSE, append = TRUE)
		write.table(data.frame(as(rownames(mat), "GRanges"))[combinedPeakScores<.1,1:3], paste0(samp,".bed"), quote = FALSE, col.names = FALSE, row.names = FALSE, append = TRUE)
	}
}
}
					     
#just write out all pseudorep peak combinations
for(cType in cTypes){
allFibroPeaks = fread(paste0("all",cType,"Peaks_intersect.bed"))
allFibroPeak_Ranges = GRanges(allFibroPeaks$V1, IRanges(allFibroPeaks$V2, allFibroPeaks$V3))
allFibroPeak_Ranges = allFibroPeak_Ranges[countOverlaps(allFibroPeak_Ranges,unlist(mmGenes))==0 & countOverlaps(allFibroPeak_Ranges,mmProm)==0]
allFibroPeaksFrame = data.frame(allFibroPeak_Ranges)[,1:3]
allFibroPeaksFrame = data.frame(loc=paste0(allFibroPeaksFrame$seqnames,":",allFibroPeaksFrame$start,"-",allFibroPeaksFrame$end))
for(cond in c("Sham","TAC","JQ1","JQ1_Withdrawn")){
    fPeaks = fread(paste0(cType,cond,"_distal_pseudoRep.bed"))
    allFibroPeaksFrame = cbind(allFibroPeaksFrame, as.numeric(countOverlaps(allFibroPeak_Ranges, GRanges(fPeaks$V1, IRanges(as.numeric(fPeaks$V2), as.numeric(fPeaks$V3))))>0))
}
colnames(allFibroPeaksFrame) = c("loc","Sham","TAC","JQ1","JQ1_Withdrawn")
test = sapply(0:1, function(Sham) sapply(0:1, function(TAC) sapply(0:1, function(JQ1) sapply(0:1, function(JQ1_Withdrawn){
	thesePeaks = allFibroPeaksFrame[allFibroPeaksFrame$Sham==Sham & allFibroPeaksFrame$TAC==TAC & allFibroPeaksFrame$JQ1==JQ1 & allFibroPeaksFrame$JQ1_Withdrawn==JQ1_Withdrawn,]
	thesePeaks = as(thesePeaks$loc, "GRanges")
	if(length(thesePeaks)!=0){
	write.table(data.frame(data.frame(thesePeaks)[,1:3],"","","*"), paste0(paste0(cType,Sham,TAC,JQ1,JQ1_Withdrawn,"_distal_pseudoRep.bed")), row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")}
	}))))

}
#or just the opening/closing
for(comp in list(c("Sham", "TAC"), c("TAC", "JQ1"), c("JQ1", "JQ1_Withdrawn"))){
    for(cType in cTypes){
        for(cond in comp){
            
            fPeaks = fread(paste0(cType,cond,"_distal_pseudoRep.bed"))
            
            assign(paste0("allCondPeaks",cond), GRanges(fPeaks$V1, IRanges(as.numeric(fPeaks$V2), as.numeric(fPeaks$V3))))
        }
        #split cell types
        file = paste(c(cType,comp,"closing"),collapse = "_")
        thesePeaks = get(paste0("allCondPeaks",comp[1]))[countOverlaps(get(paste0("allCondPeaks",comp[1])), get(paste0("allCondPeaks",comp[2])))==0]
        write.table(data.frame(data.frame(thesePeaks)[,1:3],"","","*"), paste0(file,"_distal_pseudoRep.bed"), quote = FALSE, col.names = FALSE, row.names = FALSE, sep="\t")
        file = paste(c(cType,comp,"opening"),collapse = "_")
        thesePeaks = get(paste0("allCondPeaks",comp[2]))[countOverlaps(get(paste0("allCondPeaks",comp[2])), get(paste0("allCondPeaks",comp[1])))==0]
        write.table(data.frame(data.frame(thesePeaks)[,1:3],"","","*"), paste0(file,"_distal_pseudoRep.bed"), quote = FALSE, col.names = FALSE, row.names = FALSE, sep="\t")
        file = paste(c(cType,comp,"same"),collapse = "_")
        thesePeaks = reduce(union(get(paste0("allCondPeaks",comp[2]))[countOverlaps(get(paste0("allCondPeaks",comp[2])), get(paste0("allCondPeaks",comp[1])))>0], get(paste0("allCondPeaks",comp[1]))[countOverlaps(get(paste0("allCondPeaks",comp[1])), get(paste0("allCondPeaks",comp[2])))>0]))
        write.table(data.frame(data.frame(thesePeaks)[,1:3],"","","*"), paste0(file,"_distal_pseudoRep.bed"), quote = FALSE, col.names = FALSE, row.names = FALSE, sep="\t")
        
    }
}
#and all peaks if needed
for(cType in cTypes){
    allFibroPeaks = fread(paste0("all",cType,"Peaks_intersect.bed"))
    allFibroPeak_Ranges = GRanges(allFibroPeaks$V1, IRanges(allFibroPeaks$V2, allFibroPeaks$V3))
    allFibroPeaksFrame = data.frame(allFibroPeak_Ranges)[,1:3]
    allFibroPeaksFrame = data.frame(loc=paste0(allFibroPeaksFrame$seqnames,":",allFibroPeaksFrame$start,"-",allFibroPeaksFrame$end))
    for(cond in c("Sham","TAC","JQ1","JQ1_Withdrawn")){
        fPeaks = fread(paste0(cType,cond,"_distal.bed"))
        allFibroPeaksFrame = cbind(allFibroPeaksFrame, as.numeric(countOverlaps(allFibroPeak_Ranges, GRanges(fPeaks$V1, IRanges(as.numeric(fPeaks$V2), as.numeric(fPeaks$V3))))>0))
        
        writeRanges = data.frame(allFibroPeak_Ranges[countOverlaps(allFibroPeak_Ranges, GRanges(fPeaks$V1, IRanges(as.numeric(fPeaks$V2), as.numeric(fPeaks$V3))))>0])
        write.table(data.frame(writeRanges[,1:3],"","","*"), paste0(cType,cond,"_allPeaks_pseudoRep.bed"), row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
    }
}
