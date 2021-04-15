
samp = "Sample1"

projTest = get(paste0("proj_",samp))
mat = as.matrix(get(paste0("atac_",samp)))
lsa_unit_dist = t(scale(t(projTest[,2:16])))
lsa_unit_dist_e = dist(lsa_unit_dist)
lsa_unit_dist_e_frac = 1/lsa_unit_dist_e
cells_in_markers_all = get(paste0("cells_in_markers_",samp))

projTestMapSNE = Rtsne(t(scale(t(projTest[,2:16]))), pca=FALSE, normalize=FALSE, max_iter = 5000)
tSNE_dist_e = dist(projTestMapSNE$Y[,1:2])
testDists1fabp4 = sapply(which(mat[which(countOverlaps(as(rownames(mat), "GRanges"),mmProm[mmGenes[["ENSMUSG00000062515"]]$tx_name])>0),]>0), function(x) median(rank(as.matrix(tSNE_dist_e)[x,])[which(mat[which(countOverlaps(as(rownames(mat), "GRanges"),mmProm[mmGenes[["ENSMUSG00000062515"]]$tx_name])>0),]>0)]))
testDists1dcn = sapply(which(colSums(mat[which(countOverlaps(as(rownames(mat), "GRanges"),mmProm[mmGenes[["ENSMUSG00000019929"]]$tx_name])>0),])>0), function(x) median(rank(as.matrix(tSNE_dist_e)[x,])[which(colSums(mat[which(countOverlaps(as(rownames(mat), "GRanges"),mmProm[mmGenes[["ENSMUSG00000019929"]]$tx_name])>0),])>0)]))
testDists1lyz2 = sapply(which(mat[which(countOverlaps(as(rownames(mat), "GRanges"),mmProm[mmGenes[["ENSMUSG00000069516"]]$tx_name])>0),]>0), function(x) median(rank(as.matrix(tSNE_dist_e)[x,])[which(mat[which(countOverlaps(as(rownames(mat), "GRanges"),mmProm[mmGenes[["ENSMUSG00000069516"]]$tx_name])>0),]>0)]))

scale = c(.01,.1,.5,1,2,10,20,50)
for(s in scale){

simMat = cbind(rbind(matrix(0,i,i),t(cells_in_markers_all)), rbind(s*cells_in_markers_all,as.matrix(lsa_unit_dist_e_frac)))
infMat = random.walk(simMat)

projTestMapSNE = Rtsne(t(scale(t(cbind(projTest[,2:16],as.matrix(t(1000000*infMat[1:i,-(1:i)])))))), pca=FALSE, normalize=FALSE, max_iter = 5000)
tSNE_dist_e = dist(projTestMapSNE$Y[,1:2])

assign(paste0("testDists2fabp4",s), sapply(which(mat[which(countOverlaps(as(rownames(mat), "GRanges"),mmProm[mmGenes[["ENSMUSG00000062515"]]$tx_name])>0),]>0), function(x) median(rank(as.matrix(tSNE_dist_e)[x,])[which(mat[which(countOverlaps(as(rownames(mat), "GRanges"),mmProm[mmGenes[["ENSMUSG00000062515"]]$tx_name])>0),]>0)])))

assign(paste0("testDists2dcn",s), sapply(which(colSums(mat[which(countOverlaps(as(rownames(mat), "GRanges"),mmProm[mmGenes[["ENSMUSG00000019929"]]$tx_name])>0),])>0), function(x) median(rank(as.matrix(tSNE_dist_e)[x,])[which(colSums(mat[which(countOverlaps(as(rownames(mat), "GRanges"),mmProm[mmGenes[["ENSMUSG00000019929"]]$tx_name])>0),])>0)])) )

assign(paste0("testDists2lyz2",s), sapply(which(mat[which(countOverlaps(as(rownames(mat), "GRanges"),mmProm[mmGenes[["ENSMUSG00000069516"]]$tx_name])>0),]>0), function(x) median(rank(as.matrix(tSNE_dist_e)[x,])[which(mat[which(countOverlaps(as(rownames(mat), "GRanges"),mmProm[mmGenes[["ENSMUSG00000069516"]]$tx_name])>0),]>0)])) )

}

projTestMapSNE = Rtsne(t(scale(t(cbind(projTest[,2:16],as.matrix(t(1000*cells_in_markers_all)))))), pca=FALSE, normalize=FALSE, max_iter = 5000)
tSNE_dist_e = dist(projTestMapSNE$Y[,1:2])
testDists3fabp4 = sapply(which(mat[which(countOverlaps(as(rownames(mat), "GRanges"),mmProm[mmGenes[["ENSMUSG00000062515"]]$tx_name])>0),]>0), function(x) median(rank(as.matrix(tSNE_dist_e)[x,])[which(mat[which(countOverlaps(as(rownames(mat), "GRanges"),mmProm[mmGenes[["ENSMUSG00000062515"]]$tx_name])>0),]>0)]))
testDists3dcn = sapply(which(colSums(mat[which(countOverlaps(as(rownames(mat), "GRanges"),mmProm[mmGenes[["ENSMUSG00000019929"]]$tx_name])>0),])>0), function(x) median(rank(as.matrix(tSNE_dist_e)[x,])[which(colSums(mat[which(countOverlaps(as(rownames(mat), "GRanges"),mmProm[mmGenes[["ENSMUSG00000019929"]]$tx_name])>0),])>0)]))
testDists3lyz2 = sapply(which(mat[which(countOverlaps(as(rownames(mat), "GRanges"),mmProm[mmGenes[["ENSMUSG00000069516"]]$tx_name])>0),]>0), function(x) median(rank(as.matrix(tSNE_dist_e)[x,])[which(mat[which(countOverlaps(as(rownames(mat), "GRanges"),mmProm[mmGenes[["ENSMUSG00000069516"]]$tx_name])>0),]>0)]))

ggplot() + geom_point(aes(x=c(.01,.1,.5,1,2,10,20,50), y=sapply(c(.01,.1,.5,1,2,10,20,50), function(x) median(get(paste0("testDists2dcn",x)))))) + geom_line(aes(x=c(.01,.1,.5,1,2,10,20,50), y=sapply(c(.01,.1,.5,1,2,10,20,50), function(x) median(get(paste0("testDists2dcn",x)))))) + geom_hline(aes(yintercept=median(testDists1dcn)), linetype="dashed", color="green") + geom_hline(aes(yintercept=median(testDists3dcn)), linetype="dashed", color="red") + xlab("RNA-seq edge weight") + ylab("Median tSNE Distance Rank: Dcn")

ggplot() + geom_point(aes(x=c(.01,.1,.5,1,2,10,20,50), y=sapply(c(.01,.1,.5,1,2,10,20,50), function(x) median(get(paste0("testDists2lyz2",x)))))) + geom_line(aes(x=c(.01,.1,.5,1,2,10,20,50), y=sapply(c(.01,.1,.5,1,2,10,20,50), function(x) median(get(paste0("testDists2lyz2",x)))))) + geom_hline(aes(yintercept=median(testDists1lyz2)), linetype="dashed", color="green") + geom_hline(aes(yintercept=median(testDists3lyz2)), linetype="dashed", color="red") + xlab("RNA-seq edge weight") + ylab("Median tSNE Distance Rank: Lyz2")

ggplot() + geom_point(aes(x=c(.01,.1,.5,1,2,10,20,50), y=sapply(c(.01,.1,.5,1,2,10,20,50), function(x) median(get(paste0("testDists2fabp4",x)))))) + geom_line(aes(x=c(.01,.1,.5,1,2,10,20,50), y=sapply(c(.01,.1,.5,1,2,10,20,50), function(x) median(get(paste0("testDists2fabp4",x)))))) + geom_hline(aes(yintercept=median(testDists1fabp4)), linetype="dashed", color="green") + geom_hline(aes(yintercept=median(testDists3fabp4)), linetype="dashed", color="red") + xlab("RNA-seq edge weight") + ylab("Median tSNE Distance Rank: Fabp4")

#Sample1: s=1
s=1
simMat = cbind(rbind(matrix(0,i,i),t(cells_in_markers_all)), rbind(s*cells_in_markers_all,as.matrix(lsa_unit_dist_e_frac)))
infMat = random.walk(simMat)
projTestMapSNE = Rtsne(t(scale(t(cbind(projTest[,2:16],as.matrix(t(1000000*infMat[1:i,-(1:i)])))))), pca=FALSE, normalize=FALSE, max_iter = 5000)
qplot(projTestMapSNE$Y[,1],projTestMapSNE$Y[,2]) + geom_point(aes(projTestMapSNE$Y[mat[which(countOverlaps(as(rownames(mat), "GRanges"),mmProm[mmGenes[["ENSMUSG00000062515"]]$tx_name])>0),]>0,1],projTestMapSNE$Y[mat[which(countOverlaps(as(rownames(mat), "GRanges"),mmProm[mmGenes[["ENSMUSG00000062515"]]$tx_name])>0),]>0,2]), color="red") + geom_point(aes(projTestMapSNE$Y[colSums(mat[which(countOverlaps(as(rownames(mat), "GRanges"),mmProm[mmGenes[["ENSMUSG00000019929"]]$tx_name])>0),])>0,1],projTestMapSNE$Y[colSums(mat[which(countOverlaps(as(rownames(mat), "GRanges"),mmProm[mmGenes[["ENSMUSG00000019929"]]$tx_name])>0),])>0,2]), color="green") + geom_point(aes(projTestMapSNE$Y[mat[which(countOverlaps(as(rownames(mat), "GRanges"),mmProm[mmGenes[["ENSMUSG00000069516"]]$tx_name])>0),]>0,1],projTestMapSNE$Y[mat[which(countOverlaps(as(rownames(mat), "GRanges"),mmProm[mmGenes[["ENSMUSG00000069516"]]$tx_name])>0),]>0,2]), color="blue")

clusters = hclust(dist(projTestMapSNE$Y), method = "ward.D")
clusters_split = factor(cutree(clusters, 9))
qplot(projTestMapSNE$Y[,1],projTestMapSNE$Y[,2], color=clusters_split) 
#Sample1: clust 8
assign(paste0("fibCells_",samp), clusters_split==8)
assign(paste0("atac_tNSE_",samp), projTestMapSNE)


#Sample3: s=1
s=.5
simMat = cbind(rbind(matrix(0,i,i),t(cells_in_markers_all)), rbind(s*cells_in_markers_all,as.matrix(lsa_unit_dist_e_frac)))
infMat = random.walk(simMat)
projTestMapSNE = Rtsne(t(scale(t(cbind(projTest[,2:16],as.matrix(t(1000000*infMat[1:i,-(1:i)])))))), pca=FALSE, normalize=FALSE, max_iter = 5000)
qplot(projTestMapSNE$Y[,1],projTestMapSNE$Y[,2]) + geom_point(aes(projTestMapSNE$Y[mat[which(countOverlaps(as(rownames(mat), "GRanges"),mmProm[mmGenes[["ENSMUSG00000062515"]]$tx_name])>0),]>0,1],projTestMapSNE$Y[mat[which(countOverlaps(as(rownames(mat), "GRanges"),mmProm[mmGenes[["ENSMUSG00000062515"]]$tx_name])>0),]>0,2]), color="red") + geom_point(aes(projTestMapSNE$Y[colSums(mat[which(countOverlaps(as(rownames(mat), "GRanges"),mmProm[mmGenes[["ENSMUSG00000019929"]]$tx_name])>0),])>0,1],projTestMapSNE$Y[colSums(mat[which(countOverlaps(as(rownames(mat), "GRanges"),mmProm[mmGenes[["ENSMUSG00000019929"]]$tx_name])>0),])>0,2]), color="green") + geom_point(aes(projTestMapSNE$Y[mat[which(countOverlaps(as(rownames(mat), "GRanges"),mmProm[mmGenes[["ENSMUSG00000069516"]]$tx_name])>0),]>0,1],projTestMapSNE$Y[mat[which(countOverlaps(as(rownames(mat), "GRanges"),mmProm[mmGenes[["ENSMUSG00000069516"]]$tx_name])>0),]>0,2]), color="blue")

clusters = hclust(dist(projTestMapSNE$Y), method = "ward.D")
clusters_split = factor(cutree(clusters, 9))
qplot(projTestMapSNE$Y[,1],projTestMapSNE$Y[,2], color=clusters_split)
#Sample3: clust 1
assign(paste0("fibCells_",samp), clusters_split==1)
assign(paste0("atac_tNSE_",samp), projTestMapSNE)



#Sample 2 and 4 (with s=1 for now):
for(samp in paste0("Sample", c(2,4))){
    projTest = get(paste0("proj_",samp))
    mat = as.matrix(get(paste0("atac_",samp)))
    lsa_unit_dist = t(scale(t(projTest[,2:16])))
    lsa_unit_dist_e = dist(lsa_unit_dist)
    lsa_unit_dist_e_frac = 1/lsa_unit_dist_e
    cells_in_markers_all = get(paste0("cells_in_markers_",samp))
    
    s=1
    simMat = cbind(rbind(matrix(0,i,i),t(cells_in_markers_all)), rbind(s*cells_in_markers_all,as.matrix(lsa_unit_dist_e_frac)))
    infMat = random.walk(simMat)
    projTestMapSNE = Rtsne(t(scale(t(cbind(projTest[,2:16],as.matrix(t(1000000*infMat[1:i,-(1:i)])))))), pca=FALSE, normalize=FALSE, max_iter = 5000)
    assign(paste0("atac_tNSE_",samp), projTestMapSNE)

}
#s2:
clusters = hclust(dist(projTestMapSNE$Y), method = "ward.D")
clusters_split = factor(cutree(clusters, 9))
assign(paste0("fibCells_",samp), clusters_split==7)
#s4:
clusters = hclust(dist(projTestMapSNE$Y), method = "ward.D")
clusters_split = factor(cutree(clusters, 9))
assign(paste0("fibCells_",samp), clusters_split==5)
