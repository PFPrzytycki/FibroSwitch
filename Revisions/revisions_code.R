#motif tf co-enrichment:
cellType = "fibro"
cond1 = "Sham"
cond2 = "TAC"
learningData = fread(paste0("Mouse_Stress_TF_pseudoRep/LearningData_",
                            cellType,"_",cond1,"_",cond2,"_",
                            cellType,"_",cond2,"_",cond1,"_2000_5e+05.csv"))
learningDataTrim = data.frame(learningData[,3:157])
learningDataTrim = learningDataTrim[,colSums(learningDataTrim)!=0]
correlations = cor(learningDataTrim) 

toGRanges <- function(x){GRanges(x[[1]], IRanges(x[[2]],x[[3]]))}
extend <- function(x, upstream=0, downstream=0) {
  on_plus <- strand(x) == "+" | strand(x) == "*"
  new_start <- start(x) - ifelse(on_plus, upstream, downstream)
  new_end <- end(x) + ifelse(on_plus, downstream, upstream)
  ranges(x) <- IRanges(new_start, new_end)
  x
}

#gene bodies
mmGene = genes(TxDb.Mmusculus.UCSC.mm10.ensGene)
region = keepStandardChromosomes(mmGene, pruning.mode = "coarse")

#h3k27ac regions
Ctrl2Peaks = toGRanges(fread("Ctrl_TGFb_h3k23ac/K27Ac_Unstim2.2_mm10.chr_bcpPeaks.bed"))
Ctrl1Peaks = toGRanges(fread("Ctrl_TGFb_h3k23ac/K27Ac_Unstim1.1_mm10.chr_bcpPeaks.bed"))
TGFb1Peaks = toGRanges(fread("Ctrl_TGFb_h3k23ac/K27Ac_TGFb6.3_mm10.chr_bcpPeaks.bed"))
TGFb2Peaks = toGRanges(fread("Ctrl_TGFb_h3k23ac/K27Ac_TGFb7.3_genome_bcpPeaks.bed"))
CtrlPeaks = reduce(c(Ctrl1Peaks,Ctrl2Peaks))
TGFbPeaks = reduce(c(TGFb1Peaks, TGFb2Peaks))
region = intersect(TGFbPeaks, CtrlPeaks, ignore.strand=TRUE)
#Do it useing LFC
allh3kPeaks = reduce(c(Ctrl1Peaks, Ctrl2Peaks, TGFb1Peaks, TGFb2Peaks))
allh3kPeaks = allh3kPeaks[countOverlaps(allh3kPeaks, extend(mmGene, 0,0))==0]
sigsCTRL_K27Ac = bamProfile("Ctrl_TGFb_h3k23ac/K27Ac_Unstim1.1_mm10.chr_q30.bam", allh3kPeaks, ss=FALSE, verbose=FALSE)
sigsCTRL_K27Ac = sapply(sigsCTRL_K27Ac, sum)
sigsTGFb_K27Ac = bamProfile("Ctrl_TGFb_h3k23ac/K27Ac_TGFb6.3_mm10.chr_q30.bam", allh3kPeaks, ss=FALSE, verbose=FALSE)
sigsTGFb_K27Ac = sapply(sigsTGFb_K27Ac, sum)
LFC_K27Ac = log((sigsTGFb_K27Ac+1)/(sigsCTRL_K27Ac+1))
allh3kPeaks_unchanged = allh3kPeaks[LFC_K27Ac < .5 & LFC_K27Ac > -.5]
region = allh3kPeaks_unchanged

regionWidth = width(region)

sigsCTRL1_Meox1 = bamCount("Ctrl_TGFb_MEOX1/New_Data/MEOX1_ChIPseq_antiHA.Rep1.bam", region, verbose=FALSE)
sigsCTRL1_Meox1_count = countBam("Ctrl_TGFb_MEOX1/New_Data/MEOX1_ChIPseq_antiHA.Rep1.bam")$records
sigsTGFb1_Meox1 = bamCount("Ctrl_TGFb_MEOX1/New_Data/HA_3_6.3_genome_q30.bam", region, verbose=FALSE)
sigsTGFb1_Meox1_count = countBam("Ctrl_TGFb_MEOX1/New_Data/HA_3_6.3_genome_q30.bam")$records
sigsCTRL2_Meox1 = bamCount("Ctrl_TGFb_MEOX1/New_Data/Dec19_HA_5.1_genome_q30.bam", region, verbose=FALSE)
sigsCTRL2_Meox1_count = countBam("Ctrl_TGFb_MEOX1/New_Data/Dec19_HA_5.1_genome_q30.bam")$records
sigsTGFb2_Meox1 = bamCount("Ctrl_TGFb_MEOX1/New_Data/HA_3_5.3_genome_q30.bam", region, verbose=FALSE)
sigsTGFb2_Meox1_count = countBam("Data/Ctrl_TGFb_MEOX1/New_Data/HA_3_5.3_genome_q30.bam")$records
sigsCTRL3_Meox1 = bamCount("Ctrl_TGFb_MEOX1/New_Data/Nov2_HA_1_2.1_genome_q30.bam", region, verbose=FALSE)
sigsCTRL3_Meox1_count = countBam("Ctrl_TGFb_MEOX1/New_Data/Nov2_HA_1_2.1_genome_q30.bam")$records
sigsTGFb3_Meox1 = bamCount("Ctrl_TGFb_MEOX1/New_Data/Nov2_HA_3_2.3_genome_q30.bam", region, verbose=FALSE)
sigsTGFb3_Meox1_count = countBam("Ctrl_TGFb_MEOX1/New_Data/Nov2_HA_3_2.3_genome_q30.bam")$records

ggplot() + geom_boxplot(aes("Unstim1",1000*sigsCTRL1_Meox1/regionWidth/sigsCTRL1_Meox1_count*1e6)) + 
  geom_boxplot(aes("Unstim2",1000*sigsCTRL2_Meox1/regionWidth/sigsCTRL2_Meox1_count*1e6)) + 
  geom_boxplot(aes("Unstim3",1000*sigsCTRL3_Meox1/regionWidth/sigsCTRL3_Meox1_count*1e6)) + 
  geom_boxplot(aes("TGFb1",1000*sigsTGFb1_Meox1/regionWidth/sigsTGFb1_Meox1_count*1e6)) + 
  geom_boxplot(aes("TGFb2",1000*sigsTGFb2_Meox1/regionWidth/sigsTGFb2_Meox1_count*1e6)) + 
  geom_boxplot(aes("TGFb3",1000*sigsTGFb3_Meox1/regionWidth/sigsTGFb3_Meox1_count*1e6)) + 
  scale_y_log10() + theme_classic() + ylab("Meox1 Coverage") + xlab("Condition")

sigsCTRL_Meox1 = (1000*sigsCTRL1_Meox1/regionWidth/sigsCTRL1_Meox1_count*1e6 + 
                    1000*sigsCTRL2_Meox1/regionWidth/sigsCTRL2_Meox1_count*1e6 + 
                    1000*sigsCTRL3_Meox1/regionWidth/sigsCTRL3_Meox1_count*1e6)/3
sigsTGFb_Meox1 = (1000*sigsTGFb1_Meox1/regionWidth/sigsTGFb1_Meox1_count*1e6 + 
                    1000*sigsTGFb2_Meox1/regionWidth/sigsTGFb2_Meox1_count*1e6 + 
                    1000*sigsTGFb3_Meox1/regionWidth/sigsTGFb3_Meox1_count*1e6)/3

ggplot() + geom_boxplot(aes("Unstim",(1000*sigsCTRL1_Meox1/regionWidth/sigsCTRL1_Meox1_count*1e6 + 
                              1000*sigsCTRL2_Meox1/regionWidth/sigsCTRL2_Meox1_count*1e6 + 
                              1000*sigsCTRL3_Meox1/regionWidth/sigsCTRL3_Meox1_count*1e6)/3)) + 
  geom_boxplot(aes("TGFb",(1000*sigsTGFb1_Meox1/regionWidth/sigsTGFb1_Meox1_count*1e6 + 
                     1000*sigsTGFb2_Meox1/regionWidth/sigsTGFb2_Meox1_count*1e6 + 
                     1000*sigsTGFb3_Meox1/regionWidth/sigsTGFb3_Meox1_count*1e6)/3)) + 
  scale_y_log10() + theme_classic() + ylab("Meox1 Coverage") + xlab("Condition")

pvals = sapply(1:length(region), function(x) t.test(c(1000*sigsCTRL1_Meox1[x]/regionWidth[x]/sigsCTRL1_Meox1_count*1e6,
                                                      1000*sigsCTRL2_Meox1[x]/regionWidth[x]/sigsCTRL2_Meox1_count*1e6,
                                                      1000*sigsCTRL3_Meox1[x]/regionWidth[x]/sigsCTRL3_Meox1_count*1e6),
                                                    c(1000*sigsTGFb1_Meox1[x]/regionWidth[x]/sigsTGFb1_Meox1_count*1e6,
                                                      1000*sigsTGFb2_Meox1[x]/regionWidth[x]/sigsTGFb2_Meox1_count*1e6,
                                                      1000*sigsTGFb3_Meox1[x]/regionWidth[x]/sigsTGFb3_Meox1_count*1e6))$p.value)

LFC_TGFb_Unstim = log(sigsTGFb_Meox1/sigsCTRL_Meox1)
qplot(LFC_TGFb_Unstim, -log10(pvals)) + theme_classic() + xlab("LFC TGFb vs Unstim Meox1 Coverage") + ylab("-log10 p-value") + 
  geom_hline(aes(yintercept=-log10(1/length(region))), linetype="dashed")

acPeaks = toGRanges(fread("ChIP_K27Ac_based/TGFb_AcPeaks_distal.bed"))
Meox_TGFb = toGRanges(fread("Ctrl_TGFb_MEOX1/New_Data/intersect_MEOX1_TGFb_for_GREAT_6335.bed"))
Meox_Unstim = toGRanges(fread("Ctrl_TGFb_MEOX1/New_Data/intersect_MEOX1_Unstim_for_GREAT_5845.bed"))
Ac_noMeox1 = acPeaks[countOverlaps(acPeaks, Meox_TGFb)==0 & countOverlaps(acPeaks, Meox_Unstim)==0]
Ac_Meox1 = acPeaks[countOverlaps(acPeaks, Meox_TGFb)>0 | countOverlaps(acPeaks, Meox_Unstim)>0]

write.table(data.frame(Ac_noMeox1)[,1:3], "TGFb_AcPeaks_distal_noMeox1.bed", row.names = FALSE, col.names = FALSE, quote=FALSE, sep="\t")
write.table(data.frame(Ac_Meox1)[,1:3], "TGFb_AcPeaks_distal_Meox1.bed", row.names = FALSE, col.names = FALSE, quote=FALSE, sep="\t")

sigsCTRL_K27Ac = bamProfile("Ctrl_TGFb_h3k23ac/K27Ac_Unstim1.1_mm10.chr_q30.bam", Ac_noMeox1, ss=FALSE, verbose=FALSE)
sigsCTRL_K27Ac = sapply(sigsCTRL_K27Ac, sum)
sigsTGFb_K27Ac = bamProfile("Data/Ctrl_TGFb_h3k23ac/K27Ac_TGFb6.3_mm10.chr_q30.bam", Ac_noMeox1, ss=FALSE, verbose=FALSE)
sigsTGFb_K27Ac = sapply(sigsTGFb_K27Ac, sum)
LFC_K27Ac = log((sigsTGFb_K27Ac+1)/(sigsCTRL_K27Ac+1))
Ac_noMeox1_highLFC = Ac_noMeox1[LFC_K27Ac>1]

write.table(data.frame(Ac_noMeox1_highLFC)[,1:3], "TGFb_AcPeaks_distal_highLFC_noMeox1.bed", row.names = FALSE, col.names = FALSE, quote=FALSE, sep="\t")
