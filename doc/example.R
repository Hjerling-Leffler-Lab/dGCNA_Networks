library(Matrix)
library(WGCNA)
library(ggplot2)
library(igraph)
library(dendextend)
#LOAD "Networks" packge HERE
library(colorspace)
library(lme4)
library(blme)
library(performance)
library(reshape2)
library(parallel)

setwd('/YourWorkingDirectory/')
dir<-"/DirectoryOfYourFiles/"

Disease<-readRDS(file=paste0(dir,"Disease_Annotation.rds"))

samples<-readRDS(file=paste0(dir,"SampleList.rds"))
clustersAnnotation<-readRDS(file=paste0(dir,"clustersAnnotation.rds"))

NCORES<-16
MIN_CELLS<-30

print(paste0("Number of cores detected: ",detectCores()))

#Beta cells is cluster 3, 1st level conos clusterig for merged datasets
sa<-SelectCellCommunities(samples,clustersAnnotation,c("3"),min_cells = MIN_CELLS)

#sfn<-GetFixedNumberOfCells(sa,NumberOfCells = FIXED_NUMBER_CELLS)

states<-SamplesStates(sa,Disease)
states

GetNumberOfCells(samples,clustersAnnotation,c("3"),states,save=TRUE,filename="Ncells.rds")


genes<-GetCommonExpressedGenes(sa,samp=NULL,min_expr=1,met="NumberOfCells",min_cells=5,cell_proportion=0.3,print=TRUE,FileName="BackgroundGeneSet.txt")
genes

lm<-GetLinearModel(sa,sample_groups=states,genes=genes,model="blmer",normalization=TRUE)

sm<-FindSimilaritiesForNestedData(data_lm=lm,method="pearson")

saveRDS(sm,file="similarities_matrices.rds")

mask<-GetMaskFromRandomDonors(input_data=sa,sample_groups=states,
                              SimilarityDifferences=sm[["Disease"]]-sm[["Control"]],
                              genes=genes,model="blmer",normalization=TRUE,replacement=FALSE,method="pearson",
                              T=512,Nctrl=10,Ncases=10,verbose=TRUE,ncores=NCORES)

saveRDS(mask,file="mask.rds")

connections<-GetNetworkConnectionsFromMask(SimilarityDifferences=sm[["Disease"]]-sm[["Control"]],
                                           mask=mask,TH_MASK=0.975,weighted=TRUE)

connections_after_pruning<-NetworkPruning(connections,weight_th = 0)

#clustering with WGCNA package

#WGCNA
dissTOM=TOMdist(as.matrix(connections_after_pruning),TOMType="signed")
row.names(dissTOM)<-row.names(connections_after_pruning)

hierTOM = hclust(as.dist(dissTOM),method="average");

ModulesList<-list()

#4 different cut-offs for the dendrogram
ModulesList[[1]]=labels2colors (cutreeDynamic(hierTOM,method="tree",minClusterSize = 30,cutHeight = 0.974,deepSplit=FALSE))
ModulesList[[2]]=labels2colors (cutreeDynamic(hierTOM,method="tree",minClusterSize = 30,cutHeight = 0.981,deepSplit=FALSE))
ModulesList[[3]]=labels2colors (cutreeDynamic(hierTOM,method="tree",minClusterSize = 30,cutHeight = 0.987,deepSplit=FALSE))
ModulesList[[4]]=labels2colors (cutreeDynamic(hierTOM,method="tree",minClusterSize = 30,cutHeight = 0.992,deepSplit=FALSE))

plotDendroAndColors(hierTOM, data.frame(ModulesList), c("Modules1","Modules2","Modules3","Modules4"),dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05, main="Hierarchical clustering")

