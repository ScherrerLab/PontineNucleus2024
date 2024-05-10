##Adam Hantman seq data (pontine grey)
library(Seurat)
library(ggplot2)
library(sctransform)
library(loomR)
library(hdf5r)
library(factoextra)
library(cluster)
library(glmGamPoi)
library(Matrix)
library(dplyr)
library(future)
library(UpSetR)
library(caTools)
library(ape)
library(ggtree)
library(grid)
library(dplyr)
library(patchwork)
library(UpSetR)
library(eulerr)
#BPN1
source('/proj/gs25/users/Jesse/scripts/DSSeurat_functions.R')
genes=read.table('features.tsv')
genes=genes$V2
cells=read.table('PN1.bcs.txt')
cells=cells$V1
cells2=read.table('PN2.bcs.txt')
cells2=cells2$V1
datamat=readMM('PN1.countmat.txt')
datamat2=readMM('PN2.countmat.txt')
datamat=as.matrix(datamat)
datamat2=as.matrix(datamat2)
rownames(datamat)=genes
colnames(datamat)=cells
rownames(datamat2)=genes
colnames(datamat2)=cells2
df1=rowsum(datamat,rownames(datamat))
df2=rowsum(datamat2,rownames(datamat2))


s1<-CreateSeuratObject(counts=df1,min.features=200,min.cells=5)
s2<-CreateSeuratObject(counts=df2,min.features=200,min.cells=5)

s1[['percent.mt']] <- PercentageFeatureSet(s1, pattern='^mt-') #max 9.837%
s1@meta.data$orig.ident='PN1'
s2[['percent.mt']] <- PercentageFeatureSet(s2, pattern='^mt-') # max 11.68831%
s2@meta.data$orig.ident='PN2'
dir.create('qc')
pdf('qc/pn1violins.pdf')
print(VlnPlot(s1,features=c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
dev.off()
pdf('qc/pn2violins.pdf')
print(VlnPlot(s2,features=c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
dev.off()

mat=df1
nfeatcell=colSums(mat!=0)#get num features for each cell
#calculate point with largest 'jump' in num genes/cell as feature cutoff
hh=c()
num=seq(0.1,99.9,.1)
for (ii in 1:999){
jj=num[ii]
q=quantile(nfeatcell,(jj/100))
hh[ii]=q
}
minfeat=data.frame(x=hh,y=0)
for (a in 1:999){
if (a > 1 ){
y=minfeat$x[a]-minfeat$x[a-1]
minfeat[a,2]=y
}
else {
minfeat[a,2]=0
}
}
#calculate point halfway b/w largest 'jump' in genes/cell

minfeatcutoff=round(minfeat$x[which.max(minfeat$y)])
minfeatcutoff2=round(minfeat$x[which.max(minfeat$y)-1])
minfin=minfeatcutoff-((minfeatcutoff-minfeatcutoff2)/2)


s1 <- SCTransform(s1, method = "glmGamPoi", verbose = FALSE,ncells=ncol(mat)/2,variable.features.n=quantile(nfeatcell,.5)) #changing variable features has little effect
s2 <- SCTransform(s2, method = "glmGamPoi", verbose = FALSE,ncells=ncol(mat)/2,variable.features.n=quantile(nfeatcell,.5)) #changing variable features has little effect
s.list=c(s1,s2)
features <- SelectIntegrationFeatures(object.list = s.list)
s.list<-PrepSCTIntegration(object.list=s.list,anchor.features=features)
s.10x.anchors=FindIntegrationAnchors(object.list=s.list,normalization.method='SCT',anchor.features = features)
s.10x.combined=IntegrateData(anchorset=s.10x.anchors,normalization.method='SCT')
s.10x.combined <- RunPCA(s.10x.combined, verbose = FALSE,npcs=100)

pdf('qc/Elbow.PN.10x.combined.pdf')
print(ElbowPlot(s.10x.combined,ndims=100))
dev.off()

co1=20

#Next to UMAP: for smaller datasets, decrease n.neighbors parameter
#To automate, this parameter should be a proportion of cells 
s.10x.combined <- RunUMAP(s.10x.combined, dims = 1:co1, verbose = FALSE)

pdf('qc/UMAP_bysample.pdf'),width=10,height=8)
print(DimPlot(s.10x.combined,group.by='orig.ident'))
dev.off()

saveRDS(s.10x.combined,'PN.10x.thruNorm_20PCs.rds')


minclustnum=round(ncol(s.10x.combined)/100) #smallest clust

#manually decide number of NNs
 if (ncol(s.10x.combined)>250) {
      mink=round((ncol(s.10x.combined)/100)+1)
       } else {
        mink=3
        }
if (ncol(s.10x.combined) > 2000) {
mink = 20
}
kparam=seq(mink,minclustnum)#what's the minimum cluster size? 10% of the population?
avgsils=rep(NA,minclustnum)
medsils=rep(NA,minclustnum)
badcellnum=rep(NA,minclustnum)
numclust=rep(NA,minclustnum)
dist.matrix <- dist(x = Embeddings(object = s.10x.combined[['pca']])[, 1:co1])
print('Iterative calculation of silhouettes')
for (k in kparam) {
	#print(paste0('Using ',k,' neighbors'))
	d=s.10x.combined #create a 'testing' object so you don't ruin anything
	d <- FindNeighbors(d, dims = 1:co1, verbose = FALSE,k.param=k)#first construct a shared nearest neighbors graph. K param indicates how many neighbors to use
	
	d <- FindClusters(d, verbose = FALSE,resolution=0.8)#test resolution for different sil scores?
	numclust[k]=length(unique(Idents(d)))
	#print(paste0('we get ',length(unique(Idents(d))),' clusters'))
	
	#should we optimize n neighbors + resolution by silhouette scores?
		
	clusters=d$seurat_clusters
	sil <- silhouette(x = as.numeric(x = as.factor(x = clusters)), dist = dist.matrix)
	#d$sil <- sil[, 3]
	avgsils[k]=mean(sil[,3])
	medsils[k]=median(sil[,3])
	
	badcells=length(which(sil[,3]<0))
	print(paste0(k,' neighbors; ',length(unique(Idents(d))),' clusters; ',badcells,' bad cells'))
	badcellnum[k]=badcells
}

optavgsil=which.max(avgsils)
optmedsil=which.max(medsils)
print(paste0('opt avg sil = ',optavgsil))


min.sils = round(min(avgsils,na.rm=T),2) - 0.01
max.sils = round(max(avgsils,na.rm=T),2) + 0.01
min.clust = round(min(numclust,na.rm=T),2) - 5
max.clust = round(max(numclust,na.rm=T),2) + 5
min.bad=round(min(badcellnum,na.rm=T),5)-10
max.bad = round(max(badcellnum,na.rm=T),5) + 10

pdf(paste0(qc,'/Sils_optKparams_newmeth.pdf'))
print(ggplot(as.data.frame(avgsils), aes(x=seq(1:minclustnum),y=avgsils),useDingbats=FALSE) + geom_point() + geom_smooth(span=0.5) + scale_x_continuous(breaks=c(seq(0,minclustnum,10)),limits=c(0,minclustnum),minor_breaks=NULL) + scale_y_continuous(breaks=c(seq(min.sils,max.sils,0.01)),limits=c(min.sils,max.sils),minor_breaks=NULL))
print(ggplot(as.data.frame(numclust), aes(x=seq(1:minclustnum),y=numclust),useDingbats=FALSE) + geom_point() + geom_smooth(span=0.2) + scale_x_continuous(breaks=c(seq(0,minclustnum,10)),limits=c(0,minclustnum),minor_breaks=NULL) + scale_y_continuous(breaks=c(seq(min.clust,max.clust,5)),limits=c(min.clust,max.clust),minor_breaks=NULL))
print(ggplot(as.data.frame(badcellnum),aes(x=seq(1:minclustnum),y=badcellnum),useDingbats=FALSE) + geom_point() + geom_smooth(span=0.2) + scale_x_continuous(breaks=c(seq(0,minclustnum,10)),limits=c(0,minclustnum),minor_breaks=NULL) + scale_y_continuous(breaks=c(seq(min.bad,max.bad,5)),limits=c(min.bad,max.bad),minor_breaks=NULL))
dev.off()

optavgsil=49

s.10x.combined<-FindNeighbors(s.10x.combined, dims = 1:co1, verbose = FALSE,k.param=optavgsil)
s.10x.combined <- FindClusters(s.10x.combined, verbose = FALSE,resolution=0.8)

sil <- silhouette(x = as.numeric(x = as.factor(x = clusters)), dist = dist.matrix)
s.10x.combined$sil <- sil[, 3]
pdf('qc/silhouetteplot_silopt.pdf')
print(fviz_silhouette(sil,label=TRUE,print.summary=TRUE,ggtheme=theme_classic()))
dev.off()

dir.create('umap_10x')
pdf('umap_10x/UMAP_clusters_preRefinement_optsil.pdf',width=5,height=4)
print(DimPlot(s.10x.combined,reduction="umap"))
dev.off()

#starting at cluster ID 1
current=levels(Idents(s.10x.combined))
newID=seq(1,length(unique(Idents(s.10x.combined))))
newID=as.character(newID)
names(newID)=current
s.10x.combined=RenameIdents(s.10x.combined,newID)
s.10x.combined$seurat_clusters=Idents(s.10x.combined)


#merge clusters w/o specific DE gene (might need to tweak)

data.use=s.10x.combined[['pca']][,1:co1]
data.use=t(data.use)
dist.matrix <- dist(x = Embeddings(object = s.10x.combined[['pca']])[, 1:co1])
d=s.10x.combined
d=merge.clusters.DE(object=s.10x.combined,min.de.genes=1,effect.size=2,pval.cutoff=.01,pcs.use=1:co1)
clust.dists=ComputeClusterDistances(d,reduction.use='pca',dist.type='centroid',pcs.use=1:co1)
for (i in 1:length(levels(Idents(d)))){
	if (table(Idents(d))[i]<3) {
		neigh=names(sort(clust.dists[,i],decreasing=TRUE)[2])
		d=SetIdent(object=d,cells.use=which(Idents(d)==i),ident.use=neigh)

}
	
current=levels(Idents(d))
newID=seq(1,length(unique(Idents(d))))
newID=as.character(newID)
names(newID)=current
d=RenameIdents(d,newID)
d$seurat_clusters=Idents(d)
}
pdf('umap_10x/UMAP_clusters_postRefinement.pdf',width=5,height=4)
print(DimPlot(d,reduction='umap'))
dev.off()

#works?
s.10x.combined=d
#else no


genes=c('Snap25','Rbfox3','Slc17a6','Slc17a7','Slc17a8','Gad1','Gad2','Slc32a1','Mbp','Pdgfra','Pecam1','Pdgfrb','Csf1r')
pdf('qc/celltypeViolins.pdf',width=12,height=12)
 print(VlnPlot(s.10x.combined,features=genes,assay='SCT',pt.size=0))
 dev.off()
##############
#subset neuronal and high quality clusters
#############
subs=subset(s.10x.combined,idents=c(1,4,5,6,7,8,9,11,12,13,15,18,19)) 
subs1=subset(subs,orig.ident=='PN1')
subs2=subset(subs,orig.ident=='PN2')
mat=subs1@assays$RNA@counts
nfeatcell=colSums(mat!=0)#get num features for each cell
#calculate point with largest 'jump' in num genes/cell as feature cutoff
hh=c()
num=seq(0.1,99.9,.1)
for (ii in 1:999){
jj=num[ii]
q=quantile(nfeatcell,(jj/100))
hh[ii]=q
}
minfeat=data.frame(x=hh,y=0)
for (a in 1:999){
if (a > 1 ){
y=minfeat$x[a]-minfeat$x[a-1]
minfeat[a,2]=y
}
else {
minfeat[a,2]=0
}
}
#calculate point halfway b/w largest 'jump' in genes/cell

minfeatcutoff=round(minfeat$x[which.max(minfeat$y)])
minfeatcutoff2=round(minfeat$x[which.max(minfeat$y)-1])
minfin=minfeatcutoff-((minfeatcutoff-minfeatcutoff2)/2)

subs1 <- SCTransform(subs1, method = "glmGamPoi", verbose = FALSE,ncells=ncol(subs)/2,variable.features.n=quantile(nfeatcell,.5)) #changing variable features has little effect
subs2 <- SCTransform(subs2, method = "glmGamPoi", verbose = FALSE,ncells=ncol(subs)/2,variable.features.n=quantile(nfeatcell,.5)) #changing variable features has little effect

s.list=c(subs1,subs2)
features <- SelectIntegrationFeatures(object.list = s.list)
s.list<-PrepSCTIntegration(object.list=s.list,anchor.features=features)
subs.anchors=FindIntegrationAnchors(object.list=s.list,normalization.method='SCT',anchor.features = features)
subs=IntegrateData(anchorset=subs.anchors,normalization.method='SCT')
subs <- RunPCA(subs, verbose = FALSE,npcs=100)

subs <- RunPCA(subs, verbose = FALSE,npcs=100)
pdf('qc/Elbow.PNneurons.10x.combined.pdf')
print(ElbowPlot(subs,ndims=100))
dev.off()


#manually determine elbow
co1=12
subs <- RunUMAP(subs, dims = 1:co1, verbose = FALSE)

pdf('qc/PNneurons_UMAP_bysample.pdf',width=10,height=8)
print(DimPlot(subs,group.by='orig.ident'))
dev.off()


minclustnum=round(ncol(subs)/100) #smallest clust

#manually decide number of NNs
# if (ncol(s.10x.combined)>250) {
#      mink=round((ncol(s.10x.combined)/100)+1)
#       } else {
#        mink=3
#        }
#if (ncol(s.10x.combined) > 2000) {
#mink = 20
#}
mink=20
kparam=seq(mink,minclustnum)#what's the minimum cluster size? 10% of the population?
avgsils=rep(NA,minclustnum)
medsils=rep(NA,minclustnum)
badcellnum=rep(NA,minclustnum)
numclust=rep(NA,minclustnum)
dist.matrix <- dist(x = Embeddings(object = subs[['pca']])[, 1:co1])
print('Iterative calculation of silhouettes')
for (k in kparam) {
	#print(paste0('Using ',k,' neighbors'))
	d=subs #create a 'testing' object so you don't ruin anything
	d <- FindNeighbors(d, dims = 1:co1, verbose = FALSE,k.param=k)#first construct a shared nearest neighbors graph. K param indicates how many neighbors to use
	
	d <- FindClusters(d, verbose = FALSE,resolution=0.8)#test resolution for different sil scores?
	numclust[k]=length(unique(Idents(d)))
	#print(paste0('we get ',length(unique(Idents(d))),' clusters'))
	
	#should we optimize n neighbors + resolution by silhouette scores?
		
	clusters=d$seurat_clusters
	sil <- silhouette(x = as.numeric(x = as.factor(x = clusters)), dist = dist.matrix)
	#d$sil <- sil[, 3]
	avgsils[k]=mean(sil[,3])
	medsils[k]=median(sil[,3])
	
	badcells=length(which(sil[,3]<0))
	print(paste0(k,' neighbors; ',length(unique(Idents(d))),' clusters; ',badcells,' bad cells'))
	badcellnum[k]=badcells
}

optavgsil=which.max(avgsils)
optmedsil=which.max(medsils)
print(paste0('opt avg sil = ',optavgsil))


min.sils = round(min(avgsils,na.rm=T),2) - 0.01
max.sils = round(max(avgsils,na.rm=T),2) + 0.01
min.clust = round(min(numclust,na.rm=T),2) - 5
max.clust = round(max(numclust,na.rm=T),2) + 5
min.bad=round(min(badcellnum,na.rm=T),5)-10
max.bad = round(max(badcellnum,na.rm=T),5) + 10

pdf('qc/Sils_optKparams_PNneurons.pdf')
print(ggplot(as.data.frame(avgsils), aes(x=seq(1:minclustnum),y=avgsils),useDingbats=FALSE) + geom_point() + geom_smooth(span=0.5) + scale_x_continuous(breaks=c(seq(0,minclustnum,10)),limits=c(0,minclustnum),minor_breaks=NULL) + scale_y_continuous(breaks=c(seq(min.sils,max.sils,0.01)),limits=c(min.sils,max.sils),minor_breaks=NULL))
print(ggplot(as.data.frame(numclust), aes(x=seq(1:minclustnum),y=numclust),useDingbats=FALSE) + geom_point() + geom_smooth(span=0.2) + scale_x_continuous(breaks=c(seq(0,minclustnum,10)),limits=c(0,minclustnum),minor_breaks=NULL) + scale_y_continuous(breaks=c(seq(min.clust,max.clust,5)),limits=c(min.clust,max.clust),minor_breaks=NULL))
print(ggplot(as.data.frame(badcellnum),aes(x=seq(1:minclustnum),y=badcellnum),useDingbats=FALSE) + geom_point() + geom_smooth(span=0.2) + scale_x_continuous(breaks=c(seq(0,minclustnum,10)),limits=c(0,minclustnum),minor_breaks=NULL) + scale_y_continuous(breaks=c(seq(min.bad,max.bad,5)),limits=c(min.bad,max.bad),minor_breaks=NULL))
dev.off()

optavgsil=36

subs<-FindNeighbors(subs, dims = 1:co1, verbose = FALSE,k.param=optavgsil)
subs <- FindClusters(subs, verbose = FALSE,resolution=0.8)

sil <- silhouette(x = as.numeric(x = as.factor(x = clusters)), dist = dist.matrix)
subs$sil <- sil[, 3]
pdf('qc/neuronssilhouetteplot_silopt.pdf')
print(fviz_silhouette(sil,label=TRUE,print.summary=TRUE,ggtheme=theme_classic()))
dev.off()

pdf('umap_10x/neurons_UMAP_clusters_preRefinement_optsil.pdf',width=5,height=4)
print(DimPlot(subs,reduction="umap"))
dev.off()


current=levels(Idents(subs))
newID=seq(1,length(unique(Idents(subs))))
newID=as.character(newID)
names(newID)=current
subs=RenameIdents(subs,newID)
subs$seurat_clusters=Idents(subs)


#merge clusters w/o specific DE gene (might need to tweak)

data.use=subs[['pca']][,1:co1]
data.use=t(data.use)
dist.matrix <- dist(x = Embeddings(object = subs[['pca']])[, 1:co1])
d=subs
d=merge.clusters.DE(object=subs,min.de.genes=1,effect.size=2,pval.cutoff=.01,pcs.use=1:co1)
clust.dists=ComputeClusterDistances(d,reduction.use='pca',dist.type='centroid',pcs.use=1:co1)
for (i in 1:length(levels(Idents(d)))){
	if (table(Idents(d))[i]<3) {
		neigh=names(sort(clust.dists[,i],decreasing=TRUE)[2])
		d=SetIdent(object=d,cells.use=which(Idents(d)==i),ident.use=neigh)

}

current=levels(Idents(d))
newID=seq(1,length(unique(Idents(d))))
newID=as.character(newID)
names(newID)=current
d=RenameIdents(d,newID)
d$seurat_clusters=Idents(d)
}
pdf('umap_10x/neurons_UMAP_clusters_postRefinement.pdf',width=5,height=4)
print(DimPlot(d,reduction='umap'))
dev.off()

#works?
subs=d
#else no

saveRDS(subs,'thruneuron_clusters.rds')
genes=c('Snap25','Rbfox3','Slc17a6','Slc17a7','Slc17a8','Gad1','Gad2','Slc32a1','Oprd1','Oprm1','Oprk1','Penk','Pdyn','Pomc')
pdf('qc/neuronViolins.pdf',width=12,height=12)
 print(VlnPlot(subs,features=genes,assay='SCT',pt.size=0))
 dev.off()

s.markers<-FindAllMarkers(subs,only.pos=TRUE,min.pct=.25,logfc.threshold=.25)
markername=paste0('PNneurons.10x.AllMarkerGenes.txt')
write.table(s.markers,markername,quote=F,sep='\t',col.names=NA)

markergenes1=s.markers %>% group_by(cluster) %>% top_n(n=1,wt=avg_log2FC) %>% pull()
markergenes3=s.markers %>% group_by(cluster) %>% top_n(n=3,wt=avg_log2FC) %>% pull()
if(length(unique(subs$seurat_clusters))>5){
pdf("umap_10x/neurons_UMAP_Markergenes.pdf",width=length(unique(subs$seurat_clusters)),height=length(unique(subs$seurat_clusters)))
} else{
pdf("umap_10x/neurons_UMAP_Markergenes.pdf")
}
print(FeaturePlot(subs,markergenes1))
dev.off()


pdf('markerdot.pdf',width=12,height=10)
opgenes=c('Oprm1','Oprd1','Oprk1','Penk','Pdyn','Pomc')
print(DotPlot(subs,features=unique(c(opgenes,markergenes3)),assay='SCT',scale=FALSE) + theme(axis.text.x = element_text(angle=90,vjust=.5,hjust=1)))
dev.off()

pdf('fnopdot.pdf',width=12,height=10)
print(DotPlot(subs,features=unique(c(fngenes,opgenes)),assay='SCT') + theme(axis.text.x=element_text(angle=90,vjust=.5,hjust=1)))
dev.off()

mgenes=c()
for (i in 1:length(levels(subs))){
	print(paste0('working on clust ',i))
	j=levels(subs)[i]
	a=markers.binom(subs,clust.1=j,effect.size=.5)
	a=a[order(a$pval),]
	mgenes=c(mgenes,rownames(a)[1:3])
	assign(paste0('markers.',j),a)
}

celltypes=c('Snap25','Slc4a4','Atp1a2','Aqp4','Mag','Mog','Pdgfra','Vcan','Cspg4','Siglech','Csf1r','Cx3cr1','Cst3','Hexb')
pdf('doublets.pdf',width=12,height=10)
print(DotPlot(subs,features=celltypes,assay='SCT',scale=FALSE) + theme(axis.text.x = element_text(angle=90,vjust=.5,hjust=1)))
dev.off()
####Doublet clusters: 8,9,10,11,12,16,19
pdf('violin_qc_clust.pdf')
print(VlnPlot(subs,features=c('nFeature_RNA','nCount_RNA','percent.mt'),pt.size=0))
dev.off()
####Low qc clusters based on nFeat, nCount, and mt:1,4,5,8,10,11,12,16,19ish (def overlap w/ doublet clusters)
#bad = 1,4,5,8,9,10,11,12,16,19


#do we subset by cluster and reanalyze, subset by a minimum nFeat/nCount and reanalyze, or simply discard problem clusters and analyze remaining cells?
#prob w/ 1: 'good' cells in 'bad' clusters.
#prob w/ 2: 'bad' cells w/ high # transcripts (appears uncommon)
#prob w/ 3: not getting entire picture (see 1)
pdf('scatterQC_all.pdf')
print(FeatureScatter(s.10x.combined,feature1='nCount_RNA',feature2='nFeature_RNA'))
dev.off()

pdf('scatterQC_neuro.pdf')
print(FeatureScatter(subs,feature1='nCount_RNA',feature2='nFeature_RNA'))
dev.off()

#neither plots indicate a clear separation of 'good' and 'bad' cells.
#for ease, removing clusters and moving forward may be sufficient

goodonly=subset(subs,idents=c(2,3,6,7,13,14,15,17,18,20))
current=levels(Idents(goodonly))
newID=seq(1,length(unique(Idents(goodonly))))
newID=as.character(newID)
names(newID)=current
goodonly=RenameIdents(goodonly,newID)
goodonly$seurat_clusters=Idents(goodonly)

opgenes=c('Oprm1','Oprd1','Oprk1','Oprl1','Penk','Pdyn','Pomc','Pnoc')
#fngenes=c('Snap25','Rbfox3','Slc17a7','Slc17a6','Slc17a8','Gad1','Gad2','Slc32a1','Chat','Ache','Slc5a7','Slc6a1','Slc6a4','Slc6a3','Slc6a2','Slc6a5')
fngenes=c('Snap25','Slc17a7','Slc17a6','Slc17a8','Slc32a1','Slc6a5')

pdf('fnopdot_goodonly2.pdf',width=7,height=5)
print(DotPlot(goodonly,features=unique(c(fngenes,opgenes)),scale=FALSE,assay='SCT') + theme(axis.text.x=element_text(angle=90,vjust=.5,hjust=1)))
dev.off()
#from ggenes; positive only
posggenes=c()
for (i in 1:length(levels(goodonly))){
	print(paste0('working on clust ',i))
	j=levels(goodonly)[i]
	a=markers.binom(goodonly,clust.1=j,effect.size=.5)
	a=a[which(a$log.effect>0),]
	a=a[order(a$pval),]
	posggenes=c(posggenes,rownames(a)[1:3])
	assign(paste0('goodmarkers.',j,'pos'),a)
}


pdf('fnopmarkdot_goodposonly.pdf',width=10,height=4)
print(DotPlot(goodonly,features=unique(c(fngenes,opgenes,posggenes)),scale=FALSE,assay='SCT') + theme(axis.text.x=element_text(angle=90,vjust=.5,hjust=1)))
dev.off()

finmark=c('Fam107b','Pappa2','Reln','Gria1','Maf','Calb1','Cpne4','Sox2ot','Cadps2','Rorb')
pdf('fnopmarkdot_clean.pdf',width=7,height=4)
print(DotPlot(goodonly,features=unique(c(fngenes,opgenes,finmark)),scale=FALSE,assay='SCT') + theme(axis.text.x=element_text(angle=90,vjust=.5,hjust=1)))
dev.off()

pdf('umap_goodonly.pdf')
print(DimPlot(goodonly,reduction='umap'))
dev.off()


fn2=c('Slc17a7','Slc17a8','Slc32a1','Slc6a5')
pdf('markerdot_binom.pdf',width=16,height=10)
print(DotPlot(goodonly,features=unique(c(fngenes,opgenes,mgenes)),assay='SCT',scale=FALSE) + theme(axis.text.x = element_text(angle=90,vjust=.5,hjust=1)))
dev.off()

save.image('AH_thrumarkers.RData')




########################
#integration of the two 10x datasets + AIBS smartseq data#
########################

#prep
mat=goodonly@assays$RNA@counts
ss=readRDS('/pine/scr/k/y/kylius0/scherrer/Allen_SSData/AllenPons/PONS_PG/PGthruMarkers.rds')
ss=read.csv('PG_counts.csv',header=T,row.names=1)
ss = CreateSeuratObject(ss)
ss@meta.data$orig.ident='SS'
nfeatcell=colSums(mat!=0)#get num features for each cell
#calculate point with largest 'jump' in num genes/cell as feature cutoff
hh=c()
num=seq(0.1,99.9,.1)
for (ii in 1:999){
jj=num[ii]
q=quantile(nfeatcell,(jj/100))
hh[ii]=q
}
minfeat=data.frame(x=hh,y=0)
for (a in 1:999){
if (a > 1 ){
y=minfeat$x[a]-minfeat$x[a-1]
minfeat[a,2]=y
}
else {
minfeat[a,2]=0
}
}
#calculate point halfway b/w largest 'jump' in genes/cell

minfeatcutoff=round(minfeat$x[which.max(minfeat$y)])
minfeatcutoff2=round(minfeat$x[which.max(minfeat$y)-1])
minfin=minfeatcutoff-((minfeatcutoff-minfeatcutoff2)/2)

#first breakout each 10x dataset from 'goodonly' and load AIBS rds (already filtered?); Normalize each separately and integrate
#smartseq variable features should probably differ....

#split 10x and add AIBS smartseq
sparselist<-SplitObject(goodonly,split.by='orig.ident')
#normalize

sparselist=c(sparselist$PN1,sparselist$PN2,ss)



sparselist<- lapply(X=sparselist, FUN = function(x) {
	x <- SCTransform(x, method='glmGamPoi',verbose=FALSE,ncells=ncol(x)/2,variable.features.n=quantile(nfeatcell,.5))
})

features <- SelectIntegrationFeatures(object.list = sparselist)
sparselist<-PrepSCTIntegration(object.list=sparselist,anchor.features=features)
s.int=FindIntegrationAnchors(object.list=sparselist,normalization.method='SCT',anchor.features = features)
s.int=IntegrateData(anchorset=s.int,normalization.method='SCT',features.to.integrate=rownames(all_genes)
s.int <- RunPCA(s.int, npcs=100, verbose=FALSE)

dir.create('AIBS_Janelia10x_combined')
setwd('AIBS_Janelia10x_combined')
dir.create('qc')
pdf('qc/Elbow.PN.int.030722.pdf')
print(ElbowPlot(s.int,ndims=100))
dev.off()

co1=10
co2=16

#Next to UMAP: for smaller datasets, decrease n.neighbors parameter
#To automate, this parameter should be a proportion of cells 
s.int <- RunUMAP(s.int, dims = 1:co1, verbose = FALSE)
s.int <- RunTSNE(s.int, dims=1:co1,verbose=FALSE)
pdf('qc/UMAP_bysample_10PCs_030722.pdf',width=10,height=8)
print(DimPlot(s.int,reduction='umap',group.by='orig.ident'))
dev.off()
pdf('qc/tSNE_bysample_10PCs_030722.pdf',width=10,height=8)
print(DimPlot(s.int,reduction='tsne',group.by='orig.ident'))
dev.off()


s.int <- RunUMAP(s.int, dims = 1:co2, verbose = FALSE)
s.int <- RunTSNE(s.int, dims = 1:co2, verbose=FALSE)
pdf('qc/UMAP_bysample_16PCs.pdf',width=10,height=8)
print(DimPlot(s.int,reduction='umap',group.by='orig.ident'))
dev.off()
pdf('qc/tSNE_bysample_16PCs.pdf',width=10,height=8)
print(DimPlot(s.int,reduction='tsne',group.by='orig.ident'))
dev.off()


minclustnum=round(ncol(s.int)/40) #smallest clust

#manually decide number of NNs
 if (ncol(s.int)>250) {
      mink=round((ncol(s.int)/100)+1)
       } else {
        mink=3
        }
if (ncol(s.int) > 2000) {
mink = 20
}
kparam=seq(mink,minclustnum)#what's the minimum cluster size? 10% of the population?
avgsils=rep(NA,minclustnum)
medsils=rep(NA,minclustnum)
badcellnum=rep(NA,minclustnum)
numclust=rep(NA,minclustnum)
dist.matrix <- dist(x = Embeddings(object = s.int[['pca']])[, 1:co1])
print('Iterative calculation of silhouettes')
for (k in kparam) {
	#print(paste0('Using ',k,' neighbors'))
	d=s.int #create a 'testing' object so you don't ruin anything
	d <- FindNeighbors(d, dims = 1:co1, verbose = FALSE,k.param=k)#first construct a shared nearest neighbors graph. K param indicates how many neighbors to use
	
	d <- FindClusters(d, verbose = FALSE,resolution=0.8)#test resolution for different sil scores?
	numclust[k]=length(unique(Idents(d)))
	#print(paste0('we get ',length(unique(Idents(d))),' clusters'))
	
	#should we optimize n neighbors + resolution by silhouette scores?
		
	clusters=d$seurat_clusters
	sil <- silhouette(x = as.numeric(x = as.factor(x = clusters)), dist = dist.matrix)
	#d$sil <- sil[, 3]
	avgsils[k]=mean(sil[,3])
	medsils[k]=median(sil[,3])
	
	badcells=length(which(sil[,3]<0))
	print(paste0(k,' neighbors; ',length(unique(Idents(d))),' clusters; ',badcells,' bad cells'))
	badcellnum[k]=badcells
}

optavgsil=which.max(avgsils)
optmedsil=which.max(medsils)
print(paste0('opt avg sil = ',optavgsil))


min.sils = round(min(avgsils,na.rm=T),2) - 0.01
max.sils = round(max(avgsils,na.rm=T),2) + 0.01
min.clust = round(min(numclust,na.rm=T),2) - 5
max.clust = round(max(numclust,na.rm=T),2) + 5
min.bad=round(min(badcellnum,na.rm=T),5)-10
max.bad = round(max(badcellnum,na.rm=T),5) + 10

pdf(paste0(qc,'/Sils_optKparams_int.pdf'))
print(ggplot(as.data.frame(avgsils), aes(x=seq(1:minclustnum),y=avgsils),useDingbats=FALSE) + geom_point() + geom_smooth(span=0.5) + scale_x_continuous(breaks=c(seq(0,minclustnum,10)),limits=c(0,minclustnum),minor_breaks=NULL) + scale_y_continuous(breaks=c(seq(min.sils,max.sils,0.01)),limits=c(min.sils,max.sils),minor_breaks=NULL))
print(ggplot(as.data.frame(numclust), aes(x=seq(1:minclustnum),y=numclust),useDingbats=FALSE) + geom_point() + geom_smooth(span=0.2) + scale_x_continuous(breaks=c(seq(0,minclustnum,10)),limits=c(0,minclustnum),minor_breaks=NULL) + scale_y_continuous(breaks=c(seq(min.clust,max.clust,5)),limits=c(min.clust,max.clust),minor_breaks=NULL))
print(ggplot(as.data.frame(badcellnum),aes(x=seq(1:minclustnum),y=badcellnum),useDingbats=FALSE) + geom_point() + geom_smooth(span=0.2) + scale_x_continuous(breaks=c(seq(0,minclustnum,10)),limits=c(0,minclustnum),minor_breaks=NULL) + scale_y_continuous(breaks=c(seq(min.bad,max.bad,5)),limits=c(min.bad,max.bad),minor_breaks=NULL))
dev.off()

optavgsil=68


s.int<-FindNeighbors(s.int, dims = 1:co1, verbose = FALSE,k.param=optavgsil)
s.int <- FindClusters(s.int, verbose = FALSE,resolution=0.8)
clusters=Idents(s.int)
sil <- silhouette(x = as.numeric(x = as.factor(x = clusters)), dist = dist.matrix)
s.int$sil <- sil[, 3]
pdf('qc/silhouetteplot_silopt.pdf')
print(fviz_silhouette(sil,label=TRUE,print.summary=TRUE,ggtheme=theme_classic()))
dev.off()
s.int <- RunUMAP(s.int, dims = 1:co1, verbose = FALSE)
s.int <- RunTSNE(s.int, dims=1:co1,verbose=FALSE)
dir.create('umap_int')
pdf('umap_int/UMAP_clusters_preRefinement_optsil.pdf',width=5,height=4)
print(DimPlot(s.int,reduction="umap"))
dev.off()

#starting at cluster ID 1
current=levels(Idents(s.int))
newID=seq(1,length(unique(Idents(s.int))))
newID=as.character(newID)
names(newID)=current
s.int=RenameIdents(s.int,newID)
s.int$seurat_clusters=Idents(s.int)


#merge clusters w/o specific DE gene (might need to tweak)

data.use=s.int[['pca']][,1:co1]
data.use=t(data.use)
dist.matrix <- dist(x = Embeddings(object = s.int[['pca']])[, 1:co1])
d=s.int
d=merge.clusters.DE(object=s.int,min.de.genes=1,effect.size=2,pval.cutoff=.01,pcs.use=1:co1) #good to look for DE genes, spotty on merging
clust.dists=ComputeClusterDistances(d,reduction.use='pca',dist.type='centroid',pcs.use=1:co1)
for (i in 1:length(levels(Idents(d)))){
	if (table(Idents(d))[i]<3) {
		neigh=names(sort(clust.dists[,i],decreasing=TRUE)[2])
		d=SetIdent(object=d,cells.use=which(Idents(d)==i),ident.use=neigh)

}
	
current=levels(Idents(d))
newID=seq(1,length(unique(Idents(d))))
newID=as.character(newID)
names(newID)=current
d=RenameIdents(d,newID)
d$seurat_clusters=Idents(d)
}
d=RunUMAP(d,dims=1:co1,
pdf('umap_int/UMAP_clusters_postRefinement030722.pdf',width=5,height=4)
print(DimPlot(d,reduction='umap',pt.size=1))
dev.off()

#works?
s.int=d #switch to d... both clustering approaches look decent, 
#else no

fnopgenes=c('Snap25','Slc17a7','Slc17a6','Slc6a5','Gad2','Slc32a1','Oprd1','Oprm1','Oprk1','Oprl1','Penk','Pdyn','Pomc','Pnoc')
pdf('opfnDotplot_merge_int.pdf',width=12)
print(DotPlot(d, features=fnopgenes,scale=FALSE,assay='SCT') + theme(axis.text.x = element_text(angle=90,vjust=.5,hjust=1)))
dev.off()
pdf('opfnDotplot_int.pdf',width=12)
print(DotPlot(s.int, features=fnopgenes,scale=FALSE,assay='SCT') + theme(axis.text.x = element_text(angle=90,vjust=.5,hjust=1)))
dev.off()

#split data by SS vs 10x
d@meta.data$method='10X'
d@meta.data[which(d@meta.data$orig.ident=='SS'),'method']='SS'
s.int@meta.data$method='10X'
s.int@meta.data[which(s.int@meta.data$orig.ident=='SS'),'method']='SS'

pdf('opfnDotplot_merge_int_split.pdf',width=12)
print(DotPlot(d, features=fnopgenes,scale=FALSE,assay='SCT',split.by='method') + theme(axis.text.x = element_text(angle=90,vjust=.5,hjust=1)))
dev.off()
pdf('opfnDotplot_int_split.pdf',width=12)
print(DotPlot(s.int, features=fnopgenes,scale=FALSE,assay='SCT',split.by='method') + theme(axis.text.x = element_text(angle=90,vjust=.5,hjust=1)))
dev.off()

pdf('umap_int/UMAP_clusters_postRefinement_method.pdf',width=5,height=4)
print(DimPlot(d,reduction='umap',group.by='method'))
dev.off()

pdf('VlnFnOp_merge_method.pdf',width=12,height=12)
print(VlnPlot(d,features=fnopgenes,split.by='method',assay='SCT',pt.size=0))
dev.off()
precgenes=c('Slc17a7','Slc32a1','Oprd1','Penk')
DefaultAssay(d)='SCT'
pdf('umap_int/UMAP_clusters_opfn.pdf',width=8,height=13)
print(FeaturePlot(d,features=precgenes,split.by='method',pt.size=1,cols=c('grey','blue')))
dev.off()
pdf('umap_int/UMAP_clusters_postRefinement_splitmethod.pdf',width=8,height=4)
print(DimPlot(d,reduction='umap',split.by='method'))
dev.off()



save.image('AIBS_Janelia10x_combined.RData')


#PLOTTING Oprd1/Slc17a7 CLUST-SPEC COEXPRESSION BY EULER VENN DIAGRAMS
library(eulerr)
library(grid)
library(gridExtra)

dir.create('PN_Venns')
#Subset by cluster
#get Oprm1 onlys, Oprm1UTR onlys, and Oprm1/Oprm1UTR colocalizations; should be better than upsets for only 2 genes
#also, generate a grid of plots (p1,p2 etc.) with no labels + titles?
for (f in 1:length(levels(d))){
print(paste0('Working on cluster ',f))
clust=levels(d)[f]
subs=subset(d,subset=seurat_clusters==clust)
genes=c('Slc17a7','Oprd1')
mat<-subs@assays$RNA@counts[genes,]
print(dim(mat))
opdf=data.frame(row.names=colnames(mat),a=rep(TRUE,ncol(mat)))
colnames(opdf)='Slc17a7-; Oprd1-'
opdf$Slc17a7=mat['Slc17a7',]>0
opdf$Oprd1=mat['Oprd1',]>0
plot.new()
print(dim(opdf))
set=euler(opdf,shape='ellipse',quantities=TRUE)
p=plot(set,labels=NULL,quantities=TRUE)
pdf(paste0('PN_Venns/',clust,'_Oprd1_Slc17a7.pdf'),width=7,height=4)
assign(x=paste0('p',f),value=p)
print((plot(euler(opdf,shape='ellipse',mar=c(5,5,5,5)),quantities=TRUE)))
grid.text(paste0(clust,' ',nrow(opdf),' cells'),x = 0.2, y=0.95, gp=gpar(fontsize=12))
dev.off()
}

pdf('PN_Venns/Eulergrid_Oprd1_Slc17a7_grid.pdf')
grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10, nrow=3)
dev.off()


genes=c('Slc17a7','Slc17a6','Oprd1')
mat<-d@assays$RNA@counts[genes,]
opdf=data.frame(row.names=colnames(mat),a=rep(TRUE,ncol(mat)))
colnames(opdf)='Remaining Pn'
opdf$Slc17a7=mat['Slc17a7',]>0
opdf$Slc17a6=mat['Slc17a6',]>0
opdf$Oprd1=mat['Oprd1',]>0
#opdf$cluster=d@meta.data$seurat_clusters
print(dim(opdf))
pdf('PN_Venns/PNall_Slc17a7Slc17a6Oprd1.pdf')
print(plot(euler(opdf,shape='ellipse'),quantities=TRUE))
dev.off()


genes=c('Slc17a7','Penk','Oprd1')
mat<-d@assays$RNA@counts[genes,]
opdf=data.frame(row.names=colnames(mat),a=rep(TRUE,ncol(mat)))
colnames(opdf)='Remaining Pn'
opdf$Slc17a7=mat['Slc17a7',]>0
opdf$Penk=mat['Penk',]>0
opdf$Oprd1=mat['Oprd1',]>0
#opdf$cluster=d@meta.data$seurat_clusters
print(dim(opdf))
pdf('PN_Venns/PNall_slc17a7penkoprd1.pdf')
print(plot(euler(opdf,shape='ellipse'),quantities=TRUE))
dev.off()

genes=c('Slc17a7','Slc17a6','Gad2','Penk','Oprd1')
mat<-d@assays$RNA@counts[genes,]
opdf=data.frame(row.names=colnames(mat),a=rep(TRUE,ncol(mat)))
colnames(opdf)='Remaining Pn'
opdf$Slc17a7=mat['Slc17a7',]>0
opdf$Slc17a6=mat['Slc17a6',]>0
opdf$Gad2=mat['Gad2',]>0
opdf$Penk=mat['Penk',]>0
opdf$Oprd1=mat['Oprd1',]>0
#opdf$cluster=d@meta.data$seurat_clusters
print(dim(opdf))
pdf('PN_Venns/PNall_all.pdf')
print(plot(euler(opdf,shape='ellipse'),quantities=TRUE))
dev.off()



#ss spec eulers
ssub=subset(d,orig.ident=='SS')


genes=c('Slc17a7','Penk','Oprd1')
mat<-ssub@assays$RNA@counts[genes,]
opdf=data.frame(row.names=colnames(mat),a=rep(TRUE,ncol(mat)))
colnames(opdf)='Remaining Pn'
opdf$Slc17a7=mat['Slc17a7',]>0
opdf$Penk=mat['Penk',]>0
opdf$Oprd1=mat['Oprd1',]>0
#opdf$cluster=d@meta.data$seurat_clusters
print(dim(opdf))
pdf('PN_Venns/PNss_penk.pdf')
print(plot(euler(opdf,shape='ellipse'),quantities=TRUE))
dev.off()


genes=c('Oprm1','Cplx3')
mat<-s@assays$RNA@counts[genes,]
opdf=data.frame(row.names=colnames(mat),a=rep(TRUE,ncol(mat)))
colnames(opdf)='DoubleNegative'
opdf$Oprm1=mat['Oprm1',]>0
opdf$Cplx3=mat['Cplx3',]>0
pdf('ACA_Oprm1Cplx3.pdf')
print(plot(euler(opdf,shape='ellipse'),quantities=TRUE))
dev.off()



####
#markers and dotplots
####
posggenes=c()
for (i in 1:length(levels(d))){
	print(paste0('working on clust ',i))
	j=levels(d)[i]
	a=markers.binom(d,clust.1=j,effect.size=.5)
	a=a[which(a$log.effect>0),]
	a=a[order(a$pval),]
	posggenes=c(posggenes,rownames(a)[1:3])
	assign(paste0('goodmarkers.',j,'pos'),a)
}
goodgenes=c('Snap25','Slc17a7','Slc17a6','Gad2','Slc32a1','Slc6a5','Oprd1','Oprm1','Oprk1','Oprl1','Penk','Pdyn','Pomc','Pnoc',posggenes)

pdf('posgenesDotplot.pdf',height=3,width=10)
print(DotPlot(d,features=unique(goodgenes), assay='SCT',scale=FALSE) + theme(axis.text.x = element_text(angle=90,vjust=.5,hjust=1)))
dev.off()

bettergenes=c('Snap25','Slc17a7','Slc17a6','Slc32a1','Gad2','Oprd1','Oprm1','Penk','Ntng2','Jcad','Slc7a15','P2rx6','Pappa2','Atp2b4','Adcy2','Nell1','Alcam','Tox','Gm26871','Pcsk5','Col12a1','Ebf2','Maf','Lypd6','Ndnf','Cpne7','Crhbp','Gm38505','Ankfn1','Lmx1a')
pdf('posgenesDotplot_better.pdf',height=3,width=8)
print(DotPlot(d,features=unique(bettergenes), assay='SCT',scale=FALSE) + theme(axis.text.x = element_text(angle=90,vjust=.5,hjust=1)))
dev.off()

 
#allcells upset
	dir.create('Upset')
	
	opgenes=c('Slc17a7','Slc32a1','Oprd1','Penk')
	df=d@assays$RNA@counts[opgenes,]
	df[df>0]=1
	df=t(as.matrix(df))
	df=as.data.frame(df)
	pdf('Upset/Allcells_OpFn.pdf')
	print(upset(df,order.by='freq',nsets=8,keep.order=TRUE))
	grid.text(paste0('PN combined ',nrow(df),' cells'),x = 0.65, y=0.95, gp=gpar(fontsize=20))
	dev.off()
	
#split ss and 10x
d10x=subset(x=d,subset = method=='10X')
dss=subset(x=d,subset = method=='SS')


opgenes=c('Slc17a7','Slc32a1','Oprd1','Penk')
df=dss@assays$RNA@counts[opgenes,]
df[df>0]=1
df=t(as.matrix(df))
df=as.data.frame(df)
pdf('Upset/Splitseq_OpFn.pdf')
print(upset(df,order.by='freq',nsets=4,keep.order=TRUE))
grid.text(paste0('PN SS ',nrow(df),' cells'),x = 0.65, y=0.95, gp=gpar(fontsize=20))
dev.off()

opgenes=c('Slc17a7','Slc32a1','Oprd1','Penk')
df=d10x@assays$RNA@counts[opgenes,]
df[df>0]=1
df=t(as.matrix(df))
df=as.data.frame(df)
pdf('Upset/10x_OpFn.pdf')
print(upset(df,order.by='freq',nsets=4,keep.order=TRUE))
grid.text(paste0('PN 10X ',nrow(df),' cells'),x = 0.65, y=0.95, gp=gpar(fontsize=20))
dev.off()



opgenes=c('Slc17a7','Slc32a1','Pdyn','Pnoc','Penk','Pomc','Oprm1','Oprd1','Oprk1','Oprl1')
df=d@assays$RNA@counts[opgenes,]
df[df>0]=1
df=t(as.matrix(df))
df=as.data.frame(df)
pdf('Upset052922.pdf',width=10)
print(upset(df,order.by='freq',nsets=length(opgenes),keep.order=TRUE))
grid.text(paste0('Pn total 4932 neurons'), x=0.65,y=0.95,gp=gpar(fontsize=20))
dev.off()

pdf('umap_int/umaps_5b.pdf')
print(FeaturePlot(d,features=c('sct_Snap25','sct_Gad2','sct_Slc17a7','sct_Slc32a1')))
dev.off()

pdf('umap_int/umaps_5b2.pdf')
print(FeaturePlot(d,features=c('sct_Snap25','sct_Slc6a5','sct_Slc17a7','sct_Slc32a1')))
dev.off()

pdf('umap_int/umaps_5brna.pdf')
print(FeaturePlot(d,features=c('rna_Snap25','rna_Gad2','rna_Slc17a7','rna_Slc32a1')))
dev.off()



current=levels(Idents(d))
#newID=seq(1,length(unique(Idents(d))))
newID=c('E1','E2','I1','E3','E4','I2','I3','I4','I5','E5')
newID=as.character(newID)
names(newID)=current
d=RenameIdents(d,newID)
d$seurat_clusters=Idents(d)
ordID=newID[c(1,2,4,5,10,3,6,7,8,9)]
Idents(d) <- factor(Idents(d), levels=ordID)

bestgenes=c('Snap25','Slc17a7','Slc17a6','Slc32a1','Gad1','Oprd1','Oprm1','Penk','Ntng2','Pappa2','Alcam','Pcsk5','Col12a1','Ebf2','Maf','Ndnf','Ankfn1','Lmx1a')
bestgenes=c('Snap25','Slc17a7','Slc17a6','Slc32a1','Gad1','Oprd1','Oprm1','Penk','Ntng2','Pappa2','Pcsk5','Col12a1','Lmx1a','Alcam','Ebf2','Maf','Ndnf','Ankfn1')

pdf('Dotplot030802.pdf',height=3,width=7)
print(DotPlot(d,features=unique(bestgenes),assay='SCT',scale=FALSE) + theme(axis.text.x = element_text(angle=90,vjust=.5,hjust=1)))
dev.off()

f=d #save
d=RunUMAP(d,dims=1:co1, n.neighbors=5)
pdf('umap_int/UMAP_clusters_fin030822.pdf',width=15,height=12)
print(DimPlot(d,reduction='umap',pt.size=.1))
dev.off()

pdf('umap_int/Umap_bymethod030822.pdf',width=15,height=12)
print(DimPlot(d, reduction='umap',group.by='method'))
dev.off()

pdf('umap_int/umaps_genes030822.pdf',width=10,height=8)
print(FeaturePlot(d,features=c('sct_Snap25','sct_Gad1','sct_Gad2','sct_Slc32a1','sct_Slc6a5','sct_Slc17a7','sct_Slc17a6','sct_Oprd1','sct_Penk'),min.cutoff=0,max.cutoff=2))
dev.off()

fnopgenes=c('Snap25','Slc17a7','Slc17a6','Slc32a1','Slc6a5','Gad1','Gad2','Oprd1','Oprm1','Oprk1','Oprl1','Penk','Pdyn','Pomc','Pnoc')
pdf('Supplemental_Dotplot_splitmethods.pdf',width=8,height=4)
print(DotPlot(d, features=fnopgenes,scale=FALSE,assay='SCT',split.by='method') + theme(axis.text.x = element_text(angle=90,vjust=.5,hjust=1)))
dev.off()
