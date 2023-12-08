####Binomial test to test probability of difference in prop non-zero values b/w groups
binomcount.test=function(object, cells.1,cells.2, effect.size) {
x=object@assays$RNA@counts#raw counts
            
#Test for enrichments in cluster #1
m=rowSums(x[,cells.2]>0) #Number of cells expressing marker in cluster 2
m1 = m; m1[m==0]=1; # Regularization. Add a pseudo count of 1 to unexpressed markers to avoid false positives
n=rowSums(x[,cells.1]>0) #number of cells expressing marker in cluster 1
#Find the probability of finding n or more cells +ve for marker in cluster 1 given the fraction in cluster 2
pv1 = pbinom(n, length(cells.1), m1/length(cells.2), lower.tail = FALSE) + dbinom(n, length(cells.1), m1/length(cells.2))
            
log_fold_express = log((n+1)*length(cells.2)/((m+1)*length(cells.1)),base=2) #log10 proportion of expressing cells

d1 <- data.frame(log.effect=log_fold_express,pval=pv1)
if(!is.na(effect.size)) {
d1 <- subset(d1, log.effect >= effect.size)
}
d1 <- d1[order(d1$pval,decreasing=FALSE),]
            
#Enrichments in cells.2
n1 = n; n1[n==0]=1; # Regularization.  Add a pseudo count of 1 to unexpressed markers to avoid false positives
#Find the probability of finding m or more cells +ve for marker in cluster 2 given the fraction in cluster 1
pv2 = pbinom(m, length(cells.2), n1/length(cells.1), lower.tail=FALSE) + dbinom(m, length(cells.2), n1/length(cells.1))
d2 <- data.frame(log.effect=log_fold_express,pval=pv2)
if(!is.na(effect.size)) {
d2 <- subset(d2, log.effect <= -effect.size)
}            
d2 <- d2[order(d2$pval,decreasing=FALSE),]

d = rbind(d1, d2);
d = d[order(d$pval, decreasing=FALSE),]
return(d)
} 
		  
		 


###identify marker genes based on presence/absence binomial test (can be done vs single, muliple, or all other clusters.
###Note, use of seurat_clusters as cluster IDs.
markers.binom = function(object, clust.1,clust.2=NULL,effect.size,DEAssaySlot,posFrac=0.25) {
strt=as.numeric(Sys.time())
genes.use=rownames(object)
clust.use=Idents(object)
cells.1=names(clust.use[which(clust.use%in%clust.1)])
            
if (is.null(clust.2)) {
	clust.2="rest"
	cells.2=names(clust.use)
	cells.2=cells.2[!(cells.2%in%cells.1)]
} else {
	cells.2=names(clust.use[which(clust.use%in%clust.2)])
}
            
Count.mat = object@assays$RNA@counts
TPM.mat=DEAssaySlot
result=binomcount.test(object, cells.1, cells.2, effect.size)
if (nrow(result)>1){
#ensure countmat = dataframe; throws error for only 1 DEG
#get prop posFrac per marker gene
posFrac.1 = apply(Count.mat[rownames(result),cells.1],1,function(x) round(sum(x > 0)/length(x),2)) 
posFrac.2 = apply(Count.mat[rownames(result),cells.2],1,function(x) round(sum(x > 0)/length(x),2))
} else if (nrow(result)==1) {
  posFrac.1 = length(which((Count.mat[rownames(result),cells.1]>0)))/length(cells.1)
  posFrac.2 = length(which((Count.mat[rownames(result),cells.2]>0)))/length(cells.2)
}

if (clust.2=="rest"){
genes.include = posFrac.1 >= posFrac
} else{
genes.include = (posFrac.1 >= posFrac) | (posFrac.2 >= posFrac)
}
            
result = result[genes.include,]
result = result[order(abs(result$log.effect), decreasing=TRUE),]
if (nrow(result)>0){
if (nrow(result)>1){
#Mean number of transcripts per cell
nTrans.1 = apply(TPM.mat[rownames(result), cells.1], 1, function(x) round(mean(x),3))
nTrans.2 = apply(TPM.mat[rownames(result), cells.2], 1, function(x) round(mean(x),3))
} else if (nrow(result)==1) {
  nTrans.1=mean(TPM.mat[rownames(result),cells.1])
  nTrans.2=mean(TPM.mat[rownames(result),cells.2])
}
result[,paste0("nTrans_", clust.1)] = nTrans.1
result[, paste0("nTrans_", clust.2)] = nTrans.2
#THIS IF STATEMENT NEEDS TO BE RUN, OR IT WILL FAIL WHEN NO DEGS EXIST BETWEEN CLUSTERS
result$ClusterID=clust.1 
}           
print(as.numeric(Sys.time()) - strt)
return(result)
} 



####determine cluster centroids... This should be easy. they're using PCs.
ClusterCentroids = function(object,reduction.use="pca", pcs.use = 1:10, cells.use=NULL, genes.use=NULL) {
            cells.use = colnames(object)
            group.use=Idents(object)[cells.use]
            if (reduction.use == "pca"){
              data.use = object@reductions$pca@cell.embeddings[cells.use,pcs.use]
            }
            
            centroids = c()
            for (i in levels(group.use)){
              cells.in.cluster = names(Idents(object))[Idents(object) %in% i]
              cells.in.cluster = cells.in.cluster[cells.in.cluster %in% cells.use]
              centroids = rbind(centroids, colMeans(data.use[cells.in.cluster,]))
            }
            centroids = as.data.frame(centroids)
            colnames(centroids) = colnames(data.use)
            rownames(centroids) = as.numeric(levels(Idents(object)))
            
            return(centroids)
          }
          
          


#Compute distances b/w cluster centroids... useful for merging and determining cluster differences
ComputeClusterDistances=function(object,reduction.use="pca",dist.type="nn", pcs.use = 1:10, cells.use=NULL, genes.use=NULL) {
            cells.use =  colnames(object)
            group.use=Idents(object)[cells.use]
            if (reduction.use == "pca"){
              data.use = object@reductions$pca@cell.embeddings[cells.use,pcs.use]
              centroids = ClusterCentroids(object, reduction.use="pca", pcs.use=pcs.use, cells.use=cells.use)
            }
            
            if (dist.type=="centroid"){
              clust.dists = as.matrix(dist(centroids, upper=TRUE))
              diag(clust.dists) = 1e6
            }
            
            num.clust = length(levels(group.use))
            
            
            if (dist.type == "nn"){
              clust.dists = matrix(0, nrow=num.clust, ncol=num.clust)
              diag(clust.dists) = 1e6
              rownames(clust.dists) = levels(group.use)
              colnames(clust.dists) = rownames(clust.dists)
              for (i in 1:nrow(clust.dists)){
                for(j in ((i+1):ncol(clust.dists))){
                  if (j>nrow(clust.dists)) break
                  cells.in.cluster_i = names(Idents(object))[Idents(object) %in% i]
                  cells.in.cluster_i = cells.in.cluster_i[cells.in.cluster_i %in% cells.use]
                  cells.in.cluster_j = names(Idents(object))[Idents(object) %in% j]
                  cells.in.cluster_j = cells.in.cluster_j[cells.in.cluster_j %in% cells.use]
                  
                  nnA = nn2(data.use[cells.in.cluster_i,], query = centroids[j,], k=1)
                  nnB = nn2(data.use[cells.in.cluster_j,], query = centroids[i,],k=1)
                  clust.dists[i,j] = min(c(nnA$nn.dists, nnB$nn.dists))
                  
                  clust.dists[j,i] = clust.dists[i,j]
                }
              }
            }
            
            colnames(clust.dists) = c(1:ncol(clust.dists))
            rownames(clust.dists) = colnames(clust.dists)
            return(clust.dists)
          }
          




 merge.clusters.DE=function(object,min.de.genes, effect.size,pval.cutoff, pcs.use,DEAssaySlot,posFrac=0.25) {
          
            genes.use = rownames(object@assays$RNA@meta.features)
            clust.test = as.numeric(levels(Idents(object)))
           # if (is.null(tag)){
           #   filename = "CLUSTER_PAIRWISE_MARKERS.txt"
           # } else {
           #   filename = paste0("CLUSTER_PAIRWISE_MARKERS_",tag,".txt")
           # }
           # zz = file(filename,open="wt")
            
            
            num.clust=length(clust.test) 
            print(paste0("Starting with ", num.clust, " clusters"))
            
            pass.thresh=1e6*data.frame(diag(length(levels(Idents(object))))); 
            
            for (i in setdiff(as.numeric(levels(Idents(object))), clust.test)){
              pass.thresh[i,]=1e6; pass.thresh[,i]=1e6;
            } 
            
            dist.clust = pass.thresh
            
            #Find the number of differentially expressed genes between every pair of clusters
            for(k in 1:num.clust) {
              i=clust.test[k]
              print(paste0("Testing Cluster ", i))
              for(m in ((k+1):num.clust)) {
                j=clust.test[m]
                print(j)
                if (m>num.clust) break
                if (pass.thresh[i,j]==0) {
                  marker=markers.binom(object=object, clust.1=i,clust.2=j,effect.size=effect.size,DEAssaySlot=DEAssaySlot,posFrac=posFrac)
					        paste0('found markers for',i,' vs ',j)
                 
                  
                  num.de.genes = nrow(marker)
                  pass.thresh[i,j]=num.de.genes; pass.thresh[j,i]=pass.thresh[i,j];
                  
                }
                
              }
              
              print(pass.thresh[i,])
              
            }
            
            colnames(pass.thresh) = levels(Idents(object))
            rownames(pass.thresh) = levels(Idents(object))
            
            write.table(pass.thresh, file=paste0("DE_genes_matrix_2.txt"), sep="\t", quote=FALSE)
            
            #iteratively merge clusters
	    #remove diagonal of effectively zeros; we already did this, it was positive e.
	    for (i in 1:length(pass.thresh)){
		pass.thresh[i,i]=1000
	    }
	    print(pass.thresh)
            min.val = min(pass.thresh) #will this always detect diagonal of 0.000001?; above for loop should fix
	    min.val.ind = as.data.frame(which(pass.thresh==min.val, arr.ind = TRUE))#index of minimum value (row; col; value)
            min.val.ind = min.val.ind[min.val.ind$row < min.val.ind$col,]
            min.val.ind$val = min(pass.thresh)
            rownames(min.val.ind) = 1:nrow(min.val.ind)
            
            merge.ind=-1
            while(min.val <= min.de.genes & length(unique(Idents(object)))>2) {
              merge.ind=merge.ind+1 #merge index becomes 0 
              
		#THIS IS FUCKED UP 021923
		#pass.thresh is filled with zeros; recompute pairwise markers.
		#possibilities: markers.binom throwin zeros= incorrect assignment of test.1?
		#test.1 does not equal new cluster for markers.binom. 
		#alternatively, markers.binom is fine, but rest of loop results in 0 for some reason... marker.pass?
              #In case of ties, merge clusters that are closest in PC space
	      print('New pass.threshold marker matrix')
	      print(pass.thresh)
              print('ComputingClusterDistances')
              clust.dists = ComputeClusterDistances(object, reduction.use="pca", dist.type="centroid", pcs.use=pcs.use) 
              ind.min=which.min(clust.dists[cbind(min.val.ind$row, min.val.ind$col)])#what's the distance bw clusters to be merged? using pass.thresh index (min.val.ind)
              test.1 = min.val.ind[ind.min,]$row; test.2 = min.val.ind[ind.min,]$col #cluster IDs to be tested. Note test.2 would get merged into test.1
              
              if (pass.thresh[test.1,test.2]<= min.de.genes) { #if the number of DEGs b/w to-be-merged clusters x and y are less than min.de.genes...
                Idents(object)[which(Idents(object)==test.2)]=test.1 #cells assigned cluster y are now assigned to cluster x (merged)
                pass.thresh = pass.thresh[-test.2,]; pass.thresh = pass.thresh[,-test.2]# cluster y is removed from pass.thresh matrix
                old.group.levels = as.numeric(levels(Idents(object))) #cluster names w/o merged cluster y
                old.group.levels = setdiff(old.group.levels, test.2) #this doesn't do anything...
                clust.test = setdiff(clust.test, test.2) #equivalent to old.group.levels (2 removed from list of clusters)
                
                Idents(object) = droplevels(Idents(object)) #remove cluster names
                levels(Idents(object)) = c(1:length(levels(Idents(object)))) #rename clusters such that they're sequential 1:x
                object$seurat_clusters = Idents(object) #rename seurat_clusters metadata column
                
                new.group.levels = as.numeric(levels(Idents(object)))
                names(new.group.levels) = as.character(old.group.levels)
                clust.test = new.group.levels[as.character(clust.test)] #current clusters  with old cluster name; e.g. cluster 2 merged w/ 1; values=1,2,3; names=1,3,4
                
                
                
                #Recompute pairwise markers for merged cluster
                print(paste0("Recomputing pairwise markers for new clust ", test.1))
                for (i in setdiff(clust.test, test.1)){ #all current clusters minus the cluster we generated from merging
                  print(i)
                  marker=markers.binom(object,test.1,i,effect.size=effect.size,DEAssaySlot=DEAssaySlot,posFrac=posFrac)
                  #marker= FindMarkers(object,test.1,i,only.pos=TRUE)#use wilcoxon test?
                  P2 = p.adjust(marker$pval, method="fdr")
                  marker$pval = P2
                  marker.pass=subset(marker,pval<pval.cutoff)
                  pass.thresh[test.1,i]=2*min(nrow(subset(marker.pass, log.effect>0)),nrow(subset(marker.pass, log.effect<0))); pass.thresh[i,test.1]=pass.thresh[test.1,i];
                  #pass.thresh[test.1,i]=nrow(marker.pass); 
                  pass.thresh[i,test.1]=pass.thresh[test.1,i];
                  
                }
                
              }
              colnames(pass.thresh) = 1:length(levels(Idents(object)))
              rownames(pass.thresh) = colnames(pass.thresh)
              
              min.val = min(pass.thresh)
              min.val.ind = as.data.frame(which(pass.thresh==min.val, arr.ind = TRUE))
              min.val.ind = min.val.ind[min.val.ind$row < min.val.ind$col,]
              min.val.ind$val = min(pass.thresh)
              rownames(min.val.ind) = 1:nrow(min.val.ind)
              
            }
            return(object)
          }


#Stacked Violins
modify_vlnplot<- function(obj, 
                          feature, 
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  p<- VlnPlot(obj, features = feature, pt.size = pt.size, ... )  + 
    xlab("") + ylab(feature) + ggtitle("") + 
    theme(legend.position = "none", 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.title.y = element_text(size = rel(1), angle = 0), 
          axis.text.y = element_text(size = rel(1)), 
          plot.margin = plot.margin )
  return(p)
}

## extract the max value of the y axis
extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}


## main function
StackedVlnPlot<- function(obj, features,
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  
  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(angle = 45, hjust=1,vjust=1), axis.ticks.x = element_line()) 

  
  # change the y-axis tick to only max value 
  ymaxs<- purrr::map_dbl(plot_list, extract_max)
  plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x + 
                            scale_y_continuous(breaks = c(y)) + 
                            expand_limits(y = y))

  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}



###identify marker genes based on presence/absence binomial test (can be done vs single, muliple, or all other clusters.
###Note, use of seurat_clusters as cluster IDs.
markers.binom.metagroup = function(object,metagroup, clust.1,clust.2=NULL,effect.size,DEAssaySlot,posFrac=0.25) {
genes.use=rownames(object)
clust.use=unlist(object[[metagroup]])
names(clust.use)=colnames(object)
cells.1=names(clust.use[which(clust.use%in%clust.1)])

if (is.null(clust.2)) {
        clust.2="rest"
        cells.2=names(clust.use)
        cells.2=cells.2[!(cells.2%in%cells.1)]
} else {
        cells.2=names(clust.use[which(clust.use%in%clust.2)])
}

Count.mat = object@assays$RNA@counts
TPM.mat=DEAssaySlot
result=binomcount.test(object, cells.1, cells.2, effect.size)


posFrac.1 = apply(Count.mat[rownames(result),cells.1],1,function(x) round(sum(x > 0)/length(x),2))
posFrac.2 = apply(Count.mat[rownames(result),cells.2],1,function(x) round(sum(x > 0)/length(x),2))

if (clust.2=="rest"){
genes.include = posFrac.1 >= posFrac
} else{
genes.include = (posFrac.1 >= posFrac) | (posFrac.2 >= posFrac)
}

result = result[genes.include,]
result = result[order(abs(result$log.effect), decreasing=TRUE),]

#Mean number of transcripts per cell
if (nrow(result)>0){
nTrans.1 = apply(TPM.mat[rownames(result), cells.1], 1, function(x) round(mean(x),3))
nTrans.2 = apply(TPM.mat[rownames(result), cells.2], 1, function(x) round(mean(x),3))
result[,paste0("nTrans_", clust.1)] = nTrans.1
result[, paste0("nTrans_", clust.2)] = nTrans.2
#THIS IF STATEMENT NEEDS TO BE RUN, OR IT WILL FAIL WHEN NO DEGS EXIST BETWEEN CLUSTERS
result$ClusterID=clust.1
}
return(result)
}


