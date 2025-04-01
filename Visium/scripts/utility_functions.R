#Wrapper functions used as part of the Tumor Microenvironment example for the CRI-Biocore Spatial Transcriptomics Workshop

#distfreq() takes as input a Giotto object and two sets of barcodes for different types of spots (groups A and B, say). It returns the shortest distance from any spot in A to a spot in B, and vice versa.

distfreq <- function(gObject, group1, group2, dmax, netname='Delaunay_network'){
  group1res = array(0, dim=c(length(group1),dmax+1))
  rownames(group1res) = group1
  colnames(group1res) = c(0:dmax)
  group2res = array(0, dim=c(length(group2),dmax+1))
  rownames(group2res) = group2
  colnames(group2res) = c(0:dmax)

  #Identify group1 and group2 cells in same spot
  group1res[intersect(group1, group2),1] = 1
  group2res[intersect(group1, group2),1] = 1

  group1new = group1
  group2new = group2
  for(i in 1:dmax){
    group1nb = findNetworkNeighbors(gObject, spatial_network_name = netname, source_cell_ids = group1new)
    group2nb = findNetworkNeighbors(gObject, spatial_network_name = netname, source_cell_ids = group2new)
    group1new = group1nb$cell_ID[group1nb$nb_cells!='others']
    group2new = group2nb$cell_ID[group2nb$nb_cells!='others']
    group1res[intersect(group1, group2new),i+1] = 1
    group2res[intersect(group1new, group2),i+1] = 1
  }

  group1dist = (dmax+1) - apply(group1res,1,sum)
  group2dist = (dmax+1) - apply(group2res,1,sum)
  results=list()
  results$group1=group1res
  results$group2=group2res
  results$group1dist = group1dist
  results$group2dist = group2dist

  return(results)

}

#dc1prep() takes in a seurat and giotto object of the same data, along with a collection of spots of interest ("cells"), and adds new features to the seurat object that indicate which spots are within 1, 2, or 3 spaces of the focal spots.

dc1prep <- function(scobject, giotto_object, cells, exclude=2, network = 'knn'){
  if(network=='knn'){
    gnet = 'knn_network'
  }
  if(network=='delaunay'){
    gnet = 'Delaunay_network'
  }
  dc1nb=findNetworkNeighbors(giotto_object, spatial_network_name = gnet, source_cell_ids = cells)
  dc1neighbors=dc1nb$cell_ID[dc1nb$nb_cells!='others']  #keep source + neighbors
  dc1nb_nb=findNetworkNeighbors(giotto_object, spatial_network_name = gnet, source_cell_ids = dc1neighbors)
  dc1nb_neighbors=dc1nb_nb$cell_ID[dc1nb_nb$nb_cells=='neighbor'] #keep new neighbors

  scobject$dc1nb=rep(0,nrow(scobject@meta.data))
  scobject$dc1nb[dc1neighbors]=1    #IMPORTANT: the initial NB group is both the original POS spots + Neighbors
  scobject$dc1nb_nb=rep(0,nrow(scobject@meta.data))
  scobject$dc1nb_nb[dc1nb_neighbors]=1  #second group adds the next neighbor "shell"

  dc1cd8=which((scobject$dc1==1&scobject$cd8==1)|(scobject$dc1nb==1&scobject$cd8==1))   #define CD8+ DC1+ as any spot that is itself CD8+ DC1+ or is CD8+ and neighbor to DC1+
  scobject$dc1cd8=rep(0,nrow(scobject@meta.data))
  scobject$dc1cd8[dc1cd8]=1   #Expands the possible spots from 11 to 74!

  #Next, need to define exclusion zone as any non-positive spot (by dc1cd8) within 3 spaces of a positive spot.
  dc1cd8spots=names(which(scobject$dc1cd8==1))
  dc1cd8nb=findNetworkNeighbors(giotto_object, spatial_network_name = gnet, source_cell_ids = dc1cd8spots)
  dc1cd8nb_1=dc1cd8nb$cell_ID[dc1cd8nb$nb_cells=='neighbor']    #direct neighbor to a CD8+ DC1+
  dc1cd8nb_nb=findNetworkNeighbors(giotto_object, spatial_network_name = gnet, source_cell_ids = dc1cd8nb_1)
  dc1cd8nb_2=dc1cd8nb_nb$cell_ID[dc1cd8nb_nb$nb_cells=='neighbor'&dc1cd8nb$nb_cells!='source'&dc1cd8nb$nb_cells!='both']    #Find neighbors to neighbors but make sure to exclude the original sources
  dc1cd8nb_nb_nb=findNetworkNeighbors(giotto_object, spatial_network_name = gnet, source_cell_ids = dc1cd8nb_2)
  dc1cd8nb_3=dc1cd8nb_nb_nb$cell_ID[dc1cd8nb_nb_nb$nb_cells=='neighbor'&dc1cd8nb$nb_cells!='source'&dc1cd8nb_nb$nb_cells!='source'&dc1cd8nb$nb_cells!='both'&dc1cd8nb_nb$nb_cells!='both']    #Find neighbors to neighbors but make sure to exclude the original sources

  #Note: want distance from the original DC1 not from the inferred CD8+ DC1+ spots. So we define this separately using dc1nb_nb and dc1nb
  dc1nb_2 = names(scobject$dc1nb_nb)[scobject$dc1nb_nb==1]
  dc1nb_nb_nb = findNetworkNeighbors(giotto_object, spatial_network_name = gnet, source_cell_ids = dc1nb_2)

  #Modified to allow different distance options quickly. Default is exclude=3
  if(exclude==3){
    dc1nb_ex = unique(dc1nb_nb_nb$cell_ID[dc1nb_nb_nb$nb_cells=='neighbor'|dc1nb_nb$nb_cells=='neighbor'|dc1nb_nb$nb_cells=='both'|dc1nb$nb_cells=='neighbor'])
  }
  if(exclude==2){
    dc1nb_ex = unique(dc1nb_nb$cell_ID[dc1nb_nb$nb_cells=='neighbor'|dc1nb_nb$nb_cells=='both'|dc1nb$nb_cells=='neighbor'])
  }
  if(exclude==1){
    dc1nb_ex = unique(dc1nb)
  }
  #This will include the original dc1s and some of the dc1cd8s that we want to exclude due to distances.
  keep=unique(c(cells, dc1cd8spots))
  dc1nb_ex=setdiff(dc1nb_ex, keep)
  if(exclude==3){
    scobject$dc1nb_ex=scobject$dc1nb_nb
  }
  if(exclude==2){
    scobject$dc1nb_ex=rep(0,nrow(scobject@meta.data))
  }

  scobject$dc1nb_ex[dc1nb_ex]=1
  scobject$dc1ex_vis=rep('Negative',nrow(scobject@meta.data))
  scobject$dc1ex_vis[dc1cd8spots]='Positive'
  scobject$dc1ex_vis[scobject$dc1nb_ex==1]='Excluded'
  #scobject$dc1ex_vis[scobject$dc1==1&scobject$cd8==0]='Excluded'

  scobject$dc1cd8_vis = rep('Negative',nrow(scobject@meta.data))
  scobject$dc1cd8_vis[dc1cd8spots]='Positive'
  scobject$dc1cd8_vis[dc1cd8nb_1]='Excluded'
  scobject$dc1cd8_vis[dc1cd8nb_2]='Excluded'
  scobject$dc1cd8_vis[dc1cd8nb_3]='Excluded'

  scobject$dc1_vis = rep('dc1-',nrow(scobject@meta.data))
  if(exclude==3){
    scobject$dc1_vis[scobject$dc1nb_nb==1]='Excluded'
  }
  if(exclude==2){
    scobject$dc1_vis[scobject$dc1nb==1]='Excluded'
  }

  scobject$dc1_vis[scobject$dc1nb==1]='DC1 Neighbor'
  scobject$dc1_vis[scobject$dc1==1]='DC1+'

  #Visualize the spot exclusions
  newvis = rep(0,ncol(scobject))
  newvis[scobject$dc1==1]='DC1+'
  newvis[scobject$cd8==1&(scobject$dc1==0|scobject$dc1nb==0)]='DC1-CD8+'
  newvis[scobject$cd8==1&scobject$dc1==1]='DC1+CD8+'
  newvis[scobject$cd8==1&scobject$dc1nb==1&scobject$dc1==0]='DC1+CD8nb'
  newvis[scobject$cd8==0&scobject$dc1==0]='DC1-CD8-'
  #newvis[scobject$dc1ex_vis=='Excluded'&scobject$dc1==0]='Excluded'
  newvis[scobject$dc1ex_vis=='Excluded'&scobject$dc1==0]='Excluded'
  newvis[scobject$dc1ex_vis=='Excluded'&scobject$cd8==1]='Excluded CD8+'
  newvis=factor(newvis, levels=c('DC1+','DC1-CD8+','DC1+CD8+','DC1+CD8nb','DC1-CD8-','Excluded','Excluded CD8+'))
  scobject@meta.data$Legend=newvis

  return(scobject)
}

#Spatial correlation function taking advantage of Giotto's tools
getSpatCorGeneSigValues <- function(gobject, subset_genes, order = 'original'){
  smooth_mat=do_spatial_knn_smoothing(gobject=gobject, subset_genes=subset_genes)
  smooth_mat=t_giotto(smooth_mat)
  cormat=array(0, dim=c(length(subset_genes),length(subset_genes)))
  pmat=cormat
  rownames(cormat)=rownames(pmat)=colnames(cormat)=colnames(pmat)=colnames(smooth_mat)
  for(i in 1:length(subset_genes)){for(j in 1:length(subset_genes)){
    cormat[i,j]=cor(smooth_mat[,i],smooth_mat[,j])
    pmat[i,j]=cor.test(smooth_mat[,i],smooth_mat[,j])$p.value
  }
  }
  cormat[is.na(cormat)]=0
  pmat[is.na(pmat)]=1

  if(order == 'original'){
    cormat = cormat[subset_genes, subset_genes]
    pmat = pmat[subset_genes, subset_genes]
  }

  adjpmat = p.adjust(pmat, method='BH')
  adjpmat = array(adjpmat, dim=dim(pmat))
  colnames(adjpmat) = colnames(pmat)
  rownames(adjpmat) = rownames(pmat)

  results=list()
  results$cormat=cormat
  results$pmat=pmat
  results$adjpmat = adjpmat
  return(results)
}
