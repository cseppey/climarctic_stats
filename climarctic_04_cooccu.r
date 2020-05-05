#####
# climarctic cooccu
#####

print('##### Climarctic 04 coocurrence networks #####')

rm(list=ls())

require(doSNOW)
require(foreach)
require(igraph)
require(abind)
require(plotrix)

source('~/bin/src/my_prog/R/pie_taxo.r')
source('~/bin/src/my_prog/R/legend_pie_taxo.r')

# prep cluster
cl <- makeSOCKcluster(4)

registerDoSNOW(cl)

# dir loads ####
dir_out  <- 'Projets/Climarctic/stats/MBC/out/'
dir_save <- paste0(dir_out, 'saves/') 
dir_cooc  <- paste0(dir_out, '04_cooc/')
dir.create(dir_cooc, showWarnings=F)

#---
file <- paste0(dir_save, '00_lst_comm.Rdata')
load(file)
#





#--------------------------------------------------------------
# functions
#--------------------------------------------------------------
#see http://psbweb05.psb.ugent.be/conet/documents/networksFromSimCounts.R


compute.kld=function(x, pseudocount=0.00000001){
  # diagonal is zero
  kld=matrix(data=0,nrow=nrow(x),ncol=nrow(x))
  for(i in 1:nrow(x)){
    for(j in 1:i){
      kld[i,j]=get.kld(x[i,],x[j,], pseudocount=pseudocount)   
      kld[j,i]=kld[i,j]  
    }
  }
  kld
}

get.kld=function(x,y, pseudocount=0.00000001){
  if(length(x) != length(y)){
    stop("The two vectors should have the same length!")
  }
  x[x==0]=pseudocount
  y[y==0]=pseudocount
  dis = 0
  x = x/sum(x)
  y = y/sum(y)
  for(i in 1:length(x)){
    if(!is.nan(x[i]) && !is.nan(y[i])){
      ratioxy = log(x[i]/y[i])
      ratioyx = log(y[i]/x[i])
      dis = x[i]*ratioxy+y[i]*ratioyx + dis
    }
  }
  dis
}

get.pval = function(matrix, x.index, y.index, N.rand=1000, method="spearman", renorm=F, permutandboot=F, plot=F, verbose=F) {
  x = matrix[x.index,]
  y = matrix[y.index,]
  lower.tail = TRUE
  # bray and kld are dissimilarities, so one-sided p-value needs to be computed from the upper tail
  if(method == "bray" || method == "kld"){
    lower.tail = FALSE
  }
  if(method == "spearman"){
    this.sim = cor(x, y, use="complete.obs", method="spearman")
  }else if(method == "pearson"){
    this.sim = cor(x, y, use="complete.obs", method="pearson")
  }else if(method == "bray"){
    this.sim= vegdist(rbind(x,y),method="bray")
  }else if(method == "kld"){
    this.sim=get.kld(x,y)
  }else{
    stop("Select either spearman, pearson, kld or bray as method.")
  }
  rand.sim = rep(NA, N.rand)
  boot.sim = rep(NA, N.rand)
  for (i in 1:N.rand) {
    rand = sample(x, length(x))
    if(renorm == T){
      mat.copy=matrix
      mat.copy[x.index,]=rand
      mat.copy = normalize(mat.copy)
      rand = mat.copy[x.index,]
      y = mat.copy[y.index,]
    }
    if(method == "spearman"){
      rand.sim[i] = cor(rand, y, method="spearman", use="complete.obs")
    }else if(method == "pearson"){
      rand.sim[i] = cor(rand, y, method="pearson",use="complete.obs")
    }else if(method == "bray"){
      rand.sim[i] = vegdist(rbind(rand,y),method="bray") 
    }else if(method == "kld"){
      rand.sim[i] = get.kld(rand,y)
    }
  }
  rand.sim = na.omit(rand.sim)
  if(plot == T){
    col1=rgb(0,0,1,1/3)
    col2=rgb(1,0,0,1/3)
    hist(rand.sim,col=col1)
    abline(v=mean(rand.sim),col="blue")
  }
  if(permutandboot){
    x=matrix[x.index,]
    y=matrix[y.index,]
    for (i in 1:N.rand) {
      rand.idx = sample(1:length(x),replace=TRUE)
      x.boot=x[rand.idx]
      y.boot=y[rand.idx]
      if(method == "spearman"){
        boot.sim[i] = cor(x.boot, y.boot, method="spearman", use="complete.obs")
      }else if(method == "pearson"){
        boot.sim[i] = cor(x.boot, y.boot, method="pearson",use="complete.obs")
      }else if(method == "bray"){
        boot.sim[i] = vegdist(rbind(x.boot,y.boot),method="bray") 
      }else if(method == "kld"){
        boot.sim[i] = get.kld(x.boot,y.boot)
      }
    }
    boot.sim = na.omit(boot.sim)
    if(plot == T){
      hist(boot.sim,col=col2,add=T)
      abline(v=mean(boot.sim),col="red")
      legend(x="topleft", c("Permut","Boot"), bg="white",col=c(col1,col2),lty=rep(1,2),merge=T)
    }
    # if we got enough non-NA permutation and bootstrap values, compute p-value
    if(length(rand.sim) > round(N.rand/3) && length(boot.sim) > round(N.rand/3)){
      pval = pnorm(mean(rand.sim),mean=mean(boot.sim),sd=sd(boot.sim), lower.tail=lower.tail)
    }else{
      pval = 0.5
    }
  }else{
    # if we got enough non-NA permutation values, compute p-value
    if(length(rand.sim) > round(N.rand/3)){
      if (lower.tail) {
        pval = (sum(this.sim > rand.sim) / length(rand.sim))
      } else {
        pval = (sum(this.sim < rand.sim) / length(rand.sim))
      }
    }else{
      pval = 0.5
    }
  }
  # set missing value (from constant vector) to intermediate p-value (worst possible p-value in this context)
  if(is.na(pval)){
    pval = 0.5
  }
  # p-values are one-sided, so high p-values signal mutual exclusion and are converted into low ones
  if(pval > 0.5){
    pval = 1 - pval
  }
  if(verbose == T){
    print(paste("p-value =",pval))
    print(paste("original score",this.sim))
    print(paste("mean of null distrib",mean(rand.sim)))
    print(paste("sd of null distrib",sd(rand.sim)))
    if(permutandboot == T){
      print(paste("mean of boot distrib",mean(boot.sim)))
      print(paste("sd of boot distrib",sd(boot.sim)))
    }
  }
  pval
}

#--------------------------------------------------------------
# initiation
#--------------------------------------------------------------
# library(vegan)
# library(ggplot2)
# library(cooccur)
# library(igraph)
# library(flexmix)
# 
# 
# #load data
# #ITS
# source('D://Msweetlo/Documents/POLAR_FUNGI/analyses2/polarfungi2_LoadData.R')
# phyla<-c('Rozellomycota', 'Ascomycota', 'Basidiomycota', 'Chytridiomycota')
# phyla2<-c('Fungi', 'Rozellomycota', 'Ascomycota', 'Basidiomycota', 'Chytridiomycota')
# 
# FOTUtab_prab<-FOTUtab
# FOTUtab_prab[FOTUtab_prab>1]<-1
# 
# #18S UPARSE
# crange=97;euk=TRUE;bac=FALSE;OriginalData=TRUE;prab=TRUE
# source('D://Msweetlo/Documents/GREENLAND_EXT/Analysis/bipol5_LoadData.R')
# 
# Efdata<-Edata97[Edata97$phylum=='Fungi',fnames_allsamples]
# Efdata<-Efdata[rowSums(Efdata)>0,]
# EfOTUtab<-data.frame(t(Efdata))
# Efdata<-Edata97[colnames(EfOTUtab), c('class', 'order', 'family', 'genus', 'species', fnames_allsamples)]
# 
# EfOTUtab_prab<-EfOTUtab
# EfOTUtab_prab[EfOTUtab_prab>1]<-1
# 
# fnames2<-fnames_allsamples
# 
# 
# #18S SWARM
# source('D://Msweetlo/Documents/GREENLAND_EXT/Analysis_Swarm/bipol6_LoadSwarm.R')
# setwd('D://Msweetlo/Documents/GREENLAND_EXT/Analysis_Swarm')
# 
# fnames2<-intersect(fnames_allsamples, enames_allsamples)
# 
# Efdata<-Edata97[Edata97$phylum=='Fungi',fnames2]
# Efdata<-Efdata[rowSums(Efdata)>0,]
# EfOTUtab<-data.frame(t(Efdata))
# Efdata<-Edata97[colnames(EfOTUtab), c('class', 'order', 'family', 'genus', 'species', fnames2)]
# 
# EfOTUtab_prab<-EfOTUtab
# EfOTUtab_prab[EfOTUtab_prab>1]<-1
# 
# 
# FOTUtab_prab<-FOTUtab[fnames2,]
# FOTUtab_prab[FOTUtab_prab>1]<-1
# 
# 
# netmatrix<-read.csv2(file="cor_matrix.csv", header=TRUE, dec=",", sep=";", row.names=1)
# 
# #--------------------------------------------------------------
# setwd('D://Msweetlo/Documents/POLAR_FUNGI/analyses2')
# #--------------------------------------------------------------
# # Co Occurence alanysis ITS 18S all taxa
# #--------------------------------------------------------------
# fnames2<-intersect(fnames_allsamples, enames_allsamples)
# fnames2<-intersect(c(fnames_ant, fnames_pen), enames_allsamples)
# fnames2<-intersect(fnames_north, enames_allsamples)
# 
# Sys.time()
# EOTUtab_r<-rrarefy(EOTUtab[fnames2,], 5000)
# EOTUtab_r<-EOTUtab_r[,colSums(EOTUtab_r)>0]
# Edata_r<- data.frame(Edata97[colnames(EOTUtab_r),c('phylum','class', 'order', 'family', 'genus', 'species')], data.frame(t(EOTUtab_r)))
# # 
# Efdata_r<-Edata_r#[Edata_r$phylum=='Fungi',]
# EfOTUtab_r<-data.frame(t(Efdata_r[,fnames2]))
# EfOTUtab_r<-EfOTUtab_r[,colSums(EfOTUtab_r)>0]
# Efdata_r<- data.frame(Efdata_r[colnames(EfOTUtab_r),c('class', 'order', 'family', 'genus', 'species')], data.frame(t(EfOTUtab_r)))
# 
# FOTUtab_r<-rrarefy(FOTUtab[fnames2,],5000)

#%%%%%%%%%%
pairs <- data.frame(combn(n_comm, 2))

for(h in pairs){
  hc <- as.character(h)  
  print(hc)
  
  env_1 <- lst_comm[[hc[1]]]$raw$env
  mr_1 <-  lst_comm[[hc[1]]]$raw$mr[env_1$low_seq == F,]
  
  env_2 <- lst_comm[[hc[2]]]$raw$env
  mr_2 <-  lst_comm[[hc[2]]]$raw$mr[env_2$low_seq == F,]
  
  smps_inter <- intersect(row.names(mr_1), row.names(mr_2))
  
  # I would favour relative abundance as it loose much less information
  # a rarefaction threshold at 5000 sequences would lose 80% of the eukaryotic samples
  mr_1_relabu <- decostand(mr_1[smps_inter,], 'total')
  mr_2_relabu <- decostand(mr_2[smps_inter,], 'total')
  
  #%%%%%%%%%%
  
  indata_1<-mr_1_relabu
  indata_2<-mr_2_relabu
  
  data_comb<-cbind(indata_2, indata_1)
  # data_comb<-data_comb[,colSums(data_comb>0) >= 5]#otus should be present in at least 5 samples
  data_comb<-data_comb[,colSums(data_comb > 0) >= nrow(data_comb)/5]
  indata_1<-indata_1[,intersect(colnames(indata_1), colnames(data_comb))]
  indata_2<-indata_2[,intersect(colnames(indata_2), colnames(data_comb))]
  
  
  #%%%%%%%%%%
  rng <- 1:100 #700 seq for first loops 1230 seq for snd loop
  
  # calculate network matrix only consider interactions between datasets, not within
  netmatrix <- data.frame(matrix(NA,nrow=ncol(indata_1),ncol=ncol(indata_2)))
  rownames(netmatrix)<-colnames(indata_1)
  colnames(netmatrix)<-colnames(indata_2)
  p_matrix <- netmatrix #save p-values for false discovery analysis
  kldmatrix <- netmatrix
  p_matrxkld <- netmatrix
  
  #cpute correlation ensembles
  file <- paste0(dir_save, '04_arr_cor_', paste(h, collapse='_'), '.Rdata')
  if(file.exists(file)){
    load(file)
  } else {
    system.time(for (sp1 in colnames(indata_1)){
      print(sp1)
      for (sp2 in colnames(indata_2)){
        #spearman
        cor_out<-cor.test(data_comb[,sp1], data_comb[,sp2], method="spearman") 
        p_matrix[sp1, sp2]<-cor_out$p.value 
        
        #Kullback-Leibler distance
        y<-t(as.matrix(data.frame(spp1=data_comb[,sp1]/sum(data_comb[,sp1]),spp2=data_comb[,sp2]/sum(data_comb[,sp2]))))
        kld_value<-compute.kld(y)[1,2]
        
        set.seed(0)
        pval_kld<-as.numeric(get.pval(y, 1, 2, N.rand=1000, permutandboot=F, method="kld"))
        p_matrxkld[sp1, sp2]<-pval_kld
        
        if (cor_out$p.value > 0.05 | pval_kld > 0.05){
          netmatrix[sp1,sp2] <- as.numeric(0)
          kldmatrix[sp1,sp2] <- as.numeric(0)
        }
        else{
          if (sp1!=sp2){
            netmatrix[sp1,sp2] <- cor_out$estimate
            kldmatrix[sp1,sp2] <- kld_value
            
          }
          else{
            netmatrix[sp1,sp2] <- as.numeric(0)
            kldmatrix[sp1,sp2] <- as.numeric(0)
          }
        }
      }
    })
    arr_cor <- abind(netmatrix=netmatrix, kldmatrix=kldmatrix, p_matrix=p_matrix, p_matrxkld=p_matrxkld, along=3)
    
    save(arr_cor, file=file)
  }
  
  netmatrix <- arr_cor[,,'netmatrix']
  p_matrxkld <- arr_cor[,,'p_matrxkld']
  p_matrix <- arr_cor[,,'p_matrix']
  
  #correct p-values and remove new non-significants (Benjamini-hochberg)
  netmatrix2<-netmatrix
  pvalsadj2<-p.adjust(unlist(p_matrxkld), 'BH')   
  pvalsadj2<-matrix(nrow=nrow(p_matrxkld), ncol=ncol(p_matrxkld), data=pvalsadj2)
  for(i in 1:nrow(pvalsadj2)){
    for(j in 1: ncol(pvalsadj2)){
      if(pvalsadj2[i,j]>0.05){
        netmatrix2[i,j]<-as.numeric(0)
      }
    }
  }
  
  pvalsadj2<-p.adjust(unlist(p_matrix), 'BH')   
  pvalsadj2<-matrix(nrow=nrow(p_matrix), ncol=ncol(p_matrix), data=pvalsadj2)
  for(i in 1:nrow(pvalsadj2)){
    for(j in 1: ncol(pvalsadj2)){
      if(pvalsadj2[i,j]>0.05){
        netmatrix2[i,j]<-as.numeric(0)
      }
    }
  }
  
  arr_cor <- abind(arr_cor, netmatrix2=netmatrix2, along=3)
  
  #polarfungi3_Network_ITS18S_uparseAll.csv
  #polarfungi3_Network_ITS18S_uparseAntarctica.csv
  #polarfungi3_Network_ITS18S_swarmAll.csv
  #polarfungi3_Network_ITS18S_swarmAntarctica.csv
  # write.csv2(netmatrix, "polarfungi3_Network_ITS18S_swarmAntarctica.csv")
  #netmatrix<-read.csv2(file="cor_matrix.csv", header=TRUE, dec=",", sep=";", row.names=1)
  
  #remove OTUs that have no edges 
  #%%%%%%%%%%% why not using netmatrix2?
  # # netmatrix_parse<-netmatrix[,colSums(netmatrix)>0]
  # netmatrix_parse<-netmatrix2[,colSums(netmatrix2)>0] #%%%%%%%%%%% mistake as there can be some OTU with positive correlation while the sum is negative
  # netmatrix_parse<-netmatrix_parse[rowSums(netmatrix_parse)>0,]
  # #change matrix to presence of an edege (=1) vs. no edge (=0)
  # netmatrix_parse[netmatrix_parse >= 0.60]<-1
  # netmatrix_parse[netmatrix_parse < 0.60]<-0
  # 
  # netmatrix_parse<-netmatrix_parse[,colSums(netmatrix_parse)>0]
  # netmatrix_parse<-netmatrix_parse[rowSums(netmatrix_parse)>0,]
  
  netmatrix_parse <- ifelse(netmatrix2 > 0.6, 1, 0)
  netmatrix_parse <- netmatrix_parse[rowSums(netmatrix_parse) != 0, colSums(netmatrix_parse) != 0]
  
  # #square up matrix
  # m1<-matrix(ncol=ncol(netmatrix_parse), nrow=ncol(netmatrix_parse), data=0, 
  #            dimnames = list(c(colnames(netmatrix_parse)), c(colnames(netmatrix_parse))))
  # m2<-rbind(netmatrix_parse, m1)
  # 
  # m3<-matrix(ncol=nrow(netmatrix_parse), nrow=nrow(netmatrix_parse), data=0, 
  #            dimnames = list(c(rownames(netmatrix_parse)), c(rownames(netmatrix_parse))))
  # m4<-rbind(m3, t(netmatrix_parse))
  # 
  # netmatrix_parse<-(cbind(m4, m2))
  # 
  # 
  # #make the network-object ####
  # # g <- graph.adjacency(as.matrix(netmatrix_parse),mode="undirected",weighted=NULL)
  
  #---
  g <- graph.incidence(netmatrix_parse)
  
  #---  
  vdf <- get.data.frame(g, 'vertices')
  edf <- get.data.frame(g, 'edges')
  
  # taxonomic color
  lst_pie <- list(hc1 = list(pal=colorRampPalette(c('red','green','blue')),
                             relabu=mr_1_relabu,
                             root=switch(hc[1],
                                         '01_16S_bact' = 'Bacteria',
                                         '02_18S_euk'  = 'Eukaryota',
                                         '05_ITS_fun'  = 'Fungi')),
                  hc2 = list(pal=colorRampPalette(c('cyan','yellow','magenta')),
                             relabu=mr_2_relabu,
                             root=switch(hc[2],
                                         '02_18S_euk' = 'Eukaryota',
                                         '05_ITS_fun' = 'Fungi',
                                         '08_16S_cyano' = 'Cyanobacteria')))
  
  pal <- NULL
  tax_lev <- 2:5
  for(i in seq_along(lst_pie)){
    taxo <- lst_comm[[hc[i]]]$raw$taxo
    taxo <- lst_pie[[i]][['taxo']] <- taxo[row.names(taxo) %in% unlist(dimnames(netmatrix_parse)),]
    
    mr <- lst_pie[[i]]$relabu[,row.names(taxo)]
    
    pie <- pie_taxo(mr, taxo, tax_lev, pal_ini=lst_pie[[i]]$pal,
                    thresh=0, show=F, last_tax_text=F, root=lst_pie[[i]]$root)
    
    nlt <- names(pie$lst_pal[[length(pie$lst_pal)]])
    nlt <- paste0(nlt, ' (', seq_along(nlt), ')')
    
    pal_add <- pie$lst_pal[[length(pie$lst_pal)]]
    names(pal_add) <- nlt
    
    pal <- c(pal, pal_add)
    
    lst_pie[[i]][['pie']] <- pie
  }
  
  tax_tot <- lapply(lst_pie, function(x) x$taxo[,tax_lev])
  names(tax_tot[[1]]) <- names(tax_tot[[2]]) <- paste0('tax', seq_along(tax_lev))
  tax_tot <- rbind(tax_tot[[1]], tax_tot[[2]])
  tax_tot <- tax_tot[vdf$name,]
  
  col <- NULL
  for(i in tax_tot[,ncol(tax_tot)]){
    col <- c(col, pal[i == sapply(strsplit(names(pal), ' '), '[[', 1)])
  }
  
  tax_tot <- cbind.data.frame(tax_tot, nb=sapply(strsplit(names(col), ' '), '[[', 2), col)
  
  V(g)$colour <- V(g)$color <- as.character(tax_tot$col)
  V(g)$label <- V(g)$name
  V(g)$size <- 3
  V(g)$label.cex <- 0.5
  
  #---
  # file <- paste0(dir_cooc, 'graph_bjorn.graphml')
  # write.graph(g, file, 'graphml')
  # 
  # # creation of the graf----
  # file <- paste0(dir_cooc, 'coord_bjorn')
  # system(paste0('extract_pos_gexf.pl ', dir_cooc, 'bjorn.gexf ', file))
  # 
  # coord <- read.table(file)
  # row.names(coord) <- coord$V2
  # coord <- coord[V(g)$name,]
  
  V(g)$label <- paste(row.names(tax_tot),tax_tot$nb)
  
  
  #---
  wdt <- 15
  hei <- 0.2*max(dim(netmatrix_parse)) + 2
  
  # cairo_ps(paste(dir_cooc, 'network_ITS_18S.eps'), width=12, height=8)
  pdf(paste(dir_cooc, 'network_', paste(hc, collapse='_'), '.pdf'), width=wdt, height=hei)
  layout(matrix(c(1,2, 1,3), nrow=2, byrow=T), width=c(1.5,1))
  par(mai=rep(1,4), xpd=NA)
  
  coord <- layout_as_bipartite(g)
  tax_tot <- data.frame(tax_tot, x=coord[,2], y=coord[,1])
  
  #---  
  plot(tax_tot$x, tax_tot$y, axes=F, type='n', xlab='', ylab='')
  
  seg <- sapply(edf, function(x) sapply(x, function(y) tax_tot[y,'y']))
  is_lic <- rownames(tax_tot)[apply(tax_tot[,1:4], 1, function(x) any(grepl('Cyanobacteria|Basidiomycota|Ascomycota|Chlorophyta', x)))]
  col <- ifelse(apply(edf, 1, function(x) all(x %in% is_lic)), 2, 'grey80')
  segments(0,seg[,2],1,seg[,1], col=col)
  
  points(tax_tot$x, tax_tot$y, pch=19, col=as.character(tax_tot$col))
  for(i in 0:1){
    ind <- tax_tot$x == i
    text(tax_tot$x[ind], tax_tot$y[ind], labels=paste(row.names(tax_tot), tax_tot$nb)[ind], pos=c(2,4)[i+1])
  }
  
  #---
  for(i in names(lst_pie)){
    plot.new()
    legend_pie_taxo(lst_pie[[i]]$pie, 0.5,0.5, cex=0.75, last_tax_text=F)
  }
  
  dev.off()
}
#####














