#####
# climarctic cooccu
#####

print('#####
      Climarctic 04 coocurrence networks
      #####')

rm(list=ls())

require(doSNOW)
require(foreach)
require(igraph)
require(compositions)

source('~/bin/src/my_prog/R/pie_taxo.r')
source('~/bin/src/my_prog/R/legend_pie_taxo.r')

# prep cluster
cl <- makeSOCKcluster(4)

clusterEvalQ(cl, library(igraph))
clusterEvalQ(cl, library(compositions))

registerDoSNOW(cl)

# dir loads ####
dir_out  <- 'Projets/Climarctic/stats/MBC/out/'
dir_save <- paste0(dir_out, 'saves/') 
dir_cooc  <- paste0(dir_out, '04_cooc/')
dir.create(dir_cooc, showWarnings=F)

#---
file <- paste0(dir_save, '00_lst_comm.Rdata')
load(file)

# file <- paste0(dir_save, '03_lsts.Rdata')
# load(file)
# 
file <- paste0(dir_save, '01_lst_pie.Rdata')
load(file)

# some function (https://github.com/ggloor/CoDaSeq/blob/master/chunk/phiCode.R)

sma.df <- function(df){
  df.cor <- stats::cor(df, use="pairwise.complete.obs")
  df.var <- stats::cov(df, use="pairwise.complete.obs")
  df.sd <- sqrt(diag(df.var))
  
  r.rf2 <-
    (outer(diag(df.var), diag(df.var), "-")^2 ) /
    (outer(diag(df.var), diag(df.var), "+")^2 - 4 * df.var^2 )
  
  diag(r.rf2) <- 0
  res.dof     <- nrow(df) - 2
  F           <- r.rf2/(1 - r.rf2) * res.dof
  
  list(b=sign(df.cor) * outer(df.sd, df.sd, "/"),
       p=1 - pf(F, 1, res.dof),
       r2=df.cor^2)
}

propr.phisym <- function (X){ 
  Cov    <- stats::var(X)
  tmp    <- 2 * Cov / outer(diag(Cov), diag(Cov), "+")
  return((1-tmp)/(1+tmp))
}

# loop on the communities ####

foreach(i=names(lst_comm), .verbose=T) %dopar% {

  env <- lst_comm[[i]]$raw$env
  mr <- lst_comm[[i]]$raw$mr
  taxo <- lst_comm[[i]]$raw$taxo
  
  # preparation of the correlation matrix----
  # supress low sequence sample
  ind_low <- env$low_seq
  
  mr <- mr[ind_low == F,]
  mr <- as.matrix(mr[rowSums(mr) != 0,colSums(mr) != 0])
  
  taxo <- taxo[colnames(mr),]
  
  env <- env[row.names(mr),]
  
  # loop on OTU or family
  for(j in c('OTUs','families')){
    
    print(j)
    
    fam_col <- switch(i,
                      '01_16S_bact'  = 5,
                      "02_18S_euk"   = 6,
                      "05_ITS_fun"   = 5,
                      "08_16S_cyano" = 5)
    
    root <- switch(i,
                   "01_16S_bact"   = 'Bacteria',
                   "02_18S_euk"    = 'Eukaryota',
                   "05_ITS_fun"    = 'Fungi',
                   "08_16S_cyano"  = 'Cyanobacteriia')
    
    # compositional without the low relabu
    mr_relabu <- as.matrix(decostand(mr, 'total'))
    if(j == 'families'){
      mr_relabu <- aggregate(t(mr_relabu), list(taxo[[fam_col]]), sum)
      row.names(mr_relabu) <- mr_relabu[,1]
      mr_relabu <- t(mr_relabu[,-1])
    }
    
    mr_hoc <- ifelse(mr_relabu >= 0.0005, mr, 0)
    mr_hoc <- mr_hoc[rowSums(mr_hoc) != 0, colSums(mr_hoc) != 0]
    
    mr_ra_hoc <- mr_relabu[row.names(mr_hoc),colnames(mr_hoc)]
    
    clr_hoc <- as.data.frame(clr(mr_hoc))
  
    # generate a symmetrical regression line of best fit
    sym_reg <- sma.df(clr_hoc) # RAM consuming
    sym_reg[['phi']] <- propr.phisym(clr_hoc) # RAM comsuming
    
    lt <- lapply(sym_reg, function(x) ifelse(col(x) < row(x), x, 0))
    lt <- lapply(lt, function(x) {
      y <- x
      dimnames(y) <- dimnames(sym_reg[[1]])
      return(y)
    })
    
    r2 <- ifelse(lt$phi < 0.01 & lt$p < 0.05, lt$r2, 0)
    
    # creation of the network----
    g <- graph.adjacency(r2, mode='undirected', weighted=T)
    print(length(E(g)))
    print(length(V(g)))

    df_edge <- get.data.frame(g, 'edges')
    E(g)$weight <- E(g)$weight*5
    
    # size
    sz <- log(colSums(mr_ra_hoc))
    sz <- sz-min(sz)
    sz <- (sz/max(sz)+0.1)*10000
    V(g)$size <- sz
    V(g)$label <- colnames(mr_ra_hoc)
    
    # pal taxo
    if(j == 'OTUs'){
      tax_hoc <- droplevels(taxo[colnames(mr_hoc),])
    } else {
      tax_hoc <- as.data.frame(t(as.data.frame(strsplit(unique(apply(taxo[grep(paste0(colnames(mr_hoc), collapse='|'),
                                                                               taxo[[fam_col]]),1:fam_col], 1, 
                                                                     function(x) paste(x, collapse=';'))), ';'))))
      row.names(tax_hoc) <- tax_hoc[,ncol(tax_hoc)]
      tax_hoc <- tax_hoc[colnames(mr_hoc),]
    }
    
    pie <- pie_taxo(mr_hoc, tax_hoc, 1:fam_col, show=F, thresh=0.01, root=root)
    pal_tax <- pie$lst_pal[[fam_col-1]]
    
    col_tax <- rep(NA, nrow(tax_hoc))
    new_names <- names(clr_hoc)
    
    for(k in seq_along(col_tax)){
      fam <- as.character(tax_hoc[k,fam_col])
      for(l in fam_col:1){
        
        if(fam %in% names(pal_tax)){
          new_names[k] <- paste0(new_names[k], ' (', which(names(pal_tax) == fam), ')')
          col_tax[k] <- pal_tax[fam]
          break
        
        } else {
          if(l == 1){
            new_names[k] <- paste0(new_names[k], ' (', length(pal_tax), ')')
            col_tax[k] <- pal_tax[length(pal_tax)]
            break
          }
          
          fam <- paste(c(as.character(tax_hoc[k,l-1]), '_', rep('X', fam_col-(l-1))), collapse='')
        
        }     
      }
    }
    
    V(g)$color <- col_tax
    
    # pal factors
    pal_fact <- as.data.frame(sapply(fact_3, function(w) {
      col_fact <- sapply(as.data.frame(sapply(as.data.frame(mr_hoc), function(x) rank(tapply(x, list(env[row.names(mr_hoc),w]), sum)))),
                         function(y) names(lst_palev[[w]])[which(y == length(lst_palev[[w]]))])
      return(sapply(col_fact, function(x) if(length(x)){x} else {NA}))
    }))
    
    V(g)$site <- as.character(pal_fact$site)
    V(g)$moisture <- as.character(pal_fact$moisture)
    V(g)$depth <- as.character(pal_fact$depth)
    
    # write the graf
    file <- paste0(dir_cooc, 'graph_', i, '_', j, '.graphml')
    write.graph(g, file, 'graphml')
    
    # creation of the graf----
    file <- paste0(dir_cooc, i, '_', j, '.coord')
    system(paste0('extract_pos_gexf.pl ', dir_cooc, i, '_', j, '.gexf ', file))

    coord <- read.table(file)
    row.names(coord) <- coord$V2
    coord <- coord[V(g)$label,]

    #---
    pdf(paste0(dir_cooc, i, '_', j, '.pdf'), width=30, height=20)
    par(mar=c(0,0,0,50))

    for(k in c('taxo','site','moisture','depth')){

      # arrange vertex
      V(g)$frame.color <- NA

      if(k == 'taxo'){
        V(g)$shape <- 'circle'
        V(g)$label <- new_names
        V(g)$label.cex <- 0.3

        if(nrow(df_edge)){
          E(g)$color <- apply(apply(df_edge[,1:2], 1,
                                    function(x) {
                                      apply(col2rgb(col_tax[colnames(mr_hoc) %in% x]), 1, function(y) mean(y)/255)
                                    }), 2, function(z) rgb(z[1],z[2],z[3]))
        }

      } else {

        agg <- aggregate(mr_ra_hoc, list(env[row.names(mr_ra_hoc),k]), mean)
        row.names(agg) <- agg$Group.1
        agg <- agg[,-1]

        V(g)$shape <- 'pie'
        V(g)$pie <- as.list(agg)
        V(g)$pie.color <- list(lst_palev[[k]])
        V(g)$pie.lty <- 0
        V(g)$frame.color <- NA
        V(g)$label <- NA

        t(col2rgb(lst_palev[[k]])/255)

        if(nrow(df_edge)){
          E(g)$color <- apply(df_edge[,1:2], 1,
                              function(x) {
                                r <- rowMeans(agg[,colnames(agg) %in% x])
                                r <- r/sum(r)
                                rgb <- colSums(t(col2rgb(lst_palev[[k]])/255)*r)
                                return(rgb(rgb[1],rgb[2],rgb[3]))
                              })
        }

      }

      # plot
      plot(g, layout=as.matrix(coord[,3:4]), rescale=F, xlim=range(coord[,3]), ylim=range(coord[,4]))
      usr <- par('usr')

      # legend
      x_leg <- usr[2]+diff(usr[1:2])*0.20
      y_leg <- mean(usr[3:4])
      if(k == 'taxo'){
        legend_pie_taxo(pie, x_leg, y_leg, diff(usr[1:2])*0.4, last_tax_text=F)
      } else {
        legend(x_leg, y_leg, names(lst_palev[[k]]), col=lst_palev[[k]], bty='n', pch=19)
      }

    }

    dev.off()
  }
}








#####





