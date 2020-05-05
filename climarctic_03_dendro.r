#####
# climarctic dendro
#####

print('##### Climarctic 03 dendrogram #####')

rm(list=ls())

require(dendextend)


# dir loads ####
dir_out  <- 'Projets/Climarctic/stats/MBC/out/'
dir_save <- paste0(dir_out, 'saves/') 
dir_dend  <- paste0(dir_out, '03_dendro/')
dir.create(dir_dend, showWarnings=F)

#---
file <- paste0(dir_save, '00_lst_comm.Rdata')
load(file)

file <- paste0(dir_save, '02_lst_rda.Rdata')
load(file)

file <- paste0(dir_save, '01_lst_pie.Rdata')
load(file)

# loop the taxa
for(i in n_comm) {
  
  mr_clr <- lst_comm[[i]]$clr2$mr
  mr_raw <- lst_comm[[i]]$raw$mr[row.names(mr_clr),names(mr_clr)]
  taxo <- lst_comm[[i]]$raw$taxo[names(mr_raw),]
  env <- lst_comm[[i]]$clr2$env
  
  # dendrogram ----
  
  # Aitchison dist ---
  di_ait <- dist(mr_clr)
  
  # clustering ---
  hc <- hclust(di_ait)
  dend <- as.dendrogram(hc)

  # arrange according to the factors ---
  fact1 <- 'depth'
  fact2 <- 'moisture'
  fact3 <- 'site'
  
  dend <- set(dend, 'branches_lwd', 2)
  
  # dend <- set(dend, 'by_labels_branches_col', row.names(env)[env[[fact1]] == 'top'],
  #             TF_value=c(lst_palev[[fact1]][1], lst_palev[[fact1]][2]))
  
  for(j in names(lst_palev[[fact1]])){
    for(k in row.names(env)[env[[fact1]] == j]){
      dend <- set(dend, 'by_labels_branches_col', k, TF_value=c(lst_palev[[fact1]][j], Inf))
    }
  }
  
  dend <- set(dend, 'labels_colors', lst_palev[[fact3]][env[labels(dend),fact3]])
  
  # rotate according to rh  
  dend <- rotate(dend, row.names(env)[order(env$rh)])
  
  # heatmap variables ---
  pal_var <- colorRampPalette(c('red','green'))(101)
  
  var <- rev(attributes(lst_rda$`top|deep`[[i]]$lst_rda$parsi$mod$terminfo$terms)$term.labels)
  
  e <- env[,var]
  e <- e[labels(dend),sapply(e, is.numeric)]

  # scale the variables
  p <- sapply(e, function(x){
    x <- x-min(x, na.rm=T)
    return(pal_var[round((x/max(x, na.rm=T))*100)+1])
  })
  row.names(p) <- row.names(e)
  
  # barplots ---
  # comm
  agg <- lst_pie[[i]]$abundance$per_smp$agg
  fact <- sapply(agg, is.factor)
  
  comm <- agg[,fact == F]
  names(comm) <- sapply(strsplit(names(comm), ' '), '[[', 1)
  comm <- aggregate(comm[,labels(dend)], list(agg[,2]), sum)
  
  row.names(comm) <- comm[,1]
  comm <- comm[,-1]
  root <- switch(i,
                 "01_16S_bact"   = 'Bacteria_X',
                 "02_18S_euk"    = 'Eukaryota_X',
                 "05_ITS_fun"    = 'Fungi_X',
                 "08_16S_cyano"  = 'Cyanobacteria_X')
  ind_root <- which(row.names(comm) == root)
  comm <- rbind.data.frame(comm[-ind_root,],comm[ind_root,])
  
  # pal
  pal_tax <- lst_pie[[i]]$abundance$per_smp$lst_pal[[1]]
  
  # plot ---  
  cairo_ps(paste0(dir_dend, 'dend_', i, '.eps'), width=15, height=15)
  par(mar=c(ncol(e)+30, 4, 4, 7), xpd=NA)
  
  plot(dend)
  # plot(hang.dendrogram(dend))
  usr <- par('usr')
  xrng <- diff(usr[1:2])
  yrng <- diff(usr[3:4])
  
  points(1:nrow(e), rep(0, nrow(e)), pch=19, col=lst_palev$moisture[env[labels(dend),fact2]])
  
  # variable
  x0 <- 1:ncol(comm)
  y0 <- usr[3]-yrng*0.2
  y_stretch <- yrng*0.1
  
  ys <- seq(y0, y0+y_stretch, length.out=ncol(p)+1)
  xleg <- c(usr[2]+xrng*c(0.01,0.04))
  
  for(j in 1:ncol(p)){
    y <- mean(ys[c(j,j+1)])
    rect(x0-0.5, ys[j], x0+0.5, ys[j+1],col=p[,j], border=NA)
    text(c(usr[1]+xrng*0.00,xleg), y, c(colnames(p)[j], signif(range(e[,j], na.rm=T), 2)),
         pos=c(2,2,4), cex=0.75)
    points(seq(xleg[1], xleg[2], length.out=length(pal_var)), rep(y, length(pal_var)),
           pch=19, col=pal_var)
  }
  
  # taxa
  y0 <- usr[3]-yrng*.75
  y_stretch <- yrng*0.4
  
  y <- mean(c(y0,y0+y_stretch))
  
  for(j in x0){
    cs <- cumsum(rev(comm[,j]))*y_stretch
    rect(j-0.5, c(0,cs[-length(cs)])+y0, j+0.5, cs+y0, col=rev(pal_tax), border=NA, xpd=NA)
  }
  
  text(usr[1]+xrng*0.00, y, 'taxa relative abundance', srt=90)

  legend(usr[2]+xrng*0, y, row.names(comm), pch=19,
         col=pal_tax, xpd=NA, bty='n', cex=0.75, xjust=0, yjust=0.5)
    
  
  #---
  dev.off()
  
}







#####























