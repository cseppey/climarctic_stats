#####
# climarctic comm vs env
#####

rm(list=ls())
gc()

require(foreach)
require(doSNOW) # e.g. makeSOCKcluster()

# prep cluster
cl <- makeSOCKcluster(4)

registerDoSNOW(cl)


# dir loads ####
dir_out  <- 'Projets/Climarctic/stats/MBC/out/'
dir_save <- paste0(dir_out, 'saves/') 
dir_cve  <- paste0(dir_out, '03_comm_vs_env/')
dir.create(dir_cve, showWarnings=F)

#---
file <- paste0(dir_save, '00_lst_comm.Rdata')
load(file)
file <- paste0(dir_save, '01_lst_comm.Rdata')
load(file)

fact_3 <- c('site','moisture','depth')

# loop on communities
file <- paste0(dir_cve, 'test_models.csv')
if(file.exists(file)){
  file.remove(file)
}

for(i in names(lst_comm)){
  for(j in c('raw','rrf')){

    if(is.null(lst_comm[[i]][[j]])){
      next
    }
    
    env <- lst_comm[[i]][[j]]$env
    e <- na.omit(env[,c('site','moisture','depth','pH','C_N','sand','silt','clay','rh','om')])
    ee <- e
    
    is.fact <- sapply(ee, is.factor)
    
    # transfo ####
    print(paste(i, j, 'transfo'))
    
    mr <- lst_comm[[i]][[j]]$mr[row.names(e),]
    mr <- mr[rowSums(mr) != 0,colSums(mr) != 0]
    e <- e[row.names(mr),]
    ee <- e
    
    transfo <- c('raw','hell','log')
    # lst_trsf <- foreach(k=transfo) %dopar% {
    #   if(k == 'raw'){
    #     return(list(mr=mr, bc=vegdist(mr)))
    #   } else {
    #     m <- decostand(mr, k)
    #     return(list(mr=m, bc=vegdist(m)))
    #   }
    # }
    # names(lst_trsf) <- transfo
    
    #---
    file <- paste0(dir_save, '03_transfo_', i, '_', j, '.Rdata')
    # save(lst_trsf, file=file)
    load(file)
    
    
    # model ####
    print(paste(i, j, 'model'))
    
    res_tst <- foreach(k=transfo) %dopar% {
      bc <- lst_trsf[[k]]$bc
      
      lst_mod <- NULL
      #
      
      # NMDS ----
      nmds <- metaMDS(bc)
      
      # test
      ef <- envfit(nmds, ee, na.rm=T)
      
      pvv <- ef$vectors
      pvf <- ef$factors
      
      # coord
      site <- nmds$points
      
      var <- pvv$arrows 
      row.names(var) <- paste(row.names(var), pvv$pvals)
      
      fact <- paste(names(pvf$pvals), pvf$pvals)
      
      axis_n <- c('NMDS 1','NMDS 2')
      
      lst_mod[['NMDS']] <- list(site=site, var=var, fact=fact, axis_n=axis_n)
      
      
      # RDA ----
      rda <- capscale(bc~., data=ee) 
      
      s <- summary(rda)
      
      variance <- signif(s$cont$importance[2,1:2], 2)
      
      # test
      pvs <- NULL
      for(l in names(ee)){
        set.seed(0)
        pvs <- c(pvs, signif(anova(rda, update(rda, as.formula(paste('.~.-', l))), parallel=4, permutations=1000)$`Pr(>F)`[2], 2))
      }
      pvs[is.na(pvs)] <- 1
      names(pvs) <- names(ee)
      
      # coord
      site <- s$sites[,1:2]
      
      var <- s$biplot[grep(paste(names(ee)[is.fact == F], collapse='|'), row.names(s$biplot)),1:2]
      row.names(var) <- paste(row.names(var), pvs[row.names(var)])

      fact <- paste(c(names(pvs[is.fact]), 'aic'), c(pvs[is.fact], signif(extractAIC(rda)[2], 3)))
      
      axis_n <- paste('RDA', 1:2, variance) 
      
      lst_mod[['RDA']] <- list(site=site, var=var, fact=fact, axis_n=axis_n)
      
      
      # graf ----
      for(ln in names(lst_mod)){
        l <- lst_mod[[ln]]
        
        coord <- rbind(l$site, l$var)
        
        #---
        cairo_ps(paste0(dir_cve, ln, '_', i, '_', j, '_', k, '.eps'), width=10, height=7)
        layout(matrix(c(1,3,4, 1,5,6, 0,2,2), nrow=3, byrow=T), width=c(0.2,1,1), height=c(1,1,0.3))
        par(mar=rep(0,4))
        
        for(m in 2:1){
          plot.new()
          text(0.5,0.5,labels=l$axis_n[m], srt=ifelse(m == 1, 0, 90))
        }
        
        for(m in fact_3){
          plot(NA, xlim=range(coord[,1]), ylim=range(coord[,2]))
          
          abline(v=0,h=0,lty=3)
          
          arrows(0,0,l$var[,1],l$var[,2], length=0, lty=2, col='grey70')
          text(l$var, labels=row.names(l$var))
          
          text(l$site, labels=row.names(l$site), col=lst_palev[[m]]$pal[e[[m]]])        
          
        }
        
        plot.new()
        legend(0.33,0.5,legend=c(unlist(sapply(lst_palev, '[[', 2))), bty='n', xjust=0.5, yjust=0.5,
               pch=19, col=unlist(sapply(lst_palev, '[[', 1)))
        if(length(l$fact)){legend(0.66,0.5,legend=l$fact, bty='n', xjust=0.5, yjust=0.5)}
        
        dev.off()
      }
      
      #---
      return(c(i,j,k,'NMDS', row.names(lst_mod$NMDS$var),lst_mod$NMDS$fact, 'RDA', row.names(lst_mod$RDA$var),lst_mod$RDA$fact))
    }
    
    file <- paste0(dir_cve, 'test_models.csv')
    for(k in res_tst){
      write.table(k, file=file, append=T, row.names=F, col.names=F, quote=F)      
    }
    
    
  }
}


#




















