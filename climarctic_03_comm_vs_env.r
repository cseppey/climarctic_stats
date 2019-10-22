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
# file <- paste0(dir_save, '01_lst_comm.Rdata')
# load(file)
file <- paste0(dir_save, '03_lst_comm.Rdata')
if(file.exists(file)){
  load(file)
}

#########
permu <- 1000
#########


fact_3 <- c('site','moisture','depth')

# file for the model summaries
file <- paste0(dir_cve, 'test_models.csv')
if(file.exists(file)){
  file.remove(file)
}

# file for centroids summaries
file <- paste0(dir_cve, 'centro.csv')
if(file.exists(file)){
  file.remove(file)
}

# loop on communities
for(i in names(lst_comm)){
  for(j in c('raw','rrf','rrf2')){

    if(is.null(lst_comm[[i]][[j]])){
      next
    }
    
    env <- lst_comm[[i]][[j]]$env
    e <- na.omit(env[,c('site','moisture','depth','plot_in_moist_in_site','pH','C_N','sand','silt','clay','rh','om')])
    ee <- e
    
    is.fact <- sapply(ee, is.factor)
    
    
    # transfo ####
    print(paste(i, j, 'transfo'))
    
    mr <- lst_comm[[i]][[j]]$mr[row.names(e),]
    mr <- mr[rowSums(mr) != 0,colSums(mr) != 0]
    e <- e[row.names(mr),]
    ee <- e
    
    fact_plot <- apply(sapply(ee[,4:3], as.character), 1, function(x) paste(x, collapse='_'))
    fact_plot <- factor(fact_plot, levels=unique(fact_plot))
    lfp <- levels(fact_plot)
    
    ee <- ee[,-4]
    
    transfo <- c('raw','hell','log')
    
    if('lst_trsf_bc' %in% names(lst_comm) == F){
      lst_trsf <- foreach(k=transfo) %dopar% {
        if(k == 'raw'){
          return(list(mr=mr, bc=vegdist(mr)))
        } else {
          m <- decostand(mr, k)
          return(list(mr=m, bc=vegdist(m)))
        }
      }
      names(lst_trsf) <- transfo
      
      lst_comm[[i]][[j]][['lst_trsf_bc']] <- lst_trsf
    }
    
    
    # model ####
    print(paste(i, j, 'model'))
    
    res_tst <- foreach(k=transfo, .verbose=T) %dopar% {
      bc <- lst_trsf[[k]]$bc
      
      lst_mod <- NULL
      
      var_coor_pv <- array(NA, dim=c(length(ee)+1, 3, 2),
                           dimnames=list(c(names(ee), 'stress_aic'), c('x','y','pv'), c('NMDS','RDA')))
      dn <- dimnames(var_coor_pv)
      #
      
      # NMDS ----
      nmds <- metaMDS(bc)
      
      # test + coord
      ef <- envfit(nmds, ee, na.rm=T, permutations=permu)
      
      for(l in dn[[1]]){
        if(l != ('stress_aic')){
          if(is.fact[names(is.fact) == l]) {
            var_coor_pv[l,3,1] <- ef$factors$pvals[l]
          } else {
            var_coor_pv[l,,1] <- c(ef$vectors$arrows[l,], ef$vectors$pvals[l])
          }
        } else {
          var_coor_pv[l,3,1] <- nmds$stress
        }
      }
      
      # coord
      site <- nmds$points
      
      # centroids per plot
      smp_cls <- NULL
      
      for(l in lfp){
        
        ind_plt <- which(fact_plot == l)
        
        lp <- length(ind_plt)
        rn <- row.names(site)[ind_plt]
        
        if(lp == 1){
          smp_cls <- c(smp_cls, rn)
        
        } else if (lp == 2){
          smp_cls <- c(smp_cls, paste(rn, collapse='|'))
          
        } else {
          site_plt <- site[rn,]
          cs <- colSums(apply(site_plt, 1, function(x) abs(x-apply(site_plt, 2, mean))))
          smp_cls <- c(smp_cls, names(which(cs == min(cs))))
        }
      }
      names(smp_cls) <- lfp
      
      axis_n <- c('NMDS 1','NMDS 2')
      
      lst_mod[['NMDS']] <- list(site=site, axis_n=axis_n, smp_cls=smp_cls)
      
      
      # RDA ----
      rda <- capscale(bc~., data=ee) 
      
      s <- summary(rda)
      
      variance <- signif(s$cont$importance[2,1:2], 2)
      
      # test
      pvs <- NULL
      for(l in names(ee)){
        set.seed(0)
        pvs <- c(pvs, signif(anova(rda, update(rda, as.formula(paste('.~.-', l))), parallel=4, permutations=permu)$`Pr(>F)`[2], 2))
      }
      pvs[is.na(pvs)] <- 1
      names(pvs) <- names(ee)
      
      # coord
      
      for(l in dn[[1]]){
        if(l != ('stress_aic')){
          if(is.fact[names(is.fact) == l]) {
            var_coor_pv[l,3,2] <- pvs[l]
          } else {
            if(l %in% row.names(s$biplot)){
              var_coor_pv[l,,2] <- c(s$biplot[l,1:2], pvs[l])
            }
          }
        } else {
          var_coor_pv[l,3,2] <- extractAIC(rda)[2]
        }
      }
      
      # coord
      site <- s$sites[,1:2]

      # centroids per plot
      smp_cls <- NULL
      
      for(l in lfp){
        
        ind_plt <- which(fact_plot == l)
        
        lp <- length(ind_plt)
        rn <- row.names(site)[ind_plt]
        
        if(lp == 1){
          smp_cls <- c(smp_cls, rn)
          
        } else if (lp == 2){
          smp_cls <- c(smp_cls, paste(rn, collapse='|'))
          
        } else {
          site_plt <- s$site[rn,]
          cs <- colSums(apply(site_plt, 1, function(x) abs(x-apply(site_plt, 2, mean))))
          smp_cls <- c(smp_cls, names(which(cs == min(cs))))
        }
      }
      names(smp_cls) <- lfp
      
      axis_n <- paste('RDA', 1:2, variance) 
      
      lst_mod[['RDA']] <- list(site=site, axis_n=axis_n, smp_cls=smp_cls)
      
      
      # graf ----
      for(ln in names(lst_mod)){
        l <- lst_mod[[ln]]
        v <- as.data.frame(var_coor_pv[,,ln])
        
        coord <- rbind(l$site, as.matrix(v[,1:2]))
        
        #---
        cairo_ps(paste0(dir_cve, ln, '_', i, '_', j, '_', k, '.eps'), width=10, height=7)
        layout(matrix(c(1,3,4, 1,5,6, 0,2,2), nrow=3, byrow=T), width=c(0.2,1,1), height=c(1,1,0.3))
        par(mar=rep(0,4))
        
        # axis
        for(m in 2:1){
          plot.new()
          text(0.5,0.5,labels=l$axis_n[m], srt=ifelse(m == 1, 0, 90))
        }
        
        # ordinations
        for(m in fact_3){
          plot(NA, xlim=range(coord[,1], na.rm=T), ylim=range(coord[,2], na.rm=T))
          
          abline(v=0,h=0,lty=3)
          
          arrows(0,0,v[,1],v[,2], length=0, lty=2, col='grey70')
          text(v[,1:2], labels=paste(row.names(v), signif(v$pv, 3)))
          
          ordispider(l$site, fact_plot, col='grey')
          
          text(l$site, labels=row.names(l$site), col=lst_palev[[m]]$pal[e[[m]]])     
          
        }
        
        # legend ---
        plot.new()
        
        # factors colors
        legend(0.33,0.5, legend=c(unlist(sapply(lst_palev, '[[', 2))), bty='n', xjust=0.5, yjust=0.5,
               pch=19, col=unlist(sapply(lst_palev, '[[', 1)))
        
        # factors significance
        is_fact_sa <- c(is.fact, TRUE)
        df_fct <- v[is_fact_sa,]
        
        legend(0.66,0.5, legend=paste(row.names(df_fct), signif(df_fct$pv,3)), bty='n', xjust=0.5, yjust=0.5)
        
        #---
        dev.off()
      }
      
      #---
      return(list(var=var_coor_pv[,3,], smp_cls=sapply(lst_mod, '[[', 3)))
    }
    
    names(res_tst) <- transfo
    
    # output model
    df <- cbind(res_tst$raw$var, res_tst$hell$var, res_tst$log$var) 
    colnames(df) <- apply(expand.grid(colnames(res_tst[[1]][[1]]), names(res_tst)), 1, function(x) paste(x, collapse='_'))  
    
    file <- paste0(dir_cve, 'test_models.csv')
    
    write.table(paste(i, j), file=file, append=T, row.names=F, col.names=F, quote=F)
    write.table(df, file=file, append=T, quote=F)      
    
    # output centroids
    df <- cbind(res_tst$raw$smp_cls, res_tst$hell$smp_cls, res_tst$log$smp_cls) 
    colnames(df) <- apply(expand.grid(colnames(res_tst[[1]][[1]]), names(res_tst)), 1, function(x) paste(x, collapse='_'))  
    
    file <- paste0(dir_cve, 'centro.csv')
    
    write.table(paste(i, j), file=file, append=T, row.names=F, col.names=F, quote=F)
    write.table(df, file=file, append=T, quote=F)      
    
  }
}

file <- paste0(dir_save, '03_lst_comm.Rdata')
save(lst_comm, file=file)

#




















