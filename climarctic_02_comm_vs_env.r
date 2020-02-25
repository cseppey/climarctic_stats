#####
# climarctic comm vs env
#####

rm(list=ls())

require(foreach)
require(doSNOW) # e.g. makeSOCKcluster()

# prep cluster
cl <- makeSOCKcluster(2)

clusterEvalQ(cl, library(foreach))
clusterEvalQ(cl, library(doSNOW))

registerDoSNOW(cl)


# dir loads ####
dir_out  <- 'Projets/Climarctic/stats/MBC/out/'
dir_save <- paste0(dir_out, 'saves/') 
dir_cve  <- paste0(dir_out, '03_comm_vs_env/')
dir.create(dir_cve, showWarnings=F)

#---
file <- paste0(dir_save, '00_lst_comm.Rdata')
load(file)


#########
permu <- 1000
#########


fact_3 <- c('site','moisture','depth')

# loop on communities
lsts <- foreach(h=c('top|deep','top','deep'), .verbose=T) %dopar% {
  
  cl2 <- makeSOCKcluster(2)
  registerDoSNOW(cl2)
  
  lst <- NULL
  for(i in names(lst_comm)) {
    
    mr <- lst_comm[[i]]$clr2$mr
    env <- lst_comm[[i]]$clr2$env
    
    env <- env[row.names(mr),]
    env <- env[,c('site','moisture','depth','plot_in_moist_in_site',
                  'pH','C','N','C_N','S','silt','clay','rh','om')]
    
    # ordination ####
    print(paste(i, 'model'))
    
    cond <- 'site'
    # Note: if many factor included in the condition (interaction or not), 
    # the only factor counting is the one having the highest resolution if hierarchical factors
    
    e <- na.omit(env[grepl(h, env$depth),])
    is.fact <- sapply(e, is.factor)
    m <- mr[row.names(e),]
    m <- m[,colSums(m) != 0]
    
    # RDA ----
    f <- formula(paste0(ifelse(h == 'top|deep', 'm~depth+', 'm~'), paste(names(e)[is.fact == F], collapse='+'),
                        '+Condition(', cond, ')', collapse=''))
    rda_full <- capscale(f, data=e, scale=T)
    
    lst_rda <- list(full=list(e=e, m=m, f=f, mod=rda_full))
    
    # retreive a more parsimonious model ---
    
    # Note: the parsimonious model did not depend of the order of the variable in the input model
    
    print('ordistep')
    
    file <- paste0(dir_save, '03_rda_parsi_clr2_', i, '_', h, '.Rdata')
    if(file.exists(file) == F){
      set.seed(0)
      f <- formula(paste0(ifelse(h == 'top|deep', 'm~depth+rh+Condition(', 'm~rh+Condition('), cond, ')'))
      rda_parsi <- ordistep(capscale(f, data=e, scale=T), formula(rda_full), parallel=4, trace=F)
      save(rda_parsi, file=file)
    } else {
      load(file)
    }
    
    v_parsi <- attributes(rda_parsi$terminfo$terms)$term.labels
    
    e <- na.omit(env[grepl(h, env$depth),v_parsi])
    
    m <- mr[row.names(e),]
    m <- m[,colSums(m) != 0]
    
    f <- formula(paste(c(ifelse(h == 'top|deep', 'm~depth+rh', 'm~rh'),
                         v_parsi[grepl('depth|rh|site', v_parsi) == F], 
                         paste('Condition(', cond, ')')), collapse='+'))
    
    rda_parsi <- capscale(f, data=e, scale=T)
  
    lst_rda[['parsi']] <- list(e=e, m=m, f=f, mod=rda_parsi)
    
    # retreive info on the two models ---
    
    file <- paste0(dir_save, '03_lst_ordi_clr2_', i, '_', h, '.Rdata')
    if(file.exists(file) == F){
    
      lst_ordi <- NULL
      for(j in 1:2){
        
        print(names(lst_rda)[j])
        
        e   <- lst_rda[[j]]$e
        f   <- lst_rda[[j]]$f
        mod <- lst_rda[[j]]$mod
        
        v <- rev(rev(attributes(mod$terminfo$terms)$term.labels)[-1])
        
        if(length(v) > 1){
          
          stat <- as.data.frame(foreach(k=v, .combine=rbind) %dopar% {
            
            m <- lst_rda[[j]]$m
            mod <- capscale(eval(parse(text=paste(as.character(f)[c(2,1,3)], collapse=''))), 
                            data=e, scale=T) # some fuck with the foreach environment 
                                             # but output are the same then if tested sequencially
            
            set.seed(0)
            
            upd <- update(mod, as.formula(paste('.~.-', k)))
            pv1 <- signif(anova(mod, upd, parallel=4, permutations=permu)$`Pr(>F)`[2], 2)
            va1 <- RsquareAdj(mod)[[1]] - RsquareAdj(upd)[[1]] # no R2 adj because Conditon(site)
            
            rda1 <- capscale(as.formula(paste('m~', k, '+Condition(', cond, ')')), data=e, scale=T)
            pv2 <- signif(anova(rda1, parallel=4, permutations=permu)$`Pr(>F)`[1], 2)
            va2 <- RsquareAdj(rda1)[[1]]
            
            return(c(pv1, va1, pv2, va2))
          })
          
        } else {
          pv <- signif(anova(rda_parsi, parallel=4, permutations=permu)$`Pr(>F)`[1], 2)
          va <- RsquareAdj(rda_parsi)[[1]]
          stat <- matrix(c(pv, va, pv, va), nrow=1)
        }
        
        dimnames(stat) <- list(v, c('pv1','va1','pv2','va2'))
        
        #---
        s <- summary(mod)
        
        n_fact <- row.names(s$biplot)[row.names(s$biplot) %in% names(is.fact)[is.fact == F]]
        var <- s$biplot[n_fact,1:2]
        if(is.matrix(var) == F){
          var <- matrix(var, ncol=2)
          dimnames(var) <- list(n_fact, c('CAP1','CAP2'))
        }
        
        lst_ordi[[j]] <- list(site       = s$sites[,1:2],
                              var        = var,
                              axes_names = paste('RDA; var =', signif(s$cont$importance[2,1:2], 2)),
                              stat       = stat)
      }
      
      names(lst_ordi) <- names(lst_rda)
      
      
      # NMDS ----
      e <- env[grep(h, env$depth),]
      m <- mr[row.names(e),]
      
      nmds <- metaMDS(m, dist='euc') # aitchinson distance = euclydian dist on CLR
      
      # test variables of parsi
      set.seed(0)
      
      envfit <- envfit(nmds, e[,v_parsi], na.rm=T)
      
      var_env <- lapply(envfit[1:2], function(x) matrix(unlist(x[c('pvals','r')]), ncol=2, 
                                                        dimnames=list(names(x$r), c('pv2','va2'))))
      var_env <- rbind(var_env[[2]],var_env[[1]])
      
      lst_ordi[['nmds']] <- list(site       = nmds$points,
                                 var        = envfit$vectors$arrows,
                                 axes_names = c('NMDS1','NMDS2'),
                                 stat       = var_env)
      
      save(lst_ordi, file=file)
      
    } else {
      load(file)
    }
    
    # graf ----
    print('graf')
    
    pdf(paste0(dir_cve, 'ordination_clr2_', i, '_', h, '.pdf'), width=10, height=7)
    
    for(jn in names(lst_ordi)){
      
      j <- lst_ordi[[jn]]

      # cairo_ps(paste0(dir_cve, 'ordination_', jn, '_', i, '_', h, '.eps'), width=10, height=7)
      layout(matrix(c(1,3,4, 1,5,6, 0,2,2), nrow=3, byrow=T), width=c(0.2,1,1), height=c(1,1,0.3))
      
      # axes
      par(mar=rep(0,4), oma=c(1,1,3,3))
      for(k in 2:1){
        plot.new()
        text(0.5,0.5, labels=j$axes_names[k], srt=ifelse(k == 1, 0, 90))
      }
      
      # ordinations
      site <- j$site
      var <- j$var
      
      coord <- rbind(site, var)
      
      for(k in seq_along(fact_3)){
        
        fact <- fact_3[k]
        
        env_in_ordi <- env[row.names(env) %in% row.names(j$site),]
        
        pal <- lst_palev[[fact]]$pal[env_in_ordi[[fact]]]
        
        #---
        plot(j$site, xlim=range(rbind(coord[,1])), ylim=range(coord[,2]), xaxt='n', yaxt='n',
             col=pal, pch=NA)
        
        if(k > 1){
          axis(1)
        }
        
        if(k %% 2 == 1){
          axis(2)
        }
        
        ordispider(site, env_in_ordi$plot_in_moist_in_site, col='grey80')
        
        ls <- lst_comm[[i]]$clr2$env$low_seq[row.names(lst_comm[[i]]$clr2$env) %in% row.names(site)]
        points(site[ls,], pch=19, col=2)
        
        points(site[,1], site[,2], pch=19, col=pal)
        # text(site, labels=row.names(site), col=pal)
        
        abline(v=0, h=0, lty=3)
        
        #---
        
        mult <- min(abs(sapply(as.data.frame(site), range)))
        v2 <- var*mult
        
        arrows(0,0,v2[,1],v2[,2], length=0, lty=2)

        if(k < 3){        
          axis(3, at=seq(-mult, mult, length.out=5), labels=seq(-1, 1, length.out=5), col=2)
        }
        
        if(k > 1){
          axis(4, at=seq(-mult, mult, length.out=5), labels=seq(-1, 1, length.out=5), col=2)
        }
        
        text(v2[,1:2], labels=paste(row.names(v2)))
        
      }
      
      # legend ---
      plot.new()
      
      # factors colors
      legend(0.5,0.5, legend=c(unlist(sapply(lst_palev, '[[', 2))), bty='n', xjust=0.5, yjust=0.5,
             pch=19, col=unlist(sapply(lst_palev, '[[', 1)))
      
      # dev.off()
    }
    
    dev.off()
  
    ###
    lst[[i]] <- list(lst_rda=lst_rda, lst_ordi=lst_ordi)
  }
  
  stopCluster(cl2)
  
  return(lst)
}

names(lsts) <- c('top|deep','top','deep')

file <- paste0(dir_save, '03_lsts.Rdata')
save(lsts, file=file)
load(file)

#---
file <- paste0(dir_cve, 'stat_ordi.csv')
if(file.exists(file)){
  file.remove(file)
}

for(i in names(lsts)){
  for(j in names(lsts[[i]])){
    for(k in names(lsts[[i]][[j]]$lst_ordi)){
      write.table(c(i,j,k), file, T, F, '\t')
      write.table(lsts[[i]][[j]]$lst_ordi[[k]]$stat, file, T, F, '\t')
    }
  }
}

# procrust rotation on the NMDS with common samples ####

# find the samples found in all datasets
smp_4 <- names(which(table(unlist(lapply(lst_comm, function(x) row.names(x$clr$mr)))) == 4))

# NMDS ---
lst_nmds_4 <- lapply(lst_comm, function(x) {
  m <- x$clr$mr[smp_4,]
  m <- m[,colSums(m) != 0]
  return(metaMDS(m, distance='euc'))
})

# grafs
for(i in names(lst_nmds_4)){
  
  # pdf(paste0(dir_cve, 'nmds_4_', i, '.pdf'), width=10, height=7)
  cairo_ps(paste0(dir_cve, 'nmds_4_', i, '.eps'), width=10, height=7)
  layout(matrix(c(1,3,4, 1,5,6, 0,2,2), nrow=3, byrow=T), width=c(0.2,1,1), height=c(1,1,0.3))
  
  # axes
  par(mar=rep(0,4), oma=c(1,1,3,3))
  for(j in 2:1){
    plot.new()
    text(0.5,0.5, labels=paste0('NMDS', j), srt=ifelse(j == 1, 0, 90))
  }
  
  # ordinations
  site <- lst_nmds_4[[i]]$points
  
  for(j in seq_along(fact_3)){
    
    fact <- fact_3[j]
    
    env_in_ordi <- env_tot[row.names(env_tot) %in% row.names(site),]
    
    pal <- lst_palev[[fact]]$pal[env_in_ordi[[fact]]]
    
    #---
    plot(site, xlim=range(site[,1]), ylim=range(site[,2]), xaxt='n', yaxt='n',
         col=pal, pch=NA)
    
    if(j > 1){
      axis(1)
    }
    
    if(j %% 2 == 1){
      axis(2)
    }
    
    ordispider(site, env_in_ordi$plot_in_moist_in_site, col='grey80')
    
    ls <- lst_comm[[i]]$clr$env$low_seq[row.names(lst_comm[[i]]$clr$env) %in% row.names(site)]
    points(site[ls,], pch=19, col=2)
    
    points(site[,1], site[,2], pch=19, col=pal)
    # text(site, labels=row.names(site), col=pal)
    
    abline(v=0, h=0, lty=3)
    
  }
  
  # legend ---
  plot.new()
  
  text(0.5,0.75,i)
  
  # factors colors
  legend(0.5,0.5, legend=c(unlist(sapply(lst_palev, '[[', 2))), bty='n', xjust=0.5, yjust=0.5,
         pch=19, col=unlist(sapply(lst_palev, '[[', 1)))
  
  dev.off()
  
}

# procrustes ---
lgt <- length(lst_nmds_4)

pvs <- matrix(NA, lgt, lgt, dimnames=list(names(lst_nmds_4),names(lst_nmds_4)))

cairo_ps(paste0(dir_cve, 'procrust.eps'), 11, 8)
par(mfrow=c(2,3))

set.seed(0)
for(i in 1:(lgt-1)){
  for(j in (i+1):lgt){
    pt <- protest(lst_nmds_4[[i]]$points,lst_nmds_4[[j]]$points)
    plot(procrustes(lst_nmds_4[[i]]$points,lst_nmds_4[[j]]$points),
         main=paste(c(names(lst_comm)[c(i,j)], '\np =', pt$signif, 'SS =', signif(pt$ss, 2)), collapse=' '))
    pvs[i,j] <- pt$signif
    pvs[j,i] <- pt$ss
  }
}

dev.off()






#####




















