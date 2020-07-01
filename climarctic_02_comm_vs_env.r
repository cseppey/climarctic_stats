#####
# climarctic comm vs env
#####

print('##### Climarctic 02 communities vs environment #####')

rm(list=ls())

require(foreach)
require(doSNOW) # e.g. makeSOCKcluster()
require(xtable)

# prep cluster
cl <- makeSOCKcluster(4)

registerDoSNOW(cl)


# dir loads ####
dir_out  <- 'Projets/Climarctic/stats/MBC/out/'
dir_save <- paste0(dir_out, 'saves/') 
dir_cve  <- paste0(dir_out, '02_comm_vs_env/')
dir.create(dir_cve, showWarnings=F)

#---
file <- paste0(dir_save, '00_lst_comm.Rdata')
load(file)


#########
permu <- 10000
#########

# ordinations ####
file_out <- paste0(dir_save, '02_lst_rda.Rdata')
if(file.exists(file_out)){
  load(file_out)
} else {

  # loop on communities
  
  depths <- c('top|deep','top','deep')
  out_rda <- NULL
  for(h in depths) {
    
    print(h)
      
    lst <- NULL
    for(i in n_comm) {
      
      mr <- lst_comm[[i]]$clr_nls$mr
      env <- lst_comm[[i]]$clr_nls$env
      
      # env <- env[row.names(mr),]
      env <- env[,c('site','moisture','depth','PiMiS','combi',
                    'pH','C','N','C_N','S','silt','clay','rh','om')]
      
      # ordination ####
      print(paste(i, 'model'))
      
      cond <- 'site'
      # Note: if many factor included in the condition (interaction or not), 
      # the only factor counting is the one having the highest resolution if hierarchical factors
      
      # RDA ----
      e <- na.omit(env[grepl(h, env$depth),])
      is.fact <- sapply(e, is.factor)
      m <- mr[row.names(e),]
      m <- m[,colSums(m) != 0]
      
      f <- formula(paste0(ifelse(h == 'top|deep', 'm~depth+', 'm~'), paste(names(e)[is.fact == F], collapse='+'),
                          '+Condition(', cond, ')', collapse=''))
      rda_full <- capscale(f, data=e, scale=T)
      
      lst_rda <- list(full=list(e=e, m=m, f=f, mod=rda_full))
      
      # retreive a more parsimonious model ---
      # Note: the parsimonious model did not depend of the order of the variable in the input model
      
      print('ordistep on full model')
      
      file <- paste0(dir_save, '02_rda_parsi_clr_nls_', i, '_', h, '.Rdata')
      if(file.exists(file)){
        load(file)
      } else {
        set.seed(0)
        f <- formula(paste0(ifelse(h == 'top|deep', 'm~depth+rh+Condition(', 'm~rh+Condition('), cond, ')'))
        rda_parsi <- ordistep(capscale(f, data=e, scale=T), formula(rda_full), parallel=4, trace=F)
        save(rda_parsi, file=file)
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
      
      file <- paste0(dir_save, '02_lst_ordi_clr_nls_', i, '_', h, '.Rdata')
      if(file.exists(file)){
        load(file)
      } else {
        
        lst_ordi <- NULL
        for(j in 1:2){
          
          print(paste('tests on', names(lst_rda)[j]))
          
          e   <- lst_rda[[j]]$e
          f   <- lst_rda[[j]]$f
          mod <- lst_rda[[j]]$mod
          
          v <- rev(rev(attributes(mod$terminfo$terms)$term.labels)[-1])
          
          if(length(v) > 1){
            
            stat <- data.frame(NA)
            v_out <- 'None'
            f2 <- f
            v2 <- v
            
            while(is.na(stat[1,1])){ # cyaB top samples did not have enough samples for the number of variables selected in the parimonious model
              
              nrh <- v2 != 'rh'
              print(v_out <- v2[nrh][stat[nrh,'V2'] == min(stat[nrh,'V2'])])
              ind_min_var <- which(v2 == v_out)
              
              if(length(v_out)){
                v2 <- v2[-ind_min_var]
                f2 <- update(f2, paste('.~.-',v_out))
              }
              
              #---
              stat <- foreach(k=v2) %dopar% {
                
                m <- lst_rda[[j]]$m
                mod <- capscale(eval(parse(text=paste(as.character(f2)[c(2,1,3)], collapse=''))),
                                data=e, scale=T) # some fuck with the foreach environment
                                                 # but output are the same then if tested sequencially
                
                set.seed(0)
                
                upd <- update(mod, as.formula(paste('.~.-', k)))
                pv1 <- signif(anova(mod, upd, permutations=permu)$`Pr(>F)`[2], 2)
                va1 <- RsquareAdj(mod)[[1]] - RsquareAdj(upd)[[1]] # no R2 adj because Conditon(site)
                
                rda1 <- capscale(as.formula(paste('m~', k, '+Condition(', cond, ')')), data=e, scale=T)
                pv2 <- signif(anova(rda1, permutations=permu)$`Pr(>F)`[1], 2)
                va2 <- RsquareAdj(rda1)[[1]]
                
                return(list(c(pv1, va1, pv2, va2), mod))
              }
              
              #---
              mod <- stat[[length(stat)]][[2]]
              
              stat <- as.data.frame(t(sapply(stat, '[[', 1)))
              row.names(stat) <- v2
            
            }
            
            v <- v2
            
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
        
        nmds <- metaMDS(m, dist='euc', trace=0) # aitchinson distance = euclydian dist on CLR
        
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
        
      }
      
      # graf ----
      print('graf')
      for(jn in names(lst_ordi)[2]){

        j <- lst_ordi[[jn]]

        # cairo_ps(paste0(dir_cve, 'ordination_', jn, '_', i, '_', h, '.eps'), width=10, height=7)
        svg(paste0(dir_cve, 'ordination_', jn, '_', i, '_', h, '.svg'), width=10, height=7)
        lay(j$axes_names)

        # ordinations
        smp <- j$site
        var <- j$var
        env_in_ordi <- env[row.names(env) %in% row.names(smp),]
        
        print(c(nrow(smp), nrow(env_in_ordi)))
        
        coord <- rbind(smp, var)

        for(k in seq_along(fact_3)){

          fact <- fact_3[k]

          pal <- lst_palev[[fact]][env_in_ordi[[fact]]]

          #---
          plot(smp, xlim=range(coord[,1]), ylim=range(coord[,2]), xaxt='n', yaxt='n', col=pal, pch=NA)

          abline(v=0, h=0, lty=3)

          # variables
          rng_smp <- sapply(as.data.frame(smp), range, simplify='matrix')
          rng_var <- sapply(as.data.frame(var), range, simplify='matrix')

          divis <- max(rng_var/rng_smp)
          v2 <- var/divis

          arrows(0,0,v2[,1],v2[,2], length=0, lty=2)

          if(k < 3){
            max_smp <- max(abs(rng_var[,1]/divis))
            max_var <- max(abs(rng_var[,1]))
            axis(3, at=seq(-max_smp, max_smp, length.out=9), labels=round(seq(-max_var, max_var, length.out=9), 2), col=2)
          }

          if(k > 1){
            max_smp <- max(abs(rng_var[,2]/divis))
            max_var <- max(abs(rng_var[,2]))
            axis(4, at=seq(-max_smp, max_smp, length.out=9), labels=round(seq(-max_var, max_var, length.out=9), 2), col=2)
          }

          if(nrow(v2) == 1){
            text(v2[1,1], v2[1,2], labels=paste(row.names(v2)))
          } else {
            text(v2[,1:2], labels=paste(row.names(v2)))
          }

          # samples
          if(k > 1){
            axis(1)
          }

          if(k %% 2 == 1){
            axis(2)
          }

          ordispider(smp, env_in_ordi$combi, col='grey80')

          points(smp[,1], smp[,2], pch=19, col=pal)
        }

        #---
        leg()

        dev.off()
      }
      
      ###
      lst[[i]] <- list(lst_rda=lst_rda, lst_ordi=lst_ordi)
    }
    
    out_rda[[h]] <- lst
  }
  
  save(out_rda, file=file_out)
}

#---
file <- paste0(dir_cve, 'stat_ordi.csv')
if(file.exists(file)){
  file.remove(file)
}

for(i in names(out_rda)){
  for(j in names(out_rda[[i]])){
    for(k in names(out_rda[[i]][[j]]$lst_ordi)){
      write.table(c(i,j,k), file, T, F, '\t')
      write.table(out_rda[[i]][[j]]$lst_ordi[[k]]$stat, file, T, F, '\t')
    }
  }
}

# xtable
lapply(out_rda, function(x) {
  xtable(cbind.data.frame(lapply(x, function(y) {
    m <- y$lst_ordi$parsi$stat[,1:2]
    out <- NULL
    v <- c('rh','depth','pH','om','C','N','S','C_N','clay','silt')
    for(i in v){
      if(i %in% row.names(m)){
        o <- m[i,]
      } else {
        o <- c(NA,NA)
      }
      out <- rbind(out, o)
    }
    row.names(out) <- v
    return(out)
  })), digit=3)
})

# procrust rotation on the NMDS with common samples ####
print("Procrust")

cmbn <- combn(names(lst_comm), 2)

#---
# cairo_ps(paste0(dir_cve, 'procrust.eps'), 10, 14)
svg(paste0(dir_cve, 'procrust.svg'), 10, 14)
par(mfrow=c(4,3))

pvs <- apply(cmbn, 2, function(x){
  
  print(x)
  
  lst <- list(comm1=list(mr=lst_comm[[x[1]]]$clr_nls$mr),
              comm2=list(mr=lst_comm[[x[2]]]$clr_nls$mr))
  
  # get the common samples
  smps <- intersect(row.names(lst[[1]]$mr),row.names(lst[[2]]$mr))
  
  # do the NMDS
  for(i in seq_along(lst)){
    print(x[i])
    mr <- lst[[i]]$mr[smps,]
    lst[[i]][['nmds']] <- metaMDS(mr, distance='euc', trymax=100, trace=F)
  }
  
  # test the procrustes and plot
  pt <- protest(lst[[1]]$nmds, lst[[2]]$nmds, scale=T, symmetric=T, score='sites')
  plot(procrustes(lst[[1]]$nmds, lst[[2]]$nmds, scale=T, symmetric=T, score='sites'),
       main=paste(c(x, '\np =', pt$signif, 'SS =', signif(pt$ss, 2)), collapse=' '))
  
  return(c(pt$signif, pt$ss))
  
})

dimnames(pvs) <- list(c('pvs','ss'), apply(cmbn, 2, function(x) paste(x, collapse=':')))

dev.off()

# write table of pvs

pvs_new <- matrix(NA, nrow=length(n_comm), ncol=length(n_comm), dimnames=list(n_comm, n_comm))

for(i in 4:1){
  rng <- rev(rev(cumsum(4:i))[1:2])+c(1,0)
  rng <- ifelse(is.na(rng), 1, rng)
  
  seq <- rng[1]:rng[2]
  seq1 <- (1:i)+(5-i)
  
  pvs_new[(5-i),seq1] <- pvs[1,seq]
  pvs_new[seq1,(5-i)] <- pvs[2,seq]
}

xtable(pvs_new, digit=3)




#####




















