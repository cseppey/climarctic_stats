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

# file <- paste0(dir_save, '03_lst_comm.Rdata')
# if(file.exists(file)){
#   load(file)
# }

#########
permu <- 1000
#########


fact_3 <- c('site','moisture','depth')

# loop on communities
for(i in names(lst_comm)){
  stat2 <- NULL
  for(j in 'rrf2'){#names(lst_comm[[i]])){

    if(is.null(lst_comm[[i]][[j]])){
      next
    }
    
    env <- lst_comm[[i]][[j]]$env
    env <- env[,c('site','moisture','depth','plot_in_moist_in_site',
                  'pH','C_N','S','NO3','NH4','P_labile',
                  'silt','clay','rh','om',"CH4","CO2","N2O")]
    
    # transfo ####
    print(paste(i, j, 'transfo'))
    
    mr <- lst_comm[[i]][[j]]$mr
    mr <- mr[rowSums(mr) != 0,colSums(mr) != 0]
    env <- env[row.names(mr),]
    
    # transfo <- c('raw','hell','log')
    transfo <- c('log')
    
    if('lst_trsf_bc' %in% names(lst_comm[[i]][[j]]) == F){
      lst_trsf <- foreach(k=transfo) %dopar% {
        if(k == 'raw'){
          return(list(mr=mr, bc=vegdist(mr)))
        } else {
          m <- decostand(mr, k)
          return(list(mr=m, bc=vegdist(m)))
        }
      }
      names(lst_trsf) <- transfo
      lst_trsf <- append(lst_trsf, list(env))
      names(lst_trsf)[length(lst_trsf)] <- 'env'
      
      lst_comm[[i]][[j]][['lst_trsf_bc']] <- lst_trsf
    }
    
    # model ####
    print(paste(i, j, 'model'))
    
    env <- lst_comm[[i]][[j]]$lst_trsf_bc$env
    
    res_tst <- foreach(k=transfo, .verbose=T) %dopar% {
      
      bc <- lst_comm[[i]][[j]]$lst_trsf_bc[[k]]$bc
      
      e <- na.omit(env)
      is.fact <- sapply(e, is.factor)
      b <- as.dist(as.matrix(bc)[row.names(e),row.names(e)])
      
      if(nrow(e) >= ncol(e)-1){
      
        # rda
        formu <- formula(paste0('b~depth+', paste(names(e)[is.fact == F], collapse='+'), '+Condition(site)', collapse=''))
        rda <- capscale(formu, data=e, scale=T)
        
        # retreive a more parcimonious model
        ### note: the parcimonious model did not depend of the order of the variable in the input model
        set.seed(0)
        rda_parci <- ordistep(capscale(b~depth+rh+Condition(site), data=e, scale=T), formula(rda), trace=F)
        v <- attributes(rda_parci$terminfo$terms)$term.labels
        
        e <- na.omit(env[,v])
        b <- as.dist(as.matrix(bc)[row.names(e),row.names(e)])
        rda_parci <- capscale(formula(paste(c('b~depth+rh', names(e)[grepl('depth|rh|site', names(e)) == F],
                                              'Condition(site)'), collapse='+')), data=e)
        
        # graf
        pdf(paste0(dir_cve, 'RDA_', i, '_', j, '_', k, '.pdf'), width=10, height=7)
        layout(matrix(c(1,3,4, 1,5,6, 0,2,2), nrow=3, byrow=T), width=c(0.2,1,1), height=c(1,1,0.3))
        
        stat <- NULL
        for(l in list(rda, rda_parci)){
          
          par(mar=rep(0,4), oma=rep(1,4))
          
          # coord prep
          s <- summary(l)
          
          site <- s$sites[,1:2]
          
          v <- unlist(strsplit(gsub(' ', '', as.character(s$call$formula[[3]][2])), '+', fixed=T))
          var <- s$biplot
          var <- var[row.names(var) %in% v,1:2]
          
          e <- env[row.names(site),]
          b <- as.dist(as.matrix(bc)[row.names(e),row.names(e)])
          
          coord <- rbind(site, var)
          
          # axis
          vari <- s$cont$importance[2,1:2]
          vari <- paste0(names(vari), ' ', signif(vari*100,2), '%')
          for(m in 2:1){
            plot.new()
            text(0.5,0.5, labels=vari[m], srt=ifelse(m == 1, 0, 90))
          }
          
          # ordinations
          for(m in seq_along(fact_3)){
            
            mn <- fact_3[m]
            
            plot(site, xlim=range(coord[,1], na.rm=T), ylim=range(coord[,2], na.rm=T), xaxt='n', yaxt='n',
                 col=lst_palev[[mn]]$pal[e[[mn]]], pch=NA)
            if(m != 1){
              axis(1)
            }
            if(m != 2){
              axis(2)
            }
            
            ordispider(site, factor(paste(e$depth, e$plot_in_moist_in_site)), col='grey80')
            
            text(site, labels=row.names(site), col=lst_palev[[mn]]$pal[e[[mn]]])
            
            abline(v=0,h=0,lty=3)
            
            arrows(0,0,var[,1],var[,2], length=0, lty=2, col='grey70')
            text(var[,1:2], labels=paste(row.names(var)))
            
          }
          
          # legend ---
          plot.new()
          
          # factors colors
          legend(0.5,0.5, legend=c(unlist(sapply(lst_palev, '[[', 2))), bty='n', xjust=0.5, yjust=0.5,
                 pch=19, col=unlist(sapply(lst_palev, '[[', 1)))
          
          # test each variable
          vtt <- rev(rev(attributes(l$terms)$term.labels)[-1])
          pvs <- matrix(NA, nrow=ncol(env)-3, ncol=4, dimnames=list(names(env)[-c(1,2,4)], c('pv_strict','va_strict', 'pv_all','va_all')))
          for(m in vtt){
            set.seed(0)
            
            upd <- update(l, as.formula(paste('.~.-', m)))
            pv1 <- signif(anova(l, upd, parallel=4, permutations=permu)$`Pr(>F)`[2], 2)
            va1 <- RsquareAdj(l)[[1]] - RsquareAdj(upd)[[1]]
            
            rda1 <- capscale(as.formula(paste('b~', m)), data=e, scale=T)
            pv2 <- signif(anova(rda1, permutations=permu)$`Pr(>F)`[1], 2)
            va2 <- RsquareAdj(rda1)[[1]] 
            
            pvs[m,] <- c(pv1, va1, pv2, va2)
          }
          
          ### note: the AIC did not depend on the order of the variables in the model
          ### note: the RsqrtAdj is not available because of the Condition(site)
          stat <- cbind(stat, rbind(pvs, info=c(extractAIC(l)[2], RsquareAdj(l)[1], nrow(site),
                                                length(unlist(strsplit(as.character(l$call$formula[[3]][2]), '+', fixed=T))))))
          
        }
        
        dev.off()
        
        st <- apply(stat, 2, as.numeric)
        row.names(st) <- c(row.names(pvs), 'AIC')
        
        return(st)      
      } else {
        return(dim(e))
      }
    }
    
    names(res_tst) <- transfo
    stat2[[j]] <- res_tst

  }
  
  # out stat
  file <- paste0(dir_cve, 'stat_', i, '.csv')
  if(file.exists(file)){
    file.remove(file)
  }
  
  for(k in names(stat2)){
    write.table(k, file, T)
    for(l in names(stat2[[k]])){
      write.table(l, file, T)      
      write.table(stat2[[k]][[l]], file, T)
    }
    
  }
  
}

file <- paste0(dir_save, '03_lst_comm.Rdata')
save(lst_comm, file=file)

#




















