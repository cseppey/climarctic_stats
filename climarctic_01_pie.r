#####
# climarctic pie charts
#####

rm(list=ls())
gc()

require(foreach)
require(doSNOW) # e.g. makeSOCKcluster()
require(plotrix)


# dir loads ####
dir_out  <- 'Projets/Climarctic/stats/MBC/out/'
dir_save <- paste0(dir_out, 'saves/') 
dir_pie  <- paste0(dir_out, '01_pie/')
dir.create(dir_pie, showWarnings=F)

source('bin/src/my_prog/R/pie_taxo.r')

#---
file <- paste0(dir_save, '00_lst_comm.Rdata')
load(file)

# loop on raw and rraref
lst_pie <- NULL
for(h in c('cmb','raw','rrf2','rrf')){
    
  print(h)

  # pie ####
  for(i in names(lst_comm)){

    print(i)
    
    tax_lev <- switch(i,
                      "01_16S_V1-3"   = 1:5,
                      "02_18S_V4"     = 1:5,
                      "03_pmoA_mb661" = 4:7,
                      "04_pmoA_A682"  = 4:7,
                      "05_ITS2"       = 1:5,
                      "06_phoD"       = 1:4,
                      "07_nifH"       = 1:4,
                      "08_cyaB"       = 4:6,
                      "09_nirS"       = 1:4)
    
    #---
    mr <- lst_comm[[i]][[h]]$mr
    taxo <- lst_comm[[i]][[h]]$taxo
    env <- lst_comm[[i]][[h]]$env
    
    # prep cluster
    cl <- makeSOCKcluster(2)
    
    clusterEvalQ(cl, library(plotrix))
    clusterEvalQ(cl, library(foreach))
    clusterEvalQ(cl, library(doSNOW))
    
    registerDoSNOW(cl)
    
    # loop on abundance and diversity
    lst_pie1 <- foreach(j=c('abundance','richness'), .verbose=T) %dopar% {
      
      # prep cluster
      cl2 <- makeSOCKcluster(2)
      
      clusterEvalQ(cl2, library(plotrix))
      clusterEvalQ(cl2, library(foreach))
      clusterEvalQ(cl2, library(doSNOW))
      
      registerDoSNOW(cl2)
      
      #-------------
      if(j == 'richness'){
        mr <- decostand(mr, 'pa')
      }
      
      #---
      # arrange the smp nb and seq nb for per fact and cross fact
      selec_smp1 <-list(Knud  =which(env$site     == 'Knudsenheia'),
                        Ossian=which(env$site     == 'Ossian'),
                        dry   =which(env$moisture == 'dry'),
                        medium=which(env$moisture == 'medium'),
                        wet   =which(env$moisture == 'wet'),
                        top   =which(env$depth    == 'top'),
                        deep  =which(env$depth    == 'deep'))
      
      if(j == 'abundance'){
        nb_seq_otu <- parLapply(cl2, selec_smp1, function(x, mr=mr) if(length(x)){return(round(sum(mr[x,])))}, mr)
      } else {
        nb_seq_otu <- parLapply(cl2, selec_smp1, function(x, mr=mr) if(length(x)){return(length(which(colSums(mr[x,]) != 0)))}, mr)
      }
      
      names(selec_smp1) <- paste0(names(selec_smp1), ' smp nb: ', 
                                  parLapply(cl2, selec_smp1, function(x, mr=mr) nrow(mr[x,]), mr),
                                  ifelse(j == 'abundance', '\nseq nb: ', '\nOTU nb: '), 
                                  nb_seq_otu)
      
      #---
      selec_smp2 <- factor(paste(env$moist_in_site, env$depth, sep='_'))
      lev <- apply(expand.grid(levels(env$moist_in_site), levels(env$depth)),
                                         1, function(x) paste(x, collapse='_'))
      selec_smp2 <- as.character(selec_smp2)
      
      rs <- rowSums(mr)
        
      ltot <- NULL
      for(k in lev){
        cond <- strsplit(k, '_')[[1]]
        ind <- which(env$moisture == cond[1] & env$site == cond[2] & env$depth == cond[3])
        ltot <- c(ltot, paste0(k, '\nsmp nb: ', length(ind),
                               ifelse(j == 'abundance', ' seq nb: ', ' OTU nb: '), round(sum(rs[ind]))))
        if(length(ind)){
          selec_smp2[ind] <- ltot[length(ltot)]
        }
      }
      
      selec_smp2 <- factor(selec_smp2, levels=ltot)
      
      # arrange the layout for the per smp
      lay <- NULL
      ind=1
      for(k in 1:108){
        if(k %in% as.numeric(substr(row.names(mr), 2, nchar(row.names(mr))))){
          lay <- c(lay, ind)
          ind <- ind+1
        } else {
          lay <- c(lay, 0)
        }
      }
      
      mat_per_smp <- matrix(lay, ncol=6, byrow=T)
      mat_per_smp <- rbind(mat_per_smp[,1:3],mat_per_smp[,4:6])
      
      mat_per_smp2 <- NULL
      for(k in 1:12){
        mat_per_smp2 <- cbind(mat_per_smp2, mat_per_smp[(k*3-2):(k*3),])
      }
      mat_per_smp <- NULL
      for(k in 1:4){
        mat_per_smp <- rbind(mat_per_smp, mat_per_smp2[,(k*9-8):(k*9)])
      }
      mat_per_smp <- cbind(mat_per_smp[1:6,], mat_per_smp[7:12,])
      rm(mat_per_smp2)
      
      # arg pie
      lst_arg_pie <- list(tot       =list(selec_smp=factor(rep(paste0('top\nsmp_nb: ', nrow(mr), ' seq nb: ', sum(mr)), nrow(mr))),
                                          mat_lay=matrix(c(0,1,2,0), nrow=1),
                                          wdt_lay=c(0.1,1,2.5,0.1), hei_lay=c(1.5),
                                          wdt=9, hei=6),
                          per_fact  =list(selec_smp=selec_smp1,
                                          mat_lay=matrix(c(0,1,2,8, 3:5,8, 0,6,7,8), nrow=3, byrow=T),
                                          wdt_lay=c(1,1,1,3), hei_lay=c(rep(1.1, 3)),
                                          wdt=15, hei=7),
                          cross_fact=list(selec_smp=selec_smp2,
                                          mat_lay=matrix(c(1,2,7,8,13, 3,4,9,10,13, 5,6,11,12,13), nrow=3, byrow=T),
                                          wdt_lay=c(1,1,1,1, 3), hei_lay=c(rep(1.1, 3)),
                                          wdt=15, hei=7),
                          per_smp   =list(selec_smp=factor(paste0(row.names(mr), 
                                                                  ifelse(j == 'abundance', ' nb seq: ', ' OTU nb: '),
                                                                  rowSums(mr), '\n', env$moist_in_site, '_', env$depth)),
                                          mat_lay=cbind(mat_per_smp, rep(max(lay)+1, nrow(mat_per_smp))),
                                          wdt_lay=c(rep(1, 18), 3), hei_lay=c(rep(1.1,6)),
                                          wdt=57, hei=18)
      )
      
      #---
      lst_pie2 <- foreach(k=rev(seq_along(lst_arg_pie)), .verbose=T) %dopar% {
        kn <- names(lst_arg_pie)[k]
        kl <- lst_arg_pie[[k]]

        source('bin/src/my_prog/R/pie_taxo.r')
        
        # pdf(paste0(dir_pie, 'pie_', h, '_', kn, '_', i, '_', j, '.pdf'), width=kl$wdt, height=kl$hei)
        cairo_ps(paste0(dir_pie, 'pie_', h, '_', kn, '_', i, '_', j, '.eps'), width=kl$wdt, height=kl$hei)
        
        pie <- pie_taxo(decostand(mr, 'total'), taxo, tax_lev, kl$selec_smp, mat_lay=kl$mat_lay, cex=1.2, box=F, thresh=0.02,
                        wdt_lay=kl$wdt_lay, hei_lay=kl$hei_lay, info_perc=F, rshift=0.0, last_tax_text=F)
        
        dev.off()
        
        return(pie)
      }
      
      names(lst_pie2) <- rev(names(lst_arg_pie))
      
      #---
      stopCluster(cl2)
      
      return(lst_pie2)
    }
    
    stopCluster(cl)
    
    names(lst_pie1) <- c('abundance','richness')
    
    lst_pie[[h]][[i]] <- lst_pie1
  }
  
}

#---
file <- paste0(dir_save, '01_lst_comm.Rdata')
save(lst_comm, file=file)
file <- paste0(dir_save, '01_lst_pie.Rdata')
save(lst_pie, file=file)

#























