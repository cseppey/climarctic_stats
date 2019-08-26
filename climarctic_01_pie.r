#####
# climarctic preparation
#####

rm(list=ls())
gc()

require(foreach)
require(doSNOW) # e.g. makeSOCKcluster()
require(plotrix)

# prep cluster
cl <- makeSOCKcluster(2)

clusterEvalQ(cl, library(plotrix))
clusterEvalQ(cl, library(foreach))
clusterEvalQ(cl, library(doSNOW))

registerDoSNOW(cl)


# dir loads ####
dir_out  <- 'Projets/Climarctic/stats/MBC/out/'
dir_save <- paste0(dir_out, 'saves/') 
dir_pie  <- paste0(dir_out, '01_pie/')
dir.create(dir_pie, showWarnings=F)

source('bin/src/my_prog/R/pie_taxo.r')

#---
file <- paste0(dir_save, 'lst_comm.Rdata')
load(file)

# loop on raw and rraref
for(h in c('raw','rrf')){
  
  print(h)
  
  # make a 18S version without metazoa or embryophyceae
  p18S <- lst_comm$`18S_V4`[[h]]
  
  ind_ME <- p18S$taxo$V3 == 'Metazoa' | p18S$taxo$V4 == 'Embryophyceae'
  
  p18S$mr   <- p18S$mr[,ind_ME == F]
  p18S$ass  <- p18S$ass[ind_ME == F,]
  p18S$taxo <- p18S$taxo[ind_ME == F,]
  
  lst_comm$`18S_V4_no_ME`[[h]] <- p18S
  
  # make a 18S version without metazoa or embryophyceae
  p16S <- lst_comm$`16S_V1-3`[[h]]
  
  ind_cya <- p16S$taxo$V2 == 'Cyanobacteria'
  
  p16S$mr   <- p16S$mr[,ind_cya]
  p16S$ass  <- p16S$ass[ind_cya,]
  p16S$taxo <- p16S$taxo[ind_cya,]
  
  lst_comm$`16S_V1-3_cya`[[h]] <- p16S
  
  # pie ####
  foreach(i=names(lst_comm)[4], .verbose=T) %dopar% {
    
    cl2 <- makeSOCKcluster(2)
    
    clusterEvalQ(cl2, library(plotrix))
    
    registerDoSNOW(cl2)
    
    #---
    
    mr <- lst_comm[[i]][[h]]$mr
    taxo <- lst_comm[[i]][[h]]$taxo
    env <- lst_comm[[i]][[h]]$env

    # arrange the smp nb and seq nb for per fact and cross fact
    selec_smp1 <-list(Knud  =which(env$site     == 'Knudsenheia'),
                      Ossian=which(env$site     == 'Ossian'),
                      dry   =which(env$moisture == 'dry'),
                      medium=which(env$moisture == 'medium'),
                      wet   =which(env$moisture == 'wet'),
                      top   =which(env$depth    == 'T'),
                      deep  =which(env$depth    == 'D'))
    names(selec_smp1) <- paste0(names(selec_smp1), ' smp nb: ', lapply(selec_smp1, function(x) nrow(mr[x,])),
                                '\nseq nb: ', lapply(selec_smp1, function(x) sum(mr[x,])))
    
    #---
    selec_smp2 <- factor(paste(env$moist_in_site, env$depth, sep='_'))
    lev <- levels(selec_smp2)
    selec_smp2 <- as.character(selec_smp2)
    
    for(j in lev){
      cond <- strsplit(j, '_')[[1]]
      ind <- which(env$moisture == cond[1] & env$site == cond[2] & env$depth == cond[3])
      selec_smp2[ind] <- paste0(j, '\nsmp nb: ', length(ind), ' seq nb: ', sum(mr[ind,]))
    }
    
    selec_smp2 <- as.factor(selec_smp2)
    
    # arrange the layout for the per smp
    lay <- NULL
    ind=1
    for(j in 1:108){
      if(j %in% as.numeric(substr(row.names(mr), 2, nchar(row.names(mr))))){
        lay <- c(lay, ind)
        ind <- ind+1
      } else {
        lay <- c(lay, 0)
      }
    }
    
    # arg pie
    lst_arg_pie <- list(per_fact  =list(selec_smp=selec_smp1,
                                        mat_lay=matrix(c(0,1,2,8, 3:5,8, 0,6,7,8), nrow=3, byrow=T),
                                        wdt_lay=c(1,1,1,3), hei_lay=c(rep(1.1, 3)),
                                        wdt=15, hei=7),
                        cross_fact=list(selec_smp=selec_smp2,
                                        mat_lay=matrix(c(1:4,13, 5:8,13, 9:12,13), nrow=3, byrow=T),
                                        wdt_lay=c(1,1,1,1, 3), hei_lay=c(rep(1.1, 3)),
                                        wdt=15, hei=7),
                        per_smp   =list(selec_smp=factor(paste0(row.names(mr), ' nb seq: ', rowSums(mr), '\n', env$moist_in_site, '_', env$depth)),
                                        mat_lay=cbind(matrix(lay, nrow=6, byrow=T), rep(max(lay)+1, 6))[,c(1:3,7:9,13:15, 4:6,10:12,16:18, 19)],
                                        wdt_lay=c(rep(1, 18), 3), hei_lay=c(rep(1.1,6)),
                                        wdt=57, hei=18)
                        )
    
    # loop on abundance and diversity
    foreach(j=c('abundance','richness'), .verbose=T) %dopar% {
      
      if(j == 'richness'){
        mr <- decostand(mr, 'pa')
      }
      
      for(k in seq_along(lst_arg_pie)){
        kn <- names(lst_arg_pie)[k]
        kl <- lst_arg_pie[[k]]
        
        pdf(paste0(dir_pie, 'pie_', h, '_', kn, '_', i, '_', j, '.pdf'), width=kl$wdt, height=kl$hei)
        
        pie_taxo(mr, taxo, ifelse(i == '16S_V1-3_cya', 4, 1):5, kl$selec_smp, mat_lay=kl$mat_lay, 
                 wdt_lay=kl$wdt_lay, hei_lay=kl$hei_lay, last_tax_text=F)
        
        dev.off()
      }
      
    }
    
  }
  
}



#























