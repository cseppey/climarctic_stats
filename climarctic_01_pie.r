#####
# climarctic preparation
#####

rm(list=ls())
gc()

require(foreach)
require(doSNOW) # e.g. makeSOCKcluster()
require(plotrix)

# prep cluster
cl <- makeSOCKcluster(4)

clusterEvalQ(cl, library(plotrix))

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

# make a 18S version without metazoa or embryophyceae
p18S <- lst_comm$`18S_V4`

ind_ME <- p18S$taxo$V3 == 'Metazoa' | p18S$taxo$V4 == 'Embryophyceae'

p18S$mr_sort <- p18S$mr_sort[,ind_ME == F]
p18S$ass_sort <- p18S$ass_sort[ind_ME == F,]
p18S$taxo_sort <- p18S$taxo_sort[ind_ME == F,]

lst_comm$`18S_V4_no_ME` <- p18S

# pie ####
foreach(i=names(lst_comm), .verbose=T) %dopar% {
  
  mr <- lst_comm[[i]]$mr_sort
  taxo <- lst_comm[[i]]$taxo_sort
  
  # selec smps
  e <- env_sort[row.names(mr),]
  
  selec_smp <- list(Knud  =which(e$site     == 'Knudsenheia'),
                    Ossian=which(e$site     == 'Ossian'),
                    dry   =which(e$moisture == 'dry'),
                    medium=which(e$moisture == 'medium'),
                    wet   =which(e$moisture == 'wet'),
                    top   =which(e$depth    == 'T'),
                    deep  =which(e$depth    == 'D'))
  
  # loop on abundance and diversity
  for(j in c('abundance','richness')){
    
    if(j == 'richness'){
      mr <- decostand(mr, 'pa')
    }
    
    # graf
    pdf(paste0(dir_pie, 'pie1_', i, '_', j, '.pdf'), width=15, height=7)
    
    pie <- pie_taxo(mr, taxo, 1:5, selec_smp, mat_lay=matrix(c(0,1,2,8, 3:5,8, 0,6,7,8), nrow=3, byrow=T),
                    wdt_lay=c(1,1,1,3), hei_lay=c(rep(1.1, 3)), last_tax_text=F)
    
    dev.off()
    
    #---
    selec_smp2 <- factor(apply(e[,c("moist_in_site",'depth')], 1, function(x) paste(x, collapse='_')))
    
    pdf(paste0(dir_pie, 'pie2_', i, '_', j, '.pdf'), width=15, height=7)
    
    pie <- pie_taxo(mr, taxo, 1:5, selec_smp2, mat_lay=matrix(c(1:4,13, 5:8,13, 9:12,13), nrow=3, byrow=T),
                    wdt_lay=c(1,1,1,1, 3), hei_lay=c(rep(1.1, 3)), last_tax_text=F)
    
    dev.off()
    
  }
  
}



#























