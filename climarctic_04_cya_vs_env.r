#####
# climarctic preparation
#####

rm(list=ls())
gc()

require(foreach)
require(doSNOW) # e.g. makeSOCKcluster()
require(PMCMR)
require(multcompView)

# prep cluster
cl <- makeSOCKcluster(2)

registerDoSNOW(cl)


# dir loads ####
dir_out  <- 'Projets/Climarctic/stats/MBC/out/'
dir_save <- paste0(dir_out, 'saves/') 
dir_cya  <- paste0(dir_out, '04_cya/')
dir.create(dir_cya, showWarnings=F)

#---
file <- paste0(dir_save, '00_lst_comm.Rdata')
load(file)
file <- paste0(dir_save, '01_lst_comm.Rdata')
load(file)

fact_3 <- c('depth','site','moisture')

# loop on raw and rraref
for(h in c('raw','rrf')){
  
  print(h)
  
  mr   <- lst_comm$`cyaB_cl`[[h]]$mr
  ass  <- lst_comm$`cyaB_cl`[[h]]$ass
  taxo <- lst_comm$`cyaB_cl`[[h]]$taxo
  env  <- lst_comm$`cyaB_cl`[[h]]$env
  
  # mr_16S   <- lst_comm$`16S_V1-3`[[h]]$mr
  # taxo_16S <- lst_comm$`16S_V1-3`[[h]]$taxo
  
  for(i in c('abundance','richness')){
    
    if(i == 'richness'){
      mr <- decostand(mr,'pa')
      # mr_16S <- decostand(mr_16S,'pa')
    }

    rs <- rowSums(mr)
    # rs_16S <- rowSums(mr_16S)
    
    # df frac cya vs 16S and cya grp
    cya_grp <- c('Gloeobacterales','Leptolyngbyales','Nostocales','Cyanobacteria')
    
    df_frac_sgrp <- data.frame(#frac=rowSums(mr) / rowSums(mr_16S),
                               sapply(cya_grp, function(x) rowSums(mr[,grep(x, ass$taxo)])))
    
    #---
    pdf(paste0(dir_cya, 'frac_sgrp_cyaB_', h, '_', i, '.pdf'), width=11, height=11)
    par(mfrow=c(3,2))

    for(jn in names(df_frac_sgrp)){
      j <- df_frac_sgrp[[jn]]
      
      for(k in fact_3){
        
        # make the interaction factors
        fact_inter <- fact_3[fact_3 != k]
        e <- factor(apply(env[fact_inter], 1, function(x) paste(x, collapse='_')))
        
        dfl <- data.frame(envk=env[[k]], e=e)
        for(ln in names(dfl)){
          
          l <- dfl[[ln]]
          lev <- levels(l)
          
          # test per factor and interaction      
          pvs <- posthoc.kruskal.nemenyi.test(j, l)$p.value
          pvs <- as.dist(cbind(rbind(rep(NA,ncol(pvs)), pvs), rep(NA, nrow(pvs)+1)))
          attributes(pvs)$Labels <- lev
          mcl <- multcompLetters(pvs)$Letters
          
          # measure the number of smp and seq per grp
          nb_smp <- table(l[j != 0])
          if(jn == 'frac'){
            nb_seq <- tapply(rs, list(l), sum)
          } else {
            nb_seq <- tapply(j, list(l), sum)
          }
          # nb_seq_16S <- tapply(rs_16S, list(l), sum)
          
          # graf args
          if(ln == 'e'){
            main <- paste(c('interaction', fact_inter), collapse=' ')
            names <- F
          } else {
            main <- k
            names <- lev
          }
          
          # graf
          boxplot(j~l, ylim=range(j)*c(1,1.5), xlim=c(0,length(lev)+0.5),
                  ylab=paste(i, ifelse(jn == 'frac', 'fraction cyano vs 16S', paste('of', jn))), 
                  main=main, names=names)
          usr <- par('usr')
          
          if(ln == 'e'){
            mtext(gsub('_','\n',lev), 1, 1.5, at=1:length(lev), cex=0.5)
          }
          
          # info
          x <- 1:length(lev)
          y <- usr[3]+diff(usr[3:4])*seq(0.95,0.8,length.out=4)
          text(x, y[1], labels=mcl)
          text(x, y[2], labels=nb_smp)
          text(x, y[3], labels=nb_seq)
          if(jn == 'frac'){
            text(x, y[4], labels=nb_seq_16S)
            text(usr[1]+diff(c(usr[1], 1))*0.5, y, c('grp','smp', paste(c('cya', '16S'), ifelse(i == 'abundance','seq','OTU'))), adj=c(0.5,0.5))
          } else {
            text(usr[1]+diff(c(usr[1], 1))*0.5, y[1:3], c('grp','smp', paste('nb', ifelse(i == 'abundance','seq','OTU'))), adj=c(0.5,0.5))
          }
        }
        
      }
    }
    
    dev.off()
    
  }
  
}

#####
# top10

mr  <- lst_comm$cyaB_cl$rrf$mr
ass <- lst_comm$cyaB_cl$rrf$ass
env <- lst_comm$cyaB_cl$rrf$env

print(ass[names(tail(sort(colSums(mr)), n=10)),])

lst <- list(otu_tb_cyaB_rrf          = cbind.data.frame(ass$taxo, t(mr)),
            otu_tb_cyaB_rrf_site     = cbind.data.frame(ass$taxo, t(apply(mr, 2, function(x) tapply(x, list(env$site), sum)))),
            otu_tb_cyaB_rrf_moisture = cbind.data.frame(ass$taxo, t(apply(mr, 2, function(x) tapply(x, list(env$moisture), sum)))),
            otu_tb_cyaB_rrf_depth    = cbind.data.frame(ass$taxo, t(apply(mr, 2, function(x) tapply(x, list(env$depth), sum)))),
            otu_tb_cyaB_rrf_combi    = cbind.data.frame(ass$taxo, t(apply(mr, 2, function(x) tapply(x, list(env$combi), sum)))))

file.fa <- 'Projets/Climarctic/stats/MBC/out/04_cya/cyaB_clean.fa'

if(file.exists(file.fa)){file.remove(file.fa)}
for(i in 1:nrow(ass)){
  write.table(paste0('>', row.names(ass)[i], ';', ass$taxo[i]), file.fa, T, F, row.names=F, col.names=F)
  write.table(ass$seq[i], file.fa, T, F, row.names=F, col.names=F)
}

save(lst, file='Projets/Climarctic/stats/MBC/out/04_cya/otu_tab_191004.Rdata')

for(i in names(lst)){
  write.table(lst[[i]], paste0('Projets/Climarctic/stats/MBC/out/04_cya/', i, '.csv'), quote=F, sep='\t')
}

#####























