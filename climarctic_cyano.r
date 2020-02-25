#####
# climarctic preparation
#####

rm(list=ls())

# require(foreach)
# require(doSNOW) # e.g. makeSOCKcluster()
# require(PMCMR)
# require(multcompView)
# 
# # prep cluster
# cl <- makeSOCKcluster(2)
# 
# registerDoSNOW(cl)
# 
# 
# # dir loads ####
# dir_out  <- 'Projets/Climarctic/stats/MBC/out/'
# dir_save <- paste0(dir_out, 'saves/') 
# dir_cya  <- paste0(dir_out, '04_cya/')
# dir.create(dir_cya, showWarnings=F)
# 
# #---
# file <- paste0(dir_save, '00_lst_comm.Rdata')
# load(file)
# file <- paste0(dir_save, '01_lst_comm.Rdata')
# load(file)
# 
# fact_3 <- c('depth','site','moisture')
# #
# 
# # loop on raw and rraref ####
# for(h in c('raw','rrf', 'rrf2')){
#   
#   print(h)
#   
#   mr   <- lst_comm$`cyaB_cl`[[h]]$mr
#   ass  <- lst_comm$`cyaB_cl`[[h]]$ass
#   taxo <- lst_comm$`cyaB_cl`[[h]]$taxo
#   env  <- lst_comm$`cyaB_cl`[[h]]$env
#   
#   # mr_16S   <- lst_comm$`16S_V1-3`[[h]]$mr
#   # taxo_16S <- lst_comm$`16S_V1-3`[[h]]$taxo
#   
#   for(i in c('abundance','richness')){
#     
#     if(i == 'richness'){
#       mr <- decostand(mr,'pa')
#       # mr_16S <- decostand(mr_16S,'pa')
#     }
# 
#     rs <- rowSums(mr)
#     # rs_16S <- rowSums(mr_16S)
#     
#     # df frac cya vs 16S and cya grp
#     cya_grp <- c('Gloeobacterales','Leptolyngbyales','Nostocales','Cyanobacteria')
#     
#     df_frac_sgrp <- data.frame(#frac=rowSums(mr) / rowSums(mr_16S),
#                                sapply(cya_grp, function(x) rowSums(mr[,grep(x, ass$taxo)])))
#     
#     #---
#     pdf(paste0(dir_cya, 'frac_sgrp_cyaB_', h, '_', i, '.pdf'), width=11, height=11)
#     par(mfrow=c(3,2))
# 
#     for(jn in names(df_frac_sgrp)){
#       j <- df_frac_sgrp[[jn]]
#       
#       for(k in fact_3){
#         
#         # make the interaction factors
#         fact_inter <- fact_3[fact_3 != k]
#         e <- factor(apply(env[fact_inter], 1, function(x) paste(x, collapse='_')))
#         
#         dfl <- data.frame(envk=env[[k]], e=e)
#         for(ln in names(dfl)){
#           
#           l <- dfl[[ln]]
#           lev <- levels(l)
#           
#           # test per factor and interaction      
#           pvs <- posthoc.kruskal.nemenyi.test(j, l)$p.value
#           pvs <- as.dist(cbind(rbind(rep(NA,ncol(pvs)), pvs), rep(NA, nrow(pvs)+1)))
#           attributes(pvs)$Labels <- lev
#           mcl <- multcompLetters(pvs)$Letters
#           
#           # measure the number of smp and seq per grp
#           nb_smp <- table(l[j != 0])
#           if(jn == 'frac'){
#             nb_seq <- tapply(rs, list(l), sum)
#           } else {
#             nb_seq <- tapply(j, list(l), sum)
#           }
#           # nb_seq_16S <- tapply(rs_16S, list(l), sum)
#           
#           # graf args
#           if(ln == 'e'){
#             main <- paste(c('interaction', fact_inter), collapse=' ')
#             names <- F
#           } else {
#             main <- k
#             names <- lev
#           }
#           
#           # graf
#           boxplot(j~l, ylim=range(j)*c(1,1.5), xlim=c(0,length(lev)+0.5),
#                   ylab=paste(i, ifelse(jn == 'frac', 'fraction cyano vs 16S', paste('of', jn))), 
#                   main=main, names=names)
#           usr <- par('usr')
#           
#           if(ln == 'e'){
#             mtext(gsub('_','\n',lev), 1, 1.5, at=1:length(lev), cex=0.5)
#           }
#           
#           # info
#           x <- 1:length(lev)
#           y <- usr[3]+diff(usr[3:4])*seq(0.95,0.8,length.out=4)
#           text(x, y[1], labels=mcl)
#           text(x, y[2], labels=nb_smp)
#           text(x, y[3], labels=nb_seq)
#           if(jn == 'frac'){
#             text(x, y[4], labels=nb_seq_16S)
#             text(usr[1]+diff(c(usr[1], 1))*0.5, y, c('grp','smp', paste(c('cya', '16S'), ifelse(i == 'abundance','seq','OTU'))), adj=c(0.5,0.5))
#           } else {
#             text(usr[1]+diff(c(usr[1], 1))*0.5, y[1:3], c('grp','smp', paste('nb', ifelse(i == 'abundance','seq','OTU'))), adj=c(0.5,0.5))
#           }
#         }
#         
#       }
#     }
#     
#     dev.off()
#     
#   }
#   
# }
# 
# # top10 ####
# 
# mr  <- lst_comm$cyaB_cl$rrf$mr
# ass <- lst_comm$cyaB_cl$rrf$ass
# env <- lst_comm$cyaB_cl$rrf$env
# 
# print(ass[names(tail(sort(colSums(mr)), n=10)),])
# 
# lst <- list(otu_tb_cyaB_rrf          = cbind.data.frame(ass$taxo, t(mr)),
#             otu_tb_cyaB_rrf_site     = cbind.data.frame(ass$taxo, t(apply(mr, 2, function(x) tapply(x, list(env$site), sum)))),
#             otu_tb_cyaB_rrf_moisture = cbind.data.frame(ass$taxo, t(apply(mr, 2, function(x) tapply(x, list(env$moisture), sum)))),
#             otu_tb_cyaB_rrf_depth    = cbind.data.frame(ass$taxo, t(apply(mr, 2, function(x) tapply(x, list(env$depth), sum)))),
#             otu_tb_cyaB_rrf_combi    = cbind.data.frame(ass$taxo, t(apply(mr, 2, function(x) tapply(x, list(env$combi), sum)))))
# 
# file.fa <- 'Projets/Climarctic/stats/MBC/out/04_cya/cyaB_clean.fa'
# 
# if(file.exists(file.fa)){file.remove(file.fa)}
# for(i in 1:nrow(ass)){
#   write.table(paste0('>', row.names(ass)[i], ';', ass$taxo[i]), file.fa, T, F, row.names=F, col.names=F)
#   write.table(ass$seq[i], file.fa, T, F, row.names=F, col.names=F)
# }
# 
# save(lst, file='Projets/Climarctic/stats/MBC/out/04_cya/otu_tab_191004.Rdata')
# 
# for(i in names(lst)){
#   write.table(lst[[i]], paste0('Projets/Climarctic/stats/MBC/out/04_cya/', i, '.csv'), quote=F, sep='\t')
# }
# 
# # relabu 23 oct 2019 ####
# 
# mr <- lst_comm$`08_cyaB`$raw$mr
# fa <- lst_comm$`08_cyaB`$raw$fa
# 
# relabu <- decostand(mr, 'total')
# 
# write.table(cbind.data.frame(fa=fa, t(relabu)), paste0(dir_cya, 'cyaB_seq_relabu.csv'), quote=F, sep='\t')
# 
# write.table(cbind.data.frame(fa=fa, t(mr)), paste0(dir_cya, 'cyaB_seq_raw.csv'), quote=F, sep='\t')
# 
# 
# 
# # MEETING IN ROSTOCK ####
# dominant vs other datasets ####

# tree

require(ape)

for(h in c('GTRCAT','GTRGAMMA')){
  for(i in c('_trim_V1-V3','_no_V1-V3','-antar')){
    tree <- read.tree(paste0('Projets/Climarctic/bioinfo/Katya_cyaB/01_dominant_cya_vs_other_dataset/', h, '/RAxML_bestTree.out', i))
    pdf(paste0('Projets/Climarctic/stats/MBC/out/04_cya/phylo_climarctic_vs_other/', h, '/', i, '.pdf'), width=10, height=20)
    
    plus <- ifelse(i == '_trim_V1-3', 0, 1)
    
    plot(tree, tip.color=as.numeric(factor(sapply(strsplit(tree$tip.label, '_'), '[[', 1))) + plus,
         align.tip.label=T)
    
    dev.off()
  }
}


# match with alignment

blst <- read.table('Projets/Climarctic/bioinfo/Katya_cyaB/01_dominant_cya_vs_other_dataset/03_pwaln/out.blst')

names(blst) <- c('qname','sname','pid','length','mismatch','gapopen',
                 'qstart','qend','sstart','send','evalue','bitscore')

lin <- read.table('Projets/Climarctic/bioinfo/Katya_cyaB/01_dominant_cya_vs_other_dataset/01_msaln/tot.lin')

lin$V2 <- substr(as.character(lin$V2), 2, nchar(as.character(lin$V2)))

blst_mat <- as.matrix(blst)

file <- 'Projets/Climarctic/stats/MBC/out/04_cya/top10.csv'
if(file.exists(file)){file.remove(file)}

for(i in lin$V2){
  top10 <- blst$sname[blst$qname == i]
  for(j in top10){
    write.table(t(c(as.character(lin$V1[lin$V2 == j]),
                  blst_mat[blst_mat[,1] == i & blst_mat[,2] == j,])), 
                quote=F, append=T, file=file, row.names=F, col.names=F)
  }
}


# cyaB vs 16S ####

load(paste0('Projets/Climarctic/stats/MBC/out/saves/00_lst_comm.Rdata'))

# cyaB ---
# mr and ass

mr_ass_cyaB <- read.table('Projets/Climarctic/papers/Katya_cyanobact/data/Cyano_CyaB.csv', h=T, sep='\t', row.names=1)

ass_cya <- mr_ass_cyaB[,grep('X', names(mr_ass_cyaB))]
mr_cya <- as.data.frame(t(mr_ass_cyaB[,-grep('X', names(mr_ass_cyaB))]))

# 16S ---
# mr and ass

mr_ass_16S <- read.table('Projets/Climarctic/papers/Katya_cyanobact/data/Cyano_16S.csv', h=T, sep='\t', row.names=1)

ass_cyano <- mr_ass_16S[,grep('X', names(mr_ass_16S))]
mr_cyano <- as.data.frame(t(mr_ass_16S[,-grep('X', names(mr_ass_16S))]))

names(mr_cyano) <- sub('OTU_', 'X_', names(mr_cyano))

# richness and abds
lst_cya <- list(cya  =list(nb_seq_ini=rowSums(lst_comm$`08_16S_cyano`$raw$mr)[row.names(mr_cya)],
                           rich=specnumber(mr_cya), abds=rowSums(mr_cya)),
                cyano=list(nb_seq_ini=rowSums(lst_comm$`01_16S_bact`$raw$mr)[row.names(mr_cyano)],
                           rich=specnumber(mr_cyano), abds=rowSums(mr_cyano)))

# test and plot
par(mfrow=c(2,2))
for(i in names(lst)){
  for(j in c('rich','abds')){
    plot(lst_cya[[i]][[j]]~log(lst_cya[[i]]$nb_seq_ini), main=i, xlab='nb_seq_ini', ylab=j)
    abline(coef(lm(lst_cya[[i]][[j]]~log(lst_cya[[i]]$nb_seq_ini))))
    cor <- cor.test(lst_cya[[i]][[j]],log(lst_cya[[i]]$nb_seq_ini), method='spear')
    print(c(cor$p.value, cor$estimate))
  }
}


#########
# there is samples in Katya's cyaB that are not in the initial dataset
#########













#










