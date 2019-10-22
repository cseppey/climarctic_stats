#####
# climarctic preparation
#####

rm(list=ls())
gc()

require(foreach)
require(doSNOW) # e.g. makeSOCKcluster()
require(RColorBrewer) # brewer.pal 

# prep cluster
cl <- makeSOCKcluster(4)

registerDoSNOW(cl)


# dir loads ####
dir_in <- 'Projets/Climarctic/stats/MBC/in/'
dir_out <- 'Projets/Climarctic/stats/MBC/out/'
dir_prep <- paste0(dir_out, '00_prep/') 
dir_save <- paste0(dir_out, 'saves/') 
dir.create(dir_prep, showWarnings=F, recursive=T)
dir.create(dir_save, showWarnings=F, recursive=T)

# env ----
env_ini_fact <- read.table(paste0(dir_in, 'env_clim_grad.csv'), row.names=1)
env_ini_chim <- read.table(paste0(dir_in, 'env_clim_chim.csv'), h=T)

# reorganize env fact
rn <- row.names(env_ini_fact)
row.names(env_ini_fact) <- ifelse(as.numeric(rn) < 100, ifelse(as.numeric(rn) < 10, paste0('T00', rn), paste0('T0', rn)), paste0('T', rn))

# reorganize env_chim
as_char <- as.character(env_ini_chim$sample)
env_ini_chim$site  <- substr(as_char, 1, 1)
env_ini_chim$moist <- rep(gl(3,3, labels=c('dry','medium','wet')), 4)
env_ini_chim$depth <- gl(2, 18, labels=c('top','deep'))

env_ini_chim <- env_ini_chim[c(matrix(1:36, nrow=2, byrow=T)),]
env_ini_chim <- env_ini_chim[as.numeric(gl(36,3)),]

# the the interesting variables
env_tot <- cbind.data.frame(env_ini_fact[,c(2:5,7:11)], env_ini_chim[,c('pH','N','C','sand','silt','clay')])
names(env_tot)[1:9] <- c('site','moisture','plot','depth','quadrat','empty','fresh','dry','burn')

env_tot[,c('fresh','dry','burn')] <- sapply(env_tot[,c('fresh','dry','burn')], function(x) x-env_tot$empty)
env_tot$rh <- (env_tot$fresh-env_tot$dry) / env_tot$fresh
env_tot$om <- (env_tot$dry-env_tot$burn) / env_tot$dry
env_tot$C_N <- env_tot$C / env_tot$N

env_tot <- env_tot[,-grep(paste(c('empty','fresh','dry','burn'), collapse='|'), names(env_tot))]

env_tot$site <- gl(2, 54, labels=c('Knudsenheia','Ossian'))
env_tot$moisture <- ordered(gl(3,18,108, c('dry','medium','wet')))
env_tot$plot <- as.factor(env_tot$plot)
env_tot$depth <- rep(gl(2,3, labels=c('top','deep')), 18)
env_tot$depth <- factor(env_tot$depth, levels=c('top','deep'))

env_tot$moist_in_site <- factor(apply(env_tot[,c('moisture','site')], 1, function(x) paste(x, collapse='_')))
env_tot$plot_in_moist_in_site <- factor(apply(env_tot[,c('plot','moist_in_site')], 1, function(x) paste(x, collapse='_')))
env_tot$quad_in_plot_in_moist <- factor(apply(env_tot[,c('quadrat','plot_in_moist_in_site')], 1, function(x) paste(x, collapse='_')))
env_tot$combi <- factor(paste(env_tot$moist_in_site, env_tot$depth, sep='_'))

# palette ----
lev_site <- levels(env_tot$site)
lev_mois <- levels(env_tot$moisture)
lev_dept <- levels(env_tot$depth)

pal_dark2 <- brewer.pal(8, 'Dark2')
lst_palev <- list(site    =list(pal=pal_dark2[1:2], lev=lev_site),
                  moisture=list(pal=pal_dark2[3:5], lev=lev_mois),
                  depth   =list(pal=pal_dark2[6:7], lev=lev_dept)
)

# permu and primer names ----
permu <- 10000

#---
prim_names <- c('01_16S_V1-3','02_18S_V4','03_pmoA_mb661','04_pmoA_A682','05_ITS2','06_phoD','07_nifH', '08_cyaB','09_nirS')
ind_prim <- c(1:9)

# loop the primers ####
lst_comm <- NULL
for(i in ind_prim) {
  
  print(prim_names[i])
  
  id_plate <- as.character(i)
  if(i < 10){
    id_plate <- paste0('0', id_plate)
  }
  
  #---
  file <- paste0(dir_save, '00_ini_', id_plate, '.Rdata')
  if(file.exists(file)){
    load(file)
  } else {
    mr_ini  <- read.table(paste0(dir_in, 'from_cluster/', id_plate, '/', id_plate, '_clust.mr'), h=T)
    # ass_ini <- read.table(paste0(dir_in, 'from_cluster/', id_plate, '/', id_plate, '_clust.DB.wang.taxonomy'), row.names=1)
    fa_ini  <- read.table(paste0(dir_in, 'from_cluster/', id_plate, '/', id_plate, '_clust.fa'))
    # save(mr_ini, ass_ini, fa_ini, file=file)
    save(mr_ini, fa_ini, file=file)
  }
  
  # reorganize mr
  mr_tot <- mr_ini[grep('T|B', row.names(mr_ini)),]
  
  rs <- rowSums(mr_tot)
  ord <- order(rs)
  plot(rs[ord], col=as.numeric(grepl('B', row.names(mr_tot))[ord])+1, pch=19, main=prim_names[i])
  
  mr_tot <- mr_tot[grep('T', row.names(mr_tot)),]
  mr_tot <- mr_tot[,colSums(mr_tot) != 0]
  
  # reorganize fa
  n_fa <- as.character(fa_ini[seq(1,nrow(fa_ini), by=2),])
  n_fa <- substr(n_fa, 2, nchar(n_fa))
  
  fa_tot <- fa_ini[seq(2,nrow(fa_ini), by=2),]
  names(fa_tot) <- n_fa
  
  # # reorganize ass
  # ass_tot <- data.frame(taxo=ass_ini$V2, seq=fa_tot)
  # 
  # ass_tot$taxo <- gsub('[[:punct:]][[:digit:]]{2,3}[[:punct:]]{2}|;', '|', ass_tot$taxo)
  # 
  # ass_tot <- ass_tot[names(mr_tot),]
  # 
  # taxo_tot <- strsplit(as.character(ass_tot$taxo), '|' , fixed=T)
  # nb_lev <- length(taxo_tot[[1]])
  # 
  # taxo_tot <- matrix(unlist(taxo_tot), ncol=nb_lev, byrow=T)
  # row.names(taxo_tot) <- row.names(ass_tot)
  # 
  # taxo_tot <- as.data.frame(t(apply(taxo_tot, 1, function(x) {
  #   x <- gsub('Incertae_Sedis', 'X', x)
  #   x <- gsub('Unknown', 'unknown', x)
  #   
  #   ind_unc <- grep('^[[:lower:]]', x)
  #   for(i in ind_unc){
  #     x[i] <- ifelse(grepl('_X', x[i-1]), paste0(x[i-1], 'X'), paste0(x[i-1], '_X'))
  #   }
  #   return(x)
  # })))
  # 
  # # sort taxo
  # ord_taxo <- order(ass_tot$taxo)
  # 
  # ass_sort <- ass_tot[ord_taxo,]
  # mr_sort <- mr_tot[,ord_taxo]
  # taxo_sort <- taxo_tot[ord_taxo,]
  
  mr_sort <- mr_tot
  env_sort <- env_tot[row.names(mr_sort),]
  
  # rrarefy ----
  rs <- rowSums(mr_sort)
  thresh <- seq(min(rs), max(rs), length.out=20)
  thresh <- thresh[-length(thresh)]
  
  # nb seq
  plot(sort(rs), main=prim_names[i])
  abline(h=thresh)
  
  # perc lost
  perc_lost <- foreach(j=thresh, .combine=cbind, .verbose=T) %dopar% {
    if(nrow(mr_sort[rs >= j,]) > 1){
      set.seed(0)
      mr_rrf <- rrarefy(mr_sort[rs >= j,], j)
      mr_rrf <- mr_rrf[,colSums(mr_rrf) != 0]

      return(c(sum(mr_rrf) / sum(mr_sort), ncol(mr_rrf) / ncol(mr_sort), nrow(mr_rrf) / nrow(mr_sort)))
    }
  }

  dimnames(perc_lost) <- list(c('seq','otu','smp'), paste0('t', thresh[1:ncol(perc_lost)]))
  
  l <- NULL
  op <- NULL
  for(j in c(0.5, 0.75, 1)){
    optimum <- ceiling(mean(apply(perc_lost[-3,], 1, function(x) thresh[which(x == max(x, na.rm=T))])))*j
    op <- c(op, optimum)
    
    # which smps did we loose
    set.seed(0)
    mr_rrf <- rrarefy(mr_sort[rs >= optimum,], optimum)
    mr_rrf <- mr_rrf[,colSums(mr_rrf) != 0]
    
    e <- env_tot[row.names(mr_rrf),]
    
    ll <- NULL
    for(k in c('site','moisture','depth')){
      ll[[k]] <- paste(c(j, k, signif(table(e[[k]]) / table(env_tot[[k]]), 2)), collapse=' ')
    }
    
    l <- cbind(l, ll)
  }
  
  #---
  # pdf(paste0(dir_prep, 'rraref_optimum_', prim_names[i], '.pdf'))
  cairo_ps(paste0(dir_prep, 'rraref_optimum_', prim_names[i], '.eps'))
  
  plot(NA, xlim=range(thresh), ylim=c(0,1), xlab='thresh', ylab='percentage',
       main=paste('seq ini:', sum(mr_sort), 'otu ini:', ncol(mr_sort), 'smp ini:', nrow(mr_sort)))
  
  for(j in 1:nrow(perc_lost)){
    lines(thresh[1:ncol(perc_lost)], perc_lost[j,], col=j)
  }
  
  abline(v=op, lty=3)

  legend('topright', legend=c('sequence','OTU','sample'), text.col=1:3, bty='n')
  
  dev.off()
  
  #---
  set.seed(0)
  mr_rrf <- rrarefy(mr_sort[rs >= op[3],], optimum)
  mr_rrf <- as.data.frame(mr_rrf[,colSums(mr_rrf) != 0])
  
  env_rrf <- env_tot[row.names(mr_rrf),]
  
  # ass_rrf <- ass_sort[names(mr_rrf),]
  # taxo_rrf <- taxo_sort[names(mr_rrf),]

  #---  
  set.seed(0)
  mr_rrf2 <- rrarefy(mr_sort[rs >= op[1],], optimum)
  mr_rrf2 <- as.data.frame(mr_rrf2[,colSums(mr_rrf2) != 0])
  
  env_rrf2 <- env_tot[row.names(mr_rrf2),]
  
  # ass_rrf2 <- ass_sort[names(mr_rrf2),]
  # taxo_rrf2 <- taxo_sort[names(mr_rrf2),]
  
  lst_comm[[prim_names[i]]] <- list(raw=list(env=env_sort, mr=mr_sort),#, ass=ass_sort, taxo=taxo_sort),
                                    rrf=list(env=env_rrf,  mr=mr_rrf),#,  ass=ass_rrf,  taxo=taxo_rrf),
                                    rrf2=list(env=env_rrf2,  mr=mr_rrf2))#,  ass=ass_rrf2,  taxo=taxo_rrf2))

}

#---
file <- paste0(dir_save, '00_lst_comm.Rdata')
save(lst_comm, env_tot, lst_palev, permu, file=file)

file <- paste0(dir_save, '00_env_tot.Rdata')
save(env_tot, file=file)

# check length distro ####

lst_lim <- list(`16S_V1-3`=c(450,570),
                `18S_V4`=c(380,480),
                pmoA_mb661=c(490,530),
                pmoA_A682=c(510,550),
                ITS=c(270,500),
                phoD=c(330,450),
                nifH=c(310,450),
                cyaB=c(390,470),
                nirS=c(380,460))

pdf(paste0(dir_prep, 'lgt_distro.pdf'), width=15, height=15)
par(mfrow=c(3,3))

for(i in ind_prim){
  lst_file <- list.files(paste0('Projets/Climarctic/bioinfo/archive/191019_all_clustering/0', i, '/filter_test/distro_test'), full.names=T)
  
  lgt_dis <- NULL
  for(j in lst_file){
    lgt_dis <- c(lgt_dis, unlist(read.table(j)))
  }
  
  hist(lgt_dis, breaks=50, main=names(lst_lim)[i])
  abline(v=lst_lim[[i]])
}

dev.off()

# check bioinfo output ####

pdf(paste0(dir_prep, 'bioinfo_check.pdf'), width=15, height=20)
par(mfrow=c(5,4))

ra_tot <- NULL
for(i in seq_along(prim_names[ind_prim])){
  out_bf <- read.table(paste0('Projets/Climarctic/bioinfo/archive/191019_all_clustering/filter_test/cnt_output/cnt_output', i),
                       h=T, sep='\t', row.names=1)
  
  out_bf <- out_bf[order(row.names(out_bf)),]

  B <- switch(i,
              190,
              127,
              NULL,
              NULL,
              121:129,
              107:109,
              72:74,
              82:84,
              96:99)
  
  cs <- colSums(out_bf, na.rm=T)  
  out_bf <- rbind.data.frame(out_bf, cs)
  
  lst <- list(raw=out_bf, relabu=decostand(out_bf, 'max', 1, na.rm=T))
  
  ra_tot <- cbind(ra_tot, cbind(cs, decostand(cs, 'max')))
  
  for(jn in names(lst)){
    j <- lst[[jn]]
    
    plot(NA, xlim=c(1,ncol(j)), xaxt='n', xlab='', 
         ylim=range(j, na.rm=T), ylab=ifelse(jn == 'raw', 'nb seq', '% nb_seq'),
         log=ifelse(jn == 'raw', 'y', ''), main=paste(jn, prim_names[i]))
    axis(1, at=1:ncol(j), labels=names(j))
    
    for(k in 1:nrow(j)){
      kb <- k %in% B
      if(kb) print(k)
      
      col <- ifelse(kb, 2, 1)
      col <- col2rgb(col)/255
      col <- rgb(col[1],col[2],col[3],alpha=ifelse(kb, 1, 0.1))
      
      lines(1:5, j[k,], col=col)
    }  
  }
}

# total
pal_prim <- brewer.pal(9, 'Set1') 

for(i in 1:2){
  ind <- seq(i, 18, by=2)
  
  df <- as.data.frame(t(ra_tot[,ind]))
  
  plot(NA, xlim=c(1,ncol(df)), xaxt='n', xlab='', 
       ylim=range(df, na.rm=T), ylab=ifelse(i == 1, 'nb seq', '%nb_seq'),
       log=ifelse(i == 1, '', 'y'), main=paste(c('raw','relabu')[i], 'total dataset'))
  axis(1, at=1:ncol(df), labels=names(df))
  
  for(j in 1:nrow(df)){
    lines(1:5, df[j,], col=pal_prim[j])
  }
  
  legend('bottomleft', legend=prim_names, pch=19, col=pal_prim, bty='n')
}

dev.off()



#####




















