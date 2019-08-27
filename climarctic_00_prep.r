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

env_tot <- env_tot[,-grep(paste(c('empty','fresh','dry','burn'), collapse='|'), names(env_tot))]

env_tot$site <- gl(2, 54, labels=c('Knudsenheia','Ossian'))
env_tot$moisture <- gl(3,18,108, c('dry','medium','wet'))
env_tot$plot <- as.factor(env_tot$plot)
env_tot$depth <- rep(gl(2,3, labels=c('top','deep')), 18)

env_tot$moist_in_site <- factor(apply(env_tot[,c('moisture','site')], 1, function(x) paste(x, collapse='_')))
env_tot$plot_in_moist_in_site <- factor(apply(env_tot[,c('plot','moist_in_site')], 1, function(x) paste(x, collapse='_')))
env_tot$quad_in_plot_in_moist <- factor(apply(env_tot[,c('quadrat','plot_in_moist_in_site')], 1, function(x) paste(x, collapse='_')))

# palette ----
lev_site <- levels(env_tot$site)
lev_mois <- levels(env_tot$moisture)
lev_dept <- levels(env_tot$depth)

lst_palev <- list(site    =list(pal=brewer.pal(length(lev_site), 'Dark2')[seq_along(lev_site)],
                                lev=lev_site),
                  moisture=list(pal=brewer.pal(length(lev_mois), 'Set1'),
                                lev=lev_mois),
                  depth   =list(pal=brewer.pal(length(lev_dept), 'Set2')[seq_along(lev_dept)],
                                lev=lev_dept)
)

# permu ----
permu <- 10000

# loop the primers ####
prim_names <- c('16S_V1-3','18S_V4','pmoA_mb661','pmoA_A682','ITS2','phoD','nifH')
ind_prim <- c(1,2)

lst_comm <- NULL
for(i in ind_prim) {
  
  print(i)
  
  id_plate <- as.character(i)
  if(i < 10){
    id_plate <- paste0('0', id_plate)
  }
  
  #---
  file <- paste0(dir_save, 'ini_', id_plate, '.Rdata')
  if(file.exists(file)){
    load(file)
  } else {
    mr_ini  <- read.table(paste0(dir_in, 'from_cluster/', id_plate, '/', id_plate, '_clust.mr'), h=T)
    ass_ini <- read.table(paste0(dir_in, 'from_cluster/', id_plate, '/', id_plate, '_clust.DB.wang.taxonomy'), row.names=1)
    fa_ini  <- read.table(paste0(dir_in, 'from_cluster/', id_plate, '/', id_plate, '_clust.fa'))
    save(mr_ini, ass_ini, fa_ini, file=file)
  }
  
  # reorganize mr
  n_mr <- sapply(strsplit(names(mr_ini), '_'), '[', 1:2)
  ord <- order(as.numeric(n_mr[2,]))
  
  mr_tot <- mr_ini[,ord]
  names(mr_tot) <- sapply(lapply(strsplit(names(mr_tot), '_'), '[', 1:2), function(x) paste(x, collapse='_'))
  
  # reorganize fa
  n_fa <- as.character(fa_ini[seq(1,nrow(fa_ini), by=2),])
  n_fa <- substr(n_fa, 2, nchar(n_fa))
  
  fa_tot <- fa_ini[seq(2,nrow(fa_ini), by=2),]
  names(fa_tot) <- n_fa
  
  # reorganize ass
  ass_tot <- data.frame(taxo=ass_ini$V2, seq=fa_tot)
  
  ass_tot$taxo <- gsub('[[:punct:]][[:digit:]]{2,3}[[:punct:]]{2}|;', '|', ass_tot$taxo)
  
  taxo_tot <- strsplit(as.character(ass_tot$taxo), '|' , fixed=T)
  nb_lev <- length(taxo_tot[[1]])
  
  taxo_tot <- matrix(unlist(taxo_tot), ncol=nb_lev, byrow=T)
  row.names(taxo_tot) <- row.names(ass_tot)
  
  taxo_tot <- as.data.frame(t(apply(taxo_tot, 1, function(x) {
    x <- gsub('Incertae_Sedis', 'X', x)
    x <- gsub('Unknown', 'unknown', x)
    
    ind_unc <- grep('^[[:lower:]]', x)
    for(i in ind_unc){
      x[i] <- ifelse(grepl('_X', x[i-1]), paste0(x[i-1], 'X'), paste0(x[i-1], '_X'))
    }
    return(x)
  })))
  
  # sort taxo
  ord_taxo <- order(ass_tot$taxo)
  
  ass_sort <- ass_tot[ord_taxo,]
  mr_sort <- mr_tot[,ord_taxo]
  taxo_sort <- taxo_tot[ord_taxo,]
  env_sort <- env_tot[row.names(mr_sort),]
  
  # rrarefy ----
  rs <- rowSums(mr_sort)
  thresh <- seq(min(rs), max(rs), length.out=20)
  thresh <- thresh[-length(thresh)]
  
  # nb seq
  plot(sort(rs))
  abline(h=thresh)
  
  # perc lost
  set.seed(0)
  perc_lost <- foreach(j=thresh, .combine=cbind, .verbose=T) %dopar% {
    if(nrow(mr_sort[rs >= j,]) > 1){
      mr_rrf <- rrarefy(mr_sort[rs >= j,], j)
      mr_rrf <- mr_rrf[,colSums(mr_rrf) != 0]

      return(c(sum(mr_rrf) / sum(mr_sort), ncol(mr_rrf) / ncol(mr_sort), nrow(mr_rrf) / nrow(mr_sort)))
    }
  }

  dimnames(perc_lost) <- list(c('seq','otu','smp'), paste0('t', thresh[1:ncol(perc_lost)]))
  
  l <- NULL
  for(j in c(0.5, 0.75, 1)){
    optimum <- ceiling(mean(apply(perc_lost[-3,], 1, function(x) thresh[which(x == max(x, na.rm=T))])))*j
    
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
  pdf(paste0(dir_prep, 'rraref_optimum_', prim_names[i], '.pdf'))
  
  plot(NA, xlim=range(thresh), ylim=c(0,1), xlab='thresh', ylab='percentage',
       main=paste('seq ini:', sum(mr_sort), 'otu ini:', ncol(mr_sort), 'smp ini:', nrow(mr_sort)))
  
  for(j in 1:nrow(perc_lost)){
    lines(thresh[1:ncol(perc_lost)], perc_lost[j,], col=j)
  }
  
  abline(v=c(optimum, optimum*0.75, optimum*0.5), lty=3)

  legend('topright',legend=l)
  
  dev.off()
  
  #---
  set.seed(0)
  mr_rrf <- rrarefy(mr_sort[rs >= optimum,], optimum)
  mr_rrf <- as.data.frame(mr_rrf[,colSums(mr_rrf) != 0])
  
  env_rrf <- env_tot[row.names(mr_rrf),]
  
  ass_rrf <- ass_sort[names(mr_rrf),]
  taxo_rrf <- taxo_sort[names(mr_rrf),]
  
  lst_comm[[prim_names[i]]] <- list(raw=list(env=env_sort, mr=mr_sort, ass=ass_sort, taxo=taxo_sort),
                                    rrf=list(env=env_rrf,  mr=mr_rrf,  ass=ass_rrf,  taxo=taxo_rrf))
  
}

#---
file <- paste0(dir_save, 'lst_comm.Rdata')
save(lst_comm, env_tot, lst_palev, permu, file=file)

#










