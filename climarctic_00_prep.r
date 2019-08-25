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
env_tot <- read.table(paste0(dir_in, 'env_clim_grad.csv'), row.names=1)

# reorganize env
rn <- row.names(env_tot)
row.names(env_tot) <- ifelse(as.numeric(rn) < 100, ifelse(as.numeric(rn) < 10, paste0('T00', rn), paste0('T0', rn)), paste0('T', rn))

env_sort <- env_tot[,c(2:5,7)]
names(env_sort) <- c('site','moisture','plot','depth','quadrat')

env_sort$site <- gl(2, 54, labels=c('Knudsenheia','Ossian'))
env_sort$moisture <- gl(3,18,108, c('dry','medium','wet'))

env_sort$plot <- as.factor(env_sort$plot)

env_sort$moist_in_site <- factor(apply(env_sort[,c('moisture','site')], 1, function(x) paste(x, collapse='_')))

env_sort$plot_in_moist_in_site <- factor(apply(env_sort[,c('plot','moist_in_site')], 1, function(x) paste(x, collapse='_')))

env_sort$quad_in_plot_in_moist <- factor(apply(env_sort[,c('quadrat','plot_in_moist_in_site')], 1, function(x) paste(x, collapse='_')))

# palette ----
lev_site <- levels(env_sort$site)
lev_mois <- levels(env_sort$moisture)
lev_dept <- levels(env_sort$depth)

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

lst_comm <- foreach(i=ind_prim, .verbose=T) %dopar% {
  
  id_plate <- as.character(i)
  if(i < 10){
    id_plate <- paste0('0', id_plate)
  }
  
  #---
  mr_ini  <- read.table(paste0(dir_in, 'from_cluster/', id_plate, '/', id_plate, '_clust.mr'), h=T)
  ass_ini <- read.table(paste0(dir_in, 'from_cluster/', id_plate, '/', id_plate, '_clust.DB.wang.taxonomy'), row.names=1)
  fa_ini  <- read.table(paste0(dir_in, 'from_cluster/', id_plate, '/', id_plate, '_clust.fa'))
  
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
  
  # rrarefy ----
  rs <- rowSums(mr_sort)
  thresh <- seq(min(rs), max(rs), length.out=20)
  thresh <- thresh[-length(thresh)]
  
  # nb seq
  plot(sort(rs))
  abline(h=thresh)
  
  # perc lost
  perc_lost <- matrix(NA, nrow=3, ncol=length(thresh), dimnames=list(c('seq','otu','smp'), paste0('t', thresh)))
  
  set.seed(0)
  for(j in thresh){
    if(nrow(mr_sort[rs >= j,]) > 1){
      mr_rrf <- rrarefy(mr_sort[rs >= j,], j)
      mr_rrf <- mr_rrf[,colSums(mr_rrf) != 0]
      
      n <- paste0('t', j)
      perc_lost[1,n] <- sum(mr_rrf) / sum(mr_sort)
      perc_lost[2,n] <- ncol(mr_rrf) / ncol(mr_sort)
      perc_lost[3,n] <- nrow(mr_rrf) / nrow(mr_sort)
    }
  }
  
  l <- NULL
  for(j in c(0.5, 0.75, 1)){
    optimum <- ceiling(mean(apply(perc_lost[-3,], 1, function(x) thresh[which(x == max(x, na.rm=T))])))*j
    
    # which smps did we loose
    set.seed(0)
    mr_rrf <- rrarefy(mr_sort[rs >= optimum,], optimum)
    mr_rrf <- mr_rrf[,colSums(mr_rrf) != 0]
    
    e <- env_sort[row.names(mr_rrf),]
    
    ll <- NULL
    for(k in c('site','moisture','depth')){
      ll[[k]] <- paste(c(j, k, signif(table(e[[k]]) / table(env_sort[[k]]), 2)), collapse=' ')
    }
    
    l <- cbind(l, ll)
  }
  
  #---
  pdf(paste0(dir_prep, 'rraref_optimum_', prim_names[i], '.pdf'))
  
  plot(NA, xlim=range(thresh), ylim=c(0,1), xlab='thresh', ylab='percentage',
       main=paste('seq ini:', sum(mr_sort), 'otu ini:', ncol(mr_sort), 'smp ini:', nrow(mr_sort)))
  
  for(j in 1:nrow(perc_lost)){
    lines(thresh, perc_lost[j,], col=j)
  }
  
  abline(v=c(optimum, optimum*0.75, optimum*0.5), lty=3)

  legend('topright',legend=l)
  
  dev.off()
  
  #---
  set.seed(0)
  mr_rrf <- rrarefy(mr_sort[rs >= optimum,], optimum)
  mr_rrf <- mr_rrf[,colSums(mr_rrf) != 0]
  
  env_rrf <- env_sort[row.names(mr_rrf),]
  
  ass_rrf <- ass_sort[names(mr_rrf),]
  taxo_rrf <- taxo_sort[names(mr_rrf),]
  
  return(list(mr_sort=mr_sort, ass_sort=ass_sort, taxo_sort=taxo_sort, mr_rrf=mr_rrf, ass_rrf=ass_rrf, taxo_rrf=taxo_rrf, env_rrf=env_rrf))
}

names(lst_comm) <- prim_names[ind_prim]

#---
file <- paste0(dir_save, 'lst_comm.Rdata')
save(lst_comm, env_sort, lst_palev, permu, file=file)

#










