#####
# climarctic preparation
#####

rm(list=ls())
gc()

require(RColorBrewer) # brewer.pal 

# dir loads ####
dir_in <- 'Projets/Climarctic/stats/in/'
dir_out <- 'Projets/Climarctic/stats/out_janv/'
dir_save <- paste0(dir_out, 'saves/') 
dir.create(dir_save, showWarnings=F, recursive=T)

#---
mr_tot  <- read.table(paste0(dir_in, 'from_cluster_190115/02_clust.mr'), h=T)
ass_tot <- read.table(paste0(dir_in, 'from_cluster_190115/02_clust.DB.wang.taxonomy'), row.names=1)
fa_tot  <- read.table(paste0(dir_in, 'from_cluster_190115/02_clust.fa'))
env_tot <- read.table(paste0(dir_in, 'env_clim_grad.csv'), row.names=1)

# reorganize mr
n_mr <- sapply(strsplit(names(mr_tot), '_'), '[', 1:2)
ord <- order(as.numeric(n_mr[2,]))

mr_tot <- mr_tot[,ord]
names(mr_tot) <- sapply(lapply(strsplit(names(mr_tot), '_'), '[', 1:2), function(x) paste(x, collapse='_'))

# reorganize fa
n_fa <- as.character(fa_tot[seq(1,nrow(fa_tot), by=2),])
n_fa <- substr(n_fa, 2, nchar(n_fa))

fa_tot <- fa_tot[seq(2,nrow(fa_tot), by=2),]
names(fa_tot) <- n_fa

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

# reorganize ass
ass_tot <- data.frame(taxo=ass_tot$V2, seq=fa_tot)
ass_tot$taxo <- gsub('[[:punct:]][[:digit:]]{2,3}[[:punct:]]{2}|;', '|', ass_tot$taxo)

taxo_tot <- matrix(unlist(strsplit(as.character(ass_tot$taxo), '|', fixed=T)), ncol=7, byrow=T)[ord,]
dimnames(taxo_tot) <- list(row.names(ass_tot), c('reign','phylum','division','order','family','genus','species'))

taxo_tot <- as.data.frame(t(apply(taxo_tot, 1, function(x) {
  ind_unc <- grep('^[[:lower:]]', x)
  for(i in ind_unc){
    x[i] <- ifelse(grepl('_X', x[i-1]), paste0(x[i-1], 'X'), paste0(x[i-1], '_X'))
  }
  return(x)
})))

# sort taxo
ord_taxo <- order(ass_tot$taxo)

ass_sort <- ass_tot[ord,]
mr_sort <- mr_tot[,ord]
taxo_sort <- taxo_tot[ord,]

# palette
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

# permu
permu <- 10000

#---
file <- paste0(dir_save, 'lst_comm.Rdata')
save(mr_sort, env_sort, ass_sort, taxo_sort, lst_palev, file=file)

#










