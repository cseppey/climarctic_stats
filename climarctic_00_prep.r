#####
# climarctic preparation
#####

rm(list=ls())

require(foreach)
require(doSNOW) # e.g. makeSOCKcluster()
require(RColorBrewer) # brewer.pal 
require(compositions)
require(abind)

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
env_ini_fact <- read.table(paste0(dir_in, 'env_clim_grad.csv'), row.names=1, h=T)
env_ini_chim <- read.table(paste0(dir_in, 'env_clim_chim.csv'), h=T)

# env fact %%%
rn <- row.names(env_ini_fact)
row.names(env_ini_fact) <- ifelse(as.numeric(rn) < 100, ifelse(as.numeric(rn) < 10, paste0('T00', rn), paste0('T0', rn)), paste0('T', rn))

# env_chim %%%
as_char <- as.character(env_ini_chim$sample)
env_ini_chim$site  <- substr(as_char, 1, 1)
env_ini_chim$moist <- rep(gl(3,3, labels=c('dry','intermediate','wet')), 4)
env_ini_chim$depth <- gl(2, 18, labels=c('top','deep'))

env_ini_chim <- env_ini_chim[c(matrix(1:36, nrow=2, byrow=T)),]
env_ini_chim <- env_ini_chim[as.numeric(gl(36,3)),]

# take the interesting variables
env_tot <- cbind.data.frame(env_ini_fact[,c(2:5,7:15)], env_ini_chim[,c('sand','silt','clay')])
names(env_tot)[1:9] <- c('site','moisture','plot','depth','quadrat','empty','fresh','dry','burn')

env_tot[,c('fresh','dry','burn')] <- sapply(env_tot[,c('fresh','dry','burn')], function(x) x-env_tot$empty)
env_tot$rh <- (env_tot$fresh-env_tot$dry) / env_tot$fresh
env_tot$om <- (env_tot$dry-env_tot$burn) / env_tot$dry
env_tot$C_N <- env_tot$C / env_tot$N

env_tot <- env_tot[,-grep(paste(c('empty','fresh','dry','burn'), collapse='|'), names(env_tot))]

env_tot$site <- gl(2, 54, labels=c('Knudsenheia','Ossian'))
env_tot$moisture <- ordered(gl(3,18,108, c('dry','intermediate','wet')))
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
prim_names <- c('01_16S_bact','02_18S_euk','03_pmoA_mb661','04_pmoA_A682','05_ITS_fun','06_phoD','07_nifH', '08_16S_cyano','09_nirS')
ind_prim <- c(1,2,5,8)

# loop the primers ####
# pdf(paste0(dir_prep, '00_distro.pdf'))
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
    ass_ini <- read.table(paste0(dir_in, 'from_cluster/', id_plate, '/', id_plate, '_clust.ass'), row.names=1)
    fa_ini  <- read.table(paste0(dir_in, 'from_cluster/', id_plate, '/', id_plate, '_clust.fa'))
    save(mr_ini, ass_ini, fa_ini, file=file) 
  }
  
  # cleaning ----
  # reorganize mr ---
  mr_tot <- mr_ini[grep('T|B', row.names(mr_ini)),]
  
  # remove the OTU found in the blanks
  ind_blk <- grepl('B', row.names(mr_tot))
  rs <- rowSums(mr_tot)
  ord <- order(rs)
  plot(rs[ord], col=as.numeric(ind_blk[ord])+1, pch=19, main=prim_names[i])
  
  if(length(which(ind_blk))){
    mr_blk <- ifelse(mr_tot[grep('B', row.names(mr_tot)),] == 0, 0, 1)
    ind_conta <- colSums(mr_blk) != 0
  } else {
    ind_conta <- rep(F, nrow(mr_tot))
  }
  
  mr_tot <- mr_tot[grep('T', row.names(mr_tot)),ind_conta == F]
  mr_tot <- mr_tot[rowSums(mr_tot) != 0,colSums(mr_tot) != 0]
  
  print(c(sum(mr_tot), ncol(mr_tot)))

  # find low sequence samples
  lrs <- log(sort(rowSums(mr_tot)))
  
  x <- brks <- 1:length(lrs)
  
  mse <- NULL
  for(j in 1:length(brks)){
    mod <- lm(lrs~x*(x <= brks[j]) + x*(x < brks[j]))
    mse <- c(mse, summary(mod)$sigma)
  }
  
  min_mse <- which(mse == min(mse))
  
  mod <- lm(lrs ~ x*(x < min_mse) + x*(x > min_mse))
  co <- coef(mod)
  
  #---
  cairo_ps(paste0(dir_prep, 'piecewise_', prim_names[i],'.eps'))
  
  plot(lrs~x, xlab='', ylab='log of rowSums', main=prim_names[i])
  
  curve(co[1]+co[3] + (co[2]+co[5])*x, add=T, from=1, to=min_mse+2)
  curve(co[1]+co[4] + (co[2])*x, add=T, from=min_mse+2, to=length(lrs))
  abline(v=min_mse+2, lty=3)  
  
  low_seq <- names(lrs)[1:which(mse == min(mse))]
  
  dev.off()
  
  # reorganize fa ---
  n_fa <- as.character(fa_ini[seq(1,nrow(fa_ini), by=2),])
  n_fa <- substr(n_fa, 2, nchar(n_fa))

  fa_tot <- fa_ini[seq(2,nrow(fa_ini), by=2),]
  names(fa_tot) <- n_fa
  
  fa_tot <- fa_tot[names(mr_tot)]

  # reorganize ass ---
  tax <- ass_ini$V2
  names(tax) <- row.names(ass_ini)
  tax <- tax[names(mr_tot)]
  
  ass_tot <- data.frame(taxo=tax, seq=fa_tot)
  
  ass_tot$taxo <- gsub('[[:punct:]][[:digit:]]{2,3}[[:punct:]]{2}|;', '|', ass_tot$taxo)

  # correct the taxon that are found at many tax lev
  ass_tot$taxo <- gsub('Lineage_IIa|Elusimicrobia', 'Elusimicrobia|Lineage_IIa', ass_tot$taxo, fixed=T)
  ass_tot$taxo <- gsub('Candidatus_Magasanikbacteria|Parcubacteria', 'Candidatus_Magasanikbacteria|Candidatus_Magasanikbacteria_unclassified',
                       ass_tot$taxo, fixed=T)
  ass_tot$taxo <- gsub('TM7|Candidatus_Saccharibacteria', 'TM7|u', ass_tot$taxo, fixed=T)
  ass_tot$taxo <- gsub('Elusimicrobia|Lineage_IIa', 'Elusimicrobia|u', ass_tot$taxo, fixed=T)
  
  # correct taxon with different parent taxon
  ass_tot <- ass_tot[names(mr_tot),]

  #---
  taxo_tot <- strsplit(as.character(ass_tot$taxo), '|' , fixed=T)
  nb_lev <- length(taxo_tot[[1]])

  taxo_tot <- matrix(unlist(taxo_tot), ncol=nb_lev, byrow=T)
  row.names(taxo_tot) <- row.names(ass_tot)

  taxo_tot <- as.data.frame(t(apply(taxo_tot, 1, function(x) {
    
    # scrap taxonomy correction
    ind_scrap <- which(x %in% 'Geobacter_sp.'
                       | x %in% 'Acidobacteria'
                       | x %in% 'Acidobacteria_bacterium'
                       | x %in% 'Rhodospirillales_bacterium'
                       | x %in% 'Actinobacteria_bacterium'
                       | x %in% 'Sterolibacterium'
                       | x %in% 'Sphingosinicella_sp.'
                       | x %in% 'Sphingomonas_sp.'
                       | x %in% 'Rhodoferax_sp.'
                       | x %in% 'Rhodococcus_sp.'
                       | x %in% 'Rhodobacter_sp.'
                       | x %in% 'Pseudanabaena'
                       | x %in% 'Nocardioides_sp.'
                       | x %in% 'Microbacterium_sp.'
                       | x %in% 'Mesorhizobium_sp.'
                       | x %in% 'Leptothrix_sp.'
                       | x %in% 'Kribbella_sp.'
                       | x %in% 'Devosia_sp.'
                       | x %in% 'Cellulomonas_sp.'
                       | x %in% 'Caulobacter_sp.'
                       | x %in% 'Caenimonas_sp.'
                       | x %in% 'Betaproteobacteria_bacterium'
                       | x %in% 'Bacteroidetes_bacterium'
                       | x %in% 'Afipia_sp.'
                       | x %in% 'Afipia_genosp.'
                       | x %in% 'Xanthomonadaceae_bacterium'
                       | x %in% 'Sphingobacteriaceae_bacterium'
                       | x %in% 'Rhodospirillaceae_bacterium'
                       | x %in% 'Monodopsis_sp.'
                       | x %in% 'Cytophagaceae_bacterium'
                       | x %in% 'Caulobacteraceae_bacterium'
                       | x %in% 'Candidatus'
                       | x %in% 'Sphingobacteriales_bacterium'
                       | x %in% 'Bdellovibrionales_bacterium'
                       | x %in% 'Burkholderiales_bacterium'
                       | x %in% 'Hypsibius_dujardini'
                       | x %in% 'Nostoc'
                       | x %in% 'Verrucomicrobia'
                       | x %in% 'Synechococcus'
                       | x %in% 'Sphingobacterium'
                       | x %in% 'Sorangiineae'
                       | x %in% 'Nitrospirae'
                       | x %in% 'Armatimonadetes_bacterium'
                       | x %in% 'Antarctic'
                       | x %in% 'Clostridiales_bacterium'
    )[1]
    if(is.na(ind_scrap) == F) x[ind_scrap:length(x)] <- 'u'
    
    x <- gsub('.*Incertae_.*|Unknown|.*-[pcofgs]$|.*_unclassified', 'u', x)
    x <- gsub('_Incertae|_clade|_lineage|_terrestrial|_marine|_genosp.|_sensu|_sp.|_bacterium','', x)

    for(j in 2:length(x)){
      if(x[j] %in% x[1:(j-1)]){
        x[j] <- 'u'
      }
    }
    
    ind_unc <- grep('^[[:lower:]]', x)
    if(length(ind_unc)){
      if(ind_unc[1] == 1){
        x[1] <- 'Life'
        ind_unc <- ind_unc[-1]
      }
    }
    for(j in ind_unc){
      x[j] <- ifelse(grepl('_X', x[j-1]), paste0(x[j-1], 'X'), paste0(x[j-1], '_X'))
    }
    return(x)
  })))

  # sort taxo
  ord_taxo <- order(ass_tot$taxo)

  ass_sort <- ass_tot[ord_taxo,]
  mr_sort <- mr_tot[,ord_taxo]
  mr_sort <- mr_sort[rowSums(mr_sort) != 0,]
  taxo_sort <- taxo_tot[ord_taxo,]
  
  # remove false positive
  if(i != 1 & i != 2 & i != 5 & i != 8){
    ind_true <- taxo_sort$V1 == 'TRUE'
    
    ass_sort <- droplevels(ass_sort[ind_true,])
    mr_sort <- mr_sort[,ind_true]
    taxo_sort <- droplevels(taxo_sort[ind_true,-1])
  }
  
  # taxo_clean
  taxo_false <- switch(i,
                       '1' = 'Chloroplast|Mitochondria',
                       '2' = 'Metazoa|Embryophyceae',
                       '3' = 'Life|TRUE_X',
                       '4' = 'Life|TRUE_X',
                       '5' = 'Life|Plantae',
                       '6' = 'Life|TRUE_X',
                       '7' = 'Life|TRUE_X',
                       '8' = 'Life|TRUE_X|Bacteria_X|Chloroplast|Acidobacteriota|Actinobacteriota|Chloroflexi|Firmicutes|Methylomirabilota|Nitrospinota|Patescibacteria|Planctomycetota|Verrucomicrobiota|WS4',
                       '9' = 'Life|TRUE_X')
  ind_tf <- grepl(taxo_false, taxo_sort[,1]) | grepl(taxo_false, taxo_sort[,2]) | grepl(taxo_false, ass_sort$taxo)
  
  if(length(which(ind_tf))){
    ass_sort <- droplevels(ass_sort[ind_tf == F,])
    mr_sort <- mr_sort[,ind_tf == F]
    mr_sort <- mr_sort[rowSums(mr_sort) != 0,]
    taxo_sort <- droplevels(taxo_sort[ind_tf == F,])
  }
  
  env_sort <- env_tot[row.names(mr_sort),]
  env_sort <- data.frame(env_sort, low_seq=row.names(env_sort) %in% low_seq)
  
  # communities compositional normalisation ----
  mr_clr <- clr(mr_sort)
  mr_clr <- as.data.frame(matrix(c(mr_clr), nrow=nrow(mr_clr), dimnames=dimnames(mr_clr)))
  
  env_clr <- env_tot[row.names(mr_clr),]
  env_clr <- data.frame(env_clr, low_seq=row.names(env_clr) %in% low_seq)
  
  ass_clr <- ass_sort[names(mr_clr),]
  taxo_clr <- taxo_sort[names(mr_clr),]

  # wiehout the low_sequences samples (piecewise regression)
  mr_nls <- mr_sort[env_sort$low_seq == F,]
  mr_nls <- mr_nls[,colSums(mr_nls) != 0]
  
  mr_clr2 <- clr(mr_nls)
  mr_clr2 <- as.data.frame(matrix(c(mr_clr2), nrow=nrow(mr_clr2), dimnames=dimnames(mr_clr2)))
  
  env_clr2 <- env_tot[row.names(mr_clr2),]
  
  ass_clr2 <- ass_sort[names(mr_clr2),]
  taxo_clr2 <- taxo_sort[names(mr_clr2),]
  
  lst_comm[[prim_names[i]]] <- list(raw = list(env=env_sort, mr=mr_sort, ass=ass_sort, taxo=taxo_sort),
                                    clr = list(env=env_clr , mr=mr_clr , ass=ass_clr , taxo=taxo_clr),
                                    clr2= list(env=env_clr2, mr=mr_clr2, ass=ass_clr2, taxo=taxo_clr2))

}

# dev.off()

#---
file <- paste0(dir_save, '00_lst_comm.Rdata')
save(lst_comm, env_tot, lst_palev, permu, file=file)

file <- paste0(dir_save, '00_env_tot.Rdata')
save(env_tot, file=file)

stopCluster(cl)

# check length distro ####

lst_lim <- list(`16S_bact`=c(450,570),
                `18S_euk`=c(380,480),
                pmoA_mb661=c(490,530),
                pmoA_A682=c(510,550),
                ITS_fun=c(270,500),
                phoD=c(330,450),
                nifH=c(310,450),
                `16S_cyano`=c(390,470),
                nirS=c(380,460))

pdf(paste0(dir_prep, 'lgt_distro.pdf'), width=15, height=15)
par(mfrow=c(3,3))

for(i in c(1,8,2,5,3,4,6,7,9)){
  files <- list.files(paste0('Projets/Climarctic/bioinfo/archive/191031/0', i, '/filter_test/distro_test'), full.names=T)

  lgt_dis <- NULL
  for(j in files){
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
for(i in c(1,8,2,5,3,4,6,7,9)){
  out_bf <- read.table(paste0('Projets/Climarctic/bioinfo/archive/191229/filter_test/cnt_output/cnt_output', i),
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
    axis(1, at=1:ncol(j), labels=names(j), las=2)

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
  axis(1, at=1:ncol(df), labels=names(df), las=2)

  for(j in 1:nrow(df)){
    lines(1:5, df[j,], col=pal_prim[j])
  }

  if(i == 2)legend('bottomleft', legend=prim_names, pch=19, col=pal_prim, bty='n', cex=0.5)
}

dev.off()



#####




















