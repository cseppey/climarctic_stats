#####
# climarctic preparation
#####

print('##### Climarctic 00 data prep #####')

rm(list=ls())

require(foreach)
require(doSNOW) # e.g. makeSOCKcluster()
require(RColorBrewer) # brewer.pal 
require(zCompositions)

source('~/bin/src/my_prog/R/rarecurveMPI.r')

# prep cluster
cl <- makeSOCKcluster(4)

clusterEvalQ(cl, library(zCompositions))

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

env_tot$MiS <- factor(apply(env_tot[,c('moisture','site')], 1, function(x) paste(x, collapse='_')))
env_tot$PiMiS <- factor(apply(env_tot[,c('plot','MiS')], 1, function(x) paste(x, collapse='_')))
env_tot$QiPiMiS <- factor(apply(env_tot[,c('quadrat','PiMiS')], 1, function(x) paste(x, collapse='_')))
env_tot$combi <- factor(paste(env_tot$MiS, env_tot$depth, sep='_'))

env_tot$pH[env_tot$pH > 8.5] <- NA

# palette ----
lev_site <- levels(env_tot$site)
lev_mois <- levels(env_tot$moisture)
lev_dept <- levels(env_tot$depth)

pal_dark2 <- brewer.pal(8, 'Dark2')
lst_palev <- list(site    =pal_dark2[1:2],# lev=lev_site,
                  moisture=pal_dark2[3:5],# lev=lev_mois,
                  depth   =pal_dark2[6:7]#, lev=lev_dept
)
names(lst_palev$site) <- lev_site
names(lst_palev$moisture) <- lev_mois
names(lst_palev$depth) <- lev_dept

#---
fact_3 <- c('site','moisture','depth')

# legend function
leg <- function(x = NULL){
  plot.new()
  
  text(0.5,0.85,x)
  
  # factors colors
  legend(0.5,0.5, legend=sapply(strsplit(names(unlist(lst_palev)), '.', fixed=T), '[[', 2),
         bty='n', xjust=0.5, yjust=0.5, pch=19, col=unlist(lst_palev))
} 

# layout function
lay <- function() {
  layout(matrix(c(1,3,4, 1,5,6, 0,2,2), nrow=3, byrow=T), width=c(0.2,1,1), height=c(1,1,0.3))
  par(mar=rep(0,4), oma=c(1,1,3,3))
}

# PCA ----
for(i in c('top|deep', 'top', 'deep')){
  env_sc <- env_tot[grep(i, env_tot$depth),-which(names(env_tot) == 'sand')]
  is.num <- sapply(env_sc, is.numeric)
  env_sc[,is.num] <- scale(env_sc[,is.num])
  env_sc <- na.omit(env_sc)
  
  pca <- rda(env_sc[,is.num])
  
  # plot
  cairo_ps(paste0(dir_prep, 'PCA_', i, '.eps'), width=10, height=7)
  lay()
  
  for(j in 2:1){
    plot.new()
    text(0.5,0.5, labels=paste('PCA',j), srt=ifelse(j == 1, 0, 90))
  }
  
  # ordinations
  smp <- pca$CA$u[,1:2]
  var <- pca$CA$v[,1:2]
  
  for(j in seq_along(fact_3)){
    
    fact <- fact_3[j]
    
    pal <- lst_palev[[fact]][env_sc[[fact]]]
    
    #---
    plot(smp, xlim=range(smp[,1]), ylim=range(smp[,2]), xaxt='n', yaxt='n',
         col=pal, pch=NA)
    
    abline(v=0, h=0, lty=3)
  
    # variables
    rng_smp <- sapply(as.data.frame(smp), range, simplify='matrix')
    rng_var <- sapply(as.data.frame(var), range, simplify='matrix')
    
    divis <- max(rng_var/rng_smp)
    v2 <- var/divis
    
    arrows(0,0,v2[,1],v2[,2], length=0, lty=2)
    
    if(j < 3){       
      max_smp <- max(abs(rng_var[,1]/divis))
      max_var <- max(abs(rng_var[,1]))
      axis(3, at=seq(-max_smp, max_smp, length.out=9), labels=round(seq(-max_var, max_var, length.out=9), 2), col=2)
    }
    
    if(j > 1){
      max_smp <- max(abs(rng_var[,2]/divis))
      max_var <- max(abs(rng_var[,2]))
      axis(4, at=seq(-max_smp, max_smp, length.out=9), labels=round(seq(-max_var, max_var, length.out=9), 2), col=2)
    }
    
    text(v2[,1:2], labels=paste(row.names(v2)))
    
    # samples
    if(j > 1){
      axis(1)
    }
    
    if(j %% 2 == 1){
      axis(2)
    }
    
    ordispider(smp, env_sc$combi, col='grey80')
    
    points(smp[,1], smp[,2], pch=19, col=pal)
    # text(site, labels=row.names(site), col=pal)
    
  }
  
  # legend ---
  leg()
  
  dev.off()
}


# permu and primer names ----
permu <- 10000

#---
prim_names <- c('01_16S_bact','02_18S_euk','03_pmoA_mb661','04_pmoA_A682','05_ITS_fun','06_phoD','07_nifH', '08_16S_cyano','09_nirS')
ind_prim <- c(1,2,5,8)
n_comm <- prim_names[ind_prim]

# loop the primers ####
lst_comm <- NULL
for(i in ind_prim) {
  
  # initialize ----
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
    fa_ini  <- read.table(paste0(dir_in, 'from_cluster/', id_plate, '/', id_plate, '_clust.fa'))
    lst_ass_ini <- lapply(list.files(paste0(dir_in, 'from_cluster/', id_plate), '.ass', full.names=T), read.table)
    save(mr_ini, lst_ass_ini, fa_ini, file=file) 
  }
  
  # cleaning ####
  # reorganize mr ----
  mr_tot <- mr_ini[grep('T|B', row.names(mr_ini)),]
  
  # removal of faild cyaB sample
  if(i == 8){
    ech_scrap <- c('T019','T029','T031','T033','T097','T103','T105','T004','T006','T010','T012','T016','T018',
                   'T020','T028','T022','T030','T034','T036','T040','T046','T042','T048','T056','T062','T068',
                   'T086','T098','T005','T007','T011','T079','T023','T032','T078','T052','T072','T038','T045')
    ind_scrap <- grepl(paste(ech_scrap, collapse='|'), row.names(mr_tot))
    mr_tot <- mr_tot[ind_scrap,]
  }
  
  # remove the OTU found in the blanks
  ind_blk <- grepl('B', row.names(mr_tot))

  print(c('smps', nrow(mr_tot[ind_blk == F,]), 'OTUs', ncol(mr_tot[ind_blk == F,colSums(mr_tot[ind_blk == F,]) != 0]),
          'reads', sum(mr_tot[ind_blk == F,])))
  
  rs <- rowSums(mr_tot)
  ord <- order(rs)
  plot(rs[ord], col=as.numeric(ind_blk[ord])+1, pch=19, main=prim_names[i])
  
  if(length(which(ind_blk))){
    mr_blk <- ifelse(mr_tot[grep('B', row.names(mr_tot)),] == 0, 0, 1)
    ind_conta <- colSums(mr_blk) != 0
  } else {
    ind_conta <- rep(F, ncol(mr_tot))
  }
  
  mr_tot <- mr_tot[grep('T', row.names(mr_tot)),ind_conta == F]
  mr_tot <- mr_tot[rowSums(mr_tot) != 0,colSums(mr_tot) != 0]
  
  # reorganize ass ----
  # reorganize fa
  n_fa <- as.character(fa_ini[seq(1,nrow(fa_ini), by=2),])
  n_fa <- substr(n_fa, 2, nchar(n_fa))

  fa_tot <- fa_ini[seq(2,nrow(fa_ini), by=2),]
  names(fa_tot) <- n_fa
  
  fa_tot <- fa_tot[names(mr_tot)]

  # reorganize ass
  lst_ass_tot <- NULL
  for(j in seq_along(lst_ass_ini)){
    ass_ini <- lst_ass_ini[[j]]
    
    tax <- ass_ini$V2
    names(tax) <- ass_ini$V1
    tax <- tax[names(mr_tot)]
  
    ass_tot <- data.frame(taxo=tax, seq=fa_tot)
  
    ass_tot$taxo <- gsub('[[:punct:]][[:digit:]]{2,3}[[:punct:]]{2}|;', '|', ass_tot$taxo)

    # correct the taxon that are found at many tax lev
    ass_tot$taxo <- gsub('Lineage_IIa|Elusimicrobia', 'Elusimicrobia|Lineage_IIa', ass_tot$taxo, fixed=T)
    ass_tot$taxo <- gsub('Candidatus_Magasanikbacteria|Parcubacteria', 'Candidatus_Magasanikbacteria|Candidatus_Magasanikbacteria_unclassified',
                         ass_tot$taxo, fixed=T)
    ass_tot$taxo <- gsub('TM7|Candidatus_Saccharibacteria', 'TM7|u', ass_tot$taxo, fixed=T)
    ass_tot$taxo <- gsub('Elusimicrobia|Lineage_IIa', 'Elusimicrobia|u', ass_tot$taxo, fixed=T)
    ass_tot$taxo <- gsub('_|', '|', ass_tot$taxo, fixed=T)
    ass_tot$taxo <- gsub('__', '_', ass_tot$taxo, fixed=T)
    
    # correct taxon with different parent taxon
    ass_tot <- ass_tot[names(mr_tot),]
  
    # make the taxonomy ----
    taxo_tot <- strsplit(as.character(ass_tot$taxo), '|' , fixed=T)
    nb_lev <- length(taxo_tot[[1]])
  
    taxo_tot <- matrix(unlist(taxo_tot), ncol=nb_lev, byrow=T)
    row.names(taxo_tot) <- row.names(ass_tot)
  
    taxo_tot <- as.data.frame(t(parApply(cl, taxo_tot, 1, function(x) {
      
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
  
    lst_ass_tot[[j]] <- list(ass_tot=ass_tot, taxo_tot=taxo_tot) 
  }
  
  # rename OTU according to dataset
  id <- switch(prim_names[i],
               '01_16S_bact' = 'bac',
               '02_18S_euk' = 'euk',
               '05_ITS_fun' = 'fun',
               '08_16S_cyano' = 'cya')
  
  names(mr_tot) <- sub('X_', paste0(id, '_'), names(mr_tot))
  
  lst_ass_tot <- lapply(lst_ass_tot, function(x) lapply(x, function(y){
    z <- y
    row.names(z) <- names(mr_tot)
    return(z)
  }))
  
  # sort taxo
  ass_tot <- lst_ass_tot[[1]]$ass_tot
  taxo_tot <- lst_ass_tot[[1]]$taxo_tot
  
  # improve the taxo with blastn best match of the OTUs against Unite
  if(i == 5){ 
    tt1 <- as.matrix(lst_ass_tot[[1]]$taxo_tot)
    tt2 <- as.matrix(lst_ass_tot[[2]]$taxo_tot)
    ass1 <- as.matrix(lst_ass_tot[[1]]$ass_tot)
    ass2 <- as.matrix(lst_ass_tot[[2]]$ass_tot)
    
    ind_FX <- tt1[,2] == 'Fungi_X'
    
    tt1[ind_FX,] <- tt2[ind_FX,]
    ass1[ind_FX,] <- ass2[ind_FX,]
    
    taxo_tot <- as.data.frame(tt1)
    ass_tot <- as.data.frame(ass1)
    
  }
  
  #---
  ord_taxo <- order(ass_tot$taxo)

  ass_sort <- ass_tot[ord_taxo,]
  mr_sort <- mr_tot[,ord_taxo]
  mr_sort <- mr_sort[rowSums(mr_sort) != 0,]
  taxo_sort <- taxo_tot[ord_taxo,]
  
  # taxo_clean
  taxo_false <- switch(i,
                       '1' = c('Chloroplast','Mitochondria'),
                       '2' = c('Metazoa','Embryophyceae'),
                       '3' = c('Life','TRUE_X'),
                       '4' = c('Life','TRUE_X'),
                       '5' = c('Life','Plantae','Chromista','Animalia','Protista'),
                       '6' = c('Life','TRUE_X'),
                       '7' = c('Life','TRUE_X'),
                       '8' = c('Life','TRUE_X','Bacteria_X','Chloroplast','Acidobacteriota','Actinobacteria','Actinobacteriota',
                               'Armatimonadetes','Chloroflexi','Firmicutes','GN02','Methylomirabilota','MVP-21','Nitrospinota',
                               'OD1','OP11','OP3','OP8','Patescibacteria','Planctomycetota','Planctomycetes','Proteobacteria',
                               'TM7','Verrucomicrobiota','WS4'),
                       '9' = c('Life','TRUE_X'))
  taxo_false <- paste0(taxo_false, collapse='|')
  ind_tf <- grepl(taxo_false, taxo_sort[,1]) | grepl(taxo_false, taxo_sort[,2]) | grepl(taxo_false, ass_sort$taxo)
  
  if(length(which(ind_tf))){
    ass_sort <- droplevels(ass_sort[ind_tf == F,])
    mr_sort <- mr_sort[,ind_tf == F]
    mr_sort <- mr_sort[rowSums(mr_sort) != 0,]
    taxo_sort <- droplevels(as.data.frame(taxo_sort)[ind_tf == F,])
  }
  
  env_sort <- env_tot[row.names(mr_sort),]
  
  # supress the last column of the taxo of bacteria (uninformative)
  if(i == 1){
    taxo_sort <- taxo_sort[,-ncol(taxo_sort)]
  }
  if(i == 8){
    taxo_sort <- taxo_sort[,-c(1,ncol(taxo_sort))]
  }
  
  # name tax_level
  names(taxo_sort) <- switch(as.character(i), 
                             '1' = c('reign','phylum','class','order','family','genus','species'),
                             '2' = c('reign','phylum','division','class','order','family','genus','species'),
                             '5' = c('reign','division','class','order','family','genus','species'),
                             '8' = c('phylum','class','order','family','genus'))
  
  # find low sequence samples (piecewise linear model) ----
  lrs <- log(sort(rowSums(mr_sort)))
  
  x <- brks <- seq_along(lrs)
  
  mse <- lh_steep <- NULL
  for(j in brks){
    mod <- lm(lrs~x*(x <= brks[j]) + x*(x < brks[j]))
    mse <- c(mse, summary(mod)$sigma)
    
    co <- coef(mod)
    lh_steep <- c(lh_steep, co[2]+co[5] > co[2])
  }
  
  min_mse <- which(mse == min(mse[lh_steep], na.rm=T))
  
  mod <- lm(lrs ~ x*(x < min_mse) + x*(x > min_mse))
  co <- coef(mod)
  
  low_seq <- names(lrs)[1:(min_mse-1)]
  
  env_sort$low_seq <- row.names(env_sort) %in% low_seq
  
  #---
  cairo_ps(paste0(dir_prep, 'piecewise_', prim_names[i], '.eps'))
  plot(lrs~x, xlab='', ylab='log of rowSums', main=prim_names[i])
  
  curve(co[1]+co[3] + (co[2]+co[5])*x, add=T, from=1, to=min_mse)
  curve(co[1]+co[4] + (co[2])*x, add=T, from=min_mse, to=length(lrs))
  abline(v=min_mse, lty=3)
  
  dev.off()
  
  # communities building ----
  # without the low_sequences samples (piecewise regression)
  mr_nls <- mr_sort[env_sort$low_seq == F,]
  mr_nls <- mr_nls[,colSums(mr_nls) != 0]
  
  env_nls <- env_sort[row.names(mr_nls),]
  
  ass_nls <- ass_sort[names(mr_nls),]
  taxo_nls <- taxo_sort[names(mr_nls),]
  
  # centered log ratio on nls
  mr_rcz <- mr_nls[,colSums(decostand(mr_nls, 'pa')) > 1] # take out the OTU found in only one sample otherwise, the replacment of 0 give negative values
  system.time(mr_rcz <- cmultRepl(mr_rcz, method='CZM', output='p-counts')) #30 sec
  mr_clr_nls <- as.data.frame(t(apply(mr_rcz, 1, function(x) log(x) - mean(log(x)) )))
  
  env_clr_nls <- env_tot[row.names(mr_clr_nls),]
  
  ass_clr_nls <- ass_sort[names(mr_clr_nls),]
  taxo_clr_nls <- taxo_sort[names(mr_clr_nls),]
  
  #---
  lst <- list(raw = list(env=env_sort, mr=mr_sort, ass=ass_sort, taxo=taxo_sort, lst_ass_tot=lst_ass_tot, lst_ass_ini=lst_ass_ini),
              nls = list(env=env_nls,  mr=mr_nls,  ass=ass_nls,  taxo=taxo_nls),
              clr_nls = list(env=env_clr_nls, mr=mr_clr_nls, ass=ass_clr_nls, taxo=taxo_clr_nls)
              )

  # return(lst)
  lst_comm[[prim_names[i]]] <- lst
  
  file <- paste0(dir_save, '00_lst_comm_', prim_names[i], '.Rdata')
  save(lst, file=file)
}

names(lst_comm) <- n_comm

#---
file <- paste0(dir_save, '00_lst_comm.Rdata')
save(lst_comm, env_tot, lst_palev, permu, fact_3, n_comm, leg, lay, file=file)

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

for(i in c(1,2,5,8)){
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
par(mfrow=c(2,4))

ra_tot <- NULL
for(i in c(1,2,5,8)){
  out_bf <- read.table(paste0('Projets/Climarctic/bioinfo/archive/200525/cnt_output', i),
                       h=T, sep='\t', row.names=1)

  out_bf <- out_bf[order(row.names(out_bf)),]
  
  B <- switch(as.character(i),
              '1' = 190:194,
              '2' = 127:129,
              '5' = 121:129,
              '8' = 82:84)
  
  print(colSums(out_bf[-B,]))

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

pal_prim <- brewer.pal(4, 'Set1')

for(i in 1:2){
  ind <- seq(i, 8, by=2)

  df <- as.data.frame(t(ra_tot[,ind]))

  plot(NA, xlim=c(1,ncol(df)), xaxt='n', xlab='',
       ylim=range(df, na.rm=T), ylab=ifelse(i == 1, 'nb seq', '%nb_seq'),
       log=ifelse(i == 1, '', 'y'), main=paste(c('raw','relabu')[i], 'total dataset'))
  axis(1, at=1:ncol(df), labels=names(df), las=2)

  for(j in 1:nrow(df)){
    lines(1:5, df[j,], col=pal_prim[j])
  }

  if(i == 2)legend('bottomleft', legend=prim_names[c(1,2,5,8)], pch=19, col=pal_prim, bty='n', cex=0.5)
}

dev.off()



#####




















