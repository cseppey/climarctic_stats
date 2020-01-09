#####
# climarctic preparation
#####

rm(list=ls())
gc()

require(foreach)
require(doSNOW) # e.g. makeSOCKcluster()
require(RColorBrewer) # brewer.pal 
require(zCompositions)
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
env_ini_gas_raw  <- read.table(paste0(dir_in, 'env_clim_gas_raw.csv'), h=T)
env_ini_gas  <- read.table(paste0(dir_in, 'env_clim_gas.csv'), h=T)

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

# add noise to texture
tex <- env_ini_chim[,c('sand','silt','clay')]
set.seed(0)
tex <- as.data.frame(sapply(tex, function(x) x+rnorm(nrow(tex), sd=min(tex, na.rm=T)/1000)))
names(tex) <- paste0(names(tex), '_nz')

# env_gas %%%
env_gas <- env_ini_gas_raw
env_gas$Plot <- factor(env_gas$Plot)
env_gas$Measurment <- factor(env_gas$Measurment)

#---
pdf(paste0(dir_prep, 'gas_measurements.pdf'), width=40, height=40)
par(mfrow=c(12,12), mar=c(4,5,3,8))

arr_gas <- NULL
nm <- NULL

for(e in levels(env_gas$Site)){
  ae <- NULL
  for(f in levels(env_gas$Habitat)){
    af <- NULL
    for(g in levels(env_gas$Plot)){
      ag <- NULL
      for(h in levels(env_gas$Condition)){
        ah <- NULL
        for(i in levels(env_gas$Measurment)){

          non_clim <- e == 'Knudsenheia_mette' | e == 'Sulvaten'
          
          if(non_clim){
            ind <- which(env_gas$Site == e 
                         & env_gas$Plot == g)
          } else {
            ind <- which(env_gas$Site == e 
                         & env_gas$Habitat == f 
                         & env_gas$Plot == g 
                         & env_gas$Condition == h 
                         & env_gas$Measurment == i)
          }
          
          df <- na.omit(data.frame(env_gas[ind,c('CH4','N2O','CO2')], time = env_gas$time[ind]))
          
          lms <- as.data.frame(sapply(df[,1:3], function(x) {
            lm <- lm(x~df$time)
            coef <- lm$coefficients
            R2_P <- c(R2=summary(lm)$r.squared, P=ifelse(length(x) == 3, cor.test(x,df$time)$p.value, NA))
            
            return(c(coef,R2_P))
          }))
          row.names(lms) <- c('intercept','slp','R2','P')
          
          #---
          plot(df$CH4~df$time, xlab='minute', ylab='CH4', main=paste(c(e,f,g,h,i), collapse='_'), pch=15)
          
          abline(lms$CH4[1:2])
          
          for(j in 2:3){
            par(new=T)
            plot(df[[j]]~df$time, axes=F, bty='n', xlab='', ylab='', col=j, pch=c(16,17)[j-1])
            
            axis(side=4, line=(j-2)*4, at=pretty(range(df[[j]])))
            mtext(names(df)[j], side=4, line=(j-2)*4+2, col=j, cex=0.6)
            abline(lms[[j]][1:2], col=j)
            
          }
          usr <- par('usr')
          
          ys <- usr[3]+diff(usr[3:4])*seq(0.15,0.1, length.out=3)
          text(usr[1]+diff(usr[1:2])*0.1, ys, c('slp:','R2:','P:'), cex=0.5)
          
          for(j in 1:3){
            text(usr[1]+diff(usr[1:2])*seq(0.2,0.5, length.out=3)[j], ys, signif(lms[2:4,j], 2), col=j, cex=0.5)
          }
          
          #---
          ah <- abind(ah, lms, along=3)
        }
        dimnames(ah)[[3]] <- paste0('mes_',levels(env_gas$Measurment))
        ag <- abind(ag, ah, along=4)
      }
      dimnames(ag)[[4]] <- levels(env_gas$Condition)
      af <- abind(af, ag, along=5)
    }
    dimnames(af)[[5]] <- paste0('plot_', levels(env_gas$Plot))
    ae <- abind(ae, af, along=6)
  }
  dimnames(ae)[[6]] <- levels(env_gas$Habitat)
  arr_gas <- abind(arr_gas, ae, along=7)
}
dimnames(arr_gas)[[7]] <- levels(env_gas$Site)

dev.off()

# sort a bit the fluxes
flux <- as.data.frame.table(arr_gas[c(2,4),,,'light',,,c(1,3)])
flux <- cbind.data.frame(flux[seq(1, nrow(flux), 2),-c(1,ncol(flux))], 
                         slp=flux$Freq[seq(1, nrow(flux), 2)], P=flux$Freq[seq(2, nrow(flux), 2)])
fl <- flux[seq(1, nrow(flux), 3),sapply(flux, is.factor)][,-1]
for(i in levels(flux$Var2)){
  fl <- cbind.data.frame(fl, flux[flux$Var2 == i,c('slp','P')])
}
names(fl)[sapply(fl, function(x) is.factor(x) == F)] <- apply(expand.grid(c('slp','P'), levels(flux$Var2)), 
                                                          1, function(x) paste(x[2:1], collapse='_'))

fl_best <- NULL
for(i in seq(1, nrow(fl), 2)){
  var <- NULL
  for(j in seq(6, ncol(fl), 2)){
    ind <- c(i:(i+1))
    Ps <- fl[ind,j]
    var <- c(var, fl[ind[which(Ps == min(Ps, na.rm=T))], j-1])
  }
  fl_best <- rbind(fl_best, var)
}
fl_best <- cbind.data.frame(fl[seq(1,nrow(fl),2),2:4], fl_best)
names(fl_best)[4:6] <- c('CH4','N2O','CO2')

#---
flux <- fl_best[gl(18,6),4:6]

# take the interesting variables
env_tot <- cbind.data.frame(env_ini_fact[,c(2:5,7:15)], env_ini_chim[,c('sand','silt','clay','NO3','NH4','P_H2O','P_NAHCO3','P_labile')],
                            tex, flux)
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
    ass_ini <- read.table(paste0(dir_in, 'from_cluster/', id_plate, '/', id_plate, '_clust.ass'), row.names=1)
    fa_ini  <- read.table(paste0(dir_in, 'from_cluster/', id_plate, '/', id_plate, '_clust.fa'))
    save(mr_ini, ass_ini, fa_ini, file=file) 
  }
  
  # reorganize mr
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
  mr_tot <- mr_tot[,colSums(mr_tot) != 0]
  
  print(c(sum(mr_tot), ncol(mr_tot)))

  # reorganize fa
  n_fa <- as.character(fa_ini[seq(1,nrow(fa_ini), by=2),])
  n_fa <- substr(n_fa, 2, nchar(n_fa))

  fa_tot <- fa_ini[seq(2,nrow(fa_ini), by=2),]
  names(fa_tot) <- n_fa

  # reorganize ass
  ass_tot <- data.frame(taxo=ass_ini$V2, seq=fa_tot)

  ass_tot$taxo <- gsub('[[:punct:]][[:digit:]]{2,3}[[:punct:]]{2}|;', '|', ass_tot$taxo)

  ass_tot <- ass_tot[names(mr_tot),]

  taxo_tot <- strsplit(as.character(ass_tot$taxo), '|' , fixed=T)
  nb_lev <- length(taxo_tot[[1]])

  taxo_tot <- matrix(unlist(taxo_tot), ncol=nb_lev, byrow=T)
  row.names(taxo_tot) <- row.names(ass_tot)

  taxo_tot <- as.data.frame(t(apply(taxo_tot, 1, function(x) {
    
    x <- gsub('.*Incertae_.*|Unknown|.*-[pcofgs]$|.*_unclassified', 'u', x)
    x <- gsub('_Incertae','', x)

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
  
  if(i != 1 & i != 2 & i != 5 & i != 8){
    ind_true <- taxo_sort$V1 == 'TRUE'
    
    ass_sort <- droplevels(ass_sort[ind_true,])
    mr_sort <- mr_sort[,ind_true]
    taxo_sort <- droplevels(taxo_sort[ind_true,-1])
  }
  
  # taxo_clean
  taxo_false <- switch(i,
                       '1' = 'Chloroplast',
                       '2' = 'Metazoa|Embryophyceae',
                       '3' = 'Life|TRUE_X',
                       '4' = 'Life|TRUE_X',
                       '5' = 'Plantae',
                       '6' = 'Life|TRUE_X',
                       '7' = 'Life|TRUE_X',
                       '8' = 'Life|TRUE_X|Sericytochromatia|Vampirivibrionia|Bacteria_X|Chloroplast|Acidobacteriota|Actinobacteriota|Chloroflexi|Firmicutes|Methylomirabilota|Nitrospinota|Patescibacteria|Planctomycetota|Verrucomicrobiota|WS4',
                       '9' = 'Life|TRUE_X')
  ind_tf <- grepl(taxo_false, taxo_sort[,1]) | grepl(taxo_false, taxo_sort[,2]) | grepl(taxo_false, ass_sort$taxo)
  
  if(length(which(ind_tf))){
    ass_sort <- droplevels(ass_sort[ind_tf == F,])
    mr_sort <- mr_sort[,ind_tf == F]
    taxo_sort <- droplevels(taxo_sort[ind_tf == F,])
  }
  
  env_sort <- env_tot[row.names(mr_sort),]
  
  # communities normalisation ----
  rs <- rowSums(mr_sort)
  thresh <- seq(min(rs), max(rs), length.out=20)
  thresh <- thresh[-length(thresh)]
  
  # nb seq
  plot(sort(rs), main=prim_names[i])
  abline(h=thresh)
  
  # calculation of percentage of sequence and OTU lost
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
  
  # combinatorial
  mr_cmb <- mr_sort[,colSums(decostand(mr_sort, 'pa')) > 1]
  mr_cmb <- mr_cmb[rowSums(mr_cmb) != 0,]
  mr_cmb <- cmultRepl(mr_cmb)
  
  mr_cmb <- as.data.frame(t(apply(mr_cmb, 1, function(x) {
    G <- exp(mean(log(x)))
    return(sapply(x, function(y) y/G))
  })))
  
  env_cmb <- env_tot[row.names(mr_cmb),]

  ass_cmb <- ass_sort[names(mr_cmb),]
  taxo_cmb <- taxo_sort[names(mr_cmb),]
  
  # rrf
  set.seed(0)
  mr_rrf <- rrarefy(mr_sort[rs >= op[3],], op[3])
  mr_rrf <- as.data.frame(mr_rrf[,colSums(mr_rrf) != 0])
  
  env_rrf <- env_tot[row.names(mr_rrf),]
  
  ass_rrf <- ass_sort[names(mr_rrf),]
  taxo_rrf <- taxo_sort[names(mr_rrf),]

  # rrf2 
  set.seed(0)
  mr_rrf2 <- rrarefy(mr_sort[rs >= op[1],], op[1])
  mr_rrf2 <- as.data.frame(mr_rrf2[,colSums(mr_rrf2) != 0])
  
  env_rrf2 <- env_tot[row.names(mr_rrf2),]
  
  ass_rrf2 <- ass_sort[names(mr_rrf2),]
  taxo_rrf2 <- taxo_sort[names(mr_rrf2),]
  
  lst_comm[[prim_names[i]]] <- list(raw= list(env=env_sort, mr=mr_sort, ass=ass_sort, taxo=taxo_sort),
                                    cmb= list(env=env_cmb , mr=mr_cmb , ass=ass_cmb , taxo=taxo_cmb),
                                    rrf= list(env=env_rrf , mr=mr_rrf , ass=ass_rrf , taxo=taxo_rrf),
                                    rrf2=list(env=env_rrf2, mr=mr_rrf2, ass=ass_rrf2, taxo=taxo_rrf2))

}

#---
file <- paste0(dir_save, '00_lst_comm.Rdata')
save(lst_comm, env_tot, lst_palev, permu, file=file)

file <- paste0(dir_save, '00_env_tot.Rdata')
save(env_tot, file=file)

stopCluster(cl)

stop('pouet')

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
  out_bf <- read.table(paste0('Projets/Climarctic/bioinfo/archive/next/filter_test/cnt_output/cnt_output', i),
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

  #---
  pdf(paste0(dir_prep, 'bioinfo_check_', i, '.pdf'), width=8, height=4)
  par(mfrow=c(1,2))
  
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
  
  dev.off()
}

# total

pdf(paste0(dir_prep, 'bioinfo_check_tot.pdf'), width=8, height=4)
par(mfrow=c(1,2))

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




















