#####
# climarctic environmental variables
#####

rm(list=ls())
gc()

require(randomForest)
require(dendextend)
require(PMCMR)
require(multcompView)


# dir loads ####
dir_out  <- 'Projets/Climarctic/stats/MBC/out/'
dir_save <- paste0(dir_out, 'saves/') 
dir_env  <- paste0(dir_out, '02_env/')
dir.create(dir_env, showWarnings=F)

#---
file <- paste0(dir_save, '00_lst_comm.Rdata')
load(file)


# select variables of interest
env_sub <- na.omit(env_tot[,c('site','moisture','depth',
                              'pH','C_N','sand','silt','clay','rh','om')])

is_num <- sapply(env_sub, is.numeric)

env_sc <- env_sub
env_sc[,is_num] <- scale(env_sc[,is_num])

# random forest ####
rf <- randomForest(env_sc)

clust <- hclust(as.dist(1-rf$proximity))
dend <- as.dendrogram(clust)
ord <- order.dendrogram(dend)

imp <- importance(rf)
print(imp[order(imp[,1]),])

#---
# pdf(paste0(dir_env, 'random_forest_dend.pdf'), width=17, height=7)
cairo_ps(paste0(dir_env, 'random_forest_dend.eps'), width=17, height=7)
par(mar=c(15,4,4,2))

plot(dend)
usr <- par('usr')

for(i in 1:3){
  f <- c('site','moisture','depth')[i]
  col <- lst_palev[[f]]$pal[as.numeric(env_sub[[f]])][order.dendrogram(dend)]
  points(1:nrow(env_sub), rep(usr[3]-(diff(usr[3:4])*0.05)*(i)-0.1, nrow(env_sub)), xpd=NA, pch=19, col=col)
}

legend(usr[1]+diff(usr[1:2])*0.5, usr[3]-diff(usr[3:4])*0.6, legend=c(unlist(sapply(lst_palev, '[[', 2))), bty='n', xjust=0.5, yjust=0.5,
       pch=19, col=unlist(sapply(lst_palev, '[[', 1)), xpd=NA)

dev.off()


# PCA ####
pca <- rda(env_sc[,-c(1:3)])
s <- summary(pca)

coord_smp <- s$site[,1:2]
coord_facvar <- s$species[,1:2]

perc_var <- signif(pca$CA$eig/sum(pca$CA$eig), 2)

#---
# pdf(paste0(dir_env, 'pca.pdf'), width=15, height=15)
cairo_ps(paste0(dir_env, 'pca.eps'), width=10, height=7)
layout(matrix(c(1,3,4, 1,5,6, 0,2,2), nrow=3, byrow=T), width=c(0.2,1,1), height=c(1,1,0.2))
par(mar=rep(0,4))

for(i in 2:1){
  plot.new()
  text(0.5,0.5,labels=paste0('PC', i, ' ', perc_var[i], '%'), srt=ifelse(i == 1, 0, 90))
}

for(i in c('site','moisture','depth')){
  plot(rbind(coord_smp, coord_facvar), type='n')
  abline(v=0,h=0,lty=3)
  
  arrows(0,0,coord_facvar[,1],coord_facvar[,2], length=0, col='grey70', lty=2)
  text(coord_facvar, labels=row.names(coord_facvar))
  
  text(coord_smp, labels=row.names(env_sc), col=lst_palev[[i]]$pal[as.numeric(env_sc[[i]])])
}

plot.new()
legend(0.5,0.5,legend=c(unlist(sapply(lst_palev, '[[', 2))), bty='n', xjust=0.5, yjust=0.5,
       pch=19, col=unlist(sapply(lst_palev, '[[', 1)))

dev.off()


# var vs factors ####
var <- env_sub[,is_num]
fact <- env_sub[,c('site','moisture','depth')]

# pdf(paste0(dir_env, 'var_vs_fact.pdf'), width=15, height=8)
cairo_ps(paste0(dir_env, 'var_vs_fact.eps'), width=15, height=8)
layout(matrix(1:(ncol(var)*4), byrow=T, nrow=4), width=rep(1, ncol(var)), height=c(0.2,1,1,1))
par(mar=c(0,3,0,1))

for(i in names(var)){
  plot.new()
  text(0.5,0.5, i, font=2)
}

par(mar=c(3,3,1,1))
for(i in names(fact)){
  for(j in names(var)){
    
    f <- env_sub[[i]]
    v <- env_sub[[j]]
    
    lev <- lst_palev[[i]]$lev
    
    # test per factor and interaction      
    pvs <- posthoc.kruskal.nemenyi.test(v, f)$p.value
    pvs <- as.dist(cbind(rbind(rep(NA,ncol(pvs)), pvs), rep(NA, nrow(pvs)+1)))
    attributes(pvs)$Labels <- lev
    mcl <- multcompLetters(pvs)$Letters
    
    #---
    boxplot(v~f, ylim=range(v)*c(1,1.2))
    usr <- par('usr')
    
    text(1:length(lev), usr[4]-diff(usr[3:4])*0.1, mcl)
  }
}

dev.off()
#




#























