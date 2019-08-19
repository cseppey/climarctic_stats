#####
# climarctic preparation
#####

rm(list=ls())
gc()


# dir loads ####
dir_out  <- 'Projets/Climarctic/stats/out/'
dir_save <- paste0(dir_out, 'saves/') 
dir_rda  <- paste0(dir_out, '02_rda/')
dir.create(dir_rda, showWarnings=F)

#---
file <- paste0(dir_save, 'lst_comm.Rdata')
load(file)

# prep ####
mr_abrel   <- decostand(mr_sort, method='total')
mr_ar_hell <- decostand(mr_abrel, method='hell')
  
bc      <- vegdist(mr_abrel)
bc_hell <- vegdist(mr_ar_hell)

lst_rda <- list(abrel  =list(simple=capscale(bc     ~site+moisture+depth, data=env_sort), 
                             inter =capscale(bc     ~site*moisture*depth, data=env_sort)),
                ar_hell=list(simple=capscale(bc_hell~site+moisture+depth, data=env_sort), 
                             inter =capscale(bc_hell~site*moisture*depth, data=env_sort)))

lst_pvs_terms <- lapply(lst_rda, function(x) lapply(x, function(y) anova(y, by='terms')$`Pr(>F)`))

#---
pdf(paste0(dir_rda, 'rda.pdf'), width=15, height=10)
layout(matrix(c(1,2,5,3,4,0), nrow=2, byrow=T))

for(i in seq_along(lst_rda)){
  for(j in seq_along(lst_rda[[i]])){
    
    rda <- lst_rda[[i]][[j]]
    
    pv <- signif(anova(rda)$`Pr(>F)`[1], 2)
    aic <- signif(extractAIC(rda)[2], 2)
    
    s <- summary(rda)
    
    site <- s$sites
    
    var <- signif(s$cont$importance[2,1:2]*100, 2)
      
    #---
    plot(site[,1:2], type='n', main=paste(names(lst_rda)[i], names(lst_rda[[i]])[j]),
         xlab=paste('RDA1 var:', var[1], '%'), ylab=paste('RDA2 var:', var[2], '%'))
    usr <- par('usr')
    
    ordispider(rda,env_sort$moisture, col=lst_palev$moisture$pal)
    
    points(site[,1:2], pch=c(16,17)[as.numeric(env_sort$depth)], col=lst_palev$site$pal[as.numeric(env_sort$site)])
    
    text(usr[2]-diff(usr[1:2])*0.1, usr[3]+diff(usr[3:4])*0.9, paste0('pv: ', pv, '\nAIC: ', aic), cex=0.75)
  }
}

plot.new()
legend(0.5,0.5, unlist(sapply(lst_palev, function(x) x$lev)), 
       col=c(unlist(sapply(lst_palev[1:2], function(x) x$pal)), 1, 1),
       pch=c(rep(19, 5), 16,17), bty='n', xjust=0.5, yjust=0.5)

dev.off()






