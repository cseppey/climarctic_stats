#####
# climarctic alpha div
#####

print('##### Climarctic 05 alpha #####')

rm(list=ls())

# require(betapart)
require(foreach)
require(doSNOW) # e.g. makeSOCKcluster()
require(abind)
require(fossil)
require(MuMIn)
require(lme4)
require(lmerTest)

# dir loads ####
dir_in <- 'Projets/Climarctic/stats/MBC/in/'
dir_out <- 'Projets/Climarctic/stats/MBC/out/'
dir_div <- paste0(dir_out, '05_div/') 
dir_save <- paste0(dir_out, 'saves/') 
dir.create(dir_div, showWarnings=F, recursive=T)

file <- paste0(dir_save, '00_lst_comm.Rdata')
load(file)

#---
lst_pvs_div <- NULL
for(i in n_comm){
  
  print(i)
  
  env <- lst_comm[[i]]$nls$env
  
  rh <- env$rh
  site <- env$site
  moisture <- env$moisture
  depth <- env$depth
  
  mr <- lst_comm[[i]]$nls$mr
  
  # factors
  lev_fac <- sapply(lst_palev, names)
  
  ls <- lev_fac$site
  ld <- lev_fac$depth
  lm <- lev_fac$moisture
  
  # whole dataset vs rh ----
  
  # alpha
  alpha <- apply(mr, 1, function(x) chao1(x))
  
  # correct according to the number of sequence
  log_nb_seq <- log(rowSums(mr))
  lm <- lm(alpha~log_nb_seq)
  alpha_r <- residuals(lm)
  
  env <- cbind.data.frame(env, alpha_r)
  
  print(table(env$site, env$depth))
  
  plot(alpha_r~rh, pch=c(15,16)[env$site], col=lst_palev$depth[env$depth])
  for(d  in c('top','deep')){
    for(s in c('Knudsenheia','Ossian')){
      lm <- lm(alpha_r~rh, data=env[env$depth == d & env$site == s,])
      abline(coef(lm), col=lst_palev$depth[[d]], lty=ifelse(s == 'Knudsenheia', 2, 3))
      print(anova(lm))
    }
  }
  
  # linear mixed effect model (nlme: https://www.r-bloggers.com/linear-mixed-effect-models-in-r/)
  # interesting site on mixed effect models (lme4: http://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#model-definition)
  
  for(j in c('top|deep','top','deep')){
    e <- env[grep(j, env$depth),]
    
    lst_lm <- NULL
    
    f0 <- formula(alpha_r~rh)
    lst_lm[['lm0']] <- lm0 <- lm(f0, data=e)
    
    lst_lm[['lms']] <- lms <- lm(alpha_r~rh*site, data=e)
    
    if(j == 'top|deep'){
      f1 <- update(f0, ~. * depth)
      lst_lm[['lmd']] <- lmd <- lm(f1, data=e)
    } else {
      f1 <- f0
    }
  
    ff1 <- update(f1, ~.+(rh|site))
    lst_lm[['lme1']] <- lme1 <- lmer(ff1, data=e)
    
    # test
    
    lst_pvs_div[[i]][[j]] <- pvs <- lapply(lst_lm, function(x) {
      a <- anova(x)
      p <- a$'Pr(>F)'
      names(p) <- row.names(a)
      return(c(p[is.na(p) == F], norm=shapiro.test(resid(x))$p.value))
    })
    
    title <- sapply(pvs, function(x) paste(names(x), round(x, 3), collapse=', '))
    #
    
    #---
    # cairo_ps(paste0(dir_div, 'lm_poly1_', i, '_', j, '.eps'), width=10, height=7)
    svg(paste0(dir_div, 'lm_poly1_', i, '_', j, '.svg'), width=10, height=7)
    lay(5, c('rh','residual of\nlm(Chao~log(seq nb)'), lab_crd=list(x=c(0.5, 0.33), y=c(0.33,0.5)))
    
    for(k in 1:3){
      plot(alpha_r~rh, data=e, pch=c(15,16)[e$site], yaxt=ifelse(k == 2, 'n', 's'),col=lst_palev$depth[e$depth])
    
      cm <- 0.8
      ln <- 1
      
      abline(coef(lm0))
      
      # vs depth
      if(k == 1){
        if(j == 'top|deep'){
          title(paste(title[c('lm0','lmd')], collapse='\n'), cex.main=cm, line=ln)
          
          abline(coef(lmd)[1:2], col=lst_palev$depth[1])
          abline(sum(coef(lmd)[c(1,3)]), sum(coef(lmd)[c(2,4)]), col=lst_palev$depth[2])
        } else {
          title(title['lm0'], cex.main=cm, line=ln)
        }
      }
      
      # vs site
      if(k == 2){
        title(paste(title[c('lm0','lms')], collapse='\n'), cex.main=cm, line=ln)
        
        abline(coef(lms)[1:2], lty=2, col=lst_palev$site[1])
        abline(sum(coef(lms)[c(1,3)]), lty=3, sum(coef(lms)[c(2,4)]), col=lst_palev$site[2])
      }
      
      # lme
      if(k == 3){
        title(paste(title[c('lm0','lme1')], collapse='\n'), cex.main=cm, line=ln)
        
        cfme <- as.matrix(coef(lme1)$site)
        for(k in 1:2){
          if(j == 'top|deep'){
            abline(cfme[k,1:2], lty=k+1, col=lst_palev$depth[1])
            abline(sum(cfme[k,c(1,3)]), sum(cfme[k,c(2,4)]), lty=k+1, col=lst_palev$depth[2])
          } else {
            abline(cfme[k,], lty=k+1, col=lst_palev$depth[[j]])
          }
        }
      }
    }
    
    leg(title=i, lay='mixed', bg=unlist(lst_palev), lty=c(2,3,rep(0,5)))
    
    #---
    dev.off()
    
  }
  
  
}

#---
lapply(lst_pvs_div, function(x){
  lapply(x, function(y) {
    sapply(y, function(z) {
      gds <- grepl('depth|site', names(z)) 
      grh <- grepl('rh', names(z))
      
      p1 <- z[1]
      
      p2 <- z[gds & grh == F]
      p2 <- ifelse(length(p2), p2, NA)
      
      p3 <- z[gds & grh]
      p3 <- ifelse(length(p3), p3, NA)
      
      p4 <- z[names(z) == 'norm']
      
      return(c(p1, p2, p3, p4))
    })
  })
}) 







