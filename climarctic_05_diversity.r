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
require(lme4)

# dir loads ####
dir_in <- 'Projets/Climarctic/stats/MBC/in/'
dir_out <- 'Projets/Climarctic/stats/MBC/out/'
dir_div <- paste0(dir_out, '05_div/') 
dir_save <- paste0(dir_out, 'saves/') 
dir.create(dir_div, showWarnings=F, recursive=T)

file <- paste0(dir_save, '00_lst_comm.Rdata')
load(file)

# gamma ---
#   specpool$chao
#   whole dataset
#   per site*moisture*depth

## from cluster
# lst_sp <- NULL
# for(i in n_comm){
#   lst_sp[[i]] <- specpool(lst_comm[[i]]$nls$mr) 
# }
# 
# write.tavle(t(sapply(lst_co, function(x) x)), 05_sp_csv, F, F, '\t')

specpool <- read.table(paste0(dir_div, 'from_cluster/05_sp.csv'))

# per sample to compare to the rh
#   alpha = chao
#   beta = gamma_whole_dataset/alpha
# per site*moisture*depth to know where is alpha and beta the highest
#   alpha = distro of alpha per smd
#   beta = gamma_smd / alpha

for(i in n_comm){
  
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
  # gamma 
  # (from cluster for the whole community)
  
  gamma <- specpool[i,'chao']
  
  # alpha
  # alpha <- diversity(mr)
  alpha <- apply(mr, 1, function(x) chao1(x))
  
  # correct according to the number of sequence and set between 0 and 1
  log_nb_seq <- log(rowSums(mr))
  lm <- lm(alpha~log_nb_seq)
  alpha_r <- residuals(lm)
  
  alpha_r <- alpha_r-min(alpha_r)
  alpha_r <- alpha_r/max(alpha_r)

  # beta
  beta <- gamma/alpha_r
  
  df_div <- data.frame(alpha=alpha_r, beta=beta)
  
  # models ----
  for(j in names(df_div)){
    
    metric <- df_div[[j]]
  
    # # linear mixed effect model (https://www.r-bloggers.com/linear-mixed-effect-models-in-r/)
    # 
    # # OPTIMAL RANDOM STRUCTURE
    # lst_mod <- list( glm = gls(metric ~ rh,                              data=env, method = "ML"),
    #                 lmm1 = lme(metric ~ rh, random = ~1|site           , data=env, method = "ML"),
    #                 lmm2 = lme(metric ~ rh, random = ~1|moisture       , data=env, method = "ML"),
    #                 lmm3 = lme(metric ~ rh, random = ~1|depth          , data=env, method = "ML"),
    #                 lmm4 = lme(metric ~ rh, random = ~1|combi          , data=env, method = "ML"),
    #                 lmm5 = lme(metric ~ rh, random = ~1|site/MiS       , data=env, method = "ML"),
    #                 lmm6 = lme(metric ~ rh, random = ~1|MiS/PiMiS      , data=env, method = "ML"),
    #                 lmm7 = lme(metric ~ rh, random = ~1|site/MiS/PiMiS , data=env, method = "ML"))
    # 
    # anova(lst_mod[[1]], lst_mod[[2]], lst_mod[[3]], lst_mod[[4]],
    #       lst_mod[[5]], lst_mod[[6]], lst_mod[[7]], lst_mod[[8]])
    # sapply(lst_mod, function(x) shapiro.test(resid(x))$p.value)
    # 

    
    # vs rh linear----
    cairo_ps(paste0(dir_div, 'lm_poly1_', i, '_', j, '.eps'), 11, 8)
    par(mfrow=c(2,2))
    for(k in c('site','moisture','depth')){
      
      plot(metric~rh, ylab=j, pch=19, col=lst_palev[[k]][env[[k]]])
      
      # overall
      lm <- lm(metric~rh, data=env)
      
      pvs <- anova(lm)$`Pr(>F)`
      norm <- shapiro.test(resid(lm))$p.value
      
      abline(coef(lm))
      
      # no interac
      formu <- formula(paste('metric~rh+', k))
      lm <- lm(formu, data=env)
      
      pvs <- c(pvs, anova(lm)$`Pr(>F)`)
      norm <- c(norm, shapiro.test(resid(lm))$p.value)
      
      cf <- coef(lm)
      
      for(l in seq_along(lev_fac[[k]])){
        if(l == 1){
          abline(cf[1], cf[2], col=lst_palev[[k]][l], lty=2)
        } else {
          abline(cf[1]+cf[1+l], cf[2], col=lst_palev[[k]][l], lty=2)
        }
      }
      
      # with interac
      formu <- formula(paste('metric~rh*', k))
      lm <- lm(formu, data=env)
      
      pvs <- c(pvs,anova(lm)$`Pr(>F)`)
      norm <- c(norm, shapiro.test(resid(lm))$p.value)
      
      cf <- coef(lm)
      
      for(l in seq_along(lev_fac[[k]])){
        if(l == 1){
          abline(cf[1], cf[2], col=lst_palev[[k]][l], lty=3)
        } else {
          abline(cf[1]+cf[1+l], cf[2]+cf[length(lev_fac[[k]])+l], col=lst_palev[[k]][l], lty=3)
        }
      }
      
      # pvs
      pvs <- round(na.omit(pvs),2)
      norm <- round(norm,2)
      mtext(paste0('overall: P rh=', pvs[1], ', norm:', norm[1],
                   '\n+fact: P rh=', pvs[2], ' fact=', pvs[3], ', norm:', norm[2],
                  '\n*fact: P rh=',  pvs[4], ' fact=', pvs[5], ' inter=', pvs[6], ', norm=', norm[3]), line=0.5, cex=0.75)
    }
    leg(i)
    
    #---
    dev.off()
    
    # vs rh poly2
    cairo_ps(paste0(dir_div, 'lm_poly2_', i, '_', j, '.eps'), 11, 8)
    par(mfrow=c(2,2))
    for(k in c('site','moisture','depth')){
      
      lf <- length(lev_fac[[k]])
      
      plot(metric~rh, ylab=j, pch=19, col=lst_palev[[k]][env[[k]]])
      
      pred_x <- seq(min(rh), max(rh), length.out=100)
      
      # overall
      lm <- lm(metric~poly(rh, 2, raw=T), data=env)
      pvs <- anova(lm)$`Pr(>F)`
      
      cf <- coef(lm)
      pred_y <- cf[1] + cf[2]*pred_x + cf[3]*pred_x^2
      lines(pred_y~pred_x)
      
      # no interac
      formu <- formula(paste('metric~poly(rh,2,raw=T)+', k))
      lm <- lm(formu, data=env)
      
      pvs <- c(pvs, anova(lm)$`Pr(>F)`)
      
      cf <- coef(lm)
      
      for(l in seq_along(lev_fac[[k]])){
        if(l == 1){
          pred_y <- cf[1] + cf[2]*pred_x + cf[3]*pred_x^2
        } else {
          pred_y <- cf[1]+cf[2+l] + cf[2]*pred_x + cf[3]*pred_x^2
        }
        lines(pred_y~pred_x, col=lst_palev[[k]][l], lty=2)
      }
      
      # with interac
      formu <- formula(paste('metric~poly(rh, 2, raw=T)*', k))
      lm <- lm(formu, data=env)
      
      pvs <- c(pvs,anova(lm)$`Pr(>F)`)
      
      cf <- coef(lm)
      
      for(l in seq_along(lev_fac[[k]])){
        if(l == 1){
          pred_y <- cf[1] + cf[2]*pred_x + cf[3]*pred_x^2
        } else {
          pred_y <- cf[1]+cf[2+l] + 
            (cf[2]+cf[3 + (lf-1) + 2*(l-1)-1])*pred_x +
            (cf[3]+cf[3 + (lf-1) + 2*(l-1)])*pred_x^2
        }
        lines(pred_y~pred_x, col=lst_palev[[k]][l], lty=3)
      }
        
      # pvs
      pvs <- round(na.omit(pvs),2)
      mtext(paste('overall:rh', pvs[1], '\n+fact:rh', pvs[2], 'fact', pvs[3],
                  '\n*fact:rh', pvs[4], 'fact', pvs[5], 'inter', pvs[6]), line=0.5, cex=0.75)
    }
    leg(i)
    
    #---
    dev.off()
    
  }
  
}








