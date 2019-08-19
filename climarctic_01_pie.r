#####
# climarctic preparation
#####

rm(list=ls())
gc()

require(plotrix)

# dir loads ####
dir_out  <- 'Projets/Climarctic/stats/out_janv/'
dir_save <- paste0(dir_out, 'saves/') 
dir_pie  <- paste0(dir_out, '01_pie/')
dir.create(dir_pie, showWarnings=F)

source('bin/src/my_prog/R/pie_taxo.r')

#---
file <- paste0(dir_save, 'lst_comm.Rdata')
load(file)

# pie ####

selec_smp <- list(Knud  =which(env_sort$site     == 'Knudsenheia'),
                  Ossian=which(env_sort$site     == 'Ossian'),
                  dry   =which(env_sort$moisture == 'dry'),
                  medium=which(env_sort$moisture == 'medium'),
                  wet   =which(env_sort$moisture == 'wet'),
                  top   =which(env_sort$depth    == 'T'),
                  deep  =which(env_sort$depth    == 'D'))

pdf(paste0(dir_pie, 'pie1.pdf'), width=15, height=7)

pie <- pie_taxo(mr_sort, taxo_sort, 1:5, selec_smp, mat_lay=matrix(c(0,1,2,8, 3:5,8, 0,6,7,8), nrow=3, byrow=T),
                wdt_lay=c(1,1,1,3), hei_lay=c(rep(1.1, 3)), last_tax_text=F)

dev.off()

#---
selec_smp2 <- factor(apply(env_sort[,c("moist_in_site",'depth')], 1, function(x) paste(x, collapse='_')))

pdf(paste0(dir_pie, 'pie2.pdf'), width=15, height=7)

pie <- pie_taxo(mr_sort, taxo_sort, 1:5, selec_smp2, mat_lay=matrix(c(1:4,13, 5:8,13, 9:12,13), nrow=3, byrow=T),
                wdt_lay=c(1,1,1,1, 3), hei_lay=c(rep(1.1, 3)), last_tax_text=F)

dev.off()

#























