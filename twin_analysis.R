source('twin_analysis_function.R')
###############################################################
## pipelines
# load demo data
dat <- read.delim('data.txt')
family_var <- 'FAMID';
target_var <- 'BPD';
zygo_var <- 'zygo'; # 1-Monochorionic，0-Dichorionic
r1 <- twin.chisqTest(dat,family_var,target_var,zygo_var)
r2 <- twin.ICCTest(dat,family_var,target_var,zygo_var)
r3 <- twin.ACEModel.OpenMx(dat,family_var,target_var,zygo_var)
r4 <- twin.ADEModel.OpenMx(dat,family_var,target_var,zygo_var)
