## Monochorionic vs. Dichorionic
## 单绒毛膜双胎(同卵) vs. 双绒毛膜双胎(异卵)
## functions
library(irr) # icc
library(OpenMx)
library(lavaan)
# Chisq test
twin.chisqTest <- function(dat,family_var,target_var,zygo_var,
                           Monochorionic_value=1,
                           Dichorionic_value=0){
  t1 <- table(target=dat[,target_var],zygo=dat[,zygo_var],dat[,family_var])
  t2 <- rbind(table(apply(t1,3,function(x)x['1',as.character(Monochorionic_value)])),
              table(apply(t1,3,function(x)x['1',as.character(Dichorionic_value)])))
  rownames(t2) <- c('Monochorionic','Dichorionic')
  colnames(t2) <- c('Neither','One','Both')
  t2 <- as.data.frame.matrix(t2)
  t2$Total <- rowSums(t2)
  t2$Chisq <- chisq.test(t2)$statistic
  t2$P <- chisq.test(t2)$p.value
  return(t2)
}
# ICC model: 组内相关系数（intraclass correlation, ICC）
twin.ICCTest <- function(dat,family_var,target_var,zygo_var,
                         Monochorionic_value=1,
                         Dichorionic_value=0,n_bootstrap=1000){
  w1 <- which(dat[,zygo_var]==Monochorionic_value)
  t1 <- aggregate(dat[w1,c(target_var)],list(dat[w1,family_var]),c)
  w2 <- which(dat[,zygo_var]==Dichorionic_value)
  t2 <- aggregate(dat[w2,c(target_var)],list(dat[w2,family_var]),c)
  r1 <- icc(t1$x); v1 <- sprintf('%.2f(%.2f,%.2f)',r1$value,r1$lbound,r1$ubound)
  r2 <- icc(t2$x); v2 <- sprintf('%.2f(%.2f,%.2f)',r2$value,r2$lbound,r2$ubound)
  # Fisher's z-transformation：将ICC值通过Fisher's z-transformation转换成正态分布的变量，然后用t检验或方差分析（ANOVA）来比较两组ICC值之间的差异。
  icc1_z <- atanh(r1$value) # Fisher's z-transformation
  icc2_z <- atanh(r2$value) # Fisher's z-transformation
  # 计算ICC值之间的差异
  diff_z <- icc1_z - icc2_z
  # Bootstrap方法计算差异的置信区间和p值
  diff_z_bootstrap <- numeric(n_bootstrap)
  for (i in 1:n_bootstrap) {
    twin1_boot <- t1$x[sample(1:nrow(t1$x), replace = TRUE),]
    twin2_boot <- t2$x[sample(1:nrow(t1$x), replace = TRUE),]
    icc1_boot <- icc(twin1_boot)$value
    icc2_boot <- icc(twin2_boot)$value
    diff_z_bootstrap[i] <- atanh(icc1_boot) - atanh(icc2_boot)
  }
  se_diff <- sqrt(var(diff_z_bootstrap))
  p_value <- 2 * min(sum(diff_z_bootstrap <= -diff_z), sum(diff_z_bootstrap >= diff_z))/n_bootstrap # 双侧检验
  return(list(ICC_Monochorionic=v1,ICC_Dichorionic=v2,P=p_value,
              ori_Monochorionic=r1,ori_Dichorionic=r2))
}

##
## 结构方程模型
# functions to get para for each model: ACE/ADE
# 使用函数mxExpectationNormal()来计算潜变量的方差的期望值
# expectation <- mxExpectationNormal(covariance = fit$output$cov, means = fit$output$means, dimnames = names(fit$output$means))
get_sta <- function(modelx,model=c('ACE')){
  sumACE <- summary(modelx)
  if(model=='ACE') all_var <- c('a','c','e')
  if(model=='ADE') all_var <- c('a','d','e')
  # Proportion of variance explained; standardised
  V_conf <- modelx$output$confidenceIntervals
  V_ACE_stand <- V_conf[1:3,2] #c(coef(modelx)[all_var])^2
  #V_ACE_stand <- V_ACE/sum(V_ACE,na.rm=T)
  V_ACE_conf_stand <- cbind('2.5%'=V_conf[1:3,1],'97.5%'=V_conf[1:3,3])
  names(V_ACE_stand) <- rownames(V_ACE_conf_stand) <- all_var
  #-2 log likelihood of ACE model
  LL_ACE <- as.numeric(-logLik(modelx)*2)
  output_ACE <- list(coef=V_ACE_stand,CI=V_ACE_conf_stand,
                     `-2LL`=LL_ACE,'df'=as.numeric(sumACE$degreesOfFreedom),
                     'AIC'=as.numeric(sumACE$AIC),model=modelx)
  if(model=='ACE'){
    output_ACE_dis <- c(sprintf('A(95%sCI):%.3f(%.3f-%.3f)','%',
                                V_ACE_stand['a'],V_ACE_conf_stand['a','2.5%'],V_ACE_conf_stand['a','97.5%']),
                        sprintf('C(95%sCI):%.3f(%.3f-%.3f)','%',
                                V_ACE_stand['c'],V_ACE_conf_stand['c','2.5%'],V_ACE_conf_stand['c','97.5%']),
                        sprintf('E(95%sCI):%.3f(%.3f-%.3f)','%',
                                V_ACE_stand['e'],V_ACE_conf_stand['e','2.5%'],V_ACE_conf_stand['e','97.5%']),
                        sprintf('-2LL:%.3f',LL_ACE),sprintf('df:%d',sumACE$degreesOfFreedom),
                        sprintf('AIC:%.3f',sumACE$AIC)
    )
    
  }else{
    output_ACE_dis <- c(sprintf('A(95%sCI):%.3f(%.3f-%.3f)','%',
                                V_ACE_stand['a'],V_ACE_conf_stand['a','2.5%'],V_ACE_conf_stand['a','97.5%']),
                        sprintf('D(95%sCI):%.3f(%.3f-%.3f)','%',
                                V_ACE_stand['d'],V_ACE_conf_stand['d','2.5%'],V_ACE_conf_stand['d','97.5%']),
                        sprintf('E(95%sCI):%.3f(%.3f-%.3f)','%',
                                V_ACE_stand['e'],V_ACE_conf_stand['e','2.5%'],V_ACE_conf_stand['e','97.5%']),
                        sprintf('-2LL:%.3f',LL_ACE),sprintf('df:%d',sumACE$degreesOfFreedom),
                        sprintf('AIC:%.3f',sumACE$AIC)
    )
  }
  return(list(value=output_ACE,display=output_ACE_dis))
}
# ACE model
twin.ACEModel.OpenMx <- function(dat,family_var,target_var,zygo_var,
                          #confounding_var,
                          Monochorionic_value=1,
                          Dichorionic_value=0){
  confounding_var <- NULL
  # extract value
  w1 <- which(dat[,zygo_var]==Monochorionic_value)
  t1 <- aggregate(dat[w1,c(target_var,confounding_var)],list(dat[w1,family_var]),c)
  mzData <- do.call(cbind,t1); 
  colnames(mzData) <- c(family_var,apply(expand.grid(c(1,2),c(target_var,confounding_var)),1,function(x)sprintf('%s_%s',x[2],x[1])))
  w2 <- which(dat[,zygo_var]==Dichorionic_value)
  t2 <- aggregate(dat[w2,c(target_var,confounding_var)],list(dat[w2,family_var]),c)
  dzData <- do.call(cbind,t2); 
  colnames(dzData) <- c(family_var,apply(expand.grid(c(1,2),c(target_var,confounding_var)),1,function(x)sprintf('%s_%s',x[2],x[1])))
  ## set starting value
  target_var_2 <- sprintf('%s_%s',target_var,c(1,2))
  mzData <- as.data.frame(mzData);dzData <- as.data.frame(dzData)
  use_mzData <- as.data.frame(mzData[,target_var_2])
  use_dzData <- as.data.frame(dzData[,target_var_2])
  MeanStartValue <- mean(colMeans(rbind(use_mzData,use_dzData),na.rm = T), na.rm = T)
  varMz <- cov(use_mzData,use="complete")
  varDz <- cov(use_dzData,use="complete")
  VarStartValue <- sqrt((mean(c(varMz[1],varMz[4],varDz[1],varDz[4]), na.rm = T))/3)
  aceVars <- c("A1","C1","E1","A2","C2","E2")
  selVars <- target_var_2
  # create path
  latVariances <- mxPath(from=aceVars, arrows=2, free=FALSE, values=1)
  obsMeans <- mxPath(from="one", to=selVars, arrows=1, free=TRUE, values = MeanStartValue, labels="mean")
  latMeans <- mxPath(from="one", to=aceVars, arrows=1, free=FALSE, values=0)
  pathAceT1 <- mxPath(from=c("A1","C1","E1"), to=selVars[1], arrows=1, free=TRUE, values = VarStartValue, label=c("a","c","e"))
  pathAceT2 <- mxPath(from=c("A2","C2","E2"), to=selVars[2], arrows=1, free=TRUE, values = VarStartValue, label=c("a","c","e"))
  # ACE model: A-A: rMZ=1,rDZ=0.5; C-C: rMZ=rDZ=1
  # ADE model: A-A: rMZ=1,rDZ=0.5; C-C: rMZ=1; rDZ=0.25
  covA1A2_MZ <- mxPath(from="A1", to="A2", arrows=2, free=FALSE, values=1)
  covA1A2_DZ <- mxPath(from="A1", to="A2", arrows=2, free=FALSE, values=0.5)
  covC1C2 <- mxPath(from="C1", to="C2", arrows=2, free=FALSE, values=1)
  #####
  allvars <- colnames(mzData)[-1]
  # combine all path
  paths <- list(latVariances, latMeans, obsMeans, pathAceT1, pathAceT2, covC1C2)
  modelMZ <- mxModel(model="MZ", type="RAM", manifestVars=selVars,latentVars=aceVars, paths, covA1A2_MZ, mxData(observed=use_mzData, type="raw"),
                     mxCI(c('a','c','e')))
  modelDZ <- mxModel(model="DZ", type="RAM", manifestVars=selVars, latentVars=aceVars, paths, covA1A2_DZ, mxData(observed=use_dzData, type="raw"),
                     mxCI(c('a','c','e')))
  #####
  #create the -2 loglikelihood function
  minus2ll <- mxAlgebra(expression=MZ.fitfunction + DZ.fitfunction,
                        name="minus2loglikelihood")
  #specified -2 loglikelihood and combining Dz and Mz models
  obj <- mxFitFunctionAlgebra("minus2loglikelihood")
  
  #--------------ACE model--------------
  covA      <- mxAlgebra( expression=a %*% t(a)/(a %*% t(a)+c %*% t(c)+e %*% t(e)), name="covA" )
  covC      <- mxAlgebra( expression=c %*% t(c)/(a %*% t(a)+c %*% t(c)+e %*% t(e)), name="covC" ) 
  covE      <- mxAlgebra( expression=e %*% t(e)/(a %*% t(a)+c %*% t(c)+e %*% t(e)), name="covE" )
  modelACE <- mxModel(model="ACE", modelMZ, modelDZ, minus2ll,obj,
                      covA,covC,covE,mxCI(c('covA','covC','covE')))
  fitACE <- mxRun(modelACE,intervals=T);

  #--------------AE model--------------
  modelAE <- mxModel(modelACE, name = "AE")
  modelAE <- omxSetParameters(modelAE,labels = "c", free = FALSE,values = 0)
  fitAE <- mxRun(modelAE,intervals = T);
  
  #--------------CE model--------------
  modelCE <- mxModel(modelACE, name = "CE")
  modelCE <- omxSetParameters(modelCE,labels = "a", free = FALSE,values = 0)
  fitCE <- mxRun(modelCE,intervals = T);
  
  #--------------E model--------------
  modelE <- mxModel(modelACE, name = "E")
  modelE <- omxSetParameters(modelE,labels = "a", free = FALSE,values = 0)
  modelE <- omxSetParameters(modelE,labels = "c", free = FALSE,values = 0)
  fitE <- mxRun(modelE,intervals = T);
  ## combine
  all_model <- list(ACE=fitACE,AE=fitAE,CE=fitCE,E=fitE)
  all_model_res <- lapply(all_model,function(x)get_sta(x,model='ACE'))
  output_res <- do.call(rbind,lapply(all_model_res,function(x)x$display))
  output_res_df <- as.data.frame(gsub('(.*):(.*)','\\2',output_res))
  colnames(output_res_df) <- gsub('(.*):(.*)','\\1',output_res)[1,]
  output_res_df['AE','P.Value'] <- anova(all_model$AE,all_model$ACE)[2,'p']
  output_res_df['CE','P.Value'] <- anova(all_model$CE,all_model$ACE)[2,'p']
  output_res_df['E','P.Value'] <- anova(all_model$E,all_model$ACE)[2,'p']
  print(output_res_df)
  return(list(display=output_res_df,all_model=all_model,all_model_res=all_model_res))
}
# ADE model
twin.ADEModel.OpenMx <- function(dat,family_var,target_var,zygo_var,
                          #confounding_var,
                          Monochorionic_value=1,
                          Dichorionic_value=0){
  # extract value
  confounding_var <- NULL
  w1 <- which(dat[,zygo_var]==Monochorionic_value)
  t1 <- aggregate(dat[w1,c(target_var,confounding_var)],list(dat[w1,family_var]),c)
  mzData <- do.call(cbind,t1); 
  colnames(mzData) <- c(family_var,apply(expand.grid(c(1,2),c(target_var,confounding_var)),1,function(x)sprintf('%s_%s',x[2],x[1])))
  w2 <- which(dat[,zygo_var]==Dichorionic_value)
  t2 <- aggregate(dat[w2,c(target_var,confounding_var)],list(dat[w2,family_var]),c)
  dzData <- do.call(cbind,t2); 
  colnames(dzData) <- c(family_var,apply(expand.grid(c(1,2),c(target_var,confounding_var)),1,function(x)sprintf('%s_%s',x[2],x[1])))
  ## set starting value
  target_var_2 <- sprintf('%s_%s',target_var,c(1,2))
  use_mzData <- as.data.frame(mzData[,target_var_2])
  use_dzData <- as.data.frame(dzData[,target_var_2])
  MeanStartValue <- mean(colMeans(rbind(use_mzData,use_dzData),na.rm = T), na.rm = T)
  varMz <- cov(use_mzData,use="complete")
  varDz <- cov(use_dzData,use="complete")
  VarStartValue <- sqrt((mean(c(varMz[1],varMz[4],varDz[1],varDz[4]), na.rm = T))/3)
  aceVars <- c("A1","D1","E1","A2","D2","E2")
  selVars <- target_var_2
  # create path
  latVariances <- mxPath(from=aceVars, arrows=2, free=FALSE, values=1)
  obsMeans <- mxPath(from="one", to=selVars, arrows=1, free=TRUE, values = MeanStartValue, labels="mean")
  latMeans <- mxPath(from="one", to=aceVars, arrows=1, free=FALSE, values=0)
  pathAceT1 <- mxPath(from=c("A1","D1","E1"), to=selVars[1], arrows=1, free=TRUE, values = VarStartValue, label=c("a","d","e"))
  pathAceT2 <- mxPath(from=c("A2","D2","E2"), to=selVars[2], arrows=1, free=TRUE, values = VarStartValue, label=c("a","d","e"))
  # ADE model: A-A: rMZ=1,rDZ=0.5; D-D: rMZ=1; rDZ=0.25
  covA1A2_MZ <- mxPath(from="A1", to="A2", arrows=2, free=FALSE, values=1)
  covA1A2_DZ <- mxPath(from="A1", to="A2", arrows=2, free=FALSE, values=0.5)
  covD1D2_MZ <- mxPath(from="D1", to="D2", arrows=2, free=FALSE, values=1)
  covD1D2_DZ <- mxPath(from="D1", to="D2", arrows=2, free=FALSE, values=0.25)
  # combine all path
  paths <- list(latVariances, latMeans, obsMeans, pathAceT1, pathAceT2)
  modelMZ <- mxModel(model="MZ", type="RAM", manifestVars=selVars,latentVars=aceVars, paths, covA1A2_MZ,covD1D2_MZ,mxData(observed=use_mzData, type="raw"))
  modelDZ <- mxModel(model="DZ", type="RAM", manifestVars=selVars, latentVars=aceVars, paths, covA1A2_DZ,covD1D2_DZ,mxData(observed=use_dzData, type="raw"))
  #create the -2 loglikelihood function
  minus2ll <- mxAlgebra(expression=MZ.fitfunction + DZ.fitfunction,
                        name="minus2loglikelihood")
  #specified -2 loglikelihood and combining Dz and Mz models
  obj <- mxFitFunctionAlgebra("minus2loglikelihood")
  
  #--------------ADE model--------------
  covA      <- mxAlgebra( expression=a %*% t(a)/(a %*% t(a)+d %*% t(d)+e %*% t(e)), name="covA" )
  covD      <- mxAlgebra( expression=d %*% t(d)/(a %*% t(a)+d %*% t(d)+e %*% t(e)), name="covD" ) 
  covE      <- mxAlgebra( expression=e %*% t(e)/(a %*% t(a)+d %*% t(d)+e %*% t(e)), name="covE" )
  modelADE <- mxModel(model="ADE", modelMZ, modelDZ, minus2ll,obj,
                      covA,covD,covE,mxCI(c('covA','covD','covE')))
  fitADE <- mxRun(modelADE,intervals = T);
  
  #--------------AE model--------------
  modelAE <- mxModel(modelADE, name = "AE")
  modelAE <- omxSetParameters(modelAE,labels = "d", free = FALSE,values = 0)
  fitAE <- mxRun(modelAE,intervals = T);
  
  #--------------DE model--------------
  modelDE <- mxModel(modelADE, name = "DE")
  modelDE <- omxSetParameters(modelDE,labels = "a", free = FALSE,values = 0)
  fitDE <- mxRun(modelDE,intervals = T);
  
  #--------------E model--------------
  modelE <- mxModel(modelADE, name = "E")
  modelE <- omxSetParameters(modelE,labels = "a", free = FALSE,values = 0)
  modelE <- omxSetParameters(modelE,labels = "d", free = FALSE,values = 0)
  fitE <- mxRun(modelE,intervals = T);
  ## combine
  all_model <- list(ADE=fitADE,AE=fitAE,DE=fitDE,E=fitE)
  all_model_res <- lapply(all_model,function(x)get_sta(x,model='ADE'))
  output_res <- do.call(rbind,lapply(all_model_res,function(x)x$display))
  output_res_df <- as.data.frame(gsub('(.*):(.*)','\\2',output_res))
  colnames(output_res_df) <- gsub('(.*):(.*)','\\1',output_res)[1,]
  output_res_df['AE','P.Value'] <- anova(all_model$AE,all_model$ADE)[2,'p']
  output_res_df['DE','P.Value'] <- anova(all_model$DE,all_model$ADE)[2,'p']
  output_res_df['E','P.Value'] <- anova(all_model$E,all_model$ADE)[2,'p']
  print(output_res_df)
  return(list(display=output_res_df,all_model=all_model,all_model_res=all_model_res))
}

