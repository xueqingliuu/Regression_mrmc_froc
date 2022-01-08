twopart_regression_sex <- function(dataset, nN, nD, nR, nI, nBoots) {
  
  #BIdata[which(BIdata == 0.5)] <- 1 ####This line is only for ordinal data!	
  
  #####Two part regression-Bootstrap
  
  alpha1 <- matrix(nrow = nBoots, ncol = 4)
  alpha2 <- matrix(nrow = nBoots, ncol = 4)
  alpha3 <- matrix(nrow = nBoots, ncol = 4)
  Sen.det <- matrix(nrow = nBoots, ncol = 6)
  FPR.det <- matrix(nrow = nBoots, ncol = 6)
  cAFAUC <- matrix(nrow = nBoots, ncol = 6)
  
  
  for (nB in 1:nBoots) {
    nN_sample <- sample(c(1:nN),size = nN, replace = TRUE)
    nd_sample <- sample(c(1:nD),size = nD, replace = TRUE)
    nd_new <- nd_sample + nN
    
    BIdata_boots <- c()
    tpdata_boots <- c()
    fpdata_boots <- c()
    dataset_nn <- c()
    nR_sample <-c()
    for (j in 1:nR) {
      index <- sample(c(1:nR),size = 1)
      dataset_j <- dataset[dataset[,5]==index,]
      dataset_nn <- dataset_j[nN_sample,]
      dataset_nd<-c()
      for (nd in 1:nD) {
        temp1 <- dataset_j[dataset_j[,6]==nd_new[nd],]
        dataset_nd <- rbind(dataset_nd,temp1)
      }
      temp2 <- rbind(dataset_nn, dataset_nd)
      
      BIdata <- gen.data.boots.sex(temp2, nN, nD, index)
      dataset_tp <- gen.data.tp.cov(temp2)
      dataset_fp <- gen.data.fp.cov(temp2)
      BIdata_boots <- rbind(BIdata_boots, BIdata)
      tpdata_boots <- rbind(tpdata_boots, dataset_tp)
      fpdata_boots <- rbind(fpdata_boots, dataset_fp)
    }
    
    #No.1 TP marks
    response1 <- as.numeric(c(tpdata_boots[,1], tpdata_boots[,2]))
    len1 <- length(tpdata_boots[,1])
    #sex <- as.factor(rep(tpdata_boots[,5],2))
    Z_i1 <- as.factor(c(rep(1, len1), rep(2, len1)))
    #Z_i1 <- matrix(c(rep(1, len1), rep(0, len1), rep(0, len1), rep(1, len1)), nrow = 2*len1, ncol = 2)
    mod.1 <- glm(response1~Z_i1+sex, family=binomial(link = "logit"))
    alpha11 <- mod.1$coefficients[1] #intercept
    alpha12 <- mod.1$coefficients[2] #z_i12
    alpha13 <- mod.1$coefficients[3] 
    alpha1[nB,] <- c(alpha11, alpha11+alpha12, alpha11+alpha13, alpha11+alpha13+alpha12)
    Sen.det[nB,1:4] <- sapply(alpha1[nB,], function(a){exp(a)/(1+exp(a))})
    
    
    #No.2 FP marks
    response2 <- as.numeric(c(fpdata_boots[,1], fpdata_boots[,2]))
    len2 <- length(fpdata_boots[,1])
    sex <- as.factor(rep(fpdata_boots[,5],2))
    Z_i2 <- as.factor(c(rep(1, len2), rep(2, len2)))
    mod.2 <- glm(response2 ~ Z_i2+sex, family = binomial(link = "logit"))
    alpha21 <- mod.2$coefficients[1]
    alpha22 <- mod.2$coefficients[2]
    alpha23 <- mod.2$coefficients[3]
    alpha2[nB,] <- c(alpha21, alpha21+alpha22, alpha21+alpha23, alpha21+alpha23+alpha22)
    FPR.det[nB,1:4] <- sapply(alpha2[nB,], function(a){exp(a)/(1+exp(a))})
    
    #No.3 Conditional AFROC-AUC
    # response3 <- as.numeric(c(BIdata_boots[which(BIdata_boots[,6]==1&BIdata_boots[,7]==1),1], BIdata_boots[which(BIdata_boots[,8]==1&BIdata_boots[,9]==1),2]))
    # response3[response3==0.5] <- 1
    # len3 <- length(BIdata_boots[which(BIdata_boots[,6]==1&BIdata_boots[,7]==1),1])
    # len4 <- length(BIdata_boots[which(BIdata_boots[,8]==1&BIdata_boots[,9]==1),2])
    # Z_i3 <- as.factor(c(rep(1, len3), rep(2, len4)))
    # sex_N <- as.factor(c(BIdata_boots[which(BIdata_boots[,6]==1&BIdata_boots[,7]==1),10], BIdata_boots[which(BIdata_boots[,8]==1&BIdata_boots[,9]==1),10]))
    # sex_D <- as.factor(c(BIdata_boots[which(BIdata_boots[,6]==1&BIdata_boots[,7]==1),11], BIdata_boots[which(BIdata_boots[,8]==1&BIdata_boots[,9]==1),11]))
    # mod.3 <- glm(response3 ~ Z_i3+sex_N+sex_D, family = binomial(link = "logit"))
    # alpha31 <- mod.3$coefficients[1]
    # alpha32 <- mod.3$coefficients[2]
    # alpha33 <- mod.3$coefficients[3]
    # alpha34 <- mod.3$coefficients[4]
    # alpha3[nB,] <- c(alpha31, alpha31+alpha32, alpha31+alpha33+alpha34, alpha31+alpha32+alpha33+alpha34)
    # cAFAUC[nB,1:4] <- sapply(alpha3[nB,], function(a){exp(a)/(1+exp(a))})
  }
  
  ##
  
  Sen.det[,5] <- Sen.det[,2] - Sen.det[,1]
  Sen.det[,6] <- Sen.det[,4] - Sen.det[,3]
  FPR.det[,5] <- FPR.det[,2] - FPR.det[,1]
  FPR.det[,6] <- FPR.det[,4] - FPR.det[,3]
  cAFAUC[,5] <- cAFAUC[,2] - cAFAUC[,1]
  cAFAUC[,6] <- cAFAUC[,4] - cAFAUC[,3]
  
  #Point Estimates
  Sen_all <- colMeans(Sen.det)
  FPR_all <- colMeans(FPR.det)
  cAFAUC_all <- colMeans(cAFAUC)
  
  var_alpha1_man <- diag(apply(alpha1, 2, var)[-c(3,4)])
  var_alpha1_woman <- diag(apply(alpha1, 2, var)[-c(1,2)])
  var_alpha2_man <- diag(apply(alpha2, 2, var)[-c(3,4)])
  var_alpha2_woman <- diag(apply(alpha2, 2, var)[-c(1,2)])
  var_alpha3_man <- diag(apply(alpha3, 2, var)[-c(3,4)])
  var_alpha3_woman <- diag(apply(alpha3, 2, var)[-c(1,2)])
  
  alpha1 <- colMeans(alpha1)
  alpha2 <- colMeans(alpha2)
  alpha3 <- colMeans(alpha3)
  
  ###SE
  sd_Sen <- apply(Sen.det, 2, sd)
  sd_FPR <- apply(FPR.det, 2, sd)
  sd_cAFAUC <- apply(cAFAUC, 2, sd)
  
  #CI
  CI_Sen_1_man <- c(Sen_all[1]-qnorm(0.975)*sd_Sen[1], Sen_all[1]+qnorm(0.975)*sd_Sen[1])
  CI_Sen_2_man <- c(Sen_all[2]-qnorm(0.975)*sd_Sen[2], Sen_all[2]+qnorm(0.975)*sd_Sen[2])
  CI_Sen_1_woman <- c(Sen_all[3]-qnorm(0.975)*sd_Sen[3], Sen_all[3]+qnorm(0.975)*sd_Sen[3])
  CI_Sen_2_woman <- c(Sen_all[4]-qnorm(0.975)*sd_Sen[4], Sen_all[4]+qnorm(0.975)*sd_Sen[4])
  CI_Sen_D_man <- c(Sen_all[5]-qnorm(0.975)*sd_Sen[5], Sen_all[5]+qnorm(0.975)*sd_Sen[5])
  CI_Sen_D_woman <- c(Sen_all[6]-qnorm(0.975)*sd_Sen[6], Sen_all[6]+qnorm(0.975)*sd_Sen[6])
  
  CI_FPR_1_man <- c(FPR_all[1]-qnorm(0.975)*sd_FPR[1], FPR_all[1]+qnorm(0.975)*sd_FPR[1])
  CI_FPR_2_man <- c(FPR_all[2]-qnorm(0.975)*sd_FPR[2], FPR_all[2]+qnorm(0.975)*sd_FPR[2])
  CI_FPR_1_woman <- c(FPR_all[3]-qnorm(0.975)*sd_FPR[3], FPR_all[3]+qnorm(0.975)*sd_FPR[3])
  CI_FPR_2_woman <- c(FPR_all[4]-qnorm(0.975)*sd_FPR[4], FPR_all[4]+qnorm(0.975)*sd_FPR[4])
  CI_FPR_D_man <- c(FPR_all[5]-qnorm(0.975)*sd_FPR[5], FPR_all[5]+qnorm(0.975)*sd_FPR[5])
  CI_FPR_D_woman <- c(FPR_all[6]-qnorm(0.975)*sd_FPR[6], FPR_all[6]+qnorm(0.975)*sd_FPR[6])
  
  CI_cAFAUC_1_man <- c(cAFAUC_all[1]-qnorm(0.975)*sd_cAFAUC[1], cAFAUC_all[1]+qnorm(0.975)*sd_cAFAUC[1])
  CI_cAFAUC_2_man <- c(cAFAUC_all[2]-qnorm(0.975)*sd_cAFAUC[2], cAFAUC_all[2]+qnorm(0.975)*sd_cAFAUC[2])
  CI_cAFAUC_1_woman <- c(cAFAUC_all[3]-qnorm(0.975)*sd_cAFAUC[3], cAFAUC_all[3]+qnorm(0.975)*sd_cAFAUC[3])
  CI_cAFAUC_2_woman <- c(cAFAUC_all[4]-qnorm(0.975)*sd_cAFAUC[4], cAFAUC_all[4]+qnorm(0.975)*sd_cAFAUC[4])
  CI_cAFAUC_D_man <- c(cAFAUC_all[5]-qnorm(0.975)*sd_cAFAUC[5], cAFAUC_all[5]+qnorm(0.975)*sd_cAFAUC[5])
  CI_cAFAUC_D_woman <- c(cAFAUC_all[6]-qnorm(0.975)*sd_cAFAUC[6], cAFAUC_all[6]+qnorm(0.975)*sd_cAFAUC[6])
  
  ####Hypothesis Testing 
  n_parameters <- 2
  L <- matrix(0, ncol = n_parameters, nrow = 1)
  for (i in c(1:(2-1))){
    L[i,1] = 1
    L[i,i+1] = -1
  }
  
  ###LDR
  res1_man <- WaldTest(L, alpha1[1:2], var_alpha1_man)
  res1_woman <- WaldTest(L, alpha1[3:4], var_alpha1_woman)
  ###NSMR
  res2_man <- WaldTest(L, alpha2[1:2], var_alpha2_man)
  res2_woman <- WaldTest(L, alpha2[1:2], var_alpha2_woman)
  ###cAFAUC
  res3_man <- WaldTest(L, alpha3[1:2], var_alpha3_man)
  res3_woman <- WaldTest(L, alpha3[3:4], var_alpha3_woman)
  
  # return a list
  return(list(Sen_all=Sen_all, FPR_all=FPR_all, cAFAUC=cAFAUC_all, 
              sd_Sen=sd_Sen, sd_FPR=sd_FPR, sd_cAFAUC=sd_cAFAUC,
              CI_Sen_1_man=CI_Sen_1_man, CI_Sen_1_woman=CI_Sen_1_woman,
              CI_Sen_2_man=CI_Sen_2_man, CI_Sen_2_woman=CI_Sen_2_woman, 
              CI_Sen_D_man=CI_Sen_D_man, CI_Sen_D_woman=CI_Sen_D_woman,
              CI_FPR_1_man=CI_FPR_1_man, CI_FPR_1_woman=CI_FPR_1_woman, 
              CI_FPR_2_man=CI_FPR_2_man, CI_FPR_2_woman=CI_FPR_2_woman, 
              CI_FPR_D_man=CI_FPR_D_man, CI_FPR_D_woman=CI_FPR_D_woman,
              CI_cAFAUC_1_man = CI_cAFAUC_1_man, CI_cAFAUC_1_woman = CI_cAFAUC_1_woman, 
              CI_cAFAUC_2_man = CI_cAFAUC_2_man, CI_cAFAUC_2_woman = CI_cAFAUC_2_woman, 
              CI_cAFAUC_D_man  = CI_cAFAUC_D_man, CI_cAFAUC_D_woman  = CI_cAFAUC_D_woman,
              p1.man=res1_man$mm.p, p1.woman=res1_woman$mm.p,
              p2.man=res2_man$mm.p, p2.woman=res2_woman$mm.p, 
              p3.man=res3_man$mm.p, p3.woman=res3_woman$mm.p))
}

