twopart_regression <- function(dataset, nN, nD, nR, nI, t) {
  BIdata <- gen.data(dataset, nN, nD, nR, t)
  dataset_tp <- gen.data.tp(dataset)
  dataset_fp <- gen.data.fp(dataset)
  #BIdata[which(BIdata == 0.5)] <- 1 ####This line is only for ordinal data!	
  
  #####Two part regression
  #No.1 TP marks
  response1 <- c(dataset_tp[,1], dataset_tp[,2])
  len1 <- length(dataset_tp[,1])
  Z_i1 <- matrix(c(rep(1, len1), rep(0, len1), rep(0, len1), rep(1, len1)), nrow = 2*len1, ncol = 2)  
  mod.1 <- glm(response1 ~ Z_i1-1, family = binomial(link = "logit"))
  alpha11 <- mod.1$coefficients[1]
  alpha12 <- mod.1$coefficients[2]
  alpha1 <- c(alpha11, alpha12)
  
  #No.2 FP marks
  response2 <- c(dataset_fp[,1], dataset_fp[,2])
  len2 <- length(dataset_fp[,1])
  Z_i2 <- matrix(c(rep(1, len2), rep(0, len2), rep(0, len2), rep(1, len2)), nrow = 2*len2, ncol = 2)  
  mod.2 <- glm(response2 ~ Z_i2-1, family = binomial(link = "logit"))
  alpha21 <- mod.2$coefficients[1]
  alpha22 <- mod.2$coefficients[2]
  alpha2 <- c(alpha21, alpha22)
  
  #No.3 Conditional AFROC-AUC
  response3 <- c(BIdata[which(BIdata[,6]==1&BIdata[,7]==1),1], BIdata[which(BIdata[,8]==1&BIdata[,9]==1),2])
  len3 <- length(BIdata[which(BIdata[,6]==1&BIdata[,7]==1),1])
  len4 <- length(BIdata[which(BIdata[,8]==1&BIdata[,9]==1),2])
  Z_i3 <- matrix(c(rep(1, len3), rep(0, len4), rep(0, len3), rep(1, len4)), nrow = (len3+len4), ncol = 2)
  mod.3 <- glm(response3 ~ Z_i3-1, family = binomial(link = "logit"))
  alpha31 <- mod.3$coefficients[1]
  alpha32 <- mod.3$coefficients[2]
  alpha3 <- c(alpha31, alpha32)
  
  #Point Estimates
  Sen.det <- sapply(alpha1, function(a){exp(a)/(1+exp(a))})
  FPR.det <- sapply(alpha2, function(a){exp(a)/(1+exp(a))})
  cAFAUC <- sapply(alpha3, function(a){exp(a)/(1+exp(a))})
  diff.Sen = Sen.det[2] - Sen.det[1]
  diff.FPR = FPR.det[2] - FPR.det[1]
  diff.AUC = cAFAUC[2] - cAFAUC[1]
  Sen_all <- c(Sen.det, diff.Sen)
  FPR_all <- c(FPR.det, diff.FPR)
  cAFAUC_all <- c(cAFAUC, diff.AUC)
  
  #Variance (Sparse correlation) 
  
  input_data1 <- data.frame(Phi = response1, Z_i = Z_i1,
                            i_lst = c(rep(1, nD*t*nR), rep(2, nD*t*nR)), j_lst = c(dataset_tp[,3], dataset_tp[,3]),
                            k_lst = c(dataset_tp[,4], dataset_tp[,4]))
  C_1 <- C_matrix_func(alpha1, input_data1)
  B_1 <- B_matrix1(alpha1, input_data1, nN, nD, nR, nI, t) 
  var_alpha1 <- solve(C_1)%*%B_1%*%t(solve(C_1))
  
  input_data2 <- data.frame(Phi = response2, Z_i = Z_i2,
                            i_lst = c(rep(1, nN*nR), rep(2, nN*nR)), j_lst = c(dataset_fp[,3], dataset_fp[,3]),
                            s_lst = c(dataset_fp[,4], dataset_fp[,4]))
  C_2 <- C_matrix_func(alpha2, input_data2)
  B_2 <- B_matrix2(alpha2, input_data2, nN, nR, nI) 
  var_alpha2 <- solve(C_2)%*%B_2%*%t(solve(C_2))
  
  input_data3 <- data.frame(Phi = response3, Z_i = Z_i3,
                            i_lst = c(rep(1, len3),rep(2, len4)), 
                            j_lst = c(BIdata[which(BIdata[,6]==1&BIdata[,7]==1),5], BIdata[which(BIdata[,8]==1&BIdata[,9]==1),5]), 
                            k_lst = c(BIdata[which(BIdata[,6]==1&BIdata[,7]==1),3], BIdata[which(BIdata[,8]==1&BIdata[,9]==1),3]), 
                            s_lst = c(BIdata[which(BIdata[,6]==1&BIdata[,7]==1),4], BIdata[which(BIdata[,8]==1&BIdata[,9]==1),4]))
  C_3 <- C_matrix_func(alpha3, input_data3)
  B_3 <- B_matrix3(alpha3, input_data3, nI, len3, len4) 
  var_alpha3 <- solve(C_3)%*%B_3%*%t(solve(C_3))
  
  ####Hypothesis Testing 
  n_parameters <- 2
  L <- matrix(0, ncol = n_parameters, nrow = 1)
  for (i in c(1:(2-1))){
    L[i,1] = 1
    L[i,i+1] = -1
  }
  res1 <- WaldTest(L, alpha1, var_alpha1)
  res2 <- WaldTest(L, alpha2, var_alpha2)
  res3 <- WaldTest(L, alpha3, var_alpha3)
  
  
  ####Variance of parameters of interest 
  gradient_Sen<- diag( mapply(function(a){exp(a)/(1+exp(a))**2}, alpha1) )
  var_Sen <- gradient_Sen%*%var_alpha1%*%gradient_Sen
  se_Sen <- c(sqrt(var_Sen[1,1]), sqrt(var_Sen[2,2]), sqrt(var_Sen[1,1]+var_Sen[2,2]-2*var_Sen[1,2]))
  
  gradient_FPR<- diag( mapply(function(a){exp(a)/(1+exp(a))**2}, alpha2) )
  var_FPR <- gradient_FPR%*%var_alpha2%*%gradient_FPR
  se_FPR <- c(sqrt(var_FPR[1,1]), sqrt(var_FPR[2,2]), sqrt(var_FPR[1,1]+var_FPR[2,2]-2*var_FPR[1,2]))
  
  gradient_cAFAUC <- diag( mapply(function(a){exp(a)/(1+exp(a))**2}, alpha3) )
  var_cAFAUC <- gradient_cAFAUC%*%var_alpha3%*%gradient_cAFAUC
  se_cAFAUC <- c(sqrt(var_cAFAUC[1,1]), sqrt(var_cAFAUC[2,2]), sqrt(var_cAFAUC[1,1]+var_cAFAUC[2,2]-2*var_cAFAUC[1,2]))
  
  ####Confidence interval
  CI_Sen_1 <- c(Sen_all[1]-qnorm(0.975)*se_Sen[1], Sen_all[1]+qnorm(0.975)*se_Sen[1])
  CI_Sen_2 <- c(Sen_all[2]-qnorm(0.975)*se_Sen[2], Sen_all[2]+qnorm(0.975)*se_Sen[2])
  CI_Sen_D <- c(Sen_all[3]-qnorm(0.975)*se_Sen[3], Sen_all[3]+qnorm(0.975)*se_Sen[3])
  
  CI_FPR_1 <- c(FPR_all[1]-qnorm(0.975)*se_FPR[1], FPR_all[1]+qnorm(0.975)*se_FPR[1])
  CI_FPR_2 <- c(FPR_all[2]-qnorm(0.975)*se_FPR[2], FPR_all[2]+qnorm(0.975)*se_FPR[2])
  CI_FPR_D <- c(FPR_all[3]-qnorm(0.975)*se_FPR[3], FPR_all[3]+qnorm(0.975)*se_FPR[3])
  
  CI_cAFAUC_1 <- c(cAFAUC_all[1]-qnorm(0.975)*se_cAFAUC[1], cAFAUC_all[1]+qnorm(0.975)*se_cAFAUC[1])
  CI_cAFAUC_2 <- c(cAFAUC_all[2]-qnorm(0.975)*se_cAFAUC[2], cAFAUC_all[2]+qnorm(0.975)*se_cAFAUC[2])
  CI_cAFAUC_D <- c(cAFAUC_all[3]-qnorm(0.975)*se_cAFAUC[3], cAFAUC_all[3]+qnorm(0.975)*se_cAFAUC[3])
  
  return(list(Sen=Sen_all, FPR=FPR_all, cAFAUC=cAFAUC_all, 
              se_Sen=se_Sen, se_FPR=se_FPR, se_cAFAUC=se_cAFAUC,
              CI_Sen_1 = CI_Sen_1, CI_Sen_2 = CI_Sen_2, CI_Sen_D = CI_Sen_D,
              CI_FPR_1 = CI_FPR_1, CI_FPR_2 = CI_FPR_2, CI_FPR_D = CI_FPR_D,
              CI_cAFAUC_1 = CI_cAFAUC_1, CI_cAFAUC_2 = CI_cAFAUC_2, CI_cAFAUC_D  = CI_cAFAUC_D,
              p1=res1$mm.p, p2=res2$mm.p, p3=res3$mm.p))
}

