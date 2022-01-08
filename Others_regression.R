#####Transform the dataset
gen.data.1 <- function(data, K1, K2){
  n.d <- length(data$lesions$perCase)  
  n.lesion <- sum(data$lesions$perCase)
  n.mod <- 2
  n.reader <- length(data$descriptions$readerID)
  
  y.diseased <- array(NA, dim=c(n.mod,n.reader,1)) 
  for (i in 1:n.d){
    add <- data$ratings$LL[,,i,which(data$lesions$IDs[i,]!=-Inf)]
    y.diseased <- abind(y.diseased,add)
  }
  y.diseased <- y.diseased[,,-1]
  
  y.non.diseased <- array(NA,dim=c(n.mod,n.reader,K1)) 
  y.non.diseased.0 <- data$ratings$NL[,,c(1:K1),]
  for (k in 1:K1){
    for (i in 1:n.mod){
      for (j in 1:n.reader){
        y.non.diseased[i,j,k] <- max(y.non.diseased.0[i,j,k,])
      }
    }
  }
  
  output <- matrix( nrow = (K1+K2)*n.reader, ncol = 6 )
  output[,4] <- rep(1:(K1+K2),n.reader)
  
  rr <- NULL
  for ( r in 1:n.reader ){
    rr <- c(rr,rep(r,K1+K2))
    count <- (r-1)*(K1+K2)
    output[(count+1):(count+K1+K2), 1] = c(y.non.diseased[1,r,], y.diseased[1,r,])
    output[(count+1):(count+K1+K2), 2] = c(y.non.diseased[2,r,], y.diseased[2,r,])
  }
  
  output[,5] <- rr
  
  dcase <- NULL
  for ( m in 1:n.d ) {
    temp1 <- data$lesions$perCase[m]
    dcase <- c(dcase, rep(m+K1, temp1))
  }
  
  output[,6] <- rep(c(1:K1,dcase),n.reader)
  kk <- c(rep(0, K1), rep(1, K2))
  output[,3] <- rep(kk, n.reader)
  return(output)
}

# gen.data.boots <- function(dataset, nN, nD, index, cov_r){
#   len1 <- length(dataset[dataset[,3]==1,1])
#   len2 <- length(dataset[dataset[,3]==0,1])
#   output<- matrix(nrow = len1*len2, ncol = 10)
#   output[1:(len1*len2), 5] <- index
#   output[1:(len1*len2), 10] <- cov_r
#     for (j in 1:len2){
#       count.2 <- len1*(j-1)
#       #case number k in 3rd column, s in 4th column:
#       output[(count.2+1):(count.2+ len1), 3] = j
#       #output[(count+count.2+1):(count + count.2+ num_tp), 4] = temp1
#       #fill in indicators
#       #Smth_phi: 1 /0.5 / 0
#       output[(count.2+1):(count.2+ len1), 1] = 1*(dataset[j, 1] < dataset[(len2 + 1):(len2 + len1), 1]) + .5*(dataset[j, 1] ==dataset[(len2 + 1):(len2 + len1), 1])
#       #Dir_phi
#       output[(count.2+1):(count.2+ len1), 2] = 1*(dataset[j, 2] < dataset[(len2 + 1):(len2 + len1), 2]) + .5*(dataset[j, 2] ==dataset[(len2 + 1):(len2 + len1), 2])
#       #Phi1 and Phi2 for i=1
#       output[(count.2+1):(count.2+ len1), 6] = 1*(dataset[,1]!=-Inf)[(len2+1) : (len2+len1)]
#       output[(count.2+1):(count.2+ len1), 7] = rep(1*(dataset[,1]!=-Inf)[j], len1)
#       #Phi1 and Phi2 for i=2
#       output[(count.2+1):(count.2+ len1), 8] = 1*(dataset[,2]!=-Inf)[(len2+1) : (len2+len1)]
#       output[(count.2+1):(count.2+ len1), 9] = rep(1*(dataset[,2]!=-Inf)[j], len1)
#       output[]
#     }
#   return(output)
# }

gen.data.boots.sex <- function(dataset, nN, nD, index){
  len1 <- length(dataset[dataset[,3]==1,1]) #number of lesions
  len2 <- length(dataset[dataset[,3]==0,1]) #number of normal subjects
  output<- matrix(nrow = len1*len2, ncol = 11)
  output[1:(len1*len2), 5] <- index
  for (j in 1:len2){
    count.2 <- len1*(j-1)
    #case number k in 3rd column, s in 4th column:
    output[(count.2+1):(count.2+ len1), 3] = j
    #output[(count+count.2+1):(count + count.2+ num_tp), 4] = temp1
    #fill in indicators
    #Smth_phi: 1 /0.5 / 0
    output[(count.2+1):(count.2+len1), 1] = 1*(dataset[j, 1] < dataset[(len2 + 1):(len2 + len1), 1]) + .5*(dataset[j, 1] ==dataset[(len2 + 1):(len2 + len1), 1])
    #Dir_phi
    output[(count.2+1):(count.2+len1), 2] = 1*(dataset[j, 2] < dataset[(len2 + 1):(len2 + len1), 2]) + .5*(dataset[j, 2] ==dataset[(len2 + 1):(len2 + len1), 2])
    #Phi1 and Phi2 for i=1
    output[(count.2+1):(count.2+len1), 6] = 1*(dataset[,1]!=-Inf)[(len2+1) : (len2+len1)]
    output[(count.2+1):(count.2+len1), 7] = rep(1*(dataset[,1]!=-Inf)[j], len1)
    #Phi1 and Phi2 for i=2
    output[(count.2+1):(count.2+len1), 8] = 1*(dataset[,2]!=-Inf)[(len2+1) : (len2+len1)]
    output[(count.2+1):(count.2+len1), 9] = rep(1*(dataset[,2]!=-Inf)[j], len1)
    #covariates:sex
    output[(count.2+1):(count.2+len1), 10] = rep(dataset[j,7], len1) #normal subject
    output[(count.2+1):(count.2+len1), 11] = dataset[(len2+1):(len2+len1),7] #diseased subjects
  }
  return(output)
}


gen.data.boots <- function(dataset, nN, nD, index){
  len1 <- length(dataset[dataset[,3]==1,1])
  len2 <- length(dataset[dataset[,3]==0,1])
  output<- matrix(nrow = len1*len2, ncol = 9)
  output[1:(len1*len2), 5] <- index
  for (j in 1:len2){
    count.2 <- len1*(j-1)
    #case number k in 3rd column, s in 4th column:
    output[(count.2+1):(count.2+ len1), 3] = j
    #output[(count+count.2+1):(count + count.2+ num_tp), 4] = temp1
    #fill in indicators
    #Smth_phi: 1 /0.5 / 0
    output[(count.2+1):(count.2+ len1), 1] = 1*(dataset[j, 1] < dataset[(len2 + 1):(len2 + len1), 1]) + .5*(dataset[j, 1] ==dataset[(len2 + 1):(len2 + len1), 1])
    #Dir_phi
    output[(count.2+1):(count.2+ len1), 2] = 1*(dataset[j, 2] < dataset[(len2 + 1):(len2 + len1), 2]) + .5*(dataset[j, 2] ==dataset[(len2 + 1):(len2 + len1), 2])
    #Phi1 and Phi2 for i=1
    output[(count.2+1):(count.2+ len1), 6] = 1*(dataset[,1]!=-Inf)[(len2+1) : (len2+len1)]
    output[(count.2+1):(count.2+ len1), 7] = rep(1*(dataset[,1]!=-Inf)[j], len1)
    #Phi1 and Phi2 for i=2
    output[(count.2+1):(count.2+ len1), 8] = 1*(dataset[,2]!=-Inf)[(len2+1) : (len2+len1)]
    output[(count.2+1):(count.2+ len1), 9] = rep(1*(dataset[,2]!=-Inf)[j], len1)
  }
  return(output)
}

gen.data <- function(dataset, nN, nD, nR, t){
  Phi1 <- 1*(dataset[,1] != -Inf)
  Phi2 <- 1*(dataset[,2] != -Inf)
  output<- matrix(nrow = nR*nN*nD*t, ncol = 9)
  #Make the data set
  #fill in reader, normal case, and diseased case numbers	
  for(i in 1:nR){ #go through the readers
    count = nN*nD*t*(i - 1) #for indexing
    subset <- dataset[dataset[,5]==i, ]  #data just for reader i
    num1 <- dim(subset)
    #reader number in 5th column:
    output[(count+1):(count+ nN*nD*t), 5] = i
    for (j in 1:nN){
      count.2<-(j-1)*nD*t
      #case number k in 3rd column, s in 4th column:
      output[(count+count.2+1):(count + count.2+ nD*t), 3] = j
      output[(count+count.2+1):(count + count.2+ nD*t), 4] = rep(c((nN+1):(nN+nD)),each = t)
      #fill in indicators
      #Smth_phi: 1 /0.5 / 0
      output[(count+count.2+1):(count + count.2+ nD*t), 1] = 1*(subset[j, 1] < subset[(nN + 1):(nN + nD*t), 1]) + .5*(subset[j, 1] ==subset[(nN + 1):(nN + nD*t), 1])
      #Dir_phi
      output[(count+count.2+1):(count + count.2+ nD*t), 2] = 1*(subset[j, 2] < subset[(nN + 1):(nN + nD*t), 2]) + .5*(subset[j, 2] ==subset[(nN + 1):(nN + nD*t), 2])
      #Phi1 and Phi2 for i=1
      output[(count+count.2+1):(count + count.2+ nD*t), 6] = Phi1[((i-1)*(nN+nD*t)+nN+1) : ((i-1)*(nN+nD*t)+nN+nD*t)]
      output[(count+count.2+1):(count + count.2+ nD*t), 7] = rep(Phi1[((i-1)*(nN+nD*t)+j)], nD*t)
      #Phi1 and Phi2 for i=2
      output[(count+count.2+1):(count + count.2+ nD*t), 8] = Phi2[((i-1)*(nN+nD*t)+nN+1) : ((i-1)*(nN+nD*t)+nN+nD*t)]
      output[(count+count.2+1):(count + count.2+ nD*t), 9] = rep(Phi2[((i-1)*(nN+nD*t)+j)], nD*t)
    }
  }
  return(output)
}

gen.data.varyt <- function(dataset, nN, nD, nR, L){
  len1 <- L
  len2 <- nN
  output<- matrix(nrow = nR*len1*len2, ncol = 9)
  for (i in 1:nR) {
    count = nN*L*(i - 1) #for indexing
    subset <- dataset[dataset[,5]==i, ]  #data just for reader i
    #num1 <- dim(subset)
    #reader number in 5th column:
    output[(count+1):(count+ nN*L), 5] = i
    for (j in 1:len2){
      count.2 <- len1*(j-1)
      #case number k in 3rd column, s in 4th column:
      output[(count+count.2+1):(count+count.2+ len1), 3] = j
      output[(count+count.2+1):(count+count.2+ len1), 4] = subset[subset[,3]==1,6]
      #fill in indicators
      #Smth_phi: 1 /0.5 / 0
      output[(count+count.2+1):(count+count.2+ len1), 1] = 1*(subset[j, 1] < subset[(len2 + 1):(len2 + len1), 1]) + .5*(subset[j, 1] == subset[(len2 + 1):(len2 + len1), 1])
      #Dir_phi
      output[(count+count.2+1):(count+count.2+ len1), 2] = 1*(subset[j, 2] < subset[(len2 + 1):(len2 + len1), 2]) + .5*(subset[j, 2] == subset[(len2 + 1):(len2 + len1), 2])
      #Phi1 and Phi2 for i=1
      output[(count+count.2+1):(count+count.2+ len1), 6] = 1*(subset[,1]!=-Inf)[(len2+1) : (len2+len1)]
      output[(count+count.2+1):(count+count.2+ len1), 7] = rep(1*(subset[,1]!=-Inf)[j], len1)
      #Phi1 and Phi2 for i=2
      output[(count+count.2+1):(count+count.2+ len1), 8] = 1*(subset[,2]!=-Inf)[(len2+1) : (len2+len1)]
      output[(count+count.2+1):(count+count.2+ len1), 9] = rep(1*(subset[,2]!=-Inf)[j], len1)
    }
  }
  return(output)
}


gen.data.tp <- function(dataset) {
  output <- dataset[dataset[,3]==1, c(1:2,5:6)]
  output[, 1] = 1*(output[, 1] != -Inf)
  output[, 2] = 1*(output[, 2] != -Inf)
  num1 <- length(output[,1])
  # output <- cbind(output, rep(cov_r, num1))
  return(output)
}

gen.data.tp.cov <- function(dataset) {
  output <- dataset[dataset[,3]==1, c(1:2,5:6,7)]
  output[, 1] = 1*(output[, 1] != -Inf)
  output[, 2] = 1*(output[, 2] != -Inf)
  num1 <- length(output[,1])
  # output <- cbind(output, rep(cov_r, num1))
  return(output)
}

gen.data.fp <- function(dataset) {
  output <- dataset[dataset[,3]==0, c(1:2,5:6)]
  output[, 1] = 1*(output[, 1] != -Inf)
  output[, 2] = 1*(output[, 2] != -Inf)
  num2 <- length(output[,1])
  # output <- cbind(output, rep(cov_r, num2))
  return(output)
}

gen.data.fp.cov <- function(dataset) {
  output <- dataset[dataset[,3]==0, c(1:2,5:6,7)]
  output[, 1] = 1*(output[, 1] != -Inf)
  output[, 2] = 1*(output[, 2] != -Inf)
  num2 <- length(output[,1])
  # output <- cbind(output, rep(cov_r, num2))
  return(output)
}
#---------calculate C matrix----------------#
# v1* v2^T [utils function (vec * vec^T to get an matrix)  where vec is a vector (dim:n*1)]
vecA_multiply_vecBinverse <- function(v1, v2){
  v1 <- unlist(v1); v2<-unlist(v2)
  n <- length(v1)
  out_array <- array(0,c(n,n))
  for (i in c(1:n)){
    out_array[i,] <- v1[i]*v2
  }
  return (out_array)
}

C_matrix_func <-function(alpha1, input_data1){
  W_df <- subset(input_data1, select = c(Z_i.1, Z_i.2))
  #calculate \alpha %*% W
  alpha_W <- apply(W_df, 1, (function(v)  sum(v*alpha1) ))
  #calculate -\frac{exp(alpha^T W)}{1+exp(alpha^T W)}^2
  scalar_value_lst <- sapply(alpha_W, (function(s) -exp(s)/((1+exp(s))^2) ))
  
  n <- length(W_df[1,])
  out_matrix <- array(0, c(n,n))
  
  uniqiue_w_vecs <- unique(W_df)
  unique_scalars <- unique(scalar_value_lst)
  
  if (length(uniqiue_w_vecs)==length(unique_scalars)){
    rep_time <- dim(W_df)[1]/length(uniqiue_w_vecs)
    for (t in 1:length(uniqiue_w_vecs)) {
      w_w_inverse_matrix <- vecA_multiply_vecBinverse(uniqiue_w_vecs[t,], uniqiue_w_vecs[t,]) 
      out_matrix <- out_matrix + unique_scalars[t] * w_w_inverse_matrix * rep_time
    }
  }else{
    #---initial format ---#
    for (i in c(1:dim(W_df)[1])){
      w_w_inverse_matrix <- vecA_multiply_vecBinverse(W_df[i,], W_df[i,])  #W_ijks*W_ijsk^T
      out_matrix <- out_matrix + scalar_value_lst[i] *w_w_inverse_matrix
    }
  }
  return (out_matrix)
}



#------- calculate B_matrix -----------#
# utils function : calculate U(beta) given a row_vector 
U_func <- function(hat_beta, w_vec, tag){
  if (tag){ 
    out_vec <- 1/(1+exp( sum(hat_beta*w_vec) )) * w_vec
  }else{
    out_vec <- -exp(sum(  hat_beta *w_vec )) /(1+exp( sum(hat_beta*w_vec) )) * w_vec
  }
  return (out_vec)
}


# calculate B_matrix for Phi1
B_matrix1 <- function(alpha1, input_data1, nN, nD, nR, nI, L){
  num <- nR*nI*L
  W <- subset(input_data1, select = c(-i_lst, -j_lst, -k_lst, -Phi))
  W <- matrix(as.numeric(unlist(W)), nrow = num, ncol = nI)
  
  condition <- input_data1$Phi==1   # the indicator variable
  U_matrix <- matrix(nrow = num, ncol = nI)
  for (i in 1:num) {
    U_matrix[i,] <- U_func(alpha1, W[i,], condition[i])
  }
  
  #assumption I  S_ijks = {(i',j',k',s'):k=k' or s=s'}
  out_matrix <- matrix(0, nrow = nI, ncol = nI)
  array_sum_for_ijks <- matrix(nrow = nD, ncol = nI)
  for (k in 1:nD) {
      k_new = k+nN
      vec_sum_for_u_ijks <- colSums(U_matrix[which(input_data1$k_lst==k_new),] )  
      array_sum_for_ijks[k,] <- vec_sum_for_u_ijks
  }
  
  for (ijks in 1:num){
    k <- input_data1$k_lst[ijks]
    u_ijks <- U_matrix[ijks,]
    #the vector sum \sum_{i'j'k's'}
    vec_sum_for_u_ijks <- array_sum_for_ijks[k-nN, ]
    out_matrix <- out_matrix + vecA_multiply_vecBinverse(u_ijks, vec_sum_for_u_ijks)
  }
  return (out_matrix)
}

# calculate B_matrix for Phi2
B_matrix2 <- function(alpha1, input_data1, nN, nR, nI){
  num <- nR*nI*nN
  W <- subset(input_data1, select = c(-i_lst, -j_lst, -s_lst, -Phi))
  W <- matrix(as.numeric(unlist(W)), nrow = num, ncol = nI)
  
  condition <- input_data1$Phi==1   # the indicator variable
  U_matrix <- matrix(nrow = num, ncol = nI)
  for (i in 1:num) {
    U_matrix[i,] <- U_func(alpha1, W[i,], condition[i])
  }
  
  #assumption I  S_ijks = {(i',j',k',s'):k=k' or s=s'}
  out_matrix <- matrix(0, nrow = nI, ncol = nI)
  array_sum_for_ijks <- matrix(nrow = nN, ncol = nI)
  for (s in 1:nN){
      vec_sum_for_u_ijks <- colSums(U_matrix[which(input_data1$s_lst==s),] )  
      array_sum_for_ijks[s,] <- vec_sum_for_u_ijks
  }
  
  for (ijks in 1:num){
    s <- input_data1$s_lst[ijks]
    u_ijks <- U_matrix[ijks,]
    vec_sum_for_u_ijks <- array_sum_for_ijks[s, ]
    out_matrix <- out_matrix + vecA_multiply_vecBinverse(u_ijks, vec_sum_for_u_ijks)
  }
  return (out_matrix)
}

# calculate B_matrix for Phi3
B_matrix3 <- function(alpha1, input_data3, nI, len3, len4){
  num <- len3 + len4
  W <- subset(input_data3, select = c(-i_lst, -j_lst, -k_lst, -s_lst, -Phi))
  W <- matrix(as.numeric(unlist(W)), nrow = num, ncol = nI)
  
  condition <- input_data3$Phi==1   # the indicator variable
  U_matrix <- matrix(nrow = num, ncol = nI)
  for (i in 1:num) {
    U_matrix[i,] <- U_func(alpha1, W[i,], condition[i])
  }
  
  #assumption I  S_ijks = {(i',j',k',s'):k=k' or s=s'}
  new_klst <- unique(input_data3$k_lst)
  new_slst <- unique(input_data3$s_lst)
  new_nN <- length(new_klst)
  new_nD <- length(new_slst)
  out_matrix <- matrix(0, nrow = nI, ncol = nI)
  array_sum_for_ijks <- array(dim=c(new_nD, new_nN, nI))
  for (s in 1:new_nD){
    s_new <- new_slst[s]
    for (k in 1:new_nN){
      k_new <- new_klst[k]
      vec_sum_for_u_ijks <- colSums(U_matrix[which((input_data3$s_lst==s_new)|(input_data3$k_lst==k_new)),,drop=FALSE] )  
      array_sum_for_ijks[s,k,] <- vec_sum_for_u_ijks
    }
  }
  
  for (ijks in 1:num){
    k <- input_data3$k_lst[ijks]
    s <- input_data3$s_lst[ijks]
    u_ijks <- U_matrix[ijks,]
    k_new <- which(new_klst==k)
    s_new <- which(new_slst==s)
    #the vector sum \sum_{i'j'k's'}
    vec_sum_for_u_ijks <- array_sum_for_ijks[s_new, k_new, ]
    out_matrix <- out_matrix + vecA_multiply_vecBinverse(u_ijks, vec_sum_for_u_ijks)
  }
  return (out_matrix)
}


WaldTest = function(L,thetahat,Vn, h=0) # H0: L theta = h
{
  r = dim(L)[1]
  W = t(L%*%thetahat-h) %*% solve( L%*% Vn %*%t(L) ) %*% (L%*%thetahat-h)  #m被消掉了
  W = as.numeric(W)
  pval = 1-pchisq(W,r)
  return(list(W=W,r=r,mm.p=pval))
}
