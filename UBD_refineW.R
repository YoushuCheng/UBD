library(MASS)
library(data.table)
library(LaplacesDemon)
library(tmvtnorm)
library(tmvnsim)
library(DirichletReg)
library(BayesPrism)
library(dplyr)
options(stringsAsFactors=F)
args = commandArgs(trailingOnly=TRUE)

bulkM_file = args[1]
bp.res_file = args[2]
output_file  = args[3]

############################### read data
bulkM <- get(load(bulkM_file))
bp.res <- get(load(bp.res_file))  


#get the estimated proportions
mu_w <- get.fraction (bp=bp.res,
                      which.theta="final",
                      state.or.type="type")
head(mu_w)
table(mu_w<=0)

#get the uncertainty for the estimated proportions
V_w2 <- sapply(colnames(mu_w), function(name) {
  # 1. Get the expression data for that cell type
  exp_data <- get.exp(bp = bp.res,
                      state.or.type = "type",
                      cell.name = name)
  # 2. Calculate the row sums
  rowSums(exp_data)
})
colnames(V_w2) = colnames(mu_w)

############################### default arguments (optional)
ref_file <- gsub(x = args[grep(x = args, pattern = "ref_file=")], pattern = "ref_file=", replacement = "")
scale_MH <- gsub(x = args[grep(x = args, pattern = "scale_MH=")], pattern = "scale_MH=", replacement = "")

if (length(ref_file) == 0) {
  alpha_hat <- NULL
} else {
  alpha_hat <- get(load(ref_file))  
  
  dim(alpha_hat)
  #[1] 200 5
  head(alpha_hat)
}
scale_MH = ifelse(length(scale_MH)==0, 200, as.numeric(scale_MH))


N = nrow(mu_w)
#see dimensions
dim(bulkM)
#[1]   100 200   bulk sample x gene
dim(mu_w)
#[1] 100   3
dim(V_w2)
#[1] 100   3

############################### determine scale_prior
#function to calculate the standard deviation for a Dirichlet distribution
D_sd <- function(alpha) {
  # Calculate the sum of the parameters
  alpha_0 <- sum(alpha)
  # Calculate the variance of each component
  variance <- alpha * (alpha_0 - alpha) / (alpha_0^2 * (alpha_0 + 1))
  # Return standard deviation
  return(sqrt(variance))
}

objective_function <- function(scalar, V_w2_matrix, target_sd) {
  # Scale the input matrix
  V_w2_scaled <- scalar * V_w2_matrix
  sd_matrix <- t(apply(V_w2_scaled, 1, D_sd))
  current_mean_sd <- mean(colMeans(sd_matrix))
  
  # Return the difference from the target.
  return(current_mean_sd - target_sd)
}

optim_result <- uniroot(f = objective_function,interval = c(1e-6, 1),V_w2_matrix = V_w2,target_sd = 0.22)
scale_prior <- optim_result$root
cat("scale_prior found:", scale_prior, "\n")



############################### main functions
# matrix inverse with Choleski decomposition
sinv <- function(A){
  return(chol2inv(chol((A + t(A)) / 2)))
}


#h function in MH, log(liklihood) + log(dirichlet prior)
log_h_N<-function(weight, P, a, sigma2, x, mu_w, V_w, i, kk, beta, gamma, c1, c2){
  #print(kk)
  log_lik = 0
  for (j in 1:P){
    this_mu = t(weight) %*% a[i, j, , kk + 1]
    if (hasArg("c1")) {this_mu <- this_mu + crossprod(c1[i, ], beta[j, , kk])}
    if (hasArg("c2")) {this_mu <- this_mu + t(weight) %*% gamma[j, , , kk] %*% c2[i, ]}
    
    this_log_lik = dnorm(x[i, j], mean=this_mu, sd=sqrt(sigma2[j, kk]), log=TRUE)
    log_lik = log_lik + this_log_lik
  }
  #sum up
  ans<-log_lik + DirichletReg::ddirichlet(matrix(weight, nrow = 1), alpha = scale_prior*as.vector(V_w[i, ]), log=TRUE)
  return(ans)
}



CTS <- function(x, mu_w, V_w,
                alpha_hat = NULL, v_alpha = 0.5,
                Sigma_hat = NULL, v_Sigma = 50,
                c1, c2, 
                mh_delta = 10, 
                iter = 1000){
  # sample N, probe P, cell type K.
  # dimension of covariates T1 & T2.
  
  
  # x: bulk expression / methylation. N * P.
  # mu_w: true unobnserved cell type proportion. (by BayesPrism) N * K.
  # V_w: uncertainty of the cell type proportion. (by BayesPrism, alpha of Dirichlet) N * K.
  # alpha_hat: (optional) average CTS expression obtained from scRNA-seq data. P * K.
  # v_alpha: (optional) uncertainty of CTS expression. scalar.
  # Sigma_hat: (optional) cell type covariance matrix estimated from scRNA-seq data. P * K * K.
  # v_Sigma: (optional) parameter in Inv-Wishart prior. scalar.
  # c1: covariates affecting bulk expression / methylation. N * T1.
  # c2: covariates affecting CTS expression / methylation. N * T2.
  # mh_delta: scale factor for the MH process
  # iter: iteration time for Gibbs sampling.
  
  N = nrow(x)
  P = ncol(x)
  K = ncol(mu_w)
  if (hasArg("c1")) {T1 = ncol(c1)}
  if (hasArg("c2")) {T2 = ncol(c2)}
  
  #update alpha_hat and Sigma_hat if not provided
  X_bM = t(x)
  #P * N dimension
  if(is.null(alpha_hat)) alpha_hat = matrix(1, nrow(X_bM), K) * apply(X_bM, 1, mean)
  if(is.null(Sigma_hat)) {
    Sigma_hat = array(NA, dim = c(nrow(X_bM), K, K))
    for(i in 1:nrow(X_bM)) Sigma_hat[i,,] = diag(K) * var(X_bM[i,]) / sum(colMeans(mu_w)^2)
  }
  
  if(is.null(rownames(alpha_hat))) rownames(alpha_hat) = rownames(X_bM)
  if(is.null(rownames(Sigma_hat))) rownames(Sigma_hat) = rownames(X_bM)
  if(is.null(colnames(alpha_hat))) colnames(alpha_hat) = colnames(mu_w)
  if(is.null(colnames(Sigma_hat))) colnames(Sigma_hat) = colnames(mu_w)
  
  print(head(alpha_hat))
  print(Sigma_hat[1,,])
  
  
  # Seven parameters in the model.
  a = array(NA, dim = c(N, P, K, iter))
  w = array(NA, dim = c(N, K, iter))
  if (hasArg("c1")) {beta = array(NA, dim = c(P, T1, iter))}
  if (hasArg("c2")) {gamma = array(NA, dim = c(P, K, T2, iter))}
  sigma2 = array(NA, dim = c(P, iter))
  alpha = array(NA, dim = c(P, K, iter))
  Sigma = array(NA, dim = c(P, K, K, iter))
  accept_w = rep(1, N)
  
  # Initialization
  
  # alpha_j & Sigma_j
  
  for (j in 1:P){
    alpha[j, , 1] = mvrnorm(1, alpha_hat[j, ], v_alpha * diag(rep(1, K)))
    Sigma[j, , , 1] = rinvwishart(v_Sigma, Sigma_hat[j, , ])
  }
  
  # a_ij
  for (j in 1:P){
    a[, j, , 1] <- mvrnorm(n = N, mu = alpha[j, , 1], Sigma = Sigma[j, , , 1])
  }
  
  # w_i
  #w[, , 1] <- mvrnorm(n = N, mu = mu_w[i, ], Sigma = V_w[i, , ])
  w[, , 1] <- mu_w
  
  # beta_j
  if (hasArg("c1")) {
    beta[, , 1] <- matrix(0, P, T1) #mvrnorm(n = P, mu = rep(0, T1), Sigma = diag(rep(1, T1)))
  }
  # gamma_jk
  if (hasArg("c2")) {
    for (k in 1:K){
      gamma[, k, , 1] <- matrix(0, P, T2) #mvrnorm(n = P, mu = rep(0, T2), Sigma = diag(rep(1, T2)))
    }
  }
  
  # sigma^2_j
  sigma2[, 1] <- runif(P, 0.5, 1)
  
  ##### All lack one sampling dimension. #######
  
  
  kk = 1
  # Gibbs sampling
  for (kk in 1:(iter-1)){
    
    # update a_ij
    for (i in 1:N){
      for (j in 1:P){
        V <- sinv(tcrossprod(w[i, , kk]) / sigma2[j, kk] + sinv(Sigma[j, , , kk]))
        
        t1 <- 0
        t2 <- 0
        if (hasArg("c1")) { t1 <- crossprod(c1[i, ], beta[j, , kk])}
        if (hasArg("c2")) {t(w[i, , kk]) %*% gamma[j, , , kk] %*% c2[i, ]}
        
        u <- c(x[i, j] - t1 - t2) / sigma2[j, kk] * w[i, , kk] + sinv(Sigma[j, , , kk]) %*% alpha[j, , kk]
        a[i, j, , kk + 1] <- mvrnorm(n = 1, mu = V %*% u, Sigma = V)
      }
    }
    
    # Update w_i
    for(i in 1:N){
      #propose a w using dirichlet 
      w_proposed <- DirichletReg::rdirichlet(n=1, mh_delta*as.vector(w[i, , kk]))
      w[i, , kk + 1] <- w[i, , kk]
      #avoid 0 in w_proposed
      epsilon <- 1e-20
      w_proposed[w_proposed  == 0] <- epsilon
      # Normalize to ensure that the sum is 1
      w_proposed <- w_proposed / sum(w_proposed)
      
      top_h <- log_h_N(as.vector(w_proposed), P, a, sigma2, x, mu_w, V_w, i, kk)
      bottom_h <- log_h_N(as.vector(w[i, , kk]), P, a, sigma2, x, mu_w, V_w, i, kk)
      top_q <- DirichletReg::ddirichlet(matrix(w[i, , kk], nrow = 1), mh_delta*as.vector(w_proposed), log = T)
      bottom_q <- DirichletReg::ddirichlet(w_proposed, mh_delta*as.vector(w[i, , kk]), log = T)
      ratio <- exp(top_h + top_q - bottom_h - bottom_q)
      
      if(ratio >= runif(n=1, min=0, max=1)){
        w[i, , kk + 1] <- w_proposed
        accept_w[i] <- accept_w[i] + 1
      }
      
      #print(i)
    }

    
    
    # update alpha_j
    for(j in 1:P){
      V <- sinv(N * sinv(Sigma[j, , , kk]) + diag(rep(1, K)) / v_alpha)
      u <- (alpha_hat[j, ] / v_alpha + sinv(Sigma[j, , , kk]) %*% colSums(a[, j, , kk + 1]))
      alpha[j, , kk + 1] <- mvrnorm(n = 1, mu = V %*% u, Sigma = V)
    }
    
    # Update Sigma_j
    for(j in 1:P){
      V <- Sigma_hat[j, , ] + tcrossprod(t(a[, j, , kk + 1]) - alpha[j, , kk + 1])
      Sigma[j, , , kk + 1] <- rinvwishart(v_Sigma + N, V)
    }
    

    
    # lm
    for (j in 1:P){
      y <- x[, j] - rowSums(w[, , kk + 1] * a[, j, , kk + 1])
      if (hasArg("c1")) { X <- c1[, ]}
      if (hasArg("c2")) {
        for (k in 1:K){
          X <- cbind(X, c2[, ] * w[, k, kk + 1])
        }
      }
      
      if ((!hasArg("c1")) & (!hasArg("c2"))) {
        sigma2[j, kk + 1] <- rinvwishart(N, 1 + sum(y^2))
      } else{
        f1 <- lm(y ~ X + 0)
        u <- f1$coefficients
        V <- vcov(f1) / sigma(f1)^2 * sigma2[j, kk]
        temp <- mvrnorm(n = 1, mu = u, Sigma = V)
        
        if (hasArg("c1")) {
          beta[j, , kk + 1] <- temp[1:T1]
          if (hasArg("c2")) {
            for (k in 1:K){
              gamma[j, k, , kk + 1] <- temp[(T1 + k * T2 + 1 - T2):(T1 + k * T2)]
            }
          }
        } else{
          for (k in 1:K){
            gamma[j, k, , kk + 1] <- temp[(k * T2 + 1 - T2):(k * T2)]
          }
        }
        
        nu = 1 + sigma(f1)^2 * (N - length(u))
        sigma2[j, kk + 1] <- rinvwishart(N, nu)
      }
    }
    print(kk)
    print(accept_w)
  }

  return(list(alpha = alpha, Sigma = Sigma, a = a, w = w, beta = beta, gamma = gamma, sigma2 = sigma2, accept_w = accept_w))
}


############################### Run
res <- CTS(bulkM, mu_w, V_w2,
           alpha_hat = alpha_hat, v_alpha = 0.5,
           Sigma_hat = NULL, v_Sigma = 50,
           # c1, c2,
           mh_delta = scale_MH,
           iter = 3000)
a = res$a
w = res$w
accept_w = res$accept_w

dim(w)
#[1]  100    3 3000
dim(a)
#[1]  100  200    3 3000
quantile(res$accept_w/3000)


CTS2 <- apply(a[, , ,1000:3000], c(1, 2, 3), mean)
dim(CTS2)
#[1]  100   200         3
#  sample  gene        CT  


oc = list(CTS2 = CTS2, w_iter = w,
          accept_w = accept_w)
save(oc, file  = paste0(output_file,'.RData'))


