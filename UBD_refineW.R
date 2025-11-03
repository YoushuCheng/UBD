args = commandArgs(trailingOnly=TRUE)
library(MASS)
library(data.table)
library(LaplacesDemon)
library(DirichletReg)
library(BayesPrism)
library(dplyr)
options(stringsAsFactors=F)
library(doRNG) 
library(progressr)
handlers(global = TRUE)
handlers("txtprogressbar")
library(doFuture)
registerDoFuture()

bulkM_file = args[1]
bp.res_file = args[2]
output_file  = args[3]
num_threads = as.numeric(args[4])


#reference document about doFuture: https://cran.r-project.org/web/packages/doFuture/vignettes/doFuture-2-dopar.html
plan(multicore, workers = num_threads)
options(future.globals.maxSize = 1000 * 1024^3)  # 1000 GB


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
target_coverage <- gsub(x = args[grep(x = args, pattern = "target_coverage=")], pattern = "target_coverage=", replacement = "")

if (length(ref_file) == 0) {
  alpha_hat <- NULL
} else {
  alpha_hat <- get(load(ref_file))  
  dim(alpha_hat)
  #[1] 200 5
  head(alpha_hat)
}
scale_MH = ifelse(length(scale_MH)==0, 200, as.numeric(scale_MH))
target_coverage = ifelse(length(target_coverage)==0, 0.9, as.numeric(target_coverage))

N = nrow(mu_w)
#see dimensions
dim(bulkM)
#[1]   100 200   bulk sample x gene
dim(mu_w)
#[1] 100   3
dim(V_w2)
#[1] 100   3

############################### determine scale_prior
#function to calculate the width of each 99% CI for a Dirichlet distribution
Width_CI <- function(alpha) {
  # Calculate the sum of the parameters
  alpha_0 <- sum(alpha)
  # Using the Beta approximation for each component
  # The marginal distribution of component i is Beta(alpha_i, sum(alphas) - alpha_i)
  ci_lower <- qbeta(0.005, shape1 = alpha, shape2 = alpha_0 - alpha)
  ci_upper <- qbeta(0.995, shape1 = alpha, shape2 = alpha_0 - alpha)
  # Calculate the width of each CI
  ci_widths <- ci_upper - ci_lower
  
  return(ci_widths)
}

objective_function <- function(scalar, V_w2_matrix, target_coverage) {
  # Scale the input matrix
  V_w2_scaled <- scalar * V_w2_matrix
  width_matrix <- t(apply(V_w2_scaled, 1, Width_CI))
  current_mean_width <- mean(colMeans(width_matrix))
  
  # Return the difference from the target.
  return(current_mean_width - target_coverage)
}

optim_result <- uniroot(f = objective_function,interval = c(1e-6, 1),V_w2_matrix = V_w2,target_coverage = target_coverage)
scale_prior <- optim_result$root
cat("scale_prior found:", scale_prior, "\n")



############################### main functions
# matrix inverse with Choleski decomposition
sinv <- function(A){
  return(chol2inv(chol((A + t(A)) / 2)))
}


#h function in MH, log(liklihood) + log(dirichlet prior)
log_h_N<-function(weight, P, a_i, sigma2, x, mu_w, V_w, i, kk, beta, gamma, c1=NULL, c2=NULL){
  #print(kk)
  log_lik = 0
  for (j in 1:P){
    this_mu = t(weight) %*% a_i[j, ]
    if (!is.null(c1)) {this_mu <- this_mu + crossprod(c1[i, ], beta[j, , kk])}
    if (!is.null(c2)) {this_mu <- this_mu + t(weight) %*% gamma[j, , , kk] %*% c2[i, ]}
    
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
                c1=NULL, c2=NULL, 
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
  if (!is.null(c1)) {T1 = ncol(c1)}
  if (!is.null(c2)) {T2 = ncol(c2)}
  
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
  if (!is.null(c1)) {beta = array(NA, dim = c(P, T1, iter))}
  if (!is.null(c2)) {gamma = array(NA, dim = c(P, K, T2, iter))}
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
  if (!is.null(c1)) {
    beta[, , 1] <- matrix(0, P, T1) #mvrnorm(n = P, mu = rep(0, T1), Sigma = diag(rep(1, T1)))
  }
  # gamma_jk
  if (!is.null(c2)) {
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
    
    ######### start of the parallel here #########
    update_aw_parallel <- function() {
      
      prgrs <- progressor(steps = N)
      res_list <- foreach(i = 1:N,
                          .options.future = list(seed = F),
                          .packages = c("DirichletReg")) %dorng% {
                            # call progressor inside each future
                            prgrs()
                            
                            # update a_i
                            # Each worker returns one slice a[i, , , kk + 1]
                            a_i <- array(NA, dim = c(P, K))  # store local results
                            
                            for (j in 1:P) {
                              V <- sinv(tcrossprod(w[i, , kk]) / sigma2[j, kk] + sinv(Sigma[j, , , kk]))
                              
                              t1 <- 0
                              t2 <- 0
                              if (!is.null(c1)) t1 <- crossprod(c1[i, ], beta[j, , kk])
                              if (!is.null(c2)) t2 <- t(w[i, , kk]) %*% gamma[j, , , kk] %*% c2[i, ]
                              
                              u <- c(x[i, j] - t1 - t2) / sigma2[j, kk] * w[i, , kk] + sinv(Sigma[j, , , kk]) %*% alpha[j, , kk]
                              
                              a_i[j, ] <- mvrnorm(n = 1, mu = V %*% u, Sigma = V)
                            }
                            
                            # update w_i
                            w_current <- as.vector(w[i, , kk])
                            w_new <- w_current
                            accept <- 0L
                            
                            # propose a w using dirichlet
                            w_proposed <- DirichletReg::rdirichlet(n = 1, mh_delta * w_current)
                            w_proposed <- as.vector(w_proposed)
                            
                            # avoid zeros in w_proposed
                            eps <- 1e-20
                            w_proposed[w_proposed == 0] <- eps
                            # Normalize to ensure that the sum is 1
                            w_proposed <- w_proposed / sum(w_proposed)
                            
                            
                            top_h    <- log_h_N(w_proposed, P, a_i, sigma2, x, mu_w, V_w, i, kk)
                            bottom_h <- log_h_N(w_current,  P, a_i, sigma2, x, mu_w, V_w, i, kk)
                            top_q    <- DirichletReg::ddirichlet(matrix(w_current, nrow = 1),
                                                                 mh_delta * w_proposed, log = TRUE)
                            bottom_q <- DirichletReg::ddirichlet(matrix(w_proposed, nrow = 1),
                                                                 mh_delta * w_current, log = TRUE)
                            
                            ratio <- exp(top_h + top_q - bottom_h - bottom_q)
                            if (ratio >= runif(n=1, min=0, max=1)) {
                              w_new <- w_proposed
                              accept <- 1L
                            }
                            
                            list(a_i = a_i, w = w_new, accept = accept)
                          } 
      
    }
    res_list <- update_aw_parallel()
    
    # Combine all parallel outputs into the final array
    for (i in 1:N) {
      a[i, , , kk + 1] <- res_list[[i]]$a_i
      w[i, , kk + 1] <- res_list[[i]]$w
      accept_w[i]    <- accept_w[i] + res_list[[i]]$accept
    }
    ################ end of the parallel  #########
    
    
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
      if (!is.null(c1)) { X <- c1[, ]}
      if (!is.null(c2)) {
        for (k in 1:K){
          X <- cbind(X, c2[, ] * w[, k, kk + 1])
        }
      }
      
      if ((!!is.null(c1)) & (!!is.null(c2))) {
        sigma2[j, kk + 1] <- rinvwishart(N, 1 + sum(y^2))
      } else{
        f1 <- lm(y ~ X + 0)
        u <- f1$coefficients
        V <- vcov(f1) / sigma(f1)^2 * sigma2[j, kk]
        temp <- mvrnorm(n = 1, mu = u, Sigma = V)
        
        if (!is.null(c1)) {
          beta[j, , kk + 1] <- temp[1:T1]
          if (!is.null(c2)) {
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


idx = 1000:3000
#CTS1 <- apply(a[, , ,idx], c(1, 2, 3), mean) 
dims_out <- dim(a)[1:3]
#Reshape 'a' into a 2D matrix where the 4th dimension becomes columns
a_matrix <- matrix(a, nrow = prod(dims_out), ncol = dim(a)[4])
mean_vec <- rowMeans(a_matrix[, idx])
#Reshape the resulting vector back into a 3D array
CTS2 <- array(mean_vec, dim = dims_out)

# all.equal(CTS1, CTS2)


dim(CTS2)
#  sample  gene        CT  


oc = list(CTS2 = CTS2, w_iter = w,
          accept_w = accept_w)
save(oc, file  = paste0(output_file,'.RData'))


