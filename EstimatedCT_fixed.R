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

############################### default arguments (optional)
ref_file <- gsub(x = args[grep(x = args, pattern = "ref_file=")], pattern = "ref_file=", replacement = "")

if (length(ref_file) == 0) {
  alpha_hat <- NULL
} else {
  alpha_hat <- get(load(ref_file))  
  dim(alpha_hat)
  #[1] 200 5
  head(alpha_hat)
}

N = nrow(mu_w)
#see dimensions
dim(bulkM)
#[1]   100 200   bulk sample x gene
dim(mu_w)
#[1] 100   3



############################### main functions
# matrix inverse with Choleski decomposition
sinv <- function(A){
  return(chol2inv(chol((A + t(A)) / 2)))
}



CTS <- function(x, mu_w, #V_w,
                alpha_hat = NULL, v_alpha = 0.5,
                Sigma_hat = NULL, v_Sigma = 50,
                c1=NULL, c2=NULL, 
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
                            list(a_i = a_i)
                          } 
      
    }
    res_list <- update_aw_parallel()
    
    # Combine all parallel outputs into the final array
    for (i in 1:N) {
      a[i, , , kk + 1] <- res_list[[i]]$a_i
      w[i, , kk + 1] <- mu_w[i,]
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
  }

  return(list(alpha = alpha, Sigma = Sigma, a = a, w = w, beta = beta, gamma = gamma, sigma2 = sigma2))
}


############################### Run
res <- CTS(bulkM, mu_w, 
           alpha_hat = alpha_hat, v_alpha = 0.5,
           Sigma_hat = NULL, v_Sigma = 50,
           # c1, c2,
           iter = 3000)
a = res$a
w = res$w


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


oc = list(CTS2 = CTS2, w_iter = w)
save(oc, file  = paste0(output_file,'.RData'))
