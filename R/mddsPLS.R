#' The core function of the Multi-Data-Driven sparse PLS function
#'
#' This function should not be used directly by the user.
#'
#' @param Xs A data-frame of a matrix or a list of data-frames or matrices of \emph{n} rows each, the number of individuals. Some rows must be missing. The different matrices can have different numbers of columns. The length of Xs is denoted by \emph{K}.
#' @param Y A matrix of n rows of a vector of length n detailing the response matrix. No missing values are allowed in that matrix.
#' @param lambda A real \eqn{[0,1]} where 1 means just perfect correlations will be used and 0 no regularization is used.
#' @param R A strictly positive integer detailing the number of components to build in the model.
#' @param mode A character chain. Possibilities are "\emph{reg}", which implies regression problem or anything else which means clustering is considered. Default is "\emph{reg}".
#' @param verbose Logical. If TRUE, the function cats specificities about the model. Default is FALSE.
#'
#' @return A list containing the following objects:
#' \describe{
#'   \item{u}{A list of length \emph{K}. Each element is a \emph{p_kXR} matrix : the
#'    weights per block per axis.}
#'   \item{v}{A \emph{qXR} matrix : the weights for the \emph{Y} part.}
#'   \item{ts}{A list of length \emph{R}. Each element is a \emph{nXK} matrix : the
#'    scores per axis per block.}
#'   \item{(t,s)}{Two \emph{nXR} matrices, scores of the \emph{X} and \emph{Y} parts.}
#'   \item{(t_frak,s_frak)}{Two \emph{nXR} matrices, final scores of the \emph{X} and \emph{Y} part.
#'    They correspond to \emph{PLS} scores of \emph{(t,s)} scores and so \emph{t_frak^T s_frak} is diagonal,
#'    \emph{t_frak}, respectively \emph{s_frak}, carries the same information as \emph{t}, respectively \emph{s}.}
#'   \item{B}{A list of length \emph{K}. Each element is a \emph{p_kXq} matrix : the
#'    regression matrix per block.}
#'   \item{(mu_x_s,sd_x_s)}{Two lists of length \emph{K}. Each element is a \emph{p_k} vector : the
#'    mean and standard deviation variables per block.}
#'   \item{(mu_y,sd_y)}{Two vectors of length \emph{q} : the mean and the standard deviation variables for \emph{Y} part.}
#'   \item{R}{Given as an input.}
#'   \item{q}{A non negatvie integer : the number of variables of \emph{Y} matrix. }
#'   \item{Ms}{A list of length \emph{K}. Each element is a \emph{qXp_k} matrix : the
#'    soft-thresholded empirical variance-covariance matrix \eqn{Y^TX_k/(n-1)}.}
#'   \item{lambda}{Given as an input.}
#' }
MddsPLS_core <- function(Xs,Y,lambda=0,R=1,mode="reg",verbose=FALSE){
  is.multi <- is.list(Xs)&!(is.data.frame(Xs))
  if(!is.multi){
    Xs <- list(Xs)
  }
  K <- length(Xs)
  ps <- lapply(Xs,ncol)
  ## Standardize Xs
  mu_x_s <- lapply(Xs,colMeans)
  sd_x_s <- lapply(Xs,function(X){apply(X,2,stats::sd)})
  Xs <- lapply(Xs,scale)
  pos_0 <- lapply(sd_x_s,function(sdi){which(sdi==0)})
  if(length(unlist(pos_0))>0){
    for(i_0 in which(lapply(pos_0,function(pp){length(pp)})>0)){
      Xs[[i_0]][,pos_0[[i_0]]] <- 0
    }
  }
  ## Standardize Y
  Y_0 <- Y
  if(mode=="reg"){
    if(!(is.matrix(Y)|is.data.frame(Y))){
      Y <- as.matrix(Y)
    }
  }
  else{
    Y_df <- data.frame(Y)
    Y <- scale(stats:: model.matrix( ~ Y - 1, data=Y_df))
  }
  mu_y <- colMeans(Y)
  sd_y <- apply(Y,2,stats::sd)
  for(q_j in 1:length(sd_y)){
    if(sd_y[q_j]!=0){
      Y[,q_j] <- scale(Y[,q_j])
    }
  }
  q <- ncol(Y)
  n <- nrow(Y)
  ## Create soft-thresholded matrices
  lambda_in <- lambda
  if(length(lambda_in)==1){
    lambda_in <- rep(lambda_in,K)
  }
  Ms <- lapply(1:K,function(k,Xs,Y,l,n){
    M0 <- crossprod(Y,Xs[[k]])/(n-1)
    M <- abs(M0) - l[k]
    M[which(M<0)] <- 0
    M <- sign(M0)*M
  },Xs,Y,lambda_in,n)
  if(verbose){
    N_max <- sum(unlist(lapply(Ms,function(m){length(which(colSums(abs(m))!=0))})))
    cat(paste("At most ",N_max," variable(s) can be selected",sep=""));cat("\n")
  }
  ## Solve optimization problem
  #### Inside problems
  u_t_r = u_t_r_0 <- list()
  t_r <- list()
  z_r <- list()
  # BETA_r <- list()
  for(k in 1:K){
    if(norm(Ms[[k]])==0){
      svd_k <- list(v=matrix(0,
                             nrow = ncol(Ms[[k]]),
                             ncol = R))
    }
    else{
      svd_k <- svd(Ms[[k]],nu = 0,nv = R)
    }
    u_t_r[[k]] = u_t_r_0[[k]] <- svd_k$v
    if(k==1){
      for(r in 1:R){
        t_r[[r]] <- matrix(NA,n,K)
        z_r[[r]] <- matrix(NA,q,K)
      }
    }
    for(r in 1:R){
      t_r[[r]][,k] <- Xs[[k]]%*%u_t_r[[k]][,r]
      z_r[[r]][,k] <- Ms[[k]]%*%u_t_r[[k]][,r]
    }
  }
  t <- matrix(NA,n,R)
  v <- matrix(0,q,R)

  z_all <- do.call(cbind,z_r)
  svd_all <- svd(z_all,nu=R,nv=R)
  u <- svd_all$v
  v0 <- svd_all$u
  v <- v0
  s <- Y%*%v0
  t_all <- do.call(cbind,t_r)
  t <-  t_all%*%u
  svd_frak <- svd(crossprod(t,s),nu = R,nv = R)
  u_frak <- svd_frak$u
  v_frak <- svd_frak$v
  t_frak <- t%*%u_frak
  s_frak <- s%*%v_frak
  alphas <- rep(0,R)
  for(r in 1:R){
    n_t_2<-sum(diag(crossprod(t_frak[,r])))#t[,r])))
    if(n_t_2!=0){
      alphas[r] <- sum(diag(crossprod(s_frak[,r],t_frak[,r])))/n_t_2#sum(diag(crossprod(s[,r],t[,r])))/n_t_2
    }else{
      alphas[r] <- 0
    }
  }
  if(mode=="reg"){
    B <- list()
    for(k in 1:K){
      beta_k <- u[(k-1)*R+1:R,,drop=FALSE]
      B[[k]] <- u_t_r[[k]]%*%beta_k%*%u_frak
      for(r in 1:R){
        B[[k]][,r] <- B[[k]][,r]*alphas[r]
      }
      B[[k]]  <- tcrossprod(B[[k]],v%*%v_frak)
    }
  }
  else{
    dataf <- data.frame(cbind(Y_0,t));colnames(dataf)[1]<-"Y"
    for( cc in 2:ncol(dataf)){
      dataf[,cc] <- as.numeric(levels(dataf[,cc])[dataf[,cc]])
    }
    sds <- apply(dataf[,-1,drop=FALSE],2,stats::sd)
    if(any(sds==0)){
      pos_sd0 <- as.numeric(which(sds==0))
      if(length(pos_sd0)==length(sds)){
        B <- NULL
      }else{
        dataf <- dataf[,-c(1+pos_sd0)]
        B <- MASS::lda(Y ~ ., data = dataf)
        B <- list(B=B,sds=sds)
      }
    }
    else{
      B <- MASS::lda(Y ~ ., data = dataf)
    }
  }
  if(verbose){
    a<-lapply(u_t_r,function(u){apply(u,2,function(u){length(which(abs(u)>1e-9))})})
    cat("    For each block of X, are selected in order of component:");cat("\n")
    for(k in 1:K){
      cat(paste("        @ (",paste(a[[k]],collapse = ","),") variable(s)",sep=""));cat("\n")
    }
    cat("    For the Y block, are selected in order of component:");cat("\n")
    cat(paste("        @ (",paste(apply(v,2,function(u){length(which(abs(u)>1e-9))}),
                                  collapse = ","),") variable(s)",sep=""));cat("\n")
  }
  list(u=u_t_r,v=v,ts=t_r,beta_comb=u,t=t,s=s,t_frak=t_frak,s_frak=s_frak,B=B,mu_x_s=mu_x_s,sd_x_s=sd_x_s,mu_y=mu_y,
       sd_y=sd_y,R=R,q=q,Ms=Ms,lambda=lambda)
}


#' Multi-Data-Driven sparse PLS function
#'
#' This function takes a set \eqn{X} of \eqn{K} matrices defining the same \eqn{n} individuals and a matrix \eqn{Y} defining also those individuals. According to the num-
#' ber of components \eqn{R}, the user fixes the number of components the model
#' must be built on. The coefficient lambda regularizes the quality of proximity to the data choosing to forget the least correlated bounds between
#' \eqn{X} and \eqn{Y} datasets.
#'
#' @param Xs A data-frame of a matrix or a list of data-frames or matrices of \emph{n} rows each, the number of individuals. Some rows must be missing. The different matrices can have different numbers of columns. The length of Xs is denoted by \emph{K}.
#' @param Y A matrix of \emph{n} rows of a vector of length \emph{n} detailing the response matrix. No missing values are allowed in that matrix.
#' @param lambda A real \eqn{[0,1]} where 1 means just perfect correlations will be used and 0 no regularization is used.
#' @param R A strictly positive integer detailing the number of components to build in the model.
#' @param mode A character chain. Possibilities are "\emph{reg}", which implies  regression problem or anything else which means clustering is considered.  Default is "\emph{reg}".
#' @param errMin_imput Positive real. Minimal error in the Tribe Stage of the Koh-Lanta algorithm. Default is \eqn{1e-9}.
#' @param maxIter_imput Positive integer. Maximal number of iterations in the Tribe Stage of the Koh-Lanta algorithm. If equals to \eqn{0}, mean imputation is  considered. Default is \eqn{5}.
#' @param verbose Logical. If TRUE, the function cats specificities about the model. Default is FALSE.
#'
#' @return A list containing a mddsPLS object, see \code{\link{MddsPLS_core}}.
#'
#' @export
#'
#' @seealso \code{\link{predict.mddsPLS}}, \code{\link{perf_mddsPLS}}
#'
#' @examples
#' # Single-block example :
#' ## Classification example :
#' data("penicilliumYES")
#' X <- penicilliumYES$X
#' X <- scale(X[,which(apply(X,2,stats::sd)>0)])
#' Y <- as.factor(unlist(lapply(c("Melanoconidiu","Polonicum","Venetum"),function(tt){rep(tt,12)})))
#' mddsPLS_model_class <- mddsPLS(Xs = X,Y = Y,lambda = 0.958,R = 2,mode = "clas",verbose = TRUE)
#'
#' ## Regression example :
#' data("liver.toxicity")
#' X <- scale(liver.toxicity$gene)
#' Y <- scale(liver.toxicity$clinic)
#' mddsPLS_model_reg <- mddsPLS(Xs = X,Y = Y,lambda=0.9,R = 1, mode = "reg",verbose = TRUE)
#'
#' # Multi-block example :
#' ## Classification example :
#' data("penicilliumYES")
#' X <- penicilliumYES$X
#' X <- scale(X[,which(apply(X,2,stats::sd)>0)])
#' Xs <- list(X[,1:1000],X[,-(1:1000)])
#' Xs[[1]][1:5,]=Xs[[2]][6:10,] <- NA
#' Y <- as.factor(unlist(lapply(c("Melanoconidiu","Polonicum","Venetum"),function(tt){rep(tt,12)})))
#' mddsPLS_model_class <- mddsPLS(Xs = Xs,Y = Y,lambda = 0.95,R = 2,mode = "clas",verbose = TRUE)
#'
#' ## Regression example :
#' data("liver.toxicity")
#' X <- scale(liver.toxicity$gene)
#' Xs <- list(X[,1:1910],X[,-(1:1910)])
#' Xs[[1]][1:5,]=Xs[[2]][6:10,] <- NA
#' Y <- scale(liver.toxicity$clinic)
#' mddsPLS_model_reg <- mddsPLS(Xs = Xs,Y = Y,lambda=0.9,R = 1, mode = "reg",verbose = TRUE)
mddsPLS <- function(Xs,Y,lambda=0,R=1,mode="reg",
                    errMin_imput=1e-9,maxIter_imput=50,
                    verbose=FALSE){
  is.multi <- is.list(Xs)&!(is.data.frame(Xs))
  if(!is.multi){
    Xs <- list(Xs)
  }
  K <- length(Xs)
  ps <- lapply(Xs,ncol)
  Y_0 <- Y
  if(!(is.matrix(Y)|is.data.frame(Y))){
    Y <- as.matrix(Y)
  }
  n <- nrow(Y)
  q <- ncol(Y)
  has_converged <- maxIter_imput
  id_na <- lapply(Xs,function(x){which(is.na(x[,1]),arr.ind = TRUE)})
  if(length(unlist(id_na))==0){
    ## If ther is no missing sample
    mod <- MddsPLS_core(Xs,Y,lambda=lambda,R=R,mode=mode,verbose=verbose)
  }else{
    ## If ther are some missing samples
    for(k in 1:K){## ## Imputation to mean
      if(length(id_na[[k]])>0){
        mu_k <- colMeans(Xs[[k]],na.rm = TRUE)
        for(k_ik in 1:length(id_na[[k]])){
          Xs[[k]][id_na[[k]][k_ik],] <- mu_k
        }
      }
    }
    if(K>1){
      # Xs_init <- Xs
      mod_0 <- MddsPLS_core(Xs,Y,lambda=lambda,R=R,mode=mode)
      if(sum(abs(as.vector(mod_0$s)))!=0){
        Mat_na <- matrix(0,n,K)
        for(k in 1:K){
          Mat_na[id_na[[k]],k] <- 1
        }
        err <- 2
        iter <- 0
        while(iter<maxIter_imput&err>errMin_imput){
          iter <- iter + 1
          for(k in 1:K){
            if(length(id_na[[k]])>0){
              no_k <- (1:K)[-k]
              i_k <- id_na[[k]]
              Xs_i <- mod_0$s[-i_k,,drop=FALSE]
              newX_i <- mod_0$s[i_k,,drop=FALSE]
              ## ## ## Look for selected variables
              Var_selected_k <- which(rowSums(abs(mod_0$u[[k]]))!=0)
              if(length(Var_selected_k)>0){
                ## ## ## ## Impute on the selected variables
                Y_i_k <- Xs[[k]][-i_k,Var_selected_k,drop=FALSE]
                model_here <- MddsPLS_core(Xs_i,Y_i_k,lambda=lambda)
                mod_i_k <- list(mod=model_here,R=R,mode="reg",maxIter_imput=maxIter_imput)
                class(mod_i_k) <- "mddsPLS"
                Xs[[k]][i_k,Var_selected_k] <- predict.mddsPLS(mod_i_k,newX_i)
              }
            }
          }
          mod <- MddsPLS_core(Xs,Y,lambda=lambda,R=R,mode=mode)
          if(sum(abs(mod$t_frak))*sum(abs(mod_0$t_frak))!=0){
            err <- 0
            for(r in 1:R){
              n_new <- sqrt(sum(mod$t_frak[,r]^2))
              n_0 <- sqrt(sum(mod_0$t_frak[,r]^2))
              err <- err + abs(1-as.numeric(abs(diag(crossprod(mod$t_frak[,r],mod_0$t_frak[,r]))))/(n_new*n_0))
            }
          }
          else{
            err <- 0
          }
          if(iter>=maxIter_imput){
            has_converged <- 0
          }
          if(err<errMin_imput){
            has_converged <- iter
          }
          mod_0 <- mod
        }
      }
    }
    mod <- MddsPLS_core(Xs,Y,lambda=lambda,R=R,mode=mode,verbose=verbose)
  }
  out <- list(mod=mod,Xs=Xs,Y_0=Y_0,lambda=lambda,mode=mode,
              maxIter_imput=maxIter_imput,has_converged=has_converged)
  class(out) <- "mddsPLS"
  out
}
