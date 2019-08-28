#' The core function of the Multi-Data-Driven sparse PLS function.
#'
#' This function should not be used directly by the user.
#'
#' @param Xs A matrix, if there is only one block, or a list of matrices,, if there is more than one block, of \strong{n} rows each, the number of individuals. Some rows must be missing. The different matrices can have different numbers of columns. The length of Xs is denoted by \strong{K}.
#' @param Y A matrix of n rows of a vector of length n detailing the response matrix. No missing values are allowed in that matrix.
#' @param lambda A real \eqn{[0,1]} where 1 means just perfect correlations will be used and 0 no regularization is used.
#' @param R A strictly positive integer detailing the number of components to build in the model.
#' @param L0 An integer non nul parameter giving the largest number of X variables that can be selected.
#' @param mode A character chain. Possibilities are "\strong{(reg,lda,logit)}", which implies regression problem, linear discriminant analysis (through the paclkage \code{MASS}, function \code{lda}) and logistic regression (function \code{glm}). Default is \strong{reg}.
#' @param id_na A list of na indices for each block. Initialized to NULL.
#' @param verbose Logical. If TRUE, the function cats specificities about the model. Default is FALSE.
#' @param NZV Float. The floatting value above which the weights are set to 0.
#'
#' @return A list containing the following objects:
#' \describe{
#'   \item{u}{A list of length \strong{K}. Each element is a \strong{p_kXR} matrix : the
#'    weights per block per axis.}
#'   \item{u_t_super}{A list of length \strong{K}. Each element is a \strong{p_kXR} matrix : the
#'    weights per block per axis scaled on the super description of the data set. Denoted as
#'    \strong{scaled super-weights}.}
#'   \item{v}{A \strong{qXR} matrix : the weights for the \strong{Y} part.}
#'   \item{ts}{A list of length \strong{R}. Each element is a \strong{nXK} matrix : the
#'    scores per axis per block.}
#'   \item{(t,s)}{Two \strong{nXR} matrices, super-scores of the \strong{X} and \strong{Y} parts.}
#'   \item{(t_ort,s_ort)}{Two \strong{nXR} matrices, final scores of the \strong{X} and \strong{Y} part.
#'    They correspond to \strong{PLS} scores of \strong{(t,s)} scores and so \strong{t_ort^T s_ort} is diagonal,
#'    \strong{t_ort}, respectively \strong{s_ort}, carries the same information as \strong{t}, respectively \strong{s}.}
#'   \item{B}{A list of length \strong{K}. Each element is a \strong{p_kXq} matrix : the
#'    regression matrix per block.}
#'   \item{(mu_x_s,sd_x_s)}{Two lists of length \strong{K}. Each element is a \strong{p_k} vector : the
#'    mean and standard deviation variables per block.}
#'   \item{(mu_y,sd_y)}{Two vectors of length \strong{q} : the mean and the standard deviation variables for \strong{Y} part.}
#'   \item{R}{Given as an input.}
#'   \item{q}{A non negative integer : the number of variables of \strong{Y} matrix. }
#'   \item{Ms}{A list of length \strong{K}. Each element is a \strong{qXp_k} matrix : the
#'    soft-thresholded empirical variance-covariance matrix \eqn{Y^TX_k/(n-1)}.}
#'   \item{lambda}{Given as an input.}
#' }
#'
#' @useDynLib ddsPLS
#' @importFrom Rcpp sourceCpp
#'
#' @importFrom stats sd model.matrix
#' @importFrom MASS lda
#'
MddsPLS_core <- function(Xs,Y,lambda=0,R=1,mode="reg",
                         L0=NULL,verbose=FALSE,id_na=NULL,
                         NZV=1e-9){

  my_scale <- function(a){
    if(!is.matrix(a)){
      a <- as.matrix(a,ncol=1)
    }
    if(!is.numeric(a)){
      a_ <- as.matrix(model.matrix( ~ y_obs - 1,
                                    data=data.frame(y_obs=a,ncol=1)))
      colnames(a_) <- levels(as.factor(a))
      a <- a_
    }
    scaleRcpp(a)
  }

  is.multi <- is.list(Xs)&!(is.data.frame(Xs))
  if(!is.multi){
    Xs <- list(Xs)
  }
  K <- length(Xs)
  ps <- lapply(Xs,ncol)
  ## Standardize Xs
  mu_x_s <- lapply(Xs,colMeans)
  n <- nrow(Xs[[1]])
  sd_x_s <- lapply(Xs,sdRcpp)#function(X){apply(X,2,sd)*sqrt((n-1)/n)})
  Xs <- lapply(Xs,my_scale)
  pos_0 <- lapply(sd_x_s,function(sdi){which(sdi<NZV)})
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
    if(is.data.frame(Y)){
      Y <- as.matrix(Y)
    }
  }
  else{
    Y_df <- data.frame(Y)
    momo <- model.matrix( ~ Y - 1, data=Y_df)
    attr(momo,"contrasts")=attr(momo,"assign")<-NULL
    Y <- my_scale(momo)
  }
  mu_y <- colMeans(Y)
  sd_y <- sdRcpp(Y)#apply(Y,2,function(y){sd(y)*sqrt((n-1)/n)})
  for(q_j in 1:length(sd_y)){
    if(sd_y[q_j]!=0){
      Y[,q_j] <- my_scale(Y[,q_j,drop=F])
    }
  }
  q <- ncol(Y)
  n <- nrow(Y)
  ## Create soft-thresholded matrices
  lambda_in <- lambda
  if(length(lambda_in)==1){
    lambda_in <- rep(lambda_in,K)
  }
  ps <- unlist(lapply(Xs,ncol))
  if(!is.null(L0)){
    sum_ps <- sum(ps);cum_ps <- cumsum(c(0,ps))
    all_maxs <- rep(NA,sum_ps)
    for(k in 1:K){
      ii <- cum_ps[k]
      if(is.null(id_na)){
        c_k <- suppressWarnings(abs(cor(Y,Xs[[k]])))
        c_k[which(is.na(c_k))] <- 0
      }else{
        if(length(id_na[[k]])>0){
          suppressWarnings(c_k <- abs(cor(Y[-id_na[[k]],,drop=F],Xs[[k]][-id_na[[k]],,drop=F])))
          c_k[which(is.na(c_k))] <- 0
        }else{
          suppressWarnings(c_k <- abs(cor(Y,Xs[[k]])))
          c_k[which(is.na(c_k))] <- 0
        }
      }
      all_maxs[ii+1:(ps[k])] <- apply(c_k,2,max)
      # if(is.null(id_na)){
      #   all_maxs[ii+1:(ps[k])] <- apply(abs(crossprod(Y,Xs[[k]])/(n-1)),MARGIN = 2,max)
      # }else{
      #   if(length(id_na[[k]])>0){
      #     all_maxs[ii+1:(ps[k])] <- apply(abs(crossprod(Y[-id_na[[k]],],Xs[[k]][-id_na[[k]],])/(n-1)),MARGIN = 2,max)
      #   }else{
      #     all_maxs[ii+1:(ps[k])] <- apply(abs(crossprod(Y,Xs[[k]])/(n-1)),MARGIN = 2,max)
      #   }
      # }
    }
    lambda_L0 <- sort(all_maxs,decreasing = T)[min(sum_ps,1+L0)]
    lambda_in <- rep(lambda_L0,K)
  }
  Ms <- lapply(1:K,function(k,Xs,Y,l,n){
    if(length(id_na[[k]])>0){
      M0 <- suppressWarnings(cor(Y[-id_na[[k]],],Xs[[k]][-id_na[[k]],]))
    }else{
      M0 <- suppressWarnings(cor(Y,Xs[[k]]))
    }
    M0[which(is.na(M0))] <- 0
    M <- abs(M0) - l[k]
    M[which(M<0)] <- 0
    M <- sign(M0)*M
    M
  },Xs,Y,lambda_in,n)
  if(verbose){
    N_max <- sum(unlist(lapply(Ms,function(m){length(which(colSums(abs(m))>NZV))})))
    cat(paste("At most ",N_max," variable(s) can be selected in the X part",sep=""));cat("\n")
  }
  ## Solve optimization problem
  tete <- min(R,q)#min(R,min(unlist(lapply(Ms,dim))))
  if(R<1){
    stop("Choose R superior to 1",
         call. = FALSE)
  }
  if(tete!=R){
    # warning(paste("R modified to ",tete," due to dimension consraints",sep=""),
    #         call. = FALSE)
    R <- tete
  }
  if(floor(R)!=ceiling(R)){
    R <- floor(R)
    warning(paste("R not integer and estimated to ",R,sep=""),
            call. = FALSE)
  }
  #### Inside problems
  u_t_r = u_t_r_0 <- list()
  t_r <- list()
  z_r <- list()
  z_t=t_t <- list()
  # BETA_r <- list()
  for(k in 1:K){
    if(norm(Ms[[k]])==0){
      svd_k <- list(v=matrix(0,nrow = ncol(Ms[[k]]),ncol = R),
                    d=rep(0,R))
    }
    else{
      # R_k <- min(R,min(dim(Ms[[k]])))
      # svd_k <- svd(Ms[[k]],nu = 0,nv = R)
      R_init <- min(q,sum(ps))
      svd_k_init <- svd(Ms[[k]],nu = 0,nv = R_init)
      eigen_YXk <- apply(mmultC(Ms[[k]],svd_k_init$v),2,function(t)sum(t^2))
      eigen_YXk[which(svd_k_init$d<NZV)] <- 0
      ordo_YXk <- order(eigen_YXk,decreasing = T)[1:min(R,length(eigen_YXk))]
      svd_k <- list(v=svd_k_init$v[,ordo_YXk,drop=F],
                    d=svd_k_init$d[ordo_YXk])
      ## Complete coefficients if needed
      length_val_prop <- length(svd_k$d)
      if(length_val_prop<R){
        svd_k$d <- c(svd_k$d,rep(0,R-length_val_prop))
      }
      ## Complete basis if needed
      ncol_V <- ncol(svd_k$v)
      if(ncol_V<R){
        svd_k$v <- cbind(svd_k$v,matrix(0,nrow(svd_k$v),R-ncol_V))
      }
      ## Test coefficients and put corresponding vectors to 0 if needed
      for(r in 1:R){
        if(svd_k$d[r]==0){
          svd_k$v[,r] <- 0
        }
      }
    }
    u_t_r[[k]] = u_t_r_0[[k]] <- svd_k$v
    if(k==1){
      for(r in 1:R){
        t_r[[r]] <- matrix(0,n,K)
        z_r[[r]] <- matrix(0,q,K)
      }
    }
    for(r in 1:R){
      if(svd_k$d[r]!=0){
        t_r[[r]][,k] <- mmultC(Xs[[k]],u_t_r[[k]][,r,drop=F])
        z_r[[r]][,k] <- mmultC(Ms[[k]],u_t_r[[k]][,r,drop=F])
      }
    }
    z_t[[k]] <- mmultC(Ms[[k]],u_t_r[[k]])
    t_t[[k]] <- mmultC(Xs[[k]],u_t_r[[k]])#crossprod(Y,Xs[[k]]%*%u_t_r[[k]])
  }
  ## Big SVD solution ######################### -----------------
  U_t_super <- list()
  Z <- do.call(cbind,z_t)
  R_opt <- R#min(min(dim(Z)),R)
  svd_Z <- svd(Z,nu = R_opt,nv = R_opt)
  # ## Reorder well ## -----
  # crosY_T <- crossprod(Y,do.call(cbind,t_t))
  # power <- colSums((crosY_T%*%svd_Z$v)^2)
  # orderGood <- order(power,decreasing = T)
  # svd_Z$d <- (svd_Z$d[orderGood])[1:R]
  # svd_Z$u <- (svd_Z$u[,orderGood,drop=F])[,1:R,drop=F]
  # svd_Z$v <- (svd_Z$v[,orderGood,drop=F])[,1:R,drop=F]
  ## END --- Reorder well ## -----
  beta_all <- svd_Z$v
  beta_list <- list()
  for(k in 1:K){#r in 1:R){
    #beta_list[[r]] <- beta_all[K*(r-1)+1:K,,drop=F]
    beta_list[[k]] <- beta_all[R*(k-1)+1:R,,drop=F]
  }
  V_super <- svd_Z$u
  if(length(svd_Z$d)<R){
    svd_Z$d <- c(svd_Z$d,rep(0,R-length(svd_Z$d)))
    V_super <- cbind(V_super,matrix(0,nrow=nrow(V_super),ncol=R-length(svd_Z$d)))
  }
  R_here <- ncol(beta_list[[1]])
  T_super <- matrix(0,nrow=n,ncol=R_here)
  for(k in 1:K){
    U_t_super[[k]] <- mmultC(u_t_r[[k]],beta_list[[k]])
    T_super <- T_super + mmultC(Xs[[k]],U_t_super[[k]])
  }
  ###########################################################
  vars_current <- rep(0,R_here)
  for(r in 1:R_here){
    sc_r <- T_super[,r,drop=F]#scale(x$mod$t[,r,drop=F],scale=F)
    var_t_super_r <- sum(sc_r^2)
    if(var_t_super_r!=0){
      if(n>q){
        deno_left <- norm(crossprod(Y),'f')
      }else{
        deno_left <- norm(tcrossprod(Y),'f')
      }
      deno <- deno_left*sum(sc_r^2)
      numer <- sum(mmultC(Y,crossprod(Y,sc_r))*sc_r)
      vars_current[r] <- numer/deno#sum(diag(mmultC(tcrossprod(Y),tcrossprod(sc_r))))/deno
    }
  }
  l_cur <- length(vars_current)
  if(l_cur>1){
    ord <- order(vars_current,decreasing = T)
    T_super <- T_super[,ord[1:R],drop=F]
    V_super <- V_super[,ord[1:R],drop=F]
    for(k in 1:K){
      U_t_super[[k]] <- U_t_super[[k]][,ord[1:R],drop=F]
      # u_t_r[[k]] <- u_t_r[[k]][,ord,drop=F]
      beta_list[[k]] <- beta_list[[k]][,ord[1:R],drop=F]
    }
  }
  ###########################################################
  ## -------------------------- ######################### -----------------
  S_super <- mmultC(Y,V_super)
  T_S <- crossprod(T_super,S_super)
  T_T <- crossprod(T_super)
  # svd_ort <- svd(S_T,nu = R,nv = R)
  svd_ort_T_super <- svd(T_super,nu = 0,nv = R)
  # u_ort <- svd_ort_T_super$u
  v_ort <- svd_ort_T_super$v
  Delta_ort <- svd_ort_T_super$d^2
  if(sum(Delta_ort)!=0){
    t_ort <- mmultC(T_super,v_ort)
    s_ort <- mmultC(S_super,v_ort)
    D_0_inv <- matrix(0,nrow = length(Delta_ort),ncol = length(Delta_ort))
    del_0 <- which(Delta_ort<NZV)
    if(length(del_0)>0){
      diag(D_0_inv)[-del_0] <- 1/Delta_ort[-del_0]
    }else{
      diag(D_0_inv) <- 1/Delta_ort
    }
    B_0 <- mmultC(mmultC(v_ort,tcrossprod(D_0_inv,v_ort)),T_S)
    # A <- matrix(0,R,R)
    # for(r in 1:R){
    #  A[r,r] <- as.numeric(crossprod(t_ort[,r],s_ort[,r]))/as.numeric(crossprod(t_ort[,r]))
    # }
  }else{
    t_ort=s_ort <- matrix(0,nrow = nrow(T_super),ncol=R)
    B_0 <- matrix(0,nrow = R,ncol=R)
    # A <- matrix(0,R,R)
    V_super <- matrix(0,q,R)
  }
  u <- beta_all#beta# deprecated
  if(mode=="reg"){
    B <- list()
    count <- 1
    for(k in 1:K){
      B_k <- tcrossprod(mmultC(U_t_super[[k]],B_0),V_super)
      if(anyNA(B_k)){
        B_k <- matrix(0,nrow(B_k),ncol(B_k))
      }
      B[[k]] <- B_k
    }
  }else{
    if(is.factor(Y_0)) Y_0 <- as.character(Y_0)
    dataf <- data.frame(cbind(Y_0,T_super));colnames(dataf)[1]<-"Y"
    for( cc in 2:ncol(dataf)){
      if(!is.numeric(dataf[,cc])){
        dataf[,cc] <- as.numeric(levels(dataf[,cc])[dataf[,cc]])
      }
    }
    sds <- apply(dataf[,-1,drop=FALSE],2,function(y){sd(y)*sqrt((n-1)/n)})
    if(any(abs(sds)<NZV)){
      pos_sd0 <- as.numeric(which(sds<NZV))
      if(length(pos_sd0)==length(sds)){
        B <- NULL
      }else{
        dataf <- dataf[,-c(1+pos_sd0)]
        if(mode=="lda"){
          B <- lda(Y ~ ., data = dataf)
          B <- list(B=B,sds=sds)
        }else if(mode=="logit"){
          if(!is.factor(dataf$Y)){
            if(min(dataf$Y)!=0){
              dataf$Y <- dataf$Y - 1
            }
          }
          B <- glm(Y ~ ., data = dataf,family = "binomial")
        }
      }
    }
    else{
      if(mode=="lda"){
        B <- lda(Y ~ ., data = dataf)
      }else if(mode=="logit"){
        if(!is.factor(dataf$Y)){
          dataf$Y <- factor(dataf$Y)
        }
        B <- glm(Y ~ ., data = dataf,family = "binomial")
      }
    }
  }
  if(verbose){
    a<-lapply(u_t_r,function(u){apply(u,2,function(u){length(which(abs(u)>NZV))})})
    cat("    For each block of X, are selected in order of component:");cat("\n")
    for(k in 1:K){
      cat(paste("        @ (",paste(a[[k]],collapse = ","),") variable(s)",sep=""));cat("\n")
    }
    cat("    For the Y block, are selected in order of component:");cat("\n")
    cat(paste("        @ (",paste(apply(V_super,2,function(u){length(which(abs(u)>NZV))}),
                                  collapse = ","),") variable(s)",sep=""));cat("\n")
  }
  list(u=u_t_r,u_t_super=U_t_super,V_super=V_super,ts=t_r,beta_comb=u,
       T_super=T_super,S_super=S_super,
       t_ort=t_ort,s_ort=s_ort,B=B,
       mu_x_s=mu_x_s,sd_x_s=sd_x_s,mu_y=mu_y,sd_y=sd_y,R=R,q=q,Ms=Ms,lambda=lambda_in[1])
}


#' Multi-Data-Driven sparse PLS function.
#'
#' This function takes a set \eqn{X} of \eqn{K} matrices defining the same \eqn{n} individuals and a matrix \eqn{Y} defining also those individuals. According to the num-
#' ber of components \eqn{R}, the user fixes the number of components the model
#' must be built on. The coefficient lambda regularizes the quality of proximity to the data choosing to forget the least correlated bounds between
#' \eqn{X} and \eqn{Y} data sets.
#'
#'
#' @param Xs A matrix, if there is only one block, or a list of matrices,, if there is more than one block, of \strong{n} rows each, the number of individuals. Some rows must be missing. The different matrices can have different numbers of columns. The length of Xs is denoted by \strong{K}.
#' @param Y A matrix of \strong{n} rows of a vector of length \strong{n} detailing the response matrix. No missing values are allowed in that matrix.
#' @param lambda A real \eqn{[0,1]} where 1 means just perfect correlations will be used and 0 no regularization is used.
#' @param R A strictly positive integer detailing the number of components to build in the model.
#' @param L0 An integer non nul parameter giving the largest number of X variables that can be selected.
#' @param keep_imp_mod Logical. Whether or not to keep imputation \strong{mddsPLS} models. Initialized to \code{FALSE} due to the potential size of those models.
#' @param reg_imp_model Logical. Whether or not to regularize the imputation models. Initialized to \code{TRUE}.
#' @param mode A character chain. Possibilities are "\strong{(reg,lda,logit)}", which implies regression problem, linear discriminant analysis (through the paclkage \code{MASS}, function \code{lda}) and logistic regression (function \code{glm}). Default is \strong{reg}.
#' @param errMin_imput Positive real. Minimal error in the Tribe Stage of the Koh-Lanta algorithm. Default is \eqn{1e-9}.
#' @param maxIter_imput Positive integer. Maximal number of iterations in the Tribe Stage of the Koh-Lanta algorithm. If equals to \eqn{0}, mean imputation is  considered. Default is \eqn{5}.
#' @param verbose Logical. If TRUE, the function cats specificities about the model. Default is FALSE.
#' @param NZV Float. The floatting value above which the weights are set to 0.
#' @param getVariances Logical. Whether or not to compute variances.
#'
#' @return A list containing a mddsPLS object, see \code{\link{MddsPLS_core}}. The \code{list} \code{order_values} is filled with the selected genes in each block.
#' They are oredered according to the sum of the square values of the \strong{Super-Weights} along the \code{R} dimensions.
#' The \code{rownames} give the names of the selected variables, if no name is given to the columns of \strong{Xs}, simply the indices are given.
#' Plus the \strong{Weights} and \strong{Super-Weights} are given for each of the selected variables in every \strong{R} dimension.
#' If \code{getVariances} is \code{TRUE} then the \code{Variances} is filled with two types of variances corresponding to bounds between components, or super-components and \strong{Y} vraiables, taken together or splitted.
#' Both of the types of variances are computed as follows:
#' \enumerate{
#' \item \strong{Linear}. Multivariate-linear regression matrix minimizing the Ordinary Least Squares problem is computed. Is then returned the fraction of the variance of the therefore model divide by the variance observed.
#' This represents the variance of the to be predicted parts by the predictors under a linear model.
#' \item \strong{RV}. That coefficient has permits to extend the correlation notion to matrices with the same number of rows but not necessarilye with the same number of columns \insertCite{@see @robert1976unifying}{ddsPLS}.
#' }
#'
#' @export
#' @useDynLib ddsPLS
#' @importFrom Rcpp sourceCpp
#' @importFrom stats na.omit glm
#'
#' @seealso \code{\link{summary.mddsPLS}}, \code{\link{plot.mddsPLS}}, \code{\link{predict.mddsPLS}}, \code{\link{perf_mddsPLS}}, \code{\link{summary.perf_mddsPLS}}, \code{\link{plot.perf_mddsPLS}}
#'
#' @references{
#'   \insertAllCited{}
#' }
#'
#' @examples
#' # Single-block example :
#' ## Classification example :
#' data("penicilliumYES")
#' X <- penicilliumYES$X
#' X <- scale(X[,which(apply(X,2,sd)>0)])
#' Y <- as.factor(unlist(lapply(c("Melanoconidiu","Polonicum","Venetum"),function(tt){rep(tt,12)})))
#' # mddsPLS_model_class <- mddsPLS(Xs = X,Y = Y,R = 2,L0=3,mode = "lda",verbose = TRUE)
#' # summary(mddsPLS_model_class,plot_present_indiv = FALSE)
#'
#' ## Regression example :
#' data("liverToxicity")
#' X <- scale(liverToxicity$gene)
#' Y <- scale(liverToxicity$clinic)
#' #mddsPLS_model_reg <- mddsPLS(Xs = X,Y = Y,L0=10,R = 1, mode = "reg",verbose = TRUE)
#' #summary(mddsPLS_model_reg)
#'
#' # Multi-block example :
#' ## Classification example :
#' data("penicilliumYES")
#' X <- penicilliumYES$X
#' X <- scale(X[,which(apply(X,2,sd)>0)])
#' Xs <- list(X[,1:1000],X[,-(1:1000)])
#' Xs[[1]][1:5,]=Xs[[2]][6:10,] <- NA
#' Y <- as.factor(unlist(lapply(c("Melanoconidiu","Polonicum","Venetum"),function(tt){rep(tt,12)})))
#' #mddsPLS_model_class <- mddsPLS(Xs = Xs,Y = Y,L0=3,mode = "lda",R = 2,verbose = TRUE)
#' #summary(mddsPLS_model_class)
#'
#' ## Regression example :
#' data("liverToxicity")
#' X <- scale(liverToxicity$gene)
#' Xs <- list(X[,1:1910],X[,-(1:1910)])
#' Xs[[1]][1:5,]=Xs[[2]][6:10,] <- NA
#' Y <- scale(liverToxicity$clinic)
#' #mddsPLS_model_reg <- mddsPLS(Xs = Xs,Y = Y,lambda=0.9,R = 1, mode = "reg",verbose = TRUE)
#' #summary(mddsPLS_model_reg)
mddsPLS <- function(Xs,Y,lambda=0,R=1,mode="reg",L0=NULL,
                    keep_imp_mod=FALSE,reg_imp_model=TRUE,
                    errMin_imput=1e-9,maxIter_imput=50,
                    verbose=FALSE,NZV=1E-9,getVariances=TRUE){

  my_scale <- function(a){
    if(!is.matrix(a)){
      a <- as.matrix(a,ncol=1)
    }
    if(!is.numeric(a)){
      a_ <- as.matrix(model.matrix( ~ y_obs - 1,
                                    data=data.frame(y_obs=a,ncol=1)))
      colnames(a_) <- levels(as.factor(a))
      a <- a_
    }
    scaleRcpp(a)
  }

  get_variances <- function(x,std_Y=T){
    get_var_line <- function(x,y){
      numer <- norm(mmultC(x,mmultC(solve(crossprod(x)),crossprod(x,y))),"f")
      denom <- norm(y,"f")
      numer/denom
    }
    get_rv <- function(x,y){
      numer <- sum(diag(tcrossprod(mmultC(x,crossprod(x,y)),y)))
      denom <- sum(diag(crossprod(x)))*sum(diag(crossprod(y)))
      numer/denom
    }

    Xs <- x$Xs
    K <- length(Xs)
    y_obs <- x$Y_0
    y_pred <- predict(x,Xs)
    mode <- x$mode
    if(mode=="reg"){
      if(!is.matrix(y_obs)&!is.data.frame(y_obs)){
        y_obs <- matrix(y_obs,ncol=1)
      }
      if(!is.matrix(y_pred)&!is.data.frame(y_pred)){
        y_pred <- matrix(y_pred,ncol=1)
      }
    }else{
      y_obs <- as.matrix(model.matrix( ~ y_obs - 1, data=data.frame(y_obs,ncol=1)))
      colnames(y_obs) <- levels(as.factor(x$Y_0))
      q <- ncol(y_obs)
      y_obs <- scale(y_obs,scale = F)
    }
    q <- ncol(y_obs)
    if(std_Y){
      sds <- sdRcpp(y_obs)#apply(y_obs,2,sd)
      pos_0 <- which(sds==0)
      y_obs<- my_scale(y_obs)
      if(length(pos_0)!=0){
        for(jj in pos_0){
          y_obs[,jj] <- 0
        }
      }
    }
    R <- length(x$mod$ts)
    VAR_TOT=VAR_TOT_FROB <- norm(y_obs,"f")
    VAR_GEN=VAR_GEN_FROB <- 0
    VAR_COMPS=VAR_COMPS_FROB <- matrix(0,K,R)
    VAR_SUPER_COMPS=VAR_SUPER_COMPS_FROB <- matrix(0,q,R)
    VAR_SUPER_COMPS_ALL_Y=VAR_SUPER_COMPS_ALL_Y_FROB <- matrix(0,1,R)

    t_y_obs <- tcrossprod(y_obs)

    T_GEN <- scale(x$mod$T_super)
    VAR_T_GEN <- norm(x$mod$T_super,"f")
    if(VAR_T_GEN!=0){
      # b <- mmultC(solve(crossprod(T_GEN)),crossprod(T_GEN,y_obs))
      VAR_GEN <- get_var_line(T_GEN,scale(y_obs))#(norm(mmultC(T_GEN,b),"f")/VAR_TOT)^2
      VAR_GEN_FROB <- get_rv(T_GEN,scale(y_obs))
    }

    for(r in 1:R){
      for(k in 1:K){
        t_r <- x$mod$ts[[r]]
        t_k_r <- scale(t_r[,k,drop=F],scale=F)
        var_t_k_r_all <- sum(t_k_r^2)
        if(var_t_k_r_all!=0){
          VAR_COMPS[k,r] <-  get_var_line(t_k_r,scale(y_obs))#(norm(mmultC(t_k_r,b),"f")/VAR_TOT)^2
          # deno <- sum(diag(t_y_obs))*var_t_k_r_all
          # numer <- sum(mmultC(y_obs,crossprod(y_obs,t_k_r))*t_k_r)
          # VAR_COMPS_FROB[k,r] <- numer/deno
          VAR_COMPS_FROB[k,r] <- get_rv(t_k_r,scale(y_obs))
        }
      }
    }
    for(j in 1:q){
      Y_j <- y_obs[,j,drop=F]
      var_j <- norm(Y_j,"f")
      if(var_j!=0){
        for(r in 1:R){
          t_super_r <- x$mod$T_super[,r,drop=F]
          var_t_super_r <- sum(t_super_r^2)
          if(var_t_super_r!=0){
            VAR_SUPER_COMPS[j,r] <- get_var_line(scale(t_super_r),scale(Y_j))#(norm(mmultC(t_super_r,b),"f")/var_j)^2
            # t_Y_j <- tcrossprod(Y_j)
            # t_t_super_r <- tcrossprod(t_super_r)
            # prod_num <- mmultC(t_Y_j,t_t_super_r)
            ### OO
            coef_r <- sum(Y_j*t_super_r)
            prod_num <- coef_r^2#*tcrossprod(Y_j,t_super_r)
            ###
            VAR_SUPER_COMPS_FROB[j,r] <- get_rv(scale(t_super_r),scale(Y_j))
          }
        }
      }
    }
    for(r in 1:R){
      sc_r <- scale(x$mod$T_super[,r,drop=F],scale=F)
      var_t_super_r <- sum(sc_r^2)
      if(var_t_super_r!=0){
        VAR_SUPER_COMPS_ALL_Y[r] <- get_var_line(sc_r,scale(y_obs))#(norm(mmultC(sc_r,b),"f")/VAR_TOT)^2
        VAR_SUPER_COMPS_ALL_Y_FROB[r] <- get_rv(sc_r,scale(y_obs))
      }
    }
    if(is.null(names(Xs))){
      legend_names_in <- paste("Block",1:K,sep=" ")
    }else{
      legend_names_in <- names(Xs)
      for(k in 1:K){
        if(nchar(names(Xs)[k])==0){
          legend_names_in[k] <- paste("Block",k)
        }
      }
    }
    rownames(VAR_COMPS)=rownames(VAR_COMPS_FROB) <- legend_names_in;
    colnames(VAR_COMPS)=colnames(VAR_COMPS_FROB) <- paste("Comp.",1:R)
    colnames(VAR_SUPER_COMPS)= colnames(VAR_SUPER_COMPS_ALL_Y) =
      colnames(VAR_SUPER_COMPS_FROB)= colnames(VAR_SUPER_COMPS_ALL_Y_FROB)<- paste("Super Comp.",1:R)
    rownames(VAR_SUPER_COMPS)=rownames(VAR_SUPER_COMPS_FROB) <- colnames(y_obs)
    return(list(
      Linear=list(VAR_GEN=VAR_GEN,VAR_SUPER_COMPS_ALL_Y=VAR_SUPER_COMPS_ALL_Y,
                  VAR_SUPER_COMPS=VAR_SUPER_COMPS,VAR_COMPS=VAR_COMPS),
      RV=list(VAR_GEN=VAR_GEN_FROB,VAR_SUPER_COMPS_ALL_Y=VAR_SUPER_COMPS_ALL_Y_FROB,
              VAR_SUPER_COMPS=VAR_SUPER_COMPS_FROB,VAR_COMPS=VAR_COMPS_FROB)))
  }

  if(lambda<0|lambda>1){
    stop("Choose lambda regularization parameter between 0 and 1",
         call. = FALSE)
  }
  is.multi <- is.list(Xs)&!(is.data.frame(Xs))
  if(!is.multi){
    Xs <- list(Xs)
  }
  K <- length(Xs)
  ps <- lapply(Xs,ncol)
  for(ii in 1:K){
    if(is.data.frame(Xs[[ii]])){
      Xs[[ii]] <- as.matrix(Xs[[ii]])
    }
  }
  Y_0 <- Y
  if(!(is.matrix(Y)|is.data.frame(Y))){
    Y <- as.matrix(Y)
  }
  n <- nrow(Y)
  q <- ncol(Y)
  if(keep_imp_mod) model_imputations <- list()
  has_converged <- maxIter_imput
  id_na <- lapply(Xs,function(x){which(is.na(x[,1]),arr.ind = TRUE)})
  any_na_no_all <- lapply(Xs,function(x){
    oo <- which(is.na(x),arr.ind = TRUE)[,1]
    pi <- ncol(x)
    table_o <- table(oo)
    toto <- which(table_o!=pi)
    out <- NA
    if(length(toto)!=0){
      out <- names(table_o)[toto]
    }
    as.numeric(out)
  })
  if(length(na.omit(unlist(any_na_no_all)))!=0){
    which.block <- which(unlist(lapply(any_na_no_all,function(oo){length(na.omit(oo))!=0})))
    mess1 <- "Block(s) with values missing not for all the variables:\n"
    mess2 <- paste("(",paste(which.block,collapse=","),")\n",sep="",collapse=",")
    mess3 <- "Corresponding individuals for each block:\n"
    ouou <- paste(unlist(lapply(which.block,function(i){paste(
      "(",paste(any_na_no_all[[i]],collapse=",",sep=""),
      ")",sep="")})),collapse=",")
    mess4 <- paste(ouou,"\n",sep="",collapse="")
    stop(paste(mess1,mess2,mess3,mess4),
         call. = FALSE)
  }
  mu_x_s <- lapply(Xs,colMeans,na.rm=T)
  sd_x_s <- lapply(1:length(ps),function(ii){
    p <- ps[ii]
    pos_na_ii <- which(is.na(Xs[[ii]][,1]))
    if(length(pos_na_ii)>0){
      out <- sdRcpp(na.omit(Xs[[ii]]))#apply(na.omit(Xs[[ii]]),2,sd)*sqrt((n-1-length(pos_na_ii))/(n-length(pos_na_ii)))
    }else{
      out <- sdRcpp(Xs[[ii]])
      #apply(Xs[[ii]],2,sd)*sqrt((n-1)/(n))
    }
    out
  })
  if(mode=="reg"){
    mu_y <- colMeans(Y)
    sd_y <- sdRcpp(Y)#apply(Y,2,sd)*sqrt((n-1)/n)
  }else{
    Y_class_dummies <- my_scale(model.matrix( ~ y - 1, data=data.frame(y=Y)))
    mu_y <- colMeans(Y_class_dummies)
    sd_y <- sdRcpp(Y_class_dummies)#apply(Y_class_dummies,2,sd)*sqrt((n-1)/n)
  }
  if(length(unlist(id_na))==0){
    ## If ther is no missing sample
    mod <- MddsPLS_core(Xs,Y,lambda=lambda,R=R,mode=mode,L0=L0,verbose=verbose,NZV=NZV)
  }else{
    if(!is.null(L0)){
      ps_init <- unlist(lapply(Xs,ncol))
      sum_ps_init <- sum(ps_init);cum_ps_init <- cumsum(c(0,ps_init))
      if(mode=="reg"){
        coco_i <- suppressWarnings(abs(cor(Y,do.call(cbind,Xs),use = "pairwise")))
        coco_i[which(is.na(coco_i))] <- 0
        all_maxs_init <- apply(coco_i,2,max)
      }else{
        coco_i <- suppressWarnings(abs(cor(Y_class_dummies,do.call(cbind,Xs),use = "pairwise")))
        coco_i[which(is.na(coco_i))] <- 0
        all_maxs_init <- apply(coco_i,2,max)
      }
      lambda_init <- sort(all_maxs_init,decreasing = T)[min(sum_ps_init,1+L0)]
    }else{
      lambda_init <- lambda
    }
    ## If ther are some missing samples
    for(k in 1:K){## ## Types of imputation for initialization
      if(length(id_na[[k]])>0){
        y_train <- Xs[[k]][-id_na[[k]],,drop=F]
        if(mode!="reg"){
          x_train <- Y_class_dummies[-id_na[[k]],,drop=F]
          x_test <- Y_class_dummies[id_na[[k]],,drop=F]
        }else{
          x_train <- Y[-id_na[[k]],,drop=F]
          x_test <- Y[id_na[[k]],,drop=F]
        }
        model_init <- mddsPLS(x_train,y_train,R=R,lambda = lambda_init,getVariances=F)
        y_test <- predict(model_init,x_test)
        Xs[[k]][id_na[[k]],] <- y_test
      }
    }
    if(K>1){
      mod_0 <- MddsPLS_core(Xs,Y,lambda=lambda,R=R,mode=mode,L0=L0,NZV=NZV)
      if(sum(abs(as.vector(mod_0$S_super)))!=0){
        Mat_na <- matrix(0,n,K)
        for(k in 1:K){
          Mat_na[id_na[[k]],k] <- 1
        }
        err <- 2
        iter <- 0
        ## Covariate for imputation is always the same : the projected values of Y on the initial weights
        #S_super_obj <- Y#mod_0$S_super
        if(mode!="reg"){
          S_super_obj <- Y_class_dummies
        }else{
          S_super_obj <- Y
        }
        Var_selected <- rep(NA,K)

        # ####################
        # Y_proj <- Y_class_dummies[,-1,drop=F]
        # phi <- tcrossprod(mmultC(Y_proj,solve(crossprod(Y_proj))),Y_proj)
        # ####################

        while(iter<maxIter_imput&err>errMin_imput){
          iter <- iter + 1
          for(k in 1:K){
            if(length(id_na[[k]])>0){
              no_k <- (1:K)[-k]
              i_k <- id_na[[k]]
              if(iter>1){
                # Xs_i <- mod$S_super[-i_k,,drop=FALSE]#_0$S_super[-i_k,,drop=FALSE]
                # newX_i <- mod$S_super[i_k,,drop=FALSE]#_0$S_super[i_k,,drop=FALSE]
                Var_selected_k <- which(rowSums(abs(mod$u[[k]]))>NZV)#_0$u[[k]]))!=0)
              }else{
                # Xs_i <- mod_0$S_super[-i_k,,drop=FALSE]
                # newX_i <- mod_0$S_super[i_k,,drop=FALSE]
                Var_selected_k <- which(rowSums(abs(mod_0$u[[k]]))>NZV)
              }
              # for(no_k_i in no_k){
              #   Var_selected_no_k <- which(rowSums(abs(mod_0$u[[no_k_i]]))>NZV)
              #   mat_here <- Xs[[no_k_i]][,Var_selected_no_k,drop=F]
              #   mat_here <- mat_here-mmultC(phi,mat_here)
              #   S_super_obj <- cbind(S_super_obj,mat_here)
              # }

              Xs_i <- S_super_obj[-i_k,,drop=FALSE]
              newX_i <- S_super_obj[i_k,,drop=FALSE]
              Var_selected[k] <- length(Var_selected_k)
              if(length(Var_selected_k)>0){
                ## ## ## ## Impute on the selected variables
                Y_i_k <- Xs[[k]][-i_k,Var_selected_k,drop=FALSE]
                if(reg_imp_model){
                  ## In that case the value of lambda is computed according to the initial model.
                  model_here <- MddsPLS_core(Xs_i,Y_i_k,lambda=mod_0$lambda,
                                             R=R,L0=NULL,NZV=NZV)
                }else{
                  ## In that case the value of lambda is taken equal to 0. And so
                  ## all the variables are taken in the model
                  model_here <- MddsPLS_core(Xs_i,Y_i_k,
                                             R=R,L0=NULL,NZV=NZV)
                }
                mod_i_k <- list(mod=model_here,R=R,mode="reg",maxIter_imput=maxIter_imput)
                class(mod_i_k) <- "mddsPLS"
                if(keep_imp_mod){
                  mod_i_k$Xs <- list(Xs_i)
                  mod_i_k$Y_0 <- Y_i_k
                  model_imputations[[k]] <- mod_i_k
                }
                Xs[[k]][i_k,Var_selected_k] <- predict.mddsPLS(mod_i_k,newX_i)
              }else{
                if(keep_imp_mod){
                  model_imputations[[k]] <- list()
                }
              }
            }else{
              if(keep_imp_mod){
                model_imputations[[k]] <- list()
              }
            }
          }
          if(mode!="reg"){
            if(!is.factor(Y))Y <- factor(Y)
          }
          mod <- MddsPLS_core(Xs,Y,lambda=mod_0$lambda,R=R,
                              mode=mode,L0=L0,NZV=NZV)#NULL)#######################L0)#
          if(sum(abs(mod$t_ort))*sum(abs(mod_0$t_ort))!=0){
            err <- 0
            # for(r in 1:R){
            #   n_new <- sqrt(sum(mod$t_ort[,r]^2))
            #   n_0 <- sqrt(sum(mod_0$t_ort[,r]^2))
            #   if(n_new*n_0!=0){
            #     err_i <- abs(1-as.numeric(abs(diag(crossprod(mod$t_ort[,r],
            #                                                  mod_0$t_ort[,r]))))/(n_new*n_0))
            #     err <- err + err_i
            #   }
            # }
            tsuper <- mod$T_super
            covsuper <- crossprod(tsuper)
            tsuper_0 <- mod_0$T_super
            covsuper_0 <- crossprod(tsuper_0)
            mix <- crossprod(tsuper_0,tsuper)
            numer <- sum(diag(tcrossprod(mix)))
            denom <- sqrt(sum(diag(mmultC(covsuper_0,covsuper_0)))*
                            sum(diag(mmultC(covsuper,covsuper))))
            if(denom>0){
              err <- 1-numer/denom
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
        if(keep_imp_mod){
          for(k in 1:K){
            if(length(id_na[[k]])>0 & Var_selected[k]>0){
              model_imputations[[k]]$Variances <- get_variances(model_imputations[[k]])
            }
          }
        }
      }
    }
    mod <- MddsPLS_core(Xs,Y,lambda=lambda,R=R,mode=mode,verbose=verbose,L0=L0,NZV=NZV)
  }
  mod$mu_x_s <- mu_x_s
  mod$sd_x_s <- sd_x_s
  mod$sd_y <- sd_y
  mod$mu_y <- mu_y
  var_selected <- list()
  for(k in 1:K){
    values <- rowSums(mod$u_t_super[[k]])^2
    pos <- which(values>NZV)
    if(length(pos)>0){
      order_values <- order(values[pos],decreasing = T)
      pos_ordered <- pos[order_values]
      out_k <- matrix(NA,length(pos),2*mod$R)
      coco_Xs_k <- colnames(Xs[[k]])
      if(is.null(coco_Xs_k)){
        rownames(out_k) <- pos_ordered
      }else{
        rownames(out_k) <- coco_Xs_k[pos_ordered]
      }
      colnames(out_k) <- c(paste("Weights_comp_",1:mod$R,sep=""),
                           paste("Super_Weights_comp_",1:mod$R,sep=""))
      for(r in 1:mod$R){
        out_k[,r] <- mod$u[[k]][pos_ordered,r]
        out_k[,mod$R+r] <- mod$u_t_super[[k]][pos_ordered,r]
      }
      var_selected[[k]] <- out_k
    }else{
      var_selected[[k]] <- "No variable selected"
    }
  }
  names_Xs <- names(Xs)
  if(length(names_Xs)!=0){
    names(var_selected) <- names_Xs
    names(mod$u) <- names_Xs
    names(mod$u_t_super) <- names_Xs
    if(mode=="regression")names(mod$B) <- names_Xs
    names(mod$Ms) <- names_Xs
  }
  out <- list(var_selected=var_selected,mod=mod,Xs=Xs,Y_0=Y_0,lambda=lambda,mode=mode,id_na=id_na,
              maxIter_imput=maxIter_imput,has_converged=has_converged,L0=L0,NZV=NZV)
  class(out) <- "mddsPLS"
  if(keep_imp_mod){
    if(length(names_Xs)!=0) names(model_imputations) <- names_Xs
    out$model_imputations <- model_imputations
  }
  if(getVariances){
    out$Variances <- get_variances(out)
  }
  out
}

