#' Function to compute cross-validation performances.
#'
#' That function must be applied to the given dataset and
#' the cross-validation process is made on the given set
#' of parameters.
#'
#' @param Xs A data-frame of a matrix or a list of data-frames or matrices of
#' $n$ rows each, the number of individuals. Some rows must be missing. The
#' different matrices can have different numbers of columns.
#' @param Y A matrix of n rows of a vector of length n detailing the
#' response matrix. No missing values are allowed in that matrix.
#' @param lambda_min A real in \eqn{[0,1]}. The minimum value considered.
#'  Default is \eqn{0}.
#' @param lambda_max A real in \eqn{[0,1]}. The maximum value considered.
#' Default is \eqn{NULL}, interpreted to the largest correlation between
#' \emph{X} and \emph{Y}.
#' @param n_lambda A strictly positive integer. Default to \eqn{1}.
#' @param lambdas A vector of reals in \eqn{[0,1]}. The values tested by the
#' perf process. Default is \eqn{NULL}, when that parameter is not taken into account.
#' @param R A strictly positive integer detailing the number of components to
#' build in the model.
#' @param kfolds character or integer. If equals to "loo" then a \emph{leave-one-out}
#' cross-validation is started. No other character is understood. Any strictly
#' positive integer gives the number of folds to make in the \emph{cross-validation process}
#' @param mode A character chain. Possibilities are "\emph{reg}", which implies
#'  regression problem or anything else which means clustering is considered.
#'   Default is "\emph{reg}".
#' @param fold_fixed Vector of length \eqn{n}. Each element corresponds to the
#' fold of the corresponding fold. If NULL then that argument is not considerd.
#' Default to NULL.
#' @param errMin_imput Positive real. Minimal error in the Tribe Stage of the
#' Koh-Lanta algorithm. Default is \eqn{1e-9}.
#' @param maxIter_imput Positive integer. Maximal number of iterations in the
#' Tribe Stage of the Koh-Lanta algorithm. If equals to \eqn{0}, mean imputation is
#'  considered. Default is \eqn{5}.
#' @param NCORES Integer. The number of cores. Default is \eqn{1}.
#'
#' @return A result of the perf function
#'
#' @import foreach
#'
#' @export
#'
#' @examples
#' library(doParallel)
#' # Classification example :
#' data("penicilliumYES")
#' X <- penicilliumYES$X
#' X <- scale(X[,which(apply(X,2,sd)>0)])
#' Y <- as.factor(unlist(lapply(c("Melanoconidiu","Polonicum","Venetum"),
#' function(tt){rep(tt,12)})))
#' #res_cv_class <- perf_mddsPLS(X,Y,lambda_min=0.85,n_lambda=2,R = 2,
#' #mode = "clas",NCORES = 1,fold_fixed = rep(1:12,3))
#'
#' # Regression example :
#' data("liver.toxicity")
#' X <- scale(liver.toxicity$gene)
#' Y <- scale(liver.toxicity$clinic)
#' #res_cv_reg <- perf_mddsPLS(Xs = X,Y = Y,lambda_min=0.8,n_lambda=2,R = 1,
#' # mode = "reg")
perf_mddsPLS <- function(Xs,Y,lambda_min=0,lambda_max=NULL,n_lambda=1,lambdas=NULL,R=1,kfolds="loo",
                         mode="reg",fold_fixed=NULL,maxIter_imput=5,errMin_imput=1e-9,NCORES=1){
  ## Xs shaping
  is.multi <- is.list(Xs)&!(is.data.frame(Xs))
  if(!is.multi){
    Xs <- list(Xs)
  }
  K <- length(Xs)
  ps <- lapply(Xs,ncol)
  ## Y shaping
  Y_0 <- Y
  if(mode=="reg"){
    if(!(is.matrix(Y)|is.data.frame(Y))){
      Y <- as.matrix(Y)
    }
    n<-nrow(Y);q <- ncol(Y)
  }else{
    n <- length(Y);q <- 1
    q_out <- nlevels(Y)
  }
  ## CV design
  if(kfolds=="loo"){
    kfolds <- n
    fold <- 1:n
  }
  else if(kfolds=="fixed"){
    fold <- fold_fixed
  }
  else{
    fold <- replicate(n/kfolds+1,sample(1:kfolds))[1:n]
  }
  ## Get highest Lambda
  if(is.null(lambdas)){
    if(is.null(lambda_max)){
      MMss0 <- mddsPLS(Xs,Y,lambda = 0,R = 1,mode = mode,maxIter_imput = 0)$mod$Ms
      lambda_max <- max(unlist(lapply(MMss0,
                                      function(Mi){max(abs(Mi))})))
    }
    Lambdas <- seq(lambda_min,lambda_max,length.out = n_lambda)
  }else{Lambdas <- lambdas}
  ## Write paras
  paras <- expand.grid(R,Lambdas,1:max(fold))
  if(NCORES>nrow(paras)){
    decoupe <- 1:nrow(paras)
  }else{
    decoupe <- replicate(nrow(paras)/NCORES + 1, sample(1:NCORES))[1:nrow(paras)]
  }
  NCORES_w <- min(NCORES,nrow(paras))
  if(NCORES_w!=1){
    cl <- parallel::makeCluster(NCORES_w)
    doParallel::registerDoParallel(cl)
  }
  pos_decoupe <- NULL
  ERRORS <- foreach::foreach(pos_decoupe=1:min(NCORES,nrow(paras)),
                    .combine = rbind,.packages = c("ddsPLS","MASS")) %dopar% {
                      paras_here_pos <- which(decoupe==pos_decoupe)
                      paras_here <- paras[paras_here_pos,,drop=FALSE]
                      if(mode=="reg"){
                        errors <- matrix(NA,nrow(paras_here),q)
                        select_y <- matrix(0,nrow(paras_here),q)
                      }
                      else{
                        errors <- rep(NA,nrow(paras_here))
                        select_y <- matrix(0,nrow(paras_here),nlevels(Y))
                      }
                      has_converged <- rep(0,nrow(paras_here))
                      time_build <- rep(0,nrow(paras_here))
                      for(i in 1:nrow(paras_here)){
                        R <- paras_here[i,1]
                        lambda <- paras_here[i,2]
                        i_fold <- paras_here[i,3]
                        pos_train <- which(fold!=i_fold)
                        t1 <- Sys.time()
                        X_train <- Xs
                        X_test <- Xs
                        for(k in 1:K){
                          X_train[[k]] <- X_train[[k]][pos_train,,drop=FALSE]
                          X_test[[k]] <- X_test[[k]][-pos_train,,drop=FALSE]
                        }
                        if(mode=="reg"){
                          Y_train <- Y[pos_train,,drop=FALSE]
                          Y_test <- Y[-pos_train,,drop=FALSE]
                        }
                        else{
                          Y_train <- Y[pos_train]
                          Y_test <- Y[-pos_train]
                        }
                        mod_0 <- mddsPLS(X_train,Y_train,lambda = lambda,
                                         R = R,mode = mode,errMin_imput = errMin_imput,
                                         maxIter_imput = maxIter_imput)
                        time_build[i] <- as.numeric((Sys.time()-t1))
                        has_converged[i] <- mod_0$has_converged
                        Y_est <- predict.mddsPLS(mod_0,X_test)
                        if(mode=="reg"){
                          errors_here <- Y_test-Y_est
                          errors[i,] <- sqrt(colMeans(errors_here^2))
                          v_no_null <- which(rowSums(abs(mod_0$mod$v))>1e-10)
                          select_y[i,v_no_null] <- 1
                        }
                        else{
                          errors[i] <- paste(as.character(Y_est$class),
                                             as.character(Y_test),sep="/")
                          v_no_null <- which(rowSums(abs(mod_0$mod$v))>1e-10)
                          select_y[i,v_no_null] <- 1
                        }
                      }
                      cbind(paras_here,errors,select_y,has_converged,time_build)
                    }
  if(NCORES_w!=1){
    parallel::stopCluster(cl)
  }

  paras_out <- expand.grid(R,Lambdas)
  ERRORS_OUT <- matrix(NA,nrow(paras_out),q)
  if(mode=="reg"){
    FREQ_OUT <- matrix(NA,nrow(paras_out),q)
  }
  else{
    ERRORS_OUT <- matrix(NA,nrow(paras_out),nlevels(Y))
    FREQ_OUT <- matrix(NA,nrow(paras_out),nlevels(Y))
  }
  for(i in 1:nrow(paras_out)){
    R <- paras_out[i,1]
    lambda <- paras_out[i,2]
    pos_in_errors <- intersect(which(ERRORS[,1]==R),which(ERRORS[,2]==lambda))
    if(mode=="reg"){
      ERRORS_OUT[i,] <- sqrt(colMeans(ERRORS[pos_in_errors,1:(q)+3,drop=FALSE]^2))
      FREQ_OUT[i,] <- colSums(ERRORS[pos_in_errors,1:(q)+3+q,drop=FALSE])
    }else{
      err_char <- ERRORS[pos_in_errors,1:(q)+3]
      err_char_spl <- do.call(rbind,strsplit(as.character(levels(err_char)[err_char]),
                                             split = "/",fixed = TRUE))
      colnames(err_char_spl) <- c("pred","obse")
      ERRORS_OUT[i,] <- abs(diag(table(err_char_spl[,1],err_char_spl[,2]))-table(err_char_spl[,2]))
      FREQ_OUT[i,] <- colSums(ERRORS[pos_in_errors,1:nlevels(Y)+4,drop=FALSE])
    }
  }
  out <- list(RMSEP=cbind(paras_out,ERRORS_OUT),FREQ=cbind(paras_out,FREQ_OUT),
              Conv=ERRORS[,c(1:3,ncol(ERRORS)-1)],time=ERRORS[,c(1:3,ncol(ERRORS))],mode=mode,Xs=Xs,Y=Y)
  class(out) <- "perf_mddsPLS"
  out
}
