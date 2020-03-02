#' Function to compute cross-validation performances.
#'
#' That function must be applied to the given dataset and
#' the cross-validation process is made on the given set
#' of parameters.
#'
#' @param Xs A matrix, if there is only one block, or a list of matrices,
#'  if there is more than one block, of \strong{n} rows each, the number of individuals.
#'   Some rows must be missing. The different matrices can have different numbers of columns.
#'    The length of Xs is denoted by \strong{K}.
#' @param Y A matrix of n rows of a vector of length n detailing the
#' response matrix. No missing values are allowed in that matrix.
#' @param lambda_min A real in \eqn{[0,1]}. The minimum value considered.
#'  Default is \eqn{0}.
#' @param lambda_max A real in \eqn{[0,1]}. The maximum value considered.
#' Default is \eqn{NULL}, interpreted to the largest correlation between
#' \strong{X} and \strong{Y}.
#' @param n_lambda A strictly positive integer. Default to \eqn{1}.
#' @param lambdas A vector of reals in \eqn{[0,1]}. The values tested by the
#' perf process. Default is \eqn{NULL}, when that parameter is not taken into account.
#' @param L0s A vector of non null positive integers. The values tested by the
#' perf process. Default is \eqn{NULL} and is then not taken into account.
#' @param mu A real positive. The Ridge parameter changing the bias of the regression model. If is NULL, consider the classical ddsPLS. Default to NULL.
#' @param R A strictly positive integer detailing the number of components to
#' build in the model.
#' @param deflat Logical. If TRUE, the solution uses deflations to construct the weights.
#' @param weight Logical. If TRUE, the scores are divided by the number of selected variables of their corresponding block.
#' @param kfolds character or integer. If equals to "loo" then a \strong{leave-one-out}
#' cross-validation is started. No other character is understood. Any strictly
#' positive integer gives the number of folds to make in the \strong{cross-validation process}
#' @param mode A character chain. Possibilities are "\strong{(reg,lda,logit)}", which implies regression problem, linear discriminant analysis (through the paclkage \code{MASS}, function \code{lda}) and logistic regression (function \code{glm}). Default is \strong{reg}.
#' @param fold_fixed Vector of length \eqn{n}. Each element corresponds to the
#' fold of the corresponding fold. If NULL then that argument is not considerd.
#' Default to NULL.
#' @param NCORES Integer. The number of cores. Default is \eqn{1}.
#' @param NZV Float. The floatting value above which the weights are set to 0.
#' @param plot_result Logical. Wether or not to plot the result. Initialized to \strong{TRUE}. The \strong{reg_error} argument of the \strong{plot.perf_mddsPLS} function is left to its default value.
#' @param legend_label Logical. Wether or not to add the legend names to the plot. Initialized to \strong{TRUE}.
#'
#' @return A result of the perf function
#'
#' @import foreach
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#'
#' @seealso  \code{\link{summary.perf_mddsPLS}}, \code{\link{plot.perf_mddsPLS}}, \code{\link{mddsPLS}}, \code{\link{predict.mddsPLS}},
#'
#' @export
#'
#' @examples
#' # Classification example :
#' data("penicilliumYES")
#' X <- penicilliumYES$X
#' X <- scale(X[,which(apply(X,2,sd)>0)])
#' Y <- as.factor(unlist(lapply(c("Melanoconidiu","Polonicum","Venetum"),
#' function(tt){rep(tt,12)})))
#' #res_cv_class <- perf_mddsPLS(X,Y,L0s=1:5,R = 2,
#' #mode = "lda",NCORES = 1,fold_fixed = rep(1:12,3))
#'
#' # Regression example :
#' data("liverToxicity")
#' X <- scale(liverToxicity$gene)
#' Y <- scale(liverToxicity$clinic)
#' #res_cv_reg <- perf_mddsPLS(Xs = X,Y = Y,L0s=c(1,5,10,25,50),R = 1,
#' # mode = "reg")
perf_mddsPLS <- function(Xs,Y,lambda_min=0,lambda_max=NULL,n_lambda=1,lambdas=NULL,R=1,L0s=NULL,mu=NULL,
                         deflat=FALSE,weight=FALSE,
                         kfolds="loo",mode="reg",fold_fixed=NULL,NCORES=1,
                         NZV=1e-9,plot_result=T,legend_label=T){
  ## Xs shaping ##
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
    if(!is.factor(factor(Y))) Y <- as.factor(Y)
    n <- length(Y);q <- 1
    q_out <- nlevels(Y)
  }
  ## CV design
  if(kfolds=="loo" & is.null(fold_fixed)){
    kfolds <- n
    fold <- 1:n
  }else if(!is.null(fold_fixed)){
    fold <- fold_fixed
  }else{
    fold <- replicate(n/kfolds+1,sample(1:kfolds))[1:n]
  }
  ## Get highest Lambda
  if(is.null(lambdas)&is.null(L0s)){
    if(is.null(lambda_max)){
      MMss0 <- mddsPLS(Xs,Y,lambda = 0,R = 1,mu=mu,
                       mode = mode)$mod$Ms
      lambda_max <- max(unlist(lapply(MMss0,
                                      function(Mi){max(abs(Mi))})))
    }
    Lambdas <- seq(lambda_min,lambda_max,length.out = n_lambda)
  }else{Lambdas <- lambdas}
  if(!is.null(L0s)){
    ## Write paras
    paras <- expand.grid(R,L0s,1:max(fold))
  }else{
    ## Write paras
    paras <- expand.grid(R,Lambdas,1:max(fold))
  }
  if(NCORES>nrow(paras)){
    decoupe <- 1:nrow(paras)
  }else{
    decoupe <- replicate(nrow(paras)/NCORES + 1, sample(1:NCORES))[1:nrow(paras)]
  }
  NCORES_w <- min(NCORES,nrow(paras))
  `%my_do%` <- ifelse(NCORES_w!=1,{
    out<-`%dopar%`
    cl <- makeCluster(NCORES_w)#cl <- parallel::makeCluster(NCORES_w)
    registerDoParallel(cl)#doParallel::registerDoParallel(cl)
    out},{
      out <- `%do%`
      out})
  pos_decoupe <- NULL
  ERRORS <- foreach(pos_decoupe=1:min(NCORES,nrow(paras)),
                    .combine = rbind,.packages = c("ddsPLS","MASS")) %my_do% {
                      paras_here_pos <- which(decoupe==pos_decoupe)
                      paras_here <- paras[paras_here_pos,,drop=FALSE]
                      if(mode=="reg"){
                        errors <- matrix(NA,nrow(paras_here),q)
                        select_y <- matrix(0,nrow(paras_here),q)
                      }else{
                        errors <- rep(NA,nrow(paras_here))
                        select_y <- matrix(0,nrow(paras_here),nlevels(Y))
                      }
                      number_iterations <- rep(0,nrow(paras_here))
                      time_build <- rep(0,nrow(paras_here))
                      for(i in 1:nrow(paras_here)){
                        R <- paras_here[i,1]
                        if(!is.null(L0s)){
                          L0 <- paras_here[i,2]
                        }else{
                          lambda <- paras_here[i,2]
                        }
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
                        }else{
                          Y_train <- Y[pos_train]
                          Y_test <- Y[-pos_train]
                        }
                        if(!is.null(L0s)){
                          mod_0 <- mddsPLS(X_train,Y_train,L0 = L0,mu=mu,deflat=deflat,
                                           R = R,weight = weight,
                                           mode = mode,NZV=NZV,
                                           getVariances = F)
                        }else{
                          mod_0 <- mddsPLS(X_train,Y_train,lambda = lambda,mu=mu,deflat=deflat,
                                           R = R,weight = weight,
                                           mode = mode,NZV=NZV,
                                           getVariances = F)

                        }
                        time_build[i] <- as.numeric((Sys.time()-t1))
                        number_iterations[i] <- mod_0$number_iterations
                        Y_est <- predict.mddsPLS(mod_0,X_test)$y
                        if(mode=="reg"){
                          errors_here <- Y_test-Y_est
                          errors[i,] <- sqrt(colMeans(errors_here^2))
                          v_no_null <- which(rowSums(abs(mod_0$mod$V_super))>NZV)
                          select_y[i,v_no_null] <- 1
                        }else{
                          # if(mode!="lda"){
                          # Y_est <- factor(levels(Y)[Y_est],levels=levels(Y))
                          # }else{
                          Y_est <- factor(Y_est,levels=levels(Y))
                          # }
                          errors[i] <- paste(Y_est,Y_test,sep="/",collapse = " ")
                          v_no_null <- which(rowSums(abs(mod_0$mod$V_super))>NZV)
                          select_y[i,v_no_null] <- 1
                        }
                      }
                      out <- cbind(paras_here,errors,select_y,number_iterations,time_build)
                    }
  colnames(ERRORS)[1:3] <- c("R","L0","fold")
  if(NCORES_w!=1){
    stopCluster(cl)
  }
  if(!is.null(L0s)){
    paras_out <- expand.grid(R,L0s)
    colnames(paras_out) <- c("R","L0s")
  }else{
    paras_out <- expand.grid(R,Lambdas)
    colnames(paras_out) <- c("R","Lambdas")
  }
  ERRORS_OUT=MPE  <- matrix(NA,nrow(paras_out),q)
  SDEP_OUT  <- matrix(0,nrow(paras_out),q)
  if(mode=="reg"){
    FREQ_OUT <- matrix(NA,nrow(paras_out),q)
  }
  else{
    ERRORS_OUT <- matrix(NA,nrow(paras_out),nlevels(Y))
    SDEP_OUT  <- matrix(0,nrow(paras_out),nlevels(Y))
    FREQ_OUT <- matrix(NA,nrow(paras_out),nlevels(Y))
  }
  for(i in 1:nrow(paras_out)){
    R <- paras_out[i,1]
    if(!is.null(L0s)){
      L0 <- paras_out[i,2]
      pos_in_errors <- intersect(which(ERRORS[,1]==R),which(ERRORS[,2]==L0))
    }else{
      lambda <- paras_out[i,2]
      pos_in_errors <- intersect(which(ERRORS[,1]==R),which(ERRORS[,2]==lambda))
    }
    if(mode=="reg"){
      ERRORS_OUT[i,] <- sqrt(colMeans(ERRORS[pos_in_errors,1:(q)+3,drop=FALSE]^2))
      MPE[i,] <- colMeans(abs(ERRORS[pos_in_errors,1:(q)+3,drop=FALSE]))
      SDEP_OUT[i,] <- apply(ERRORS[pos_in_errors,1:(q)+3,drop=FALSE]^2,2,
                            function(y){sd(y)*sqrt((n-1)/n)})
      FREQ_OUT[i,] <- colSums(ERRORS[pos_in_errors,1:(q)+3+q,drop=FALSE])
    }else{
      ## err_char is structured as "Y_est/Y_observed"
      err_char <- ERRORS[pos_in_errors,1:(q)+3]
      mat_errors <- matrix(NA,length(fold),2)
      lev_Y <- levels(Y)
      for(ii in 1:max(fold)){
        i_fold_ii <- which(fold==ii)
        hihi <- unlist(strsplit(as.character(err_char[ii]),split = " ",fixed = TRUE))
        for(jj in 1:length (i_fold_ii)){
          mat_errors[i_fold_ii[jj],] <- unlist(strsplit(hihi[jj],split = "/",fixed = TRUE))
        }
      }
      classes <- lev_Y
      q_err <- length(classes)
      OUT_VEC <- rep(NA,q_err)
      for(i_q in 1:q_err){
        cla <- lev_Y[i_q]
        pos <- which(mat_errors[,2]==cla)
        OUT_VEC_i_q <- length(which(mat_errors[pos,1]==cla))
        ERRORS_OUT[i,i_q] <- length(pos)-OUT_VEC_i_q
      }
      colnames(ERRORS_OUT) <- classes
      FREQ_OUT[i,] <- colSums(ERRORS[pos_in_errors,1:nlevels(Y)+4,drop=FALSE])
    }
  }
  if(mode=="reg"){
    out <- list(RMSEP=cbind(paras_out,ERRORS_OUT),SDEP=cbind(paras_out,SDEP_OUT),
                FREQ=cbind(paras_out,FREQ_OUT),
                MPE=cbind(paras_out,MPE),
                Conv=ERRORS[,c(1:3,ncol(ERRORS)-1)],time=ERRORS[,c(1:3,ncol(ERRORS))],
                mode=mode,Xs=Xs,Y=Y,kfolds=kfolds,fold=fold,BackUp=ERRORS)
  }else{
    TAB <- table(Y)
    Precision <- cbind(paras_out,1-rowSums(ERRORS_OUT)/sum(TAB))
    colnames(Precision)[3] <- "Mean Precision"
    out <- list(ERROR=cbind(paras_out,ERRORS_OUT),
                SDEP=cbind(paras_out,SDEP_OUT),
                Precision=Precision,
                FREQ=cbind(paras_out,FREQ_OUT),
                Conv=ERRORS[,c(1:3,ncol(ERRORS)-1)],time=ERRORS[,c(1:3,ncol(ERRORS))],
                mode=mode,Xs=Xs,Y=Y,kfolds=kfolds,fold=fold,BackUp=ERRORS)
  }
  class(out) <- "perf_mddsPLS"
  if(plot_result){
    if(legend_label){
      plot(out,no_occurence=T,plot_mean = T,legend_names=colnames(out$ERROR)[-c(1:2)])
    }else{
      plot(out,no_occurence=T,plot_mean = T)
    }
  }
  res_plot_no_plot <- plot(out,no_plot=T)
  out$Optim <- res_plot_no_plot[1:2]
  out
}
