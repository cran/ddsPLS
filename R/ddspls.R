#' C++ wrapper for bootstrap function
#'
#' The wrapper used to start the bootstrap commands. Not to be used by the user.
#'
#' @param U matrix, weights X
#' @param V matrix, weights Y
#' @param X matrix
#' @param Y matrix
#' @param lambdas vector, the to be tested values for lambda
#' @param lambda_prev vector, the previous selected values for lambda
#' @param R integer, the desired number of components
#' @param n_B integer, the number of bootstrap samples required
#' @param doBoot boolean, whether or not perform bootstrap. Used to build the
#' final model (FALSE)
#' @param n integer, the number of observations
#' @param p integer, the number of covariates
#' @param q integer, the number of response variables
#' @param n_lambdas integer, the number of to be tested lambdas
#' @param lambda0. the vector of lambda0
#'
#' @return A list
bootstrapWrap <- function(U,V,X,Y,lambdas,lambda_prev,
                          R,n_B,doBoot=TRUE,n,p,q,n_lambdas,lambda0.){
  res <- bootstrap_Rcpp(U,V,X,Y,lambdas,lambda_prev,
                        R,n_B,doBoot,n,p,q,n_lambdas,lambda0.)
  res
}

#' Data-Driven Sparse Partial Least Squares
#'
#' The main function of the package. It does both start the ddsPLS algorithm,
#' using bootstrap analysis. Also it estimates automatically the number of
#' components and the regularization coefficients.
#' One regularization parameter per component only is needed to select both in
#' \code{x} and in \code{y}. Build the optimal model, of the class
#' \code{ddsPLS}.
#' Among the different parameters, the \code{lambda} is the vector of parameters that are
#' tested by the algorithm along each component for each bootstrap sample. The total number
#' of bootstrap samples is fixed by the parameter \code{n_B}, for this parameter, the more
#'  the merrier, even if costs more in computation time.
#' This gives access to 3 S3 methods (\code{\link{summary.ddsPLS}}, \code{\link{plot.ddsPLS}} and \code{\link{predict.ddsPLS}}).
#'
#' @param X matrix, the covariate matrix (n,p).
#' @param Y matrix, the response matrix (n,q).
#' @param criterion character, whether \code{diffR2Q2} to be minimized, default,
#'  or \code{Q2} to be maximized.
#' @param doBoot logical, whether performing bootstrap operations, default to
#' \code{TRUE}. If equal to
#' \code{FALSE}, a model with is built on the parameters \code{lambda} and the
#'  number of components is the length of this vector.
#'   In that context, the parameter \code{n_B} is ignored. If equal to
#'    \code{TRUE}, the ddsPLS algorithm, through bootstrap validation,
#'     is started using \code{lambda} as a grid and \code{n_B} as
#'   the total number of bootstrap samples to simulate per component.
#' @param LD Boolean, wether or not consider Low-Dimensional dataset.
#' @param lambdas vector, the to be tested values for \code{lambda}.
#' Each value for \code{lambda} can be interpreted in terms of correlation
#'  allowed in the model.
#' More precisely, a covariate `x[j]` is not selected if its empirical
#' correlation with all the response variables `y[1..q]` is below \code{lambda}.
#'  A response variable `y[k]` is not selected if its empirical correlation
#'  with all the covariates `x[1..p]` is below \code{lambda}.
#' Default to \code{seq(0,1,length.out = 30)}.
#' @param n_B integer, the number of to be simulated bootstrap samples.
#' Default to \code{50}.
#' @param n_lambdas integer, the number of lambda values. Taken into account
#' only if \code{lambdas} is \code{NULL}. Default to 100.
#' @param lambda_roof limit value to be considered in the optimization.
#' @param lowQ2  real, the minimum value of Q^2_B to accept the
#' current lambda value. Default to \code{0.0}.
#' @param NCORES integer, the number of cores used. Default to \code{1}.
#' @param verbose boolean, whether to print current results. Defaut to
#' \code{FALSE}.
#' @param errorMin real, not to be used.
#'
#' @return A list with different interesting output describing the built model
#' @export
#' @importFrom foreach %dopar%
#' @importFrom foreach foreach
#' @importFrom foreach %do%
#' @importFrom parallel makeCluster
#' @importFrom parallel stopCluster
#' @importFrom doParallel registerDoParallel
#'
#' @examples
#' # n <- 100 ; d <- 2 ; p <- 20 ; q <- 2
#' # phi <- matrix(rnorm(n*d),n,d)
#' # a <- rep(1,p/4) ; b <- rep(1,p/2)
#' # X <- phi%*%matrix(c(1*a,0*a,0*b,
#' #                     1*a,3*b,0*a),nrow = d,byrow = TRUE) + matrix(rnorm(n*p),n,p)
#' # Y <- phi%*%matrix(c(1,0,
#' #                     0,0),nrow = d,byrow = TRUE) + matrix(rnorm(n*q),n,q)
#' # model_ddsPLS <- ddsPLS(X,Y,verbose=TRUE)
#'
#' @seealso \code{\link{summary.ddsPLS}}, \code{\link{plot.ddsPLS}}, \code{\link{predict.ddsPLS}}
#'
#' @useDynLib ddsPLS
ddsPLS <- function(X,Y,criterion="diffR2Q2",
                   doBoot=TRUE,LD=FALSE,
                   lambdas=NULL,n_B=50,n_lambdas=100,lambda_roof=NULL,
                   lowQ2=0.0,NCORES=1,errorMin=1e-9,verbose=FALSE){
  get_fused <- function(sigma_k,lambda1,lambda2)
  {
    # m2 <- c(flsa(sigma_k,lambda1=0,lambda2 = lambda2))
    m2 <- sigma_k
    m1 <- abs(m2)-lambda1
    m1[which(m1<0)] <- 0
    m1 <- m1*sign(m2)
    m1
  }
  getLambdas <- function(xSC,ySC,n,p,q){
    getLambda0 <- function(xSC,ySC,n,p,q){
      Sig_est <- matrix(rep(cov(xSC,ySC),n),nrow = n,byrow = TRUE)
      theta <- colMeans((xSC*matrix(rep(ySC,p),n,byrow = FALSE)-Sig_est)^2)
      mean(sqrt(log(max(p,q))*theta/n))
    }
    mean(unlist(lapply(1:q,function(j){getLambda0(xSC,ySC[,j],n,p,q)})))
  }
  gamma <- 0
  if(criterion %in% c("diffR2Q2","Q2")){
    call <- match.call()
    n <- nrow(Y)
    p <- ncol(X)
    q <- ncol(Y)
    # Standardize X and Y train and test.
    sdY <- apply(Y,2,sd)
    muY <- apply(Y,2,mean)
    id_na_y_mean <- which(is.na(sdY))
    if(length(id_na_y_mean)>0){
      muY[id_na_y_mean] = sdY[id_na_y_mean] <- 0
    }
    Y_init = scale(Y);
    Y_init[which(is.na(Y_init))] <- 0
    RSS0 <- sum(scale(Y,scale = FALSE)^2)
    RSS0_y <- apply(Y,2,function(yy){sum(scale(yy,scale = FALSE)^2)})
    sdX <- apply(X,2,sd)
    muX <- apply(X,2,mean)
    id_na_x_mean <- which(is.na(sdX))
    if(length(id_na_x_mean)>0){
      muX[id_na_x_mean] = sdX[id_na_x_mean] <- 0
    }
    X_init = scale(X);
    X_init[which(is.na(X_init))] <- 0
    sd_y_x_inv <- matrix(0,p,q)
    for(j in 1:q){
      sd_y_x_inv[,j] <- sdY[j]
    }
    for(i in 1:p){
      if(sdX[i]!=0){
        sd_y_x_inv[i,] <- 1/sdX[i]
      }else{
        sd_y_x_inv[i,] <- 0
      }
    }
    # Create lambda values
    useL0 <- F
    lambda0 <- NULL
    if(is.null(lambdas)){
      if(is.null(lambda_roof)){
        lambdas <- seq(0,max(abs(cov(X_init,Y_init))),length.out = n_lambdas)
      }else{
        lambdas <- seq(0,lambda_roof,length.out = n_lambdas+1)[-n_lambdas-1]
      }
      if(!LD){
        lambda0 <- c(lambda0,getLambdas(X_init,Y_init,n,p,q))
      }else{
        lambda0 <- 0
      }
      useL0 <- T
    }else{
      lambda0 <- 0
    }
    n_lambdas <- length(lambdas)
    # First covariance work
    COVInit = crossprod(Y_init,X_init)/(n-1);
    maxCOVInit = max(abs(COVInit))
    lambda_prev <- rep(0,n)
    TEST <- rep(0,n_lambdas)
    test_lambdas <- list()
    nb_ValsOk <- 0
    if(doBoot){
      U_out <- matrix(0,p,n); V0 <- matrix(0,q,n)
      B_previous <- matrix(0,p,q)
      R <- 1
      test <- TRUE
      R2Best <- rep(0,n)
      R2hBest <- rep(0,n)
      Q2Best <- rep(0,n)
      Q2hBest <- rep(0,n)
      explainedVar <- rep(0,n)
      Results <- list()
      varExplained = varExplainedTot = varExplained_y = varExplainedTot_y <- NULL
      Results$R2 = Results$R2h = Results$Q2 = Results$Q2h = Results$PropQ2hPos  <- list()
      Results$R2mean = Results$R2hmean = Results$Q2mean = Results$Q2hmean =
        Results$R2mean_diff_Q2mean = Results$t = Results$P = Results$C =
        Results$U_star = Results$U = Results$V = Results$B <- list()
      h <- 0 ; bestID <- 0;
      Q2_previous <- -1e9 ; bestVal <- -1e9
      if (verbose) {
        cat("                      ______________\n");
        cat("                     |    ddsPLS    |\n");
        cat("=====================----------------=====================\n");
      }
      while (test){
        if (verbose) {
          cat(paste("Should we build component " ,h+1 , " ? Bootstrap pending...\n",sep=""))
        }
        if(h>0 & useL0){
          ## Look at values for lambda_0
          X_here <- X_init
          Y_here <- Y_init
          for(hr in 1:h){
            tr <- resOUT$t[,hr]
            pr <- resOUT$P[,hr]
            cr <- resOUT$C[,hr]
            X_here <- X_here - tcrossprod(tr,pr)
            Y_here <- Y_here - tcrossprod(tr,cr)
          }
          lambda0 <- c(lambda0,getLambdas(X_here,Y_here,n,p,q))
          useL0 <- T
        }
        if(!useL0) lambda0 <- rep(0,h+1)
        NCORES_w <- min(NCORES,n_B)
        n_B_i <- ceiling(n_B/NCORES)
        `%my_do%` <- ifelse(NCORES_w!=1,{
          out<-`%dopar%`;cl <- makeCluster(NCORES_w)
          registerDoParallel(cl);out},{out <- `%do%`;out})
        res <- foreach(i_B=1:NCORES_w,
                       .combine='c',.multicombine=TRUE) %my_do% {
                         bootstrapWrap(U_out,V0,X_init,Y_init,lambdas,lambda_prev,
                                       R=h+1,n_B_i,doBoot=TRUE,n,p,q,n_lambdas,
                                       lambda0.=lambda0)
                       }
        if(NCORES_w!=1) stopCluster(cl)
        ## Criterion qulity descriptors
        Results$R2[[h+1]] <- do.call(rbind,res[which(names(res)=="R2")])
        Results$R2h[[h+1]] <- do.call(rbind,res[which(names(res)=="R2h")])
        Results$Q2[[h+1]] <- do.call(rbind,res[which(names(res)=="Q2")])
        Results$Q2h[[h+1]] <- do.call(rbind,res[which(names(res)=="Q2h")])
        Results$PropQ2hPos[[h+1]] <- apply(Results$Q2h[[h+1]],2,function(v){sum(v>=0)/length(v)})
        Results$R2mean[[h+1]] <- colMeans(Results$R2[[h+1]])
        Results$R2hmean[[h+1]] <- colMeans(Results$R2h[[h+1]])
        Results$Q2mean[[h+1]] <- colMeans(Results$Q2[[h+1]])
        Results$Q2hmean[[h+1]] <- colMeans(Results$Q2h[[h+1]])
        Results$R2mean_diff_Q2mean[[h+1]] <- Results$R2mean[[h+1]]-Results$Q2mean[[h+1]]
        Results$R2sd[[h+1]] <- apply(Results$R2[[h+1]],2,sd)
        Results$R2hsd[[h+1]] <- apply(Results$R2h[[h+1]],2,sd)
        Results$Q2sd[[h+1]] <- apply(Results$Q2[[h+1]],2,sd)
        Results$Q2hsd[[h+1]] <- apply(Results$Q2h[[h+1]],2,sd)
        Results$R2_diff_Q2sd[[h+1]] <- apply(Results$R2[[h+1]]-Results$Q2[[h+1]],2,sd)
        TEST <- (Results$Q2hmean[[h+1]]>lowQ2)*
          (Results$Q2mean[[h+1]]>Q2_previous)*
          (lambdas>=lambda0[h+1])==1#*(Results$PropQ2hPos[[h+1]]==1)
        nb_ValsOk = sum(TEST)
        test_lambdas[[h+1]] <- TEST
        # resOUT <- NULL
        test_t2 <- F
        if (nb_ValsOk>0){
          if(criterion=="diffR2Q2" & !LD){
            bestVal = min(Results$R2mean_diff_Q2mean[[h+1]][TEST])
            bestID = which(Results$R2mean_diff_Q2mean[[h+1]]==bestVal)[1]
          }else if(criterion=="Q2" | LD){
            bestVal = max(Results$Q2hmean[[h+1]][TEST])
            bestID = which(Results$Q2hmean[[h+1]]==bestVal)[1]
          }

          test_total <- TRUE
          while(test_total){
            lambda_prev[h+1] = lambdas[bestID]
            resMozna <- modelddsPLSCpp_Rcpp(U_out,V0,X_init,Y_init,
                                            lambda_prev,R=h+1,n,p,q,lambda0)
            test_t2 <- sum((resMozna$t[,h+1])^2)>errorMin
            if(test_t2){
              test_total <- FALSE
            }else{
              bestIDNext <- bestID - 1
              test_neighboor <- bestIDNext<=1
              if(!test_neighboor){
                bestID <- bestIDNext
                test_total <- TRUE
              }else{
                test_total <- FALSE
              }
            }
          }
          if(test_t2){
            resMozna -> resOUT
            resMozna <- NULL
            h <- h + 1
            for(hh in 1:h){
              P_test <- resOUT$P[,hh]
              if(which.max(P_test)!=which.max(abs(P_test))){
                resOUT$P[,hh] <- -resOUT$P[,hh]
                resOUT$C[,hh] <- -resOUT$C[,hh]
                resOUT$t[,hh] <- -resOUT$t[,hh]
                resOUT$V[,hh] <- -resOUT$V[,hh]
                resOUT$U[,hh] <- -resOUT$U[,hh]
                resOUT$U_star[,hh] <- -resOUT$U_star[,hh]
              }
            }
            Q2Best[h] = Results$Q2mean[[h]][bestID]
            Q2hBest[h] = Results$Q2hmean[[h]][bestID]
            R2Best[h] = Results$R2mean[[h]][bestID]
            R2hBest[h] = Results$R2hmean[[h]][bestID]
            Q2_previous = Q2Best[h]
            U_out[,1:h] = resOUT$U[,1:h]
            V0[,1:h] = resOUT$V[,1:h]
            ## Look at variabilities in P
            TT <- do.call(rbind,res[which(names(res)=="t")])[,c(1:n)+(bestID-1)*n,drop=F]
            PP <- do.call(rbind,res[which(names(res)=="P")])[,c(1:p)+(bestID-1)*p,drop=F]
            CC <- do.call(rbind,res[which(names(res)=="C")])[,c(1:q)+(bestID-1)*q,drop=F]
            UUSSTTAARR <- do.call(rbind,res[which(names(res)=="U_star")])[,c(1:p)+(bestID-1)*p,drop=F]
            UU <- do.call(rbind,res[which(names(res)=="U")])[,c(1:p)+(bestID-1)*p,drop=F]
            VV <- do.call(rbind,res[which(names(res)=="V")])[,c(1:q)+(bestID-1)*q,drop=F]
            for(i_B in 1:n_B){
              sign_i_B <- sign(VV[i_B,which.max(abs(VV[i_B,]))])
              if(sign_i_B<0){
                VV[i_B,] <- VV[i_B,]*(-1)
                PP[i_B,] <- PP[i_B,]*(-1)
                CC[i_B,] <- CC[i_B,]*(-1)
                UUSSTTAARR[i_B,] <- UUSSTTAARR[i_B,]*(-1)
                UU[i_B,] <- UU[i_B,]*(-1)
                TT[i_B,] <- TT[i_B,]*(-1)
              }
            }
            TT_boot_mean <- colMeans(TT)
            TT_boot_var <- apply(TT,2,var)
            VV_boot_mean <- colMeans(VV)
            VV_boot_var <- apply(VV,2,var)
            CC_boot_mean <- colMeans(CC)
            CC_boot_var <- apply(CC,2,var)
            PP_boot_mean <- colMeans(PP)
            PP_boot_var <- apply(PP,2,var)
            UU_boot_mean <- colMeans(UU)
            UU_boot_var <- apply(UU,2,var)
            UUSSTTAARR_boot_mean <- colMeans(UUSSTTAARR)
            UUSSTTAARR_boot_var <- apply(UUSSTTAARR,2,var)
            rm(VV,CC,PP,UU,UUSSTTAARR)
            Results$t[[h]] <- list(mean=TT_boot_mean,var=TT_boot_var)
            Results$U[[h]] <- list(mean=UU_boot_mean,var=UU_boot_var)
            Results$V[[h]] <- list(mean=VV_boot_mean,var=VV_boot_var)
            Results$C[[h]] <- list(mean=CC_boot_mean,var=CC_boot_var)
            Results$P[[h]] <- list(mean=PP_boot_mean,var=PP_boot_var)
            Results$U_star[[h]] <- list(mean=UUSSTTAARR_boot_mean,var=UUSSTTAARR_boot_var)
            ## Next
            B_out <- resOUT$B
            for (i in 1:p){
              if(sdX[i]>errorMin){
                B_out[i,] <- B_out[i,]/sdX[i]
              }
            }
            for (j in 1:q){
              B_out[,j] <- B_out[,j]*sdY[j]
            }
            B_out -> resOUT$B
            out0 <- list(model=list(muX=muX,muY=muY,B=B_previous),R=1);class(out0)="ddsPLS"
            out1 <- list(model=list(muX=muX,muY=muY,B=B_out),R=1);class(out1)="ddsPLS"
            Y_est_0 <- predict(out0,X,doDiagnosis=FALSE)$y_est
            Y_est_1 <- predict(out1,X,doDiagnosis=FALSE)$y_est
            cor2_0 <- unlist(lapply(1:q,function(j){
              1-sum((Y_est_0[,j]-Y[,j])^2)/sum((muY[j]-Y[,j])^2)
            }))
            cor2_1 <- unlist(lapply(1:q,function(j){
              1-sum((Y_est_1[,j]-Y[,j])^2)/sum((muY[j]-Y[,j])^2)
            }))
            ## Compute regression on current component only
            t_r <- resOUT$t[,h]
            Pi_r <- (abs(resOUT$V[,h])>1e-20)*1
            Y_est_r <- tcrossprod(t_r)%*%Y/sum(t_r^2)
            for(j in 1:q){
              if(Pi_r[j]==0){
                Y_est_r[,j] <- muY[j]
              }else{
                Y_est_r[,j] <- Y_est_r[,j] + muY[j]
              }
            }
            cor2_r <- unlist(lapply(1:q,function(j){
              1-sum((Y_est_r[,j]-Y[,j])^2)/sum((muY[j]-Y[,j])^2)
            }))
            ##
            B_previous <- B_out
            varExplained <- c(varExplained,mean(cor2_r)*100)#cor2_1-cor2_0)*100)
            varExplainedTot <- c(varExplainedTot,mean(cor2_1)*100)
            varExplained_y <- rbind(varExplained_y,cor2_r*100)#(cor2_1-cor2_0)*100)
            varExplainedTot_y <- rbind(varExplainedTot_y,(cor2_1)*100)
            if (verbose) {
              ress <- data.frame(
                list(
                  "  "="   ",
                  "lambda"=round(lambdas[bestID],2),
                  "R2"=round(Results$R2mean[[h]][bestID],2),
                  "R2h"=round(Results$R2hmean[[h]][bestID],2),
                  "Q2"=round(Results$Q2mean[[h]][bestID],2),
                  "Q2h"=round(Results$Q2hmean[[h]][bestID],2),
                  "VarExpl"=paste(round(varExplained[h]),"%",sep=""),
                  "VarExpl Tot"=paste(round(varExplainedTot[h]),"%",sep="")
                )
              )
              rownames(ress) <- ""
              colnames(ress)[1] <- "   "
              print(ress)
              cat(paste("                                     ...component ",h," built!","\n",sep=""))
            }
          }
        }
        if (!test_t2 | nb_ValsOk<=0){
          if(verbose) cat(paste("                                 ...component ",h+1," not built!","\n",sep=""))
          test = F;
          if(h==0){
            if(verbose){
              if(sum(Results$Q2hmean[[h+1]]>lowQ2)==0){
                cat("             ...no Q2r large enough for tested lambda.\n")
              }
            }
          }
        }
      }
      if(verbose) {
        cat("=====================                =====================\n");
        cat("                     ================\n");
      }
      lambda_sol=R2Sol=R2hSol=Q2Sol=Q2hSol <- rep(0,h)
      for (r in 0:h) {
        lambda_sol[r] = lambda_prev[r];
        R2Sol[r] = R2Best[r];
        R2hSol[r] = R2hBest[r];
        Q2Sol[r] = Q2Best[r];
        Q2hSol[r] = Q2hBest[r];
      }
      out <- list()
      if(h>0){
        out$model <- resOUT
        out$model$muY <- muY
        out$model$muX <- muX
        out$model$sdY <- sdY
        out$model$sdX <- sdX
        Results$lambdas <- lambdas
        out$results <- Results
        out$varExplained_in_X <- list()
        norm_X_tot <- (n-1)*p
        norm_comp <- colSums(out$model$t^2)
        var_comp <- norm_comp/norm_X_tot
        out$varExplained_in_X$Comp <- var_comp*100
        out$varExplained_in_X$Cumu <- cumsum(var_comp)*100
        out$varExplained <- list()
        out$varExplained$Comp <- varExplained
        out$varExplained$Cumu <- varExplainedTot
      }else{
        out$model = NULL
        # out$results <- NULL
        # out$resultsNotBuilt <- Results
        # out$resultsNotBuilt$lambdas <- lambdas
        out$results <- Results
        out$results$lambdas <- lambdas
        out$model$muY <- muY
        out$model$muX <- muX
        out$model$sdY <- sdY
        out$model$sdX <- sdX
        out$results$Q2hmean = out$results$Q2mean =
          out$results$R2hmean = out$results$R2mean = 0
      }
      out$R = h
      out$lambda = lambda_sol
      out$lambda_optim <- test_lambdas
      out$Q2 = Q2Sol
      out$Q2h = Q2hSol
      out$R2 = R2Sol
      out$R2h = R2hSol
      out$lowQ2=lowQ2
      out$X <- X
      out$doBoot <- doBoot
      class(out) <- "ddsPLS"
      if(h>0){
        out$Y_est <- predict(out,X,doDiagnosis=FALSE)$y_est
        out$Y_obs <- Y
        colnames(out$Y_est) = colnames(Y)
        rownames(out$model$V) = colnames(Y)
        rownames(out$model$U) = colnames(X)
        rownames(out$model$C) = colnames(Y)
        rownames(out$model$P) = colnames(X)
        rownames(out$model$U_star) = colnames(X)
        rownames(out$model$B) = colnames(X)
        colnames(out$model$B) = colnames(Y)
        out$varExplained_in_X <- list()
        norm_X_tot <- (n-1)*p
        norm_comp <- colSums(out$model$t^2)
        var_comp <- norm_comp/norm_X_tot
        out$varExplained_in_X$Comp <- var_comp*100
        out$varExplained_in_X$Cumu <- cumsum(var_comp)*100
        out$varExplained$PerY <- varExplainedTot_y[h,,drop=F]
        out$varExplained$PerYPerComp <- list()
        out$varExplained$PerYPerComp$Comp <- varExplained_y
        out$varExplained$PerYPerComp$Cumu <- varExplainedTot_y
        selX <- (rowSums(abs(out$model$B))>1e-9)*1
        selY <- (colSums(abs(out$model$B))>1e-9)*1
        out$Selection <- list(X=which(selX==1),Y=which(selY==1))
        out$call <- call
      }else{
        out$Y_est <- predict(out,X,doDiagnosis=FALSE)$y_est
        out$Y_obs <- Y
        colnames(out$Y_est) = colnames(Y)
      }
      out$criterion=criterion
      if (verbose & h>0) {
        plot(out)
      }
    }else{
      R <- length(lambdas)
      B_total <- matrix(0,p,q)
      if(R!=0){
        U_out <- matrix(0,p,R); V0 <- matrix(0,q,R)
        varExplained=varExplainedTot <- rep(0,R)
        varExplained_y=varExplainedTot_y <- matrix(0,R,q)
        lambda0 <- rep(0,R)
        if(!is.null(gamma))
        {
          x <- X_init
          y <- Y_init
          resr <- list(U=matrix(0,p,R),
                       V=matrix(0,q,R),
                       C=matrix(0,q,R),
                       P=matrix(0,p,R),
                       U_star=matrix(0,p,R),
                       t=matrix(0,n,R),
                       B=matrix(0,p,q))
        }
        for(r in 1:R){
          if (is.null(gamma))
          {
            resr <- modelddsPLSCpp_Rcpp(U_out,V0,X_init,Y_init,lambdas,
                                        R=r,n,p,q,lambda0)
            U_out[,r] = resr$U[,r]
            V0[,r] = resr$V[,r]
            # Compute regressions
            ## Compute regression on current component only
            t_r <- resr$t[,r]
            Pi_r <- (abs(resr$V[,r])>1e-20)*1
            Y_est_r <- tcrossprod(t_r)%*%Y/sum(t_r^2)
            for(j in 1:q){
              if(Pi_r[j]==0){
                Y_est_r[,j] <- muY[j]
              }else{
                Y_est_r[,j] <- Y_est_r[,j] + muY[j]
              }
            }
            cor2_r <- unlist(lapply(1:q,function(j){
              1-sum((Y_est_r[,j]-Y[,j])^2)/sum((muY[j]-Y[,j])^2)
            }))
            B_1 <- resr$B
            for (i in 1:p){
              if(sdX[i]>errorMin){
                B_1[i,] <- B_1[i,]/sdX[i]
              }
            }
            for (j in 1:q){
              B_1[,j] <- B_1[,j]*sdY[j]
            }
            out1 <- list(model=list(muX=muX,muY=muY,B=B_1),R=1);class(out1)="ddsPLS"
            Y_est_1 <- predict(out1,X,doDiagnosis=FALSE)$y_est
            cor2_1 <- unlist(lapply(1:q,function(j){
              1-sum((Y_est_1[,j]-Y[,j])^2)/sum((muY[j]-Y[,j])^2)
            }))
            # Compute explained variances
            varExplained[r] <- mean(cor2_r)*100
            varExplainedTot[r] <- mean(cor2_1)*100
            varExplained_y[r,] <- cor2_r*100
            varExplainedTot_y[r,] <- cor2_1*100
            B_total <- B_total + B_1
          }
          else
          {
            coco <- NULL
            for(k in 1:q)
            {
              coco <- cbind(coco,get_fused(c(crossprod(x,y[,k])/(n-1)),
                                           lambdas[r],gamma))
            }
            svdd <- svd(coco,nu=r,nv=r)
            ur <- svdd$u[,1]
            vr <- svdd$v[,1]
            resr$U[,r] <- ur
            resr$V[,r] <- vr
            U_out[,r] <- ur
            V0[,r] <- vr
            # Compute regressions
            ## Compute regression on current component only
            t_r <- x%*%ur#resr$t[,r]
            resr$t[,r] <- t_r
            normtr <- sum(t_r^2)
            Pi_r <- (abs(vr)>1e-20)*1
            if(normtr>1e-9)
            {
              pr <- crossprod(x,t_r)/normtr
              cr <- crossprod(y,t_r)/normtr
              for(j in 1:q){
                if(Pi_r[j]==0){
                  cr[j,] <- 0
                }
              }
            }
            else
            {
              pr <- rep(0,p)
              cr <- rep(0,q)
            }
            resr$P[,r] <- pr
            resr$U_star[,1:r] <- resr$U[,1:r]%*%solve(crossprod(resr$P[,1:r],resr$U[,1:r]))
            resr$C[,r] <- cr
            Y_est_r <- tcrossprod(t_r)%*%Y/sum(t_r^2)
            for(j in 1:q){
              if(Pi_r[j]==0){
                Y_est_r[,j] <- muY[j]
              }else{
                Y_est_r[,j] <- Y_est_r[,j] + muY[j]
              }
            }
            cor2_r <- unlist(lapply(1:q,function(j){
              1-sum((Y_est_r[,j]-Y[,j])^2)/sum((muY[j]-Y[,j])^2)
            }))
            B_1 <- ur%*%solve(crossprod(pr,ur))%*%t(cr)
            for (i in 1:p){
              if(sdX[i]>errorMin){
                B_1[i,] <- B_1[i,]/sdX[i]
              }
            }
            for (j in 1:q){
              B_1[,j] <- B_1[,j]*sdY[j]
            }
            out1 <- list(model=list(muX=muX,muY=muY,B=B_1),R=1)
            class(out1)="ddsPLS"
            Y_est_1 <- predict(out1,X,doDiagnosis=FALSE)$y_est
            cor2_1 <- unlist(lapply(1:q,function(j){
              1-sum((Y_est_1[,j]-Y[,j])^2)/sum((muY[j]-Y[,j])^2)
            }))
            # Compute explained variances
            varExplained[r] <- mean(cor2_r)*100
            varExplainedTot[r] <- mean(cor2_1)*100
            varExplained_y[r,] <- cor2_r*100
            varExplainedTot_y[r,] <- cor2_1*100
            B_total <- B_total + B_1
            x <- x - tcrossprod(t_r,pr)
            y <- y - tcrossprod(t_r,cr)
          }
        }
        resr$B <- B_total
        idBad <- which(sqrt(colSums(resr$t^2))<1e-9)
        if(length(idBad)>0){
          resr$U[,idBad] <- 0
          resr$V[,idBad] <- 0
        }
      }else{
        resr <- list(U=NULL,V=NULL,t=NULL,B=B_total)
        varExplained <- NULL
        varExplainedTot <- NULL
        varExplained_y <- NULL
        varExplainedTot_y <- NULL
      }
      out <- list()
      out$model <- resr
      out$criterion=criterion
      out$model$muY <- muY
      out$model$muX <- muX
      out$model$sdY <- sdY
      out$model$sdX <- sdX
      out$R <- R
      out$lambda = lambdas
      class(out) <- "ddsPLS"
      out$varExplained_in_X <- list()
      norm_X_tot <- (n-1)*p
      norm_comp <- colSums(out$model$t^2)
      var_comp <- norm_comp/norm_X_tot
      out$varExplained_in_X$Comp <- var_comp*100
      out$varExplained_in_X$Cumu <- cumsum(var_comp)*100
      out$varExplained <- list()
      out$varExplained$PerY <- varExplainedTot_y
      out$varExplained$Comp <- varExplained
      out$varExplained$Cumu <- varExplainedTot
      out$varExplained$PerYPerComp <- list()
      out$varExplained$PerYPerComp$Comp <- varExplained_y
      out$varExplained$PerYPerComp$Cumu <- varExplainedTot_y
      out$X <- X
      out$Y_est <- predict(out,X,doDiagnosis=FALSE)$y_est
      out$Y_obs <- Y
      colnames(out$Y_est) = colnames(Y)
      rownames(out$model$V) = colnames(Y)
      rownames(out$model$U) = colnames(X)
      rownames(out$model$C) = colnames(Y)
      rownames(out$model$P) = colnames(X)
      rownames(out$model$U_star) = colnames(X)
      rownames(out$model$B) = colnames(X)
      colnames(out$model$B) = colnames(Y)
      selX <- (rowSums(abs(out$model$B))>1e-9)*1
      selY <- (colSums(abs(out$model$B))>1e-9)*1
      out$Selection <- list(X=which(selX==1),Y=which(selY==1))
      out$call <- call
    }
  }else{
    out <- NULL
    cat("Please select a correct type of criterion among `diffR2Q2` (default), `Q2`.")
  }
  out$doBoot <- doBoot
  return(out)
}
