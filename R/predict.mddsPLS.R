#' The predict method associated to the \strong{mddsPLS} class.
#'
#' @param object A mdd-sPLS object, output from the mddsPLS function.
#' @param newdata A data-set where individuals are described by the same as for mod_0
#' @param type charcter. It can be \strong{y} to return Y estimated value of \strong{x} for the completed values of newdata. \emph{both} for both \emph{y} and \emph{x}.
#' @param ... Other plotting parameters to affect the plot.
#'
#' @return Requested predicted values. In the case of classification, object \emph{probY} gives the probability per individual and per class.
#'
#' @importFrom stats predict
#'
#' @export
#'
#' @examples
#' data("liverToxicity")
#' X <- scale(liverToxicity$gene)
#' Y <- scale(liverToxicity$clinic)
#' mod_0 <- mddsPLS(X,Y)
#' Y_test <- predict(mod_0,X)
predict.mddsPLS  <- function(object,newdata,type="y",...){
  mod_0 <- object
  newX <- newdata
  object$L0 <- NULL
  #### FUNCTION
  fill_X_test <- function(mod_0,X_test){
    lambda <- mod_0$lambda
    R <- mod_0$mod$R
    id_na_test <- unlist(lapply(X_test,function(x){anyNA(x)}))
    mod <- mod_0$mod
    if(any(id_na_test)){
      ## Create covariable matrix train
      pos_ok <- which(!id_na_test)
      t_X_here <- do.call(cbind,lapply(1:R,function(ii,ti){
        ti[[ii]][,pos_ok]
      },mod$ts))
      u_X_here <- mod$u[pos_ok]
      mu_x_here <- mod$mu_x_s[pos_ok]
      sd_x_0 <- mod$sd_x_s[pos_ok]
      ## Create to be predicted matrix train
      pos_no_ok <- (1:K)[-pos_ok]
      pos_vars_Y_here <- lapply(mod$u[pos_no_ok],function(u){which(rowSums(abs(u))!=0)})
      if(sum(unlist(pos_vars_Y_here))!=0){
        nvars_Y_here_TOTAL <- length(unlist(pos_vars_Y_here))
        vars_Y_here <- matrix(0,nrow(t_X_here),nvars_Y_here_TOTAL)
        C_pos <- 1
        for(k_id in 1:length(pos_no_ok)){
          vars_k_id <- pos_vars_Y_here[[k_id]]
          if(length(vars_k_id)>0){
            to_use <- mod_0$Xs[[pos_no_ok[k_id]]][,vars_k_id,drop=FALSE]
            if(!is.matrix(to_use)){
              to_use <- as.matrix(to_use)
            }
            vars_Y_here[,C_pos+(0:(length(vars_k_id)-1))] <- to_use
            C_pos <- C_pos + length(vars_k_id)
          }
        }
      }
      else{
        vars_Y_here <- matrix(0,nrow(t_X_here),R)
      }
      ## Generate model
      model_impute_test <- mddsPLS(t_X_here,vars_Y_here,lambda = lambda,
                                   R = R,
                                   NZV=mod_0$NZV)
      ## Create test dataset
      n_test <- nrow(X_test[[1]])
      t_X_test <- matrix(NA,n_test,ncol(t_X_here))
      K_h <- sum(1-id_na_test)
      for(r_j in 1:R){
        for(k_j in 1:K_h){
          kk <- pos_ok[k_j]
          pos_col <- (r_j-1)*K_h+k_j
          xx <- X_test[[kk]]
          for(id_xx in 1:n_test){
            xx[id_xx,] <- xx[id_xx,]-mu_x_here[[k_j]]
            pos_sd_no_nul <- which(sd_x_0[[k_j]]>1e-10)
            if(length(pos_sd_no_nul)!=0){
              xx[id_xx,pos_sd_no_nul] <-
                xx[id_xx,pos_sd_no_nul]/sd_x_0[[k_j]][pos_sd_no_nul]
            }else{
              xx[id_xx,] <- 0
            }
          }
          t_X_test[,pos_col] <- xx%*%u_X_here[[k_j]][,r_j]
        }
      }
      ## Estimate missing values
      res <- predict.mddsPLS(model_impute_test,t_X_test)$y
      ## Put results inside Xs
      C_pos <- 1
      for(k_id in 1:length(pos_no_ok)){
        vars_k_id <- pos_vars_Y_here[[k_id]]
        X_test[[pos_no_ok[k_id]]] <- matrix(mod$mu_x_s[[pos_no_ok[k_id]]],nrow = 1)
        if(length(vars_k_id)>0){
          X_test[[pos_no_ok[k_id]]][1,vars_k_id] <- res[C_pos+(0:(length(vars_k_id)-1))]
          C_pos <- C_pos + length(vars_k_id)
        }
      }
    }
    X_test
  }
  #### END FUNCTION
  is.multi <- is.list(newX)&!(is.data.frame(newX))
  if(!is.multi){
    newX <- list(newX)
  }
  for(k in 1:length(newX)){
    if(is.data.frame(newX[[k]])){
      newX[[k]] <- as.matrix(newX[[k]])
    }
  }
  n_new <- nrow(newX[[1]])
  mod <- mod_0$mod
  q <- mod$q
  if(n_new==1){
    K <- length(newX)
    id_na_test <- unlist(lapply(newX,function(x){anyNA(x)}))
    if(any(id_na_test)){
      if(K>1){ ###  & mod_0$maxIter_imput>0
        newX <- fill_X_test(mod_0,newX)
      }
      else{
        for(k in 1:K){
          if(id_na_test[k]){
            newX[[k]][1,] <- mod_0$mod$mu_x_s[[k]]
          }
        }
      }
    }
    mode <- mod_0$mode
    Y_0 <- mod_0$Y_0
    mu_x_s <- mod$mu_x_s
    sd_x_s <- mod$sd_x_s
    mu_y <- mod$mu_y
    sd_y <- mod$sd_y
    R <- mod$R
    K <- length(mu_x_s)
    for(k in 1:K){
      newX[[k]][1,]<-(newX[[k]][1,]-mu_x_s[[k]])
      ok_sd <- which(sd_x_s[[k]]!=0)
      if(length(ok_sd)>0){
        newX[[k]][1,ok_sd] <- newX[[k]][1,ok_sd]/sd_x_s[[k]][ok_sd]
      }else{
        newX[[k]][1,] <- 0
      }
    }
    if(mode=="reg"){
      probability <- NULL
      newY <- matrix(0,n_new,q)
      for(k in 1:K){
        newY <- newY + newX[[k]]%*%mod$B[[k]]
      }
      for(i in 1:n_new){
        newY[i,]<-newY[i,] + mu_y
      }
    }
    else{
      T_super_new <- matrix(0,nrow=n_new,ncol=ncol(mod$T_super))
      for(k in 1:K){
          T_super_new <- T_super_new + newX[[k]]%*%mod_0$mod$u_t_super[[k]]
      }
      df_new <- data.frame(T_super_new)# df_new <- data.frame(do.call(cbind,T_super_new))#%*%mod_0$mod$beta_comb)
      colnames(df_new) <- paste("X",2:(ncol(df_new)+1),sep="")
      if(mod_0$mode=="lda"){
        probability <- rep(0,nlevels(Y_0))
        names(probability) <- levels(Y_0)
        if(class(mod_0$mod$B)=="list"){
          colnames(df_new) <- colnames(mod_0$mod$B$B$means)
        }else{
          colnames(df_new) <- colnames(mod_0$mod$B$means)
        }
        if(is.null(mod_0$mod$B)){
          newY <- list(class=sample(1:nlevels(mod_0$Y_0),size = 1,
                                    prob = table(mod_0$Y_0)/sum(table(mod_0$Y_0))))$'class'
        }
        else if(!is.null(mod_0$mod$B$sds)){
          pos_sds_no_0 <- which(mod_0$mod$B$sds!=0)
          prediction <- predict(mod_0$mod$B$B,df_new[,pos_sds_no_0,drop=F])
          newY <- prediction$'class'
          probability <- prediction$posterior
        }
        else{
          prediction <- predict(mod_0$mod$B,df_new)
          newY <- prediction$'class'
          probability <- prediction$posterior
        }
      }
      else if(mod_0$mode=="logit"){
        probability <- rep(0,nlevels(Y_0))
        names(probability) <- levels(Y_0)
        if(is.null(mod_0$mod$B)){
          newY <- list(class=sample(1:nlevels(mod_0$Y_0),size = 1,
                                    prob = table(mod_0$Y_0)/sum(table(mod_0$Y_0))))$'class'
        }
        else{
          colnames(df_new) <- names(mod_0$mod$B$coefficients)[-1]
          newY <- predict(mod_0$mod$B,df_new)
          prob <- suppressWarnings(predict(mod_0$mod$B,df_new,type="response"))
          probability <- c(1-prob,prob)
          names(probability) <- levels(Y_0)
          if(newY<0){
            newY <- levels(mod_0$Y_0)[1]
          }else{
            newY <- levels(mod_0$Y_0)[2]
          }
        }
      }
    }
    for(k in 1:K){
      newX[[k]][1,] <- newX[[k]][1,]*sd_x_s[[k]]
      newX[[k]][1,]<- newX[[k]][1,]+mu_x_s[[k]]
    }
  }
  else{
    if(mod_0$mode=="reg"){
      newY <- matrix(NA,n_new,q)
      probability <- NULL
    }else{
      newY <- rep(NA,n_new)
      probability <- matrix(0,n_new,nlevels(mod_0$Y_0))
      colnames(probability) <- levels(mod_0$Y_0)
    }
    for(i_new in 1:n_new){
      # Solved by Soso
      RES <- predict.mddsPLS(mod_0,lapply(newX,
                                          function(nx,ix){
                                            nx[ix,,drop=FALSE]
                                          },i_new),type="both")
      if(mod_0$mode=="reg"){
        newY[i_new,] <- RES$y
      }else{
        newY[i_new] <- as.character(RES$y)
        pro_i <- RES$probY
        for(jj in 1:length(pro_i)){
          n_i_jj <- names(pro_i)[jj]
          if(is.null(n_i_jj)) n_i_jj <- colnames(pro_i)[jj]
          pos_i_jj <- which(colnames(probability)==n_i_jj)
          probability[i_new,pos_i_jj] <- pro_i[jj]
        }
      }
      if(type=="x"|type=="both"){
        for(k in 1:length(newX)){
          if(anyNA(newX[[k]][i_new,])){
            newX[[k]][i_new,] <- RES$x[[k]][1,]
          }
        }
      }
    }
  }
  if(type=="y"){
    out <- list(y=newY,probY=probability)
  }else if(type=="x"){
    out <- newX
  }else{
    out <- list(x=newX,y=newY,probY=probability)
  }
  out
}
