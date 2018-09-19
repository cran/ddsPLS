#' The predict function of a mdd-sPLS model
#'
#' @param object A mdd-sPLS object, output from the mddsPLS function.
#' @param newdata A data-set where individuals are described by the same as for mod_0
#' @param ... Other plotting parameters to affect the plot.
#'
#' @return A matrix of estimated \emph{Y_test} values.
#'
#' @importFrom stats predict
#'
#' @export
#'
#' @examples
#' data("liver.toxicity")
#' X <- scale(liver.toxicity$gene)
#' Y <- scale(liver.toxicity$clinic)
#' mod_0 <- mddsPLS(X,Y)
#' Y_test <- predict(mod_0,X)
predict.mddsPLS  <- function(object,newdata,...){
  mod_0 <- object
  newX <- newdata
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
            vars_Y_here[,C_pos+(0:(length(vars_k_id)-1))] <- mod_0$Xs[[pos_no_ok[k_id]]][,vars_k_id,drop=FALSE]
            C_pos <- C_pos + length(vars_k_id)
          }
        }
      }
      else{
        vars_Y_here <- matrix(0,nrow(t_X_here),1)
      }
      ## Generate model
      model_impute_test <- mddsPLS(t_X_here,vars_Y_here,lambda = lambda,R = R,maxIter_imput = mod_0$maxIter_imput)
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
            xx[id_xx,which(sd_x_0[[k_j]]!=0)] <-
              xx[id_xx,which(sd_x_0[[k_j]]!=0)]/sd_x_0[[k_j]][which(sd_x_0[[k_j]]!=0)]
          }
          t_X_test[,pos_col] <- xx%*%u_X_here[[k_j]][,r_j]
        }
      }
      ## Estimate missing values
      res <- predict.mddsPLS(model_impute_test,t_X_test)
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
      if(K>1 & mod_0$maxIter_imput>0){
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
      for(i in 1:n_new){
        newX[[k]][i,]<-(newX[[k]][i,]-mu_x_s[[k]])
        ok_sd <- which(sd_x_s[[k]]!=0)
        newX[[k]][i,ok_sd] <- newX[[k]][i,ok_sd]/sd_x_s[[k]][ok_sd]
      }
    }
    if(mode=="reg"){
      newY <- matrix(0,n_new,q)
      for(k in 1:K){
        newY <- newY + newX[[k]]%*%mod$B[[k]]
      }
      for(i in 1:n_new){
        newY[i,]<-newY[i,]*sd_y+mu_y
      }
    }
    else{
      t_r_new <- list()
      for(k in 1:K){
        if(k==1){
          for(r in 1:R){
            t_r_new[[r]] <- matrix(NA,n_new,K)
          }
        }
        for(r in 1:R){
          t_r_new[[r]][,k] <- newX[[k]]%*%mod_0$mod$u[[k]][,r]
        }
      }
      df_new <- data.frame(do.call(cbind,t_r_new)%*%mod_0$mod$beta_comb)
      colnames(df_new) <- paste("X",2:(ncol(df_new)+1),sep="")
      if(is.null(mod_0$mod$B)){
        newY <- list(class=sample(levels(mod_0$Y_0),size = 1,
                                  prob = table(mod_0$Y_0)/sum(table(mod_0$Y_0))))
      }
      else if(!is.null(mod_0$mod$B$sds)){
        pos_sds_0 <- 1+which(mod_0$mod$B$sds)
        newY <- predict(mod_0$mod$B,df_new[,c(1,pos_sds_0)])
      }else{
        newY <- predict(mod_0$mod$B,df_new)
      }
    }
  }
  else{
    newY <- matrix(NA,n_new,q)
    for(i_new in 1:n_new){
      newY[i_new,] <- predict.mddsPLS(mod_0,lapply(newX,
                                           function(nx,ix){
                                             nx[ix,,drop=FALSE]
                                           },i_new))
    }
  }
  newY
}
