#' The summary method of the \strong{perf_mddsPLS} function.
#'
#' This function is easy to use and gives information about the dataset and the model.
#'
#' @param object The object of class mddsPLS
#' @param plot_res_cv logical. If \strong{TRUE}, plots the results of the cross-validation
#' @param ... Other parameters.
#'
#' @importFrom graphics plot
#' @importFrom stats sd
#'
#' @export
#'
#' @seealso  \code{\link{perf_mddsPLS}}, \code{\link{plot.perf_mddsPLS}}
#'
#' @examples
#' library(ddsPLS)
#' data("liverToxicity")
#' X <- scale(liverToxicity$gene)
#' Y <- scale(liverToxicity$clinic)
#' X1<-X[,1:10];X1[1,]<-NA
#' X2<-X[,11:20];X2[2:5,]<-NA
#' X3<-X[,21:30];X3[4:20,]<-NA
#' X4<-X[,31:40]
#' Xs <- list(x1=X1,x2=X2,aaaa=X3,X4)
#' # object <- perf_mddsPLS(Xs = Xs,Y = Y[,1], lambdas=c(0.1,0.2,0.3),R = 1,
#' # mode = "reg",kfolds=5)
#' # summary(object)
summary.perf_mddsPLS <- function (object,plot_res_cv=T,
                                  ...)
{
  if(object$mode!="reg"){
    names(object)[which(names(object)=="ERROR")] <- "RMSEP"
  }
  is_L0 <- names(object$RMSEP)[2]
  K <- length(object$Xs);    sent_K <- paste("Number of blocks:",K)
  n <- nrow(object$Xs[[1]]);    sent_n <- paste("Number of individuals:",n)
  kfolds <- length(unique(object$fold))
  if(kfolds==0){
    if(object$kfolds=="loo"){
      kfolds <- n
    }else{
      kfolds <- object$kfolds
    }
  }
  if(kfolds==n){
    sent_kfolds <- "Performed Leave-One-Out cross-validation"
  }else{
    sent_kfolds <- paste("Performed ",kfolds,"-fold cross-validation",sep="")
  }
  ps <- unlist(lapply(object$Xs,ncol))
  id_na <- lapply(object$Xs,function(x){which(is.na(x[,1]))})
  na_x <- unlist(lapply(id_na,function(oo){if(length(oo)==0){out <- 0}else{out <- length(oo)};out}))
  prop_na <- signif(na_x/n*100,3)
  names_X_block <- names(object$Xs)
  if(length(names_X_block)==0){
    names_X_block <- 1:K
  }else{
    names_X_block <- unlist(lapply(1:K,function(k){
      out <- names_X_block[k]
      if(out==""){
        out <- paste("Block",k)
      }
      out
    }))
  }
  mat_miss <- matrix(NA,3,1+K)
  colnames(mat_miss) <- c(names_X_block,"Total")
  rownames(mat_miss) <- c("Number of variables","Number of missing samples","Proportion of missing samples (%)")
  mat_miss[1,] <- c(ps,sum(ps))
  mat_miss[2,] <- c(na_x,sum(na_x))
  mat_miss[3,] <- c(prop_na,signif(sum(prop_na)/K,3))
  df_miss <- data.frame(matrix(1,n,K))
  for(k in 1:K){popo <- as.numeric(object$id_na[[k]]);if(length(popo)>0){df_miss[popo,k]<-0};
  df_miss[,k] <- factor(df_miss[,k],levels = c(0,1))}
  names(df_miss) <- names_X_block
  q <- ncol(object$RMSEP)-2; sent_q <- paste("Number of variables in Y part:",q)
  mode <- object$mode;if(mode=="reg"){
    mode <- "regression"
  }else{
    mode <- "classification"
  }
  sent_mode <- paste("Model built in mode",mode)
  FREQ <- object$FREQ
  Conv <- object$Conv
  time <- object$time
  RMSEP <- object$RMSEP
  unik_paras <- RMSEP[,1:2]
  MAT_FINAL_RES <- data.frame(matrix(NA,nrow(unik_paras),5))
  MAT_FINAL_RES[,1:2] <- as.matrix(unik_paras)
  names(MAT_FINAL_RES) <- c("R","lambda","Error","Nb of convergences VS nb of fold","Mean(sd) time of computation")
  MAT_FINAL_RES$Error <- RMSEP$ERRORS_OUT
  Rs <- unique(unik_paras$R)
  ls <- unique(unik_paras$Lambdas)
  for(i in 1:nrow(MAT_FINAL_RES)){
    r <- MAT_FINAL_RES$R[i]
    lambda <- MAT_FINAL_RES$lambda[i]
    id <- which(Conv[1]==r&Conv[2]==lambda )
    ## Convergence
    num_conv <- Conv$number_iterations[id]
    num_conv <- paste(sum(num_conv!=0),kfolds,sep="/")
    MAT_FINAL_RES$`Nb of convergences VS nb of fold`[i] <- num_conv
    ## Time of computation
    time_i <- time$time_build[id]
    MAT_FINAL_RES$`Mean(sd) time of computation`[i] <- paste(signif(mean(time_i),2),
                                                             "(",signif(sd(time_i),2),")",sep="")
    ##
  }
  if(is_L0=="L0s"){
    names(MAT_FINAL_RES)[2] <- "L0"
  }
  cat("======================================================");cat("\n")
  cat("      Cross-Validation ddsPLS object description      ");cat("\n")
  cat("======================================================");cat("\n")
  cat("\n")
  cat(sent_K);cat("\n")
  cat(sent_n);cat("\n")
  cat(sent_q);cat("\n")
  cat(sent_mode);cat("\n")
  cat("\n")
  cat("\n")
  cat("    Missing value information    ");cat("\n")
  cat("---------------------------------");cat("\n")
  cat("\n")
  print(data.frame(mat_miss));cat("\n")
  cat("\n")
  cat("     Cross-Validation results    ");cat("\n")
  cat("---------------------------------");cat("\n")
  cat("\n")
  colnames(MAT_FINAL_RES)[4] <- "Mean AND sd time of computation"
  print(data.frame(MAT_FINAL_RES))
  cat("\n")
  cat("\n")
  cat("------------------------------------------------------");cat("\n")
  cat("======================================================");cat("\n")
  if(plot_res_cv){
    plot(object)
  }
}
