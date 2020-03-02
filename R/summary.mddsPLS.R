#' The summary method of the \strong{mddsPLS} function.
#'
#' This function is easy to use and gives information about the dataset and the model.
#'
#' @param object The object of class mddsPLS
#' @param main_plot_indiv character. Main of the Venn diagram. Initialized to NULL.
#' @param fontsize interger. The size of the text, initialized to 10.
#' @param alpha real between 0 and 1. The transparency parameter.
#' @param ... Other parameters.
#'
#' @importFrom graphics plot
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
#'
#' @seealso  \code{\link{mddsPLS}}
#'
#' @export
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
#' # object <- mddsPLS(Xs = Xs,Y = Y[,1],lambda=0.1,R = 1, mode = "reg",verbose = TRUE)
#' # summary(object)
summary.mddsPLS <- function (object,
                             main_plot_indiv=NULL,
                             fontsize=10,alpha=0.7,
                             ...)
{
  NZV <- object$NZV
  K <- length(object$Xs);    sent_K <- paste("Number of blocks:",K)
  R <- object$mod$R;    sent_R <- paste("Number of dimensions:",R)
  if(!is.null(object$L0)){
    sent_lambda <- paste("Regularization coefficient:",object$L0)
  }else{
    sent_lambda <- paste("Regularization coefficient:",object$lambda)
  }
  n <- nrow(object$Xs[[1]])
  sent_n <- paste("Number of individuals:",n)
  ps <- unlist(lapply(object$Xs,ncol))
  na_x <- unlist(lapply(object$id_na,function(oo){if(length(oo)==0){out <- 0}else{out <- length(oo)};out}))
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
  names(df_miss) <- paste(names_X_block," (",unlist(lapply(object$id_na,length)),")",sep="")
  q <- ncol(object$Y_0)
  if(is.null(q)){q <- 1}
  sent_q <- paste("Number of variables in Y part:",q)
  mode <- object$mode;if(mode=="reg"){mode <- "regression"}else{mode <- "classification"}
  sent_mode <- paste("Model built in mode",mode)
  number_its <- object$number_iterations
  sent_its <- paste("Koh-Lanta process realized in",number_its,"iterations")

  df_num_var_sel <- data.frame(matrix(NA,K,R))
  rownames(df_num_var_sel) <- names_X_block
  colnames(df_num_var_sel) <- paste("Super Comp.",1:R)
  for(r in 1:R){
    for(k in 1:K){
      df_num_var_sel[k,r] <- length(which(abs(object$mod$u_t_super[[k]][,r])>NZV))
    }
  }


  cat("===============================================================");cat("\n")
  cat("                     ddsPLS object description              ");cat("\n")
  cat("===============================================================");cat("\n")
  cat("\n")
  cat(sent_K);cat("\n")
  cat(sent_R);cat("\n")
  cat(sent_lambda);cat("\n")
  cat(sent_n);cat("\n")
  cat(sent_q);cat("\n")
  cat(sent_mode);cat("\n")
  cat(sent_its);cat("\n")
  cat("\n")
  cat("\n")
  if(mode=="regression"){
    cat("          Variance Explained (%)    ");cat("\n")
    cat("-------------------------------------------");cat("\n")
    cat("Total Y variance explained by all the Super Components");cat("\n")
    oo <- as.numeric(object$Variances$Linear$VAR_GEN)
    print(signif(oo,2)*100)
    cat("Total Y variance explained by each Super Component");cat("\n")
    oo <- as.numeric(object$Variances$Linear$VAR_SUPER_COMPS_ALL_Y)
    names(oo) <- names(object$Variances$Linear$VAR_SUPER_COMPS_ALL_Y)
    print(signif(oo,2)*100)
    cat("\n")
    cat("Marginal Y variable variance explained by each Super Component");cat("\n")
    print(signif(object$Variances$Linear$VAR_SUPER_COMPS,2)*100)
    cat("\n")
    cat("Total Y variance explained by each component of each block");cat("\n")
    print(signif(object$Variances$Linear$VAR_COMPS,2)*100)
    cat("\n")
    cat("\n")
    cat("                    RV coefficient    ");cat("\n")
    cat("-------------------------------------------------------");cat("\n")
    cat("Total Y with all the Super Components");cat("\n")
    oo <- as.numeric(object$Variances$RV$VAR_GEN)
    print(signif(oo,2))
    cat("\n")
    cat("Total Y with each Super Component");cat("\n")
    oo <- as.numeric(object$Variances$RV$VAR_SUPER_COMPS_ALL_Y)
    names(oo) <- names(object$Variances$RV$VAR_SUPER_COMPS_ALL_Y)
    print(signif(oo,2))
    cat("\n")
    cat("Each Y variable with each Super Component");cat("\n")
    print(signif(object$Variances$RV$VAR_SUPER_COMPS,2))
    cat("\n")
    cat("Total Y with each component of each block");cat("\n")
    print(signif(object$Variances$RV$VAR_COMPS,2))
    cat("\n")
    cat("\n")
  }
  cat("         Missing value information    ");cat("\n")
  cat("-------------------------------------------");cat("\n")
  cat("\n")
  print(data.frame(mat_miss));cat("\n")
  cat("\n")
  cat("              mddsPLS results         ");cat("\n")
  cat("-------------------------------------------");cat("\n")
  cat("\n")
  N_max <- sum(unlist(lapply(object$mod$Ms,function(m){length(which(colSums(abs(m))!=0))})))
  cat(paste("At most ",N_max," variable(s) can be selected in the X part",sep=""));cat("\n")
  cat("\n")
  cat("\n")
  a<-lapply(object$mod$u,function(u){apply(u,2,function(u){length(which(abs(u)>NZV))})})
  cat(" ---- For each block of X, are selected");cat("\n")
  print(df_num_var_sel)
  if(mode=="regression"){
    cat(" ---- For the Y block, are selected");cat("\n")
    cat(paste("        @ (",paste(apply(object$mod$V_super,2,function(u){length(which(abs(u)>NZV))}),
                                  collapse = ","),") variable(s)",sep=""));cat("\n")
  }
  cat("\n")
  cat("\n")
  cat("---------------------------------------------------------------");cat("\n")
  cat("===============================================================");cat("\n")
}
