#' Function to predict from ddsPLS objects
#'
#' @param object A ddsPLS object.
#' @param X_test matrix, a test data-set. If is "NULL", the default value,
#' the predicted values for the train test are returned.
#' @param toPlot boolean, wether or not to plot the extreme value test plot.
#' Default to `TRUE`.
#' @param legend.position character. Where to put the legend.
#' @param cex float positive. Number indicating the amount by which plotting
#' symbols should be scaled relative to the default.
#' @param cex.text float positive. Number indicating the amount by which
#' plotting text elements should be scaled relative to the default.
#' @param doDiagnosis Yes or no to perform diagnosis.
#' @param ... Other parameters
#'
#' @export
#'
#' @seealso \code{\link{ddsPLS}}, \code{\link{plot.ddsPLS}}, \code{\link{summary.ddsPLS}}
#'
#' @useDynLib ddsPLS
predict.ddsPLS <- function(object,X_test=NULL,toPlot=FALSE,doDiagnosis=T,
                           legend.position="topright",cex=1,cex.text=1,...){
  x <- object
  getDiagnoses <- function(x,X_test,y_test_est){
    n_test <- nrow(X_test)
    n_train <- nrow(x$X)
    X_test_center <- X_test-matrix(rep(x$model$muX,n_test),nrow = n_test,byrow = T)
    X_train_center <- x$X-matrix(rep(x$model$muX,n_train),nrow = n_train,byrow = T)
    ## On the predictions
    y_train_est <- x$Y_est
    sigma_y_inv <- unlist(lapply(x$model$sdY,function(ss){
      if(ss>1e-9){out <- 1/ss}else{out <- 0};out}))
    y_mean_train <- matrix(rep(colMeans(y_train_est),n_train),nrow = n_train,byrow = T)
    y_mean_test <- matrix(rep(colMeans(y_train_est),n_test),nrow = n_test,byrow = T)
    dist_y_Mahalanobis_train <- (rowMeans(
      (y_train_est-y_mean_train)^2*matrix(rep(sigma_y_inv,n_train),nrow = n_train,byrow = T)^2))
    dist_y_Mahalanobis_test <- (rowMeans(
      (y_test_est-y_mean_test)^2*matrix(rep(sigma_y_inv,n_test),nrow = n_test,byrow = T)^2))
    ## On the scores
    t_test <- X_test_center%*%x$model$U_star
    t_train <- X_train_center%*%x$model$U_star
    sigma_t_inv <- unlist(lapply(sqrt(colMeans(t_train^2)),function(ss){
      if(ss>1e-9){out <- 1/ss}else{out <- 0};out}))
    dist_score_Mahalanobis_train <- (rowMeans(
      t_train^2*matrix(rep(sigma_t_inv,n_train),nrow = n_train,byrow = T)^2))
    dist_score_Mahalanobis_test <- (rowMeans(
      t_test^2*matrix(rep(sigma_t_inv,n_test),nrow = n_test,byrow = T)^2))
    ## On the reconstruction of X
    X_test_center_est <- tcrossprod(t_test,x$model$P)
    X_train_center_est <- tcrossprod(t_train,x$model$P)
    epsilon_X_test <- X_test_center_est-X_test_center
    epsilon_X_train <- X_train_center_est-X_train_center
    # sigmaInv <- unlist(lapply(x$model$sdX,function(ss){
    #   if(ss>1e-9){out <- 1/ss}else{out <- 0};out}))
    sigma_x_Inv <- unlist(lapply(sqrt(colMeans(epsilon_X_train^2)),function(ss){
      if(ss>1e-9){out <- 1/ss}else{out <- 0};out}))
    dist_Mahalanobis_train <- (rowMeans(
      epsilon_X_train^2*matrix(rep(sigma_x_Inv,n_train),nrow = n_train,byrow = T)^2))
    dist_Mahalanobis_test <- (rowMeans(
      epsilon_X_test^2*matrix(rep(sigma_x_Inv,n_test),nrow = n_test,byrow = T)^2))
    list(# y=list(train=dist_y_Mahalanobis_train,test=dist_y_Mahalanobis_test),
      t=list(train=dist_score_Mahalanobis_train,test=dist_score_Mahalanobis_test),
      x=list(train=dist_Mahalanobis_train,test=dist_Mahalanobis_test))
  }
  diagnoses <- NULL
  if(is.null(X_test)){
    y_est <- x$Y_est
  }else{
    n_test <- nrow(X_test)
    if(x$R==0){
      y_est <- matrix(rep(x$model$muY,n_test),nrow = n_test,byrow = T)
    }else{
      y_est <- (X_test-matrix(rep(x$model$muX,n_test),nrow = n_test,byrow = T))%*%x$model$B
      y_est <- y_est + matrix(rep(x$model$muY,n_test),nrow = n_test,byrow = T)
    }
    if(doDiagnosis & x$R>0){
      diagnoses <- getDiagnoses(x,X_test,y_est)
      if(toPlot){
        xlim <- range(c(diagnoses$t$train,diagnoses$t$test))
        ylim <- range(c(diagnoses$x$train,diagnoses$x$test))
        N <- 4
        par(mar=c(5, 4, 1, 0.5) + 0.1)
        layout(matrix(c(rep(2,N),4,rep(c(rep(1,N),3),N)),nrow = N+1,byrow = T))
        plot(diagnoses$t$train,diagnoses$x$train,col=1,pch=1,
             cex=cex,xlab="",ylab="",
             xlim=xlim,ylim=ylim)
        title(xlab = expression(d[bold(t)](bold(t))),line = 2.5,cex.main=cex.text)
        title(ylab = expression(d[bold(x)](hat(bold(x)),bold(x))),
              line = 2.5,cex.main=cex.text)
        points(diagnoses$t$test,diagnoses$x$test,col="red",pch=16,cex=cex)
        legend(legend.position,c("Train","Test"),
               col=c(1,2),pch=c(1,16),cex=cex.text,pt.cex = cex)
        ##
        par(mar=c(0,4.1,1,2.1))
        plot(density(diagnoses$t$train),xlim=range(unlist(diagnoses$t)),
             col=1,bty="n",xaxt="n",yaxt="n",xlab="",ylab="",main="")
        points(density(diagnoses$t$test),type="l",col="red")
        title(main = expression("Density of"~d[bold(t)](bold(t))),
              line = -2,cex.main=cex.text)
        ##
        par(mar=c(5.1,0,4.1,0))
        dd <- density(diagnoses$x$train)
        plot(dd$y,dd$x,type="l",bty="n",xaxt="n",yaxt="n",xlab="",ylab="",main="",
             ylim=range(unlist(diagnoses$x)),col=1)
        title(main = expression("Density of"~d[bold(x)](hat(bold(x)),bold(x))),
              line = -2,cex.main=cex.text)
        dd <- density(diagnoses$x$test)
        points(dd$y,dd$x,type="l",col="red")
        ##
        plot(0,0,col="white",bty="n",xaxt="n",yaxt="n",xlab="",ylab="")
      }
    }
  }
  list(y_est=y_est,
       diagnoses=diagnoses)
}
