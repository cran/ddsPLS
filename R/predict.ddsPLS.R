#' Function to predict from ddsPLS objects
#'
#' @param object ddsPLS object.
#' @param X_test matrix, a test data-set. If is "NULL", the default value,
#' the predicted values for the train test are returned.
#' @param toPlot boolean, wether or not to plot the extreme value test plot.
#' Default to `TRUE`.
#' @param doDiagnosis boolean, wether or not to perform the diagnoses
#' operations. See \code{Value} section for more details.
#' @param legend.position character. Where to put the legend.
#' @param cex float positive. Number indicating the amount by which plotting
#' symbols should be scaled relative to the default.
#' @param cex.text float positive. Number indicating the amount by which
#' plotting text elements should be scaled relative to the default.
#' @param \ldots arguments to be passed to methods, such as graphical parameters.
#'
#' @details The diagnostic descriptors are usefull to detect potential outliers
#' in the \code{train} or in the \code{test} datasets. This can be appreciated
#' putting the parameter \code{toPlot} to \code{TRUE}. Thus a graph is created
#' projecting, the observations \eqn{i} of the \code{train} and the \code{test}
#' datasets, in the \bold{standardized} subspaces spanned by the \code{ddsPLS} model:
#'     \deqn{
#'     \left(
#'     \epsilon_x(i),\epsilon_t(i)
#'     \right)=
#'     \left(
#'         \dfrac{1}{\sqrt{p}}\sqrt{\sum_{j=1}^p \left(\dfrac{\left[\hat{\mathbf{x}}_i-
#'         \mathbf{x}_i\right]_{(j)}}{\hat{\sigma}_j^{(x)}}\right)^2},
#'         \dfrac{1}{\sqrt{R}}\sqrt{\sum_{r=1}^R \left(\dfrac{\hat{{t}}_i^{(r)}}
#'         {\hat{\sigma}_r}\right)^2}
#'     \right),
#'     }
#'     where \eqn{[\cdot]_{(j)}} takes the \eqn{j^{th}} coordinate of its argument.
#'     The different estimators are
#'     \deqn{
#'         \hat{{t}}_i^{(r)} = (\mathbf{x}_i-\hat{\boldsymbol{\mu}}_\mathbf{x})
#'         \mathbf{u}_r,
#'     }
#'     \deqn{
#'         \hat{\mathbf{x}}_i = \dfrac{1}{R}(\mathbf{x}_i-\hat{\boldsymbol{\mu}}_\mathbf{x})
#'         \sum_{r=1}^R\mathbf{u}_r\mathbf{p}_r^\top,
#'     }
#'     plus \eqn{\forall j\in[\![1,p]\!]}, \eqn{\hat{\sigma}_j^{(x)2}} is the
#'     estimated empirical variance of the \eqn{j^{th}} variable of \code{X}
#'     estimated on the \code{train} set of size \eqn{n}.
#'     Also, \eqn{\hat{\sigma}_r^2=\dfrac{1}{n}\sum_{j=1}^n\left((\mathbf{x}_j-\hat{\boldsymbol{\mu}}_\mathbf{x})
#'     \mathbf{u}_r\right)^2} is thus the estimated empirical variance of
#'     the \eqn{r^{th}}-component. Further, \eqn{R} is the approximated
#'     number of components of the \code{ddsPLS} model, \eqn{\hat{\boldsymbol{\mu}}_\mathbf{x}}
#'     is the empirical mean of \code{X}, \eqn{\mathbf{u}_r} is the weight of
#'     \code{X} along the \eqn{r^{th}} component, \eqn{\mathbf{p}_r} is the
#'     loading for \code{X}.
#'     The \code{diagnoses} object of the output is filled with two lists:
#'     \itemize{
#'        \item{\code{$object}}{ the first coordinate of the previous bivariate
#'        description corresponding to the reconstruction by the \code{ddsPLS}
#'        model.}
#'        \item{\code{$t}}{ the second coordinate of the previous bivariate
#'        description, corresponding to the score.}
#'     }
#'
#' @return A list of two objects:
#' \itemize{
#'     \item{\code{Y_est}}{ the estimated values for the response variable.}
#'     \item{\code{diagnoses}}{ the results of diagnostic operations, useful to detect
#'     potential outliers in the dataset.}}
#'
#' @export
#' @rdname predict.ddsPLS
#' @name predict.ddsPLS
#'
#' @examples
#' n <- 100 ; d <- 2 ; p <- 20 ; q <- 2 ; n_test <- 1000
#' phi <- matrix(rnorm(n*d),n,d)
#' phi_test <- matrix(rnorm(n_test*d),n_test,d)
#' a <- rep(1,p/4) ; b <- rep(1,p/2)
#' X <- phi%*%matrix(c(1*a,0*a,0*b,1*a,3*b,0*a),nrow = d,byrow = TRUE) +
#' matrix(rnorm(n*p,sd = 1/4),n,p)
#' X_test <- phi_test%*%matrix(c(1*a,0*a,0*b,1*a,3*b,0*a),nrow = d,byrow=TRUE) +
#' matrix(rnorm(n_test*p,sd = 1/4),n_test,p)
#' Y <- phi%*%matrix(c(1,0,0,0),nrow = d,byrow = TRUE) +
#' matrix(rnorm(n*q,sd = 1/4),n,q)
#' res <- ddsPLS(X,Y,verbose=FALSE)
#' pre <- predict(res,X_test = X_test,toPlot = TRUE,doDiagnosis = TRUE)
#'
#' @seealso \code{\link{ddsPLS}}, \code{\link{plot.ddsPLS}}, \code{\link{summary.ddsPLS}}
#'
#' @useDynLib ddsPLS
predict.ddsPLS <- function(object,X_test=NULL,toPlot=FALSE,doDiagnosis=TRUE,
                           legend.position="topright",cex=1,cex.text=1,...){
  getDiagnoses <- function(object,X_test,y_test_est){
    n_test <- nrow(X_test)
    n_train <- nrow(object$X)
    X_test_center <- X_test-matrix(rep(object$model$muX,n_test),nrow = n_test,byrow = TRUE)
    X_train_center <- object$X-matrix(rep(object$model$muX,n_train),nrow = n_train,byrow = TRUE)
    ## On the predictions
    # y_train_est <- object$Y_est
    # sigma_y_inv <- unlist(lapply(object$model$sdY,function(ss){
    #   if(ss>1e-9){out <- 1/ss}else{out <- 0};out}))
    # means_y_train_est <- colMeans(y_train_est)
    # y_mean_train <- matrix(rep(means_y_train_est,n_train),nrow = n_train,
    #                        byrow = TRUE)
    # y_mean_test <- matrix(rep(means_y_train_est,n_test),nrow = n_test,
    #                       byrow = TRUE)
    # dist_y_Mahalanobis_train <- (rowMeans(
    #   (y_train_est-y_mean_train)^2*matrix(rep(sigma_y_inv,n_train),
    #                                       nrow = n_train,byrow = TRUE)^2))
    # dist_y_Mahalanobis_test <- (rowMeans(
    #   (y_test_est-y_mean_test)^2*matrix(rep(sigma_y_inv,n_test),
    #                                     nrow = n_test,byrow = TRUE)^2))
    ## On the scores
    t_test <- X_test_center%*%object$model$U_star
    t_train <- X_train_center%*%object$model$U_star
    sigma_t_inv_2 <- unlist(lapply((colMeans(t_train^2)),function(ss){
      if(ss>1e-9){out <- 1/ss}else{out <- 0};out}))
    dist_score_Mahalanobis_train <- sqrt(rowMeans(
      t_train^2*matrix(rep(sigma_t_inv_2,n_train),nrow = n_train,byrow = TRUE)^2))
    dist_score_Mahalanobis_test <- sqrt(rowMeans(
      t_test^2*matrix(rep(sigma_t_inv_2,n_test),nrow = n_test,byrow = TRUE)^2))
    ## On the reconstruction of X
    X_test_center_est <- tcrossprod(t_test,object$model$P)
    X_train_center_est <- tcrossprod(t_train,object$model$P)
    epsilon_X_test <- X_test_center_est-X_test_center
    epsilon_X_train <- X_train_center_est-X_train_center
    sigma_x_Inv_2 <- unlist(lapply((colMeans(epsilon_X_train^2)),function(ss){
      if(ss>1e-9){out <- 1/ss}else{out <- 0};out}))
    dist_Mahalanobis_train <- sqrt(rowMeans(
      epsilon_X_train^2*matrix(rep(sigma_x_Inv_2,n_train),nrow = n_train,byrow = TRUE)^2))
    dist_Mahalanobis_test <- sqrt(rowMeans(
      epsilon_X_test^2*matrix(rep(sigma_x_Inv_2,n_test),nrow = n_test,byrow = TRUE)^2))
    list(t=list(train=dist_score_Mahalanobis_train,test=dist_score_Mahalanobis_test),
         x=list(train=dist_Mahalanobis_train,test=dist_Mahalanobis_test))
  }
  diagnoses <- NULL
  if(is.null(X_test)){
    X_test <- object$X
  }
  n_test <- nrow(X_test)
  if(object$R==0){
    Y_est <- matrix(rep(object$model$muY,n_test),nrow = n_test,byrow = TRUE)
  }else{
    Y_est <- (X_test-matrix(rep(object$model$muX,n_test),nrow = n_test,byrow = TRUE))%*%object$model$B
    Y_est <- Y_est + matrix(rep(object$model$muY,n_test),nrow = n_test,byrow = TRUE)
  }
  if(doDiagnosis & object$R>0){
    diagnoses <- getDiagnoses(object,X_test,Y_est)
    if(toPlot){
      #####################################
      oldpar <- par(no.readonly = TRUE)
      on.exit(par(oldpar))
      #####################################
      xlim <- range(c(diagnoses$t$train,diagnoses$t$test))
      ylim <- range(c(diagnoses$x$train,diagnoses$x$test))
      N <- 4
      par(mar=c(5, 4, 1, 0.5) + 0.1)
      layout(matrix(c(rep(2,N),4,rep(c(rep(1,N),3),N)),nrow = N+1,byrow = TRUE))
      plot(diagnoses$t$train,diagnoses$x$train,col=1,pch=1,
           cex=cex,xlab="",ylab="",
           xlim=xlim,ylim=ylim)
      title(xlab = expression(epsilon[t]),line = 2.5,cex.main=cex.text)
      title(ylab = expression(epsilon[x]),line = 2.5,cex.main=cex.text)
      points(diagnoses$t$test,diagnoses$x$test,col="red",pch=16,cex=cex)
      legend(legend.position,c("Train","Test"),
             col=c(1,2),pch=c(1,16),cex=cex.text,pt.cex = cex)
      ##
      par(mar=c(0,4.1,1,2.1))
      plot(density(diagnoses$t$train),xlim=range(unlist(diagnoses$t)),
           col=1,bty="n",xaxt="n",yaxt="n",xlab="",ylab="",main="")
      points(density(diagnoses$t$test),type="l",col="red")
      # title(main = expression("Kern. dens. est. of"~epsilon[t]),
      #       line = -2,cex.main=cex.text)
      ##
      par(mar=c(5.1,0,4.1,0))
      dd <- density(diagnoses$x$train)
      plot(dd$y,dd$x,type="l",bty="n",xaxt="n",yaxt="n",xlab="",ylab="",main="",
           ylim=range(unlist(diagnoses$x)),col=1)
      # title(main = expression("Kern. dens. est. of"~epsilon[x]),
      #       line = -2,cex.main=cex.text,srt=45)
      dd <- density(diagnoses$x$test)
      points(dd$y,dd$x,type="l",col="red")
      ##
      plot(0,0,col="white",bty="n",xaxt="n",yaxt="n",xlab="",ylab="")
    }
  }
  list(Y_est=Y_est,
       diagnoses=diagnoses)
}
