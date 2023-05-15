#' Function to sum up bootstrap performance results of the ddsPLS algorithm
#'
#' @param x A ddsPLS object.
#' @param \ldots arguments to be passed to methods, such as graphical parameters.
#'
#' @examples
#' n <- 100 ; d <- 2 ; p <- 20 ; q <- 2
#' phi <- matrix(rnorm(n*d),n,d)
#' a <- rep(1,p/4) ; b <- rep(1,p/2)
#' X <- phi%*%matrix(c(1*a,0*a,0*b,1*a,3*b,0*a),nrow = d,byrow = TRUE) +
#' matrix(rnorm(n*p,sd = 1/4),n,p)
#' Y <- phi%*%matrix(c(1,0,0,0),nrow = d,byrow = TRUE) +
#' matrix(rnorm(n*q,sd = 1/4),n,q)
#' res <- ddsPLS(X,Y,verbose=FALSE)
#' print(res)
#'
#' @return  No return value, called for side effects
#'
#' @export
#' @rdname print.ddsPLS
#' @name print.ddsPLS
#'
#' @seealso \code{\link{ddsPLS}}, \code{\link{plot.ddsPLS}}, \code{\link{predict.ddsPLS}}
#'
#' @useDynLib ddsPLS
print.ddsPLS <- function(x,...)
{
  cat("\nCall:\n",
      paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  h_opt <- x$R
  if(h_opt>0){
    cat(paste("ddsPLS model built on",h_opt,"components.\n\n"))
  }else{
    cat(paste("No ddsPLS model built.\n\n"))
  }
  invisible(x)
}

#' Function to sum up bootstrap performance results of the ddsPLS algorithm
#'
#' @param object A ddsPLS object.
#' @param returnValues boolean. Wether or not to return the printed values, default to FALSE.
#' @param digits integer indicating the number of decimal places (round) to be used.
#' @param \ldots arguments to be passed to methods, such as graphical parameters.
#'
#' @examples
#' n <- 100 ; d <- 2 ; p <- 20 ; q <- 2
#' phi <- matrix(rnorm(n*d),n,d)
#' a <- rep(1,p/4) ; b <- rep(1,p/2)
#' X <- phi%*%matrix(c(1*a,0*a,0*b,1*a,3*b,0*a),nrow = d,byrow = TRUE) +
#' matrix(rnorm(n*p,sd = 1/4),n,p)
#' Y <- phi%*%matrix(c(1,0,0,0),nrow = d,byrow = TRUE) +
#' matrix(rnorm(n*q,sd = 1/4),n,q)
#' res <- ddsPLS(X,Y,verbose=FALSE)
#' summary(res,digits=5)
#'
#' @return  No return value, called for side effects
#'
#' @export
#' @rdname summary.ddsPLS
#' @name summary.ddsPLS
#'
#' @seealso \code{\link{ddsPLS}}, \code{\link{plot.ddsPLS}}, \code{\link{predict.ddsPLS}}
#'
#' @useDynLib ddsPLS
summary.ddsPLS <- function(object,returnValues=FALSE,
                           digits=2,...){
  cat("                      ______________\n");
  cat("                     |    ddsPLS    |\n");
  cat("=====================----------------=====================\n");
  h_opt <- object$R
  if(h_opt>0){
    q <- ncol(object$Y_obs)
    p <- nrow(object$model$U)
    if(object$doBoot){
      ## R2 and Q2 final values
      R2Q2 <- matrix(0,h_opt,5)
      for(h in 1:h_opt){
        best <- object$lambda[h]
        bestID <- which(object$results$lambdas==best)
        R2Q2[h,] <- c(
          best,
          object$results$R2mean[[h]][bestID],
          object$results$R2hmean[[h]][bestID],
          object$results$Q2mean[[h]][bestID],
          object$results$Q2hmean[[h]][bestID]
        )
      }
      colnames(R2Q2) <- c("lambda","R2","R2_r","Q2","Q2_r")
      rownames(R2Q2) <- paste("Comp.",1:h_opt)
    }
    ## Explained variances for the X part
    VARX <- do.call(rbind,object$varExplained_in_X)
    rownames(VARX) <- c("Per component ", "Cumulative ")
    colnames(VARX) <- paste("Comp.",1:h_opt)
    ## Explained variances for the components
    VAR <- do.call(rbind,object$varExplained[c("Comp","Cumu")])
    rownames(VAR) <- c("Per component ", "Cumulative ")
    colnames(VAR) <- paste("Comp.",1:h_opt)
    ## Explained variances for the response variables
    VARY <- object$varExplained$PerY
    if(is.null(colnames(object$Y_obs))){
      colnames(VARY) <- paste("Y",1:q,sep="")
    } else {
      colnames(VARY) <- colnames(object$Y_obs)
    }
    rownames(VARY) <- "All_comps"#paste("Comp.",1:h_opt)
    ## Explained variances per response variables per component marginally
    VARYCompMarg <- matrix(object$varExplained$PerYPerComp$Comp,ncol=q)
    if(is.null(colnames(object$Y_obs))){
      colnames(VARYCompMarg) <- paste("Y",1:q,sep="")
    } else {
      colnames(VARYCompMarg) <- colnames(object$Y_obs)
    }
    rownames(VARYCompMarg) <- paste("Comp.",1:h_opt)
    ## ... in total
    VARYCompTotal <- matrix(object$varExplained$PerYPerComp$Cum,ncol=q)
    if(is.null(colnames(object$Y_obs))){
      colnames(VARYCompTotal) <- paste("Y",1:q,sep="")
    } else {
      colnames(VARYCompTotal) <- colnames(object$Y_obs)
    }
    rownames(VARYCompTotal) <- paste("Comp.",1:h_opt)
    ######
    ######
    if(object$doBoot){
      cat(paste("The optimal ddsPLS model is built on",h_opt,"component(s)\n\n")
          )
    }else{
      cat(paste("The chosen ddsPLS model is built on",h_opt,"component(s)\n\n"))
    }
    ######
    if(object$doBoot){
      cat(paste("The bootstrap quality statistics:\n",
                "---------------------------------\n",sep=""))
      print(round(R2Q2,digits))
    }
    ######
    cat(paste("\n\nThe explained variance (in %):\n",
              "-----------------------\n",sep=""))
    cat(paste("\nIn total: ",round(object$varExplained$Cumu[h_opt],digits),
              "\n-  -  -  \n",sep=""))

    cat(paste("\nPer component or cumulated:\n",
              "-  -  -  -  -  -  -  -  -  \n",sep=""))
    print(round(VAR,digits))
    ######
    cat(paste("\nPer response variable:\n",
              "-  -  -  -  -  -  -  -\n",sep=""))
    print(round(VARY,digits))
    ######
    cat(paste("\nPer response variable per component:\n",
              "-  -  -  -  -  -  -  -  -  -  -  -  \n",
              sep=""))
    print(round(VARYCompMarg,digits))
    ######
    cat(paste("\n...and cumulated to:\n",
              "-  -  -  -  -  -  - \n",
              sep=""))
    print(round(VARYCompTotal,digits))

    cat(paste("\nFor the X block:\n",
              "-  -  -  -  -  -  -  -  -  \n",sep=""))
    print(round(VARX,digits))
  }else{
    cat(paste("The optimal ddsPLS model is empty.\n"))
  }
  cat("=====================                =====================\n");
  cat("                     ================\n");
  if(returnValues & h_opt>0){
    if(object$doBoot){
      out <- list(R2Q2=R2Q2,VAR=VAR,VARY=VARY,
                  VARYperComp=list(marginal=VARYCompMarg,total=VARYCompTotal))
    }else{
      out <- list(VAR=VAR,VARY=VARY,
                  VARYperComp=list(marginal=VARYCompMarg,total=VARYCompTotal))
    }
    return(out)
  }
}
