#' Function to plot bootstrap performance results of the ddsPLS algorithm
#'
#' @param x A ddsPLS object
#' @param type The type of graphics. One of "criterion" (default),
#' "prop", "predict", "Q2r", "Q2", "R2r", "R2", "weightX", "weightY",
#'  "loadingX" or "loadingY". The type "prop" corresponds to the proportion of
#'  models with positive \eqn{Q_{r}^2} among the bootstrapped evaluations for
#'  each of the to be tested values for \code{lambda}.
#' @param digits double. Rounding of the written explained variance.
#' @param legend.position character. Where to put the legend.
#' @param horiz boolean. Whether to plot horizontally.
#' @param biPlot boolean. Whether to plot one component versus another one.
#' @param comp vector of two integers. Taken into account only if
#' \code{biplot=TRUE}.
#' @param col vector. Mainly to modify bars in weight plots.
#' @param mar vector. The margins for the plot.
#' @param cex.names double. Size factor for variable names.
#' @param \ldots arguments to be passed to methods, such as graphical parameters.
#'
#' @return No return value, called for side effects.
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
#'
#' ## Plot criterion
#' plot(res,type = "criterion")
#' ## Univariate weights of X, same with Y
#' plot(res,type = "weightX")
#' ## Bivariate weights of X, same with Y
#' plot(res,type = "weightX",biPlot = TRUE)
#' ## Modify margins to fit the window
#' plot(res,type = "weightY",mar = c(3,7,3,3))
#'
#' @importFrom graphics layout
#'
#' @export
#' @rdname plot.ddsPLS
#' @name plot.ddsPLS
#'
#' @seealso \code{\link{ddsPLS}}, \code{\link{predict.ddsPLS}}, \code{\link{summary.ddsPLS}}
#'
#' @useDynLib ddsPLS
plot.ddsPLS <- function(x,type="criterion",
                        digits=1,
                        legend.position="topright",
                        horiz=TRUE,biPlot=FALSE,comp=c(1,2),
                        col=NULL,
                        cex.names=1,mar=c(5, 4, 4, 2) + 0.1,...){
  ## Reset personnal plot par() settings
  opar <- par(no.readonly =TRUE)
  on.exit(par(opar))
  ## -----------------------------------
  h_opt <- x$R
  got_inside <- FALSE
  if(h_opt>0){
    q <- ncol(x$Y_obs)
    p <- nrow(x$model$U)
    lambdas <- x$results$lambdas
    par(mar=mar)
    if(x$doBoot){
      if(h_opt>1){
        R2mean_diff_Q2mean <- matrix(do.call(cbind,x$results$R2mean_diff_Q2mean)[,1:h_opt],ncol = h_opt)
        Q2hmean <- matrix(do.call(cbind,x$results$Q2hmean)[,1:h_opt],ncol = h_opt)
        Q2mean <- matrix(do.call(cbind,x$results$Q2mean)[,1:h_opt],ncol = h_opt)
        R2hmean <- matrix(do.call(cbind,x$results$R2hmean)[,1:h_opt],ncol = h_opt)
        R2mean <- matrix(do.call(cbind,x$results$R2mean)[,1:h_opt],ncol = h_opt)
        PropQ2hPos <- matrix(do.call(cbind,x$results$PropQ2hPos)[,1:h_opt],ncol = h_opt)
      }else{
        R2mean_diff_Q2mean <- matrix(x$results$R2mean_diff_Q2mean[[1]],ncol = 1)
        Q2hmean <- matrix(x$results$Q2hmean[[1]],ncol = 1)
        Q2mean <- matrix(x$results$Q2mean[[1]],ncol = 1)
        R2hmean <- matrix(x$results$R2hmean[[1]],ncol = 1)
        R2mean <- matrix(x$results$R2mean[[1]],ncol = 1)
        PropQ2hPos <- matrix(x$results$PropQ2hPos[[1]],ncol = 1)
      }
      if(type %in% c("criterion")){
        if(x$criterion=="diffR2Q2"){
          # Plot of R2-Q2
          matplot(lambdas,R2mean_diff_Q2mean,type = "l",ylab="",xlab=expression(lambda),
                  main=bquote(bar(R)[B]^2-bar(Q)[B]^2))
          for(s in 1:h_opt){
            points(lambdas,R2mean_diff_Q2mean[,s],type = "p",pch=16,cex=x$lambda_optim[[s]],col=s)
            points(x$lambda[s],R2mean_diff_Q2mean[which(lambdas==x$lambda[s]),s],pch=1,cex=2,col=s)
          }
          legend(legend.position,paste("Comp.",1:h_opt," (",round(x$varExplained$Comp),"%)",sep=""),
                 col = 1:h_opt,pch=16,bty = "n",
                 title = paste("Total explained variance ",round(x$varExplained$Cumu)[h_opt],"%",sep=""))
        }else{
          # Plot of Q2
          matplot(lambdas,Q2hmean,type = "l",ylab="",xlab=expression(lambda),
                  main=bquote(bar(Q)[B]^2))
          for(s in 1:h_opt){
            points(lambdas,Q2hmean[,s],type = "p",pch=16,cex=x$lambda_optim[[s]],col=s)
            points(x$lambda[s],Q2hmean[which(lambdas==x$lambda[s]),s],pch=1,cex=2,col=s)
          }
          legend(legend.position,paste("Comp.",1:h_opt," (",round(x$varExplained$Comp),"%)",sep=""),
                 col = 1:h_opt,pch=16,bty = "n",
                 title = paste("Total explained variance ",round(x$varExplained$Cumu)[h_opt],"%",sep=""))
        }
        got_inside <- TRUE
      }
      if(type %in% c("prop")){
        # Plot of Prop of positive Q2h
        matplot(lambdas,PropQ2hPos,type = "l",ylab="",xlab=expression(lambda),
                main=bquote("Proportion of models with positive"~Q["b,r"]^2))
        abline(h=((1:10)/10)[-5],col="gray80",lwd=0.5,lty=3)
        abline(h=5/10,col="gray60",lwd=0.7,lty=1)
        text(min(lambdas),1/2,labels = "1/2",pos = 4,col="gray40")
        for(s in 1:h_opt){
          points(lambdas,PropQ2hPos[,s],type = "p",pch=16,cex=x$lambda_optim[[s]],col=s)
          points(x$lambda[s],PropQ2hPos[which(lambdas==x$lambda[s]),s],pch=1,cex=2,col=s)
        }
        legend(legend.position,paste("Comp.",1:h_opt," (",round(x$varExplained$Comp),"%)",sep=""),
               col = 1:h_opt,pch=16,bty = "n",
               title = paste("Total explained variance ",round(x$varExplained$Cumu)[h_opt],"%",sep=""))
        q <- ncol(x$Y_obs)
        ;got_inside <- TRUE}
    }
    if(type %in% c("predict")){
      # Predicted versus observed
      matplot(x$Y_obs,x$Y_est,pch=1:q,type="p",xlab="Observed",ylab="Predicted",
              main="Predicted versus observed")
      abline(0,1)
      if(is.null(colnames(x$Y_obs))){
        nono <- paste("Y",1:q," (",round(x$varExplained$PerY,digits = digits),"%)",sep="")
      } else {
        nono <- paste(colnames(x$Y_obs)," (",round(x$varExplained$PerY,digits = digits),"%)",sep="")
      }
      legend(legend.position,nono,col = 1:q,pch=1:q,bty = "n")
      ;got_inside <- TRUE}
    if(type %in% c("selection")){
      allX <- rep(0,p)
      allY <- rep(0,q)
      allX[x$Selection$X] <- 1
      allY[x$Selection$Y] <- 1
      if(is.null(rownames(x$model$U))){
        names(allX) <- paste("X",1:p,sep="")
      }
      plot(allX,ylab="",
           main=paste("Which variables are selected in block X ?"),
           yaxt="n",xaxt="n",ylim=c(0,1),col=1+allX,pch=16+allX)
      axis(2,at = c(0,1),labels = c("No","Yes"),las=2)
      axis(1,at = 1:length(allX),labels = names(allX),
           las=2,cex.axis=cex.names)
      p_app <- 5*((p-p%%5)/5+1)
      ltys <- rep(2,p_app)
      ltys[5*(0:(p_app/5))] <- 1
      abline(v=1:p_app,col="gray",lwd=1/2,lty=ltys)
      if(is.null(rownames(x$model$V))){
        names(allY) <- paste("Y",1:q,sep="")
      }
      selY <- length(x$Selection$Y)
      plot(allY,ylab="",
           main=paste("Which variables are selected in block Y ?"),
           yaxt="n",xaxt="n",ylim=c(0,1),col=1+allY,pch=16+allY)
      axis(2,at = c(0,1),labels = c("No","Yes"),las=2)
      axis(1,at = 1:length(selY),labels = names(selY),
           las=2,cex.axis=cex.names)
      q_app <- 5*((q-q%%5)/5+1)
      ltys <- rep(2,q_app)
      ltys[5*(0:(q_app/5))] <- 1
      abline(v=1:q_app,col="gray",lwd=1/2,lty=ltys)
      ;got_inside <- TRUE}
    if(x$doBoot){
      if(type %in% c("Q2r")){
        # PLot of Q2_h
        matplot(lambdas,Q2hmean,type = "l",ylab="",xlab=expression(lambda),
                main=bquote(bar(Q)["B,r"]^2))
        for(s in 1:h_opt){  for(s in 1:h_opt){
          points(lambdas,Q2hmean[,s],type = "p",pch=16,cex=x$lambda_optim[[s]],col=s)
          points(x$lambda[s],Q2hmean[which(lambdas==x$lambda[s]),s],pch=1,cex=2,col=s)
        }
          points(x$lambda[s],Q2hmean[which(lambdas==x$lambda[s]),s],pch=16,col=s)
        }
        abline(h=x$lowQ2,lwd=2,lty=2)
        legend(legend.position,paste("Comp.",1:h_opt," (",round(x$varExplained$Comp),"%)",sep=""),
               col = 1:h_opt,pch=16,bty = "n",
               title = paste("Total explained variance ",round(x$varExplained$Cumu)[h_opt],"%",sep=""))
        ;got_inside <- TRUE}
      if(type %in% c("R2r")){
        # PLot of R2_h
        matplot(lambdas,R2hmean,type = "l",ylab="",xlab=expression(lambda),
                main=bquote(bar(R)["B,r"]^2))
        for(s in 1:h_opt){  for(s in 1:h_opt){
          points(lambdas,R2hmean[,s],type = "p",pch=16,cex=x$lambda_optim[[s]],col=s)
          points(x$lambda[s],R2hmean[which(lambdas==x$lambda[s]),s],pch=1,cex=2,col=s)
        }
          points(x$lambda[s],R2hmean[which(lambdas==x$lambda[s]),s],pch=16,col=s)
        }
        legend(legend.position,paste("Comp.",1:h_opt," (",round(x$varExplained$Comp),"%)",sep=""),
               col = 1:h_opt,pch=16,bty = "n",
               title = paste("Total explained variance ",round(x$varExplained$Cumu)[h_opt],"%",sep=""))
        ;got_inside <- TRUE}
      if(type %in% c("Q2")){
        # PLot of Q2
        matplot(lambdas,Q2mean,type = "l",ylab="",xlab=expression(lambda),
                main=bquote(bar(Q)["B"]^2))
        for(s in 1:h_opt){  for(s in 1:h_opt){
          points(lambdas,Q2mean[,s],type = "p",pch=16,cex=x$lambda_optim[[s]],col=s)
          points(x$lambda[s],Q2mean[which(lambdas==x$lambda[s]),s],pch=1,cex=2,col=s)
        }
          points(x$lambda[s],Q2mean[which(lambdas==x$lambda[s]),s],pch=16,col=s)
        }
        abline(h=0,lwd=2,lty=2)
        legend(legend.position,paste("Comp.",1:h_opt," (",round(x$varExplained$Comp),"%)",sep=""),
               col = 1:h_opt,pch=16,bty = "n",
               title = paste("Total explained variance ",round(x$varExplained$Cumu)[h_opt],"%",sep=""))
        ;got_inside <- TRUE}
      if(type %in% c("R2")){
        # PLot of R2
        matplot(lambdas,R2mean,type = "l",ylab="",xlab=expression(lambda),
                main=bquote(bar(R)["B"]^2))
        for(s in 1:h_opt){  for(s in 1:h_opt){
          points(lambdas,R2mean[,s],type = "p",pch=16,cex=x$lambda_optim[[s]],col=s)
          points(x$lambda[s],R2mean[which(lambdas==x$lambda[s]),s],pch=1,cex=2,col=s)
        }
          points(x$lambda[s],R2mean[which(lambdas==x$lambda[s]),s],pch=16,col=s)
        }
        abline(h=0,lwd=2,lty=2)
        legend(legend.position,paste("Comp.",1:h_opt," (",round(x$varExplained$Comp),"%)",sep=""),
               col = 1:h_opt,pch=16,bty = "n",
               title = paste("Total explained variance ",round(x$varExplained$Cumu)[h_opt],"%",sep=""))
        ;got_inside <- TRUE}
    }
    if(type == "weightX" | type == "weightY"){
      got_inside <- TRUE
      q <- ncol(x$Y_obs)
      if(is.null(colnames(x$Y_obs))){
        colnames(x$Y_obs) <- paste("Y",1:q,sep="")
      }
      if(is.null(rownames(x$model$U))){
        p <- nrow(x$model$U)
        rownames(x$model$U) <- paste("X",1:p,sep="")
      }
      if(!biPlot){
        par(mfrow=c(1,h_opt),mar=mar)
        for(s in 1:h_opt){
          if(type == "weightX"){
            popo <- t(x$model$U)[s,,drop=FALSE]
            if(is.null(col)){
              col <- rep(1,q)
            }
            colo <- col
            barplot(as.vector(popo),
                    xlim=c(-1,1)*max(abs(popo)),
                    horiz = horiz,axis.lty = 1,las=2,
                    names.arg=colnames(popo),
                    main=paste("X part, Comp.",s),
                    col=colo,border=colo,
                    cex.names = cex.names)
          }
          else{
            popo <- t(x$model$V)[s,,drop=FALSE]
            if(is.null(colnames(popo))){
              colnames(popo) <- paste("Y",
                                      1:q,
                                      sep="")
            }
            colnames(popo) <- paste(
              colnames(popo)," (",
              round(x$varExplained$PerYPerComp$Comp[s,],
                    digits = digits ),
              "%)",sep="")
            if(is.null(col)){
              col <- 1:q
            }
            barplot(as.vector(popo),
                    xlim=c(-1,1)*max(abs(popo)),
                    cex.names = cex.names,horiz = horiz,
                    names=colnames(popo),col=col,border=col,
                    axis.lty = 1,las=2,
                    main=paste("Y part, Comp.",s))
          }
          abline(v = c(-10:10)/10,lty=2,col="gray")
        }
      }else{
        if(type == "weightX" | h_opt>1){
          plot(x$model$U[,comp[1]],x$model$U[,comp[2]],
               main="Weights for X",xlim=c(-1,1),ylim=c(-1,1),
               xlab=paste("Comp. ",comp[1]," (",round(x$varExplained_in_X$Comp[comp[1]]),"%)",sep="" ),
               ylab=paste("Comp. ",comp[2]," (",round(x$varExplained_in_X$Comp[comp[2]]),"%)",sep="" ))
          abline(h=0,v=0,lty=2,col="gray")
          thetas <- seq(0,2*pi,length.out=1e3)
          points(cos(thetas),sin(thetas),type="l",lty=2,col="gray")
        }
      }
    }
    if(type == "loadingX" | type == "loadingY"){
      got_inside <- TRUE
      q <- ncol(x$Y_obs)
      if(is.null(colnames(x$Y_obs))){
        colnames(x$Y_obs) <- paste("Y",1:q,sep="")
      }
      if(is.null(rownames(x$model$P))){
        p <- nrow(x$model$P)
        rownames(x$model$P) <- paste("X",1:p,sep="")
      }
      if(!biPlot){
        par(mfrow=c(1,h_opt),mar=mar)
        for(s in 1:h_opt){
          if(type == "loadingX"){
            popo <- t(x$model$P)[s,,drop=FALSE]
            if(is.null(col)){
              col <- rep(1,q)
            }
            colo <- col
            barplot(as.vector(popo),
                    xlim=c(-1,1)*max(abs(popo)),
                    horiz = horiz,axis.lty = 1,las=2,
                    names.arg=colnames(popo),
                    main=paste("X part, Comp.",s),
                    col=colo,border=colo,
                    cex.names = cex.names)
          }else{
            popo <- t(x$model$C)[s,,drop=FALSE]
            if(is.null(colnames(popo))){
              colnames(popo) <- paste("Y",
                                      1:q,
                                      sep="")
            }
            colnames(popo) <- paste(
              colnames(popo)," (",
              round(x$varExplained$PerYPerComp$Comp[s,],
                    digits = digits ),
              "%)",sep="")
            if(is.null(col)){
              col <- 1:q
            }
            barplot(as.vector(popo),
                    xlim=c(-1,1)*max(abs(popo)),
                    cex.names = cex.names,horiz = horiz,
                    names=colnames(popo),col=col,border=col,
                    axis.lty = 1,las=2,
                    main=paste("Y part, Comp.",s))
          }
          abline(v = c(-10:10)/10,lty=2,col="gray")
        }
      }else{
        if(type == "loadingX" | h_opt>1){
          loading1 <- x$model$P[,comp[1]]/sqrt(sum(x$model$P[,comp[1]]^2))
          loading2 <- x$model$P[,comp[2]]/sqrt(sum(x$model$P[,comp[2]]^2))
          plot(loading1,loading2,xlim=c(-1,1),ylim=c(-1,1),
               main="Normalized loadings for X",
               xlab=paste("Comp. ",comp[1]," (",round(x$varExplained_in_X$Comp[comp[1]]),"%)",sep="" ),
               ylab=paste("Comp. ",comp[2]," (",round(x$varExplained_in_X$Comp[comp[2]]),"%)",sep="" ))
          abline(h=0,v=0,lty=2,col="gray")
          thetas <- seq(0,2*pi,length.out=1e3)
          points(cos(thetas),sin(thetas),type="l",lty=2,col="gray")
        }
      }

    }
    if(!got_inside){
      cat(
        "Please select a correct type of vizu among `criterion` (default),
      `prop`, `predict`, `Q2r`, `Q2`, `R2r`, `R2`, `weightX`, `weightY`,
      `loadingX` or `loadingY`.")
    }
  }
}
