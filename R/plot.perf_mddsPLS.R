#' Function to plot cross-validation performance results.
#'
#' That function must be applied to a perf_mddsPLS object. Extra parameters are
#'  avalaible to control the plot quality.
#'
#' @param x The perf_mddsPLS object.
#' @param plot_mean logical. Whether or not to plot the mean curve.
#' @param reg_error character. One of "MSEP" (Mean Squared Error in Prediction) or "MPE" (Mean Prediction Error). Default is "MSEP".
#' @param pos_legend character. One of "bottomleft", "topright",....
#' @param legend_names vector of characters. Each element is the name of one of the q response variables.
#' @param which_sd_plot vector of integers of length the number of columns in Y. Indicates which area of standard error must be drawn.
#' @param alpha.f factor modifying the opacity alpha; typically in [0,1]. Used by \strong{adjustcolor}
#' @param ylim numeric vectors of length 2, giving the error plot range.
#' @param no_occurence logical. Whether or not to plot the occurence plot of the \strong{Y} variables. Initialized to \strong{TRUE}.
#' @param no_plot logical. Whether or not to plot error. Initialized to \strong{FALSE}.
#' @param main character of \strong{NULL}. If null the title is given to the willing of the software. If \strong{""}, no title is given. Else is what the user wants.
#' @param ... Other plotting parameters to affect the plot.
#'
#' @return The plot visualisation
#'
#' @seealso  \code{\link{perf_mddsPLS}}, \code{\link{summary.perf_mddsPLS}}
#'
#' @importFrom graphics par axis matplot mtext points
#' @importFrom stats cor
#'
#' @export
#'
#' @examples
#' library(doParallel)
#' # Classification example :
#' data("penicilliumYES")
#' X <- penicilliumYES$X
#' X <- scale(X[,which(apply(X,2,sd)>0)])
#' Y <- as.factor(unlist(lapply(c("Melanoconidiu","Polonicum","Venetum"),
#' function(tt){rep(tt,12)})))
#' #res_cv_class <- perf_mddsPLS(X,Y,L0s=1:5,R = 2,
#' #mode = "lda",NCORES = 1,fold_fixed = rep(1:12,3))
#' #plot(res_cv_class)
#'
#' # Regression example :
#' data("liverToxicity")
#' X <- scale(liverToxicity$gene)
#' Y <- scale(liverToxicity$clinic)
#' #res_cv_reg <- perf_mddsPLS(Xs = X,Y = Y,L0s=c(1,5,10,15,20),R = 1,
#' # mode = "reg")
#' #plot(res_cv_reg)
plot.perf_mddsPLS <- function(x,plot_mean=FALSE,
                              reg_error="MSEP",
                              legend_names=NULL,
                              pos_legend="bottomleft",
                              which_sd_plot=NULL,
                              ylim=NULL,alpha.f=0.4,
                              no_occurence=T,
                              main=NULL,
                              no_plot=F,
                              ...){
  if(!no_plot){
    ## Reset personnal plot par() settings
    opar <- par(no.readonly =TRUE)
    on.exit(par(opar))
    ## -----------------------------------
  }
  res_perf_mdd <- x
  is_L0 <- names(res_perf_mdd[[1]])[2]
  X_all <- scale(do.call(cbind,res_perf_mdd$Xs))
  if(res_perf_mdd$mode=="reg"){
    if(reg_error=="MPE"){
      names(res_perf_mdd)[1] <- "OUT"
      names(res_perf_mdd)[4] <- "RMSEP"
      ylab1<-"MPE"
      main1 <- "MPE versus regularization coefficient mdd-sPLS"
      legend_0 <- "Mean MPE"
    }else{
      names(res_perf_mdd)[1] <- "RMSEP"
      ylab1<-"MSEP"
      main1 <- "MSEP versus regularization coefficient mdd-sPLS"
      legend_0 <- "Mean MSEP"
    }
    cc <- matrix(NA,nrow = ncol(res_perf_mdd$Y),ncol = ncol(X_all))
    for(j in 1:ncol(X_all)){
      cc[,j] <- abs(cor(res_perf_mdd$Y,X_all[,j],use = "pairwise.complete.obs"))
    }
    # cc <- abs(crossprod(scale(res_perf_mdd$Y),X_all)/(nrow(res_perf_mdd$Y)-1))
    col_na <- which(is.na(colSums(cc)))
    if(length(col_na)>0){
      cc <- cc[,-col_na,drop=F]
    }
  }else{
    names(res_perf_mdd)[1] <- "RMSEP"
    Y_df <- data.frame(res_perf_mdd$Y)
    colnames(Y_df) <- "Y"
    Y <- scale(model.matrix( ~ Y - 1, data=Y_df))
    cc <- abs(crossprod(Y,X_all)/(nrow(Y)-1))
  }
  ranges_0 <- sort(apply(cc,2,max))
  l_lambdas <- length(unique(res_perf_mdd$RMSEP[,2]))
  if(is_L0!="L0s"){
    xlab <- expression(lambda)
    if(ncol(cc)>1){
      if(l_lambdas>1){
        # ranges <- sort(ranges_0)
        l_in_min <- ranges_0>=min(res_perf_mdd$RMSEP[,2])
        l_in_max <- ranges_0<=max(res_perf_mdd$RMSEP[,2])
        inter_in <- l_in_min&l_in_max
        l_in <- length(ranges_0)-which(inter_in)
        ranges <- ranges_0[length(ranges_0)-l_in]
        card_ranges <- l_in
      }else{
        l_in_min <- which(ranges_0>=min(res_perf_mdd$RMSEP[,2]))
        ranges <- ranges_0[l_in_min]
        card_ranges <- rev(0:(length(ranges_0)-1))
      }
    }else{
      card_ranges <- 1
    }
  }else{
    xlab <- expression(L[0])
  }
  ERRORS <- res_perf_mdd
  FREQ <- ERRORS$FREQ
  RMSEP <- ERRORS$RMSEP
  SDEP <- ERRORS$SDEP
  q <- ncol(ERRORS$RMSEP)-2
  if(q<3){
    colors <- 1:q
  }
  else if(q>8){
    colors <- brewer.pal(8, "Dark2")
    pal <- colorRampPalette(colors)
    colors <- pal(q)
  }
  else{
    colors <- brewer.pal(q, "Dark2")
  }
  if(res_perf_mdd$mod=="reg"){
    ylab2<-"Occurences per variable (%)"
    ylim1 <- range(abs(RMSEP[,3:ncol(RMSEP)]))^2
    y1 <- RMSEP[order(RMSEP[,2,drop=FALSE]),3:ncol(RMSEP),drop=FALSE]^2
    y_mean <- rowMeans(RMSEP[order(RMSEP[,2,drop=FALSE]),3:ncol(RMSEP),drop=FALSE]^2)
    if(is.null(main)){
      main2 <- "Occurences per variable versus regularization coefficient mdd-sPLS"
    }else{
      main1 <- main
      main2 <- NULL
    }
    if(!no_occurence){
      par(mfrow=c(2,1),ann=T)
    }else{
      par(mfrow=c(1,1),ann=T)
    }
  }
  else{
    ylab1<-"#Good Classif Rate"
    ylab2<-"Occurences per class"
    ylim1<- c(0,1)
    y1 <- RMSEP[order(RMSEP[,2,drop=FALSE]),3:ncol(RMSEP),drop=FALSE]
    TAB <- table(res_perf_mdd$Y)
    for(r in 1:nlevels(res_perf_mdd$Y)){
      r_here <- which(names(TAB)==colnames(y1)[r])
      y1[,r] <- 1-y1[,r]/TAB[r_here]
    }
    y_mean <- 1-rowSums(RMSEP[order(RMSEP[,2,drop=FALSE]),3:ncol(RMSEP),drop=FALSE])/sum(TAB)
    main1 <- "Good classification rate versus regularization coefficient mdd-sPLS"
    main2 <- "Occurences per class versus regularization coefficient mdd-sPLS"
    if(!no_plot)par(mfrow=c(1,1),ann=T)
  }
  # graphics::matplot(sort(RMSEP[,2]),y1,type="l",lwd=4,lty=1,
  #                   ylim=ylim1,col=colors,
  #                   xlab=expression(lambda),
  #                   ylab=ylab1,
  #                   main=main1)
  lam_plot <- sort(RMSEP[,2])
  ord <- order(RMSEP[,2])
  delta <- rep(0,q)
  if(!is.null(which_sd_plot)){
    delta[which_sd_plot] <- 1
  }
  DELTAS <- matrix(rep(delta,nrow(y1)),ncol=length(delta),byrow = T)
  if(is.null(ylim)){
    min_y <- min(y1-SDEP[,-c(1:2)]*DELTAS)
    max_y <- max(y1+SDEP[,-c(1:2)]*DELTAS)
    if(is.na(max_y)|is.na(min_y)){
      ylim <- c(0,1)
    }else{
      ylim <- c(min_y,max_y)
    }
  }
  for(jq in 1:q){
    dat <- data.frame(list(lambda=lam_plot,MSEP=y1[ord,jq],sd=SDEP[ord,2+jq]))
    ses <- dat$MSEP + outer(dat$sd, c(1,-1)*delta[jq])
    if(jq==1){
      if(res_perf_mdd$mod=="reg"){
        if(!no_plot)with(dat,
                         plot(
                           lambda, MSEP, type="l",xlab=xlab,ylab=ylab1,
                           ylim=ylim,col=colors[jq],lwd=3,
                           main=main1,
                           panel.first=polygon(c(lambda,rev(lambda)), c(ses[,1],rev(ses[,2])),
                                               border=NA,
                                               col=adjustcolor(colors[jq],alpha.f = alpha.f))
                         )
        )
      }else{
        if(!no_plot)with(dat,
                         plot(
                           lambda, MSEP, type="l",xlab=xlab,ylab=ylab1,
                           ylim=ylim,col=colors[jq],lwd=3,
                           main=main1,xaxt="n",
                           panel.first=polygon(c(lambda,rev(lambda)), c(ses[,1],rev(ses[,2])),
                                               border=NA,
                                               col=adjustcolor(colors[jq],alpha.f = alpha.f))
                         )
        )
        if(!no_plot)axis(1,at=dat$lambda)
      }
    }else{
      if(!no_plot)with(dat,
                       points(
                         lambda, MSEP, type="l",col=colors[jq],lwd=3,
                         panel.first=polygon(c(lambda,rev(lambda)), c(ses[,1],rev(ses[,2])),
                                             border=NA,
                                             col=adjustcolor(colors[jq],alpha.f = alpha.f))
                       )
      )
    }
    if(delta[jq]==1){
      for(jj in 1:nrow(dat)){
        x <- dat$lambda
        CI.dn <- ses[,2]
        CI.up <- ses[,1]
        if(!no_plot)arrows(x,CI.dn,x,CI.up,code=3,length=0.05,angle=90,col=colors[jq])
      }
    }
  }
  if(res_perf_mdd$mod!="reg"){
    if(is_L0!="L0s"){
      pos_all <- max(which(y_mean==max(y_mean)))
      pos_one <- max(which(y1==max(y1),arr.ind = T)[,1])
    }else{
      pos_all <- min(which(y_mean==max(y_mean)))
      pos_one <- min(which(y1==max(y1),arr.ind = T)[,1])
    }
    if(!no_plot)abline(v=c(lam_plot[pos_all],
                           lam_plot[pos_one]),lty=4,lwd=2)
    lam_all <- lam_plot[pos_all]
    lam_one <- lam_plot[pos_one]
  }else{
    if(!no_plot)abline(v=c(lam_plot[which(y_mean==min(y_mean))],
                           lam_plot[which(y1==min(y1),arr.ind = T)[,1]]),lty=4,lwd=2)
    lam_all <- lam_plot[which(y_mean==min(y_mean))]
    lam_one <- lam_plot[which(y1==min(y1),arr.ind = T)[,1]]
  }
  if(res_perf_mdd$mod!="reg"){
    if(!no_plot)points(sort(RMSEP[,2]),y_mean,type = "l",lwd=4,lty=1,
                       col=adjustcolor(1,alpha.f = 0.2))
    if(!no_plot)points(sort(RMSEP[,2]),y_mean,type = "l",lwd=2,lty=3,
                       col=1)
  }
  if(!is.null(legend_names)){
    if(res_perf_mdd$mod!="reg"){
      if(!no_plot)legend(pos_legend,
                         legend = c(paste(legend_names,
                                          paste(" (",TAB," indiv.)",sep=""),
                                          sep=""),
                                    "Mean good classif rate"),
                         col = c(colors,1),
                         lty = c(rep(1,length(colors)),3),
                         lwd=c(rep(2,length(colors),1.5)),bty = "n")
    }else{
      if(!is.null(plot_mean)){
        legend_names <- c(legend_names, legend_0)
        col <- c(colors,"black")
        lty <- c(rep(1,length(colors)),3)
        lwd <- c(rep(2,length(colors)),3)
      }else{
        col <- colors
        lty <- rep(1,length(colors))
        lwd <- rep(2,length(colors))
      }
      if(!no_plot)legend(pos_legend,legend = legend_names,
                         col = col,lty = lty,lwd=lwd,bty = "n")
    }
  }

  if(plot_mean){
    if(!no_plot)points(sort(RMSEP[,2]),y_mean,type="l",
                       lty=3,lwd=2)
    # graphics::points(sort(RMSEP[,2]),y_mean,type="l",
    # col=grDevices::adjustcolor("black",alpha.f = 0.2),
    # lty=1,lwd=4)
  }

  if(is_L0!="L0s"){
    y_card <- card_ranges*diff(range(y1))/diff(range(card_ranges))
    y_card <- y_card - min(y_card) + min(y1)
    if(!no_plot)par(new = TRUE)
    # graphics::plot(ranges,card_ranges, type = "l", xaxt = "n", yaxt = "n",
    # ylab = "", xlab = "", col = grDevices::adjustcolor("red",0),
    # lty = 1,lwd=5)
    if(!no_plot)axis(side = 3,at=ranges,labels=card_ranges, col="red",
                     col.axis="red")
    if(!no_plot)mtext("", side = 3, line = 3, col = "red")
  }
  if(res_perf_mdd$mod=="reg"){
    if(is_L0!="L0s"){
      ranges_y_0 <- sort(apply(cc,1,max))
      if(l_lambdas>1){
        # ranges_y <- sort(ranges_y[intersect(which(ranges_y>=min(res_perf_mdd$RMSEP[,2])),
        # which(ranges_y<=max(res_perf_mdd$RMSEP[,2])))])
        # card_ranges_y <- rev(0:(length(ranges_y)-1))
        l_in_min <- ranges_y_0>=min(res_perf_mdd$RMSEP[,2])
        l_in_max <- ranges_y_0<=max(res_perf_mdd$RMSEP[,2])
        inter_in <- l_in_min&l_in_max
        l_in <- length(ranges_y_0)-which(inter_in)
        ranges_y <- ranges_y_0[length(ranges_y_0)-l_in]
        card_ranges_y <- l_in
      }else{
        ranges_y <- sort(ranges_y[which(ranges_y>=min(res_perf_mdd$RMSEP[,2]))])
        card_ranges_y <- rev(0:(length(ranges_y)-1))
      }
    }
    if(!no_occurence){
      if(!no_plot)matplot(FREQ[order(FREQ[2]),2],
                          FREQ[order(FREQ[2]),-c(1:2)]/max(FREQ[order(FREQ[2]),-c(1:2)])*100,
                          type="l",lwd=4,col=colors,lty=1,
                          xlab=expression(lambda),ylab=ylab1,
                          main=main2)
      if(is_L0!="L0s"){
        pos_y <- unique(seq(1,length(ranges_y),length.out = 15))
        pos_y[length(pos_y)] <- min(max(pos_y),length(ranges_y))
        axis(side = 3,at=ranges_y,labels=card_ranges_y, col="blue",col.axis="blue")
        mtext("", side = 3, line = 3, col = "blue")
      }
    }
  }
  out <- list(optim_para_all = lam_all, optim_para_one = lam_one)
  class(out) <- "plot.perf_mddsPLS"
  out
}
