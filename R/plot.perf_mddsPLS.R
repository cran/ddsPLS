#' Function to plot cross-validation performance results.
#'
#' That function must be applied to a perf_mddsPLS object. Extra parameters are
#'  avalaible to control the plot quality.
#'
#' @param x The perf_mddsPLS object.
#' @param plot_mean logical. Whether or not to plot the mean curve.
#' @param pos_legend character. One of "bottomleft", "topright",....
#' @param legend_names vector of characters. Each element is the name of one of the q response variables.
#' @param ... Other plotting parameters to affect the plot.
#'
#' @return The plot visualisation
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
#' #res_cv_class <- perf_mddsPLS(X,Y,lambda_min=0.85,n_lambda=2,R = 2,
#' #mode = "clas",NCORES = 1,fold_fixed = rep(1:12,3))
#' #plot(res_cv_class)
#'
#' # Regression example :
#' data("liver.toxicity")
#' X <- scale(liver.toxicity$gene)
#' Y <- scale(liver.toxicity$clinic)
#' #res_cv_reg <- perf_mddsPLS(Xs = X,Y = Y,lambda_min=0.8,n_lambda=2,R = 1,
#' # mode = "reg")
#' #plot(res_cv_reg)
plot.perf_mddsPLS <- function(x,plot_mean=FALSE,legend_names=NULL,
                              pos_legend="bottomleft",...){
  res_perf_mdd <- x
  names(res_perf_mdd)[1] <- "RMSEP"
  X_all <- scale(do.call(cbind,res_perf_mdd$Xs))
  if(res_perf_mdd$mode=="reg"){
    cc <- matrix(NA,nrow = ncol(res_perf_mdd$Y),ncol = ncol(X_all))
    for(j in 1:ncol(X_all)){
      cc[,j] <- abs(stats::cor(res_perf_mdd$Y,X_all[,j],use = "pairwise.complete.obs"))
    }
    # cc <- abs(crossprod(scale(res_perf_mdd$Y),X_all)/(nrow(res_perf_mdd$Y)-1))
    col_na <- which(is.na(colSums(cc)))
    if(length(col_na)>0){
      cc <- cc[,-col_na,drop=F]
    }
  }
  else{
    Y_df <- data.frame(res_perf_mdd$Y)
    colnames(Y_df) <- "Y"
    Y <- scale(stats::model.matrix( ~ Y - 1, data=Y_df))
    cc <- abs(crossprod(Y,X_all)/(nrow(Y)-1))
  }
  ranges <- apply(cc,2,max)
  l_lambdas <- length(unique(res_perf_mdd$RMSEP[,2]))
  if(ncol(cc)>1){
    if(l_lambdas>1){
      ranges <- sort(ranges[intersect(which(ranges>=min(res_perf_mdd$RMSEP[,2])),
                                      which(ranges<=max(res_perf_mdd$RMSEP[,2])))])
      card_ranges <- rev(0:(length(ranges)-1))
    }else{
      ranges <- sort(ranges[which(ranges>=min(res_perf_mdd$RMSEP[,2]))])
      card_ranges <- rev(0:(length(ranges)-1))
    }
  }else{
    card_ranges <- 1
  }
  ERRORS <- res_perf_mdd
  FREQ <- ERRORS$FREQ
  RMSEP <- ERRORS$RMSEP
  q <- ncol(ERRORS$RMSEP)-2
  if(q<3){
    colors <- 1:q
  }
  else if(q>8){
    colors <- RColorBrewer::brewer.pal(8, "Dark2")
    pal <- grDevices::colorRampPalette(colors)
    colors <- pal(q)
  }
  else{
    colors <- RColorBrewer::brewer.pal(q, "Dark2")
  }
  if(res_perf_mdd$mod=="reg"){
    ylab1<-"MSEP"
    ylab2<-"Occurences per variable (%)"
    ylim1 <- range(abs(RMSEP[,3:ncol(RMSEP)]))^2
    y1 <- RMSEP[order(RMSEP[,2,drop=FALSE]),3:ncol(RMSEP),drop=FALSE]^2
    y_mean <- rowMeans(RMSEP[order(RMSEP[,2,drop=FALSE]),3:ncol(RMSEP),drop=FALSE]^2)
    main1 <- "MSEP versus regularization coefficient\n mdd-sPLS"
    main2 <- "Occurences per variable versus regularization coefficient\n mdd-sPLS"
    graphics::par(mar=c(5,5,7,5),mfrow=c(2,1))
  }
  else{
    ylab1<-"#Good Classif Rate"
    ylab2<-"Occurences per class"
    ylim1<- c(0,1)
    y1 <- RMSEP[order(RMSEP[,2,drop=FALSE]),3:ncol(RMSEP),drop=FALSE]
    TAB <- table(res_perf_mdd$Y)
    for(r in 1:nlevels(res_perf_mdd$Y)){
      y1[,r] <- 1-y1[,r]/TAB[r]
    }
    y_mean <- 1-rowSums(RMSEP[order(RMSEP[,2,drop=FALSE]),3:ncol(RMSEP),drop=FALSE])/sum(TAB)
    main1 <- "Good classification rate versus regularization coefficient\n mdd-sPLS"
    main2 <- "Occurences per class versus regularization coefficient\n mdd-sPLS"
    graphics::par(mar=c(5,5,6,5),mfrow=c(1,1))
  }
  graphics::matplot(sort(RMSEP[,2]),y1,type="l",lwd=4,lty=1,
                    ylim=ylim1,col=colors,
                    xlab=expression(lambda),
                    ylab=ylab1,
                    main=main1)
  if(res_perf_mdd$mod!="reg"){
    graphics::points(sort(RMSEP[,2]),y_mean,type = "l",lwd=4,lty=1,
                     col=grDevices::adjustcolor(1,alpha.f = 0.2))
    graphics::points(sort(RMSEP[,2]),y_mean,type = "l",lwd=2,lty=3,
                     col=1)
  }
  if(!is.null(legend_names)){
    if(res_perf_mdd$mod!="reg"){
      graphics::legend(pos_legend,
                       legend = c(paste(legend_names,paste(" (",TAB," indiv.)",sep=""),sep=""),
                                  "Mean good classif rate"),
                       col = c(colors,1),lty = c(rep(1,length(colors)),3),
                       lwd=c(rep(2,length(colors),1.5)))
    }else{
      graphics::legend(pos_legend,legend = legend_names,col = colors,lty = 1,lwd=2)
    }
  }
  if(plot_mean){
    graphics::points(sort(RMSEP[,2]),y_mean,type="l",lty=3)
    graphics::points(sort(RMSEP[,2]),y_mean,type="l",col=grDevices::adjustcolor("black",alpha.f = 0.2),lty=1,lwd=4)
  }
  y_card <- card_ranges*diff(range(y1))/diff(range(card_ranges))
  y_card <- y_card - min(y_card) + min(y1)
  graphics::par(new = TRUE)
  graphics::plot(ranges,card_ranges, type = "l", xaxt = "n", yaxt = "n",
                 ylab = "", xlab = "", col = grDevices::adjustcolor("red",0), lty = 1,lwd=5)
  graphics::axis(side = 3,at=ranges,labels=card_ranges, col="red",col.axis="red")
  graphics::mtext("", side = 3, line = 3, col = "red")
  if(res_perf_mdd$mod=="reg"){
    ranges_y <- apply(cc,1,max)
    if(l_lambdas>1){
      ranges_y <- sort(ranges_y[intersect(which(ranges_y>=min(res_perf_mdd$RMSEP[,2])),
                                          which(ranges_y<=max(res_perf_mdd$RMSEP[,2])))])
      card_ranges_y <- rev(0:(length(ranges_y)-1))
    }else{
      ranges_y <- sort(ranges_y[which(ranges_y>=min(res_perf_mdd$RMSEP[,2]))])
      card_ranges_y <- rev(0:(length(ranges_y)-1))
    }

    graphics::matplot(FREQ[order(FREQ[2]),2],
                      FREQ[order(FREQ[2]),-c(1:2)]/max(FREQ[order(FREQ[2]),-c(1:2)])*100,type="l",lwd=4,col=colors,lty=1,
                      xlab=expression(lambda),
                      ylab=ylab2,
                      main=main2)
    pos_y <- unique(seq(1,length(ranges_y),length.out = 15))
    pos_y[length(pos_y)] <- min(max(pos_y),length(ranges_y))
    graphics::axis(side = 3,at=ranges_y,labels=card_ranges_y, col="blue",col.axis="blue")
    graphics::mtext("", side = 3, line = 3, col = "blue")
  }
}
