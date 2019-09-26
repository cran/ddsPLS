#' Function to plot \strong{mddsPLS}
#'
#' That function must be applied to a \strong{mddsPLS} object. Extra parameters are
#'  avalaible to control the plot quality. Keep in mind that if a lot of covariates
#'  are selected, their names might not all fit the plot window, only the names of
#'  the most important covariates are present. To provide the names of all the
#'  covariates, the user can modify concerned parameters of the \strong{barplot}
#'  function (for example the \emph{cex.names} parameter).
#'
#' @param x The perf_mddsPLS object.
#' @param vizu character. One of \strong{weights}, \strong{coeffs}, \strong{heatmap}, \strong{correlogram}. \strong{coeffs} does not work in the case of classification (\strong{lda} or \strong{logit}). If \strong{heatmap} is selected, light colors correspond to low expressions.
#' @param super logical. If \strong{TRUE} barplots are filled with **Super-Weights** in the case of \strong{vizu=weights} of with général **X** and **Y** components else.
#' @param addY logical. Whether or not to plot **Block Y**. Initialized to \strong{FALSE}.
#' @param block vector of intergers indicating which components must be plotted. If equals \strong{NULL} then all the components are plotted. Initialized to \strong{NULL}.
#' @param comp vector of intergers indicating which blocks must be plotted. If equals \strong{NULL} then all the blocks are plotted. Initialized to \strong{NULL}.
#' @param variance character. One of \strong{Linear}, \strong{RV}. Explains the type of variance shown in the graphics.
#' @param mar_left positive float. Extra lines to add to the left margins, where the variable names are written.
#' @param mar_bottom positive float. Extra lines to add to the bottom margins. Useful when \strong{addY}=TRUE.
#' @param margins_heatmap vector of 2 positive floats. The \strong{margins} argument of the \strong{heatmap} function. Margins to the bottom and to the right of the heatmap, if plotted. Useful if samples and covariates have particularly long names. Default to c(5,5).
#' @param pos_legend Initialized to "topright". If equals NULL, then no legend is given.
#' @param legend.cex positive float. character expansion factor relative to current par("cex") for \strong{legend} function.
#' @param legend_names vector of character. Indicates the names of the blocks. Initialized to NULL and in this case just gets positions in the Xs list.
#' @param block_Y_name character. Initialized to "Block Y".
#' @param alpha.Y_sel positive float. factor modifying the opacity alpha; typically in \code{[0,1]} from \code{adjustcolor} function.
#' @param reorder_Y logical. In case \strong{addY}=TRUE. Order the \strong{Y} variances according to proportion of varaince explained on the first component.
#' @param values_corr logical. Wether of noth to write the correlation calues in the correlogram. Initialized to FALSE
#' @param ... Other plotting parameters to affect the plot.
#'
#' @return The plot visualisation
#'
#' @importFrom graphics abline arrows barplot legend par text
#' @importFrom stats heatmap model.matrix
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette adjustcolor
#' @importFrom corrplot corrplot
#'
#' @seealso  \code{\link{mddsPLS}}, \code{\link{summary.mddsPLS}}
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
#' # x <- mddsPLS(Xs = X,Y = Y,R = 3, mode = "lda",L0=20)
#' # plot(x)
#'
#' # Regression example :
#' data("liverToxicity")
#' X <- scale(liverToxicity$gene)
#' Y <- scale(liverToxicity$clinic)
#' # mod_reg <- mddsPLS(Xs = X,Y = Y,L0=10,R = 2)
#' # plot(mod_reg,addY = T,mar_left = 3)
#' # plot(mod_reg,addY = T,mar_left = 3,super = T)
plot.mddsPLS <- function(x,vizu="weights",super=FALSE,addY=FALSE,
                         block=NULL,comp=NULL,variance="Linear",
                         mar_left=2,mar_bottom=2,margins_heatmap=c(5,5),
                         pos_legend="topright",legend_names=NULL,legend.cex=1,
                         values_corr=F,block_Y_name="Y",alpha.Y_sel=0.4,
                         reorder_Y=F,
                         ...){
  ## Reset personnal plot par() settings
  opar <- par(no.readonly =TRUE)
  on.exit(par(opar))
  ## -----------------------------------
  NZV <- x$NZV
  ## Functions ---------------
  ##### HEATMAP FUNCTION .....
  plot_heatmap <- function(x,comp=NULL,out=F,position="topright",marginsRight=5){
    Xs <- x$Xs
    K <- length(Xs)
    if(!is.null(names(Xs))){
      for(k in 1:K){
        if(nchar(names(Xs)[k])==0){
          names(Xs)[k] <- paste("Block",k)
        }
      }
    }else{
      for(k in 1:K){
        names(Xs)[k] <- paste("Block",k)
      }
    }
    Y_1 <- x$Y_0
    R <- length(x$mod$T_super)
    comp_in<-comp
    if(is.null(comp_in)){
      comp_in <- 1
    }
    if(x$mode!="reg"){
      Y_1 <- as.matrix(model.matrix( ~ Y - 1, data=data.frame(x$Y_0,ncol=1)))
      colnames(Y_1) <- levels(as.factor(x$Y_0))
      q <- ncol(Y_1)
    }else if(!is.matrix(Y_1)){
      Y_1 <- matrix(Y_1,ncol=1)
    }
    n <- nrow(Y_1)
    q <- ncol(Y_1)
    which_sel <- lapply(x$mod$u_t_super,function(u){which(abs(u[,comp_in])>NZV)})
    p_sel <- sum(unlist(lapply(which_sel,length)))
    coco_imputed <- data.frame(matrix(NA,n,p_sel+q))
    coeffs <- matrix(rep(0,p_sel),nrow = 1)
    colnames(coeffs) <- rep("OOO",p_sel)
    count <- 0
    my_group <- rep(block_Y_name,p_sel+q)
    for(k in 1:K){
      Xs_k <- Xs[[k]]
      if(!is.matrix(Xs_k)){
        Xs_k <- matrix(Xs_k,nrow=1)
      }
      ls <- which_sel[[k]]
      if(length(ls)!=0){
        my_group[count + 1:length(ls)] <- names(Xs)[k]
        coco_imputed[,count + 1:length(ls)] <- Xs[[k]][,ls,drop=F]
        coeffs[1,count + 1:length(ls)] <- x$mod$u_t_super[[k]][ls,1]
        colNames <- colnames(Xs[[k]])[ls]
        if(!is.null(colNames)){
          names(coco_imputed)[count + 1:length(ls)] <- colNames
        }else{
          names(coco_imputed)[count + 1:length(ls)] <- paste(names(Xs)[k],", Var. ",ls,sep="")
        }
      }
      count <- count + length(ls)
    }
    coco_imputed[,count + 1:q] <- Y_1
    names_Y <- colnames(Y_1)
    if(is.null(colnames(Y_1))){
      names_Y <- paste("Y, Var.",1:q)
    }
    colnames(coco_imputed)[count +1:q] <- names_Y
    my_group_factor <- factor(my_group)
    if(K+1<3){
      colors <- c("black","red","blue")[1:(K+1)]
    }else if(K+1>8){
      colors <- brewer.pal(8, "Dark2")
      pal <- colorRampPalette(colors)
      colors <- pal(K+1)
    }else{
      colors <- brewer.pal(K+1, "Dark2")
    }
    my_group_factor_plot <- factor(as.character(my_group_factor),levels=c(names(Xs),block_Y_name))
    my_col <- matrix(colors[my_group_factor_plot],nrow=1)
    rownames(my_col) <- "Block"
    main <- paste("Heatmap for component",comp_in)
    if(variance=="Linear"){
      var_here <- signif(x$Variances$Linear$VAR_SUPER_COMPS_ALL_Y[comp_in],2)*100
      main <- paste(main," (",var_here,"% var. expl.)",sep="")
    }else{
      var_here <- signif(x$Variances$RV$VAR_SUPER_COMPS_ALL_Y[comp_in],2)*100
      main <- paste(main," (RV=",var_here/100,")",sep="")
    }
    level_present <- as.character(unique(my_group_factor))
    if(!out){
      heatmap(t(as.matrix(coco_imputed)),scale="row",labCol = "",
              xlab = "Individuals",RowSideColors=my_col,
              main=main,margins = margins_heatmap)
      legend(position,legend=level_present,
             fill=colors[sort(unique(as.numeric(my_group_factor_plot)))],
             border=FALSE, bty="n",cex=legend.cex)
    }
    output <- coco_imputed
    if(out){
      return(output)
    }
  }
  ##### CORRPLOT FUNCTION ....
  plot_corcor <- function(x,comp=NULL,values=F){
    Xs <- x$Xs
    K <- length(Xs)
    if(!is.null(names(Xs))){
      for(k in 1:K){
        if(nchar(names(Xs)[k])==0){
          names(Xs)[k] <- paste("Block",k)
        }
      }
    }else{
      for(k in 1:K){
        names(Xs)[k] <- paste("Block",k)
      }
    }
    Y_1 <- x$Y_0
    R <- length(x$mod$T_super)
    comp_in<-comp
    if(is.null(comp_in)){
      comp_in <- 1
    }
    if(x$mode!="reg"){
      Y_1 <- as.matrix(model.matrix( ~ Y - 1, data=data.frame(x$Y_0,ncol=1)))
      colnames(Y_1) <- levels(as.factor(x$Y_0))
      q <- ncol(Y_1)
    }else if(!is.matrix(Y_1)){
      Y_1 <- matrix(Y_1,ncol=1)
    }
    n <- nrow(Y_1)
    q <- ncol(Y_1)
    which_sel <- lapply(x$mod$u_t_super,function(u){which(abs(u[,comp_in])>NZV)})
    p_sel <- sum(unlist(lapply(which_sel,length)))
    coco_imputed <- data.frame(matrix(NA,n,p_sel+q))
    coeffs <- matrix(rep(0,p_sel),nrow = 1)
    colnames(coeffs) <- rep("OOO",p_sel)
    count <- 0
    my_group <- rep(block_Y_name,p_sel+q)
    for(k in 1:K){
      Xs_k <- Xs[[k]]
      if(!is.matrix(Xs_k)){
        Xs_k <- matrix(Xs_k,nrow=1)
      }
      ls <- which_sel[[k]]
      if(length(ls)!=0){
        my_group[count + 1:length(ls)] <- names(Xs)[k]
        coco_imputed[,count + 1:length(ls)] <- Xs[[k]][,ls,drop=F]
        coeffs[1,count + 1:length(ls)] <- x$mod$u_t_super[[k]][ls,1]
        colNames <- colnames(Xs[[k]])[ls]
        if(!is.null(colNames)){
          names(coco_imputed)[count + 1:length(ls)] <- colNames
        }else{
          names(coco_imputed)[count + 1:length(ls)] <- paste(names(Xs)[k],", Var. ",ls,sep="")
        }
      }
      count <- count + length(ls)
    }
    coco_imputed[,count + 1:q] <- Y_1
    names_Y <- colnames(Y_1)
    if(is.null(colnames(Y_1))){
      names_Y <- paste("Y, Var.",1:q)
    }
    colnames(coco_imputed)[count +1:q] <- names_Y
    my_group_factor <- factor(my_group)
    if(K+1<3){
      colors <- c("black","red","blue")[1:(K+1)]
    }else if(K+1>8){
      colors <- brewer.pal(8, "Dark2")
      pal <- colorRampPalette(colors)
      colors <- pal(K+1)
    }else{
      colors <- brewer.pal(K+1, "Dark2")
    }
    my_group_factor_plot <- factor(as.character(my_group_factor),levels=c(names(Xs),block_Y_name))
    my_col <- matrix(colors[my_group_factor_plot],nrow=1)
    rownames(my_col) <- "Block"
    level_present <- as.character(unique(my_group_factor))
    if(values){
      corrplot(cor(plot_heatmap(x,comp=comp,out=T)),method="number",tl.col = my_col)
    }else{
      corrplot(cor(plot_heatmap(x,comp=comp,out=T)),tl.col = my_col)
    }
    legend("topleft",legend=level_present,
           fill=colors[sort(unique(as.numeric(my_group_factor_plot)))],
           border=FALSE, bty="n",cex=legend.cex)
  }
  ## END FUNCTIONS -----------
  R <- x$mod$R
  if(is.null(comp)){
    comp_in <- 1:R
  }else{
    comp_in <- comp
  }
  R_in <- length(comp_in)
  K <- length(x$Xs)
  block_in <- block
  if(is.null(block_in) | super){
    block_in <- 1:K
  }
  if (any(!comp_in %in% 1:R)) {
    stop("One of the asked components does not exist",
         call. = FALSE)
  }else if(any(!block_in %in% 1:K)){
    stop("One of the asked block does not exist",
         call. = FALSE)
  }
  isReg <- x$mode=="reg"
  if(!isReg){
    if(vizu %in% c("coeffs")){
      stop("Cannot be performed for classification. \n
           Consider changing 'vizu' to something different from 'coeffs'.",
           call. = FALSE)
    }
    if(addY){
      stop("Do not print Y variance explanation for classification. \n
           Consider putting 'addY' to FALSE for example.",
           call. = FALSE)
    }
  }
  Y_in <- x$Y_0
  if(!isReg){
    names_Y <- levels(Y_in)
    q <- nlevels(Y_in)
  }else{
    if(!(is.matrix(Y_in)|is.data.frame(Y_in))){
      Y_in <- as.matrix(Y_in)
    }
    if(is.data.frame(Y_in)){
      Y_in <- as.matrix(Y_in)
    }
    q <- ncol(Y_in)
    names_Y <- colnames(Y_in)
    if(is.null(names_Y)){
      names_Y <- 1:q
    }
  }
  legends_names_y <- colnames(Y_in)
  legend_names_in <- names(x$Xs)
  if(is.null(legend_names) & is.null(legend_names_in)){
    legend_names_in <- paste("Block",block_in,sep=" ")
  }else{
    if(K!=1){
      for(k in 1:K){
        if(nchar(names(x$Xs)[k])==0){
          legend_names_in[k] <- paste("Block",k)
        }
      }
    }else{
      legend_names_in <- legend_names
    }
  }
  l_bl <- K+1
  if(l_bl<3){
    colors <- 1:l_bl
  }else if(l_bl>8){
    colors <- brewer.pal(8, "Dark2")
    pal <- colorRampPalette(colors)
    colors <- pal(l_bl)
  }else{
    colors <- brewer.pal(l_bl, "Dark2")
  }
  if(vizu=="weights"){
    viz <- x$mod$u
    viz_y <- x$mod$V_super
    if(super){
      viz <- x$mod$u_t_super
      block_in <- 1:K
    }
    toplot <- list()
    ind_1 <- 1:length(block_in)
    ind_2 <- 1:R_in
    if(!super){
      if(addY){
        par(mfrow=c(length(block_in),R_in+1),
            mar=c(5+mar_bottom,4+mar_left,4,2)+0.1)
      }else{
        par(mfrow=c(length(block_in),R_in),
            mar=c(5,4+mar_left,4,2)+0.1)
      }
    }else{
      if(addY){
        par(mfrow=c(R_in,1+1),
            mar=c(5+mar_bottom,4+mar_left,4,2)+0.1)
      }else{
        par(mfrow=c(R_in,1),
            mar=c(5,4+mar_left,4,2)+0.1)
      }
    }
    for(i_k in ind_1){
      for(i_r in ind_2){
        r <- comp_in[i_r]
        k <- block_in[i_k]
        if(i_r==1){
          toplot[[k]] <- list()
        }
        viz_k <- viz[[k]]
        viz_k_r <- viz_k[,r]
        pos_no_nul <- which(abs(viz_k_r)>NZV)
        main <- paste(legend_names_in[i_k],", component ",r,sep="")
        if(isReg){
          if(variance=="Linear"){
            var_here <- signif(x$Variances$Linear$VAR_COMPS[k,r],2)*100
            main <- paste(main," (",var_here,"% var. expl.)",sep="")
          }else{
            var_here <- signif(x$Variances$RV$VAR_COMPS[k,r],2)*100
            main <- paste(main," (RV=",var_here/100,")",sep="")
          }
        }
        if(length(pos_no_nul)>0){
          toplot[[k]][[r]] <- viz_k[pos_no_nul,r]
          names(toplot[[k]][[r]]) <- colnames(x$Xs[[k]])[pos_no_nul]
          if(is.null(names(toplot[[k]][[r]]))){
            names(toplot[[k]][[r]]) <- pos_no_nul
          }
          toplot[[k]][[r]] <- toplot[[k]][[r]][order(abs(toplot[[k]][[r]]),
                                                     decreasing = T)]

          if(!super){
            barplot(toplot[[k]][[r]],horiz = T,las=2,
                    col=colors[k],xlim = c(-1,1),
                    main=main,xlab="Coefficient",...)
            abline(v=c(0.5,-0.5,-1,1),lty=c(2,2,1,1),col=adjustcolor("black",alpha.f = 0.2))
          }
        }else{
          toplot[[k]][[r]] <- 0
          if(!super){
            plot(1, type="n", axes=F, xlab="", ylab="",main=main)
          }
        }
      }
      if(addY & !super){
        y_como <- viz_y[,r]
        pos_no_nul <- which(abs(y_como)>NZV)
        if(length(pos_no_nul)>0){
          y_como <- y_como[pos_no_nul]
          names(y_como) <- names_Y[pos_no_nul]
        }
        toplot_y <- y_como[order(abs(y_como),decreasing = T)]
        barplot(toplot_y,horiz = T,las=2,col=colors[K+1],xlim = c(-1,1),
                main=paste(block_Y_name," component ",r,sep=""),xlab="Coefficient",...)
        abline(v=c(0.5,-0.5,-1,1),lty=c(2,2,1,1),col=adjustcolor("black",alpha.f = 0.2))
        legeds <- c(legend_names_in[block_in],block_Y_name)
        colOut <- colors[c(block_in,K+1)]
      }else{
        legeds <- legend_names_in
        colOut <- colors[block_in]
      }
    }
    if(super){
      for(i_r in 1:R_in){
        r <- comp_in[i_r]
        plotR <- NULL
        cols <- NULL
        for(k in block_in){
          if(toplot[[k]][[r]][1]!=0){
            plotR <- c(plotR,toplot[[k]][[r]])
            cols <- c(cols,rep(colors[k],length(toplot[[k]][[r]])))
          }
        }
        main <- paste("Block Xs, Super Component ",r,sep="")
        if(isReg){
          if(variance=="Linear"){
            var_here <- signif(x$Variances$Linear$VAR_SUPER_COMPS_ALL_Y[r],2)*100
            main <- paste(main," (",var_here,"% var. expl. total Y)",sep="")
          }else{
            var_here <- signif(x$Variances$RV$VAR_SUPER_COMPS_ALL_Y[r],2)*100
            main <- paste(main," (RV=",var_here/100,")",sep="")
          }
        }
        if(is.null(plotR)){
          plot(1, type="n", axes=F, xlab="", ylab="",main=main)
        }else{
          oo <- order(abs(plotR),decreasing = T)
          barplot(plotR[oo],horiz = T,las=2,col=cols[oo],
                  main=main,xlim=c(-1,1),xlab="Coefficient",...)
          abline(v=c(0.5,-0.5,-1,1),lty=c(2,2,1,1),col=adjustcolor("black",alpha.f = 0.2))
        }
        if(addY){
          y_como <- viz_y[,r]
          pos_no_nul <- which(abs(y_como)>NZV)
          if(length(pos_no_nul)>0){
            y_como <- y_como[pos_no_nul]
            names(y_como) <- names_Y[pos_no_nul]
          }
          if(variance=="Linear"){
            if(!is.null(dim(x$Variances$Linear$VAR_SUPER_COMPS))){
              toplot_y <- x$Variances$Linear$VAR_SUPER_COMPS[,r]*100
            }else{
              toplot_y <- x$Variances$Linear$VAR_SUPER_COMPS[r]*100
            }
            xlab <- "Variance Explained (%)"
          }else{
            if(!is.null(dim(x$Variances$RV$VAR_SUPER_COMPS))){
              toplot_y <- x$Variances$RV$VAR_SUPER_COMPS[,r]*100
            }else{
              toplot_y <- x$Variances$RV$VAR_SUPER_COMPS[r]*100
            }
            xlab <- "Variance in Common"
          }
          y_selected <- which(abs(x$mod$V_super[,r])>NZV)
          fonts <- rep(1,length(toplot_y))
          fonts[y_selected] <- 2
          if(reorder_Y){
            new_order <- order(x$Variances$RV$VAR_SUPER_COMPS[,1],decreasing=T)
            toplot_y <- toplot_y[new_order]
            y_selected <- which(new_order %in% y_selected)
            fonts <- fonts[new_order]
          }
          if(length(y_selected)>0){
            cols_y <- rep(colors[K+1],q)
            cols_y[-y_selected] <- adjustcolor(colors[K+1],alpha.f = alpha.Y_sel)
            legeds <- c(legend_names_in,paste(block_Y_name,c("selected","not selected")))
            colOut <- c(colors[c(block_in,K+1)],adjustcolor(colors[K+1],alpha.f = alpha.Y_sel))

          }else{
            cols_y <- rep(adjustcolor(colors[K+1],alpha.f = alpha.Y_sel),q)
            legeds <- c(legend_names_in,block_Y_name)
            colOut <- colors[c(block_in,K+1)]
          }
          xx<-barplot(toplot_y,horiz = F,las=2,col=cols_y,ylim = c(0,119),
                      main=paste("Block Y, component ",r,sep=""),
                      ylab=xlab,...)
          abline(h=c(0.25,0.5,0.75,1)*100,lty=c(3,3,2,1),lwd=c(0.5,1,1,1),
                 col=adjustcolor("black",alpha.f = 0.2))
          text(xx,toplot_y,labels=round(toplot_y,0),pos=3,font=fonts)
        }else{
          legeds <- legend_names_in
          colOut <- colors[block_in]
        }
      }
    }
    if(!is.null(pos_legend)){
    legend(pos_legend,legend = legeds,
           fill = colOut,bty="n",cex=legend.cex)
    }

  }else if(vizu=="heatmap"){
    plot_heatmap(x,comp,position = pos_legend,marginsRight=margins_heatmap)
  }else if(vizu=="correlogram"){
    plot_corcor(x,comp,values=values_corr)
  }else if(vizu=="coeffs"){
    y_selected <- which(colSums(abs(do.call(rbind,x$mod$B)))>NZV)
    t_selected <- which(unlist(lapply(x$mod$B,norm))>NZV)
    par(mfrow=c(length(y_selected),length(t_selected)),mar = c(5,4 + mar_left, 4, 2) + 0.1)
    for(j in y_selected){
      for(i_t in 1:length(t_selected)){
        t <- t_selected[i_t]
        oo <- x$mod$B[[t]][,j]
        pos_ok <- which(abs(oo)>NZV)
        if(length(pos_ok)>0){
          ok <- x$mod$B[[t]][pos_ok,j,drop=F]
          rownames(ok) <- colnames(x$Xs[[t]])[pos_ok]
          barplot(t(ok),las = 2,horiz = T,col=colors[t],
                  main=paste("Block X : ",legend_names_in[t],", variable Y : ",
                             names_Y[j],sep=""),...)
        }else{
          plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
        }
      }
    }
    if(!is.null(pos_legend)){
    legend(pos_legend,legend = legend_names_in[t_selected],
           fill = colors[t_selected],bty="n",cex=legend.cex)
    }
  }
}
