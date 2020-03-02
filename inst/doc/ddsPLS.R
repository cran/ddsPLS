## ----include=FALSE-------------------------------------------------------
DPI=300
out.width="1000px"
out.height="1000px"
out.height_2="500px"
library(htmltools)
tagList(rmarkdown::html_dependency_font_awesome())

## ---- eval=F-------------------------------------------------------------
#  # For no named matrices
#  Xs <- list(X1,X2,X3)
#  
#  # For names matrices
#  Xs <- list(X_dna=X1,X_proteo=X2,X_RNA=X3)

## ---- fig.show='hold',message=FALSE--------------------------------------
library(ddsPLS)

## ----eval=F--------------------------------------------------------------
#  mode

## ---- eval=F-------------------------------------------------------------
#  data("penicilliumYES")
#  X <- penicilliumYES$X
#  X <- scale(X[,which(apply(X,2,sd)>0)])
#  Xs <- list(X[,1:1000],X[,-(1:1000)])
#  Xs[[1]][1:5,]=Xs[[2]][6:10,] <- NA
#  Y <- as.factor(unlist(lapply(c("Melanoconidiu","Polonicum","Venetum"),function(tt){rep(tt,12)})))
#  
#  mddsPLS_model_class <- mddsPLS(Xs = Xs,Y = Y,L0=3,R = 2,mode = "lda")

## ---- fig.show='hold',message=FALSE,eval=T-------------------------------
data("liverToxicity")
X <- scale(liverToxicity$gene)
Y <- scale(liverToxicity$clinic)
mddsPLS_model_reg <- mddsPLS(Xs = X,Y = Y,L0=10,R = 3,
                             mode = "reg")

## ----fig.height=10,fig.width=10,echo=T,dpi=DPI,out.width= out.width,out.height=out.height----
plot(mddsPLS_model_reg,vizu = "weights",variance = "Linear",
     super = T,comp = c(1,2),addY = T,mar_left = 3,mar_bottom = 3,
     pos_legend = "topright",reorder_Y = T,legend.cex = 1.5)

## ----fig.height=7,fig.width=10,dpi=DPI,out.width= out.width,out.height=out.height_2----
plot(mddsPLS_model_reg,vizu = "heatmap",comp = 1)

## ---- fig.show='hold',message=FALSE,eval=T-------------------------------
summary(mddsPLS_model_reg)

## ------------------------------------------------------------------------
mddsPLS_model_reg$var_selected[[1]]

## ----message=FALSE,eval=F------------------------------------------------
#  n_lambda <- 50
#  NCORES <- 7
#  
#  res_cv_reg_L0 <- perf_mddsPLS(Xs = X,Y = Y,
#                             R = 1,L0s=1:50,
#                             mode = "reg",NCORES = NCORES,kfolds = "loo")
#  
#  res_cv_reg <- perf_mddsPLS(Xs = X,Y = Y,
#                             R = 1,lambda_min=0.5,n_lambda=n_lambda,
#                             mode = "reg",NCORES = NCORES,kfolds = "loo")

## ----echo=F--------------------------------------------------------------
# res_cv_reg$Xs <- NULL
# res_cv_reg_L0$Xs <- NULL
# save(res_cv_reg,file="res_cv_reg_noXs.RData")
# save(res_cv_reg_L0,file="res_cv_reg_L0_noXs.RData")
load("res_cv_reg_noXs.RData")
load("res_cv_reg_L0_noXs.RData")
res_cv_reg$Xs <- list(X)
res_cv_reg_L0$Xs <- list(X)

## ---- fig.show='hold',message=FALSE,eval=T,fig.width=10, fig.height=6,dpi=DPI,out.width= out.width,out.height=out.height----
plot(res_cv_reg_L0,which_sd_plot = c(5,7),ylim=c(0,1.1),alpha.f = 0.4,
     plot_mean = T,legend_names = colnames(Y),pos_legend = "bottomright",no_occurence = T)

## ---- fig.show='hold',message=FALSE,eval=T,fig.width=10, fig.height=6,dpi=DPI,out.width= out.width,out.height=out.height----
plot(res_cv_reg,which_sd_plot = c(5,7),ylim=c(0,1.1),
     alpha.f = 0.4,plot_mean = T,legend_names = colnames(Y),
     no_occurence = T)

## ---- fig.show='hold',message=FALSE,eval=T-------------------------------
summary(res_cv_reg,plot_res_cv =F)

## ------------------------------------------------------------------------
data("liverToxicity")
X <- scale(liverToxicity$gene)
Y <- scale(liverToxicity$clinic)
p1=p2 <- 1000
p3 <- ncol(X)-p1-p2
Xs <- list(Block1=X[,1:p1],Matrix2=X[,p1+1:p2],Other_Variables=X[,p1+p2+1:p3])
Xs$Block1[1:10,] <- NA
Xs$Matrix2[5+1:10,] <- NA
Xs$Other_Variables[10+1:20,] <- NA

model_multi_vizu <- mddsPLS(Xs,Y,lambda = 0.8,R = 3)

## ----fig.height=10,fig.width=13,dpi=DPI,out.width= out.width,out.height=out.height_2----
plot(model_multi_vizu,vizu = "weights",super = T,comp=1 ,addY = T,
     mar_left = 5,mar_bottom = 3,reorder_Y = T)

## ----dpi=DPI,out.width= out.width,out.height=out.height_2----------------
plot(model_multi_vizu,vizu = "heatmap",comp = 1)

## ---- fig.show='hold',message=FALSE,eval=T-------------------------------
summary(model_multi_vizu)

## ------------------------------------------------------------------------
model_multi_vizu$var_selected[[1]]

