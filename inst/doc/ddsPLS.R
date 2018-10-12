## ---- fig.show='hold',message=FALSE--------------------------------------
library(ddsPLS)
library(doParallel)
library(RColorBrewer)
data("liver.toxicity")
X <- scale(liver.toxicity$gene)
Y <- scale(liver.toxicity$clinic)
mddsPLS_model_reg <- mddsPLS(Xs = X,Y = Y,lambda=0.9,R = 1,
                             mode = "reg",verbose = TRUE)

## ----fig.width=7, fig.height=10,message=FALSE----------------------------
res_cv_reg <- perf_mddsPLS(Xs = X,Y = Y,
                           R = 1,lambda_min=0.7,n_lambda=5,
                           mode = "reg",NCORES = 1,kfolds = "loo")
plot(res_cv_reg,legend_names=colnames(Y))

## ---- fig.show='hold',message=FALSE--------------------------------------
data("penicilliumYES")
X <- penicilliumYES$X
X <- scale(X[,which(apply(X,2,sd)>0)])
classes <- c("Melanoconidium","Polonicum","Venetum")
Y <- as.factor(unlist(lapply(classes,
                             function(tt){rep(tt,12)})))
mddsPLS_model_class <- mddsPLS(Xs = X,Y = Y,lambda = 0.958,R = 2,
                               mode = "clas",verbose = TRUE)
## Plot the two first axes 
plot(mddsPLS_model_class$mod$t,col=Y,pch=as.numeric(Y)+15,cex=2,
     xlab="1st variable, 2 var. selected",
     ylab="2nd variable, 1 var. selected")
legend(-2,0,legend=classes,col=1:3,pch=15+(1:3),box.lty=0,y.intersp=2)

## ----fig.width=7, fig.height=6,message=FALSE-----------------------------
res_cv_class <- perf_mddsPLS(X,Y,R = 2,lambda_min=0.94,n_lambda=5,
                             mode = "clas",NCORES = 1,
                             fold_fixed = rep(1:12,3))
plot(res_cv_class,legend_names = levels(Y))

