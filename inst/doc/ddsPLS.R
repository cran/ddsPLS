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
                           R = 1,lambda_min=0.4,n_lambda=2,
                           mode = "reg",NCORES = 1,
                           kfolds = 2)
plot(res_cv_reg,legend_names=colnames(Y))

## ---- fig.show='hold',message=FALSE--------------------------------------
data("penicilliumYES")
X <- penicilliumYES$X
X <- scale(X[,which(apply(X,2,sd)>0)])
Y <- as.factor(unlist(lapply(c("Melanoconidiu","Polonicum","Venetum"),
                             function(tt){rep(tt,12)})))
mddsPLS_model_class <- mddsPLS(Xs = X,Y = Y,lambda = 0.958,R = 2,
                               mode = "clas",verbose = TRUE)


## ----fig.width=7, fig.height=6,message=FALSE-----------------------------
res_cv_class <- perf_mddsPLS(X,Y,R = 2,lambda_min=0.94,n_lambda=2,
                             mode = "clas",NCORES = 1,
                             fold_fixed = rep(1:12,3))
plot(res_cv_class,legend_names = levels(Y))

