## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

getData <- function(n=100,alpha=0.4,beta_0=0.2,sigma=0.3,
                    p1=50,p2=25,p3=25,p=1000){
  R1 = R2 = R3 <- 1
  d <- R1+R2+R3
  A0 <- matrix(c(
    c(rep(1/sqrt(R1),p1),rep(sqrt(alpha),p2),rep(0,p3),rep(0,p-p1-p2-p3)),
    c(rep(0,p1),rep(sqrt(1-alpha),p2),rep(0,p3),rep(0,p-p1-p2-p3)),
    c(rep(0,p1),rep(0,p2),rep(1,p3),rep(0,p-p1-p2-p3))
  ),nrow = d,byrow = TRUE)
  A <- eps*A0
  D0 <- matrix(c(1,0,0,
                 sqrt(beta_0),sqrt(1-beta_0),0,
                 0,0,0),nrow = d,byrow = FALSE)
  D <- eps*D0
  q <- ncol(D)
  L_total <- q+p
  psi <- MASS::mvrnorm(n,mu = rep(0,d+L_total),Sigma = diag(d+L_total))
  phi <- psi[,1:d,drop=F]
  errorX <- matrix(rep(sqrt(1-apply(A^2,2,sum)),n),n,byrow = TRUE)
  errorY <- matrix(rep(sqrt(1-apply(D^2,2,sum)),n),n,byrow = TRUE)
  X <- phi%*%A + errorX*psi[,d+1:p,drop=F]
  Y <- phi%*%D + errorY*psi[,d+p+1:q,drop=F]
  list(X=X,Y=Y)
}

R2mean_diff_Q2mean <- cbind(
  c(0.0390018049710666,0.0341751037010943,0.0300566345943933,0.0270372843280886,0.0251730242086019,0.0242227185400858,0.0234891342998435,0.0193133989423563,0.0169097166961483,0.0171263425383856,0.0170385680914704,0.0163279588077112,0.0152494718863702,0.0144263961042266,0.0139846253273588,0.0138793727811233,0.0138875727882883,0.0139142088056178,0.013947769332828,0.0151373262113938,0.0151964897734,0.0152786585690076,0.0184760931873834,0.0245523530157278,0.0329589020372058,0.0414793872444706,0.0487733189333913,0.0513953941897703,0.0513953941897703,0.0513953941897703)
  ,
  c(0.0511147118390687,0.0366502939988904,0.0244192885247774,0.0157596067948765,0.0103693740213229,0.00774690072471196,0.00654474100829194,0.00602821715754076,0.00456584615102629,0.000647908840480049,-0.0012117730626553,-0.00324779272702536,-0.00339597202321107,-0.00349873504246012,-0.00347050350187528,-0.00343289518836143,-0.0033818961424037,-0.00330807936244493,-0.00318820602484982,-0.00296287672757123,-0.00264380624349336,-0.00161555612407471,-0.00511044138700312,-0.00517941187733331,-0.0027506954114751,-0.000122694727787143,0,0,0,0)
)

Q2mean <- cbind(
  c(
    0.377432184277843,0.380140639849587,0.381440573733071,0.381211078689433,0.379661476366436,0.377071669682569,0.373718997129306,0.371737203943991,0.366583285363937,0.357192628930804,0.346720293574788,0.337697081091153,0.332215307215219,0.330014817204357,0.329463268875761,0.329394099789224,0.329391348298594,0.32939602992879,0.329401599141234,0.32892978442128,0.328934064910596,0.328947551675465,0.327555460720922,0.325057907988492,0.321266212859659,0.318924669092771,0.315393736421526,0.295485665498064,0.295485665498064,0.295485665498064
  ),
  c(
    0.550769148647778,0.564051120615041,0.574049282463525,0.579753590218078,0.582232070714579,0.582828526524413,0.582822316455892,0.582729371555363,0.583602822903569,0.586270800202963,0.587558962483909,0.588999196398028,0.589131477829779,0.589208302519195,0.589207789772327,0.58920691960863,0.589205165067086,0.589201245713807,0.589190361787023,0.589150058389219,0.581193572554757,0.540787461471241,0.287353467160313,-0.174809002088594,-0.730811514935191,-0.992123542028864,-1,-1,-1,-1
  )
)

PropQ2hPos <- cbind(
  c(
    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1
  ),
  c(
    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0.995,0.97,0.81,0.52,0.17,0.005,0,0,0,0
  )
)



lambda_optim <- list(
  c(TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE),
  c(TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE),
  c(FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE)
)

lambdaSparse <- c(0.5172414,0.4482759)



## ----loadPAckage,eval=FALSE---------------------------------------------------
#  devtools::install_github("hlorenzo/ddsPLS")
#  library(ddsPLS)

## ----loadPAckage2,echo=FALSE--------------------------------------------------
library(ddsPLS)

## ----noCor, eval=TRUE---------------------------------------------------------
eps <- 0.95
n <- 200
alpha <- 0.1
datas <- getData(n=n,alpha=alpha,beta_0=0.1,sigma=sqrt(1-eps^2),
                 p1=50,p2=25,p3=25,p=1000)

## ----noCor2, eval=TRUE--------------------------------------------------------
lambdas <- seq(0,1,length.out = 30)
n_B <- 200

## ----noCor3, eval=FALSE-------------------------------------------------------
#  NCORES <- 4
#  mo <- ddsPLS( datas$X, datas$Y,lambdas = lambdas,
#                n_B=n_B,NCORES=NCORES,verbose = FALSE)

## ----noCor3hide, echo=FALSE---------------------------------------------------
mo <- ddsPLS( datas$X, datas$Y,lambdas = lambdaSparse,
              n_B=n_B,NCORES=1,verbose = FALSE)

## ----summ---------------------------------------------------------------------
sum_up <- summary(mo,return = TRUE)

## ----selX---------------------------------------------------------------------
setdiff(1:75,mo$Selection$X)

## ----noCor30, eval=FALSE------------------------------------------------------
#  mo0 <- ddsPLS( datas$X, datas$Y,lambdas = 0,
#                 n_B=n_B,NCORES=NCORES,verbose = FALSE)
#  sum0 <- summary(mo0,return = TRUE)
#  print(sum0$R2Q2[,c(1,4)])

## -----------------------------------------------------------------------------
mo0 <- ddsPLS( datas$X, datas$Y,lambdas = 0,
               n_B=n_B,NCORES=1,verbose = FALSE)
sum0 <- summary(mo0,return = TRUE)

## -----------------------------------------------------------------------------
print(sum0$R2Q2[,c(1,4)])
print(sum_up$R2Q2[,c(1,4)])

## ----est,fig.width=5,fig.height=5,fig.align="center"--------------------------
plot(mo,type="predict",legend.position = "topleft")

## ----criterion,fig.width=7,fig.height=5,fig.align="center",eval=FALSE---------
#  plot(mo,type="criterion",legend.position = "top")

## ----criterionhide,fig.width=7,fig.height=5,fig.align="center",echo=FALSE-----
h_opt <- 2
matplot(lambdas,R2mean_diff_Q2mean,type = "l",ylab="",xlab=expression(lambda),
              main=bquote(bar(R)[B]^2-bar(Q)[B]^2))
for(s in 1:h_opt){
  points(lambdas,R2mean_diff_Q2mean[,s],type = "p",pch=16,cex=lambda_optim[[s]],col=s)
  points(lambdaSparse[s],R2mean_diff_Q2mean[which.min(abs(lambdas-lambdaSparse[s])),s],pch=1,cex=2,col=s)
}
legend("top",paste("Comp.",1:h_opt," (",round(mo$varExplained$Comp),"%)",sep=""),
       col = 1:h_opt,pch=16,bty = "n",
       title = paste("Total explained variance ",round(mo$varExplained$Cumu)[h_opt],"%",sep=""))

## ----Q2r,fig.width=7,fig.height=5,eval=FALSE----------------------------------
#  plot(mo,type="Q2",legend.position = "bottomleft")

## ----Q2rhide,fig.width=7,fig.height=5,echo=FALSE------------------------------
h_opt <- 2
matplot(lambdas,Q2mean,type = "l",ylab="",xlab=expression(lambda),
              main=bquote(bar(R)[B]^2-bar(Q)[B]^2))
for(s in 1:h_opt){
  points(lambdas,Q2mean[,s],type = "p",pch=16,cex=lambda_optim[[s]],col=s)
  points(lambdaSparse[s],Q2mean[which.min(abs(lambdas-lambdaSparse[s])),s],pch=1,cex=2,col=s)
}
abline(h=0,lwd=2,lty=2)
legend("bottomleft",paste("Comp.",1:h_opt," (",round(mo$varExplained$Comp),"%)",sep=""),
       col = 1:h_opt,pch=16,bty = "n",
       title = paste("Total explained variance ",round(mo$varExplained$Cumu)[h_opt],"%",sep=""))

## ----prop,fig.width=7,fig.height=5,fig.align="center",eval=FALSE--------------
#  plot(mo,type="prop",legend.position = "bottomleft")

## ----prophide,fig.width=7,fig.height=5,echo=FALSE-----------------------------
h_opt <- 2
# Plot of Prop of positive Q2h
matplot(lambdas,PropQ2hPos,type = "l",ylab="",xlab=expression(lambda),
        main=bquote("Proportion of models with positive"~Q["b,r"]^2))
abline(h=((1:10)/10)[-5],col="gray80",lwd=0.5,lty=3)
abline(h=5/10,col="gray60",lwd=0.7,lty=1)
text(min(lambdas),1/2,labels = "1/2",pos = 4,col="gray40")
for(s in 1:h_opt){
  points(lambdas,PropQ2hPos[,s],type = "p",pch=16,cex=mo$lambda_optim[[s]],col=s)
  points(lambdaSparse[s],PropQ2hPos[which.min(abs(lambdas-lambdaSparse[s])),s],pch=1,cex=2,col=s)
}
legend("bottomleft",paste("Comp.",1:h_opt," (",round(mo$varExplained$Comp),"%)",sep=""),
       col = 1:h_opt,pch=16,bty = "n",
       title = paste("Total explained variance ",round(mo$varExplained$Cumu)[h_opt],"%",sep=""))

## ----wy,fig.width=7,fig.height=3,fig.align="center"---------------------------
plot(mo,type="weightY",mar=c(4,7,2,1))

## ----wx2,fig.width=7,fig.height=3,fig.align="center"--------------------------
plot(mo,type="weightX",cex.names = 0.5 )

## ----getData2,echo=TRUE-------------------------------------------------------
getData

