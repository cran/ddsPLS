#include "ddsplsCore.h"
#include <RcppEigen.h>
#include <random>
#include <chrono>
#include <thread>

using namespace Eigen;
using namespace std;

stdMatTrainTest standardizeCpp(Eigen::MatrixXd MIn,Eigen::MatrixXd MOutTestStd,
                               const int p,const int n,const int countOOB){
  Eigen::VectorXd colD(n);
  Eigen::VectorXd MU(n);
  Eigen::VectorXd sdOut(p);
  double mui=0.0;
  double sdi=0.0;
  for (int i = 0u; i < p; ++i) {
    colD = MIn.block(0,i,n,1);
    mui = colD.mean();
    for(int k = 0; k < n; k++){
      MU(k)     = mui;
    }
    sdi = sqrt((colD - MU).squaredNorm()/(n-1));
    sdOut(i) = sdi;
    for(int k = 0u; k < n; ++k) {
      MIn(k,i) -= mui;
      if(sdi>0){
        MIn(k,i) /= sdi;
      }else{
        MIn(k,i) = 0;
      }
    }
    for(int k = 0u; k < countOOB; ++k) {
      MOutTestStd(k,i) -= mui;
      if(sdi>0){
        MOutTestStd(k,i) /= sdi;
      }else{
        MOutTestStd(k,i) = 0;
      }
    }
  }
  stdMatTrainTest out;
  out.MTrain = MIn;
  out.MTest = MOutTestStd;
  out.sd = sdOut;
  return out;
}

stdMatTrainTest standardizeCppNoTest(Eigen::MatrixXd MIn,
                                     const int p,const int n){
  Eigen::VectorXd colD(n);
  Eigen::VectorXd MU(n);
  Eigen::VectorXd sdOut(p),muOut(p);
  double mui=0.0;
  double sdi=0.0;
  for (int i = 0u; i < p; ++i) {
    colD = MIn.block(0,i,n,1);
    mui = colD.mean();
    for(int k = 0; k < n; k++){
      MU(k)     = mui;
    }
    sdi = sqrt((colD - MU).squaredNorm()/(n-1));
    sdOut(i) = sdi;
    muOut(i) = mui;
    for(int k = 0u; k < n; ++k) {
      MIn(k,i) -= mui;
      if(sdi>0){
        MIn(k,i) /= sdi;
      }else{
        MIn(k,i) = 0;
      }
    }
  }
  stdMatTrainTest out;
  out.MTrain = MIn;
  out.sd = sdOut;
  out.mu = muOut;
  return out;
}

oneComponent do_one_componentCpp(const Eigen::MatrixXd x0,const Eigen::MatrixXd y0,
                                 const Eigen::MatrixXd COV,
                                 const int n,const int p,const int q,
                                 const double lam,const double errorMin=1e-9){
  // Create random generator
  std::random_device rd2;
  std::mt19937 generator2(rd2());
  std::normal_distribution<> d2(0, 1);
  double valueRd=1.0;

  Eigen::VectorXd max_cov_y = Eigen::VectorXd::Zero(q);
  Eigen::VectorXd max_cov_x = Eigen::VectorXd::Zero(p);
  VectorXi id_y_high = VectorXi::Zero(q);
  int countNoNullY=0;
  VectorXi id_x_high = VectorXi::Zero(p);
  int countNoNullX=0;
  Eigen::MatrixXd COVIn = COV;
  // Get max values per column and per row
  for (int i = 0u; i < q; ++i) {
    max_cov_y(i) = COV.block(i,0,1,p).lpNorm<Infinity>();
    if (max_cov_y(i)>lam){
      id_y_high(i) = 1;
      countNoNullY += 1;
    }
  }
  for (int j = 0u; j < p; ++j) {
    max_cov_x(j) = COV.block(0,j,q,1).lpNorm<Infinity>();
    if (max_cov_x(j)>lam){
      id_x_high(j) = 1;
      countNoNullX += 1;
    }
  }
  // Get reduced covariance matrix
  Eigen::MatrixXd COV_high = Eigen::MatrixXd::Zero(countNoNullY,countNoNullX);
  Eigen::VectorXd U0 = Eigen::VectorXd::Zero(p);
  Eigen::VectorXd V0 = Eigen::VectorXd::Zero(q);
  Eigen::VectorXd V_svd = Eigen::VectorXd::Zero(q);
  Eigen::VectorXd t = Eigen::VectorXd::Zero(n);
  double coefIJ=0.0,value=0.0;
  int countY = 0u;
  int countX = 0u;
  if (countNoNullY>0){
    for (int i = 0u; i < q; ++i) {
      countX = 0u;
      if (id_y_high(i)==1){
        countX = 0u;
        for (int j = 0u; j < p; ++j) {
          if (id_x_high(j)==1){
            coefIJ = COV(i,j);
            value = abs(coefIJ)-lam;
            if (value>0) {
              if (coefIJ>0){
                COV_high(countY,countX) = value;
              } else {
                COV_high(countY,countX) = value*(-1.0);
              }
            }
            else {
              COV_high(countY,countX) = 0;
            }
            countX += 1;
          }
        }
        countY += 1;
      }
    }
    double error = 2.0;
    Eigen::VectorXd u0In = Eigen::VectorXd::Zero(countNoNullX);
    Eigen::VectorXd u0In2 = Eigen::VectorXd::Zero(countNoNullX);
    Eigen::VectorXd v0In = Eigen::VectorXd::Zero(countNoNullY);
    for (int j = 0u; j < countNoNullX; ++j) {
      valueRd = d2(generator2);
      u0In(j) = valueRd;
    }
    u0In /= sqrt(u0In.squaredNorm());

    while(error>errorMin){
      v0In = COV_high*u0In;
      u0In2 = COV_high.transpose()*v0In;
      u0In2 /= sqrt(u0In2.squaredNorm());
      error = (u0In2-u0In).squaredNorm();
      u0In = u0In2;
    }
    v0In  /= sqrt(v0In.squaredNorm());
    countY = 0u;
    for (int i = 0u; i < q; ++i) {
      if (id_y_high(i)==1){
        V0(i) = v0In(countY);
        countY += 1;
      }
    }
    countX = 0u;
    for (int j = 0u; j < p; ++j) {
      if (id_x_high(j)==1){
        U0(j) = u0In(countX);
        countX += 1;
      }
    }
    // Build score
    t = x0*U0;
    // Build y0 masked
    for (int i = 0u; i < q; ++i) {
      if (id_y_high(i)==1){
        V_svd(i) = (y0.block(0,i,n,1).transpose()*t).sum();
      }
    }
    V_svd /= sqrt(V_svd.squaredNorm());
  }
  // Generate output
  oneComponent out;
  out.t = t;
  out.U0 = U0;
  out.V_svd = V_svd;
  out.V0 = V0;
  return(out);
}

multiComponent modelddsPLSCpp(Eigen::MatrixXd U_out,Eigen::MatrixXd V0,
                              const Eigen::MatrixXd COVInit,const double maxCOVInit,
                              const Eigen::MatrixXd x, const Eigen::MatrixXd y,
                              const int n,const int p,const int q,const Eigen::VectorXd lam,
                              const int R=1,double errorMin=1e-9){
  double RSS0=0.0;
  double normU02;
  double RSSr;
  double maxCOV=0.0;
  Eigen::MatrixXd muX = Eigen::MatrixXd::Zero(n,p);
  Eigen::MatrixXd sdXInvMat = Eigen::MatrixXd::Zero(p,q);
  Eigen::MatrixXd muY = Eigen::MatrixXd::Zero(n,q);
  Eigen::MatrixXd sdYMat = Eigen::MatrixXd::Zero(p,q);
  Eigen::MatrixXd sdYXInvMat = Eigen::MatrixXd::Zero(p,q);
  Eigen::MatrixXd y_plus_un = Eigen::MatrixXd::Zero(n,q);
  Eigen::MatrixXd x_plus_un = Eigen::MatrixXd::Zero(n,p);
  Eigen::MatrixXd U_star = Eigen::MatrixXd::Zero(p,R);
  Eigen::MatrixXd V_out = Eigen::MatrixXd::Zero(q,R);
  Eigen::MatrixXd bXr = Eigen::MatrixXd::Zero(R,p);
  Eigen::MatrixXd bYr = Eigen::MatrixXd::Zero(R,q);
  Eigen::MatrixXd y_est = Eigen::MatrixXd::Zero(n,q);
  Eigen::MatrixXd B_r_0 = Eigen::MatrixXd::Zero(p,q);
  Eigen::MatrixXd y0 = y;
  Eigen::MatrixXd x0 = x;
  Eigen::MatrixXd COV = Eigen::MatrixXd::Zero(q,p);
  Eigen::VectorXd vectHere = Eigen::VectorXd::Zero(n);
  Eigen::VectorXd var_expl = Eigen::VectorXd::Zero(R);
  Eigen::VectorXd t_r = Eigen::VectorXd::Zero(n);
  Eigen::VectorXd U0 = Eigen::VectorXd::Zero(p);
  Eigen::VectorXd V_svd = Eigen::VectorXd::Zero(q);
  Eigen::VectorXd V0_r = Eigen::VectorXd::Zero(q);
  Eigen::VectorXd bt = Eigen::VectorXd::Zero(p);
  Eigen::VectorXd deltaU = Eigen::VectorXd::Zero(p);
  // Compute initial residual sum of squares
  RSS0 = y0.squaredNorm();
  double lambda_h;
  // Begin to build subspaces
  for (int r = 0u; r < R; ++r) {
    // Build empirical covariance matrix
    if (r==0){
      COV = COVInit;
      maxCOV = maxCOVInit;
    }else{
      COV = y0.transpose()*x0/double(n-1.0);
      maxCOV = COV.lpNorm<Infinity>();
    }
    lambda_h = lam[r];
    if(maxCOV>lambda_h){
      oneComponent c_h = do_one_componentCpp(x0,y0,COV,n,p,q,lambda_h,errorMin);
      t_r = c_h.t;
      U0 = c_h.U0;
      U_out.block(0,r,p,1) = U0;
      V_svd = c_h.V_svd;
      V0_r = c_h.V0;
      V0.block(0,r,q,1) = V0_r;
      // Build regression matrices
      normU02 = U0.squaredNorm();
      if(normU02>errorMin){
        bt = x0.transpose()*t_r/normU02;
        x_plus_un = t_r*bt.transpose();
        U_star.block(0,r,p,1) = U0;
        V_out.block(0,r,q,1) = V_svd;
        B_r_0 = U0*V_svd.transpose();
        y_plus_un = t_r*V_svd.transpose();
        y_est += y_plus_un;
        bXr.block(r,0,1,p) = bt.transpose();
        if(r>0){
          for (int s_r = r-1; s_r > 0; --s_r){
            deltaU = U_star.block(0,s_r,p,1)*(bXr.block(s_r,0,1,p)*U_star.block(0,s_r,p,1)).sum();
            U_star.block(0,r,p,1) = U_star.block(0,r,p,1) - deltaU;
          }
        }
        // Computation of explained variance
        Eigen::MatrixXd diff_here = y-y_est;
        RSSr = diff_here.squaredNorm();
        var_expl(r) = 1.0-RSSr/RSS0;
        y0 -= y_plus_un;
        x0 -= x_plus_un;
      }
    }
  }
  multiComponent out;
  out.U_out = U_out;
  out.V0 = V0;
  out.explainedVar = var_expl;
  return out;
}

ddsPLSCpp bootstrap_pls_CT_Cpp(const Eigen::MatrixXd X_init,const Eigen::MatrixXd Y_init,
                               const Eigen::VectorXd lambdas,const Eigen::VectorXd lambda_prev,
                               Eigen::MatrixXd uIN, Eigen::MatrixXd vIN,
                               const int n,const int p,const int q,const int N_lambdas,
                               const Eigen::VectorXd lambda0,
                               const bool doBoot = true,
                               const int h=1,
                               const double errorMin=1.0e-9){
  Eigen::VectorXd idIB = Eigen::VectorXd::Zero(n);
  Eigen::VectorXd idIBChosen = Eigen::VectorXd::Zero(n);
  Eigen::VectorXd idOOBChosen(n);
  for (int i=0; i<n; ++i) {
    idOOBChosen(i) = 1;
  }
  Eigen::MatrixXd X_train = Eigen::MatrixXd::Zero(n,p);
  Eigen::MatrixXd Y_train = Eigen::MatrixXd::Zero(n,q);
  int countOOB = 0;
  int nbOOB=0;
  bool test = true;
  Eigen::MatrixXd u = Eigen::MatrixXd::Zero(p,h);
  Eigen::MatrixXd v = Eigen::MatrixXd::Zero(q,h);
  Eigen::MatrixXd V_reconstruct = Eigen::MatrixXd::Zero(q,h);
  Eigen::MatrixXd t_all = Eigen::MatrixXd::Zero(n,h);
  Eigen::MatrixXd X_r = Eigen::MatrixXd::Zero(n,p);
  Eigen::MatrixXd Y_r = Eigen::MatrixXd::Zero(n,q);
  Eigen::MatrixXd Y_r_mask = Eigen::MatrixXd::Zero(n,q);
  Eigen::MatrixXd x_r = Eigen::MatrixXd::Zero(n,p);
  Eigen::MatrixXd y_r = Eigen::MatrixXd::Zero(n,q);
  Eigen::MatrixXd B_youyou = Eigen::MatrixXd::Zero(p,q);
  ddsPLSCpp out;
  out.B = Eigen::MatrixXd::Zero(N_lambdas*p,q);
  out.U = Eigen::MatrixXd::Zero(1,N_lambdas*p);
  out.U_star = Eigen::MatrixXd::Zero(1,N_lambdas*p);
  out.V = Eigen::MatrixXd::Zero(1,N_lambdas*q);
  out.P = Eigen::MatrixXd::Zero(1,N_lambdas*p);
  out.C = Eigen::MatrixXd::Zero(1,N_lambdas*q);
  out.t = Eigen::MatrixXd::Zero(1,N_lambdas*n);
  int r=0;
  double sdyi=0.0,sdxj=0.0;
  Eigen::MatrixXd diff_B(p,q);
  Eigen::MatrixXd y_train_pred(n,q);
  Eigen::MatrixXd y_train_pred_next(n,q);
  double n_t_p_i,n_t_p_n_i,d_t_i;
  Eigen::VectorXd sdY(q),sdX(p),muY(q),muX(p);
  Eigen::VectorXd idOOB;
  // Build bootstrap indexes, IB and OOB
  if(doBoot==true){
    while (test){
      std::random_device rd;
      std::mt19937 generator(rd());
      std::uniform_int_distribution<> distribution(0, n-1);
      for (int i=0; i<n; ++i) {
        int number = distribution(generator);
        idIB(i) = number;
        idIBChosen(number) = 1;
        idOOBChosen(number) = 0;
      }
      nbOOB = idOOBChosen.sum();
      if(nbOOB>0){
        test = false;
      }
    }
  } else {
    nbOOB = 1;
  }
  idOOB = Eigen::VectorXd(nbOOB);
  if(doBoot==true){
    countOOB = 0;
    for (int k = 0u; k < n; ++k) {
      if (idOOBChosen(k)==1){
        idOOB(countOOB) = k;
        countOOB += 1;
      }
    }
    // Build train matrices
    for (int k = 0u; k < n; ++k) {
      X_train.block(k,0,1,p) = X_init.block(idIB(k),0,1,p);
      Y_train.block(k,0,1,q) = Y_init.block(idIB(k),0,1,q);
    }
    // Build test matrices
  } else {
    countOOB = 1;
  }
  Eigen::MatrixXd X_test_normalize(countOOB,p);
  Eigen::MatrixXd Y_test_normalize(countOOB,q);
  Eigen::MatrixXd y_test_pred(countOOB,q);
  Eigen::MatrixXd y_test_pred_RSS(countOOB,q);
  double n_Q2_i,d_Q2_i,d_Q2_a_i;
  // Standardize X and Y train and test.
  if(doBoot==true){
    for (int k = 0u; k < countOOB; ++k) {
      X_test_normalize.block(k,0,1,p) = X_init.block(idOOB(k),0,1,p);
      Y_test_normalize.block(k,0,1,q) = Y_init.block(idOOB(k),0,1,q);
    }
    stdMatTrainTest stdY = standardizeCpp(Y_train,Y_test_normalize,q,n,countOOB);
    Y_r = stdY.MTrain;
    sdY = stdY.sd;
    stdY = stdMatTrainTest();
    stdMatTrainTest stdX = standardizeCpp(X_train,X_test_normalize,p,n,countOOB);
    X_train = stdX.MTrain;
    X_r = stdX.MTrain;
    X_test_normalize = stdX.MTest;
    sdX = stdX.sd;
    stdX = stdMatTrainTest();
  } else {
    X_train = X_init;
    Y_train = Y_init;
    stdMatTrainTest stdY = standardizeCppNoTest(Y_train,q,n);
    Y_train = stdY.MTrain;
    Y_r = stdY.MTrain;
    muY = stdY.mu;
    sdY = stdY.sd;
    stdY = stdMatTrainTest();
    stdMatTrainTest stdX = standardizeCppNoTest(X_train,p,n);
    X_train = stdX.MTrain;
    X_r = stdX.MTrain;
    muX = stdX.mu;
    sdX = stdX.sd;
    stdX = stdMatTrainTest();
  }

  // Build model on past components if needed
  Eigen::MatrixXd u_r = Eigen::MatrixXd::Zero(p,1);
  Eigen::MatrixXd v_r = Eigen::MatrixXd::Zero(q,1);
  Eigen::VectorXd t_r = Eigen::VectorXd::Zero(n);
  Eigen::VectorXd bt = Eigen::VectorXd::Zero(p);
  Eigen::MatrixXd t_out = Eigen::MatrixXd::Zero(n,n);
  Eigen::MatrixXd C = Eigen::MatrixXd::Zero(p,h);
  Eigen::MatrixXd D = Eigen::MatrixXd::Zero(q,h);
  Eigen::MatrixXd COV(q,p);
  double norm2t=0.0;
  bool test_previous_ok = true;
  bool test_no_null=true;
  double maxCOV;
  Eigen::VectorXd lam_r = Eigen::VectorXd(1);
  // Need to build the h-1 components for the current bootstrap dataset
  if (h>1){
    for (int s = 0u; s < h-1; ++s) {
      u.block(0,s,p,1) = uIN.block(0,s,p,1);
      v.block(0,s,q,1) = vIN.block(0,s,q,1);
    }
    r=0;
    test=true;
    while (test){
      Eigen::MatrixXd onlyToInv(r+1,r+1);
      Eigen::MatrixXd InversedMat(r+1,r+1);
      Eigen::MatrixXd U_star_cl(p,r);
      COV = Y_r.transpose()*X_r/double(n-1.0);
      maxCOV = COV.lpNorm<Infinity>();
      lam_r(0) = lambda_prev(r);
      multiComponent res1 = modelddsPLSCpp(u,v,COV,maxCOV,X_r,Y_r,n,p,q,lam_r);
      v_r = res1.V0;
      u_r = res1.U_out;
      u.block(0,r,p,1) = u_r.block(0,0,p,1);
      v.block(0,r,q,1) = v_r.block(0,0,q,1);
      t_r = X_r*u_r;
      t_out.block(0,r,n,1) = t_r;
      norm2t = (t_r.transpose()*t_r).sum();
      if (norm2t>errorMin) {
        // Build X regressors and estimated matrix
        bt = Eigen::VectorXd(p);
        bt = X_r.transpose()*t_r/norm2t;
        C.block(0,r,p,1) = bt;
        x_r = t_r*bt.transpose();
        // Build Y regressors and estimated matrix
        for (int i = 0u; i < q; ++i) {
          if (abs(v_r(i,0))<errorMin){
            D(i,r) = 0.0;
          } else {
            D(i,r) = (Y_r.block(0,i,n,1).transpose()*t_r).sum()/norm2t;
          }
          y_r.block(0,i,n,1) = t_r*D(i,r);
        }
        // Only matrix inversion of the problem
        onlyToInv = C.block(0,0,p,r+1).transpose()*u.block(0,0,p,r+1);
        InversedMat = onlyToInv.inverse();
        U_star_cl = u.block(0,0,p,r+1)*InversedMat;
        B_youyou = U_star_cl*D.block(0,0,q,r+1).transpose();
        for (int i = 0u; i < q; ++i) {
          sdyi=sdY(i);
          if(sdyi>errorMin){
            for (int j = 0u; j < p; ++j) {
              sdxj=sdX(j);
              B_youyou(j,i) *= sdyi;
              if(sdxj>errorMin){
                B_youyou(j,i) /= sdxj;
              }else{
                B_youyou(j,i) = 0.0;
              }
            }
          }else{
            for (int j = 0u; j < p; ++j) {
              B_youyou(j,i) = 0.0;
            }
          }
        }
        X_r -= x_r;
        Y_r -= y_r;
      } else {
        test_no_null = false;
        test_previous_ok = false;
      }
      if (r==h-2) {
        test = false;
      }
      // Update r
      r += 1;
    }
  } else {
    norm2t = 1.0;
    test_previous_ok = true;
  }
  // Begin to test each lambda
  r = h-1;
  Eigen::MatrixXd onlyToInv(h,h);
  Eigen::MatrixXd InversedMat(h,h);
  Eigen::MatrixXd B_all = Eigen::MatrixXd::Zero(p,q);
  Eigen::MatrixXd U_star_cl(p,h);
  Eigen::MatrixXd u_il = Eigen::MatrixXd::Zero(p,1);
  Eigen::MatrixXd V_il = Eigen::MatrixXd::Zero(q,1);
  double coeffTest = 0.0;
  bool testGood = false;
  COV = Y_r.transpose()*X_r/double(n-1.0);
  maxCOV = COV.lpNorm<Infinity>();
  int N_simu_lams = N_lambdas;
  if (doBoot == false){
    int N_simu_lams = 1;
  }

  Eigen::VectorXd vars_expl(N_simu_lams), vars_expl_h(N_simu_lams), Q2(N_simu_lams), Q2_all(N_simu_lams);
  for (int iLam = 0u; iLam < N_simu_lams; ++iLam){
    if ( (test_previous_ok==true) & (lambdas(iLam)>=lambda0(r)) ) {
      lam_r(0) = lambdas(iLam);
      if (doBoot == false){
        lam_r(0) = lambdas(0);
      }
      multiComponent res1 = modelddsPLSCpp(u,v,COV,maxCOV,X_r,Y_r,n,p,q,lam_r);
      u_il = res1.U_out;
      V_il = res1.V0;
      u.block(0,r,p,1) = u_il.block(0,0,p,1);
      v.block(0,r,q,1) = V_il.block(0,0,q,1);
      t_r = X_r*u_il.block(0,0,p,1);
      t_out.block(0,r,n,1) = t_r;
      norm2t = (t_r.transpose()*t_r).sum();
      if (norm2t>errorMin) {
        testGood = true;
        // Build X regressors and estimated matrix
        bt = X_r.transpose()*t_r/norm2t;
        C.block(0,r,p,1) = bt;
        x_r = t_r*bt.transpose();
        // Build Y regressors and estimated matrix
        for (int i = 0u; i < q; ++i) {
          coeffTest = abs(V_il(i,0));
          if (coeffTest<errorMin){
            D(i,r) = 0.0;
          } else {
            D(i,r) = (Y_r.block(0,i,n,1).transpose()*t_r).sum()/norm2t;
          }
        }
        onlyToInv = C.transpose()*u;
        InversedMat = onlyToInv.inverse();
        U_star_cl = u*InversedMat;
        B_all = U_star_cl*D.transpose();
        for (int i = 0u; i < q; ++i) {
          sdyi=sdY(i);
          if(sdyi>errorMin){
            for (int j = 0u; j < p; ++j) {
              sdxj=sdX(j);
              B_all(j,i) *= sdyi;
              if(sdxj>errorMin){
                B_all(j,i) /= sdxj;
              }else{
                B_all(j,i) = 0.0;
              }
            }
          }else{
            for (int j = 0u; j < p; ++j) {
              B_all(j,i) = 0.0;
            }
          }
        }
        out.B.block(iLam*p,0,p,q) = B_all;
        out.t.block(0,iLam*n,1,n) = t_r.block(0,0,n,1).transpose();
        out.U.block(0,iLam*p,1,p) = u_il.block(0,0,p,1).transpose();
        out.U_star.block(0,iLam*p,1,p) = U_star_cl.transpose();
        out.V.block(0,iLam*q,1,q) = V_il.block(0,0,q,1).transpose();
        out.P.block(0,iLam*p,1,p) = bt.transpose();
        out.C.block(0,iLam*q,1,q) = D.block(0,r,q,1).transpose();
        // Create the prediction matrices and metrics
        diff_B = B_all - B_youyou;
        y_train_pred = X_train*B_all;
        y_train_pred_next = X_train*diff_B;
        n_t_p_i = (y_train_pred-Y_train).squaredNorm();
        n_t_p_n_i = (y_train_pred_next-Y_train).squaredNorm();
        d_t_i = (Y_train).squaredNorm();
        vars_expl(iLam) = 1.0-n_t_p_i/d_t_i;
        vars_expl_h(iLam) = 1.0-n_t_p_n_i/d_t_i;
        if(doBoot==true){
          y_test_pred = X_test_normalize*B_all;
          y_test_pred_RSS = X_test_normalize*B_youyou;
          n_Q2_i = (y_test_pred-Y_test_normalize).squaredNorm();
          d_Q2_i = (y_test_pred_RSS-Y_test_normalize).squaredNorm();
          d_Q2_a_i = (Y_test_normalize).squaredNorm();
          Q2(iLam) = 1.0-n_Q2_i/d_Q2_i;
          Q2_all(iLam) = 1.0-n_Q2_i/d_Q2_a_i;
        }
      }
    }
    if (testGood == false){
      vars_expl(iLam) = -1.0;
      vars_expl_h(iLam) = -1.0;
      if(doBoot==true){
        Q2(iLam) = -1.0;
        Q2_all(iLam) = -1.0;
      }
    }
    testGood = false;
  }
  out.R2 = vars_expl; // R^2
  out.R2h = vars_expl_h; // R^2_h
  if (doBoot==true) {
    out.Q2 = Q2_all; // Q^2
    out.Q2h = Q2; // Q^2_h
  } else {
    out.t = t_out.block(0,0,n,h);
    out.P = C.block(0,0,p,h);
    out.C = D.block(0,0,q,h);
    out.U_star = U_star_cl;
    out.U = u;
    out.V = v;
    out.B = B_all;
  }
  return out;
}


//' @title C++ code to build models, internal function
//' @description
//' Build a ddsPLS model once the bootstrap operations has allowed to find a correct lambda.
//' @param U The weights for X part.
//' @param V The weights for Y part.
//' @param X The matrix of X part.
//' @param Y The matrix of X part.
//' @param lambdas The to be tested values for lambda.
//' @param R The number of components to build.
//' @param n The number of observations.
//' @param p The number of variables of X part.
//' @param q The number of variables of Y part.
//' @param lambda0 The lowest value to be tested for lambda.
//'
//' @return A list containing the PLS parameters:
//' \itemize{
//'   \item \code{$P}: Loadings for \code{X}.
//'   \item \code{$C}: Loadings for \code{Y}.
//'   \item \code{$t}: Scores.
//'   \item \code{$V}: Weights for \code{Y}.
//'   \item \code{$U}: Loadings for \code{X}.
//'   \item \code{$U_star}: Loadings for \code{X} in original base: $U_star=U(P'U)^{-1}$.
//'   \item \code{$B}: Regression matrix of \code{Y} on \code{X}.
//'   \item \code{$muY}: Empirical mean of \code{Y}.
//'   \item \code{$muX}: Empirical mean of \code{X}.
//'   \item \code{$sdY}: Empirical standard deviation of \code{Y}.
//'   \item \code{$sdX}: Empirical standard deviation of \code{X}.
//'}
//'
// [[Rcpp::export]]
Rcpp::List  modelddsPLSCpp_Rcpp(const Eigen::MatrixXd U,const Eigen::MatrixXd V,
                                const Eigen::MatrixXd X, const Eigen::MatrixXd Y,
                                const Eigen::VectorXd lambdas,const int R,
                                const int n,const int p,const int q,
                                const Eigen::VectorXd lambda0){
  Eigen::VectorXd lambda_prev(R-1);
  for (int r = 0u; r < R-1; ++r) {
    lambda_prev(r) = lambdas(r);
  }
  Eigen::VectorXd lambda_next(1);
  lambda_next(0) = lambdas(R-1);
  ddsPLSCpp res = bootstrap_pls_CT_Cpp(X,Y,lambda_next,lambda_prev,
                                       U,V,n,p,q,1,lambda0,false,R);
  Rcpp::List out;
  out["P"] = res.P;
  out["C"] = res.C;
  out["t"] = res.t;
  out["V"] = res.V;
  out["U"] = res.U;
  out["U_star"] = res.U_star;
  out["B"] = res.B;
  return out;
}

//' @title C++ implementation of the bootstrap operations
//' @description
//' Start the bootstrap operations.
//' Should not be used by user.
//' @param U The weights for X part.
//' @param V The weights for Y part.
//' @param X The matrix of X part.
//' @param Y The matrix of X part.
//' @param lambdas The to be tested values for lambda.
//' @param lambda_prev The previously selected values for lambda.
//' @param R The number of components to build.
//' @param n_B The number of bootstrap samples to generate and analyse.
//' @param doBoot Wheteher do bootstrap operations.
//' @param n The number of observations.
//' @param p The number of variables of X part.
//' @param q The number of variables of Y part.
//' @param N_lambdas The number of to be tested values for lambda.
//' @param lambda0 The minimum value to be checked in lambdas.
//'
//' @return The bootstrapped statistics
//'
// [[Rcpp::export]]
Rcpp::List  bootstrap_Rcpp(const Eigen::MatrixXd U,const Eigen::MatrixXd V,
                           const Eigen::MatrixXd X,const Eigen::MatrixXd Y,
                           const Eigen::VectorXd lambdas,const Eigen::VectorXd lambda_prev,
                           const int R,const int n_B,const bool doBoot,
                           const int n,const int p,const int q,const int N_lambdas,
                           const Eigen::VectorXd lambda0){
  Eigen::MatrixXd R2 = Eigen::MatrixXd(n_B,N_lambdas);
  Eigen::MatrixXd R2h = Eigen::MatrixXd(n_B,N_lambdas);
  Eigen::MatrixXd Q2 = Eigen::MatrixXd(n_B,N_lambdas);
  Eigen::MatrixXd Q2h = Eigen::MatrixXd(n_B,N_lambdas);
  Eigen::MatrixXd Uout = Eigen::MatrixXd::Zero(n_B,N_lambdas*p);
  Eigen::MatrixXd tt = Eigen::MatrixXd::Zero(n_B,N_lambdas*n);
  Eigen::MatrixXd Ustarout = Eigen::MatrixXd::Zero(n_B,N_lambdas*p);
  Eigen::MatrixXd Vout = Eigen::MatrixXd::Zero(n_B,N_lambdas*q);
  Eigen::MatrixXd P = Eigen::MatrixXd::Zero(n_B,N_lambdas*p);
  Eigen::MatrixXd C = Eigen::MatrixXd::Zero(n_B,N_lambdas*q);
  Eigen::MatrixXd B = Eigen::MatrixXd::Zero(N_lambdas*p,q*n_B);
  // Rcpp::List B;
  ddsPLSCpp res;
  for (int i = 0u; i < n_B; ++i) {
    res = bootstrap_pls_CT_Cpp(X,Y,lambdas,lambda_prev,U,V,n,p,q,N_lambdas,lambda0,doBoot,R);
    R2.block(i,0,1,N_lambdas) = res.R2.transpose();
    R2h.block(i,0,1,N_lambdas) = res.R2h.transpose();
    Q2.block(i,0,1,N_lambdas) = res.Q2.transpose();
    Q2h.block(i,0,1,N_lambdas) = res.Q2h.transpose();
    tt.block(i,0,1,n*N_lambdas) = res.t.transpose();
    Uout.block(i,0,1,p*N_lambdas) = res.U.transpose();
    Ustarout.block(i,0,1,p*N_lambdas) = res.U_star.transpose();
    Vout.block(i,0,1,q*N_lambdas) = res.V.transpose();
    P.block(i,0,1,p*N_lambdas) = res.P.transpose();
    C.block(i,0,1,q*N_lambdas) = res.C.transpose();
    B.block(0,i*q,p*N_lambdas,q) = res.B;
  }
  Rcpp::List out;
  out["R2"] = R2;
  out["R2h"] = R2h;
  out["Q2"] = Q2;
  out["Q2h"] = Q2h;
  out["lambdas"] = lambdas;
  out["t"] = tt;
  out["U"] = Uout;
  out["U_star"] = Ustarout;
  out["V"] = Vout;
  out["P"] = P;
  out["C"] = C;
  out["B"] = B;
  out["lambdas"] = lambdas;
  return out;
}
