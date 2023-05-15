#ifndef DDSPLS_H
#define DDSPLS_H

#include <RcppEigen.h>

using namespace Eigen;
using namespace std;

class stdMatTrainTest {
public:
  Eigen::MatrixXd MTrain;
  Eigen::MatrixXd MTest;
  Eigen::VectorXd sd;
  Eigen::VectorXd mu;
};

class oneComponent {
public:
  Eigen::VectorXd t;
  Eigen::VectorXd U0;
  Eigen::VectorXd V0;
  Eigen::VectorXd V_svd;
};

class multiComponent {
public:
  Eigen::MatrixXd U_out;
  Eigen::MatrixXd V0;
  Eigen::VectorXd explainedVar;
};

class ddsPLSCpp {
public:
  Eigen::VectorXd R2;
  Eigen::VectorXd R2h;
  Eigen::VectorXd Q2;
  Eigen::VectorXd Q2h;
  Eigen::VectorXd idIB;
  Eigen::VectorXd idOOB;
  Eigen::MatrixXd t;
  Eigen::MatrixXd U_star;
  Eigen::MatrixXd U;
  Eigen::MatrixXd V;
  Eigen::MatrixXd P;
  Eigen::MatrixXd C;
  Eigen::MatrixXd B;
  double lambda0;
};

class ddsPLSCpp_B {
public:
  Eigen::MatrixXd R2;
  Eigen::MatrixXd R2h;
  Eigen::MatrixXd Q2;
  Eigen::MatrixXd Q2h;
};

class megaResults {
public:
  Eigen::MatrixXd R2;
  Eigen::MatrixXd R2h;
  Eigen::MatrixXd Q2;
  Eigen::MatrixXd Q2h;
  Eigen::VectorXd R2mean;
  Eigen::VectorXd R2hmean;
  Eigen::VectorXd Q2mean;
  Eigen::VectorXd Q2hmean;
  Eigen::VectorXd R2mean_diff_Q2mean;
};

class ddsPLSBoot {
public:
  multiComponent model;
  int R_sol;
  Eigen::VectorXd explainedVar;
  Eigen::VectorXd lambda_sol;
  Eigen::VectorXd Q2;
  Eigen::VectorXd Q2h;
  Eigen::VectorXd R2;
  Eigen::VectorXd R2h;
  megaResults results;
};

#endif
