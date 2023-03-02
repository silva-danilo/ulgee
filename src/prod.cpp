// packages
#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;
using namespace Eigen;
using namespace std;

//' @export
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Eigen::VectorXd prod_1(const Eigen::MatrixXd M, Eigen::VectorXd m){
  return M*m;
}

// [[Rcpp::export]]
Eigen::MatrixXd prod_2(const Eigen::MatrixXd R, Eigen::VectorXd sigma, 
                       Eigen::VectorXd D){
  MatrixXd Om = sigma.asDiagonal()*R*sigma.asDiagonal();
  return D.asDiagonal()*Om.inverse()*D.asDiagonal();
}

// [[Rcpp::export]]
Eigen::MatrixXd prod_3(const Eigen::MatrixXd X, Eigen::MatrixXd W){
  return X.transpose()*W*X;
}

// [[Rcpp::export]]
Eigen::VectorXd prod_4(const Eigen::MatrixXd S1, Eigen::VectorXd s2){
  return S1.inverse()*s2;
}

// [[Rcpp::export]]
Eigen::MatrixXd prod_5(const Eigen::MatrixXd S1, Eigen::MatrixXd S3){
  return S1.inverse()*S3*S1.inverse();
}
