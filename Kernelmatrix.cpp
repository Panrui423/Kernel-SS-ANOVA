#include<Rmath.h>
#include<RcppCommon.h>
#include<RcppArmadillo.h>
#include<RcppEigen.h>

// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]
using namespace Eigen;
using namespace arma;
using Eigen::Map;
using Eigen::MatrixXd;
using Rcpp::as;


// [[Rcpp::export]]
VectorXd Map_Expball(const double d1, const VectorXd d2){
  const int n = d2.rows();
  VectorXd d = (d2.array() - d1).array().exp();
  return ( 1 - ( 1 + d.array() ).array().inverse());  
}

// [[Rcpp::export]]
extern "C" SEXP MetricKern(const MatrixXd& dx){
  const int n = dx.rows();
  MatrixXd out = MatrixXd::Zero(n, n);
  
  for (int i = 0; i < n; i++){
    for (int j = 0; j <= i; j++){
      double element = 0.0;
      
      for (int k = 0; k < n; k++){
        VectorXd dx_row = dx.row(k);
        VectorXd vec1 = Map_Expball( dx(i, k), dx_row);
        VectorXd vec2 = Map_Expball( dx(j, k), dx_row);
        element += (vec1.transpose() * vec2);
      }
      
      out(i, j) = element;
      out(j, i) = out(i, j);
    }
  }
  
  return Rcpp::wrap(out / (n*n));
}

// [[Rcpp::export]]
extern "C" SEXP MetricKernSim(const MatrixXd& dx, const MatrixXd& dx1, const VectorXd loc){
  const int n = dx.rows();
  const int n1 = dx1.rows();
  Rcpp::IntegerVector  lab(n1);
  int m;
  MatrixXd out = MatrixXd::Zero(n, n);
  
  lab = loc.cast<int>();
  
  for (int i = 0; i < n; i++){
    for (int j = 0; j <= i; j++){
      double element = 0.0;
      
      for (int k = 0; k < n1; k++){
        VectorXd dx1_row = dx1.row(k);
        m = lab[k];
        VectorXd vec1 = Map_Expball( dx(i, m), dx1_row);
        VectorXd vec2 = Map_Expball( dx(j, m), dx1_row);
        element += (vec1.transpose() * vec2);
      }
      
      out(i, j) = element;
      out(j, i) = out(i, j);
    }
  }
  
  return Rcpp::wrap(out / (n1 * n1));
}


// [[Rcpp::export]]
extern "C" Rcpp::IntegerVector test(const MatrixXd& dx, const MatrixXd& dx1, const VectorXd id){
  const int n1 = dx1.rows();
  Rcpp::IntegerVector  lab(n1);
  int m;
  
  lab = id.cast<int>();
  return lab;
}


// [[Rcpp::export]]
extern "C" SEXP MetricKtest(const MatrixXd& dxtest, int train_num, int test_num) {
  MatrixXd out = MatrixXd::Zero(test_num, train_num);
  
  for(int i = 0; i < test_num; i++) {
    for(int j = 0; j < train_num; j++) {
      double element = 0.0;
      
      for(int k = 0; k < train_num; k++){
        VectorXd dxtest_row = dxtest.row(k);
        VectorXd vec1 = Map_Expball( dxtest(i + train_num, k), dxtest_row);
        VectorXd vec2 = Map_Expball( dxtest(j, k), dxtest_row);
        element += (vec1.transpose() * vec2);
      } 
      
      out(i, j) = element;
    }
  }
  
  return Rcpp::wrap(out / (train_num * train_num));
}


// [[Rcpp::export]]
extern "C" SEXP MetricKtest_rest(const MatrixXd& dxtest, int train_num, int return_num, int test_num) {
  MatrixXd out = MatrixXd::Zero(test_num, return_num);
  
  for(int i = 0; i < test_num; i++) {
    for(int j = 0; j < return_num; j++) {
      double element = 0.0;
      
      for(int k = 0; k < train_num; k++){
        VectorXd dxtest_row = dxtest.row(k);
        VectorXd vec1 = Map_Expball( dxtest(i + train_num, k), dxtest_row.head(train_num));
        VectorXd vec2 = Map_Expball( dxtest(j + train_num, k), dxtest_row.head(train_num));
        element += (vec1.transpose() * vec2);
      } 
      
      out(i, j) = element;
    }
  }
  
  return Rcpp::wrap(out / (train_num * train_num));
}

// [[Rcpp::export]]
extern "C" SEXP MetricKtestSim(const MatrixXd& dx, const MatrixXd& dx1, const VectorXd loc, int train_num, int test_num) {
  
  const int n1 = dx1.rows();
  Rcpp::IntegerVector  lab(n1);
  int m;
  MatrixXd out = MatrixXd::Zero(test_num, train_num);
  
  lab = loc.cast<int>();
  
  for(int i = 0; i < test_num; i++) {
    for(int j = 0; j < train_num; j++) {
      double element = 0.0;
      
      for(int k = 0; k < n1; k++){
        VectorXd dxtest_row = dx1.row(k);
        m = lab[k];
        VectorXd vec1 = Map_Expball( dx(i + train_num, m), dxtest_row);
        VectorXd vec2 = Map_Expball( dx(j, m), dxtest_row);
        element += (vec1.transpose() * vec2);
      } 
      
      out(i, j) = element;
    }
  }
  
  return Rcpp::wrap(out / (n1 * n1));
}

// [[Rcpp::export]]
VectorXd Map_Maxout(const double d1, const VectorXd d2){
   const int n = d2.rows();
   VectorXd d;
   d.setZero(n, 1);
   return (d2.array() - d1).array().max(d.array()) / (d2.array() + 0.0001).array();
}

// [[Rcpp::export]]
extern "C" SEXP Maxout(const MatrixXd& dx){
   const int n = dx.rows();
   MatrixXd out = MatrixXd::Zero(n, n);
   
   for (int i = 0; i < n; i++){
      for (int j = 0; j <= i; j++){
         double element = 0.0;
         
         for (int k = 0; k < n; k++){
          //  VectorXd dx_row = dx.row(k).t();
            VectorXd dx_row = dx.row(k);
            VectorXd vec1 = Map_Maxout( dx(i, k), dx_row );
            VectorXd vec2 = Map_Maxout( dx(j, k), dx_row );
            element += (vec1.transpose() * vec2);
           // element = vec1.array() * vec2.array();
         }
         
         out(i, j) = element;
         out(j, i) = out(i, j);
      }
   }
   
   return Rcpp::wrap(out / (n * n));
}


// [[Rcpp::export]]
extern "C" SEXP Maxout_dxtest(const MatrixXd& dxtest, int train_num, int test_num) {
   MatrixXd out = MatrixXd::Zero(test_num, train_num);
   
   for(int i = 0; i < test_num; i++) {
      for(int j = 0; j < train_num; j++) {
         double element = 0.0;
         
         for(int k = 0; k < train_num; k++){
            VectorXd dxtest_row = dxtest.row(k);
            VectorXd vec1 = Map_Maxout( dxtest(i + train_num, k), dxtest_row );
            VectorXd vec2 = Map_Maxout( dxtest(j, k), dxtest_row );
            element += vec1.transpose() * vec2;
         } 
         
         out(i, j) = element;
      }
   }
   
   return Rcpp::wrap(out / (train_num * train_num));
}

