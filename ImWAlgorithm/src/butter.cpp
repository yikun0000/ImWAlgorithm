#include <RcppArmadillo.h>
#include <cmath>

// 使用arma命名空间
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]
template<typename T>
T fftshift(const T& mat) {
  size_t N1 = mat.n_rows;
  if (mat.n_rows%2!=0){
    N1 = mat.n_rows-1;
  }
  else{
    N1 = mat.n_rows;
  }
  size_t N2 = mat.n_cols;
  if (mat.n_cols%2!=0){
    N2 = mat.n_cols-1;
  }
  else{
    N2 = mat.n_cols;
  }
  T result(mat);
  // 重新排列四个象限
  result.submat(0, 0, N1 / 2 - 1, N2 / 2 - 1) = mat.submat(N1 / 2, N2 / 2, N1 - 1, N2 - 1);
  result.submat(N1 / 2, N2 / 2, N1 - 1, N2 - 1) = mat.submat(0, 0, N1 / 2 - 1, N2 / 2 - 1);
  result.submat(0, N2 / 2, N1 / 2 - 1, N2 - 1) = mat.submat(N1 / 2, 0, N1 - 1, N2 / 2 - 1);
  result.submat(N1 / 2, 0, N1 - 1, N2 / 2 - 1) = mat.submat(0, N2 / 2, N1 / 2 - 1, N2 - 1);
  return result;
}
template<typename T>
T ifftshift(const T& mat) {
  size_t N1 = mat.n_rows;
  if (mat.n_rows%2!=0){
    N1 = mat.n_rows-1;
  }
  else{
    N1 = mat.n_rows;
  }
  size_t N2 = mat.n_cols;
  if (mat.n_cols%2!=0){
    N2 = mat.n_cols-1;
  }
  else{
    N2 = mat.n_cols;
  }
  T result(mat);
  // 重新排列四个象限
  result.submat(0, 0, N1 / 2 - 1, N2 / 2 - 1) = mat.submat((N1 + 1) / 2, (N2 + 1) / 2, N1 - 1, N2 - 1);
  result.submat((N1 + 1) / 2, (N2 + 1) / 2, N1 - 1, N2 - 1) = mat.submat(0, 0, N1 / 2 - 1, N2 / 2 - 1);
  result.submat(0, (N2 + 1) / 2, N1 / 2 - 1, N2 - 1) = mat.submat((N1 + 1) / 2, 0, N1 - 1, N2 / 2 - 1);
  result.submat((N1 + 1) / 2, 0, N1 - 1, N2 / 2 - 1) = mat.submat(0, (N2 + 1) / 2, N1 / 2 - 1, N2 - 1);
  return result;
}
// [[Rcpp::export]]
arma::mat butterworthLowPassFilter(arma::mat image, double d0, int n) {
  arma::cx_mat s = arma::fft2(image); // 二维傅里叶变换
  s = fftshift(s); // 频谱平移
  int N1 = s.n_rows;
  int N2 = s.n_cols;
  int n1 = N1 / 2;
  int n2 = N2 / 2;
  for (int i = 0; i < N1; ++i) {
    for (int j = 0; j < N2; ++j) {
      double distance = std::sqrt(std::pow(i - n1, 2) + std::pow(j - n2, 2));
      double h = 0;
      if (distance != 0) {
        h = 1 / (1 + std::pow((distance / d0), 2 * n));
      }
      s(i, j) *= h; // 应用滤波器
    }
  }
  s = ifftshift(s); // 频谱反平移
  arma::mat result = arma::real(arma::ifft2(s)); // 二维傅里叶反变换
  return result;
}
