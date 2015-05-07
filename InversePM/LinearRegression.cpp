#include "LinearRegression.h"

using namespace std;

LinearRegression::LinearRegression() {
}


void LinearRegression::train(const cv::Mat_<double>& X, const cv::Mat_<double>& Y) {
	W = X.inv(cv::DECOMP_SVD) * Y;
}

cv::Mat LinearRegression::predict(const cv::Mat_<double>& x) {
	return x * W;
}

double LinearRegression::conditionNumber() {
	cv::Mat_<double> w1, w2, u, vt;
	cv::SVD::compute(W, w1, u, vt);

	cv::SVD::compute(W.inv(cv::DECOMP_SVD), w2, u, vt);

	return w1(0,0) * w2(0, 0);
}

