#include "LinearRegressionRegularization.h"

using namespace std;

LinearRegressionRegularization::LinearRegressionRegularization() {
}


void LinearRegressionRegularization::train(const cv::Mat_<double>& X, const cv::Mat_<double>& Y) {
	//W = X.inv(cv::DECOMP_SVD) * Y;
	W = cv::Mat_<double>::zeros(X.cols, Y.cols);

	float lambda = 0.01f;
	float alpha = 0.1f;
	int maxIter = 100;

	for (int iter = 0; iter < maxIter; ++iter) {
		cv::Mat_<double> dW = X.t() * (X * W - Y) / X.rows + lambda * W;

		W -= dW * alpha;
	}
}

cv::Mat LinearRegressionRegularization::predict(const cv::Mat_<double>& x) {
	return x * W;
}

double LinearRegressionRegularization::conditionNumber() {
	cv::Mat_<double> w1, w2, u, vt;
	cv::SVD::compute(W, w1, u, vt);

	cv::SVD::compute(W.inv(cv::DECOMP_SVD), w2, u, vt);

	return w1(0,0) * w2(0, 0);
}

