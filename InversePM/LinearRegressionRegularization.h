#pragma once

#include <opencv/cv.h>
#include <opencv/highgui.h>

class LinearRegressionRegularization {
public:
	cv::Mat_<double> W;

public:
	LinearRegressionRegularization();

	void train(const cv::Mat_<double>& X, const cv::Mat_<double>& Y);
	cv::Mat predict(const cv::Mat_<double>& x);
	double conditionNumber();
};

