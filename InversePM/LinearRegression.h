#pragma once

#include <opencv/cv.h>
#include <opencv/highgui.h>

class LinearRegression {
public:
	cv::Mat_<double> W;

public:
	LinearRegression();

	void train(const cv::Mat_<double>& X, const cv::Mat_<double>& Y);
	cv::Mat predict(const cv::Mat_<double>& x);
	double conditionNumber();
};

