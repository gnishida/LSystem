#pragma once

#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <opencv/cv.h>
#include <opencv/highgui.h>

namespace ml {

std::vector<std::string> splitDataset(const std::string &str, char sep);
void splitDataset(const cv::Mat_<double>& data, float ratio1, cv::Mat_<double>& data1, cv::Mat_<double>& data2);
void splitDataset(const cv::Mat_<double>& data, float ratio1, float ratio2, cv::Mat_<double>& data1, cv::Mat_<double>& data2, cv::Mat_<double>& data3);
bool loadDataset(char* filename, cv::Mat_<double>& X, cv::Mat_<double>& Y);
void saveDataset(char* filename, const cv::Mat_<double>& X, const cv::Mat_<double>& Y);
void normalizeDataset(cv::Mat_<double> mat, cv::Mat_<double>& normalized_mat, cv::Mat_<double>& mu, cv::Mat_<double>& abs_max);
void addBias(cv::Mat_<double>& data);
cv::Mat mat_square(const cv::Mat& m);


}