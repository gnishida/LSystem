﻿#pragma once

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/string_cast.hpp>
#include <vector>
#include <map>
#include <random>

using namespace std;

#include <opencv/cv.h>
#include <opencv/highgui.h>

namespace greedylsystem {

// 統計情報（high-level indicators）
class Stats {
public:
	cv::Mat_<int> coverage;
	cv::Mat_<double> density;
};

class GreedyLSystem {
public:
	static double GRID_SIZE;
	static int MIN_DELTA;
	static int MAX_DELTA;
	static int MIN_LEVEL;
	static int MAX_LEVEL;
	static int MIN_LENGTH;
	static int MAX_LENGTH;

private:
	int NUM_GRID;
	int NUM_STAT_GRID;

	std::mt19937 mt;

	cv::Mat_<float> deltas;
	cv::Mat_<int> levels;
	cv::Mat_<float> lengths;

public:
	Stats stats;

public:
	GreedyLSystem();
	void draw(bool clearUnusedParams);
	void meanInit(int num_grid, int num_stat_grid);
	void randomInit(int num_grid, int num_stat_grid, int seed);
	void nearestSample(int num_grid, int num_stat_grid, const cv::Mat_<double>& target_density, const cv::Mat_<double>& sampleX, const cv::Mat_<double>& sampleY);
	void inverse(const cv::Mat_<double>& target_density);
	void setParams(int num_grid, int num_stat_grid, const cv::Mat_<double>& mat);
	cv::Mat_<double> getParams() const;
	cv::Mat_<double> getStatistics() const;

private:
	void drawSegment(glm::mat4 modelMat, int level);
	void drawCylinder(const glm::mat4& modelMat, float top_radius, float base_radius, float height, const glm::vec3& color);
	void drawCircle(const glm::mat4& modelMat, float length, float width, const glm::vec3& color);
	pair<int, int> XYtoUV(float x, float y);
	float deg2rad(float deg);
	float genRand();
	float genRand(float a, float b);
};

}
