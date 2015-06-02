#pragma once

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

private:
	int NUM_GRID;
	int NUM_STAT_GRID;

	char axiom;
	map<char, vector<pair<double, string> > > rules;
	string rule;
	std::mt19937 mt;

	cv::Mat_<float> deltas;
	cv::Mat_<int> levels;
	cv::Mat_<float> lengths;

public:
	Stats stats;

public:
	GreedyLSystem();
	void draw();
	void meanInit(int num_grid, int num_stat_grid);
	void randomInit(int num_grid, int num_stat_grid, int seed);
	void generate(int num_grid, int num_stat_grid, const cv::Mat_<double>& target_density, const cv::Mat_<double>& sampleX, const cv::Mat_<double>& sampleY);
	void setParams(int num_grid, int num_stat_grid, const cv::Mat_<double>& mat);
	cv::Mat_<double> getParams() const;
	cv::Mat_<double> getStatistics() const;

private:
	void drawSegment(glm::mat4& modelMat, int level, string rule);
	void drawCylinder(const glm::mat4& modelMat, float top_radius, float base_radius, float height, const glm::vec3& color);
	void drawCircle(const glm::mat4& modelMat, float length, float width, const glm::vec3& color);
	string chooseRule(const vector<pair<double, string> >& rules);
	pair<int, int> XYtoUV(float x, float y);
	float deg2rad(float deg);
	float genRand();
	float genRand(float a, float b);
};

}
