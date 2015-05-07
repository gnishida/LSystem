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

namespace lsystem {

// 統計情報（high-level indicators）
class Stats {
public:
	cv::Mat_<float> density;
};

class LSystem {
public:
	static double GRID_SIZE;
	static int NUM_GRID;
	static double CELL_SIZE;
	static int NUM_STATS_GRID;

public:
	int N;
	double delta;
	char axiom;
	map<char, vector<pair<double, string> > > rules;
	string rule;
	std::mt19937 mt;
	double segment_length;

	cv::Mat_<float> deltas;
	cv::Mat_<int> levels;
	cv::Mat_<float> lengths;

	Stats stats;

public:
	LSystem();
	void draw();
	void randomInit(int seed);
	void setParams(const cv::Mat_<float>& mat);
	vector<float> getParams();
	vector<float> getStatistics();


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
