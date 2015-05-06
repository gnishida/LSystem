#pragma once

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/string_cast.hpp>
#include <vector>
#include <map>
#include <random>

using namespace std;

namespace lsystem {

class State {
public:
	double radius;
	glm::mat4 modelMat;
	glm::vec3 color;

public:
	State() : radius(1.0), color(0, 0.7, 0) {};
};

class LSystem {
public:
	int N;
	double delta;
	char axiom;
	map<char, vector<pair<double, string> > > rules;
	string rule;
	std::mt19937 mt;

public:
	LSystem();
	string derive();
	void draw();

private:
	void drawSegment(string rule);
	void drawCylinder(const glm::mat4& modelMat, float top_radius, float base_radius, float height, const glm::vec3& color);
	void drawCircle(const glm::mat4& modelMat, float length, float width, const glm::vec3& color);
	string chooseRule(const vector<pair<double, string> >& rules);
	float deg2rad(float deg);
};

}
