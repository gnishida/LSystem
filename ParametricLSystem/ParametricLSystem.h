#pragma once

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/string_cast.hpp>
#include <vector>
#include <map>

using namespace std;

namespace parametriclsystem {

class State {
public:
	double radius;
	glm::mat4 modelMat;
	glm::vec3 color;

public:
	State() : radius(0.3), color(0, 0.7, 0) {};
};

class Rule {
public:
	string left_hand;
	vector<string> variables;
	string condition;
	string right_hand;

public:
	Rule() {}
	Rule(const string& left_hand, const string& condition, const string& right_hand);

	bool isTrue(const vector<double>& values);
	string derive(const string arg);
	string derive(const vector<double>& values);
};

class ParametricLSystem {
public:
	int N;
	double delta;
	string axiom;
	map<char, Rule> rules;
	string rule;

public:
	ParametricLSystem();
	string derive();
	void draw();

private:
	void drawSegment(string rule);
	void drawCylinder(const glm::mat4& modelMat, float top_radius, float base_radius, float height, const glm::vec3& color);
	void drawCircle(const glm::mat4& modelMat, float length, float width, const glm::vec3& color);
	string chooseRule(const vector<pair<double, string> >& rules);
};

void replaceAll(string& str, const string& from, const string& to);
float deg2rad(float deg);
vector<string> split(const string& str, char delim);
double extractArgument(const string& str, int startIndex, int& nextIndex);

}
