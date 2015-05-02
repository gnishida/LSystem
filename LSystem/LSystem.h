#pragma once

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/string_cast.hpp>
#include <vector>
#include <map>
#include <QColor>

using namespace std;

namespace lsystem {

class State {
public:
	double radius;
	glm::mat4 modelMat;
	glm::vec3 color;

public:
	State() : radius(1.0), color(0, 1, 0) {};
};

class LSystem {
public:
	int N;
	double delta;
	char axiom;
	map<char, string> rules;

public:
	LSystem();
	void draw();

private:
	void drawSegment(State& state, int level, char left_hand);
	void drawCylinder(const glm::mat4& modelMat, float top_radius, float base_radius, float height, const glm::vec3& color);
	void drawCircle(const glm::mat4& modelMat, float length, float width, const glm::vec3& color);
	float deg2rad(float deg);
};

}
