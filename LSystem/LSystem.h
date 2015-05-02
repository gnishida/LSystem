#pragma once

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/string_cast.hpp>
#include <vector>
#include <map>
#include <QColor>

using namespace std;

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
	void drawSegment(glm::mat4& modelMat, int level, char left_hand, double scale);
	void drawCylinder(const glm::mat4& modelMat, float top_radius, float base_radius, float height, const QColor& color);
	float deg2rad(float deg);
};

