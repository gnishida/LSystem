#pragma once

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/string_cast.hpp>
#include <vector>
#include <QColor>

using namespace std;

class LSystem {
public:
	int N;
	double delta;
	string rule;

public:
	LSystem();
	void draw();

private:
	void drawSegment(glm::mat4& modelMat, int level);
	void drawQuad(const glm::mat4& modelMat, float top, float base, float height, const QColor& color);
	float deg2rad(float deg);
};

