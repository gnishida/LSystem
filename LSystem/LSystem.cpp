#include "LSystem.h"
#include <QGLWidget>
#include <iostream>

#define M_PI		3.141592653592

LSystem::LSystem() {
	N = 5;
	delta = 25.7;
	rule = "F[+F]F[-F]F";

}

void LSystem::draw() {
	glm::mat4 modelMat;
	drawSegment(modelMat, 0);
}

void LSystem::drawSegment(glm::mat4& modelMat, int level) {
	std::list<glm::mat4> listMat;

	if (N == level) {
		drawQuad(modelMat, 1, 1, 5, QColor(0, 255, 0));
		modelMat = glm::translate(modelMat, glm::vec3(0, 5, 0));
	} else {
		for (int i = 0; i < rule.length(); ++i) {
			if (rule[i] == 'F') {
				drawSegment(modelMat, level + 1);
			} else if (rule[i] == '[') {
				listMat.push_back(modelMat);
			} else if (rule[i] == ']') {
				modelMat = listMat.back();
				listMat.pop_back();
			} else if (rule[i] == '+') {
				modelMat = glm::rotate(modelMat, deg2rad(delta), glm::vec3(0, 0, 1));
			} else if (rule[i] == '-') {
				modelMat = glm::rotate(modelMat, deg2rad(-delta), glm::vec3(0, 0, 1));
			}
		}
	}
}

void LSystem::drawQuad(const glm::mat4& modelMat, float top, float base, float height, const QColor& color) {
	glm::vec4 p1(-base * 0.5, 0, 0, 1);
	glm::vec4 p2(base * 0.5, 0, 0, 1);
	glm::vec4 p3(top * 0.5, height, 0, 1);
	glm::vec4 p4(-top * 0.5, height, 0, 1);

	p1 = modelMat * p1;
	p2 = modelMat * p2;
	p3 = modelMat * p3;
	p4 = modelMat * p4;

	glBegin(GL_QUADS);
	glNormal3f(0, 0, 1);
	glColor3f(color.redF(), color.greenF(), color.blueF());
	glVertex3f(p1.x, p1.y, p1.z);
	glVertex3f(p2.x, p2.y, p2.z);
	glVertex3f(p3.x, p3.y, p3.z);
	glVertex3f(p4.x, p4.y, p4.z);
	glEnd();
}

float LSystem::deg2rad(float deg) {
	return deg * M_PI / 180.0;
}