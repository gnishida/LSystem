#include "LSystem.h"
#include <QGLWidget>
#include <iostream>

#define M_PI		3.141592653592

LSystem::LSystem() {
	/*
	N = 5;
	delta = 25.7;
	axiom = 'F';
	rules['F'] = "F[+F]F[-F]F";
	*/

	/*
	N = 5;
	delta = 20;
	axiom = 'F';
	rules['F'] = "F[+F]F[-F][F]";
	*/

	/*
	N = 4;
	delta = 22.5;
	axiom = 'F';
	rules['F'] = "FF-[-F+F+F]+[+F-F-F]";
	*/

	/*
	N = 7;
	delta = 20;
	axiom = 'X';
	rules['X'] = "F[+X]F[-X]+X";
	rules['F'] = "FF";
	*/

	/*
	N = 7;
	delta = 25.7;
	axiom = 'X';
	rules['X'] = "F[+X][-X]FX";
	rules['F'] = "FF";
	*/

	/*
	N = 5;
	delta = 22.5;
	axiom = 'X';
	rules['X'] = "F-[[X]+X]+F[+FX]-X";
	rules['F'] = "FF";
	*/

	N = 7;
	delta = 22.5;
	axiom = 'A';
	rules['A'] = "[&F!A]/////'[&F!A]///////'[&F!A]";
	rules['F'] = "S/////F";
	rules['S'] = "F";
}

void LSystem::draw() {
	glm::mat4 modelMat;
	drawSegment(modelMat, 0, axiom, 1.0);
}

void LSystem::drawSegment(glm::mat4& modelMat, int level, char left_hand, double scale) {
	std::list<glm::mat4> listMat;

	string rule = rules.find(left_hand)->second;

	if (N == level) {
		if (left_hand == 'F') {
			drawCylinder(modelMat, scale, scale, 5.0 * scale, QColor(0, 255, 0));
			modelMat = glm::translate(modelMat, glm::vec3(0, 5.0 * scale, 0));
		}
	} else {
		for (int i = 0; i < rule.length(); ++i) {
			if (rule[i] == '[') {
				listMat.push_back(modelMat);
			} else if (rule[i] == ']') {
				modelMat = listMat.back();
				listMat.pop_back();
			} else if (rule[i] == '+') {
				modelMat = glm::rotate(modelMat, deg2rad(delta), glm::vec3(0, 0, 1));
			} else if (rule[i] == '-') {
				modelMat = glm::rotate(modelMat, deg2rad(-delta), glm::vec3(0, 0, 1));
			} else if (rule[i] == '\\') {
				modelMat = glm::rotate(modelMat, deg2rad(delta), glm::vec3(0, 1, 0));
			} else if (rule[i] == '/') {
				modelMat = glm::rotate(modelMat, deg2rad(-delta), glm::vec3(0, 1, 0));
			} else if (rule[i] == '&') {
				modelMat = glm::rotate(modelMat, deg2rad(-delta), glm::vec3(1, 0, 0));
			} else if (rule[i] == '^') {
				modelMat = glm::rotate(modelMat, deg2rad(delta), glm::vec3(1, 0, 0));
			} else if (rule[i] == '!') {
				scale *= 0.9;
			} else if (rule[i] == '\'') {

			} else {
				drawSegment(modelMat, level + 1, rule[i], scale);
			}
		}
	}
}

void LSystem::drawCylinder(const glm::mat4& modelMat, float top_radius, float base_radius, float height, const QColor& color) {
	int slices = 12;

	glBegin(GL_QUADS);
	for (int i = 0; i < slices; ++i) {
		float theta1 = (float)i / slices * 2 * M_PI;
		float theta2 = (float)(i + 1) / slices * 2 * M_PI;

		glm::vec4 p1(cosf(theta1) * base_radius, 0, -sinf(theta1) * base_radius, 1);
		glm::vec4 p2(cosf(theta2) * base_radius, 0, -sinf(theta2) * base_radius, 1);
		glm::vec4 p3(cosf(theta2) * base_radius, height, -sinf(theta2) * base_radius, 1);
		glm::vec4 p4(cosf(theta1) * base_radius, height, -sinf(theta1) * base_radius, 1);

		glm::vec4 n1(cosf(theta1) * base_radius, 0, -sinf(theta1) * base_radius, 1);
		glm::vec4 n2(cosf(theta2) * base_radius, 0, -sinf(theta2) * base_radius, 1);
		glm::vec4 n3(cosf(theta2) * base_radius, 0, -sinf(theta2) * base_radius, 1);
		glm::vec4 n4(cosf(theta1) * base_radius, 0, -sinf(theta1) * base_radius, 1);
		
		p1 = modelMat * p1;
		p2 = modelMat * p2;
		p3 = modelMat * p3;
		p4 = modelMat * p4;

		glColor3f(color.redF(), color.greenF(), color.blueF());
		glNormal3f(n1.x, n1.y, n1.z);
		glVertex3f(p1.x, p1.y, p1.z);
		glNormal3f(n2.x, n2.y, n2.z);
		glVertex3f(p2.x, p2.y, p2.z);
		glNormal3f(n3.x, n3.y, n3.z);
		glVertex3f(p3.x, p3.y, p3.z);
		glNormal3f(n4.x, n4.y, n4.z);
		glVertex3f(p4.x, p4.y, p4.z);
	}
	glEnd();


}

float LSystem::deg2rad(float deg) {
	return deg * M_PI / 180.0;
}