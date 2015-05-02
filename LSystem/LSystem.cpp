#include "LSystem.h"
#include <QGLWidget>
#include <iostream>
#include <time.h>

namespace lsystem {

const double M_PI = 3.141592653592;

LSystem::LSystem() {
	/*
	N = 5;
	delta = 25.7;
	axiom = 'F';
	rules['F'].push_back(pair<double, string>(1.0, "F[+F]F[-F]F"));
	*/

	/*
	N = 5;
	delta = 20;
	axiom = 'F';
	rules['F'].push_back(pair<double, string>(1.0, "F[+F]F[-F][F]"));
	*/

	/*
	N = 4;
	delta = 22.5;
	axiom = 'F';
	rules['F'].push_back(pair<double, string>(1.0, "FF-[-F+F+F]+[+F-F-F]"));
	*/

	/*
	N = 7;
	delta = 20;
	axiom = 'X';
	rules['X'].push_back(pair<double, string>(1.0, "F[+X]F[-X]+X"));
	rules['F'].push_back(pair<double, string>(1.0, "FF"));
	*/

	/*
	N = 7;
	delta = 25.7;
	axiom = 'X';
	rules['X'].push_back(pair<double, string>(1.0, "F[+X][-X]FX"));
	rules['F'].push_back(pair<double, string>(1.0, "FF"));
	*/

	/*
	N = 5;
	delta = 22.5;
	axiom = 'X';
	rules['X'].push_back(pair<double, string>(1.0, "F-[[X]+X]+F[+FX]-X"));
	rules['F'].push_back(pair<double, string>(1.0, "FF"));
	*/

	/*
	N = 7;
	delta = 22.5;
	axiom = 'A';
	rules['A'].push_back(pair<double, string>(1.0, "[&FL!A]/////'[&FL!A]///////'[&FL!A]"));
	rules['F'].push_back(pair<double, string>(1.0, "S/////F"));
	rules['S'].push_back(pair<double, string>(1.0, "FL"));
	rules['L'].push_back(pair<double, string>(1.0, "['''^^f]"));
	*/

	N = 5;
	delta = 25.7;
	axiom = 'F';
	rules['F'].push_back(pair<double, string>(0.33, "F[+F]F[-F]F"));
	rules['F'].push_back(pair<double, string>(0.33, "F[+F]F"));
	rules['F'].push_back(pair<double, string>(0.34, "F[-F]F"));

	random_seed = time(NULL);
}

void LSystem::draw() {
	srand(random_seed);

	State state;
	drawSegment(state, 0, axiom);
}

void LSystem::drawSegment(State& state, int level, char left_hand) {
	std::list<State> listState;

	if (rules.find(left_hand) == rules.end()) return;

	string rule = chooseRule(rules[left_hand]);

	for (int i = 0; i < rule.length(); ++i) {
		if (rule[i] == '[') {
			listState.push_back(state);
		} else if (rule[i] == ']') {
			state = listState.back();
			listState.pop_back();
		} else if (rule[i] == '+') {
			state.modelMat = glm::rotate(state.modelMat, deg2rad(delta), glm::vec3(0, 0, 1));
		} else if (rule[i] == '-') {
			state.modelMat = glm::rotate(state.modelMat, deg2rad(-delta), glm::vec3(0, 0, 1));
		} else if (rule[i] == '\\') {
			state.modelMat = glm::rotate(state.modelMat, deg2rad(delta), glm::vec3(0, 1, 0));
		} else if (rule[i] == '/') {
			state.modelMat = glm::rotate(state.modelMat, deg2rad(-delta), glm::vec3(0, 1, 0));
		} else if (rule[i] == '&') {
			state.modelMat = glm::rotate(state.modelMat, deg2rad(delta), glm::vec3(1, 0, 0));
		} else if (rule[i] == '^') {
			state.modelMat = glm::rotate(state.modelMat, deg2rad(-delta), glm::vec3(1, 0, 0));
		} else if (rule[i] == '|') {
			state.modelMat = glm::rotate(state.modelMat, deg2rad(180), glm::vec3(0, 0, 1));
		} else if (rule[i] == '!') {
			state.radius *= 0.7;
		} else if (rule[i] == '\'') {
			state.color.g = min(1, state.color.g + 0.2);
		} else if (rule[i] == 'f') {
			drawCircle(state.modelMat, 20.0, 4.0, state.color);
		} else {
			if (level < N && rules.find(rule[i]) != rules.end()) {
				drawSegment(state, level + 1, rule[i]);
			} else {
				if (rule[i] == 'F') {
					drawCylinder(state.modelMat, state.radius, state.radius, 5.0, state.color);
					state.modelMat = glm::translate(state.modelMat, glm::vec3(0, 5.0, 0));
				}
			}
		}
	}
}

/**
 * XZ面を底面とし、Y軸に伸びる円筒形を描画する。
 *
 * @modelMat		モデル行列
 * @top_radius		円筒形の上面の半径
 * @base_radius		円筒形の底面の半径
 * @height			円筒形の高さ
 * @color			色
 */
void LSystem::drawCylinder(const glm::mat4& modelMat, float top_radius, float base_radius, float height, const glm::vec3& color) {
	int slices = 12;

	glBegin(GL_TRIANGLES);
	glColor3f(color.r, color.g, color.b);
	for (int i = 0; i < slices; ++i) {
		float theta1 = (float)i / slices * 2 * M_PI;
		float theta2 = (float)(i + 1) / slices * 2 * M_PI;

		glm::vec4 p1(cosf(theta1) * base_radius, 0, sinf(theta1) * base_radius, 1);
		glm::vec4 p2(cosf(theta2) * base_radius, 0, sinf(theta2) * base_radius, 1);
		glm::vec4 p3(cosf(theta2) * top_radius, height, sinf(theta2) * top_radius, 1);
		glm::vec4 p4(cosf(theta1) * top_radius, height, sinf(theta1) * top_radius, 1);

		glm::vec4 n1(cosf(theta1) * base_radius, 0, sinf(theta1) * base_radius, 0);
		glm::vec4 n2(cosf(theta2) * base_radius, 0, sinf(theta2) * base_radius, 0);
		glm::vec4 n3(cosf(theta2) * top_radius, 0, sinf(theta2) * top_radius, 0);
		glm::vec4 n4(cosf(theta1) * top_radius, 0, sinf(theta1) * top_radius, 0);
		
		p1 = modelMat * p1;
		p2 = modelMat * p2;
		p3 = modelMat * p3;
		p4 = modelMat * p4;

		n1 = modelMat * n1;
		n2 = modelMat * n2;
		n3 = modelMat * n3;
		n4 = modelMat * n4;

		glNormal3f(n1.x, n1.y, n1.z);
		glVertex3f(p1.x, p1.y, p1.z);
		glNormal3f(n2.x, n2.y, n2.z);
		glVertex3f(p2.x, p2.y, p2.z);
		glNormal3f(n3.x, n3.y, n3.z);
		glVertex3f(p3.x, p3.y, p3.z);

		glNormal3f(n1.x, n1.y, n1.z);
		glVertex3f(p1.x, p1.y, p1.z);
		glNormal3f(n3.x, n3.y, n3.z);
		glVertex3f(p3.x, p3.y, p3.z);
		glNormal3f(n4.x, n4.y, n4.z);
		glVertex3f(p4.x, p4.y, p4.z);
	}
	glEnd();
}

/**
 * XY平面上に、原点を中心とする楕円を描画する。
 *
 * @modelMat		モデル行列
 * @length			楕円のY軸方向の直径
 * @width			楕円のX軸方向の直径
 * @color			色
 */
void LSystem::drawCircle(const glm::mat4& modelMat, float length, float width, const glm::vec3& color) {
	int slices = 12;

	glm::mat4 mat = glm::translate(modelMat, glm::vec3(0, length * 0.5, 0));

	glm::vec4 p0(0, 0, 0, 1);
	p0 = mat * p0;
	glm::vec4 n(0, 0, 1, 0);
	n = mat * n;

	glBegin(GL_TRIANGLES);
	glNormal3f(n.x, n.y, n.z);
	glColor3f(color.r, color.g, color.b);
	for (int i = 0; i < slices; ++i) {
		float theta1 = (float)i / slices * 2 * M_PI;
		float theta2 = (float)(i + 1) / slices * 2 * M_PI;

		glm::vec4 p1(cosf(theta1) * width * 0.5, sinf(theta1) * length * 0.5, 0, 1);
		glm::vec4 p2(cosf(theta2) * width * 0.5, sinf(theta2) * length * 0.5, 0, 1);

		p1 = mat * p1;
		p2 = mat * p2;

		glVertex3f(p0.x, p0.y, p0.z);
		glVertex3f(p1.x, p1.y, p1.z);
		glVertex3f(p2.x, p2.y, p2.z);
	}
	glEnd();
}

/**
 * ルールリストから、確率に基づいて１つのルールを選択する。
 * リストの各要素は、<確率、ルール>のペアとなっている。
 *
 * @param rules		ルールリスト
 * @reutnr			選択されたルール
 */
string LSystem::chooseRule(const vector<pair<double, string> >& rules) {
	vector<double> cdf;
	{
		double total = 0.0;
		for (int i = 0; i < rules.size(); ++i) {
			total += rules[i].first;
			cdf.push_back(total);
		}
	}

	double rnd = (double)rand() / (RAND_MAX + 1) * cdf.back();
	for (int i = 0; i < cdf.size(); ++i) {
		if (rnd <= cdf[i]) {
			return rules[i].second;
		}
	}

	return rules.back().second;
}

float LSystem::deg2rad(float deg) {
	return deg * M_PI / 180.0;
}

}
