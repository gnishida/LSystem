#include "ParametricLSystem.h"
#include <QGLWidget>
#include <iostream>
#include <time.h>
#include "Eval.h"
#include <algorithm>

namespace parametriclsystem {

const double M_PI = 3.141592653592;

Rule::Rule(const string& left_hand, const string& condition, const string& right_hand) {
	int index1 = left_hand.find('(');
	int index2 = left_hand.find(')');

	string arg = left_hand.substr(index1 + 1, index2 - index1 - 1);
	vector<string> arg_list = split(arg, ',');
	for (int i = 0; i < arg_list.size(); ++i) {
		variables.push_back(arg_list[i]);
	}

	this->left_hand = left_hand.substr(0, index1);
	this->condition = condition;
	this->right_hand = right_hand;
}

bool Rule::isTrue(const vector<double>& values) {
	string cond = this->condition;
	for (int i = 0; i < variables.size(); ++i) {
		replaceAll(cond, variables[i], to_string((long double)values[i]));
	}
	return eval::compare(cond);
}

string Rule::derive(const string arg) {
	vector<string> arg_list = split(arg, ',');
	vector<double> values;
	for (int i = 0; i < arg_list.size(); ++i) {
		values.push_back(stof(arg_list[i]));
	}

	return derive(values);
}

string Rule::derive(const vector<double>& values) {
	string r = right_hand;
	for (int i = 0; i < variables.size(); ++i) {
		replaceAll(r, variables[i], to_string((long double)values[i]));
	}

	string r2;
	char c;
	for (int i = 0; i < r.length(); ) {
		if (r[i] == '(') {
			r2 += '(';
			int index = r.find(')', i + 1);
			string arg = r.substr(i + 1, index - i - 1);
			vector<string> vec_arg = split(arg, ',');
			for (int j = 0; j < vec_arg.size(); ++j) {
				if (j > 0) {
					r2 += ",";
				}
				r2 += to_string((long double)eval::calculate(vec_arg[j]));
			}
			r2 += ')';
			i = index + 1;
		} else {
			r2 += r[i++];
		}
	}

	return r2;
}

ParametricLSystem::ParametricLSystem() {
	axiom = "A(20,5)";
	N = 10;
	rules['A'] = Rule("A(l,w)", "*", "!(w)F(l)[&(45)B(l*0.6,w*0.707)]/(137.5)A(l*0.9,w*0.707)");
	rules['B'] = Rule("B(l,w)", "*", "!(w)F(l)[-(45)$C(l*0.6,w*0.707)]C(l*0.9,w*0.707)");
	rules['C'] = Rule("C(l,w)", "*", "!(w)F(l)[+(45)$B(l*0.6,w*0.707)]B(l*0.9,w*0.707)");

	/*
	axiom = "A(20,5)";
	N = 10;
	rules['A'] = Rule("A(l,w)", "*", "!(w)F(l)[&(5)B(l*0.7,w*0.707)]/(137.5)A(l*0.9,w*0.707)");
	rules['B'] = Rule("B(l,w)", "*", "!(w)F(l)[-(65)$C(l*0.7,w*0.707)]C(l*0.9,w*0.707)");
	rules['C'] = Rule("C(l,w)", "*", "!(w)F(l)[+(65)$B(l*0.7,w*0.707)]B(l*0.9,w*0.707)");
	*/

	srand(time(NULL));
	rule = derive();
	//cout << rule << endl;
}

string ParametricLSystem::derive() {
	string result = axiom;
	for (int n = 0; n < N; ++n) {
		string next;
		for (int i = 0; i < result.length(); ) {
			if (result[i] >= 'A' && result[i] <= 'Z' && rules.find(result[i]) != rules.end()) {
				int index1 = result.find('(', i + 1);
				int index2 = result.find(')', i + 1);

				string arg = result.substr(index1 + 1, index2 - index1 - 1);

				next += rules[result[i]].derive(arg);

				i = index2 + 1;
			} else {
				next += result[i++];
			}
		}

		result = next;
	}

	return result;
}

void ParametricLSystem::draw() {
	drawSegment(rule);
}

void ParametricLSystem::drawSegment(string rule) {
	std::list<State> listState;

	State state;
	for (int i = 0; i < rule.length(); ) {
		if (rule[i] == '[') {
			listState.push_back(state);
			i++;
		} else if (rule[i] == ']') {
			state = listState.back();
			listState.pop_back();
			i++;
		} else if (rule[i] == '+') {
			state.modelMat = glm::rotate(state.modelMat, deg2rad(extractArgument(rule, i + 1, i)), glm::vec3(0, 0, 1));
		} else if (rule[i] == '-') {
			state.modelMat = glm::rotate(state.modelMat, deg2rad(-extractArgument(rule, i + 1, i)), glm::vec3(0, 0, 1));
		} else if (rule[i] == '\\') {
			state.modelMat = glm::rotate(state.modelMat, deg2rad(extractArgument(rule, i + 1, i)), glm::vec3(0, 1, 0));
		} else if (rule[i] == '/') {
			state.modelMat = glm::rotate(state.modelMat, deg2rad(-extractArgument(rule, i + 1, i)), glm::vec3(0, 1, 0));
		} else if (rule[i] == '&') {
			state.modelMat = glm::rotate(state.modelMat, deg2rad(extractArgument(rule, i + 1, i)), glm::vec3(1, 0, 0));
		} else if (rule[i] == '^') {
			state.modelMat = glm::rotate(state.modelMat, deg2rad(-extractArgument(rule, i + 1, i)), glm::vec3(1, 0, 0));
		} else if (rule[i] == '|') {
			state.modelMat = glm::rotate(state.modelMat, deg2rad(180), glm::vec3(0, 0, 1));
		} else if (rule[i] == '!') {
			state.diameter = extractArgument(rule, i + 1, i);
		} else if (rule[i] == '\'') {
			//state.color.g = min(1, state.color.g + 0.2);
		} else if (rule[i] == '$') {
			glm::vec4 y(0, 1, 0, 0);
			y = state.modelMat * y;
			glm::vec3 v(0, 1, 0);
			glm::vec3 x = glm::normalize(glm::cross(glm::vec3(y), v));
			glm::vec3 z = glm::cross(x, glm::vec3(y));

			state.modelMat[0] = glm::vec4(x, 0);
			state.modelMat[1] = y;
			state.modelMat[2] = glm::vec4(z, 0);

			i++;
		} else if (rule[i] == 'f') {
			drawCircle(state.modelMat, 20.0, 4.0, state.color);
		} else if (rule[i] == 'F') {
			double length = extractArgument(rule, i + 1, i);
			drawCylinder(state.modelMat, state.diameter * 0.5 * 0.707, state.diameter * 0.5, length, state.color);
			state.modelMat = glm::translate(state.modelMat, glm::vec3(0, length, 0));
		} else {
			i++;
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
void ParametricLSystem::drawCylinder(const glm::mat4& modelMat, float top_radius, float base_radius, float height, const glm::vec3& color) {
	int slices = 22;

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
void ParametricLSystem::drawCircle(const glm::mat4& modelMat, float length, float width, const glm::vec3& color) {
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
string ParametricLSystem::chooseRule(const vector<pair<double, string> >& rules) {
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

void replaceAll(std::string& str, const std::string& from, const std::string& to) {
    std::string::size_type pos = 0;
    while(pos = str.find(from, pos), pos != std::string::npos) {
        str.replace(pos, from.length(), to);
        pos += to.length();
    }
}

float deg2rad(float deg) {
	return deg * M_PI / 180.0;
}

std::vector<std::string> split(const string& str, char delim) {
	std::vector<std::string> res;
	size_t current = 0, found;
	while ((found = str.find_first_of(delim, current)) != string::npos) {
		res.push_back(string(str, current, found - current));
		current = found + 1;
	}
	res.push_back(string(str, current, str.size() - current));
	return res;
}

double extractArgument(const string& str, int startIndex, int& nextIndex) {
	int index1 = str.find('(', startIndex);
	int index2 = str.find(')', startIndex);

	nextIndex = index2 + 1;

	return stof(str.substr(index1 + 1, index2 - index1 - 1));
}

}
