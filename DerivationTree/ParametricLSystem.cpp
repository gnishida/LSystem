#include "ParametricLSystem.h"
#include <iostream>
#include <time.h>
#include "Eval.h"
#include <algorithm>
#include "MLUtils.h"

namespace parametriclsystem {

const double M_PI = 3.141592653592;

Literal::Literal(char c) {
	this->c = c;
	this->param_value = 0.0;
	this->param_defined = false;
}

Literal::Literal(char c, double param_value) {
	this->c = c;
	this->param_value = param_value;
	this->param_defined = true;
}

String::String(string str) {
	for (int i = 0; i < str.length(); ++i) {
		this->str.push_back(Literal(str[i]));
	}
}

void String::operator+=(const string& str) {
	for (int i = 0; i < str.length(); ++i) {
		this->str.push_back(Literal(str[i]));
	}
}

ostream& operator<<(ostream& os, const String& str) {
	os << setprecision(1);
	for (int i = 0; i < str.length(); ++i) {
		os << str[i].c;
		if (str[i].param_defined) {
			os << "(" << fixed << str[i].param_value << ")";
		}
	}

    return os;
}


Rule::Rule(const string& right_hand) {
	this->right_hand = right_hand;
}

TreeNode::TreeNode(double val, TreeNode* parent) {
	this->val = val;
	this->parent = parent;
}

TreeNode* TreeNode::getChild(double value) {
	if (children.find(value) == children.end()) {
		children[value] = new TreeNode(value, this);
	}

	return children[value];
}

ParametricLSystem::ParametricLSystem() {
	axiom = "X";
	rules['X'].push_back(Rule("F"));
	rules['X'].push_back(Rule("F[-X][+X]"));
}

/**
 * 適当な乱数シードに基づいて、ランダムにgenerateする。
 */
String ParametricLSystem::derive(int random_seed) {
	srand(random_seed);

	String result(axiom);
	TreeNode* node = root;

	for (int iter = 0; iter < 20; ++iter) {
		String next;
		bool found = false;
		for (int i = 0; i < result.length(); ++i) {
			if (found) {
				next += result[i];
			} else if (rules.find(result[i].c) != rules.end()) {
				int index = chooseRule(result[i].c);
				next += rules[result[i].c][index].right_hand;

				node = root->getChild(index);
				found = true;
			} else if (result[i].c == 'F' && !result[i].param_defined) {
				double val = 0.1 * (int)ml::genRand(1, 5);
				next += Literal(result[i].c, val);

				node = root->getChild(val);
				found = true;
			} else if ((result[i].c == '-' || result[i].c == '+') && !result[i].param_defined) {
				double val = 10 * (int)ml::genRand(1, 5);
				next += Literal(result[i].c, val);
				node = root->getChild(val);
				found = true;
			} else {
				next += result[i];
			}
		}

		// 新たなderivationがないなら、終了
		if (!found) break;

		result = next;
	}

	return result;
}

void ParametricLSystem::draw() {
	//drawSegment(rule);
}

void ParametricLSystem::drawTree() {
	drawSubTree(root, 0);
}

void ParametricLSystem::drawSubTree(TreeNode* node, int indent) {
	for (int i = 0; i < indent; ++i) cout << " ";
	cout << "+" << node->val << endl;
}

/*
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
			drawCylinder(state.modelMat, state.diameter * 0.5, state.diameter * 0.5, length, state.color);
			state.modelMat = glm::translate(state.modelMat, glm::vec3(0, length, 0));
		} else {
			i++;
		}
	}
}
*/

/**
 * XZ面を底面とし、Y軸に伸びる円筒形を描画する。
 *
 * @modelMat		モデル行列
 * @top_radius		円筒形の上面の半径
 * @base_radius		円筒形の底面の半径
 * @height			円筒形の高さ
 * @color			色
 */
/*
void ParametricLSystem::drawCylinder(const glm::mat4& modelMat, float top_radius, float base_radius, float height, const glm::vec3& color) {
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
*/

/**
 * XY平面上に、原点を中心とする楕円を描画する。
 *
 * @modelMat		モデル行列
 * @length			楕円のY軸方向の直径
 * @width			楕円のX軸方向の直径
 * @color			色
 */
/*
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
*/

/**
 * ルールリストから、確率に基づいて１つのルールを選択する。
 * リストの各要素は、<確率、ルール>のペアとなっている。
 *
 * @param rules		ルールリスト
 * @reutnr			選択されたルール
 */
int ParametricLSystem::chooseRule(char non_terminal) {
	return rand() % rules[non_terminal].size();
}

float deg2rad(float deg) {
	return deg * M_PI / 180.0;
}

}
