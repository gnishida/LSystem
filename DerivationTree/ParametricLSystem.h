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
	double diameter;
	glm::mat4 modelMat;
	glm::vec3 color;

public:
	State() : diameter(0.0), color(0, 0.7, 0) {};
};

class Literal {
public:
	char c;
	double param_value;
	bool param_defined;

public:
	Literal(char c);
	Literal(char c, double param_value);
};

class String {
public:
	vector<Literal> str;

public:
	String() {}
	String(string str);

	int length() const { return str.size(); }
	Literal operator[](int index) const { return str[index]; }
	Literal& operator[](int index) { return str[index]; }
	void operator+=(const Literal& l) { str.push_back(l); }
	void operator+=(const string& str);
};

ostream& operator<<(ostream& os, const String& dt);

class Rule {
public:
	string right_hand;

public:
	Rule() {}
	Rule(const string& right_hand);
};

class TreeNode {
public:
	double val;
	TreeNode* parent;
	map<int, TreeNode*> children;

public:
	TreeNode(double val, TreeNode* parent);
	TreeNode* getChild(double value);
};


class ParametricLSystem {
public:
	string axiom;
	map<char, vector<Rule> > rules;
	TreeNode* root;

public:
	ParametricLSystem();
	String derive(int random_seed);
	void draw();
	void drawTree();
	void drawSubTree(TreeNode* node, int indent);

private:
	/*
	void drawSegment(string rule);
	void drawCylinder(const glm::mat4& modelMat, float top_radius, float base_radius, float height, const glm::vec3& color);
	void drawCircle(const glm::mat4& modelMat, float length, float width, const glm::vec3& color);
	*/
	int chooseRule(char non_terminal);
};

float deg2rad(float deg);

}
