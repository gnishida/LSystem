#include "LSystem.h"
#include <QGLWidget>
#include <iostream>
#include <time.h>

namespace lsystem {

const double M_PI = 3.141592653592;
double LSystem::GRID_SIZE = 300.0;
int LSystem:: NUM_GRID = 5;
double LSystem::CELL_SIZE = GRID_SIZE / NUM_GRID;


LSystem::LSystem() {
	N = 5;
	delta = 25.7;
	axiom = 'F';
	segment_length = GRID_SIZE / pow(3.0, N);
	rules['F'].push_back(pair<double, string>(1.0, "F[+F]F[-F]F"));

	randomInit(0);
}

void LSystem::draw() {
	// 統計情報をクリア
	stats.density = cv::Mat_<int>::zeros(NUM_GRID, NUM_GRID);
	
	drawSegment(glm::mat4(), 0, string(1, axiom));

	// densityが0のセルについて、PMパラメータを真ん中の値にする
	// * densityが0の場合、PMパラメータが全く寄与しないため、
	// * indicatorからPMパラメータを予測することは不可能だから。
	for (int r = 0; r < stats.density.rows; ++r) {
		for (int c = 0; c < stats.density.cols; ++c) {
			if (stats.density(r, c) == 0) {
				deltas(r, c) = 30;
				levels(r, c) = 5;
			}
		}
	}
}

void LSystem::randomInit(int seed) {
	vector<unsigned> seeds(1);
	seeds[0] = seed;
	std::seed_seq seq(seeds.begin(), seeds.end());
	mt.seed(seq);

	// グローバルコントロールパラメータをランダムにセット
	deltas = cv::Mat_<float>::zeros(NUM_GRID, NUM_GRID);
	levels = cv::Mat_<int>::zeros(NUM_GRID, NUM_GRID);
	for (int r = 0; r < deltas.rows; ++r) {
		for (int c = 0; c < deltas.cols; ++c) {
			deltas(r, c) = genRand(10, 50);
		}
	}
	for (int r = 0; r < levels.rows; ++r) {
		for (int c = 0; c < levels.cols; ++c) {
			levels(r, c) = genRand(3, 7);
		}
	}
}

/**
 * 指定された行列に格納されたパラメータ値をセットする。
 *
 * @param mat		パラメータ値が格納された行列
 */
void LSystem::setParams(const cv::Mat_<float>& mat) {
	cv::Mat_<float> m;
	if (mat.rows == 1) {
		m = mat.t();
	} else {
		m = mat.clone();
	}

	int index = 0;
	for (int r = 0; r < deltas.rows; ++r) {
		for (int c = 0; c < deltas.cols; ++c) {
			deltas(r, c) = m(index++, 0);
		}
	}
	for (int r = 0; r < levels.rows; ++r) {
		for (int c = 0; c < levels.cols; ++c) {
			levels(r, c) = m(index++, 0);
		}
	}

	// Hard constraintsに従って、値を修正する
}

/**
 * パラメータの配列を返却する。
 *
 * @return		パラメータの配列
 */
vector<float> LSystem::getParams() {
	vector<float> ret(deltas.rows * deltas.cols + levels.rows * levels.cols);

	int index = 0;
	for (int r = 0; r < deltas.rows; ++r) {
		for (int c = 0; c < deltas.cols; ++c) {
			ret[index++] = deltas(r, c);
		}
	}
	for (int r = 0; r < levels.rows; ++r) {
		for (int c = 0; c < levels.cols; ++c) {
			ret[index++] = levels(r, c);
		}
	}
	return ret;
}

vector<float> LSystem::getStatistics() {
	vector<float> ret(stats.density.rows * stats.density.cols);
	int index = 0;
	for (int r = 0; r < stats.density.rows; ++r) {
		for (int c = 0; c < stats.density.cols; ++c) {
			ret[index++] = stats.density(r, c);
		}
	}

	return ret;
}

void LSystem::drawSegment(glm::mat4& modelMat, int level, string rule) {
	std::list<glm::mat4> stack;

	for (int i = 0; i < rule.length(); ++i) {
		glm::vec4 p(0, 0, 0, 1);
		p = modelMat * p;
		pair<int, int> uv = XYtoUV(p.x, p.y);

		if (rule[i] == '[') {
			stack.push_back(modelMat);
		} else if (rule[i] == ']') {
			modelMat = stack.back();
			stack.pop_back();
		} else if (rule[i] == '+') {
			if (uv.first >= 0 && uv.first < NUM_GRID && uv.second >= 0 && uv.second < NUM_GRID) {
				modelMat = glm::rotate(modelMat, deg2rad(deltas(uv.second, uv.first)), glm::vec3(0, 0, 1));
			}
		} else if (rule[i] == '-') {
			if (uv.first >= 0 && uv.first < NUM_GRID && uv.second >= 0 && uv.second < NUM_GRID) {
				modelMat = glm::rotate(modelMat, deg2rad(-deltas(uv.second, uv.first)), glm::vec3(0, 0, 1));
			}
		} else if (rule[i] == '\\') {
			if (uv.first >= 0 && uv.first < NUM_GRID && uv.second >= 0 && uv.second < NUM_GRID) {
				modelMat = glm::rotate(modelMat, deg2rad(deltas(uv.second, uv.first)), glm::vec3(0, 1, 0));
			}
		} else if (rule[i] == '/') {
			if (uv.first >= 0 && uv.first < NUM_GRID && uv.second >= 0 && uv.second < NUM_GRID) {
				modelMat = glm::rotate(modelMat, deg2rad(-deltas(uv.second, uv.first)), glm::vec3(0, 1, 0));
			}
		} else if (rule[i] == '&') {
			if (uv.first >= 0 && uv.first < NUM_GRID && uv.second >= 0 && uv.second < NUM_GRID) {
				modelMat = glm::rotate(modelMat, deg2rad(deltas(uv.second, uv.first)), glm::vec3(1, 0, 0));
			}
		} else if (rule[i] == '^') {
			if (uv.first >= 0 && uv.first < NUM_GRID && uv.second >= 0 && uv.second < NUM_GRID) {
				modelMat = glm::rotate(modelMat, deg2rad(-deltas(uv.second, uv.first)), glm::vec3(1, 0, 0));
			}
		} else if (rule[i] == '|') {
			if (uv.first >= 0 && uv.first < NUM_GRID && uv.second >= 0 && uv.second < NUM_GRID) {
				modelMat = glm::rotate(modelMat, deg2rad(180), glm::vec3(0, 0, 1));
			}
		} else {
			if (rules.find(rule[i]) != rules.end()) {								
				if (uv.first >= 0 && uv.first < NUM_GRID && uv.second >= 0 && uv.second < NUM_GRID && level < levels(uv.second, uv.first)) {
					drawSegment(modelMat, level + 1, chooseRule(rules[rule[i]]));
				} else {
					drawCylinder(modelMat, 1, 1, segment_length, glm::vec3(1.0, 1.0, 1.0));
					modelMat = glm::translate(modelMat, glm::vec3(0, segment_length, 0));
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

	// 統計情報を更新
	{
		glm::vec4 p1(0, 0, 0, 1);
		glm::vec4 p2(0, height, 0, 1);
		p1 = modelMat * p1;
		p2 = modelMat * p2;
		pair<int, int> uv = XYtoUV((p1.x + p2.x) * 0.5, (p1.y + p2.y) * 0.5);
		if (uv.first >= 0 && uv.first < NUM_GRID && uv.second >= 0 && uv.second < NUM_GRID) {
			stats.density(uv.second, uv.first)++;
		}
	}
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

	std::uniform_real_distribution<> randu(0, cdf.back());

	double rnd = randu(mt);
	for (int i = 0; i < cdf.size(); ++i) {
		if (rnd <= cdf[i]) {
			return rules[i].second;
		}
	}

	return rules.back().second;
}

pair<int, int> LSystem::XYtoUV(float x, float y) {
	int u = (x + GRID_SIZE * 0.5) / CELL_SIZE;
	int v = y / CELL_SIZE;

	return make_pair(u, v);
}

float LSystem::deg2rad(float deg) {
	return deg * M_PI / 180.0;
}

float LSystem::genRand() {
	std::uniform_real_distribution<> randu(0, 1);
	return randu(mt);
}

float LSystem::genRand(float a, float b) {
	std::uniform_real_distribution<> r(a, b); 
	return r(mt);
}

}
