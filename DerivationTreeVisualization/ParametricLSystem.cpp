#include "ParametricLSystem.h"
#include <iostream>
#include <time.h>
#include <algorithm>
#include "MLUtils.h"
#include <list>

#define MAX_ITERATIONS						300//200
#define MAX_LEVEL							6

namespace parametriclsystem {

const double M_PI = 3.141592653592;

Literal::Literal(char c, int level) {
	if (level > 10) {
		int xx =0;
	}
	this->c = c;
	this->level = level;
	this->param_value = 0.0;
	this->param_defined = false;
}

Literal::Literal(char c, int level, double param_value) {
	this->c = c;
	this->level = level;
	this->param_value = param_value;
	this->param_defined = true;
}

String::String(string str, int level) {
	for (int i = 0; i < str.length(); ++i) {
		this->str.push_back(Literal(str[i], level));
	}
}

String::String(Literal l) {
	this->str.push_back(l);
}

void String::operator+=(const String& str) {
	for (int i = 0; i < str.length(); ++i) {
		this->str.push_back(str[i]);
	}
}

String String::operator+(const String& str) const {
	String new_str = *this;

	for (int i = 0; i < str.length(); ++i) {
		new_str.str.push_back(str[i]);
	}

	return new_str;
}

void String::replace(int index, const String& str) {
	this->str.erase(this->str.begin() + index);
	this->str.insert(this->str.begin() + index, str.str.begin(), str.str.end());
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

int TreeNode::seq = 0;

TreeNode::TreeNode(const Literal& l, TreeNode* parent) {
	this->id = seq++;
	this->l = l;
	this->parent = parent;
}

TreeNode* TreeNode::getChild(char c, int level, double value) {
	if (children.find(value) == children.end()) {
		children[value] = new TreeNode(Literal(c, level, value), this);
	}

	return children[value];
}

ParametricLSystem::ParametricLSystem(int grid_size, int indicator_data_type, float scale) {
	this->grid_size = grid_size;
	this->indicator_data_type = indicator_data_type;
	this->scale = scale;

	axiom = "X";
	rules['X'].push_back("F");
	rules['X'].push_back("X[-X]X[+X]X");
	//rules['X'].push_back(Rule("F[-F]"));

	// ルートノードを作成
	root = new TreeNode(Literal('#', 0), NULL);
	root->model = String(axiom, 0);
}

/**
 * 適当な乱数シードに基づいて、ランダムにgenerateする。
 *
 * @param random_seed		乱数シード
 * @return					生成されたモデル
 */
String ParametricLSystem::derive(int random_seed) {
	cv::Mat indicator;
	return derive(String(axiom, 0), random_seed, MAX_ITERATIONS, true, indicator);
}

void ParametricLSystem::generateDerivationStateTree(int num_samples) {
	for (int i = 0; i < num_samples; ++i) {
		cv::Mat indicator;
		String model = derive(root->model, i, MAX_ITERATIONS, true, indicator);

		/*
		// sampleの画像を保存する
		char filename[256];
		sprintf(filename, "images/sample%d.png", i);

		cv::Mat img;
		computeIndicator(model, 1.0f, img);
		ml::mat_save(filename, img);
		*/
	}

	cout << "sampling done." << endl;

	FILE* fp = fopen("tree.txt", "w");
	gatherIndicators(root, fp);
	fclose(fp);
}

/**
 * 指定された開始モデルからスタートし、適当な乱数シードに基づいて、ランダムにgenerateする。
 *
 * @param start_model		開始モデル
 * @param random_seed		乱数シード
 * @param max_iterations	繰り返し数
 * @param build_tree		木を生成するか？
 * @param indicator [OUT]	生成されたモデルのindicator
 * @return					生成されたモデル
 */
String ParametricLSystem::derive(const String& start_model, int random_seed, int max_iterations, bool build_tree, cv::Mat& indicator) {
	srand(random_seed);

	String result = start_model;
	TreeNode* node = root;

	//for (int iter = 0; iter < max_iterations; ++iter) {
	for (int iter = 0; ; ++iter) {
		// 展開するパラメータを決定
		int i = findNextLiteralToDefineValue(result);

		// 新たなderivationがないなら、終了
		if (i == -1) break;

		if (rules.find(result[i].c) != rules.end()) {
			int index = chooseRule(result[i]);

			// 最大繰り返し数を超えたら、X->Fにする
			if (iter > max_iterations) {
				index = 0;
			}

			if (build_tree) {
				node = node->getChild(result[i].c, result[i].level + 1, index);
			}

			result.replace(i, String(rules[result[i].c][index], result[i].level + 1));
		} else if (result[i].c == 'F') {
			//double val = grid_size / pow(3.0, result[i].level - 1);
			double val = 25;
			result[i] = Literal(result[i].c, result[i].level, val);
			if (build_tree) {
				node = node->getChild(result[i].c, result[i].level, val);
			}
		} else if (result[i].c == '-' || result[i].c == '+') {
			double val = ml::genRandInt(10, 60, 3);
			result[i] = Literal(result[i].c, result[i].level, val);
			if (build_tree) {
				node = node->getChild(result[i].c, result[i].level, val);
			}
		}

		node->model = result;
	}

	// indicatorを計算する
	if (build_tree) {
		computeIndicator(result, scale, node->indicator);
		indicator = node->indicator.clone();
	} else {
		computeIndicator(result, scale, indicator);
	}

	return result;
}

/**
 * ツリーを表示する。
 *
 * @param max_depth		最大、この深さまで表示する
 */
void ParametricLSystem::drawTree(int max_depth) {
	drawSubTree(root, 0, max_depth);
}

/**
 * 指定された文字列に基づいてモデルを生成し、indicatorを計算して返却する。
 * 
 * @param rule				モデルを表す文字列
 * @param scale				grid_size * scaleのサイズでindicatorを計算する
 * @param indicator [OUT]	indicator
 */
void ParametricLSystem::computeIndicator(String rule, float scale, cv::Mat& indicator) {
	int size = grid_size * scale;

	indicator = cv::Mat::zeros(size, size, indicator_data_type);

	std::list<glm::mat4> stack;

	glm::mat4 modelMat;

	for (int i = 0; i < rule.length(); ++i) {
		if (rule[i].c == '[') {
			stack.push_back(modelMat);
		} else if (rule[i].c == ']') {
			modelMat = stack.back();
			stack.pop_back();
		} else if (rule[i].c == '+') {
			modelMat = glm::rotate(modelMat, deg2rad(rule[i].param_value), glm::vec3(0, 0, 1));
		} else if (rule[i].c == '-') {
			modelMat = glm::rotate(modelMat, deg2rad(-rule[i].param_value), glm::vec3(0, 0, 1));
		} else if (rule[i].c == 'F') {
			double length = rule[i].param_value * scale;
			
			// 線を描画する代わりに、indicatorを更新する
			glm::vec4 p1(0, 0, 0, 1);
			glm::vec4 p2(0, length, 0, 1);
			p1 = modelMat * p1;
			p2 = modelMat * p2;
			int u1 = p1.x + size * 0.5 + 0.5;
			int v1 = p1.y + 0.5;
			int u2 = p2.x + size * 0.5 + 0.5;
			int v2 = p2.y + 0.5;
			int thickness = max(1.0, 3.0 * scale);
			cv::line(indicator, cv::Point(u1, v1), cv::Point(u2, v2), cv::Scalar(1), thickness);

			/*
			for (float y = 0.0f; y <= length; y += 0.1f) {
				glm::vec4 p(0, y, 0, 1);
				p = modelMat * p;

				int c = p.x + size * 0.5 + 0.5;
				int r = p.y + 0.5;
				if (c >= 0 && c < size && r >= 0 && r < size) {
					ml::mat_set_value(indicator, r, c, 1);
				}
			}
			*/

			modelMat = glm::translate(modelMat, glm::vec3(0, length, 0));
		} else {
		}
	}
}

/**
 * 現在のモデルについて、値が未設定のパラメータで最もレベルが低く、最も左側に位置するindexを返却する。
 * もし未設定のパラメータがない場合は、-1を返却する。
 *
 * @param str		現在のモデル
 * @return			インデックス
 */
int ParametricLSystem::findNextLiteralToDefineValue(const String& str) {
	int min_level1 = std::numeric_limits<int>::max();
	int min_level2 = std::numeric_limits<int>::max();
	int min_index1 = -1;
	int min_index2 = -1;

	for (int i = 0; i < str.length(); ++i) {
		if (rules.find(str[i].c) != rules.end()) {
			if (str[i].level < min_level1) {
				min_level1 = str[i].level;
				min_index1 = i;
			}
		} else if ((str[i].c == 'F' || str[i].c == '+' || str[i].c == '-') && !str[i].param_defined) {
			if (str[i].level < min_level2) {
				min_level2 = str[i].level;
				min_index2 = i;
			}
		}
	}

	if (min_level1 < min_level2) {
		return min_index1;
	} else {
		return min_index2;
	}
}

/**
 * サブツリーを表示する。
 *
 * @param node		このノード以下のサブツリーを表示する。
 * @param depth		現在の深さ
 * @param max_depth この深さまで表示する
 */
void ParametricLSystem::drawSubTree(TreeNode* node, int depth, int max_depth) {
	if (depth >= max_depth - 1) return;

	for (auto it = node->children.begin(); it != node->children.end(); ++it) {
		double val = it->first;
		TreeNode* child = it->second;
		for (int i = 0; i < depth; ++i) cout << " ";
		cout << "+" << child->l.c << "(" << child->l.param_value << ") (*" << child->id << ")" << endl;
		drawSubTree(child, depth + 1, max_depth);
	}
}

/**
 * このノード以下のindicatorを加重平均して、返却する。
 *
 * @param node			ノード
 * @param gather_type	0 -- 加重平均 / 1 -- Max / 2 -- Min
 * @return				加重平均indicator
 */
cv::Mat ParametricLSystem::gatherIndicators(TreeNode* node, FILE* fp) {
	int num = 0;

	cv::Mat indicator;
	computeIndicator(node->model, scale, indicator);
	indicator = indicator.reshape(1, 1);

	for (auto it = node->children.begin(); it != node->children.end(); ++it) {
		TreeNode* child = it->second;

		cv::Mat indicator2 = gatherIndicators(child, fp);

		cv::vconcat(indicator, indicator2, indicator);
	}

	if (node->id > 0) {
		fprintf(fp, "%d\t%d\t%lf\n", node->id, node->parent->id, ml::mat_variance(indicator));
	}

	return indicator;
}

/**
 * ルールリストから、確率に基づいて１つのルールを選択する。
 * リストの各要素は、<確率、ルール>のペアとなっている。
 *
 * @param rules		ルールリスト
 * @reutnr			選択されたルール
 */
int ParametricLSystem::chooseRule(const Literal& non_terminal) {
	// ハードコーディング
	// 深さMAX_LEVELを超えたら、X->Fとする
	if (non_terminal.level > MAX_LEVEL) return 0;

	// ハードコーディング
	// 深さ/MAX_LEVEL の確率で、X->Fとする
	if (rand() % MAX_LEVEL <= non_terminal.level) {
		return 0;
	} else {
		return 1;
	}



	//return rand() % rules[non_terminal.c].size();
}

float deg2rad(float deg) {
	return deg * M_PI / 180.0;
}

}
