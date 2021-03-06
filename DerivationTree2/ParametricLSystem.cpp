﻿#include "ParametricLSystem.h"
#include <iostream>
#include <time.h>
#include <algorithm>
#include "MLUtils.h"
#include <list>
#include <boost/filesystem.hpp>

#define MAX_ADDITIONAL_DEPTH_FOR_ESTIMATE	5
#define NUM_RANDOM_GENERATION_FOR_ESTIMATE	200
#define MAX_DEPTH							6

namespace parametriclsystem {

const double M_PI = 3.141592653592;

Literal::Literal(char c, int level) {
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
	rules['X'].push_back("F[-X][+X]");
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
	return derive(String(axiom, 0), random_seed, MAX_DEPTH, true, indicator);
}

/**
 * 指定された開始モデルからスタートし、適当な乱数シードに基づいて、ランダムにgenerateする。
 *
 * @param start_model		開始モデル
 * @param random_seed		乱数シード
 * @param max_depth			最大深さ
 * @param build_tree		木を生成するか？
 * @param indicator [OUT]	生成されたモデルのindicator
 * @return					生成されたモデル
 */
String ParametricLSystem::derive(const String& start_model, int random_seed, int max_depth, bool build_tree, cv::Mat& indicator) {
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

			// もしmax_depth以上ら、X->Fにする
			if (result[i].level >= max_depth) {
				index = 0;
			}

			if (build_tree) {
				node = node->getChild(result[i].c, result[i].level + 1, index);
			}

			result.replace(i, String(rules[result[i].c][index], result[i].level + 1));
		} else if (result[i].c == 'F') {
			//double mean_val = grid_size * 0.5 / (result[i].level + 1);
			double mean_val = grid_size * 0.5 / pow(1.5, result[i].level);
			//double mean_val = grid_size * 0.15;
			double val = ml::genRandInt(mean_val * 0.8, mean_val * 1.2, 3);
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
 * 指定されたモデルからスタートし、ターゲットに近づくようモデルをgenerateする。
 *
 * @param start_model		初期モデル
 * @param max_depth			最大深さ
 * @target					ターゲット
 * @indicator [OUT]			生成されたモデルのindicator
 * @return					生成されたモデル
 */
String ParametricLSystem::derive(const String& start_model, int max_depth, const cv::Mat& target, cv::Mat& indicator) {
	String result = start_model;

	for (int iter = 0; ; ++iter) {
		String min_next;

		// 展開するパラメータを決定
		int i = findNextLiteralToDefineValue(result);

		// 新たなderivationがないなら、終了
		if (i == -1) break;

		if (rules.find(result[i].c) != rules.end()) {
			double min_dist = std::numeric_limits<double>::max();

			for (int k = 0; k < rules[result[i].c].size(); ++k) {
				// この値を選択した時のモデルを作成
				String next = result;
				next.replace(i, String(rules[result[i].c][k], result[i].level + 1));

				// indicatorを推定
				cv::Mat indicator;
				estimateIndicator(next, scale, result[i].level + 1 + MAX_ADDITIONAL_DEPTH_FOR_ESTIMATE, target, indicator);

				// distanceを計算
				ml::mat_save("hoge.png", indicator);
				double dist = distance(indicator, target);

				if (dist < min_dist) {
					min_dist = dist;
					min_next = next;
				}

				// もしmax_depth以上ら、X->Fにする
				if (result[i].level >= max_depth) break;
			}
		} else if (result[i].c == 'F') {
			double min_dist = std::numeric_limits<double>::max();
					 
			for (int k = 0; k < 3; ++k) {
				//double mean_val = grid_size * 0.5 / (result[i].level + 1);
				double mean_val = grid_size * 0.5 / pow(1.5, result[i].level);
				//double mean_val = grid_size * 0.15;
				//double val = grid_size * 0.5 / (result[i].level + 1) * (0.8 + 0.2 * k);
				double val = mean_val * (0.8 + 0.2 * k);

				// この値を選択した時のモデルを生成
				String next = result;
				next[i] = Literal(result[i].c, result[i].level, val);

				// indicatorを計算
				cv::Mat indicator;
				estimateIndicator(next, scale, result[i].level + MAX_ADDITIONAL_DEPTH_FOR_ESTIMATE, target, indicator);

				// distanceを計算
				double dist = distance(indicator, target);

				if (dist < min_dist) {
					min_dist = dist;
					min_next = next;
				}
			}
		} else if (result[i].c == '-' || result[i].c == '+') {
			double min_dist = std::numeric_limits<double>::max();
					 
			for (int k = 0; k < 3; ++k) {
				double val = k * 25.0 + 10.0;

				// この値を選択した時のモデルを生成
				String next = result;
				next[i] = Literal(result[i].c, result[i].level, val);

				// indicatorを計算
				cv::Mat indicator;
				estimateIndicator(next, scale, result[i].level + MAX_ADDITIONAL_DEPTH_FOR_ESTIMATE, target, indicator);

				// distanceを計算
				double dist = distance(indicator, target);

				if (dist < min_dist) {
					min_dist = dist;
					min_next = next;
				}
			}
		}

		result = min_next;
	}

	computeIndicator(result, scale, indicator);

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
 * ツリーの各ノードの平均indicatorを計算する。
 */
void ParametricLSystem::gatherIndicators(int gather_type) {
	countLeaves(root);
	gatherSubIndicators(root, gather_type);
}

/**
 * 各ノードの加重平均indicatorをファイルに保存する
 *
 * @param max_depth			この深さを超えると、ファイルに保存しない
 * @param min_threshold		0より大きい値を指定した場合、この値より小さい値は0、以上の値は1とする
 */
void ParametricLSystem::saveIndicatorImages(int max_depth, float min_threshold) {
	boost::filesystem::path dir("images");
	boost::filesystem::create_directory(dir);

	saveSubIndicatorImages(root, 0, max_depth, min_threshold);
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
 * 指定されたモデルからスタートして、ランダムにいくつかサンプルを生成し、ターゲットに最も近いindicatorを返却する。
 *
 * @param start_model		開始モデル
 * @param scale				grid_size * scaleのサイズで、indicatorを計算する
 * @param target			ターゲットindicator
 * @param indicator [OUT]	indicator
 */
void ParametricLSystem::estimateIndicator(const String start_model, float scale, int max_depth, const cv::Mat& target, cv::Mat& indicator) {
	int size = grid_size * scale;

	int N = NUM_RANDOM_GENERATION_FOR_ESTIMATE;

	double min_dist = std::numeric_limits<double>::max();
	for (int i = 0; i < N; ++i) {
		cv::Mat ind;
		derive(start_model, i, max_depth, false, ind);

		double dist = distance(ind, target);
		if (dist < min_dist) {
			min_dist = dist;
			indicator = ind.clone();
		}
	}
}

/**
 * 指定されたターゲットindicatorに近いモデルを生成する。
 *
 * @param target			ターゲットindicator
 * @param threshold			しきい値
 * @param indicator [OUT]	生成されたモデルのindicator
 * @return					生成されたモデル
 */
String ParametricLSystem::inverse(const cv::Mat& target, cv::Mat& indicator) {
	TreeNode* node = traverseTree(root, target);

	computeIndicator(node->model, scale, indicator);
	ml::mat_save("start_model.png", indicator);
	cout << "===================================" << endl;
	cout << "Start model:" << endl;
	cout << node->model << endl;
		
	return derive(node->model, MAX_DEPTH, target, indicator);
}

/**
 * indicatorとターゲットとの距離を計算して返却する。
 *
 * @param indicator		indicator
 * @param target		ターゲットindicator
 * @param threshold		しきい値
 * @return				距離
 */
double ParametricLSystem::distance(const cv::Mat& indicator, const cv::Mat& target) {
	double dist = 0.0;

	for (int r = 0; r < target.rows; ++r) {
		for (int c = 0; c < target.cols; ++c) {
			double target_val = ml::mat_get_value(target, r, c);
			double val = ml::mat_get_value(indicator, r, c);

			if (target_val > 0) {
				dist -= val * 100;
			} else {
				dist += val * 10;
			}
		}
	}

	return dist;

	/*
	// indicatorをblurする
	cv::Mat bluredIndicator;
	cv::GaussianBlur(indicator, bluredIndicator, cv::Size(5, 5), 3, 3);

	// maskを作成
	cv::Mat mask = (bluredIndicator >= threshold) / 255;
	cv::Mat mask2(bluredIndicator.size(), bluredIndicator.type());
	mask.convertTo(mask2, bluredIndicator.type());

	cout << bluredIndicator << endl;
	cout << mask2 << endl;
	cout << bluredIndicator.mul(mask2) << endl;
	cout << target << endl;
	cout << target.mul(mask2) << endl;

	double dist = ml::mat_squared_sum(bluredIndicator.mul(mask2) - target.mul(mask2));

	return dist;
	*/
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
 * このノード以下にある葉ノードの数を更新し、返却する。
 *
 * @param node		ノード
 * @return			葉ノードの数を返却する
 */
int ParametricLSystem::countLeaves(TreeNode* node) {
	if (node->children.size() == 0) {
		node->num_leaves = 1;
	} else {
		node->num_leaves = 0;
		for (auto it = node->children.begin(); it != node->children.end(); ++it) {
			TreeNode* child = it->second;
			node->num_leaves += countLeaves(child);
		}
	}

	return node->num_leaves;
}

/**
 * このノード以下のindicatorを加重平均して、返却する。
 *
 * @param node			ノード
 * @param gather_type	0 -- 加重平均 / 1 -- Max / 2 -- Min
 * @return				加重平均indicator
 */
cv::Mat ParametricLSystem::gatherSubIndicators(TreeNode* node, int gather_type) {
	int num = 0;

	for (auto it = node->children.begin(); it != node->children.end(); ++it) {
		TreeNode* child = it->second;

		if (num == 0) {
			node->indicator = gatherSubIndicators(child, gather_type) * child->num_leaves;
		} else {
			if (gather_type == 0) {
				node->indicator += gatherSubIndicators(child, gather_type) * child->num_leaves;
			} else {
				if (gather_type == 1) {
					node->indicator = ml::mat_max(node->indicator, gatherSubIndicators(child, gather_type));
				} else {
					node->indicator = ml::mat_min(node->indicator, gatherSubIndicators(child, gather_type));
				}

				/*
				cv::Mat indicator2 = gatherSubIndicators(child, gather_type);
				for (int r = 0; r < node->indicator.rows; ++r) {
					for (int c = 0; c < node->indicator.cols; ++c) {
						double v1 = ml::mat_get_value(node->indicator, r, c);
						double v2 = ml::mat_get_value(indicator2, r, c);
						if (gather_type == 1 && v2 > v1) {
							ml::mat_set_value(node->indicator, r, c, v2);
						} else if (gather_type == 2 && v1 > v2) {
							ml::mat_set_value(node->indicator, r, c, v2);
						}
					}
				}
				*/
			}
		}
		num += child->num_leaves;
	}

	if (num > 0 && gather_type == 0) {
		node->indicator /= (double)num;
	}

	return node->indicator;
}

/**
 * このノードの加重平均indicatorをファイルに保存する。
 *
 * @param node				ノード
 * @param depth				深さ
 * @param max_depth			この深さを超えたら、ファイルに保存しない
 * @param min_threshold		0より大きい値を指定した場合、この値より小さい値は0、以上の値は1とする
 */
void ParametricLSystem::saveSubIndicatorImages(TreeNode* node, int depth, int max_depth, float min_threshold) {
	if (depth >= max_depth) return;

	for (auto it = node->children.begin(); it != node->children.end(); ++it) {
		TreeNode* child = it->second;

		saveSubIndicatorImages(child, depth + 1, max_depth, min_threshold);
	}

	// このノードのindicatorをファイルに保存する
	char filename[256];
	sprintf(filename, "images/indicator_%d.png", node->id);

	if (min_threshold == 0) {
		ml::mat_save(filename, node->indicator);
	} else {
		ml::mat_save(filename, ml::mat_threhold(node->indicator, min_threshold));
	}

}

/**
 * 現在のノードから下に辿り、指定されたターゲットにもっとも近いノードを返却する。
 *
 * @param node			現在のノード
 * @param target		ターゲットindicator
 * @param threshold		しきい値
 * @return				もっともターゲットに近いノード
 */
TreeNode* ParametricLSystem::traverseTree(TreeNode* node, const cv::Mat& target) {
	int node_type = 0;
	int count = 0;
	for (auto it = node->children.begin(); it != node->children.end(); ++it) {
		TreeNode* child = it->second;

		if (child->l.c == 'F' || child->l.c == '-' || child->l.c == '+') {
			node_type = 1;
		}

		count++;
	}

	// 全てのオプション値に対応する子ノードがなければ、このノードでtraverseを終了する
	if ((node_type == 0 && count < 2) || (node_type == 1 && count < 3)) {
		return node;
	}

	// 子ノードの中から、indicatorが最もtargetに近いものを選び、traverseする
	double min_dist = std::numeric_limits<double>::max();
	TreeNode* min_child;
	for (auto it = node->children.begin(); it != node->children.end(); ++it) {
		TreeNode* child = it->second;

		double dist = distance(child->indicator, target);
		if (dist < min_dist) {
			min_dist = dist;
			min_child = child;
		}
	}

	return traverseTree(min_child, target);
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
	// 深さ6を超えたら、X->Fとする
	if (non_terminal.level > 6) return 0;

	// ハードコーディング
	// 深さ/6 の確率で、X->Fとする
	if (rand() % 6 <= non_terminal.level) {
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
