#include "LSystem.h"
#include <QGLWidget>
#include <iostream>
#include <time.h>

namespace lsystem {

const double M_PI = 3.141592653592;
double LSystem::GRID_SIZE = 300.0;
int LSystem::NUM_GRID = 5;
int LSystem::NUM_STAT_GRID = 10;

LSystem::LSystem() {
	N = 5;
	delta = 25.7;
	axiom = 'F';
	segment_length = GRID_SIZE / pow(3.0, N);
	rules['F'].push_back(pair<double, string>(1.0, "F[+F]F[-F]F"));

	randomInit(0);

	//cv::Mat_<double> m = (cv::Mat_<double>(1, 300) << 40.303252957292, 45.42651318637027, 45.85918142739293, 45.31599437049659, 40.20883135983765, 43.36990718178729, 45.00000000000006, 45.08233711354288, 44.83639136421227, 44.94791645589248, 44.61873500626289, 44.41451733791551, 43.09047344975868, 47.01686452876735, 44.17499784944837, 55.34026399731, 43.89837940446548, 45.622238327642, 44.45726496311552, 45.5537313618275, 45.93829035909113, 49.19912489317675, 47.9387493283275, 44.70746843059937, 48.94580814400437, 51.7937023905646, 43.5094938738287, 45.13197340305391, 45.51916463162067, 43.30459563515249, 45.47082864566377, 39.74088075316392, 41.77084861913104, 43.74361971072989, 44.35213657009608, 54.74233889740216, 50.37209044274015, 44.14010989177534, 44.45922744111154, 44.6007598932827, 39.57727893709945, 48.55143411264957, 42.60862442880845, 41.77261346680339, 48.05783589413782, 43.4571117924953, 42.97352676167057, 43.34972869856996, 48.42908771283724, 41.52376646926615, 41.33682139254617, 48.5843560589435, 42.19915035251733, 48.74092408324188, 39.72358270364049, 37.85830128240387, 43.69346620897605, 45.59578259444214, 46.94899580859627, 46.35493802300608, 49.95519006486925, 41.80684073205814, 55.41105127246206, 47.36774592869517, 45.76458995915002, 55.44330931950103, 42.77594512177315, 43.39125394005441, 45.0101739547796, 45.98074661103055, 52.18301924602761, 49.62083303280698, 48.4699517010096, 46.16721241915786, 40.80171988722839, 37.99557087268596, 41.76099617450774, 45.73870116056646, 47.84439284976989, 44.80083389898825, 43.33143504210584, 46.47930420056611, 43.91354461270019, 44.70858219675651, 39.8732388664087, 38.47036280123584, 43.7885515878206, 43.1294869301937, 45.67588479469301, 44.79225663563654, 47.20320950868727, 45.36916187347739, 47.69789401246765, 44.48656539034306, 50.44145486254065, 42.91776484625974, 45.72600189175505, 45.7096166599132, 45.35132168389525, 44.69070682435217, 5.276616700771267, 5.257641890381032, 4.9521929090733, 4.892976875316679, 5.184863056873061, 3.829792125839229, 5.000000000000015, 4.98878966146818, 4.987960533607775, 4.996315042259306, 5.158956953051689, 5.18247605971221, 5.132772940244786, 4.438008581883933, 5.703940234843035, 3.95732643170182, 4.998799622641927, 4.900376402977019, 4.983898725323124, 4.975621968167972, 5.178671031640198, 5.063736027454675, 4.821149345369247, 5.054635072878669, 5.529748307732377, 5.283587912365543, 5.003626677131027, 5.058319294131555, 5.084608316460255, 4.968580078828498, 5.110495908979962, 5.255359883584155, 5.162035024788262, 4.389101406050467, 5.115450162440548, 4.745798589733946, 5.000610552631832, 4.888987854673905, 5.176238243034651, 4.955969088226083, 5.069007842032237, 4.769351136942664, 5.120051081569637, 4.701095697280931, 4.428302989071524, 5.477083694828093, 5.185945916935429, 5.107727409294572, 5.02928645728968, 4.71270181596928, 5.02241797560756, 5.374772984658089, 5.01568441954427, 5.122235364490984, 4.970428717663235, 4.707290531852314, 4.847263022311885, 4.970399082212614, 5.003964329847607, 5.075538808977687, 5.12765920909692, 5.098014004639166, 5.374293166630107, 5.312073534985141, 4.922383871912405, 4.475868110724254, 4.961579428884218, 5.148898027409823, 4.946032965508302, 5.053474717089823, 5.172051002531645, 5.131720746177544, 5.137165258691042, 5.178841752343646, 5.481856614109793, 5.208297847713949, 5.012205777360402, 5.012151727685247, 4.926789913669775, 4.92408094633354, 4.802293477932388, 5.184739083526072, 4.95682729497978, 4.759107213964757, 4.648630490018539, 4.888194865028303, 4.927556078693026, 4.951795345908745, 5.00107551146813, 4.929365796846735, 4.916883888405719, 5.138656499948044, 4.964212695583082, 5.433121436270493, 4.238624197419311, 4.843950833178054, 5.157440006586765, 4.868629764603242, 4.969583937021509, 5.001084701680412, 27.50099806104303, 31.65878856178328, 27.29195716928815, 29.12053854632098, 30.52765952562965, 33.57505028711945, 30.00000000000012, 29.96226248154375, 29.9445936281242, 30.15637074482526, 30.92478914392544, 28.31402540268407, 31.53684514790622, 28.38208857676883, 29.44161606534021, 28.71018125410702, 29.97488541561197, 29.19671952636205, 29.78171673493658, 30.12879410572699, 30.08063466123329, 31.21040503234276, 28.16822589516432, 32.57744184354358, 38.76109490149969, 30.77541467145607, 30.79056512569276, 30.29828585807098, 30.54260128577718, 31.04423194402666, 27.26612384789891, 29.25292346828954, 29.93682727944105, 26.8080702086794, 28.4072409388417, 25.79700194192017, 30.47414207298604, 28.57716528969637, 29.79221216461831, 30.85545216067885, 23.24064442512047, 30.15638953215357, 34.74513288600978, 31.20559292767926, 29.79693717960585, 35.41992808824129, 32.02976429532768, 31.56900960837012, 30.36492302091277, 28.4959629400539, 25.30901009127507, 31.84878627556879, 31.3807995846957, 31.52172489116739, 30.67193471234086, 19.7660869897124, 31.00839511549184, 27.592430960884, 29.54658074594938, 30.29421671798915, 31.89774936632108, 28.65359788394929, 29.944066270218, 28.94024482537492, 27.13191294598704, 17.03583128251178, 31.29552044744389, 30.30034317206152, 31.61659485190562, 30.26196840822653, 30.64895568720903, 28.09002567472689, 31.05592414810993, 29.52223802697058, 30.64691389322694, 21.98408814962, 32.69761971601398, 29.99706500014407, 29.78179246546126, 29.55790942177752, 32.70883452788847, 29.49635153723926, 28.07635090751703, 27.77952916106336, 24.783512861126, 24.56354241655591, 32.01719286688631, 30.49803558311307, 29.35015511383196, 31.60396270286034, 29.46835038870734, 30.79368965070502, 27.86722492120627, 27.31984876995413, 29.14386727883832, 32.58951499976944, 32.58747032881468, 28.80126157971769, 30.79591467461194, 30.003122119396);
	//setParams(m);
}

void LSystem::draw() {
	// 統計情報をクリア
	stats.density = cv::Mat_<int>::zeros(NUM_GRID, NUM_GRID);
	stats.coverage = cv::Mat_<float>::zeros(NUM_STAT_GRID, NUM_STAT_GRID);

	drawSegment(glm::mat4(), 0, string(1, axiom));

	// densityが0のセルについて、PMパラメータを真ん中の値にする
	// * densityが0の場合、PMパラメータが全く寄与しないため、
	// * indicatorからPMパラメータを予測することは不可能だから。
	for (int r = 0; r < stats.density.rows; ++r) {
		for (int c = 0; c < stats.density.cols; ++c) {
			if (stats.density(r, c) == 0) {
				float x = ((float)c + 0.5f) / NUM_GRID * GRID_SIZE - GRID_SIZE * 0.5f;
				float y = ((float)r + 0.5f) / NUM_GRID * GRID_SIZE;

				pair<int, int> uv = XYtoUV(x, y);

				deltas(uv.second, uv.first) = 45;
				levels(uv.second, uv.first) = 5;
				lengths(uv.second, uv.first) = 30;
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
	lengths = cv::Mat_<float>::zeros(NUM_GRID, NUM_GRID);
	for (int r = 0; r < deltas.rows; ++r) {
		for (int c = 0; c < deltas.cols; ++c) {
			deltas(r, c) = genRand(10, 80);
		}
	}
	for (int r = 0; r < levels.rows; ++r) {
		for (int c = 0; c < levels.cols; ++c) {
			levels(r, c) = genRand(3, 8);
		}
	}
	for (int r = 0; r < lengths.rows; ++r) {
		for (int c = 0; c < lengths.cols; ++c) {
			lengths(r, c) = genRand(10, 50);
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

	deltas = cv::Mat_<double>(NUM_GRID, NUM_GRID);
	levels = cv::Mat_<double>(NUM_GRID, NUM_GRID);
	lengths = cv::Mat_<double>(NUM_GRID, NUM_GRID);

	int index = 0;
	for (int r = 0; r < deltas.rows; ++r) {
		for (int c = 0; c < deltas.cols; ++c) {
			deltas(r, c) = m(index++, 0);
		}
	}
	for (int r = 0; r < levels.rows; ++r) {
		for (int c = 0; c < levels.cols; ++c) {
			levels(r, c) = (int)(m(index++, 0) + 0.5);	// 四捨五入させる！
		}
	}
	for (int r = 0; r < lengths.rows; ++r) {
		for (int c = 0; c < lengths.cols; ++c) {
			lengths(r, c) = m(index++, 0);
		}
	}

	// Hard constraintsに従って、値を修正する
	/*
	for (int r = 0; r < deltas.rows; ++r) {
		for (int c = 0; c < deltas.cols; ++c) {
			deltas(r, c) = glm::clamp(deltas(r, c), 10.0f, 80.0f);
		}
	}
	for (int r = 0; r < levels.rows; ++r) {
		for (int c = 0; c < levels.cols; ++c) {
			levels(r, c) = glm::clamp(levels(r, c), 3, 8);
		}
	}
	for (int r = 0; r < lengths.rows; ++r) {
		for (int c = 0; c < lengths.cols; ++c) {
			lengths(r, c) = glm::clamp(lengths(r, c), 10.0f, 50.0f);
		}
	}
	*/
}

/**
 * パラメータの配列を返却する。
 *
 * @return		パラメータの配列
 */
vector<float> LSystem::getParams() {
	vector<float> ret(deltas.rows * deltas.cols + levels.rows * levels.cols + lengths.rows * lengths.cols);

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
	for (int r = 0; r < lengths.rows; ++r) {
		for (int c = 0; c < lengths.cols; ++c) {
			ret[index++] = lengths(r, c);
		}
	}
	return ret;
}

vector<float> LSystem::getStatistics() {
	vector<float> ret(stats.coverage.rows * stats.coverage.cols);
	int index = 0;
	for (int r = 0; r < stats.coverage.rows; ++r) {
		for (int c = 0; c < stats.coverage.cols; ++c) {
			ret[index++] = stats.coverage(r, c);
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
					if (uv.first >= 0 && uv.first < NUM_GRID && uv.second >= 0 && uv.second < NUM_GRID) {
						float l = lengths(uv.second, uv.first);
						drawCylinder(modelMat, 1, 1, l, glm::vec3(1.0, 1.0, 1.0));
						modelMat = glm::translate(modelMat, glm::vec3(0, l, 0));
					}
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

		int u = ((p1.x + p2.x) * 0.5 + GRID_SIZE * 0.5) / (GRID_SIZE / NUM_GRID);
		int v = (p1.y + p2.y) * 0.5 / (GRID_SIZE / NUM_GRID);
		if (u >= 0 && u < NUM_GRID && v >= 0 && v < NUM_GRID) {
			stats.density(v, u)++;
		}
	}

	// 統計情報を更新
	{
		for (float h = 0.0; h <= height; h += 1.0) {
			glm::vec4 p(0, h, 0, 1);
			p = modelMat * p;
			int u = (p.x + GRID_SIZE * 0.5) / (GRID_SIZE / NUM_STAT_GRID);
			int v = p.y / (GRID_SIZE / NUM_STAT_GRID);
			if (u >= 0 && u < NUM_STAT_GRID && v >= 0 && v < NUM_STAT_GRID) {
				stats.coverage(v, u) = 1;
			}
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
	double cell_size = (double)GRID_SIZE / NUM_GRID;
	int u = (x + GRID_SIZE * 0.5) / cell_size;
	int v = y / cell_size;

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
