#include "ParametricLSystem.h"
#include <iostream>
#include "MLUtils.h"

int main() {
	parametriclsystem::ParametricLSystem pls(300, CV_32F, 0.1);
	for (int i = 0; i < 1000; ++i) {
	//for (int i = 0; i < 10; ++i) {
		cv::Mat indicator;
		parametriclsystem::String model = pls.derive(i);
		//cout << model << endl;

		// sampleの画像を保存する
		/*
		char filename[256];
		sprintf(filename, "images/sample%d.png", i);

		cv::Mat img;
		pls.computeIndicator(model, 1.0f, img);
		ml::mat_save(filename, img);
		*/
	}

	pls.drawTree(4);

	// ツリーの各ノードのindicatorをestimateする
	pls.gatherIndicators(0);
	//pls.saveIndicatorImages(4, 0.5);


	// ターゲットindicatorを読み込む
	cv::Mat target_density = cv::imread("target_indicator_30.png", 0);
	target_density.convertTo(target_density, CV_32F, 1.0/255.0);
	cv::flip(target_density, target_density, 0);

	// 白黒を反転させる
	target_density = 1 - target_density;

	// ターゲットに近いモデルを生成する
	cv::Mat indicator;
	parametriclsystem::String model = pls.inverse(target_density, 0.1, indicator);
	ml::mat_save("result.png", indicator);
	cout << "Resulting model:" << endl;
	cout << model << endl;

	// 生成したモデルの画像を保存する
	cv::Mat img;
	pls.computeIndicator(model, 1.0f, img);
	ml::mat_save("result2.png", img);

	return 0;
}