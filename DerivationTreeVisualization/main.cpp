#include "ParametricLSystem.h"
#include <iostream>
#include "MLUtils.h"

int main(int argc, char *argv[]) {
	if (argc < 2) {
		cout << "Usage: " << argv[0] << " <# samples>" << endl;
		return -1;
	}

	int num_samples = atoi(argv[1]);

	parametriclsystem::ParametricLSystem pls(300, CV_32F, 0.1);
	pls.generateDerivationStateTree();
	/*
	for (int i = 0; i < 10; ++i) {
		cv::Mat indicator;
		parametriclsystem::String model = pls.derive(i);
		cout << model << endl;

		// sampleの画像を保存する
		char filename[256];
		sprintf(filename, "images/sample%d.png", i);

		cv::Mat img;
		pls.computeIndicator(model, 1.0f, img);
		ml::mat_save(filename, img);
	}
	*/
	pls.drawTree(4);

	// ツリーの各ノードのindicatorをestimateする
	//pls.gatherIndicators(0);

	return 0;
}