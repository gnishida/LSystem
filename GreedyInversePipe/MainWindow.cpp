#include "MainWindow.h"

MainWindow::MainWindow(QWidget *parent, Qt::WFlags flags) : QMainWindow(parent, flags) {
	ui.setupUi(this);

	connect(ui.actionExit, SIGNAL(triggered()), this, SLOT(close()));
	connect(ui.actionRandomGeneration, SIGNAL(triggered()), this, SLOT(onRandomGeneration()));
	connect(ui.actionGreedyInverse, SIGNAL(triggered()), this, SLOT(onGreedyInverse()));

	glWidget = new GLWidget3D(this);
	setCentralWidget(glWidget);
}

MainWindow::~MainWindow() {
}

void MainWindow::onRandomGeneration() {
	glWidget->model = glWidget->lsystem.derive(0);

	cout << glWidget->model << endl;

	glWidget->updateGL();
}

void MainWindow::onGreedyInverse() {
	// ターゲットindicatorを読み込む
	cv::Mat target_indicator = cv::imread("target_indicator.png", 0);
	target_indicator.convertTo(target_indicator, CV_8U, 1.0/255.0);
	cv::flip(target_indicator, target_indicator, 0);

	// 白黒を反転させる
	target_indicator = 1 - target_indicator;

	// ターゲットに近いモデルを生成する
	cv::Mat indicator;
	glWidget->model = glWidget->lsystem.inverse(target_indicator, indicator);

	cout << glWidget->model << endl;

	glWidget->updateGL();
}