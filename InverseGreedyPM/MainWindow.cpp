#include "MainWindow.h"
#include <QDir>
#include <fstream>
#include "MLUtils.h"
#include "GreedyLSystem.h"
#include <QFileDialog>

MainWindow::MainWindow(QWidget *parent, Qt::WFlags flags) : QMainWindow(parent, flags) {
	ui.setupUi(this);

	connect(ui.actionExit, SIGNAL(triggered()), this, SLOT(close()));

	connect(ui.actionGenerate, SIGNAL(triggered()), this, SLOT(onGenerate()));
	connect(ui.actionInverse, SIGNAL(triggered()), this, SLOT(onInverse()));

	glWidget = new GLWidget3D(this);
	setCentralWidget(glWidget);
}

MainWindow::~MainWindow() {
}

/**
 * N個のデータをサンプリングする。
 *
 * @param N			サンプリング数
 * @param dataX		PM parameter
 * @param dataY		high-level indicator
 */
void MainWindow::sample(int N, int num_grid, int num_stat_grid, cv::Mat_<double>& dataX, cv::Mat_<double>& dataY) {
	bool initialized = false;

	int seed_count = 0;
	int check_point_interval = N / 10;
	for (int iter = 0; iter < N; ++iter) {
		// show progress
		if (iter > 0 && iter % check_point_interval == 0) {
			cout << iter / check_point_interval * 10 << "% >>";
		}

		glWidget->lsystem.randomInit(num_grid, num_stat_grid, seed_count++);
		glWidget->lsystem.draw();

		cv::Mat_<double> params = glWidget->lsystem.getParams();
		cv::Mat_<double> statistics = glWidget->lsystem.getStatistics();

		if (!initialized) {
			dataX = cv::Mat_<double>(N, params.cols);
			dataY = cv::Mat_<double>(N, statistics.cols);
			initialized = true;
		}

		params.copyTo(dataX.row(iter));
		statistics.copyTo(dataY.row(iter));
	}
}

void MainWindow::onGenerate() {
	int num_grid = 5;
	int num_stat_grid = 5;
	int N = 10000;

	if (!QDir("samples").exists()) QDir().mkdir("samples");

	cout << "Generating " << N << " samples (" << num_grid << "x" << num_grid << ")..." << endl;
	
	cv::Mat_<double> dataX;
	cv::Mat_<double> dataY;
	sample(N, num_grid, num_stat_grid, dataX, dataY);

	QString filenameX = "samples/samplesX_" + QString::number(N) + "_" + QString::number(num_grid) + "_" + QString::number(num_stat_grid) + ".txt";
	QString filenameY = "samples/samplesY_" + QString::number(N) + "_" + QString::number(num_grid) + "_" + QString::number(num_stat_grid) + ".txt";
	ml::saveDataset(filenameX.toUtf8().data(), dataX);
	ml::saveDataset(filenameY.toUtf8().data(), dataY);

	cout << "... done." << endl;

	glWidget->update();
}

void MainWindow::onInverse() {
	cv::Mat_<double> dataX;
	cv::Mat_<double> dataY;
	ml::loadDataset("samples/samplesX_10000_5_5.txt", dataX);
	ml::loadDataset("samples/samplesY_10000_5_5.txt", dataY);

	int num_grid = 5;
	int num_stat_grid = 5;
	int N = dataX.rows;

	cv::Mat target_density = cv::imread("target_density.png", 0);
	target_density.convertTo(target_density, CV_64F, 1.0/255.0);
	cv::flip(target_density, target_density, 0);

	// 白黒を反転させる
	target_density = 1 - target_density;

	// num_stat_gridにリサイズする
	cv::resize(target_density, target_density, cv::Size(num_stat_grid, num_stat_grid));

	// target_densityを、1xDの行列にする
	target_density = target_density.reshape(1, 1);

	// target_densityに基づいて生成する
	glWidget->lsystem.generate(num_grid, num_stat_grid, target_density, dataX, dataY);


	glWidget->update();

}
