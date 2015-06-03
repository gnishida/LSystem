#include "MainWindow.h"
#include <QDir>
#include <fstream>
#include "MLUtils.h"
#include "GreedyLSystem.h"
#include "GenerateSamplesWidget.h"
#include <QFileDialog>

MainWindow::MainWindow(int num_grid, int num_stat_grid, QWidget *parent, Qt::WFlags flags) : QMainWindow(parent, flags) {
	this->num_grid = num_grid;
	this->num_stat_grid = num_stat_grid;

	ui.setupUi(this);

	connect(ui.actionExit, SIGNAL(triggered()), this, SLOT(close()));

	connect(ui.actionGenerate, SIGNAL(triggered()), this, SLOT(onGenerate()));
	connect(ui.actionNearestSample, SIGNAL(triggered()), this, SLOT(onNearestSample()));
	connect(ui.actionInverseFromNearestSample, SIGNAL(triggered()), this, SLOT(onInverseFromNearestSample()));
	connect(ui.actionInverseFromMedian, SIGNAL(triggered()), this, SLOT(onInverseFromMedian()));
	connect(ui.actionInverseFromRandom, SIGNAL(triggered()), this, SLOT(onInverseFromRandom()));

	glWidget = new GLWidget3D(this);
	setCentralWidget(glWidget);

	data_sampled = false;
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
	int check_point_interval = N / 10 + 1;
	for (int iter = 0; iter < N; ++iter) {
		// show progress
		if (iter > 0 && iter % check_point_interval == 0) {
			cout << iter / check_point_interval * 10 << "% >>";
		}

		glWidget->lsystem.randomInit(num_grid, num_stat_grid, seed_count++);
		glWidget->lsystem.draw(true);

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
	GenerateSamplesWidget dlg(this);
	if (dlg.exec() != QDialog::Accepted) {
		return;
	}

	int N = dlg.ui.lineEditNumSamples->text().toInt();

	if (!QDir("samples").exists()) QDir().mkdir("samples");

	cout << "Generating " << N << " samples (" << num_grid << "x" << num_grid << ")..." << endl;
	
	sample(N, num_grid, num_stat_grid, dataX, dataY);
	data_sampled = true;

	QString filenameX = "samples/samplesX.bin";
	QString filenameY = "samples/samplesY.bin";
	ml::saveDataset(filenameX.toUtf8().data(), dataX, true);
	ml::saveDataset(filenameY.toUtf8().data(), dataY, true);

	if (dlg.ui.checkBoxSaveImages->isChecked()) {
		for (int iter = 0; iter < N; ++iter) {
			glWidget->lsystem.setParams(num_grid, num_stat_grid, dataX.row(iter));
			glWidget->updateGL();
			QString fileName = "samples/" + QString::number(iter) + ".png";
			glWidget->grabFrameBuffer().save(fileName);
		}
	}

	cout << "... done." << endl;

	glWidget->update();
}

void MainWindow::onNearestSample() {
	// read target indicator
	QString filename = QFileDialog::getOpenFileName(this, tr("Open target indicator..."), "", tr("Indicator Files (*.png)"));
	if (filename.isEmpty()) return;
	
	cv::Mat target_density = cv::imread(filename.toUtf8().data(), 0);
	target_density.convertTo(target_density, CV_64F, 1.0/255.0);
	cv::flip(target_density, target_density, 0);

	if (!data_sampled) {
		ml::loadDataset("samples/samplesX.bin", dataX, true);
		ml::loadDataset("samples/samplesY.bin", dataY, true);
		data_sampled = true;
	}

	int N = dataX.rows;
	
	// 白黒を反転させる
	target_density = 1 - target_density;

	// num_stat_gridにリサイズする
	cv::resize(target_density, target_density, cv::Size(num_stat_grid, num_stat_grid));

	// target_densityを、1xDの行列にする
	target_density = target_density.reshape(1, 1);

	// target_densityに基づいて生成する
	glWidget->lsystem.nearestSample(num_grid, num_stat_grid, target_density, dataX, dataY);

	glWidget->update();
}

void MainWindow::onInverseFromNearestSample() {
	// read target indicator
	QString filename = QFileDialog::getOpenFileName(this, tr("Open target indicator..."), "", tr("Indicator Files (*.png)"));
	if (filename.isEmpty()) return;

	cv::Mat target_density = cv::imread(filename.toUtf8().data(), 0);
	target_density.convertTo(target_density, CV_64F, 1.0/255.0);
	cv::flip(target_density, target_density, 0);

	if (!data_sampled) {
		ml::loadDataset("samples/samplesX.bin", dataX, true);
		ml::loadDataset("samples/samplesY.bin", dataY, true);
		data_sampled = true;
	}

	int N = dataX.rows;

	// 白黒を反転させる
	target_density = 1 - target_density;

	// num_stat_gridにリサイズする
	cv::resize(target_density, target_density, cv::Size(num_stat_grid, num_stat_grid));

	// target_densityを、1xDの行列にする
	target_density = target_density.reshape(1, 1);

	// target_densityに基づいて生成する
	glWidget->lsystem.nearestSample(num_grid, num_stat_grid, target_density, dataX, dataY);
	glWidget->lsystem.inverse(target_density);

	glWidget->update();
}

void MainWindow::onInverseFromMedian() {
	// read target indicator
	QString filename = QFileDialog::getOpenFileName(this, tr("Open target indicator..."), "", tr("Indicator Files (*.png)"));
	if (filename.isEmpty()) return;

	cv::Mat target_density = cv::imread(filename.toUtf8().data(), 0);
	target_density.convertTo(target_density, CV_64F, 1.0/255.0);
	cv::flip(target_density, target_density, 0);

	if (!data_sampled) {
		ml::loadDataset("samples/samplesX.bin", dataX, true);
		ml::loadDataset("samples/samplesY.bin", dataY, true);
		data_sampled = true;
	}

	int N = dataX.rows;

	// 白黒を反転させる
	target_density = 1 - target_density;

	// num_stat_gridにリサイズする
	cv::resize(target_density, target_density, cv::Size(num_stat_grid, num_stat_grid));

	// target_densityを、1xDの行列にする
	target_density = target_density.reshape(1, 1);

	// target_densityに基づいて生成する
	glWidget->lsystem.meanInit(num_grid, num_stat_grid);
	glWidget->lsystem.inverse(target_density);

	glWidget->update();
}


void MainWindow::onInverseFromRandom() {
	// read target indicator
	QString filename = QFileDialog::getOpenFileName(this, tr("Open target indicator..."), "", tr("Indicator Files (*.png)"));
	if (filename.isEmpty()) return;

	cv::Mat target_density = cv::imread(filename.toUtf8().data(), 0);
	target_density.convertTo(target_density, CV_64F, 1.0/255.0);
	cv::flip(target_density, target_density, 0);

	if (!data_sampled) {
		ml::loadDataset("samples/samplesX.bin", dataX, true);
		ml::loadDataset("samples/samplesY.bin", dataY, true);
		data_sampled = true;
	}

	int N = dataX.rows;

	// 白黒を反転させる
	target_density = 1 - target_density;

	// num_stat_gridにリサイズする
	cv::resize(target_density, target_density, cv::Size(num_stat_grid, num_stat_grid));

	// target_densityを、1xDの行列にする
	target_density = target_density.reshape(1, 1);

	// target_densityに基づいて生成する
	glWidget->lsystem.randomInit(num_grid, num_stat_grid, 0);
	glWidget->lsystem.inverse(target_density);

	glWidget->update();
}
