﻿#include "MainWindow.h"
#include <QDir>
#include <fstream>
#include "MLUtils.h"
#include "LinearRegression.h"

MainWindow::MainWindow(QWidget *parent, Qt::WFlags flags) : QMainWindow(parent, flags) {
	ui.setupUi(this);

	connect(ui.actionExit, SIGNAL(triggered()), this, SLOT(close()));
	connect(ui.actionGenerateSamples, SIGNAL(triggered()), this, SLOT(onGenerateSamples()));
	connect(ui.actionGenerateSampleFiles, SIGNAL(triggered()), this, SLOT(onGenerateSampleFiles()));
	connect(ui.actionLinearRegression, SIGNAL(triggered()), this, SLOT(onLinearRegression()));

	glWidget = new GLWidget3D(this);
	setCentralWidget(glWidget);
}

MainWindow::~MainWindow() {
}

void MainWindow::sample(int N, cv::Mat_<double>& dataX, cv::Mat_<double>& dataY) {
	dataX = cv::Mat_<double>(N, lsystem::LSystem::NUM_GRID * lsystem::LSystem::NUM_GRID * 3);
	dataY = cv::Mat_<double>(N, lsystem::LSystem::NUM_GRID * lsystem::LSystem::NUM_GRID);

	int seed_count = 0;
	for (int iter = 0; iter < N; ++iter) {
		cout << iter << endl;

		glWidget->lsystem.randomInit(seed_count++);
		glWidget->lsystem.draw();

		vector<float> params = glWidget->lsystem.getParams();
		for (int col = 0; col < dataX.cols; ++col) {
			dataX(iter, col) = params[col];
		}

		vector<float> statistics = glWidget->lsystem.getStatistics();
		for (int col = 0; col < dataY.cols; ++col) {
			dataY(iter, col) = statistics[col];
		}
	}
}

void MainWindow::onGenerateSamples() {
	const int N = 100;

	if (!QDir("samples").exists()) QDir().mkdir("samples");

	cout << "Generating samples..." << endl;

	cv::Mat_<double> dataX;
	cv::Mat_<double> dataY;
	cv::Mat_<double> normalized_dataX;
	cv::Mat_<double> normalized_dataY;
	cv::Mat_<double> muX, muY;
	cv::Mat_<double> maxX, maxY;
	sample(N, dataX, dataY);

	for (int iter = 0; iter < N; ++iter) {
		glWidget->lsystem.setParams(dataX.row(iter));
		glWidget->updateGL();
		QString fileName = "samples/" + QString::number(iter) + ".png";
		glWidget->grabFrameBuffer().save(fileName);
	}

	glWidget->update();
}

void MainWindow::onGenerateSampleFiles() {
	const int N = 2000;

	if (!QDir("samples").exists()) QDir().mkdir("samples");

	cout << "Generating samples..." << endl;
	
	cv::Mat_<double> dataX;
	cv::Mat_<double> dataY;
	
	sample(N, dataX, dataY);

	ml::saveDataset("samples/samples.txt", dataY, dataX);

	glWidget->update();
}

void MainWindow::onLinearRegression() {
	const int N = 2000;

	if (!QDir("samples").exists()) QDir().mkdir("samples");

	cout << "Generating samples..." << endl;

	cv::Mat_<double> dataX;
	cv::Mat_<double> dataY;
	cv::Mat_<double> normalized_dataX, normalized_dataY;
	cv::Mat_<double> muX, muY;
	cv::Mat_<double> maxX, maxY;
	cv::Mat_<double> train_dataX, train_dataY, test_dataX, test_dataY;
	cv::Mat_<double> train_normalized_dataX, train_normalized_dataY, test_normalized_dataX, test_normalized_dataY;

	sample(N, dataX, dataY);
	ml::saveDataset("samples/samples.txt", dataY, dataX);

	ml::loadDataset("samples/samples.txt", dataY, dataX);

	ml::normalizeDataset(dataX, normalized_dataX, muX, maxX);
	ml::normalizeDataset(dataY, normalized_dataY, muY, maxY);
	ml::addBias(normalized_dataY);

	ml::splitDataset(dataX, 0.9, train_dataX, test_dataX);
	ml::splitDataset(dataY, 0.9, train_dataY, test_dataY);
	ml::splitDataset(normalized_dataX, 0.9, train_normalized_dataX, test_normalized_dataX);
	ml::splitDataset(normalized_dataY, 0.9, train_normalized_dataY, test_normalized_dataY);

	// Linear regressionにより、Wを求める（yW = x より、W = y^+ x)
	LinearRegression lr;
	lr.train(train_normalized_dataY, train_normalized_dataX);
	cout << "condition number: " << lr.conditionNumber() << endl;

	// reverseで木を生成する
	cv::Mat_<double> error = cv::Mat_<double>::zeros(1, dataX.cols);
	cv::Mat_<double> error2 = cv::Mat_<double>::zeros(1, dataX.cols);
	for (int iter = 0; iter < test_normalized_dataY.rows; ++iter) {
		cv::Mat normalized_x_hat = lr.predict(test_normalized_dataY.row(iter));
		cv::Mat x_hat = normalized_x_hat.mul(maxX) + muX;

		error += (test_normalized_dataX.row(iter) - normalized_x_hat).mul(test_normalized_dataX.row(iter) - normalized_x_hat);
		error2 += (test_dataX.row(iter) - x_hat).mul(test_dataX.row(iter) - x_hat);

		glWidget->lsystem.setParams(test_dataX.row(iter));
		glWidget->updateGL();
		QString fileName = "samples/" + QString::number(iter) + ".png";
		glWidget->grabFrameBuffer().save(fileName);

		glWidget->lsystem.setParams(x_hat);
		glWidget->updateGL();
		fileName = "samples/reversed_" + QString::number(iter) + ".png";
		glWidget->grabFrameBuffer().save(fileName);
	}

	error /= test_normalized_dataY.rows;
	error2 /= test_normalized_dataY.rows;
	cv::sqrt(error, error);
	cv::sqrt(error2, error2);
	cv::reduce(error, error, 1, CV_REDUCE_AVG);
	cv::reduce(error2, error2, 1, CV_REDUCE_AVG);

	cout << "Prediction error (normalized):" << endl;
	cout << error << endl;
	cout << "Prediction error:" << endl;
	cout << error2 << endl;
}