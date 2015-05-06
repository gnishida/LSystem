#include "MainWindow.h"
#include <QDir>
#include <fstream>

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
	dataX = cv::Mat_<double>(N, lsystem::LSystem::NUM_GRID * lsystem::LSystem::NUM_GRID * 2);
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

void MainWindow::normalizeData(cv::Mat_<double>& data, cv::Mat_<double>& normalized_data, cv::Mat_<double>& mu, cv::Mat_<double>& maxVal) {
	cv::reduce(data, mu, 0, CV_REDUCE_AVG);
	normalized_data = data - cv::repeat(mu, data.rows, 1);

	// [-1, 1]にする
	cv::reduce(cv::abs(normalized_data), maxVal, 0, CV_REDUCE_MAX);
	for (int r = 0; r < data.rows; ++r) {
		for (int c = 0; c < normalized_data.cols; ++c) {
			if (maxVal(0, c) != 0) {
				normalized_data(r, c) /= maxVal(0, c);
			}
		}
	}
}

/**
 * 一番右の列に1を追加する。
 */
void MainWindow::addBias(cv::Mat_<double>& data) {
	cv::Mat_<double> tmp = data.clone();
	data = cv::Mat_<double>(tmp.rows, tmp.cols + 1);
	for (int r = 0; r < tmp.rows; ++r) {
		for (int c = 0; c < tmp.cols; ++c) {
			data(r, c) = tmp(r, c);
		}
		data(r, tmp.cols) = 1;
	}
}

/**
 * データを、指定された比率に従い、trainingデータとtestデータに分割する。
 *
 */
void MainWindow::split(const cv::Mat_<double>& data, float train_ratio, float test_ratio, cv::Mat_<double>& train_data, cv::Mat_<double>& test_data) {
	int train_rows = data.rows * train_ratio;
	int test_rows = data.rows - train_rows;

	train_data = cv::Mat_<double>(train_rows, data.cols);
	test_data = cv::Mat_<double>(test_rows, data.cols);

	for (int r = 0; r < data.rows; ++r) {
		for (int c = 0; c < data.cols; ++c) {
			if (r < train_rows) {
				train_data(r, c) = data(r, c);
			} else {
				test_data(r - train_rows, c) = data(r, c);
			}
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
	
	ofstream ofs("samples/samples.txt");

	cv::Mat_<double> dataX;
	cv::Mat_<double> dataY;
	cv::Mat_<double> normalized_dataX;
	cv::Mat_<double> normalized_dataY;
	cv::Mat_<double> muX, muY;
	cv::Mat_<double> maxX, maxY;
	sample(N, dataX, dataY);

	for (int iter = 0; iter < N; ++iter) {
		ofs << "[";
		for (int c = 0; c < dataY.cols; ++c) {
			if (c > 0) ofs << ",";
			ofs << dataY(iter, c);
		}
		ofs << "],[";
		for (int c = 0; c < dataX.cols; ++c) {
			if (c > 0) ofs << ",";
			ofs << dataX(iter, c);
		}
		ofs << "]" << endl;
	}
	ofs.close();

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
	normalizeData(dataX, normalized_dataX, muX, maxX);
	normalizeData(dataY, normalized_dataY, muY, maxY);
	addBias(normalized_dataY);
	split(dataX, 0.9, 0.1, train_dataX, test_dataX);
	split(dataY, 0.9, 0.1, train_dataY, test_dataY);
	split(normalized_dataX, 0.9, 0.1, train_normalized_dataX, test_normalized_dataX);
	split(normalized_dataY, 0.9, 0.1, train_normalized_dataY, test_normalized_dataY);

	// Linear regressionにより、Wを求める（yW = x より、W = y^+ x)
	cv::Mat_<double> W = train_normalized_dataY.inv(cv::DECOMP_SVD) * train_normalized_dataX;

	// reverseで木を生成する
	cv::Mat_<double> error = cv::Mat_<double>::zeros(1, dataX.cols);
	cv::Mat_<double> error2 = cv::Mat_<double>::zeros(1, dataX.cols);
	for (int iter = 0; iter < test_normalized_dataY.rows; ++iter) {
		cv::Mat normalized_x_hat = test_normalized_dataY.row(iter) * W;
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
