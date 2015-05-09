#include "MainWindow.h"
#include <QDir>
#include <fstream>
#include "MLUtils.h"
#include "LSystem.h"
#include "LinearRegression.h"
#include "LinearRegressionRegularization.h"
#include "GenerateSamplesWidget.h"

MainWindow::MainWindow(QWidget *parent, Qt::WFlags flags) : QMainWindow(parent, flags) {
	ui.setupUi(this);

	connect(ui.actionExit, SIGNAL(triggered()), this, SLOT(close()));
	connect(ui.actionGenerateSamples, SIGNAL(triggered()), this, SLOT(onGenerateSamples()));
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
	GenerateSamplesWidget dlg(this);
	if (dlg.exec() != QDialog::Accepted) {
		return;
	}

	lsystem::LSystem::NUM_GRID = dlg.ui.lineEditNumGrid->text().toInt();
	int N = dlg.ui.lineEditNumSamples->text().toInt();

	if (!QDir("samples").exists()) QDir().mkdir("samples");

	cout << "Generating samples..." << endl;
	
	cv::Mat_<double> dataX;
	cv::Mat_<double> dataY;
	sample(N, dataX, dataY);

	ml::saveDataset("samples/samplesX.txt", dataY);
	ml::saveDataset("samples/samplesY.txt", dataX);

	if (dlg.ui.checkBoxSaveImages->isChecked()) {
		for (int iter = 0; iter < N; ++iter) {
			glWidget->lsystem.setParams(dataX.row(iter));
			glWidget->updateGL();
			QString fileName = "samples/" + QString::number(iter) + ".png";
			glWidget->grabFrameBuffer().save(fileName);
		}
	}

	glWidget->update();
}

void MainWindow::onLinearRegression() {
	GenerateSamplesWidget dlg(this);
	if (dlg.exec() != QDialog::Accepted) {
		return;
	}

	lsystem::LSystem::NUM_GRID = dlg.ui.lineEditNumGrid->text().toInt();
	int N = dlg.ui.lineEditNumSamples->text().toInt();

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
	ml::saveDataset("samples/samplesX.txt", dataY);
	ml::saveDataset("samples/samplesY.txt", dataX);

	ml::loadDataset("samples/samplesX.txt", dataY);
	ml::loadDataset("samples/samplesY.txt", dataX);

	// densityをvisualize
	{
		cv::Mat_<double> sum_Y;
		cv::reduce(dataY, sum_Y, 0, CV_REDUCE_SUM);
		cv::Mat density_img(lsystem::LSystem::NUM_GRID, lsystem::LSystem::NUM_GRID, CV_8U);
		cout << sum_Y << endl;
		int di = 0;
		for (int r = 0; r < density_img.rows; ++r) {
			for (int c = 0; c < density_img.cols; ++c) {
				density_img.at<uchar>(r, c) = (int)min(255.0, sum_Y(0, di++));//255.0);
			}
		}
		flip(density_img, density_img, 0);
		imwrite("samples/density.png", density_img);
	}

	if (dlg.ui.checkBoxNormalizeData->isChecked()) {
		// データをnormalizeしてテスト
		ml::normalizeDataset(dataX, normalized_dataX, muX, maxX);
		ml::normalizeDataset(dataY, normalized_dataY, muY, maxY);
		ml::addBias(normalized_dataY);

		ml::splitDataset(dataX, 0.9, train_dataX, test_dataX);
		ml::splitDataset(dataY, 0.9, train_dataY, test_dataY);
		ml::splitDataset(normalized_dataX, 0.9, train_normalized_dataX, test_normalized_dataX);
		ml::splitDataset(normalized_dataY, 0.9, train_normalized_dataY, test_normalized_dataY);


		test_dataX = train_dataX.clone();
		test_dataY = train_dataY.clone();
		test_normalized_dataX = train_normalized_dataX.clone();
		test_normalized_dataY = train_normalized_dataY.clone();





		// Linear regressionにより、Wを求める（yW = x より、W = y^+ x)
		LinearRegression lr;
		double residue = lr.train(train_normalized_dataY, train_normalized_dataX);
		//LinearRegressionRegularization lr;
		//double residue = lr.train(train_normalized_dataY, train_normalized_dataX, 0, 0.01, 100);
		cout << "residue: " << residue << endl;
		cout << "condition number: " << lr.conditionNumber() << endl;

		// reverseで木を生成する
		cv::Mat_<double> error = cv::Mat_<double>::zeros(1, test_normalized_dataY.cols - 1);
		for (int iter = 0; iter < test_normalized_dataY.rows; ++iter) {
			cv::Mat normalized_x_hat = lr.predict(test_normalized_dataY.row(iter));
			cv::Mat x_hat = normalized_x_hat.mul(maxX) + muX;

			glWidget->lsystem.setParams(test_dataX.row(iter));
			glWidget->updateGL();

			if (dlg.ui.checkBoxSaveImages->isChecked()) {
				QString fileName = "samples/" + QString::number(iter) + ".png";
				glWidget->grabFrameBuffer().save(fileName);
			}
			vector<float> stats_hat = glWidget->lsystem.getStatistics();

			glWidget->lsystem.setParams(x_hat);
			glWidget->updateGL();
			if (dlg.ui.checkBoxSaveImages->isChecked()) {
				QString fileName = "samples/reversed_" + QString::number(iter) + ".png";
				glWidget->grabFrameBuffer().save(fileName);
			}
			vector<float> stats = glWidget->lsystem.getStatistics();
			
			// compute error
			error += ml::mat_square(cv::Mat_<double>(cv::Mat(stats)) - cv::Mat_<double>(cv::Mat(stats_hat))).t();
		}

		error /= (double)test_normalized_dataY.rows;
		cv::sqrt(error, error);

		cout << "Prediction error (normalized):" << endl;
		cout << error << endl;

		cv::reduce(error, error, 1, CV_REDUCE_SUM);
		cout << "Total error: " << error(0, 0) << endl;
	} else {
		// データをnormalizeしないでテスト
		ml::normalizeDataset(dataY, normalized_dataY, muY, maxY);

		ml::addBias(dataY);

		ml::splitDataset(dataX, 0.9, train_dataX, test_dataX);
		ml::splitDataset(dataY, 0.9, train_dataY, test_dataY);

		// Linear regressionにより、Wを求める（yW = x より、W = y^+ x)
		LinearRegression lr;
		lr.train(train_dataY, train_dataX);
		cout << "condition number: " << lr.conditionNumber() << endl;

		// reverseで木を生成する
		cv::Mat_<double> error = cv::Mat_<double>::zeros(1, test_dataY.cols - 1);
		for (int iter = 0; iter < test_dataY.rows; ++iter) {
			cv::Mat x_hat = lr.predict(test_dataY.row(iter));

			glWidget->lsystem.setParams(test_dataX.row(iter));
			glWidget->updateGL();
			if (dlg.ui.checkBoxSaveImages->isChecked()) {
				QString fileName = "samples/" + QString::number(iter) + ".png";
				glWidget->grabFrameBuffer().save(fileName);
			}
			vector<float> stats_hat = glWidget->lsystem.getStatistics();

			glWidget->lsystem.setParams(x_hat);
			glWidget->updateGL();
			if (dlg.ui.checkBoxSaveImages->isChecked()) {
				QString fileName = "samples/reversed_" + QString::number(iter) + ".png";
				glWidget->grabFrameBuffer().save(fileName);
			}
			vector<float> stats = glWidget->lsystem.getStatistics();
			
			// compute error		
			error += ml::mat_square(cv::Mat_<double>(cv::Mat(stats)) - cv::Mat_<double>(cv::Mat(stats_hat))).t();
		}

		error /= test_dataY.rows;
		cv::sqrt(error, error);

		cout << "Prediction error (normalized):" << endl;
		cout << error << endl;
		
		cv::reduce(error, error, 1, CV_REDUCE_SUM);
		cout << "Total error: " << error(0, 0) << endl;
	}
}
