#include "MainWindow.h"
#include <QDir>
#include <fstream>
#include "MLUtils.h"
#include "LSystem.h"
#include "LinearRegression.h"
#include "LinearRegressionRegularization.h"
#include "NearestNeighborRegression.h"
#include "LocalLinearRegression.h"
#include "GenerateSamplesWidget.h"
#include <QFileDialog>

MainWindow::MainWindow(QWidget *parent, Qt::WFlags flags) : QMainWindow(parent, flags) {
	ui.setupUi(this);

	connect(ui.actionExit, SIGNAL(triggered()), this, SLOT(close()));
	connect(ui.actionGenerateSamples, SIGNAL(triggered()), this, SLOT(onGenerateSamples()));
	connect(ui.actionLinearRegression, SIGNAL(triggered()), this, SLOT(onLinearRegression()));
	connect(ui.actionNearestNeighbor, SIGNAL(triggered()), this, SLOT(onNearestNeighbor()));
	connect(ui.actionLocalRegression, SIGNAL(triggered()), this, SLOT(onLocalRegression()));

	glWidget = new GLWidget3D(this);
	setCentralWidget(glWidget);
}

MainWindow::~MainWindow() {
}

void MainWindow::sample(int N, cv::Mat_<double>& dataX, cv::Mat_<double>& dataY) {
	dataX = cv::Mat_<double>(N, lsystem::LSystem::NUM_GRID * lsystem::LSystem::NUM_GRID * 3);
	dataY = cv::Mat_<double>(N, lsystem::LSystem::NUM_STAT_GRID * lsystem::LSystem::NUM_STAT_GRID);

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
	cv::Mat_<double> dataX;
	cv::Mat_<double> dataY;
	cv::Mat_<double> normalized_dataX, normalized_dataY;
	cv::Mat_<double> muX, muY;
	cv::Mat_<double> maxX, maxY;
	cv::Mat_<double> train_dataX, train_dataY, test_dataX, test_dataY;
	cv::Mat_<double> train_normalized_dataX, train_normalized_dataY, test_normalized_dataX, test_normalized_dataY;

	// PM parameterデータを読み込む
	{
		QString filename = QFileDialog::getOpenFileName(this, tr("Open PM parameter dataset..."), "", tr("PM parameter dataset (*.txt)"));
		if (filename.isEmpty()) return;

		ml::loadDataset(filename.toUtf8().data(), dataX);
	}

	// high-level indicatorデータを読み込む
	{
		QString filename = QFileDialog::getOpenFileName(this, tr("Open high-level indicator dataset..."), "", tr("High-level indicator dataset (*.txt)"));
		if (filename.isEmpty()) return;

		ml::loadDataset(filename.toUtf8().data(), dataY);
	}

	if (!QDir("samples").exists()) QDir().mkdir("samples");

	// densityをvisualize
	{
		cv::Mat_<double> sum_Y;
		cv::reduce(dataY, sum_Y, 0, CV_REDUCE_SUM);
		cv::Mat density_img(lsystem::LSystem::NUM_STAT_GRID, lsystem::LSystem::NUM_STAT_GRID, CV_8U);
		
		int di = 0;
		for (int r = 0; r < density_img.rows; ++r) {
			for (int c = 0; c < density_img.cols; ++c) {
				density_img.at<uchar>(r, c) = (int)min(255.0, sum_Y(0, di++));//255.0);
			}
		}
		flip(density_img, density_img, 0);
		imwrite("samples/density.png", density_img);
	}

#if 1
	// データをnormalizeしてテスト
	ml::normalizeDataset(dataX, normalized_dataX, muX, maxX);
	ml::normalizeDataset(dataY, normalized_dataY, muY, maxY);
	//ml::addBias(normalized_dataY);

	float ratio = 0.9f;
	if (dataX.rows * (1.0 - ratio) > 300) {
		ratio = (float)(dataX.rows - 300) / dataX.rows;
	}

	ml::splitDataset(dataX, ratio, train_dataX, test_dataX);
	ml::splitDataset(dataY, ratio, train_dataY, test_dataY);
	ml::splitDataset(normalized_dataX, ratio, train_normalized_dataX, test_normalized_dataX);
	ml::splitDataset(normalized_dataY, ratio, train_normalized_dataY, test_normalized_dataY);

	// Linear regressionにより、Wを求める（yW = x より、W = y^+ x)
	LinearRegression lr;
	double residue = lr.train(train_normalized_dataY, train_normalized_dataX);
	//LinearRegressionRegularization lr;
	//double residue = lr.train(train_normalized_dataY, train_normalized_dataX, 0, 0.01, 100);
	cout << "residue: " << residue << endl;
	cout << "condition number: " << lr.conditionNumber() << endl;

	// reverseで木を生成する
	//cv::Mat_<double> error = cv::Mat_<double>::zeros(1, test_normalized_dataY.cols - 1);
	cv::Mat_<double> trueStats(test_normalized_dataY.rows, test_normalized_dataY.cols);
	cv::Mat_<double> predStats(test_normalized_dataY.rows, test_normalized_dataY.cols);
	for (int iter = 0; iter < test_normalized_dataY.rows; ++iter) {
		cv::Mat normalized_x_hat = lr.predict(test_normalized_dataY.row(iter));
		cv::Mat x_hat = normalized_x_hat.mul(maxX) + muX;

		glWidget->lsystem.setParams(test_dataX.row(iter));
		glWidget->updateGL();

		//if (dlg.ui.checkBoxSaveImages->isChecked()) {
		{
			QString fileName = "samples/" + QString::number(iter) + ".png";
			glWidget->grabFrameBuffer().save(fileName);
		}
		vector<float> stats = glWidget->lsystem.getStatistics();
		cv::Mat_<double> stats_m = cv::Mat_<double>(cv::Mat(stats)).t();
		stats_m.copyTo(trueStats.row(iter));

		glWidget->lsystem.setParams(x_hat);
		glWidget->updateGL();
		//if (dlg.ui.checkBoxSaveImages->isChecked()) {
		{
			QString fileName = "samples/reversed_" + QString::number(iter) + ".png";
			glWidget->grabFrameBuffer().save(fileName);
		}
		vector<float> stats_hat = glWidget->lsystem.getStatistics();
		cv::Mat_<double> stats_hat_m = cv::Mat_<double>(cv::Mat(stats_hat)).t();
		stats_hat_m.copyTo(predStats.row(iter));	
	}

	cout << "RMSE: " << ml::rmse(trueStats, predStats) << endl;
#endif

#if 0
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
#endif
}

void MainWindow::onNearestNeighbor() {
	cv::Mat_<double> dataX;
	cv::Mat_<double> dataY;
	cv::Mat_<double> normalized_dataX, normalized_dataY;
	cv::Mat_<double> muX, muY;
	cv::Mat_<double> maxX, maxY;
	cv::Mat_<double> train_dataX, train_dataY, test_dataX, test_dataY;
	cv::Mat_<double> train_normalized_dataX, train_normalized_dataY, test_normalized_dataX, test_normalized_dataY;


	// PM parameterデータを読み込む
	{
		QString filename = QFileDialog::getOpenFileName(this, tr("Open PM parameter dataset..."), "", tr("PM parameter dataset (*.txt)"));
		if (filename.isEmpty()) return;

		ml::loadDataset(filename.toUtf8().data(), dataX);
	}

	// high-level indicatorデータを読み込む
	{
		QString filename = QFileDialog::getOpenFileName(this, tr("Open high-level indicator dataset..."), "", tr("High-level indicator dataset (*.txt)"));
		if (filename.isEmpty()) return;

		ml::loadDataset(filename.toUtf8().data(), dataY);
	}

	if (!QDir("samples").exists()) QDir().mkdir("samples");

	// データをnormalizeしてテスト
	ml::normalizeDataset(dataX, normalized_dataX, muX, maxX);
	ml::normalizeDataset(dataY, normalized_dataY, muY, maxY);
	//ml::addBias(normalized_dataY);

	float ratio = 0.9f;
	if (dataX.rows * (1.0 - ratio) > 300) {
		ratio = (float)(dataX.rows - 300) / dataX.rows;
	}
	ml::splitDataset(dataX, ratio, train_dataX, test_dataX);
	ml::splitDataset(dataY, ratio, train_dataY, test_dataY);
	ml::splitDataset(normalized_dataX, ratio, train_normalized_dataX, test_normalized_dataX);
	ml::splitDataset(normalized_dataY, ratio, train_normalized_dataY, test_normalized_dataY);

	// Nearest neighborを使って、yを予想する
	NearestNeighborRegression nnr;
	cv::Mat_<double> trueStats(test_normalized_dataY.rows, test_normalized_dataY.cols);
	cv::Mat_<double> predStats(test_normalized_dataY.rows, test_normalized_dataY.cols);
	for (int iter = 0; iter < test_normalized_dataY.rows; ++iter) {
		double dist;
		cv::Mat normalized_x_hat = nnr.predict(train_normalized_dataY, train_normalized_dataX, test_normalized_dataY.row(iter), dist);
		cv::Mat x_hat = normalized_x_hat.mul(maxX) + muX;

		glWidget->lsystem.setParams(test_dataX.row(iter));
		glWidget->updateGL();

		//if (dlg.ui.checkBoxSaveImages->isChecked()) {
		{
			QString fileName = "samples/" + QString::number(iter) + ".png";
			glWidget->grabFrameBuffer().save(fileName);
		}
		vector<float> stats = glWidget->lsystem.getStatistics();
		cv::Mat_<double> stats_m = cv::Mat_<double>(cv::Mat(stats)).t();
		stats_m.copyTo(trueStats.row(iter));

		glWidget->lsystem.setParams(x_hat);
		glWidget->updateGL();
		//if (dlg.ui.checkBoxSaveImages->isChecked()) {
		{
			QString fileName = "samples/reversed_" + QString::number(iter) + ".png";
			glWidget->grabFrameBuffer().save(fileName);
		}
		vector<float> stats_hat = glWidget->lsystem.getStatistics();
		cv::Mat_<double> stats_hat_m = cv::Mat_<double>(cv::Mat(stats_hat)).t();
		stats_hat_m.copyTo(predStats.row(iter));			
	}

	cout << "RMSE: " << ml::rmse(trueStats, predStats) << endl;
}

void MainWindow::onLocalRegression() {
	cv::Mat_<double> dataX;
	cv::Mat_<double> dataY;
	cv::Mat_<double> normalized_dataX, normalized_dataY;
	cv::Mat_<double> muX, muY;
	cv::Mat_<double> maxX, maxY;
	cv::Mat_<double> train_dataX, train_dataY, test_dataX, test_dataY;
	cv::Mat_<double> train_normalized_dataX, train_normalized_dataY, test_normalized_dataX, test_normalized_dataY;


	// PM parameterデータを読み込む
	{
		QString filename = QFileDialog::getOpenFileName(this, tr("Open PM parameter dataset..."), "", tr("PM parameter dataset (*.txt)"));
		if (filename.isEmpty()) return;

		ml::loadDataset(filename.toUtf8().data(), dataX);
	}

	// high-level indicatorデータを読み込む
	{
		QString filename = QFileDialog::getOpenFileName(this, tr("Open high-level indicator dataset..."), "", tr("High-level indicator dataset (*.txt)"));
		if (filename.isEmpty()) return;

		ml::loadDataset(filename.toUtf8().data(), dataY);
	}

	if (!QDir("samples").exists()) QDir().mkdir("samples");

	// データをnormalizeしてテスト
	ml::normalizeDataset(dataX, normalized_dataX, muX, maxX);
	ml::normalizeDataset(dataY, normalized_dataY, muY, maxY);
	//ml::addBias(normalized_dataY);

	float ratio = 0.9f;
	if (dataX.rows * (1.0 - ratio) > 300) {
		ratio = (float)(dataX.rows - 300) / dataX.rows;
	}
	ml::splitDataset(dataX, ratio, train_dataX, test_dataX);
	ml::splitDataset(dataY, ratio, train_dataY, test_dataY);
	ml::splitDataset(normalized_dataX, ratio, train_normalized_dataX, test_normalized_dataX);
	ml::splitDataset(normalized_dataY, ratio, train_normalized_dataY, test_normalized_dataY);

	// Local regressionを使って、yを予想する
	LocalLinearRegression llr;
	cv::Mat_<double> trueStats(test_normalized_dataY.rows, test_normalized_dataY.cols);
	cv::Mat_<double> predStats(test_normalized_dataY.rows, test_normalized_dataY.cols);
	for (int iter = 0; iter < test_normalized_dataY.rows; ++iter) {
		double dist;
		cv::Mat normalized_x_hat = llr.predict(train_normalized_dataY, train_normalized_dataX, test_normalized_dataY.row(iter), 30.0);

		cv::Mat x_hat = normalized_x_hat.mul(maxX) + muX;

		glWidget->lsystem.setParams(test_dataX.row(iter));
		glWidget->updateGL();

		//if (dlg.ui.checkBoxSaveImages->isChecked()) {
		{
			QString fileName = "samples/" + QString::number(iter) + ".png";
			glWidget->grabFrameBuffer().save(fileName);
		}
		vector<float> stats = glWidget->lsystem.getStatistics();
		cv::Mat_<double> stats_m = cv::Mat_<double>(cv::Mat(stats)).t();
		stats_m.copyTo(trueStats.row(iter));

		glWidget->lsystem.setParams(x_hat);
		glWidget->updateGL();
		//if (dlg.ui.checkBoxSaveImages->isChecked()) {
		{
			QString fileName = "samples/reversed_" + QString::number(iter) + ".png";
			glWidget->grabFrameBuffer().save(fileName);
		}
		vector<float> stats_hat = glWidget->lsystem.getStatistics();
		cv::Mat_<double> stats_hat_m = cv::Mat_<double>(cv::Mat(stats_hat)).t();
		stats_hat_m.copyTo(predStats.row(iter));			
	}

	cout << "RMSE: " << ml::rmse(trueStats, predStats) << endl;
}
