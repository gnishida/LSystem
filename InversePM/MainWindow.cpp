#include "MainWindow.h"
#include <QDir>
#include <fstream>
#include "MLUtils.h"
#include "LSystem.h"
#include "LinearRegression.h"
#include "LinearRegressionRegularization.h"
#include "NearestNeighborRegression.h"
#include "LocalLinearRegression.h"
#include "ClusteredLinearRegression.h"
#include "GenerateSamplesWidget.h"
#include "asmOpenCV.h"
#include <QFileDialog>

MainWindow::MainWindow(QWidget *parent, Qt::WFlags flags) : QMainWindow(parent, flags) {
	ui.setupUi(this);

	connect(ui.actionExit, SIGNAL(triggered()), this, SLOT(close()));

	connect(ui.actionGenerateSamples, SIGNAL(triggered()), this, SLOT(onGenerateSamples()));
	connect(ui.actionBaseline, SIGNAL(triggered()), this, SLOT(onBaseline()));
	connect(ui.actionLinearRegression, SIGNAL(triggered()), this, SLOT(onLinearRegression()));
	connect(ui.actionNearestNeighbor, SIGNAL(triggered()), this, SLOT(onNearestNeighbor()));
	connect(ui.actionLocalRegression, SIGNAL(triggered()), this, SLOT(onLocalRegression()));
	connect(ui.actionClusteredLR, SIGNAL(triggered()), this, SLOT(onClusteredLR()));
	connect(ui.actionMCMC, SIGNAL(triggered()), this, SLOT(onMCMC()));

	connect(ui.actionFindLinearity, SIGNAL(triggered()), this, SLOT(onFindLinearity()));
	connect(ui.actionAngleEffect, SIGNAL(triggered()), this, SLOT(onAngleEffect()));

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

void MainWindow::onGenerateSamples() {
	GenerateSamplesWidget dlg(this);
	if (dlg.exec() != QDialog::Accepted) {
		return;
	}

	int num_grid = dlg.ui.lineEditNumGrid->text().toInt();
	int num_stat_grid = dlg.ui.lineEditNumStatGrid->text().toInt();
	int N = dlg.ui.lineEditNumSamples->text().toInt();

	if (!QDir("samples").exists()) QDir().mkdir("samples");

	cout << "Generating " << N << " samples (" << num_grid << "x" << num_grid << ")..." << endl;
	
	cv::Mat_<double> dataX;
	cv::Mat_<double> dataY;
	sample(N, num_grid, num_stat_grid, dataX, dataY);

	QString filenameX = "samples/samplesX_" + QString::number(N) + "_" + QString::number(num_grid) + "_" + QString::number(num_stat_grid) + ".txt";
	QString filenameY = "samples/samplesY_" + QString::number(N) + "_" + QString::number(num_grid) + "_" + QString::number(num_stat_grid) + ".txt";
	ml::saveDataset(filenameX.toUtf8().data(), dataX);
	ml::saveDataset(filenameY.toUtf8().data(), dataY);

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

void MainWindow::onBaseline() {
	cv::Mat_<double> dataX;
	cv::Mat_<double> dataY;
	cv::Mat_<double> normalized_dataX, normalized_dataY;
	cv::Mat_<double> muX, muY;
	cv::Mat_<double> maxX, maxY;
	cv::Mat_<double> train_dataX, train_dataY, test_dataX, test_dataY;
	cv::Mat_<double> train_normalized_dataX, train_normalized_dataY, test_normalized_dataX, test_normalized_dataY;

	int num_grid;
	int num_stat_grid;

	// PM parameterデータを読み込む
	{
		QString filename = QFileDialog::getOpenFileName(this, tr("Open PM parameter dataset..."), "", tr("PM parameter dataset (*.txt)"));
		if (filename.isEmpty()) return;

		ml::loadDataset(filename.toUtf8().data(), dataX);
		num_grid = sqrt(dataX.cols / 3.0);
	}

	// high-level indicatorデータを読み込む
	{
		QString filename = QFileDialog::getOpenFileName(this, tr("Open high-level indicator dataset..."), "", tr("High-level indicator dataset (*.txt)"));
		if (filename.isEmpty()) return;

		ml::loadDataset(filename.toUtf8().data(), dataY);
		num_stat_grid = sqrtf(dataY.cols);
	}

	if (!QDir("samples").exists()) QDir().mkdir("samples");

	// データをnormalizeしてテスト
	ml::normalizeDataset(dataX, normalized_dataX, muX, maxX);
	ml::normalizeDataset(dataY, normalized_dataY, muY, maxY);

	float ratio = 0.9f;
	if (dataX.rows * (1.0 - ratio) > 200) {
		ratio = (float)(dataX.rows - 200) / dataX.rows;
	}
	ml::splitDataset(dataX, ratio, train_dataX, test_dataX);
	ml::splitDataset(dataY, ratio, train_dataY, test_dataY);
	ml::splitDataset(normalized_dataX, ratio, train_normalized_dataX, test_normalized_dataX);
	ml::splitDataset(normalized_dataY, ratio, train_normalized_dataY, test_normalized_dataY);
	
	// 平均値を使って、yを予想する
	cv::Mat_<double> normalized_x_hat;
	cv::reduce(train_normalized_dataX, normalized_x_hat, 0, CV_REDUCE_AVG);
	cv::Mat_<double> trueStats(test_normalized_dataY.rows, num_stat_grid * num_stat_grid);
	cv::Mat_<double> predStats(test_normalized_dataY.rows, num_stat_grid * num_stat_grid);
	for (int iter = 0; iter < test_normalized_dataY.rows; ++iter) {
		cv::Mat x_hat = normalized_x_hat.mul(maxX) + muX;

		glWidget->lsystem.setParams(num_grid, num_stat_grid, test_dataX.row(iter));
		glWidget->updateGL();
		//if (dlg.ui.checkBoxSaveImages->isChecked()) {
		{
			QString fileName = "samples/" + QString::number(iter) + ".png";
			glWidget->grabFrameBuffer().save(fileName);
		}
		cv::Mat_<double> stats = glWidget->lsystem.getStatistics();
		stats.copyTo(trueStats.row(iter));

		glWidget->lsystem.setParams(num_grid, num_stat_grid, x_hat);
		glWidget->updateGL();
		//if (dlg.ui.checkBoxSaveImages->isChecked()) {
		{
			QString fileName = "samples/reversed_" + QString::number(iter) + ".png";
			glWidget->grabFrameBuffer().save(fileName);
		}
		cv::Mat_<double> stats2 = glWidget->lsystem.getStatistics();
		stats2.copyTo(predStats.row(iter));			
	}

	cout << "Baseline RMSE: " << ml::rmse(trueStats, predStats, false) << " / (Per cell) " << ml::rmse(trueStats, predStats, true) << endl;
}

void MainWindow::onLinearRegression() {
	cv::Mat_<double> dataX;
	cv::Mat_<double> dataY;
	cv::Mat_<double> normalized_dataX, normalized_dataY;
	cv::Mat_<double> muX, muY;
	cv::Mat_<double> maxX, maxY;
	cv::Mat_<double> train_dataX, train_dataY, test_dataX, test_dataY;
	cv::Mat_<double> train_normalized_dataX, train_normalized_dataY, test_normalized_dataX, test_normalized_dataY;

	int num_grid;
	int num_stat_grid;

	// PM parameterデータを読み込む
	{
		QString filename = QFileDialog::getOpenFileName(this, tr("Open PM parameter dataset..."), "", tr("PM parameter dataset (*.txt)"));
		if (filename.isEmpty()) return;

		ml::loadDataset(filename.toUtf8().data(), dataX);
		num_grid = sqrt(dataX.cols / 3.0);
	}

	// high-level indicatorデータを読み込む
	{
		QString filename = QFileDialog::getOpenFileName(this, tr("Open high-level indicator dataset..."), "", tr("High-level indicator dataset (*.txt)"));
		if (filename.isEmpty()) return;

		ml::loadDataset(filename.toUtf8().data(), dataY);
		num_stat_grid = sqrtf(dataY.cols);
	}

	if (!QDir("samples").exists()) QDir().mkdir("samples");

	// densityをvisualize
	{
		cv::Mat_<double> sum_Y;
		cv::reduce(dataY, sum_Y, 0, CV_REDUCE_SUM);
		cv::Mat density_img(num_stat_grid, num_stat_grid, CV_8U);
		
		int di = 0;
		for (int r = 0; r < density_img.rows; ++r) {
			for (int c = 0; c < density_img.cols; ++c) {
				density_img.at<uchar>(r, c) = (int)min(255.0, sum_Y(0, di++));//255.0);
			}
		}
		flip(density_img, density_img, 0);
		imwrite("samples/density.png", density_img);
	}

	// データをnormalizeしてテスト
	ml::normalizeDataset(dataX, normalized_dataX, muX, maxX);
	ml::normalizeDataset(dataY, normalized_dataY, muY, maxY);

	float ratio = 0.9f;
	if (dataX.rows * (1.0 - ratio) > 200) {
		ratio = (float)(dataX.rows - 200) / dataX.rows;
	}

	ml::splitDataset(dataX, ratio, train_dataX, test_dataX);
	ml::splitDataset(dataY, ratio, train_dataY, test_dataY);
	ml::splitDataset(normalized_dataX, ratio, train_normalized_dataX, test_normalized_dataX);
	ml::splitDataset(normalized_dataY, ratio, train_normalized_dataY, test_normalized_dataY);

	// Linear regressionにより、Wを求める（yW = x より、W = y^+ x)
	LinearRegression lr;
	double residue = lr.train(train_normalized_dataY, train_normalized_dataX);
	cout << "residue: " << residue << endl;
	cout << "condition number: " << lr.conditionNumber() << endl;

	// reverseで木を生成する
	cv::Mat_<double> trueStats(test_normalized_dataY.rows, num_stat_grid * num_stat_grid);
	cv::Mat_<double> predStats(test_normalized_dataY.rows, num_stat_grid * num_stat_grid);
	for (int iter = 0; iter < test_normalized_dataY.rows; ++iter) {
		cv::Mat normalized_x_hat = lr.predict(test_normalized_dataY.row(iter));
		cv::Mat x_hat = normalized_x_hat.mul(maxX) + muX;

		glWidget->lsystem.setParams(num_grid, num_stat_grid, test_dataX.row(iter));
		glWidget->updateGL();

		//if (dlg.ui.checkBoxSaveImages->isChecked()) {
		{
			QString fileName = "samples/" + QString::number(iter) + ".png";
			glWidget->grabFrameBuffer().save(fileName);
		}
		cv::Mat_<double> stats = glWidget->lsystem.getStatistics();
		stats.copyTo(trueStats.row(iter));

		glWidget->lsystem.setParams(num_grid, num_stat_grid, x_hat);
		glWidget->updateGL();
		//if (dlg.ui.checkBoxSaveImages->isChecked()) {
		{
			QString fileName = "samples/reversed_" + QString::number(iter) + ".png";
			glWidget->grabFrameBuffer().save(fileName);
		}
		cv::Mat_<double> stats2 = glWidget->lsystem.getStatistics();
		stats2.copyTo(predStats.row(iter));
	}

	cout << "LR RMSE: " << ml::rmse(trueStats, predStats, false) << " / (Per cell) " << ml::rmse(trueStats, predStats, true) << endl;
}

void MainWindow::onNearestNeighbor() {
	cv::Mat_<double> dataX;
	cv::Mat_<double> dataY;
	cv::Mat_<double> train_dataX, train_dataY, test_dataX, test_dataY;

	int num_grid;
	int num_stat_grid;

	// PM parameterデータを読み込む
	{
		QString filename = QFileDialog::getOpenFileName(this, tr("Open PM parameter dataset..."), "", tr("PM parameter dataset (*.txt)"));
		if (filename.isEmpty()) return;

		ml::loadDataset(filename.toUtf8().data(), dataX);
		num_grid = sqrt(dataX.cols / 3.0);
	}

	// high-level indicatorデータを読み込む
	{
		QString filename = QFileDialog::getOpenFileName(this, tr("Open high-level indicator dataset..."), "", tr("High-level indicator dataset (*.txt)"));
		if (filename.isEmpty()) return;

		ml::loadDataset(filename.toUtf8().data(), dataY);
		num_stat_grid = sqrtf(dataY.cols);
	}

	if (!QDir("samples").exists()) QDir().mkdir("samples");

	float ratio = 0.9f;
	if (dataX.rows * (1.0 - ratio) > 200) {
		ratio = (float)(dataX.rows - 200) / dataX.rows;
	}
	ml::splitDataset(dataX, ratio, train_dataX, test_dataX);
	ml::splitDataset(dataY, ratio, train_dataY, test_dataY);

	// Nearest neighborを使って、yを予想する
	NearestNeighborRegression nnr;
	cv::Mat_<double> trueStats(test_dataY.rows, num_stat_grid * num_stat_grid);
	cv::Mat_<double> predStats(test_dataY.rows, num_stat_grid * num_stat_grid);
	for (int iter = 0; iter < test_dataY.rows; ++iter) {
		double dist;
		cv::Mat x_hat = nnr.predict(train_dataY, train_dataX, test_dataY.row(iter), dist);
		
		glWidget->lsystem.setParams(num_grid, num_stat_grid, test_dataX.row(iter));
		glWidget->updateGL();
		//if (dlg.ui.checkBoxSaveImages->isChecked()) {
		{
			QString fileName = "samples/" + QString::number(iter) + ".png";
			glWidget->grabFrameBuffer().save(fileName);
		}
		cv::Mat_<double> stats = glWidget->lsystem.getStatistics();
		stats.copyTo(trueStats.row(iter));

		glWidget->lsystem.setParams(num_grid, num_stat_grid, x_hat);
		glWidget->updateGL();
		//if (dlg.ui.checkBoxSaveImages->isChecked()) {
		{
			QString fileName = "samples/reversed_" + QString::number(iter) + ".png";
			glWidget->grabFrameBuffer().save(fileName);
		}
		cv::Mat_<double> stats2 = glWidget->lsystem.getStatistics();
		stats2.copyTo(predStats.row(iter));

	}

	cout << "NN RMSE: " << ml::rmse(trueStats, predStats, false) << " / (Per cell) " << ml::rmse(trueStats, predStats, true) << endl;
}

void MainWindow::onLocalRegression() {
	cv::Mat_<double> dataX;
	cv::Mat_<double> dataY;
	cv::Mat_<double> normalized_dataX, normalized_dataY;
	cv::Mat_<double> muX, muY;
	cv::Mat_<double> maxX, maxY;
	cv::Mat_<double> train_dataX, train_dataY, test_dataX, test_dataY;
	cv::Mat_<double> train_normalized_dataX, train_normalized_dataY, test_normalized_dataX, test_normalized_dataY;

	int num_grid;
	int num_stat_grid;

	// PM parameterデータを読み込む
	{
		QString filename = QFileDialog::getOpenFileName(this, tr("Open PM parameter dataset..."), "", tr("PM parameter dataset (*.txt)"));
		if (filename.isEmpty()) return;

		ml::loadDataset(filename.toUtf8().data(), dataX);
		num_grid = sqrt(dataX.cols / 3.0);
	}

	// high-level indicatorデータを読み込む
	{
		QString filename = QFileDialog::getOpenFileName(this, tr("Open high-level indicator dataset..."), "", tr("High-level indicator dataset (*.txt)"));
		if (filename.isEmpty()) return;

		ml::loadDataset(filename.toUtf8().data(), dataY);
		num_stat_grid = sqrtf(dataY.cols);
	}

	if (!QDir("samples").exists()) QDir().mkdir("samples");

	// データをnormalizeしてテスト
	ml::normalizeDataset(dataX, normalized_dataX, muX, maxX);
	ml::normalizeDataset(dataY, normalized_dataY, muY, maxY);

	float ratio = 0.9f;
	if (dataX.rows * (1.0 - ratio) > 200) {
		ratio = (float)(dataX.rows - 200) / dataX.rows;
	}
	ml::splitDataset(dataX, ratio, train_dataX, test_dataX);
	ml::splitDataset(dataY, ratio, train_dataY, test_dataY);
	ml::splitDataset(normalized_dataX, ratio, train_normalized_dataX, test_normalized_dataX);
	ml::splitDataset(normalized_dataY, ratio, train_normalized_dataY, test_normalized_dataY);

	// Local regressionを使って、yを予想する
	LocalLinearRegression llr;
	cv::Mat_<double> trueStats(test_normalized_dataY.rows, num_stat_grid * num_stat_grid);
	cv::Mat_<double> predStats(test_normalized_dataY.rows, num_stat_grid * num_stat_grid);
	for (int iter = 0; iter < test_normalized_dataY.rows; ++iter) {
		double dist;
		cv::Mat normalized_x_hat = llr.predict(train_normalized_dataY, train_normalized_dataX, test_normalized_dataY.row(iter), 30);

		cv::Mat x_hat = normalized_x_hat.mul(maxX) + muX;

		glWidget->lsystem.setParams(num_grid, num_stat_grid, test_dataX.row(iter));
		glWidget->updateGL();
		//if (dlg.ui.checkBoxSaveImages->isChecked()) {
		{
			QString fileName = "samples/" + QString::number(iter) + ".png";
			glWidget->grabFrameBuffer().save(fileName);
		}
		cv::Mat_<double> stats = glWidget->lsystem.getStatistics();
		stats.copyTo(trueStats.row(iter));

		glWidget->lsystem.setParams(num_grid, num_stat_grid, x_hat);
		glWidget->updateGL();
		//if (dlg.ui.checkBoxSaveImages->isChecked()) {
		{
			QString fileName = "samples/reversed_" + QString::number(iter) + ".png";
			glWidget->grabFrameBuffer().save(fileName);
		}
		cv::Mat_<double> stats2 = glWidget->lsystem.getStatistics();
		stats2.copyTo(predStats.row(iter));		
	}

	cout << "LLR RMSE: " << ml::rmse(trueStats, predStats, false) << " / (Per cell) " << ml::rmse(trueStats, predStats, true) << endl;
}

void MainWindow::onClusteredLR() {
	cv::Mat_<double> dataX;
	cv::Mat_<double> dataY;
	cv::Mat_<double> normalized_dataX, normalized_dataY;
	cv::Mat_<double> muX, muY;
	cv::Mat_<double> maxX, maxY;
	cv::Mat_<double> train_dataX, train_dataY, test_dataX, test_dataY;
	cv::Mat_<double> train_normalized_dataX, train_normalized_dataY, test_normalized_dataX, test_normalized_dataY;

	int num_grid;
	int num_stat_grid;

	// PM parameterデータを読み込む
	{
		QString filename = QFileDialog::getOpenFileName(this, tr("Open PM parameter dataset..."), "", tr("PM parameter dataset (*.txt)"));
		if (filename.isEmpty()) return;

		ml::loadDataset(filename.toUtf8().data(), dataX);
		num_grid = sqrt(dataX.cols / 3.0);
	}

	// high-level indicatorデータを読み込む
	{
		QString filename = QFileDialog::getOpenFileName(this, tr("Open high-level indicator dataset..."), "", tr("High-level indicator dataset (*.txt)"));
		if (filename.isEmpty()) return;

		ml::loadDataset(filename.toUtf8().data(), dataY);
		num_stat_grid = sqrtf(dataY.cols);
	}

	if (!QDir("samples").exists()) QDir().mkdir("samples");

	// データをnormalizeしてテスト
	ml::normalizeDataset(dataX, normalized_dataX, muX, maxX);
	ml::normalizeDataset(dataY, normalized_dataY, muY, maxY);

	float ratio = 0.9f;
	if (dataX.rows * (1.0 - ratio) > 200) {
		ratio = (float)(dataX.rows - 200) / dataX.rows;
	}
	ml::splitDataset(dataX, ratio, train_dataX, test_dataX);
	ml::splitDataset(dataY, ratio, train_dataY, test_dataY);
	ml::splitDataset(normalized_dataX, ratio, train_normalized_dataX, test_normalized_dataX);
	ml::splitDataset(normalized_dataY, ratio, train_normalized_dataY, test_normalized_dataY);

	// clustered linear regressionを使って、yを予想する
	ClusteredLinearRegression clr(train_normalized_dataY, train_normalized_dataX, train_normalized_dataY.cols * 2);
	cv::Mat_<double> trueStats(test_normalized_dataY.rows, num_stat_grid * num_stat_grid);
	cv::Mat_<double> predStats(test_normalized_dataY.rows, num_stat_grid * num_stat_grid);
	for (int iter = 0; iter < test_normalized_dataY.rows; ++iter) {
		double dist;
		cv::Mat normalized_x_hat = clr.predict(test_normalized_dataY.row(iter));

		cv::Mat x_hat = normalized_x_hat.mul(maxX) + muX;

		glWidget->lsystem.setParams(num_grid, num_stat_grid, test_dataX.row(iter));
		glWidget->updateGL();
		//if (dlg.ui.checkBoxSaveImages->isChecked()) {
		{
			QString fileName = "samples/" + QString::number(iter) + ".png";
			glWidget->grabFrameBuffer().save(fileName);
		}
		cv::Mat_<double> stats = glWidget->lsystem.getStatistics();
		stats.copyTo(trueStats.row(iter));

		glWidget->lsystem.setParams(num_grid, num_stat_grid, x_hat);
		glWidget->updateGL();
		//if (dlg.ui.checkBoxSaveImages->isChecked()) {
		{
			QString fileName = "samples/reversed_" + QString::number(iter) + ".png";
			glWidget->grabFrameBuffer().save(fileName);
		}
		cv::Mat_<double> stats2 = glWidget->lsystem.getStatistics();
		stats2.copyTo(predStats.row(iter));		
	}

	cout << "CLR RMSE: " << ml::rmse(trueStats, predStats, false) << " / (Per cell) " << ml::rmse(trueStats, predStats, true) << endl;
}

void MainWindow::onMCMC() {
}

void MainWindow::onFindLinearity() {
	cv::Mat_<double> dataX;
	cv::Mat_<double> dataY;
	cv::Mat_<double> normalized_dataX, normalized_dataY;
	cv::Mat_<double> muX, muY;
	cv::Mat_<double> maxX, maxY;
	cv::Mat_<double> train_dataX, train_dataY, test_dataX, test_dataY;
	cv::Mat_<double> train_normalized_dataX, train_normalized_dataY, test_normalized_dataX, test_normalized_dataY;

	int num_grid;
	int num_stat_grid;

	// PM parameterデータを読み込む
	{
		QString filename = QFileDialog::getOpenFileName(this, tr("Open PM parameter dataset..."), "", tr("PM parameter dataset (*.txt)"));
		if (filename.isEmpty()) return;

		ml::loadDataset(filename.toUtf8().data(), dataX);
		num_grid = sqrt(dataX.cols / 3.0);
	}

	// high-level indicatorデータを読み込む
	{
		QString filename = QFileDialog::getOpenFileName(this, tr("Open high-level indicator dataset..."), "", tr("High-level indicator dataset (*.txt)"));
		if (filename.isEmpty()) return;

		ml::loadDataset(filename.toUtf8().data(), dataY);
		num_stat_grid = sqrtf(dataY.cols);
	}

	// データをnormalizeしてテスト
	ml::normalizeDataset(dataX, normalized_dataX, muX, maxX);
	ml::normalizeDataset(dataY, normalized_dataY, muY, maxY);

	LinearRegression lr;

	// PM parameterの各パラメータについて、線形性をチェック
	{
		ofstream ofs("pm_linearity.txt");
		cout << "For each PM parameter, check linearity..." << endl;
		double min_residue = std::numeric_limits<double>::max();
		double avg_residue = 0.0;
		for (int c = 0; c < normalized_dataX.cols; ++c) {
			double residue = lr.train(normalized_dataY, normalized_dataX.col(c));
			//cout << lr.predict(normalized_dataY.row(0)) << endl;
			//cout << normalized_dataX(0, 0) << endl;
			ofs << residue << endl;

			avg_residue += residue;
			if (residue < min_residue) {
				min_residue = residue;
			}
		}
		ofs.close();
		avg_residue /= normalized_dataX.cols;

		cout << "Min residue: " << min_residue << endl;
		cout << "Avg residue: " << avg_residue << endl;
	}

	// Indicatorの各要素について、線形性をチェック
	{
		ofstream ofs("indicator_linearity.txt");
		cout << "For each indicator, check linearity..." << endl;
		double min_residue = std::numeric_limits<double>::max();
		double avg_residue = 0.0;
		for (int c = 0; c < normalized_dataY.cols; ++c) {
			double residue = lr.train(normalized_dataX, normalized_dataY.col(c));
			//cout << lr.predict(normalized_dataX.row(0)) << endl;
			//cout << normalized_dataY(0, 0) << endl;
			ofs << residue << endl;

			avg_residue += residue;
			if (residue < min_residue) {
				min_residue = residue;
			}
		}
		ofs.close();
		avg_residue /= normalized_dataY.cols;

		cout << "Min residue: " << min_residue << endl;
		cout << "Avg residue: " << avg_residue << endl;
	}
}

void MainWindow::onAngleEffect() {
	int num_grid = 5;
	int num_stat_grid = 5;

	glWidget->lsystem.randomInit(num_grid, num_stat_grid, 18);
	glWidget->updateGL();
	cv::Mat_<double> params = glWidget->lsystem.getParams();

	ofstream ofs("angle_effect.txt");

	for (double angle = 10.0; angle <= 80.0; angle += 0.1) {
		ofs << angle;
		// 中心のセルのangleパラメータだけを変更する
		params(0, num_grid * 2.5) = angle;

		// 真ん中、下から2番目のangleパラメータだけを変更する
		//params(0, num_grid * 1.5) = angle;

		glWidget->lsystem.setParams(num_grid, num_stat_grid, params);
		glWidget->updateGL();
		cv::Mat_<double> stats = glWidget->lsystem.getStatistics();

		for (int c = 0; c < stats.cols; ++c) {
			ofs << " " << stats(0, c);
		}
		ofs << endl;
	}
	ofs.close();

	// 最後に、平均値での画像を表示する
	{
		// 中心のセルのangleパラメータだけを変更する
		params(0, num_grid * 2.5) = 45;

		// 真ん中、下から2番目のangleパラメータだけを変更する
		//params(0, num_grid * 1.5) = 45;

		glWidget->lsystem.setParams(num_grid, num_stat_grid, params);
		glWidget->updateGL();
	}
}
