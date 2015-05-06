#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QtGui/QMainWindow>
#include "ui_MainWindow.h"
#include "GLWidget3D.h"

class MainWindow : public QMainWindow {
	Q_OBJECT

private:
	Ui::MainWindowClass ui;
	GLWidget3D* glWidget;

public:
	MainWindow(QWidget *parent = 0, Qt::WFlags flags = 0);
	~MainWindow();

	void sample(int N, cv::Mat_<double>& dataX, cv::Mat_<double>& dataY);
	void normalizeData(cv::Mat_<double>& data, cv::Mat_<double>& normalized_data, cv::Mat_<double>& mu, cv::Mat_<double>& maxVal);
	void addBias(cv::Mat_<double>& data);
	void split(const cv::Mat_<double>& data, float train_ratio, float test_ratio, cv::Mat_<double>& train_data, cv::Mat_<double>& test_data);

public slots:
	void onGenerateSamples();
	void onGenerateSampleFiles();
	void onLinearRegression();
};

#endif // MAINWINDOW_H
