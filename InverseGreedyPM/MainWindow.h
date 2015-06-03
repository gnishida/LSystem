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

	int num_grid;
	int num_stat_grid;
	cv::Mat_<double> dataX;
	cv::Mat_<double> dataY;
	bool data_sampled;

public:
	MainWindow(int num_grid, int num_stat_grid, QWidget *parent = 0, Qt::WFlags flags = 0);
	~MainWindow();

	void sample(int N, int num_grid, int num_stat_grid, cv::Mat_<double>& dataX, cv::Mat_<double>& dataY);

public slots:
	void onGenerate();
	void onNearestSample();
	void onInverseFromNearestSample();
	void onInverseFromMedian();
	void onInverseFromRandom();
};

#endif // MAINWINDOW_H
