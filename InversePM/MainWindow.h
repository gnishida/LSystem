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

	void sample(int N, int num_grid, int num_stat_grid, cv::Mat_<double>& dataX, cv::Mat_<double>& dataY);

public slots:
	void onGenerateSamples();
	void onBaseline();
	void onLinearRegression();
	void onNearestNeighbor();
	void onLocalRegression();
	void onClusteredLR();
	void onMCMC();

	void onFindLinearity();
	void onAngleEffect();
};

#endif // MAINWINDOW_H
