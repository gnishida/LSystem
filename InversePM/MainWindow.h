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

public slots:
	void onGenerateSamples();
	void onLinearRegression();
};

#endif // MAINWINDOW_H
