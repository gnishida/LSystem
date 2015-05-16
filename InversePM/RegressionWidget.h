#pragma once

#include <QWidget>
#include "ui_RegressionWidget.h"

class RegressionWidget : public QDialog {
Q_OBJECT

public:
	Ui::RegressionWidget ui;

public:
	RegressionWidget(QWidget* parent);

public slots:
	//void onOK();
};

