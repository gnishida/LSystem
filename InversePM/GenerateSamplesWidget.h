#pragma once

#include <QWidget>
#include "ui_GenerateSamplesWidget.h"

class GenerateSamplesWidget : public QDialog {
Q_OBJECT

public:
	Ui::GenerateSamplesWidget ui;

public:
	GenerateSamplesWidget(QWidget* parent);

public slots:
	//void onOK();
};

