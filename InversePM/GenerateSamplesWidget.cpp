#include "GenerateSamplesWidget.h"

GenerateSamplesWidget::GenerateSamplesWidget(QWidget* parent) : QDialog((QWidget*)parent) {
	// set up the UI
	ui.setupUi(this);

	ui.lineEditNumSamples->setText("2000");
	ui.radioButtonNormalizeData->setChecked(false);

	connect(ui.okButton, SIGNAL(clicked()), this, SLOT(accept()));
	connect(ui.cancelButton, SIGNAL(clicked()), this, SLOT(reject()));
}

/*void GenerateSamplesWidget::onOK() {
	max_round = ui.lineEditMaxRound->text().toInt();
	max_step = ui.lineEditMaxStep->text().toInt();

	this->accept();
}*/
