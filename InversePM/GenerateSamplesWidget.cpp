#include "GenerateSamplesWidget.h"

GenerateSamplesWidget::GenerateSamplesWidget(QWidget* parent) : QDialog((QWidget*)parent) {
	// set up the UI
	ui.setupUi(this);

	ui.lineEditNumGrid->setText("5");
	ui.lineEditNumStatGrid->setText("5");
	ui.lineEditNumSamples->setText("2000");
	ui.checkBoxSaveImages->setChecked(false);

	connect(ui.okButton, SIGNAL(clicked()), this, SLOT(accept()));
	connect(ui.cancelButton, SIGNAL(clicked()), this, SLOT(reject()));
}

/*void GenerateSamplesWidget::onOK() {
	max_round = ui.lineEditMaxRound->text().toInt();
	max_step = ui.lineEditMaxStep->text().toInt();

	this->accept();
}*/
