#include "GenerateSamplesWidget.h"

GenerateSamplesWidget::GenerateSamplesWidget(QWidget* parent) : QDialog((QWidget*)parent) {
	// set up the UI
	ui.setupUi(this);

	ui.lineEditNumSamples->setText("10000");
	ui.checkBoxSaveImages->setChecked(false);

	connect(ui.okButton, SIGNAL(clicked()), this, SLOT(accept()));
	connect(ui.cancelButton, SIGNAL(clicked()), this, SLOT(reject()));
}
