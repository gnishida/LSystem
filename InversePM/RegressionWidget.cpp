#include "RegressionWidget.h"

RegressionWidget::RegressionWidget(QWidget* parent) : QDialog((QWidget*)parent) {
	// set up the UI
	ui.setupUi(this);

	ui.lineEditStatsGridSize->setText("5");
	ui.checkBoxSaveImages->setChecked(false);

	connect(ui.okButton, SIGNAL(clicked()), this, SLOT(accept()));
	connect(ui.cancelButton, SIGNAL(clicked()), this, SLOT(reject()));
}
