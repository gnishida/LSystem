#include "MainWindow.h"
#include <QtGui/QApplication>

int main(int argc, char *argv[]) {
	if (argc < 3) {
		std::cout << "Usage: " << argv[0] << " <grid size> <grid stat size>" << std::endl;
		return -1;
	}

	QApplication a(argc, argv);
	MainWindow w(atoi(argv[1]), atoi(argv[2]));
	w.show();
	return a.exec();
}
