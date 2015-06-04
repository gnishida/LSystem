#include "ParametricLSystem.h"
#include <iostream>

int main() {
	parametriclsystem::ParametricLSystem pls;
	for (int i = 0; i < 10; ++i) {
		std::cout << pls.derive(i) << endl;
	}

	pls.drawTree();

	return 0;
}