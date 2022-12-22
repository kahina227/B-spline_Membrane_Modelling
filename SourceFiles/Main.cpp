#include "BspMembraneSurf.h"
#include <iostream>
#include <string>
#include <fstream>
#include <vector>

int main() {
	// User inputs parameters: order, length and number of segments
	int order;
	std::cout << "Please input the order of the B-spline curve: ";
	std::cin >> order;
	std::cout << std::endl;

	double length;
	std::cout << "Please input the the length of the surface: ";
	std::cin >> length;
	std::cout << std::endl;

	int numSegts;
	std::cout << "Please input the number of segments: ";
	std::cin >> numSegts;
	std::cout << std::endl;

	// choosing initial control point coord
	double x, y;
	std::cout << "Please input the x and y coordinate seperated by space: ";
	std::cin >> x >> y;
	std::cout << std::endl;


	// Choosing physical properties
	double adensity, adamping, atension;
	std::cout << "Please input the density: ";
	std::cin >> adensity;
	std::cout << "Please input the damping: ";
	std::cin >> adamping;
	std::cout << "Please input the tension: ";
	std::cin >> atension;


	BspMembraneSurf* membraneSurfaces = new BspMembraneSurf(length, order, numSegts, adensity, adamping, atension, x, y);
	membraneSurfaces->computeSurfacesFDM(length, order, numSegts, adensity, adamping, atension, x, y);
	delete membraneSurfaces;
	membraneSurfaces = nullptr;

	return 0;
}
