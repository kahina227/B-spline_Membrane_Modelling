#ifndef BSPMEMBRANESURF
#define BSPMEMBRANESURF

#include <iostream>
#include <string>
#include <fstream>
#include <vector>

#include "PBspCurv.h"
#include "PBspSurf.h"
#include "BspSurf.h"
#include "MathFunctions.h"

class BspMembraneSurf {
	// B-spline surface parameters
	int bspOrder;
	double membraneLength;
	int numSegments;

	// Initial condition control point coordinates
	double xCoord, yCoord;

	// Physical properties
	double density, damping, tension;
	double dt = 0.000001;
public:
	// Default constructor
	BspMembraneSurf();

	// Alternate Constructor with parameters
	BspMembraneSurf(const int order, const double length, const int numSegts, const double x, const double y,
										const double adensity, const double adamping, const double atension);

	// Accessor methods
	void setBspOrder(const int& order);
	void setMembraneLength(const double& length);
	void setNumSegments(const int& numSegts);
	void setXYcoordinate(const double& x, const double& y);
	void setDensity(const double& adensity);
	void setDamping(const double& adamping);
	void setTension(const double& atension);

	/*
	* Compute the knot vectors in u and v directions of the membrane
	*/
	std::vector<double> ComputeKnotVector(int bspOrder, double membraneLength, int numSegments);

	/*
	* Set the initial surface condition - user input of the plucked point coordinates in u and v
	*/
	MMatrix<Point1D> setICSurfaceMatrix(double membraneLength, int bspOrder, int numSegments, double xCoord, double yCoord); 
	
	/*
	* Compute the mass matrix
	*/
	MMatrix<double> computeMassMatrix(int bspOrder, int numSegments, double membraneLength, double density); 

	/*
	* Compute damping matrix
	*/
	MMatrix<double> computeDampingMatrix(int bspOrder, int numSegments, double membraneLength, double damping);
	
	/*
	* Compute force matrix
	*/
	MMatrix<double> computeForceMatrix(int bspOrder, int numSegments, double membraneLength, double tension); 

	/*
	* Set boundary condition of the membrane
	*/
	MMatrix<double> setBoundaryConditionMatrix(int bspOrder, int numSegments, double membraneLength); 

	/*
	* Compute M
	*/
	MMatrix<double> computeMatrix_M(int bspOrder, int numSegments, double membraneLength, double density); 

	/*
	* Compute D
	*/
	MMatrix<double> computeMatrix_D(int bspOrder, int numSegments, double membraneLength, double damping); 

	/*
	* Compute F
	*/
	MMatrix<double> computeMatrix_F(int bspOrder, int numSegments, double membraneLength, double tension); 

	// COMPUTE THE MATRICES A, B AND C
	/*
	* Compute the matrix C
	*/
	MMatrix<double> computeMatrix_C(int bspOrder, int numSegments, double membraneLength, double density, double damping); 

	/*
	* Compute the matrix A
	*/
	MMatrix<double> computeMatrix_A(int bspOrder, int numSegments, double membraneLength, double density, double tension, double damping); 

	/*
	* Compute the matrix B
	*/
	MMatrix<double> computeMatrix_B(int bspOrder, int numSegments, double membraneLength, double damping, double density); 

	/*
	* Create the vector d0
	*/
	Vector<Point1D> compute_d0(double membraneLength, int bspOrder, int numSegments, double density, double damping, double tension,
																										double xCoord, double yCoord); 
	/*
	* Create the matrix (I-B)
	*/
	MMatrix<Point1D> computeMatrix_IsubB(double membraneLength, int bspOrder, int numSegments, double damping, double density); 

	// Compute d0_1
	/*
	* Create the vector d0_1
	*/
	Vector<Point1D> compute_d0_1(double membraneLength, int bspOrder, int numSegments, double density, double damping, double tension,
																										double xCoord, double yCoord); 

	/*
	* Create the vector d1
	*/
	Vector<Point1D> compute_d1(double membraneLength, int bspOrder, int numSegments, double density, double damping, double tension,
																									double xCoord, double yCoord);
		
	// COMPUTE THE d_{i} - THE CONTROL POINTS PARTIAL DIFFERENCTIAL EQUATION USING THE FINITE DIFFERENCE SCHEME
	/*
	* Compute d_{i} from d1 and d0
	*/
	void computeSurfacesFDM(double membraneLength, int bspOrder, int numSegments, double density, double damping, double tension,
																								double xCoord, double yCoord);
};
#endif
