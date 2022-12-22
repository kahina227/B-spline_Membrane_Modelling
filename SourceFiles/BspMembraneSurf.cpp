#include "BspMembraneSurf.h"

// Default constructor, with default values of the data members
BspMembraneSurf::BspMembraneSurf() {
	bspOrder = 4;
	membraneLength = 20;
	numSegments = 7;
	xCoord = 3;
	yCoord = 3;
	density = 0.01;
	damping = 0.02;
	tension = 3625000.0;
};

// Alternate Constructor with parameters, and default values of the data members
BspMembraneSurf::BspMembraneSurf(const int order, const double length, const int numSegts, const double x, const double y,
	const double adensity, const double adamping, const double atension) : bspOrder{order}, membraneLength{length}, 
	numSegments{numSegts}, xCoord{x}, yCoord{y}, density{adensity}, damping{adamping}, tension{atension} {
	bspOrder = 4;
	membraneLength = 20;
	numSegments = 7;
	xCoord = 3;
	yCoord = 3;
	density = 0.01;
	damping = 0.02;
	tension = 3625000.0;
}

// ACCESSORS METHODS
// Set the B-spline order
void BspMembraneSurf::setBspOrder(const int& order) {
	bspOrder = order;
}

// Set the length of the membrane
void BspMembraneSurf::setMembraneLength(const double& length) {
	membraneLength = length;
}

// Set the number of segments
void BspMembraneSurf::setNumSegments(const int& numSegts) {
	numSegments = numSegts;
}

// Set the initial control point coordinates in u and v directions
void BspMembraneSurf::setXYcoordinate(const double& x, const double& y) {
	xCoord = x;
	yCoord = y;
}

// Set the density
void BspMembraneSurf::setDensity(const double& adensity) {
	density = adensity;
}

// Set the damping
void BspMembraneSurf::setDamping(const double& adamping) {
	damping = adamping;
}

// Set the tension
void BspMembraneSurf::setTension(const double& atension) {
	tension = atension;
}


// COMPUTE THE B-SPLINE KNOT SET, ALSO KNOWN AS THE KNOT VECTOR
/*
* Compute the knot vectors in u and v directions of the membrane
* @result return the knot vector in u and the knot vector in v
* @see report for detail explanation
*/
std::vector<double> BspMembraneSurf::ComputeKnotVector(int bspOrder, double membraneLength, int numSegments) {

	double knotSpan = membraneLength / (double)numSegments;

	// Compute number of knots from user input: number of segements and B-spline order
	int numKnots = numSegments - 1 + 2 * bspOrder;

	// Define and initialise the knots vector
	std::vector<double> knotsVect(numKnots);

	// Calculate the values of the knots vector
	int count = 1;
	for (int i = 0; i <= numKnots - 1; ++i) {
		if (i <= (bspOrder - 1)) {
			knotsVect[i] = 0.0;
		}
		if (i >= bspOrder && i < numKnots - bspOrder) {
			knotsVect[i] = knotSpan * count;
			++count;
		}
		if (i >= numKnots - bspOrder) {
			knotsVect[i] = membraneLength;
		}
	}
	return knotsVect;
}


// SET THE INITIAL CONIDTION OF THE SURFACE OF THE MEMBRANE
/*
* Set the initial surface condition - user input of the plucked point coordinates
* @see report for detail explanation
*/
MMatrix<Point1D> BspMembraneSurf::setICSurfaceMatrix(double membraneLength, int bspOrder, int numSegments, double xCoord, double yCoord) {

	int numControlPts = (numSegments - 1 + 2 * bspOrder) - bspOrder;
	std::vector<double> kts = BspMembraneSurf::ComputeKnotVector(bspOrder, membraneLength, numSegments);

	// order in u and order in v directions are equals both to B-spline order bspOrder
	// since the membrane is square, bspOrder used in both u and v
	// create and initialise the elements of the control points matrix to zero
	MMatrix<Point1D> icControlPtsMatrix(numControlPts, numControlPts);
	for (int i = 0; i < numControlPts; ++i) {
		for (int j = 0; j < numControlPts; ++j) {
			icControlPtsMatrix[i][j] = 0;
		}	
	}

	// compute the fisrt possible control point that can be plucked
	int interiorControlPts = numControlPts - 2;
	double controlPtsSpan = membraneLength / interiorControlPts;
	double firstPluckedPt = 1.0; 
	double calculPosition = controlPtsSpan; 
	int firstPoint = 1;
	while (calculPosition <= xCoord && calculPosition <= yCoord) {
		calculPosition = controlPtsSpan * (double)firstPoint;
		if (calculPosition >= xCoord && calculPosition >= yCoord) {
			firstPluckedPt = firstPoint; 
		}
		firstPoint++;
	}

	// Set z coordinate of the plucked control point and is fixed to the quarter of the surface length
	//std::ofstream ofs_ICmatrix("ResultFiles/icMatrix.dat");
	double zCoord = membraneLength / 4.0;
	if (xCoord - firstPluckedPt >= 2.0) {
		for (int i = firstPluckedPt; i <= firstPluckedPt + 2.0; ++i) {
			for (int j = firstPluckedPt; j <= firstPluckedPt + 2.0; ++j) {
				icControlPtsMatrix[i][j] = zCoord;
			}
		}
	}
	if (xCoord - firstPluckedPt < 2.0) {
		for (int i = firstPluckedPt; i <= firstPluckedPt + 2.0; ++i) {
			for (int j = firstPluckedPt; j <= firstPluckedPt + 2.0; ++j) {
				icControlPtsMatrix[i][j] = zCoord;
			}
		}
	}

	// Create the initial surface in 3D
	std::ofstream ofs_ICsurface("ResultFiles/icSurf_d0.dat");
	PBspSurf3D icSurf_d0 = PBspSurf3D(FBspSurf(icControlPtsMatrix, kts, kts, bspOrder, bspOrder, numControlPts, numControlPts));
	ofs_ICsurface << icSurf_d0;

	return icControlPtsMatrix;
}



// COMPUTE THE MASS, DAMPING AND FORCE MATRICES
/**
* Compute mass matrix
*/
MMatrix<double> BspMembraneSurf::computeMassMatrix(int bspOrder, int numSegments, double membraneLength, double density) {

	int numKnots = numSegments - 1 + 2 * bspOrder;
	std::vector<double> kts = BspMembraneSurf::ComputeKnotVector(bspOrder, membraneLength, numSegments);

	BspCurvBasisFuncSet b(kts, bspOrder, numKnots);
	MMatrix<double> a1_mat = b.CreateMatrixMinimisation(0, 0.0, membraneLength);
	MMatrix<double> JTJ = Math::kronecker(a1_mat, a1_mat);
	MMatrix<double> massMatrix = Math::mmult(density, JTJ);
	return massMatrix;
}

/**
* Compute damping matrix
*/
MMatrix<double> BspMembraneSurf::computeDampingMatrix(int bspOrder, int numSegments, double membraneLength, double damping) {
	
	int numKnots = numSegments - 1 + 2 * bspOrder;
	std::vector<double> kts = BspMembraneSurf::ComputeKnotVector(bspOrder, membraneLength, numSegments);

	BspCurvBasisFuncSet b(kts, bspOrder, numKnots);
	MMatrix<double> a1_mat = b.CreateMatrixMinimisation(0, 0.0, membraneLength);
	MMatrix<double> JTJ = Math::kronecker(a1_mat, a1_mat);
	MMatrix<double> dampingMatrix = Math::mmult(damping, JTJ);
	return dampingMatrix;
}

/**
* Compute force matrix
*/
MMatrix<double> BspMembraneSurf::computeForceMatrix(int bspOrder, int numSegments, double membraneLength, double tension) {
	
	int numKnots = numSegments - 1 + 2 * bspOrder;
	std::vector<double> kts = BspMembraneSurf::ComputeKnotVector(bspOrder, membraneLength, numSegments);

	BspCurvBasisFuncSet b(kts, bspOrder, numKnots);
	MMatrix<double> a1_mat = b.CreateMatrixMinimisation(0, 0.0, membraneLength);
	MMatrix<double> JTJ = Math::kronecker(a1_mat, a1_mat);

	// Force matrix
	MMatrix<double> a2_mat = b.CreateMatrixMinimisation(1, 0.0, membraneLength);
	MMatrix<double> JuTJu = Math::kronecker(a2_mat, a1_mat);
	MMatrix<double> JvTJv = Math::kronecker(a1_mat, a2_mat);
	MMatrix<double> forceMatrix = Math::mmult(tension, Math::add(JuTJu, JvTJv));
	return forceMatrix;
}


// COMPUTE THE BOUNDARY CONDITION
/**
* Compute the reduced matrix system, which represent the boundary condition
* @see report for detail explanation
*/
MMatrix<double> BspMembraneSurf::setBoundaryConditionMatrix(int bspOrder, int numSegments, double membraneLength) {
	
	// There are 4 constraints, one for each control points governing the edges of the membrane
	int sideDim = (numSegments - 1 + 2 * bspOrder) - bspOrder; // number of control points
	int redSideDim = sideDim - 2;
	int dim = membraneLength;
	int numConstraints = 4 * sideDim - 4;

	/**
	* Implement the boundary conditions: create a constraint matrix that has 1 in the governing control points 
	* in the perimeter of square, which has zero deflection.
	*/
	// RHS zero vector
	Vector<double> bcRhs(numConstraints, 0.0);
	
	// Construct full constraint matrix
	MMatrix<double> bcMatrix(numConstraints, sideDim * sideDim + 1);

	// Arrange 1s in the appropriate places
	Vector<double> zero(sideDim * sideDim, 0.0);
	Vector<double> temp(sideDim * sideDim);

	// top side...
	for (int i = 0; i < sideDim; i++) {
		temp = zero;
		temp[i] = 1.0;
		bcMatrix.InsertRow(temp, i);
	}
	// bottom side
	for (int i = 0; i < sideDim; i++) {
		temp = zero;
		temp[i + (sideDim - 1) * (sideDim)] = 1.0;
		bcMatrix.InsertRow(temp, i + (sideDim + 2 * redSideDim));
	}
	// left and right sides
	for (int i = 1, k = sideDim; i < sideDim - 1; i++) {
		// left side
		temp = zero;
		temp[i * sideDim] = 1.0;
		bcMatrix.InsertRow(temp, k);
		k++;

		// right side
		temp = zero;
		temp[i * sideDim + (sideDim - 1)] = 1.0;
		bcMatrix.InsertRow(temp, k);
		k++;
	}

	/**
	* Insert the rhs vector into the constraint matrix, which completes
	* the formation of the constraint system matrix bcMatrix
	*/
	bcMatrix.InsertCol(bcRhs, dim);

	return bcMatrix;
}

// GET APPROPRIATE MASS, DAMPING AND FORCE MATRICES
// @see report for details
// Compute M
MMatrix<double> BspMembraneSurf::computeMatrix_M(int bspOrder, int numSegments, double membraneLength, double density) {
	MMatrix<double> massMatrix = BspMembraneSurf::computeMassMatrix(bspOrder, numSegments, membraneLength, density);
	MMatrix<double> bcMatrix = BspMembraneSurf::setBoundaryConditionMatrix(bspOrder, numSegments, membraneLength);
	Vector<int> MassElim = Math::GetEliminateIndices(massMatrix, bcMatrix);
	MMatrix<double> M = Math::EliminateMVariables(massMatrix, bcMatrix, MassElim);
	return M;
}

// Compute D
MMatrix<double> BspMembraneSurf::computeMatrix_D(int bspOrder, int numSegments, double membraneLength, double damping) {
	MMatrix<double> dampingMatrix = BspMembraneSurf::computeDampingMatrix(bspOrder, numSegments, membraneLength, damping);
	MMatrix<double> bcMatrix = BspMembraneSurf::setBoundaryConditionMatrix(bspOrder, numSegments, membraneLength);
	Vector<int> DampingElim = Math::GetEliminateIndices(dampingMatrix, bcMatrix);
	MMatrix<double> D = Math::EliminateMVariables(dampingMatrix, bcMatrix, DampingElim);
	return D;
}

// Compute F
MMatrix<double> BspMembraneSurf::computeMatrix_F(int bspOrder, int numSegments, double membraneLength, double tension) {
	MMatrix<double> forceMatrix = BspMembraneSurf::computeForceMatrix(bspOrder, numSegments, membraneLength, tension);
	MMatrix<double> bcMatrix = BspMembraneSurf::setBoundaryConditionMatrix(bspOrder, numSegments, membraneLength);
	Vector<int> ForceElim = Math::GetEliminateIndices(forceMatrix, bcMatrix);
	MMatrix<double> F = Math::EliminateMVariables(forceMatrix, bcMatrix, ForceElim);
	return F;
}

// COMPUTE THE MATRICES A, B AND C
// @see report for details
// Compute Matrix C
MMatrix<double> BspMembraneSurf::computeMatrix_C(int bspOrder, int numSegments, double membraneLength, double density, double damping) {
	MMatrix<double> M = BspMembraneSurf::computeMatrix_M(bspOrder, numSegments, membraneLength, density);
	MMatrix<double> D = BspMembraneSurf::computeMatrix_D(bspOrder, numSegments, membraneLength, damping);
	MMatrix<double> C = Math::Inverse(Math::add(Math::mmult(2.0, M), Math::mmult(dt, D)));
	return C;
}

// Compute Matrix A
MMatrix<double> BspMembraneSurf::computeMatrix_A(int bspOrder, int numSegments, double membraneLength, double density, double tension, double damping) {
	MMatrix<double> C = BspMembraneSurf::computeMatrix_C(bspOrder, numSegments, membraneLength, density, damping);
	MMatrix<double> M = BspMembraneSurf::computeMatrix_M(bspOrder, numSegments, membraneLength, density);
	MMatrix<double> F = BspMembraneSurf::computeMatrix_F(bspOrder, numSegments, membraneLength, tension);
	MMatrix<double> A = Math::mult2(C, (Math::subtract(Math::mmult(4.0, M), Math::mmult(2.0 * dt * dt, F))));
	return A;
}

// Compute Matrix B
MMatrix<double> BspMembraneSurf::computeMatrix_B(int bspOrder, int numSegments, double membraneLength, double damping, double density) {
	MMatrix<double> C = BspMembraneSurf::computeMatrix_C(bspOrder, numSegments, membraneLength, density, damping);
	MMatrix<double> D = BspMembraneSurf::computeMatrix_D(bspOrder, numSegments, membraneLength, damping);
	MMatrix<double> M = BspMembraneSurf::computeMatrix_M(bspOrder, numSegments, membraneLength, density);
	MMatrix<double> B = Math::mult2(C, (Math::subtract(Math::mmult(dt, D), Math::mmult(2.0, M))));
	return B;
}

// COMPUTE d0, (I-B) and d1
/**
* Compute d0
*/

Vector<Point1D> BspMembraneSurf::compute_d0(double membraneLength, int bspOrder, int numSegments, double density, double damping, 
																					double tension,	double xCoord, double yCoord) {
	
	int sideDim = (numSegments - 1 + 2 * bspOrder) - bspOrder; 
	int redSideDim = sideDim - 2;
	int dim = membraneLength;
	
	// Get values of C, A and B matrices
	MMatrix<double> C = BspMembraneSurf::computeMatrix_C(bspOrder, numSegments, membraneLength, density, damping);
	MMatrix<double> A = BspMembraneSurf::computeMatrix_A(bspOrder, numSegments, membraneLength, density, tension, damping);
	MMatrix<double> B = BspMembraneSurf::computeMatrix_B(bspOrder, numSegments, membraneLength, damping, density);

	MMatrix<Point1D> icSurfMatrix = BspMembraneSurf::setICSurfaceMatrix(membraneLength, bspOrder, numSegments, xCoord, yCoord);

	// Assume the plucked control points of the surface are in the matrix d - the control points matrix
	MMatrix<Point1D> d(sideDim, sideDim);
	
	// The initial condition, initial plucked control points matrix
	for (int i = 0; i < sideDim; ++i) {
		for (int j = 0; j < sideDim; ++j) {
			d[i][j] = icSurfMatrix[i][j];
		}
	}

	/**
	* Store the control points matrix d as a vector d0
	* and converting the initial condition matrix to a vector
	*/
	Vector<Point1D> d0 = Math::CreateKroneckerVector(d);

	return d0;
}

// Compute (I-B)
MMatrix<Point1D> BspMembraneSurf::computeMatrix_IsubB(double membraneLength, int bspOrder, int numSegments, double damping, double density) {
	
	int sideDim = (numSegments - 1 + 2 * bspOrder) - bspOrder; 
	int redSideDim = sideDim - 2;

	// Get values of B
	MMatrix<double> B = BspMembraneSurf::computeMatrix_B(bspOrder, numSegments, membraneLength, damping, density);
	int redDim = (sideDim - 2) * (sideDim - 2);

	MMatrix<double> I(redDim, redDim, 0.0);

	// The diagonal elements of the matrix equal to 1 and the upper and lower are set to 0
	for (int i = 0; i < redDim; i++)
		I[i][i] = 1.0;

	MMatrix<double> IsubB = Math::subtract(I, B); // I - B
	return IsubB;
}

// Compute d0_1
Vector<Point1D> BspMembraneSurf::compute_d0_1(double membraneLength, int bspOrder, int numSegments, double density, double damping, double tension,
																										double xCoord, double yCoord) {
	
	int sideDim = (numSegments - 1 + 2 * bspOrder) - bspOrder;
	int redSideDim = sideDim - 2;
	int redDim = (sideDim - 2) * (sideDim - 2);

	Vector<Point1D> d0 = BspMembraneSurf::compute_d0(membraneLength, bspOrder, numSegments, density, damping, tension, xCoord, yCoord);
	// Create rhs vector
	Vector<Point1D> d0_1(redDim, 0.0);

	// Write d0 into this rhs vector d0_1
	for (int i = 0; i < sideDim - 2; i++) {
		for (int j = 0; j < sideDim - 2; j++) {
			d0_1[i * (sideDim - 2) + j] = d0[(i + 1) * sideDim + (j + 1)];
		}
	}
	return d0_1;
}

/**
* Compute d1 from d0 (control points from initial plucked condition)
* Obtained from equation (I-B)*d1 = A*d0
* d1 = (I-B)^{-1} * Ad0
*/
Vector<Point1D> BspMembraneSurf::compute_d1(double membraneLength, int bspOrder, int numSegments, double density, double damping, 
																				double tension, double xCoord, double yCoord) {
	
	int sideDim = (numSegments - 1 + 2 * bspOrder) - bspOrder;
	
	//MMatrix<double> A = computeMatrixA(bspOrder, numSegments, membraneLength, density, tension, damping);
	MMatrix<Point1D> A = BspMembraneSurf::computeMatrix_A(bspOrder, numSegments, membraneLength, density, tension, damping);
	Vector<Point1D> d0_1 = BspMembraneSurf::compute_d0_1(membraneLength, bspOrder, numSegments, density, damping, tension, xCoord, yCoord);
	
	// multiply d0_1 by A
	Vector<Point1D> Ad = Math::mult3(A, d0_1);

	MMatrix<Point1D> IsubB = BspMembraneSurf::computeMatrix_IsubB(membraneLength, bspOrder, numSegments, damping, density);

	// calculate d1
	// Solve (I - B)^{-1} * A d0_1
	Vector<Point1D> d_1 = Math::gauss(IsubB, Ad, Ad.GetNum());

	// re-create the full sized d1
	Vector<Point1D> d1(sideDim * sideDim, 0.0);
	
	// re-create the full sized d1
	for (int i = 1; i < sideDim - 1; i++) {
		for (int j = 1; j < sideDim - 1; j++) {
			d1[i * sideDim + j] = d_1[(i - 1) * (sideDim - 2) + (j - 1)];
		}
	}

	// convert d1 to a Matrix<Point1D> object
	Math MatFromVect;
	MMatrix<Point1D> d1Matrix = MatFromVect.MatrixFromVector(d1, sideDim);

	// create a functional b-spline surface with d1
	std::vector<double> kts = BspMembraneSurf::ComputeKnotVector(bspOrder, membraneLength, numSegments);
	int numControlPts = (numSegments - 1 + 2 * bspOrder) - bspOrder;
	FBspSurf d1Surf(d1Matrix, kts, kts, bspOrder, bspOrder, numControlPts, numControlPts);
	PBspSurf3D surf1(d1Surf);

	// Print the d1 surface to .dat file
	std::ofstream ofs_surf_1("ResultFiles/surf_d1.dat");
	ofs_surf_1 << surf1;

	return d1;
}

// COMPUTE d2 FROM d1 AND d0 AND ALL OTHER DERIVATIVES OF THE FDM
/**
* Compute d2 from d1 and d0
* d2 = A*d1 + B*d0, where:
* C = (2*M + dt*D)^{-1}, A = C*(4*M - 2*(dt^2)*F), and B = C*(dt*D - 2*M)
*/
void BspMembraneSurf::computeSurfacesFDM(double membraneLength, int bspOrder, int numSegments, double density, double damping, 
																				double tension,	double xCoord, double yCoord) {
	int sideDim = (numSegments - 1 + 2 * bspOrder) - bspOrder;
	int redSideDim = sideDim - 2;
	int redDim = (sideDim - 2) * (sideDim - 2);

	Vector<Point1D> d1_1(redDim);
	
	// Get the knot vector
	std::vector<double> kts = BspMembraneSurf::ComputeKnotVector(bspOrder, membraneLength, numSegments);
	
	// Get d0
	Vector<Point1D> d0 = BspMembraneSurf::compute_d0(membraneLength, bspOrder, numSegments, density, damping, tension, xCoord, yCoord);
	
	// Get d0_1
	Vector<Point1D> d0_1 = BspMembraneSurf::compute_d0_1(membraneLength, bspOrder, numSegments, density, damping, tension, xCoord, yCoord);
	
	// Get d1
	Vector<Point1D> d1 = BspMembraneSurf::compute_d1(membraneLength, bspOrder, numSegments, density, damping, tension, xCoord, yCoord);
	
	// Get A and B
	MMatrix<double> A = BspMembraneSurf::computeMatrix_A(bspOrder, numSegments, membraneLength, density, tension, damping);
	MMatrix<double> B = BspMembraneSurf::computeMatrix_B(bspOrder, numSegments, membraneLength, damping, density);

	// Set the Timer
	clock_t time_exec;
	time_exec = clock();
	std::cerr << time_exec;
	Vector<Point1D> d2(sideDim * sideDim, 0.0);
	float t = 0.0;
	do {
		for (int i = 0; i < sideDim - 2; i++)
			for (int j = 0; j < sideDim - 2; j++)
				d0_1[i * (sideDim - 2) + j] = d0[(i + 1) * sideDim + (j + 1)];

		for (int i = 0; i < sideDim - 2; i++)
			for (int j = 0; j < sideDim - 2; j++)
				d1_1[i * (sideDim - 2) + j] = d1[(i + 1) * sideDim + (j + 1)];

		// Compute the new control points for this timeIndex
		Vector<Point1D> d_2 = Math::vadd(Math::mult3(A, d1_1), Math::mult3(B, d0_1));

		// Re-create full sized d_2

		for (int i = 1; i < sideDim - 1; i++)
			for (int j = 1; j < sideDim - 1; j++)
				d2[i * sideDim + j] = d_2[(i - 1) * (sideDim - 2) + (j - 1)];

		// Rack the control points vectors along one, ready for next iteration
		d0 = d1;
		d1 = d2;
		clock_t time = clock() - time_exec;
		t = (float)(time / CLOCKS_PER_SEC);
		std::cerr << t;

		// Re-create the matrix Md from the vector d_{i}
		Math MatFromVect;
		MMatrix<Point1D> Md = MatFromVect.MatrixFromVector(d2, sideDim);

		// Create a functional b-spline with d_{i}
		FBspSurf dSurf(Md, kts, kts, bspOrder, bspOrder, sideDim, sideDim);

		// Create a parametric b-spline with d_{i}
		PBspSurf3D surf(dSurf);

		// Print the surfaces to .dat files
		std::ofstream ofs_surfs;
		ofs_surfs.open("ResultFiles/surf_d" + std::to_string(t) + ".dat");
		ofs_surfs << surf;
		ofs_surfs << std::endl;
		ofs_surfs << std::endl;
	} while (t < 5.0);
}
