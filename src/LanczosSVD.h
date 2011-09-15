#ifndef LANCZOSSVD_H
#define LANCZOSSVD_H

#include "Vector.h"
#include "Matrix.h"

class LanczosSVD
{
public:
	LanczosSVD(Matrix *inputMatrix, int desiredRank);
	void solve();

private:
	void lanczosStep(int i);
	void init();

	Vector *currentVector;
	Vector *previousVector;
	Matrix *basisMatrix;
	double beta;
	double alpha;
	Matrix *basis;
	Matrix *triDiagonalMatrix;
	long desiredRank;
	Matrix *inputMatrix;
};

#endif // LANCZOSSVD_H
