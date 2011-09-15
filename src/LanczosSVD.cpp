#include "LanczosSVD.h"
#include <cmath>

LanczosSVD::LanczosSVD(Matrix *inputMatrix, int desiredRank)
{
	this->desiredRank = desiredRank;
	this->inputMatrix = inputMatrix;
}

void LanczosSVD::solve()
{
	this->init();

	for (int i = 0 ; i < desiredRank ; i++)
	{
		this->lanczosStep(i);
	}
	// find eigen values and vectors of triDiagonalMatrix
}

void LanczosSVD::lanczosStep(int i)
{
	Vector nextVector = inputMatrix->timesSquared(currentVector);
	// calc scale factor of next vector
	// scale next vector
	// nextVector = nextVector - beta*PrevVector
	double alpha = currentVector->dotProduct(&nextVector);
	// nextVector = nextVector - alpha*betacurrentVector
	// orthoganalize Against All But Last ( nextVector and basisMatrix )
	beta = nextVector.norm2();
	// scale nextVector to 1/beta
	basisMatrix->assignRow(i, &nextVector);
	previousVector->reassign(currentVector);
	currentVector->reassign(&nextVector);

	triDiagonalMatrix->values[i-1][i-1] = alpha;
	if (i < desiredRank)
	{
		triDiagonalMatrix->values[i-1][i] = beta;
		triDiagonalMatrix->values[i][i-1] = beta;
	}
}

void LanczosSVD::init()
{
	this->currentVector = new Vector(inputMatrix->cols);
	this->currentVector->assignValuesTo(1.0 / sqrt(inputMatrix->cols));
}
