#ifndef MATRIX_H
#define MATRIX_H

#include "Vector.h"

class Matrix
{
public:
	Matrix(long rows, long cols);
	~Matrix();

	void assignRow(long row, Vector *vector);
	Matrix times(Matrix *matrix);
	Vector times(Vector *vector);
	Vector timesSquared(Vector *vector);

	double getValue(long row, long col);

	long cols;
	long rows;
	double **values;
};

#endif // MATRIX_H
