#ifndef VECTOR_H
#define VECTOR_H

class Vector
{
public:
	Vector(long size);
	Vector(Vector *other);
	void reassign(Vector *newVector);
	void normalize();
	double norm2();
	void assignValuesTo(double newValue);
	double dotProduct(Vector *other);
	~Vector();

	long size;
	double *values;
};

#endif // VECTOR_H
