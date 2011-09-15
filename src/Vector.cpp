#include "Vector.h"
#include <cmath>

using namespace std;

Vector::Vector(long size)
{
	this->size = size;
	this->values = new double[size];
}

Vector::~Vector()
{
	delete[] values;
}

void Vector::assignValuesTo(double newValue)
{
    for (int i=0 ; i < this->size ; i++)
	values[i] = newValue;
}

void Vector::reassign(Vector *newVector)
{
    delete[] values;
    this->size = newVector->size;
    this->values = new double[size];
    for (int i=0 ; i < this->size ; i++)
	values[i] = newVector->values[i];
}

double Vector::norm2()
{
    double total = 0.0;
    for ( int i=0; i < this->size; i++ )
    {
	total += this->values[i] * this->values[i];
    }

    return sqrt(total);
}

void Vector::normalize()
{

    // Unit the vector based on Norm-2 value
    double normValue = this->norm2();

    for ( int i = 0 ; i<this->size ; i++)
	this->values[i] /=normValue;

/*
    double maxValue = 0.0;

    double tmp;
    for ( int i=0; i < this->size; i++ )
    {
	tmp = abs(this->values[i]);
	if ( tmp > maxValue )
	    maxValue = tmp;
    }

    for ( int i=0; i < this->size; i++ )
    {
	this->values[i] = this->values[i] / maxValue;
    }
*/
}

double Vector::dotProduct(Vector *other)
{
    //! Check for same size
    if(this->size != other->size )
	return 0.0;

    double total = 0.0;
    for ( int i=0; i < this->size; i++ )
    {
	total += this->values[i] * other->values[i];
    }

    return total;
}
