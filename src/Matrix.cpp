#include "Matrix.h"
#include <mpi.h>

Matrix::Matrix(long rows, long cols)
{
    this->rows = rows;
    this->cols = cols;
    values = new double*[rows];
    for (int r = 0; r < rows ; r++)
    {
	values[r] = new double[cols];
    }
}

Matrix::~Matrix()
{
    for (int r = 0; r < rows ; r++)
    {
	delete[] values[r];
    }
    delete[] values;
}

void Matrix::assignRow(long row, Vector *vector)
{
    for (int i = 0; i < vector->size ; i++)
    {
	values[row][i] = vector->values[i];
    }
}

Matrix Matrix::times(Matrix *matrix)
{

    const int MASTER_RANK = 0;
    Matrix result(this->rows, matrix->cols);

    int rank, numOfProc;

    MPI_Status status;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numOfProc);

    long averageRows = this->rows / numOfProc;
    long remainedRows = this->rows / numOfProc;

    long rowOffset = 0;
    long rowsToSend = 0;
    long count = 0;

    //! Distribute matrices values to slaves
    for (int r = 0 ; r< numOfProc ; r++)
    {
	if (r != MASTER_RANK)
	    continue;

	if (r <= remainedRows)
	    rowsToSend = averageRows + 1;
	else
	    rowsToSend = averageRows;

	MPI_Send(&rowOffset, 1, MPI_LONG, r, 0, MPI_COMM_WORLD);
	MPI_Send(&rowsToSend, 1, MPI_LONG, r, 0, MPI_COMM_WORLD);

	count = rowsToSend * this->cols;

	MPI_Send(values[rowOffset], count, MPI_DOUBLE, r, 0, MPI_COMM_WORLD);

	count = this->cols * matrix->cols;
	MPI_Send(matrix,count, MPI_DOUBLE, r, 0, MPI_COMM_WORLD);

	rowOffset += rowsToSend;
    }

    //! Collect results
    for (int r = 0 ; r< numOfProc ; r++)
    {
	if (r != MASTER_RANK)
	    continue;

	MPI_Recv(&rowOffset, 1, MPI_LONG, r, 0, MPI_COMM_WORLD, &status);
	MPI_Recv(&rowsToSend, 1, MPI_LONG, r, 0, MPI_COMM_WORLD, &status);

	count = rowsToSend * matrix->cols;

	MPI_Recv(result.values[rowOffset], count, MPI_DOUBLE, r, 0, MPI_COMM_WORLD, &status);
    }


    //! Slave Nodes
    if(rank != MASTER_RANK)
    {
	MPI_Recv(&rowOffset, 1, MPI_LONG, MASTER_RANK, 0, MPI_COMM_WORLD, &status);
	MPI_Recv(&rowsToSend, 1, MPI_LONG, MASTER_RANK, 0, MPI_COMM_WORLD, &status);

	count = rowsToSend * this->cols;
	MPI_Recv(this->values, count, MPI_DOUBLE, MASTER_RANK, 0, MPI_COMM_WORLD, &status);

	count = this->cols * matrix->cols;
	MPI_Recv(matrix->values, count, MPI_DOUBLE, MASTER_RANK, 0, MPI_COMM_WORLD, &status);

	for ( long i=0 ; i < matrix->cols ; i++)
	{
	    for ( long j=0 ; j < rowsToSend ; j++)
	    {
		result.values[i][j] = 0;
		for ( long k = 0 ; k< this->cols ; k++)
		    result.values[j][i] += this->values[j][k] * matrix->values[k][i];
	    }
	}

	MPI_Send(&rowOffset, 1, MPI_LONG, MASTER_RANK, 0, MPI_COMM_WORLD);
	MPI_Send(&rowsToSend, 1, MPI_LONG, MASTER_RANK, 0, MPI_COMM_WORLD);

	count = rowsToSend * matrix->cols;
	MPI_Send(result.values, count, MPI_DOUBLE, MASTER_RANK, 0, MPI_COMM_WORLD);
    }

    return result;
}

Vector Matrix::times(Vector *vector)
{

}

Vector Matrix::timesSquared(Vector *vector)
{

}

double Matrix::getValue(long row, long col)
{
    return this->values[row][col];
}
