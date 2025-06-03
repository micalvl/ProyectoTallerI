/**
 *  @file   Matrix.cpp
 *  @brief  Matrix's methods
 *  @author [Original Author]
 *  @transcribed by Miguel Calvo León
 *  @date   2025-03-27
 ***********************************************/

#include "../include/Matrix.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
using namespace std;


Matrix::Matrix(int fil, int col) : fil(fil), col(col)
{
    initMatrix();
}

Matrix::Matrix(int fil, int col, double v[], int n): fil(fil), col(col)
{
    initMatrix();

    int k = 0;

    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++){
            if (k < n)
                matrix[i][j] = v[k++];
            else
                matrix[i][j] = 0;
        }
}

Matrix::Matrix() : fil(0), col(0), matrix(nullptr) {}

Matrix::Matrix(const Matrix& m): fil(m.fil), col(m.col)
{
    initMatrix();
    for (int i = 0; i < fil; ++i)
        for (int j = 0; j < col; ++j)
            matrix[i][j] = m.matrix[i][j];
}


void Matrix::initMatrix()
{
    matrix = new double*[fil];
    for (int i = 0; i < fil; i++)
        matrix[i] = new double[col];

    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++)
            matrix[i][j] = 0.0;
}

Matrix::~Matrix()
{
    for (int i = 0; i < fil; i++)
        delete[] matrix[i];

    delete[] matrix;
}

Matrix& Matrix::operator=(const Matrix& matrix2)
{
    if (this != &matrix2) {
        if (fil != matrix2.fil || col != matrix2.col) {
            for (int i = 0; i < fil; ++i) delete[] matrix[i];
            delete[] matrix;
            fil = matrix2.fil;
            col = matrix2.col;

            initMatrix();
        }

        for (int i = 0; i < fil; ++i)
            for (int j = 0; j < col; ++j)
                matrix[i][j] = matrix2.matrix[i][j];
    }
    return *this;
}

Matrix Matrix::operator+(const Matrix& matrix2) const
{
    Matrix result(fil, col);

    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++)
            result.matrix[i][j] = matrix[i][j] + matrix2.matrix[i][j];

    return result;
}

Matrix Matrix::operator-(const Matrix& matrix2) const
{
    Matrix result(fil, col);

    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++)
            result.matrix[i][j] = matrix[i][j] - matrix2.matrix[i][j];

    return result;
}

Matrix Matrix::operator*(const Matrix& matrix2) const
{
    Matrix result(fil, matrix2.col);

    for (int i = 0; i < this->fil ; i++){
        for (int j = 0; j < matrix2.col; j++){
            result.matrix[i][j] = 0;
            for (int k = 0; k < this->col; k++){
                result.matrix[i][j] = result.matrix[i][j] + this->matrix[i][k] * matrix2.matrix[k][j];
            }
        }
    }

    return result;
}


double& Matrix::operator()(const int i, const int j) const
{
    // I created that to see where it be my error
    if (i < 1 || i > fil || j < 1 || j > col) {
        std::cerr
                << "ERROR: Matrix index out of range ("
                << i << "," << j << ") for Matrix of size "
                << fil << "x" << col << "\n";
        std::abort();
    }
    return matrix[i-1][j-1];
}

void Matrix::print()
{
    for (int i = 0; i < fil; i++){
        for (int j = 0; j < col; j++){
            std::cout << std::fixed << std::setprecision(14) << matrix[i][j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

double Matrix::norm() const {
    double suma = 0.0;
    for (int i = 0; i < fil; ++i) {
        for (int j = 0; j < col; ++j) {
            suma += matrix[i][j] * matrix[i][j];
        }
    }
    return sqrt(suma);
}

void Matrix::printMatrixValues(float** arr, int n, int m){
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            cout<<arr[i][j]<<"\t";
        }
        cout<<endl;
    }
    return;
}

void Matrix::printInverseMatrix(float** arr, int n, int m){
    for (int i = 0; i < n; i++) {
        for (int j = n; j < m; j++) {
            printf("%.3f\t", arr[i][j]);
        }
        cout<<endl;
    }
    return;
}

void Matrix::findInvMatGaussJordan(float** mat, int order){
    float temp;
    printf("The inverse of matrix : A = \n");
    printMatrixValues(mat, order, order);
    for (int i = 0; i < order; i++) {
        for (int j = 0; j < 2 * order; j++) {
            if (j == (i + order))
                mat[i][j] = 1;
        }
    }
    for (int i = order - 1; i > 0; i--) {
        if (mat[i - 1][0] < mat[i][0]) {
            float* temp = mat[i];
            mat[i] = mat[i - 1];
            mat[i - 1] = temp;
        }
    }
    for (int i = 0; i < order; i++) {
        for (int j = 0; j < order; j++) {
            if (j != i) {
                temp = mat[j][i] / mat[i][i];
                for (int k = 0; k < 2 * order; k++) {
                    mat[j][k] -= mat[i][k] * temp;
                }
            }
        }
    }
    for (int i = 0; i < order; i++) {
        temp = mat[i][i];
        for (int j = 0; j < 2 * order; j++) {
            mat[i][j] = mat[i][j] / temp;
        }
    }
    cout<<"A' =\n";
    printInverseMatrix(mat, order, 2 * order);
    return;
}

Matrix Matrix::inverse(){
    int n = fil;
    double** A = new double*[n];
    for (int i = 0; i < n; ++i) {
        A[i] = new double[2*n];
        for (int j = 0; j < n; ++j) {
            A[i][j] = (*this)(i+1, j+1);
        }
        for (int j = 0; j < n; ++j) {
            A[i][n+j] = (i == j ? 1.0 : 0.0);
        }
    }

    for (int p = 0; p < n; ++p) {
        int maxr = p;
        double maxv = fabs(A[p][p]);
        for (int r = p+1; r < n; ++r) {
            double v = fabs(A[r][p]);
            if (v > maxv) {
                maxv = v;
                maxr = r;
            }
        }

        if (maxr != p) {
            std::swap(A[p], A[maxr]);
        }

        double diag = A[p][p];
        for (int c = 0; c < 2*n; ++c)
            A[p][c] /= diag;

        for (int r = 0; r < n; ++r) {
            if (r == p) continue;
            double factor = A[r][p];
            for (int c = 0; c < 2*n; ++c)
                A[r][c] -= factor * A[p][c];
        }
    }

    Matrix inv(n, n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            inv(i+1, j+1) = A[i][n + j];
        }
    }

    for (int i = 0; i < n; ++i)
        delete[] A[i];
    delete[] A;

    return inv;
}

int Matrix::getFilas() const { return fil; }
int Matrix::getColumnas() const { return col; }


Matrix Matrix::transpose() const {
    Matrix result(col, fil);
    for (int i = 1; i <= fil; ++i)
        for (int j = 1; j <= col; ++j)
            result(j, i) = (*this)(i, j);
    return result;
}

Matrix Matrix::identity(int size) {
    Matrix I(size, size);
    for (int i = 1; i <= size; ++i)
        I(i, i) = 1.0;
    return I;
}

Matrix Matrix::opsc(double scalar) const {
    Matrix result(fil, col);
    for (int i = 0; i < fil; ++i) {
        for (int j = 0; j < col; ++j) {
            result.matrix[i][j] = this->matrix[i][j] * scalar;
        }
    }
    return result;
}

Matrix Matrix::divsc(double scalar) const {
    Matrix result(fil, col);
    for (int i = 0; i < fil; ++i) {
        for (int j = 0; j < col; ++j) {
            result.matrix[i][j] = this->matrix[i][j] / scalar;
        }
    }
    return result;
}

Matrix Matrix::zeros(int rows, int cols) {
    return Matrix(rows, cols);
}

Matrix Matrix::cross(const Matrix& a, const Matrix& b) {
    Matrix c(3,1);
    c(1,1) = a(2,1)*b(3,1) - a(3,1)*b(2,1);
    c(2,1) = a(3,1)*b(1,1) - a(1,1)*b(3,1);
    c(3,1) = a(1,1)*b(2,1) - a(2,1)*b(1,1);
    return c;
}


double Matrix::dot(const Matrix& a, const Matrix& b) {
    return a(1,1)*b(1,1)
           + a(2,1)*b(2,1)
           + a(3,1)*b(3,1);
}


/*
int main(){
    int order = 3;
    float** mat = new float*[20];
    for (int i = 0; i < 20; i++)
        mat[i] = new float[20];
    mat[0][0] = 6; mat[0][1] = 9; mat[0][2] = 5;
    mat[1][0] = 8; mat[1][1] = 3; mat[1][2] = 2;
    mat[2][0] = 1; mat[2][1] = 4; mat[2][2] = 7;
    findInvMatGaussJordan(mat, order);
    return 0;
}
*/

