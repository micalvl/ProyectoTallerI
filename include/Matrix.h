//
// Created by micalvl on 20/03/2025.
//

#ifndef PROYECTOTALLERI_MATRIX_H
#define PROYECTOTALLERI_MATRIX_H


class Matrix {
public:
    Matrix(int fil, int col);
    Matrix(int fil, int col, double v[], int n);
    Matrix(const Matrix& m);
    ~Matrix();

    Matrix& operator=(const Matrix& matrix2);
    Matrix operator+(const Matrix& matrix2);
    Matrix operator-(const Matrix& matrix2);
    Matrix operator*(const Matrix& matrix2);
    double& operator()(const int i, const int j) const;

    void print();


private:
    void initMatrix();

private:
    int fil, col, n;
    double** matrix;
};


#endif //PROYECTOTALLERI_MATRIX_H
