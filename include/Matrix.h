/**
 *  @file   Matrix.h
 *  @brief  Matrix's methods
 *  @author [Original Author]
 *  @transcribed by Miguel Calvo León
 *  @date   2025-03-27
 ***********************************************/


#ifndef PROYECTOTALLERI_MATRIX_H
#define PROYECTOTALLERI_MATRIX_H


class Matrix {
public:
    /**
     * @brief Constructs a matrix with the given dimensions.
     *
     * @param fil Number of rows.
     * @param col Number of columns.
     */
    Matrix(int fil, int col);

    /**
     * @brief Constructs a matrix with the given dimensions and initializes it with values from an array.
     *
     * @param fil Number of rows.
     * @param col Number of columns.
     * @param v Array containing the values to initialize the matrix.
     * @param n Number of elements in the array.
     */
    Matrix(int fil, int col, double v[], int n);

    /**
     * @brief Copy constructor. Creates a new matrix as a copy of an existing one.
     *
     * @param m The matrix to copy.
     */
    Matrix(const Matrix& m);

    /**
     * @brief Destructor.
     */
    ~Matrix();

    /**
     * @brief Assigns the contents of another matrix.
     *
     * @param matrix2 The matrix to assign.
     * @return Reference to the modified matrix.
     */
    Matrix& operator=(const Matrix& matrix2);

    /**
     * @brief Sums two matrices.
     *
     * @param matrix2 The matrix to add.
     * @return The resulting matrix.
     */
    Matrix operator+(const Matrix& matrix2) const;

    /**
     * @brief Subtracts another matrix from this matrix.
     *
     * @param matrix2 The matrix to subtract.
     * @return The resulting matrix.
     */
    Matrix operator-(const Matrix& matrix2) const;

    /**
     * @brief Multiplies two matrices.
     *
     * @param matrix2 The matrix to multiply.
     * @return The resulting matrix.
     */
    Matrix operator*(const Matrix& matrix2) const;

    /**
     * @brief Accesses an element of the matrix.
     *
     * @param i Row index.
     * @param j Column index.
     * @return Reference to the matrix element at (i, j).
     */
    double& operator()(const int i, const int j) const;

    /**
     * @brief Prints the matrix to the console.
     */
    void print();

    double norm() const;

    int getFilas() const;
    int getColumnas() const;
    Matrix transpose() const;

    static Matrix identity(int size);

    Matrix opsc(double scalar) const;
    void printMatrixValues(float** arr, int n, int m);
    void printInverseMatrix(float** arr, int n, int m);
    void findInvMatGaussJordan(float** mat, int order);
    Matrix inverse();


private:
    /**
     * @brief Initializes the matrix with default values.
     */
    void initMatrix();

private:
    int fil, col; // Number of rows and columns.
    double** matrix; // Pointer to the matrix
};


#endif //PROYECTOTALLERI_MATRIX_H
