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


    Matrix();


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


    /**
    * @brief Computes Frobenius norm (square root of sum of squares of all elements).
    * @return Frobenius norm.
    */
    double norm() const;


    /**
    * @brief Returns number of rows.
    * @return Number of rows (fil).
    */
    int getFilas() const;


    /**
    * @brief Returns number of columns.
    * @return Number of columns (col).
    */
    int getColumnas() const;


    /**
    * @brief Transposes the matrix (swap rows and columns).
    * @return Transposed matrix of dimensions (col×fil).
    */
    Matrix transpose() const;


    /**
    * @brief Creates an identity matrix of given size.
    * @param[in] size  Number of rows and columns.
    * @return Identity matrix (size×size).
    */
    static Matrix identity(int size);


    /**
    * @brief Multiplies every element by a scalar.
    * @param[in] scalar  Value to multiply each element.
    * @return Scaled matrix.
    */
    Matrix opsc(double scalar) const;


    /**
    * @brief Divides every element by a scalar.
    * @param[in] scalar  Value to divide each element.
    * @return Scaled matrix.
    */
    Matrix divsc(double scalar) const;


    /**
    * @brief Prints a raw float** array of size n×m (for debugging).
    * @param[in] arr  2D float array.
    * @param[in] n    Number of rows.
    * @param[in] m    Number of columns.
    */
    void printMatrixValues(float** arr, int n, int m);


    /**
    * @brief Prints the right half (inverse) of a Gauss-Jordan augmented matrix.
    * @param[in] arr  2D float array (augmented 2n×n).
    * @param[in] n    Order of original matrix.
    * @param[in] m    Total number of columns (2n).
    */
    void printInverseMatrix(float** arr, int n, int m);

    void findInvMatGaussJordan(float** mat, int order);


    /**
    * @brief Computes the inverse of this matrix using Gauss-Jordan elimination.
    * @return Inverted matrix of same dimensions.
    */
    Matrix inverse();


    /**
    * @brief Creates a zero matrix of given dimensions.
    * @param[in] rows  Number of rows.
    * @param[in] cols  Number of columns.
    * @return Zero matrix (rows×cols).
    */
    static Matrix zeros(int rows, int cols);


    /**
    * @brief Computes the 3D cross product of two 3×1 vectors.
    * @param[in] a  First 3×1 vector.
    * @param[in] b  Second 3×1 vector.
    * @return 3×1 cross product vector.
    */
    static Matrix cross(const Matrix& a, const Matrix& b);


    /**
    * @brief Computes the 3D dot product of two 3×1 vectors.
    * @param[in] a  First 3×1 vector.
    * @param[in] b  Second 3×1 vector.
    * @return Scalar dot product.
    */
    static double dot(const Matrix& a, const Matrix& b);



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
