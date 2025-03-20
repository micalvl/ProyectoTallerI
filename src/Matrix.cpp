//
// Created by micalvl on 20/03/2025.
//

#include "../include/Matrix.h"
#include <iostream>
#include <iomanip>

Matrix::Matrix(int fil, int col) : fil(fil), col(col){
    initMatrix();
}

void Matrix::initMatrix(){
    matrix = new double*[fil];
    for (int i = 0; i< fil; i++){
        matrix[i] = new double[col];
    }

    for(int i = 0; i < fil; i++){
        for(int j = 0; j < col; j++){
            matrix[i][j] = 0.0;
        }
    }
}

Matrix::Matrix(int fil, int col, double v[], int n) : fil(fil), col(col){
    initMatrix();

    int k = 0;

    for(int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++) {
            if (k < n)
                matrix[i][j] = v[k++];
            else
                matrix[i][j] = 0;
        }

}

