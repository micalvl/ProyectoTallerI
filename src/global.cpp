#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <sstream>
#include "..\include\global.h"
#include <cctype>

using namespace std;

Param AuxParam;
Matrix eopdata;
Matrix Cnm;
Matrix Snm;
Matrix PC;

void eop19620101(int c){
    eopdata=Matrix::zeros(13,c);

    FILE *fp = fopen("../data/eop19620101.txt","r");
    if(fp == NULL){
        printf("Fail open eop19620101.txt file\n");
        exit(EXIT_FAILURE);
    }

    for(int j=1;j<=c;j++){
        fscanf(fp,"%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",&(eopdata(1,j)),
               &(eopdata(2,j)),&(eopdata(3,j)),&(eopdata(4,j)),&(eopdata(5,j)),
               &(eopdata(6,j)),&(eopdata(7,j)),&(eopdata(8,j)),&(eopdata(9,j)),
               &(eopdata(10,j)),&(eopdata(11,j)),&(eopdata(12,j)),&(eopdata(13,j)));
    }

    fclose(fp);
}

void GGM03S(int n){
    Cnm=Matrix::zeros(n,n);
    Snm=Matrix::zeros(n,n);

    FILE *fp = fopen("../data/GGM03S.txt","r");
    if(fp == NULL){
        printf("Fail open GGM03S.txt file\n");
        exit(EXIT_FAILURE);
    }

    double aux;
    for(int i=1;i<=n;i++){
        for (int j=1;j<=i;j++){
            fscanf(fp,"%lf%lf%lf%lf%lf%lf",&aux,&aux,&(Cnm(i,j)),&(Snm(i,j)),&aux,&aux);
        }
    }

    fclose(fp);
}


static bool esValidoEnNumero(char c) {
    return std::isdigit(c) || c=='+'||c=='-'||c=='.'||c=='E'||c=='e';
}

void DE430Coeff(int expectedRows, int expectedCols) {
    PC = Matrix(expectedRows, expectedCols);
    ifstream file("../data/DE430Coeff.txt");
    if (!file.is_open()) {
        throw runtime_error("No pude abrir DE430Coeff.txt");
    }

    string line;
    int row = 1;
    while (row <= expectedRows && getline(file, line)) {
        if (line.empty()) continue;
        istringstream iss(line);
        for (int col = 1; col <= expectedCols; ++col) {
            double v;
            if (!(iss >> v)) {
                throw runtime_error(
                        "Error leyendo valor en fila " + to_string(row) + ", columna " + to_string(col));
            }
            PC(row, col) = v;
        }
        ++row;
    }
}
