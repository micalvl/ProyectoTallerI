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

void eop19620101(int nfilas, int ncols) {
    eopdata = Matrix(nfilas, ncols);

    FILE *fp = fopen("../data/eop19620101.txt", "r");
    if (fp == NULL) {
        printf("Fail open eop19620101.txt file\n");
        exit(EXIT_FAILURE);
    }

    // Leer línea a línea (por filas)
    for (int i = 1; i <= nfilas; i++) {
        for (int j = 1; j <= ncols; j++) {
            if (fscanf(fp, "%lf", &eopdata(i, j)) != 1) {
                printf("Error leyendo dato (%d,%d) del fichero eop19620101.txt\n", i, j);
                fclose(fp);
                exit(EXIT_FAILURE);
            }
        }
    }

    fclose(fp);
    printf("[OK] eopdata cargado: %d filas, %d columnas\n", nfilas, ncols);
}

void GGM03S(int nmax) {
    Cnm = Matrix::zeros(nmax+1, nmax+1);
    Snm = Matrix::zeros(nmax+1, nmax+1);

    FILE *fp = fopen("../data/GGM03S.txt","r");
    if (fp == NULL) {
        printf("Fail open GGM03S.txt file\n");
        exit(EXIT_FAILURE);
    }

    double aux;
    int n, m;
    double cnm, snm;


    while (fscanf(fp, "%d%d%lf%lf%lf%lf", &n, &m, &cnm, &snm, &aux, &aux) == 6) {

        Cnm(n+1, m+1) = cnm;
        Snm(n+1, m+1) = snm;
    }

    fclose(fp);
}


static bool esValidoEnNumero(char c) {
    return std::isdigit(c) || c=='+'||c=='-'||c=='.'||c=='E'||c=='e';
}

void DE430Coeff(int expectedRows, int expectedCols) {

    PC = Matrix(expectedRows, expectedCols);

    std::ifstream file("../data/DE430Coeff.txt");
    if (!file.is_open()) {
        throw std::runtime_error("No pude abrir DE430Coeff.txt");
    }

    std::string line;
    int row = 1;
    while (row <= expectedRows && std::getline(file, line)) {
        line.erase(0, line.find_first_not_of(" \t\r\n"));
        line.erase(line.find_last_not_of(" \t\r\n") + 1);

        if (line.empty()) continue;

        std::istringstream iss(line);
        for (int col = 1; col <= expectedCols; ++col) {
            double v;
            if (!(iss >> v)) {
                std::cerr << "Error leyendo valor en fila " << row << ", columna " << col << "\n";
                throw std::runtime_error(
                        "Error leyendo valor en fila " + std::to_string(row) + ", columna " + std::to_string(col));
            }
            PC(row, col) = v;
        }
        ++row;
    }
    file.close();

    if (row <= expectedRows) {
        std::cerr << "[ADVERTENCIA] El archivo tiene menos filas (" << row-1 << ") de las esperadas (" << expectedRows << ")\n";
    }

    if (PC.getFilas() == 0 || PC.getColumnas() == 0) {
        throw std::runtime_error("DE430Coeff: Matriz PC quedó vacía.");
    }
}
