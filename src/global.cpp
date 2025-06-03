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

void eop19620101(int rows, int cols) {
    eopdata = Matrix(rows, cols);

    FILE *fp = fopen("../data/eop19620101.txt", "r");
    if (fp == NULL) {
        printf("Error opening eop19620101.txt file\n");
        exit(EXIT_FAILURE);
    }

    for (int i = 1; i <= rows; i++) {
        for (int j = 1; j <= cols; j++) {
            if (fscanf(fp, "%lf", &eopdata(i, j)) != 1) {
                printf("Error\n");
                fclose(fp);
                exit(EXIT_FAILURE);
            }
        }
    }

    fclose(fp);
    printf("eopdata charged: %d rows, %d columns\n", rows, cols);
}

void GGM03S(int n) {
    Cnm = Matrix::zeros(n+1, n+1);
    Snm = Matrix::zeros(n+1, n+1);

    FILE *fp = fopen("../data/GGM03S.txt","r");
    if (fp == NULL) {
        printf("Fail open GGM03S.txt file\n");
        exit(EXIT_FAILURE);
    }

    double aux;
    int a, b;
    double cnm, snm;


    while (fscanf(fp, "%d%d%lf%lf%lf%lf", &a, &b, &cnm, &snm, &aux, &aux) == 6) {

        Cnm(a+1, b+1) = cnm;
        Snm(a+1, b+1) = snm;
    }

    fclose(fp);
}


void DE430Coeff(int rows, int cols) {

    PC = Matrix(rows, cols);

    std::ifstream file("../data/DE430Coeff.txt");
    if (!file.is_open()) {
        throw std::runtime_error("Cant open DE430Coeff.txt");
    }

    std::string line;
    int row = 1;
    while (row <= rows && std::getline(file, line)) {
        line.erase(0, line.find_first_not_of(" \t\r\n"));
        line.erase(line.find_last_not_of(" \t\r\n") + 1);

        if (line.empty()) continue;

        std::istringstream iss(line);
        for (int col = 1; col <= cols; ++col) {
            double v;
            if (!(iss >> v)) {
                std::cerr << "Error reading: " << row << ",  " << col << "\n";
                throw std::runtime_error(
                        "Error reading:  " + std::to_string(row) + ", " + std::to_string(col));
            }
            PC(row, col) = v;
        }
        ++row;
    }
    file.close();

    if (PC.getFilas() == 0 || PC.getColumnas() == 0) {
        throw std::runtime_error("DE430Coeff: Matriz PC is empty.");
    }
}
