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

void DE430Coeff(int rows, int cols) {
    PC = Matrix::zeros(rows, cols);
    std::ifstream file("../data/DE430Coeff.txt");
    if (!file.is_open()) {
        std::cerr << "Error abriendo DE430Coeff.txt\n";
        std::exit(EXIT_FAILURE);
    }

    std::string line;
    for (int i = 1; i <= rows; ++i) {
        if (!std::getline(file, line)) {
            std::cerr << "Línea " << i << " no encontrada\n";
            std::exit(EXIT_FAILURE);
        }
        // sustituir comas por puntos
        std::replace(line.begin(), line.end(), ',', '.');

        // tokenizar por espacios/tabs
        std::istringstream ss(line);
        for (int j = 1; j <= cols; ++j) {
            std::string tok;
            if (!(ss >> tok)) {
                std::cerr << "Token vacío en PC("<<i<<","<<j<<")\n";
                std::exit(EXIT_FAILURE);
            }

            // filtrar caracteres inválidos
            std::string f;
            for (char c : tok) {
                if (esValidoEnNumero(c)) f.push_back(c);
            }

            // separar mantisa/exponente
            auto pE = f.find_first_of("Ee");
            std::string mant = f, expn;
            if (pE != std::string::npos) {
                expn = f.substr(pE);
                mant = f.substr(0, pE);
            }

            // si la mantisa tiene >1 '.', borrar todos (eran miles)
            if (std::count(mant.begin(), mant.end(), '.') > 1)
                mant.erase(std::remove(mant.begin(), mant.end(), '.'), mant.end());

            std::string clean = mant + expn;
            double val;
            try {
                val = std::stod(clean);
            } catch (...) {
                std::cerr << "stod fallo en PC("<<i<<","<<j<<") con '"<<clean<<"'\n";
                std::exit(EXIT_FAILURE);
            }

            // —> Escalar solo las dos primeras columnas a JD:
            if (j == 1 || j == 2) {
                val /= 1e8;  // 243326450000000 → 2433264.5
            }

            PC(i, j) = val;
        }
    }
}
