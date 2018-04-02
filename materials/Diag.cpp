//
// Created by Utena on 2018/3/29.
//
#include <fstream>
#include "cmath"
//#include <Eigen/Dense>

int main(){
    double Mass[3]={1.008,15.999,1.008};
    double Hessian_Raw[9][9];
    double Hessian[9][9];
    std::ifstream infile;
    infile.open("hess.dat");
    int i=0;int j=0;
    for (int i=0;i<9;i++) {
        for (int j=0; j < 9; j++){
            infile >> Hessian_Raw[i][j];}
    }
    for (i=0; i<9;i++) {
        for (j=0; j < 9; j++){
            Hessian[i][j]=Hessian_Raw[i][j]/(sqrt(Mass[i/3]*Mass[j/3]));}
    }
    std::ofstream outfile;
    outfile.open("Hessian.dat");
    for (i=0; i<9;i++){
    for (j=0;j<9;j++){
        outfile<<Hessian[i][j]<<"  ";}
        outfile<<std::endl;
    }
};
