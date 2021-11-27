#include <Dense>
#include<tuple>
#include<vector>
#include<iostream>
#include<Sparse>

void GradPasOptimal(const Eigen::SparseMatrix<double> A, const Eigen::VectorXd b, const Eigen::VectorXd x0, const double epsilon,const int kmax, Eigen:: VectorXd & x);

void ResMin(const Eigen::SparseMatrix<double> A, const Eigen::VectorXd b, const Eigen::VectorXd x0, const double epsilon, const int kmax, Eigen::VectorXd & x);

/*Eigen::MatrixXd ArnoldiH(const Eigen::MatrixXd A, Eigen::VectorXd & v);
  std::vector<Eigen::VectorXd> ArnoldiV(const Eigen::MatrixXd A, Eigen::VectorXd & v);*/

std::vector<Eigen::MatrixXd> Arnoldi(const Eigen::SparseMatrix<double> A, Eigen::VectorXd & v, const int m);


//Eigen::MatrixXd Givens_Rotation(const Eigen::MatrixXd M, int i, int j);

//void Givens(const Eigen::MatrixXd A, Eigen::MatrixXd & Qm, Eigen::MatrixXd & Rm);

void GivensOpt(const Eigen::MatrixXd A, Eigen::MatrixXd & Q, Eigen::MatrixXd & R);

void GMRes(const Eigen::SparseMatrix<double> A, const Eigen::VectorXd b, const Eigen::VectorXd x0, const double epsilon, const int kmax, Eigen::VectorXd & x, const int m);

void resol_syst_triang_sup(const Eigen::MatrixXd A, Eigen::VectorXd & y, const Eigen::VectorXd b);

std::tuple<std::vector<double>,std::vector<int>,std::vector<int>> stockageCSR(const Eigen::MatrixXd A);

void MatVecCSR(const std::tuple<std::vector<double>, std::vector<int>, std::vector<int>> Tab, const Eigen::VectorXd v);


Eigen::SparseMatrix<double> Lecture_Matrice_A(std::string fichier);
Eigen::VectorXd Lecture_Matrice_b(std::string fichier);

void ResMin_cond_gauche(Eigen::SparseMatrix<double> A, const Eigen::VectorXd b, const Eigen::VectorXd x0, const double epsilon, const int kmax, Eigen::VectorXd & x);

Eigen::VectorXd Resol_LU(Eigen::SparseMatrix<double> M, Eigen::VectorXd b);
