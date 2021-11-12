#include <Dense>
#include<tuple>
#include<vector>
#include<iostream>

void GradPasOptimal(const Eigen::MatrixXd A, const Eigen::VectorXd b, const Eigen::VectorXd x0, const double epsilon,const int kmax, Eigen:: VectorXd & x);

void ResMin(const Eigen::MatrixXd A, const Eigen::VectorXd b, const Eigen::VectorXd x0, const double epsilon, const int kmax, Eigen::VectorXd & x);

/*Eigen::MatrixXd ArnoldiH(const Eigen::MatrixXd A, Eigen::VectorXd & v);
  std::vector<Eigen::VectorXd> ArnoldiV(const Eigen::MatrixXd A, Eigen::VectorXd & v);*/

std::vector<Eigen::MatrixXd> Arnoldi(const Eigen::MatrixXd A, Eigen::VectorXd & v, const int m);

void GMRes(const Eigen::MatrixXd A, const Eigen::VectorXd b, const Eigen::VectorXd x0, const double epsilon, const int kmax, Eigen::VectorXd & x, const int m);


