#include "algo.h"
#include <iostream>
#include <string>
#include <random>
#include <stdio.h>
#include <time.h>
using namespace std;
using namespace Eigen;

int main()
{
  //Verification GPO, ResMin
  // MatrixXd A(4,4);
  // A<<1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.,15.,16.;
  // VectorXd b(4);
  // b<<0.,1.,2.,3.;
  // VectorXd x0(4);
  // x0<<0.,1.,2.,3.;
  // VectorXd x(4);


  // GradPasOptimal(A,b,x0,0.001,1000,x);
  // cout << "Gradient pas optimal donne x="<< endl;
  // cout << x << endl;


  // ResMin(A,b,x0,0.001,100,x);
  // cout<<"Résidu minimium donne x="<<endl;
  // cout<<x<<endl;

  //Verification Arnoldi
/*
  int n; n=100;
  int m;m=1;
  SparseMatrix<double> An(n,n);
  MatrixXd Bn(n,n);
  MatrixXd Bnn(n,n);
  MatrixXd In(n,n);
  MatrixXd Qm(n,n);
  MatrixXd Rm(n,n);
  VectorXd v(n);
  VectorXd x0(n);
  VectorXd x(n);
  VectorXd y(m);
  VectorXd b(n);
  //b=VectorXd::Random(n);
  for (int i=0;i<n;i++)
  {
    x0(i)=i+1;b(i)=i;
  }


  v=VectorXd::LinSpaced(n,0.,(double)(n-1));
  An=sparseMatrix<double>::rsparseMatrix(n,n,nnz=100);
  Bn=(MatrixXd::Random(n,n)+MatrixXd::Constant(n,n,1.))/2.;
  In=MatrixXd::Identity(n,n);
  An=In+1./pow(n,2)*Bn.transpose()*Bn;*/

  //vector<Eigen::MatrixXd> HmVm;
  //MatrixXd Vm_1(n,m+1);MatrixXd Hm_barre(m+1,m);MatrixXd Vm(n,m);MatrixXd Hm(m,m);

  //HmVm=Arnoldi(An,v,m);
  //Vm_1=HmVm[3]; Hm_barre=HmVm[1]; Vm=HmVm[2]; Hm=HmVm[0];

  //cout<<"Norme des vecteurs vi"<<endl;
  //for(int i=0;i<m;i++)
  //{
  //  cout<<Vm.col(i).norm()<<endl;
  //}

  //Givens(Hm_barre, Qm, Rm)

  /*
  cout<<"Am*Vm="<<endl;
  cout<<An*Vm<<endl;
  cout<<"Vm+1*Hm_barre="<<endl;
  cout<<Vm_1*Hm_barre<<endl;
  cout<<"Hm_barre="<<endl;
  cout<<Hm_barre<<endl;
  cout<<"VmT*A*Vm="<<endl;
  cout<<Vm.transpose()*An*Vm<<endl;
  cout<<"Hm="<<endl;
  cout<<Hm<<endl;
  cout<<"----------------------------------------------------------"<<endl;
  cout<<"Givens : "<< endl;
  cout<<"Qm*Rm="<<endl;
  cout<<Qm*Rm <<endl;
  cout<<"Hm"<<endl;
  cout<<Hm<<endl;
  cout <<"QmT * Qm = "<< endl;
  cout<<Qm.transpose()*Qm<<endl;
  cout <<" Rm " << endl;
  cout<< Rm << endl;
  cout<<"Q"<<endl;
  cout<<Qm<<endl;
  cout<<"Rm"<<endl;
  cout<<Rm<<endl;*/



  int n; n=100;
  int m;m=100;
  SparseMatrix<double> An(n,n);
  SparseMatrix<double> C(n,n);
  MatrixXd Bn(n,n);
  MatrixXd Bnn(n,n);
  SparseMatrix <double> In(n,n);
  VectorXd x0(n);
  VectorXd x(n);
  VectorXd b(n);

  for (int i=0;i<n;i++)
  {
    x0(i)=1.;
  }
  b=x0;

  Bnn=MatrixXd::Random(n,n);
  Bn=(Bnn+MatrixXd::Constant(n,n,1.))/2.;
  C=Bn.sparseView();
  In.setIdentity();
  An=In+1./pow(n,2)*C.transpose()*C;


  ResMin(An,b,x0,0.00001,1000,x);
  cout<<"Résidu minimium donne r="<<endl;
  cout<<(An*x-b).norm()<<endl;

/*
  GradPasOptimal(An,b,x0,0.000001,1000,x);
  cout << "Gradient pas optimal donne r="<< endl;
  cout << An*x-b << endl;

  cout<<"Gmres"<<endl;
  GMRes(An,b,x0,0.000001,1000,x,m);
  cout<<"GMRes donne r="<<endl;
  cout<<An*x-b<<endl;*/

  cout<<"ResMin_cond_gauche"<<endl;
  ResMin_cond_gauche(An,b,x0,0.00001,1000,x);
  cout<<"Résidu minimium preconditionné à gauche donne r="<<endl;
  cout<< (An*x-b).norm() <<endl;

  cout<<"ResMin_cond_droite"<<endl;
  ResMin_cond_droite(An,b,x0,0.00001,1000,x);
  cout<<"Résidu minimium preconditionné à droite donne r="<<endl;
  cout<< (An*x-b).norm() <<endl;

  //    Definition    //

  /*SparseMatrix<double> An;
  MatrixXd bn;
  VectorXd x0;
  VectorXd x;
  VectorXd residu;

  int m=1000;
  double t1,t2;


  //    Initialisation    //

  An=Lecture_Matrice_A("/home/armand/DM/DM-main/DM_Turpault/smt/smt.mtx");
  bn=Lecture_Matrice_b("/home/armand/DM/DM-main/DM_Turpault/smt/smt_b.mtx");
  x0.resize(bn.size());
  x.resize(bn.size());
  residu.resize(bn.size());

  for (int i=0; i<bn.size();i++)
    {
      x0(i)=i;
    }


  //    Calcul    //

  cout <<"Calcul ..."<<endl;
  t1 = clock();
  //GradPasOptimal(An,bn,x0,0.00001,1000,x);
  //ResMin(An,bn,x0,0.00001,1000,x);
  ResMin_cond_gauche(An,bn,x0,0.00001,1000,x);
  //GMRes(An,bn,x0,0.00001,1000,x,m);
  t2=clock();

  residu = An*x-bn;
  cout << endl << "Norme du résidu = "<< endl;
  cout << residu.norm() << endl;
  //cout << "GMRes donne x="<<endl;
  //cout << x(25709) << endl;
  cout << "Temps d'execution : "<< (t2 - t1) / CLOCKS_PER_SEC <<" en seconde" <<endl;*/
/*
  //Test Resol_LU
  MatrixXd L(3,3); SparseMatrix<double> L_(3,3);MatrixXd U(3,3); SparseMatrix<double> U_(3,3);
  L<<1.,0.,0.,-0.5,1.,0.,0.,-0.666,1.;
  U<<2.,-1.,0.,0.,1.5,-1.,0.,0.,1.333;
  L_=L.sparseView();  U_=U.sparseView();
  cout<<"L="<<endl<<L_<<endl;cout<<"U="<<endl<<U_<<endl;
  VectorXd x(3);
  VectorXd b(3);b<<1.,2.,3.;cout<<"b="<<b<<endl;
  x=Resol_LU(L_,U_,b);
  cout<<L*U*x-b<<endl;*/
  /*
  SparseMatrix<double> E(4,4); SparseMatrix<double> F(4,4); SparseMatrix<double> D(4,4); SparseMatrix<double> D_1(4,4);

  for (int i=0; i<C.outerSize(); ++i)
  {
    for (SparseMatrix<double>::InnerIterator it(C,i); it; ++it)
    {
      //it.value();
      //it.row();   // row index
      //it.col();   // col index (here it is equal to i)
      //it.index(); // inner index, here it is equal to it.row()

      if (it.row()==it.col())
      {
        D.coeffRef(it.row(),it.col()) = it.value();
        D_1.coeffRef(it.row(),it.col()) = 1./it.value();
      }
      if (it.row()>it.col())
      {
        E.coeffRef(it.row(),it.col()) = it.value();
      }

      if (it.row()<it.col())
      {
        F.coeffRef(it.row(),it.col()) = it.value();
      }
    }
  }
    cout<<"D="<<endl<<D<<endl;cout<<"E="<<endl<<E<<endl;cout<<"F="<<endl<<F<<endl;cout<<"D-1="<<endl<<D_1<<endl;*/

  return 0;

}
