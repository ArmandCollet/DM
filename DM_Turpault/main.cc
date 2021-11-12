#include "algo.h"
#include <iostream>
#include <string>
#include <random>
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

  int n; n=10;
  int m;m=5;
  MatrixXd An(n,n);
  MatrixXd Bn(n,n);
  MatrixXd Bnn(n,n);
  MatrixXd In(n,n);
  VectorXd v(n);
  
  v=VectorXd::LinSpaced(n,0.,(double)(n-1));
  Bnn=MatrixXd::Random(n,n);
  Bn=(Bnn+MatrixXd::Constant(n,n,1.))/2.;
  In=MatrixXd::Identity(n,n);
  An=In+1./pow(n,2)*Bn.transpose()*Bn;
  
  vector<Eigen::MatrixXd> HmVm;
  MatrixXd Vm_1(n,m+1);MatrixXd Hm_barre(m+1,m);MatrixXd Vm(n,m);MatrixXd Hm(m,m);
   
  HmVm=Arnoldi(An,v,m);
  Vm_1=HmVm[3]; Hm_barre=HmVm[1]; Vm=HmVm[2]; Hm=HmVm[0];

  cout<<"Norme des vecteurs vi"<<endl;
  for(int i=0;i<m;i++)
    {
      cout<<Vm.col(i).norm()<<endl;
    }

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
	
  
    
  return 0;
  
}
