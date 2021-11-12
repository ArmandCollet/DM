#include <iostream>
#include "algo.h"
#include <vector>
#include <tuple>
using namespace std;
using namespace Eigen;

void GradPasOptimal (const MatrixXd A, const VectorXd b, const VectorXd x0, const double epsilon, const int kmax, VectorXd & x)
{
  //Initialisation
  int k=0;
  VectorXd r(b.size()); r=b-A*x0;
  VectorXd z(b.size());
  double alpha;
  x=x0;
  //Boucle
  while((r.norm()>epsilon) && (k<=kmax))
    {
      z =A*r;
      alpha = r.dot(r) / z.dot(r);
      x = x + alpha*r;
      r = r - alpha*z;
      k+=1;
    }

  if (k>kmax)
    {
      cout<<"Tolérance non atteinte: "<<endl;
    }
  
}

void ResMin (const MatrixXd A,const VectorXd b, const VectorXd x0,const double epsilon,const int kmax, VectorXd & x)

{
  //Initialisation
  int k=0;
  VectorXd r(b.size()); r=b-A*x0;
  VectorXd z(b.size());
  double alpha;
  x=x0;
  //Boucle
  while ((r.norm()>epsilon) && (k<=kmax))
    {
      z = A*r;
      alpha = r.dot(z) / z.dot(z);
      x = x + alpha*r;
      r = r - alpha*z;
      k+=1;
    }

  if (k>kmax)
    {
      cout<<"Tolérance non atteinte: "<<endl;
      }
}

/*
MatrixXd ArnoldiH(const Eigen::MatrixXd A, Eigen::VectorXd & v)
{
  int m;
  m=v.size();
  //vector <Eigen::VectorXd> Vm(m);
  MatrixXd Vm(m,m);
  MatrixXd Hm(m,m);
  //vector<Eigen::VectorXd> w(m);
  MatrixXd w(m,m);
  Vm.col(0)=v/v.norm();

  for (int j=0; j<m;j++)
    {
      w.col(j)=A*Vm.col(j);
      for (int i=0;i<j;i++)
	{
	  Hm(i,j)=w.col(j).dot(Vm.col(i));
	  w.col(j)=w.col(j)-Hm(i,j)*Vm.col(i);
	}
      Hm(j+1,j)=(w.col(j)).norm();
      if (Hm(j+1,j)==0.)
	break;
      Vm.col(j+1)=w.col(j)/Hm(j+1,j);
	
    }    
  return Hm;
  }

vector<Eigen::VectorXd> ArnoldiV(const Eigen::MatrixXd A, Eigen::VectorXd & v)
{
  int m;
  m=v.size();
  vector< Eigen::VectorXd> Vm(m);
  MatrixXd Hm(m,m);
  vector<Eigen::VectorXd> w(m); 
  Vm[0]=v/v.norm();

  for (int j=0; j<m;j++)
    {
      w[j]=A*Vm[j];
      for (int i=0;i<j;i++)
	{
	  Hm(i,j)=w[j].dot(Vm[i]);
	  w[j]=w[j]-Hm(i,j)*Vm[i];
	}
      Hm(j+1,j)=w[j].norm();
      if (Hm(j+1,j)==0.)
	break;
      Vm[j+1]=w[j]/Hm(j+1,j);
	
    }    
  return Vm;
  }*/



std::vector<Eigen::MatrixXd> Arnoldi(const Eigen::MatrixXd A, Eigen::VectorXd & v, const int m)
{
  //Déclaration des variables
  
  MatrixXd Vm_1(v.size(),m+1);
  MatrixXd Vm(v.size(),m);
  MatrixXd Hm(m,m);
  MatrixXd Hm_barre(m+1,m);
  MatrixXd w(v.size(),m);
  vector<Eigen::MatrixXd> HmVm;

  //Initialisation
  Vm_1.col(0)=v/v.norm();

  //Boucles
  for (int j=0; j<m;j++)
    {
      w.col(j)=A*Vm_1.col(j);
      for (int i=0;i<=j;i++)
	{
	  Hm_barre(i,j)=w.col(j).dot(Vm_1.col(i));
	  w.col(j)=w.col(j)-Hm_barre(i,j)*Vm_1.col(i);
	}
      Hm_barre(j+1,j)=(w.col(j)).norm();
      if (Hm_barre(j+1,j)==0.)
	break;
      Vm_1.col(j+1)=w.col(j)/Hm_barre(j+1,j);
	
    }

  //Extraction de Vm et Hm à partie de Vm+1 et Hm_barre
  for(int i=0;i<m;i++)
    {
      Vm.col(i)=Vm_1.col(i);Hm.row(i)=Hm_barre.row(i);
    }

  //Remplissage du vecteur contenant Hm, Hm_barre, Vm et Vm+1
  HmVm.push_back(Hm);
  HmVm.push_back(Hm_barre);
  HmVm.push_back(Vm);
  HmVm.push_back(Vm_1);

  return HmVm;
}


void GMRes(const Eigen::MatrixXd A, const Eigen::VectorXd b, const Eigen::VectorXd x0, const double epsilon, const int kmax, Eigen::VectorXd & x, const int m)
{
  //Initialisation
  int k=0;
  VectorXd r(b.size()); r=b-A*x0;
  double beta; beta=r.norm();
  MatrixXd Hm(m+1,m); MatrixXd Vm(A.cols(),m);
  vector<Eigen::MatrixXd> HmVm;
  
  while ((beta>epsilon) && (k<=kmax))
    {
      HmVm=Arnoldi(A,r,m);
      Hm=HmVm[0];Vm=HmVm[1];
      
     
     
    }
}
	       


