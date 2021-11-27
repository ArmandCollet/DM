#include <iostream>
#include <fstream>
#include "algo.h"
#include <vector>
#include <tuple>
using namespace std;
using namespace Eigen;

//CSR (inutile à priori)
std::tuple<vector<double>,vector<int>,vector<int>> stockageCSR(const MatrixXd A)
{
  vector<double> AA; vector<int> JA; vector<int> IA;
  int k;k=0;
  IA.push_back(1.);

  for (int i=0;i<A.rows();i++)
  {
    for (int j=0;j<A.cols();j++)
    {
      if (A(i,j)!=0.)
      {
        k+=1;
        AA.push_back(A(i,j));
        JA.push_back(j+1);
      }
    }
    IA.push_back(k+1);
  }
  return std::make_tuple(AA,JA,IA);
}
void MatVecCSR(const std::tuple<std::vector<double>, std::vector<int>, std::vector<int>> Tab, const Eigen::VectorXd v)
{
  int n; n=v.size();
  VectorXd y; y.resize(n); y=VectorXd::Constant(n,0.);
  vector<double> AA; AA=get<0>(Tab);
  vector<int> JA; JA=get<1>(Tab);
  vector<int> IA; IA=get<2>(Tab);

  for (int i=0;i<n;i++)
  {
    for (int k=IA[i]; k<IA[i+1];k++)
    {
      y(i)+=AA[k-1]*v(JA[k-1]-1);
    }
    cout<<y(i)<<endl;
  }
}
//Fin CSR

// Gradient à pas optimal
void GradPasOptimal (const SparseMatrix<double> A, const VectorXd b, const VectorXd x0, const double epsilon, const int kmax, VectorXd & x)
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
    //cout << r.norm() << endl;
  }
  cout<<"Nombre d'itérations ="<<k<<endl;
  if (k>kmax)
  {
    cout<<"Tolérance non atteinte: "<<endl;
  }

}

//Résidu minimium
void ResMin (const SparseMatrix<double> A,const VectorXd b, const VectorXd x0,const double epsilon,const int kmax, VectorXd & x)

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
  cout<<"Nombre d'itérations ="<<k<<endl;
  if (k>kmax)
  {
    cout<<"Tolérance non atteinte: "<<endl;
  }
}

//Résolution à partir d'une Décomposition LU
VectorXd Resol_LU(SparseMatrix <double> M, VectorXd b)
{
  int n=M.rows();
  VectorXd x(n), y(n);
  double s1, s2;

  y(0) = b(0);
  for(int i=1; i<n; i++)
  {
    s1=0.;cout<<"test"<<endl;
    for(int k=0; k<i; k++)
    {
      s1 += M.coeffRef(i,k)*y(k);
    }
    y(i)=b(i)-s1; cout<<y(i)<<endl;
  }
  x(n-1)=y(n-1)/M.coeffRef(n-1,n-1);


  for (int i=n-2; i>=0; i--)
  {
    s2=0.;  cout<<"test1"<<endl;

    for (int k=i+1;k<n;k++)
    {
      s2 += M.coeffRef(i,k)*x(k);
    }
    x(i) = (y(i)-s2)/M.coeffRef(i,i);
  }
  return x;
}

//Résidu minimum préconditionné à gauche
void ResMin_cond_gauche(Eigen::SparseMatrix<double> A, const Eigen::VectorXd b, const Eigen::VectorXd x0, const double epsilon, const int kmax, Eigen::VectorXd & x)
{
  int k=0; int m=A.rows();
  VectorXd r(b.size()); r=b-A*x0;
  VectorXd z(b.size()); VectorXd q(m); VectorXd w(m);
  double alpha;
  SparseMatrix<double> M(m,m); SparseMatrix<double> E(m,m); SparseMatrix<double> F(m,m); SparseMatrix<double> D(m,m); SparseMatrix<double> D_1(m,m);

  for (int i=0;i<m;i++)
  {
    for (SparseMatrix<double>::InnerIterator it(A,i);it;++it)
    {
      D.coeffRef(i,i) = it.value(); D_1.coeffRef(i,i)=1./it.value();

      for (int j=0;j<i;j++)
      {
        E.coeffRef(i,j)=it.value();
      }

      for (int j=i+1;j<m;j++)
      {
        F.coeffRef(i,j)=it.value();
      }
    }

  }

  M=(D-E)*D_1*(D-F);

  x=x0;

  q=Resol_LU(M,r);

  while ((r.norm()>epsilon) && (k<=kmax))
  {
    w=A*q;
    z=Resol_LU(M,w);
    alpha=q.dot(z)/z.dot(z);
    x = x +alpha*q;
    r = r - alpha*w;
    q=q-alpha*z;
    k+=1;
  }
  cout<<"Nombre d'itérations ="<<k<<endl;
  if (k>kmax)
  {
    cout<<"Tolérance non atteinte: "<<endl;
  }
}



//Arnoldi
std::vector<Eigen::MatrixXd> Arnoldi(const Eigen::SparseMatrix<double> A, Eigen::VectorXd & v, const int m)
{
  //Déclaration des variables

  MatrixXd Vm_1(v.size(),m+1); Vm_1 = MatrixXd::Constant(v.size(),m+1);
  MatrixXd Vm(v.size(),m);
  MatrixXd Hm(m,m);
  MatrixXd Hm_barre(m+1,m); Hm_barre = MatrixXd::Constant(m+1,m);
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

/*
//Givens avec produit de matrices honteux
Eigen::MatrixXd  Givens_Rotation(const Eigen::MatrixXd M, int i, int j)
{
  // Constitution de la matrice de rotation

  Eigen::MatrixXd Rij(max(M.rows(), M.cols()),max(M.rows(), M.cols()));
  MatrixXd In(max(M.rows(), M.cols()),max(M.rows(), M.cols()));

  In  = MatrixXd::Identity(max(M.cols(), M.rows()),max(M.cols(), M.rows()));
  Rij = In;

  double a = sqrt( pow(M(j,j),2) + pow(M(i,j),2) );
  double c = M(j,j)/a;
  double s = -M(i,j)/a;
  Rij(j,j) = c;
  Rij(i,i) = c;
  Rij(i,j) = s;
  Rij(j,i) = -s;

  return Rij;
}

void Givens(const Eigen::MatrixXd A, Eigen::MatrixXd & Q, Eigen::MatrixXd & R)

//Décomposition QR Givens
{
  Eigen::MatrixXd Rij(max(A.rows(), A.cols()),max(A.cols(), A.rows()));
  MatrixXd In(max(A.rows(), A.cols()),max(A.rows(), A.cols()));

  In = MatrixXd::Identity(max(A.cols(), A.rows()),max(A.cols(), A.rows()));
  Q = In;
  R = A;

  for (int i = 0; i<A.rows(); i++)
  {
    for (int j = 0; j<i; j++)
    {
      Rij = Givens_Rotation(R, i, j);
      Q=Q*Rij.transpose();
      R=Rij*R;
    }
  }
}*/

//Givens honnête
void GivensOpt(const Eigen::MatrixXd A, Eigen::MatrixXd & Q, Eigen::MatrixXd & R)
{
  int m; m=A.rows();int n=A.cols();
  R=A;Q=MatrixXd::Identity(m,m);

  for (int i=0;i<m;i++)
  {

    for (int j=0; j<i;j++)
    {
      double a = sqrt( pow(R(j,j),2) + pow(R(i,j),2) );
      double c = R(j,j)/a;
      double s = -R(i,j)/a;

      for (int l=0; l<m ; l++)
      {
        if (l<n)
        {
          double Ril_inter = c*R(i,l) + s*R(j,l);
          R(j,l) = c*R(j,l) - s*R(i,l);
          R(i,l) = Ril_inter;
        }
      }

      for (int l=0;l<m;l++)
      {
          double Qil_inter = s*Q(l,j)+c*Q(l,i);
          Q(l,j) = -s*Q(l,i)+c*Q(l,j);
          Q(l,i) = Qil_inter;
      }
    }
  }
}

//Résolution système triangulaire supérieur
void resol_syst_triang_sup(const Eigen::MatrixXd A, Eigen::VectorXd & y, const Eigen::VectorXd b)
{
  double S;
  int m; m=A.rows();

  for (int i=m-1; i>=0;i--)
  {
    S=0.;
    for (int j=m-1;j>i;j--)
    {
      S+=A(i,j)*y(j);
    }
    y(i)=(b(i)-S)/A(i,i);
  }

}

//GMRes
void GMRes(const Eigen::SparseMatrix<double> A, const Eigen::VectorXd b, const Eigen::VectorXd x0, const double epsilon, const int kmax, Eigen::VectorXd & x, const int m)
{
  //Initialisation
  int k=0;
  VectorXd r(b.size()); r=b-A*x0;
  VectorXd y(m);
  double beta; beta=r.norm(); double gamma;
  MatrixXd Hm(m,m); MatrixXd Vm(b.size(),m);MatrixXd Hm_barre(m+1,m);MatrixXd Vm_1(b.size(),m+1);
  vector<Eigen::MatrixXd> HmVm;
  Eigen::MatrixXd Qm(m,m); Eigen::MatrixXd Rm(m,m); Eigen::MatrixXd Qm_barre(m+1,m+1);Eigen::MatrixXd Rm_barre(m+1,m);
  Eigen::VectorXd gm_barre(m+1);Eigen::VectorXd gm(m);

  x=x0;

  while ((beta>epsilon) && (k<=kmax))
  {
    HmVm = Arnoldi(A,r,m);
    Hm = HmVm[0]; Vm = HmVm[2]; Hm_barre = HmVm[1]; Vm_1 = HmVm[3];
    //cout<<"VmT*Vm="<<endl;
    //cout<<Vm.transpose()*Vm<<endl;

    //Givens(Hm,Qm,Rm);
    GivensOpt(Hm_barre,Qm_barre,Rm_barre);

    gm_barre = beta*(Qm_barre.transpose()).col(0);
    gamma = gm_barre[m];

    for(int i=0;i<m;i++)
    {
      Rm.row(i) = Rm_barre.row(i); gm.row(i) = gm_barre.row(i);
    }

    resol_syst_triang_sup(Rm, y, gm);

    x += Vm * y;

    r = b-A*x;
    //r = gamma*Vm_1*((Qm_barre.transpose()).col(m));

    //beta = abs(gamma);
    beta = r.norm();

    k+=1;
  }
  cout<<"Nombre d'itérations ="<<k<<endl;
  if (k>kmax)
  {
    cout<<"Tolérance non atteinte: "<<endl;
  }

}



/*--------------------LECTURE MATRICE-----------------------*/
SparseMatrix <double>  Lecture_Matrice_A (string fichier)// pour les matrice du 1)
{
  cout << "Lecture du fichier : "<<fichier<<endl;
  ifstream flux_A(fichier);
  SparseMatrix <double> A;
  int n,m,kmax;
  string ligne;
  if (flux_A)
    {
      for (int i=0; i<34 ; i++)
	{
	  getline(flux_A,ligne);
	}
      flux_A >> n;
      flux_A >> m;
      flux_A >> kmax;
      cout << "n, m, imax : "<<n<<" "<<m<<" "<<kmax<<endl<<endl;
      A.resize(n,m);
      vector<Triplet<double>> triplets;
      double a(0);
      int i(0);
      int j(0);
      for (int k=0; k<kmax; k++)
	{

	  flux_A >> i >> j >> a ;
	  triplets.push_back({i-1,j-1,a});
	  triplets.push_back({j-1,i-1,a});
	}
      A.setFromTriplets(triplets.begin(), triplets.end());
    }
  else
    {
      cout << "Ce fichier n'existe pas"<<endl;
    }
  flux_A.close();
  return A;
}



VectorXd Lecture_Matrice_b(string fichier)
{
  cout << "Lecture du fichier : "<<fichier<<endl;
  ifstream flux_b(fichier); //../smt/smt_b.mtx
  VectorXd b;
  int kmax;
  double a;
  string ligne;
  if (flux_b)
    {
      for (int i=0; i<6 ; i++)
	{
	  getline(flux_b,ligne);
	}
      flux_b >> kmax >> a; //a va être réécrit
      cout << "kmax : "<<kmax<<endl<<endl;
      b.resize(kmax);
      for (int k=1; k<kmax+1 ; k++)
	{
	  flux_b >> a;
	  b[k-1] = a;
	}
    }
  else
    {
      cout<<"Ce fichier n'existe pas"<<endl;
    }
  flux_b.close();
  return b;
}