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
  //Choix de la matrice à utiliser

  int userChoiceMat;
  cout << "------------------------------------" << endl;
  cout << "Choississez le matrice : " << endl;
  cout << "1) Matrice aléatoire de la question 1"<< endl;
  cout << "2) Contraintes thermiques pour un transistor" << endl;
  cout << "3) Calcul météo avec le modèle de St Venant" << endl;
  cin >> userChoiceMat;


  int userChoiceSolv; // Choix du solveur à utiliser

// Choix du solveur
	cout << "------------------------------------" << endl;
	cout << "Choississez le solveur : " << endl;
	cout << "1) GPO"<< endl;
	cout << "2) Résidu minimum" << endl;
	cout << "3) GMRes" << endl;
	cout << "4) Résidu minimum conditionné à gauche" << endl;
	cout << "5) Résidu minimum conditionné à droite"<< endl;
	cout << "6) Résidu minimum conditionné à droite flexible"<< endl;
	cout << "7) ..."<< endl;
	cin >> userChoiceSolv;

  switch(userChoiceMat)
	{
    case 1:
    { //Matrice aléatoire de la question 1
      int n;
      cout << "Choisissez la taille de la matrice" << endl;
      cin >> n;
      SparseMatrix<double> An(n,n), C(n,n), In(n,n);
      MatrixXd Bn(n,n), Bnn(n,n);
      VectorXd x0(n), x(n), b(n);

      for (int i=0;i<n;i++)
      {
        x0(i)=i;
      }
      b=x0;

      Bnn=MatrixXd::Random(n,n);
      Bn=(Bnn+MatrixXd::Constant(n,n,1.))/2.;
      C=Bn.sparseView();
      In.setIdentity();
      An=In+1./pow(n,2)*C.transpose()*C;

      switch(userChoiceSolv)
      {
        case 1: //GPO
        GradPasOptimal(An,b,x0,0.000001,1000,x);
        break;

        case 2: //Résidu minimum
        ResMin(An,b,x0,0.00001,1000,x);
        break;

        case 3: //GMRes
        int m;
        cout << "Chosissez m"<< endl;
        cin >> m;
        GMRes(An,b,x0,0.000001,1000,x,m);
        break;

        case 4: //Résidu minimum conditionné à gauche
        ResMin_cond_gauche(An,b,x0,0.00001,1000,x);
        break;

        case 5: //Résidu minimum conditionné à droite
        ResMin_cond_droite(An,b,x0,0.00001,1000,x);
        break;

        default:
        cout << "Ce choix n’est pas possible ! Veuillez recommencer !" << endl;
        exit(0);

      }
      break;
    }

    case 2: // Contraintes thermiques pour un transistor
    {

      //    Définition    //

      SparseMatrix<double> An;
      MatrixXd bn;
      VectorXd x0, x;
      double t1,t2;

      //   Lecture de la matrice   //

      An=Lecture_Matrice_A("/home/valentin/DM/DM-main/DM_Turpault/smt/smt.mtx");
      bn=Lecture_Matrice_b("/home/valentin/DM/DM-main/DM_Turpault/smt/smt_b.mtx");
      x0.resize(bn.size());
      x.resize(bn.size());

      for (int i=0; i<bn.size();i++)
      {
        x0(i)=i;
      }

      cout <<"Calcul ..."<<endl;
      t1 = clock();

      switch(userChoiceSolv)
      {
        case 1: //GPO
        GradPasOptimal(An,bn,x0,0.000001,1000,x);
        break;

        case 2: //Résidu minimum
        ResMin(An,bn,x0,0.00001,1000,x);
        break;

        case 3: //GMRes
        int m;
        cout << "Chosissez m"<< endl;
        cin >> m;
        GMRes(An,bn,x0,0.000001,1000,x,m);
        break;

        case 4: //Résidu minimum conditionné à gauche
        ResMin_cond_gauche(An,bn,x0,0.00001,1000,x);
        break;

        case 5: //Résidu minimum conditionné à droite
        ResMin_cond_droite(An,bn,x0,0.00001,1000,x);
        break;

        default:
        cout << "Ce choix n’est pas possible ! Veuillez recommencer !" << endl;
        exit(0);

      }

      t2=clock();
      cout << "Temps d'execution : "<< (t2 - t1) / CLOCKS_PER_SEC <<" en seconde" <<endl;
      break;
    }

    case 3: // Calcul météo avec le modèle de St Venant
    {

      //    Définition    //

      SparseMatrix<double> An;
      VectorXd bn;
      VectorXd x0, x;
      double t1,t2;

      //   Lecture de la matrice   //

      An=Lecture_Matrice_A_2("/home/valentin/DM/DM-main/DM_Turpault/shallow_water1/shallow_water1.mtx");
      bn.resize(An.rows());
      x0.resize(An.rows());
      x.resize(An.rows());

      for (int i=0; i<An.rows();i++)
      {
        x0(i)=i;
      }

      bn=x0;

      cout <<"Calcul ..."<<endl;
      t1 = clock();

      switch(userChoiceSolv)
      {
        case 1: //GPO
        GradPasOptimal(An,bn,x0,0.000001,1000,x);
        break;

        case 2: //Résidu minimum
        ResMin(An,bn,x0,0.00001,1000,x);
        break;

        case 3: //GMRes
        int m;
        cout << "Chosissez m"<< endl;
        cin >> m;
        GMRes(An,bn,x0,0.000001,1000,x,m);
        break;

        case 4: //Résidu minimum conditionné à gauche
        ResMin_cond_gauche(An,bn,x0,0.00001,1000,x);
        break;

        case 5: //Résidu minimum conditionné à droite
        ResMin_cond_droite(An,bn,x0,0.00001,1000,x);
        break;

        default:
        cout << "Ce choix n’est pas possible ! Veuillez recommencer !" << endl;
        exit(0);

      }

      t2=clock();
      cout << "Temps d'execution : "<< (t2 - t1) / CLOCKS_PER_SEC <<" en seconde" <<endl;
      break;
    }

    default:
    cout << "Ce choix n’est pas possible ! Veuillez recommencer !" << endl;
    exit(0);
  }



  return 0;

}
