//Code inutile mais à garder au cas ou


//-------------------------test matrice aléatoire


// int n; n=100;
// int m; m=100;
// SparseMatrix<double> An(n,n), C(n,n), In(n,n);
// MatrixXd Bn(n,n), Bnn(n,n);
// VectorXd x0(n), x(n), b(n);
//
// for (int i=0;i<n;i++)
// {
//   x0(i)=i;
// }
// b=x0;
//
// Bnn=MatrixXd::Random(n,n);
// Bn=(Bnn+MatrixXd::Constant(n,n,1.))/2.;
// C=Bn.sparseView();
// In.setIdentity();
// An=In+1./pow(n,2)*C.transpose()*C;
//
//
// GradPasOptimal(An,b,x0,0.000001,1000,x);
// cout << "Gradient pas optimal donne r="<< endl;
// cout << b-An*x << endl;
//
// ResMin(An,b,x0,0.00001,1000,x);
// cout<<"Résidu minimium donne r="<<endl;
// double res = b-An*x;
// cout<<"r="<<res.norm()<<endl;
//
//
// cout<<"ResMin_cond_gauche"<<endl;
// ResMin_cond_gauche(An,b,x0,0.00001,1000,x);
// cout<<"Résidu minimium preconditionné à gauche donne r="<<endl;
// cout<< (An*x-b).norm() <<endl;
//
// cout<<"ResMin_cond_droite"<<endl;
// ResMin_cond_droite(An,b,x0,0.00001,1000,x);
// cout<<"Résidu minimium preconditionné à droite donne r="<<endl;
// cout<< (An*x-b).norm() <<endl;
//
// cout<<"Gmres"<<endl;
// GMRes(An,b,x0,0.000001,1000,x,m);
// cout<<"GMRes donne r="<<endl;
// cout<<An*x-b<<endl;


//---------------------matrice 25000*25000

// //    Definition    //
//
// SparseMatrix<double> An;
// MatrixXd bn;
// VectorXd x0, x, residu;
// double t1,t2;
//
//
// //    Initialisation    //
//
// An=Lecture_Matrice_A("/home/valentin/DM/DM-main/DM_Turpault/smt/smt.mtx");
// bn=Lecture_Matrice_b("/home/valentin/DM/DM-main/DM_Turpault/smt/smt_b.mtx");
// x0.resize(bn.size());
// x.resize(bn.size());
// residu.resize(bn.size());
//
// for (int i=0; i<bn.size();i++)
//   {
//     x0(i)=i;
//   }
//
//
// //    Calcul    //
//
// cout <<"Calcul ..."<<endl;
// t1 = clock();
// //GradPasOptimal(An,bn,x0,0.00001,1000,x);
// //ResMin(An,bn,x0,0.00001,1000,x);
// ResMin_cond_gauche(An,bn,x0,0.00001,1000,x);
// //GMRes(An,bn,x0,0.00001,1000,x,m);
// t2=clock();
//
// residu = An*x-bn;
// cout << endl << "Norme du résidu = "<< endl;
// cout << residu.norm() << endl;
// //cout << "GMRes donne x="<<endl;
// //cout << x(25709) << endl;
// cout << "Temps d'execution : "<< (t2 - t1) / CLOCKS_PER_SEC <<" en seconde" <<endl;
