#include <Implicit.h>
#include <iostream>
#include <math.h>
#include <armadillo>
#include <tridiag.h>
using namespace arma;
using namespace std;

Implicit::Implicit()
{
}

void Implicit::Implicit_Scheme (vec &V, double alpha, int n, int m, double dx)
{
    mat A(m,m);
    A.zeros();
    A.diag().fill(1+2*alpha);
    A.diag(1).fill(-alpha);
    A.diag(-1).fill(-alpha);
   // A.print();
   // A=inv(A);

  /*  vec V(m);
    for (int i=0; i<m; i++)
    V(i)=U(0,i);    //initial condition, t=0;
   // V.print();*/

    int j=1;
    while (j<n){
        //V=inv(A)*V;
        /*V(0)=1.0;    //boundary condition for x=0;
        V(1)+=alpha/2.0;
        V(m-1)=0.0;  //boundary condition for x=L=1;
        tridiag solver;
        solver.tridiag_solver(A, V, m);*/

          V(0)=1.0;
          vec U(m-2);
          for (int k=0; k<m-2; k++)
                   {U(k)=V(k+1);}

                   U(0)+=alpha;
                   tridiag solver;
                   solver.tridiag_solver(A, U, m-2);

                   for (int k=0; k<m-2; k++)
                   {V(k+1)=U(k);}

           V(m-1)=0.0;


        j++;}

  //  V.print("Implicit=");

    ofstream myfile;
    myfile.open ("Implicit.txt");
    for (int i=0; i<m; i++)
       myfile <<i*dx<<" "<<V(i)<<endl;
       myfile.close();

   return;
}

