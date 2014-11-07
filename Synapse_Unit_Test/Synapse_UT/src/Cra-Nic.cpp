#include <iostream>
#include <math.h>
#include <Cra-Nic.h>
#include <tridiag.h>
#include <armadillo>
using namespace std;
using namespace arma;

Crank_Nicolson::Crank_Nicolson(){}

void Crank_Nicolson::Crank_Nicolson_Scheme(vec &V, double alpha, int n, int m, double dx)
{ /*vec V(m);
    for (int i=0; i<m; i++)
    V(i)=U(0,i);    //initial condition, t=0;*/

  mat I(m,m);
    I.eye();

  mat B(m,m);
    B.diag(-1).fill(-1);
    B.diag(1).fill(-1);
    B.diag().fill(2);

    mat A(m,m);
    A=2*I+alpha*B;

    vec V_new(m);
   // vec V1=V;

    int j=1;
    while (j<n)

   /* { V=inv(2*I+alpha*B)*(2*I-alpha*B)*V;
      V(0)=1.0; //boundary condition for x=0;
      j++;}*/

    {  // V=(2*I-alpha*B)*V;

        for (int i=0; i<m; i++)
        {if (i==0) //V_new(i)=(2.0 - 2.0*alpha)*V(i) + alpha*V(i+1);
                   V_new(i)=1.0;
            else{ if (i<(m-1)) V_new(i)=alpha*V(i-1) + (2.0 - 2.0*alpha)*V(i) + alpha*V(i+1);
                 else //V_new(i)=alpha*V(i-1) + (2.0 - 2.0*alpha)*V(i);
                      V_new(i)=0.0;
            }
            } V=V_new;

    /* tridiag solve;
     solve.tridiag_solver(A, V, m);
        V(0)=1.0;   //boundary condition for x=0;
        V(m-1)=0.0; //boundary condition for x=L=1;*/

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

   // V.print("Crank-Nicolson=");

    ofstream myfile;
    myfile.open ("Cra-Nic.txt");
    for (int i=0; i<m; i++)
       myfile <<i*dx<<" "<<V(i)<<endl;
       myfile.close();

return;}


