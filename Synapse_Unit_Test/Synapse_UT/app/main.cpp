#include <iostream>
#include <armadillo>
#include <Explicit.h>
#include <Implicit.h>
#include <Closed_form.h>
#include <Cra-Nic.h>
using namespace arma;
using namespace std;

int main()
{
    int n, m;
    double a=0.0, b=1.0, t_begin=0.0, t_final;
    cout<<"Insert the number of steps in time"<<endl<<"n=";
    cin>>n;
    cout<<"Insert the end of the interval in time"<<endl<<"t_final=";
    cin>>t_final;
    cout<<"Insert the number of steps for position"<<endl<<"m=";
    cin>>m;

    //define steps for t and x:
    double dt=(t_final-t_begin)/(n-1);
    double dx=(b-a)/(m-1);
    double alpha=dt/(dx*dx);
    cout<<"alpha="<<alpha<<endl;
    cout<<"dt="<<dt<<endl;
    cout<<"dx="<<dx<<endl;

    //create array U
    mat U(n,m);

    //boundary conditions
    U.col(0).ones();
    U.col(m-1).zeros();
    //initial conditions
    U.row(0).zeros();

    vec V_Ex(m), V_Im(m), V_CN(m);
    for (int i=0; i<m; i++)
    V_Ex(i)=U(0,i);    //initial condition, t=0;
    V_Im.zeros();
    V_CN.zeros();

    Explicit method;
    method.Explicit_Scheme(V_Ex, alpha, n, m, dx);

    Implicit solve;
    solve.Implicit_Scheme(V_Im, alpha, n, m, dx);

    Crank_Nicolson result;
    result.Crank_Nicolson_Scheme(V_CN, alpha, n, m, dx);

    Closed_form test;
    test.Closed_form_solution(n, m, t_final, dx);

    cout<<"The end"<<endl;

    return 0;
}
