
#include <cstdio>
#include <functional>
#include <iostream>
#include <vector>
#include <ginac/ginac.h>
using namespace std;
using namespace GiNaC;

#include "algebra.h"

int main()
{
    int N=3;
    int P=2;

    
    symbol A("A"), X("X");
    vector<symbol> va; va.reserve(N);
    ex irr = pow(X,3) + pow(X,2) + 1;
    ex f = A + pow(A,3) + pow(A,5);
    ex aa;
    
    for(int i=0; i<N; i++) {
        char buf[8];
        sprintf(buf,"%d",i);
        va.push_back(symbol("a_"+string(buf)));
        
        aa += va[i]*pow(X,i);
    }

    // fp is multivariate polynomial
    // fp = X^(N-1)*A_(N_1) + ... + A_0

    ex fp = rem(f.subs(A==aa).expand(), irr, X);
    
    // doing reduction
    for(int i=0; i<N; i++) {
        fp = frobenius(fp, va[i], P);
    }
    fp = modular(fp, P);

    

    vector<ex> vb; vb.reserve(N);
    
    cout << fp << endl;

    {
        ex fpt = fp.collect(X);
        for(int i=0; i<N; i++) {
            vb.push_back(fpt.coeff(X,i));
            cout << "b_" << i << " = " << vb[i] << endl;
        }
        cout << fpt << endl;
    }
    
    
    return 0;
}
