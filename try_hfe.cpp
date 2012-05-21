
#include <cstdio>
#include <functional>
#include <iostream>
#include <vector>
#include <ginac/ginac.h>
using namespace std;
using namespace GiNaC;

// modular reduction
struct fmodnum : public map_function {
public:
    numeric P;
    fmodnum(const numeric& P) : P(P) { };
    ex operator()(const ex &e) {
        if(is_a<mul>(e)) {
            ex m=1;
            for(size_t i=0; i!=e.nops(); i++) {
                const ex& ei = e.op(i);
                if(is_a<numeric>(ei))
                    m *= mod(ex_to<numeric>(ei), P);
                else
                    m *= ei.map(*this);
            }
            return m;
        }
        if(is_a<add>(e)) {
            ex m=0;
            for(size_t i=0; i!=e.nops(); i++) {
                const ex& ei = e.op(i);
                if(is_a<numeric>(ei))
                    m += mod(ex_to<numeric>(ei), P);
                else
                    m += ei.map(*this);
            }
            return m;
        }
        else {
            return e.map(*this);
        }
    }
};

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
    
    
    
    
    
    
    ex fp = f
            .subs(A==aa)
            .expand();
    
    // poly reduction
    fp = rem(fp, irr, X);
    
    // doing frobenius
    for(int i=0; i<N; i++) {
        ex frob = pow(va[i],P);
        while(fp.has(frob, has_options::algebraic)) {
            fp = fp.subs(frob==va[i], subs_options::algebraic);
        }
    }
    
    // doing modular
    {
        fmodnum map1(P);
        fp = fp.map(map1).collect(X);
    }
    
    {
        ex fpt = fp/*.expand()*/;
        for(int i=0; i<N; i++) {
            cout << "b_" << i << " = " << fpt.coeff(X,i) << endl;
        }
    }
    
    cout << fp;
    
    return 0;
}