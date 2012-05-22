
#include "algebra.h"

namespace GiNaC {

ex modular(const ex& f, const numeric& n) {
    struct fmodular : public map_function {
    public:
        numeric n;
        fmodular(const numeric& n) : n(n) { };
        ex operator()(const ex &e) {
            if(is_a<mul>(e)) {
                ex m=1;
                for(size_t i=0; i!=e.nops(); i++) {
                    const ex& ei = e.op(i);
                    if(is_a<numeric>(ei))
                        m *= mod(ex_to<numeric>(ei), n);
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
                        m += mod(ex_to<numeric>(ei), n);
                    else
                        m += ei.map(*this);
                }
                return m;
            }
            else {
                return e.map(*this);
            }
        }
    } m(n);
    return f.map(m);
}

ex frobenius(const ex& f, const ex& a, const numeric& p) {
    ex fn = f;
    ex frob = pow(a,p);
    while(fn.has(frob, has_options::algebraic)) {
        fn = fn.subs(frob==a, subs_options::algebraic);
    }
    return fn;
}

}; // namespace GiNaC

