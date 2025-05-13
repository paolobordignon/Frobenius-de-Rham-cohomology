

def horizontal_basis(N,N_pow_ser = False):
    R2.<t> = LaurentSeriesRing(QQ)
    R1.<a,b,c,d> = PolynomialRing(QQ)
    R.<t> = LaurentSeriesRing(R1,'t',default_prec = 20)
    N2 = matrix([[R(N[0,0]),R(N[0,1])],[R(N[1,0]),R(N[1,1])]])
    U = identity_matrix(R,2)
    U1 = matrix([[a,b],[c,d]])
    
    for i in range(1,R.default_prec()):
        U += U1*t^i
        N1 = ((N2*U+U.derivative(t))/t^(i-1))(t=0)
        J = R1.ideal(N1.list())
        U = U(a=J.reduce(a),b=J.reduce(b),c=J.reduce(c),d=J.reduce(d))

    U2 = matrix([[R2(U[0,0]),R2(U[0,1])],[R2(U[1,0]),R2(U[1,1])]])

    if N_pow_ser:
        return (U2,N2)
    else:
        return U2