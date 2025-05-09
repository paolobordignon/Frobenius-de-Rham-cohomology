import sys
from sage.all import *

p = 5
Z5 = Zp(p, prec = 10, type = 'capped-abs', print_mode = 'series')
R = PowerSeriesRing(Z5,'x')
x = R.gen()
S = PowerSeriesRing(R,'z')
z = S.gen()

E = EllipticCurve([26,1])
f = x**3+E.a4()*x+E.a6()
frob_x = x**p
#frob_y = y^p*(1+(f(x^p)-f(x)^p)*1/f(x)^p)^(1/2)
frob_z = z**p*(1+(f(x**p)-f(x)**p)*z^(2*p))**(-1/2)

frob_omega = p*x**(p-1)*frob_z
frob_eta = p*x**(2*p-1)*frob_z

V = LaurentSeriesRing(S)

from sage.all import *

def reduction_coeff(pow_ser):
    pow_ser_2 = pow_ser
    for i in range(len(pow_ser_2.coefficients())):
        pow_ser += -z**pow_ser_2.exponents()[i]*pow_ser_2.coefficients()[i]+z**pow_ser_2.exponents()[i]*reduction_deg_2g(f,pow_ser_2.coefficients()[i])
    return pow_ser

def reduction_deg_2g(f,pol):
    #f=x^3+ax+b
    if pol(z=1).degree()<3:
        return pol
    else: 
        n = pol(z=1).degree()
        pol += - pol(z=1).coefficients()[len(pol(z=1).coefficients())-1]*(f**(floor(n/3))-1/z**(2*floor(n/3)))*x**(n%3)
        return reduction_deg_2g(f,pol)

def differential_zn(f,i,n):
    """
    Compute the exact differential d(x^iz^n) with i=0,1,2 and n integer
    return p(x,z) where p(x,z)dx=d(x^i z^n) and p(x,z)=\sum p_i(x)z^i with degree(p_i(x))<3 
    """
    if n == 0:
        return i*x^(i-1)*z^0
    return (i*x^(i-1)*z^n+reduction_deg_2g(f,-1/2*n*f.derivative(x)*x^i)*z^(n+2))


def matrix_coeff_differentials_pos(n):
    matrix_diff = matrix(QQ,3)
    for i in range(3):
        if differential_zn(f,i,n) == 0:
            pol_coef = 0
        else:
            pol_coef = differential_zn(f,i,n).coefficients()[-1]
        for j in range(3):
            if pol_coef == 0:
                matrix_diff[j,i] = 0
            else:
                matrix_diff[j,i] = pol_coef[j]
    return matrix_diff

def matrix_coeff_differentials_neg(n):
    matrix_diff = matrix(QQ,3)
    for i in range(3):
        if differential_zn(f,i,n) == 0:
            pol_coef = 0
        else:
            if i == 0:
                pol_coef = differential_zn(f,i,n-2).coefficients()[0]
            else:
                pol_coef = differential_zn(f,i,n).coefficients()[0]
        for j in range(3):
            if pol_coef == 0:
                matrix_diff[j,i] = 0
            else:
                matrix_diff[j,i] = pol_coef[j]
    return matrix_diff

def vector_coeff(pol):
    return vector([pol[0],pol[1],pol[2]])

def reduction_z_pos(pow_ser,exact_form=0):
    if pow_ser.degree()<3:
        return (pow_ser,exact_form)
    else:
        k = pow_ser.degree()
        M = matrix_coeff_differentials_pos(k-2)
        V = vector_coeff(pow_ser.coefficients()[-1])
        X = M.solve_right(V)
        exact_form += X.dot_product(vector([x^0*z^(k-2),x^1*z^(k-2),x^2*z^(k-2)]))
        return reduction_z_pos(pow_ser-X.dot_product(vector([differential_zn(f,0,k-2),differential_zn(f,1,k-2),differential_zn(f,2,k-2)])),exact_form) 

def reduction_z_neg(pow_ser,exact_form=0):
    if pow_ser.exponents()[0]>-1:
        return (pow_ser,exact_form)
    else:
        k = pow_ser.exponents()[0]
        M = matrix_coeff_differentials_neg(k)
        V = vector_coeff(pow_ser.coefficients()[0])
        X = M.solve_right(V)
        exact_form += X.dot_product(vector([x^0*z^(k-2),x^1*z^k,x^2*z^k]))
        return reduction_z_neg(pow_ser-X.dot_product(vector([differential_zn(f,0,k-2),differential_zn(f,1,k),differential_zn(f,2,k)])),exact_form)

def matrix_frobenius(p,frob_omega,frob_eta,prec=10, exact_dif = False):

    Z5 = Zp(p, prec, type = 'capped-abs', print_mode = 'series')

    red_frob_eta = reduction_coeff(frob_eta)
    red_frob_omega = reduction_coeff(frob_omega)

    red_frob_eta_pos, exact_pos= reduction_z_pos(red_frob_eta)
    red_frob_eta_fin, exact_form_eta= reduction_z_neg(red_frob_eta_pos,exact_pos)
    red_frob_omega_pos, exact_pos= reduction_z_pos(red_frob_omega)
    red_frob_omega_fin, exact_form_omega= reduction_z_neg(red_frob_omega_pos,exact_pos)

    matrix_frobenius = matrix(Z5,2)
    matrix_frobenius[0,0] = red_frob_omega_fin[1][0]
    matrix_frobenius[0,1] = red_frob_omega_fin[1][1]
    matrix_frobenius[1,0] = red_frob_eta_fin[1][0]
    matrix_frobenius[1,1] = red_frob_eta_fin[1][1]

    if exact_dif:
        return (matrix_frobenius, [exact_form_omega,exact_form_eta])
    else:
        return matrix_frobenius

print(matrix_frobenius(p,frob_omega,frob_eta))