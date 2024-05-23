%Gauss-Legendre method for integral calculation
%Author: You Xiran    Time: October 2023
function[I]=IntGauss_Legendre(f,a,b)
%Gauss Legendre integration(n=5):
t_start=(b-a)/2.;
t_end=(a+b)/2.;
I=t_start*(0.2369268851*f(t_end+t_start*(-0.9061798459))+0.4786286705*...
    f(t_end+t_start*(-0.5384693101))+ 0.5688888889*f(t_end+t_start*0.0)+0.4786286705*...
    f(t_end+t_start*0.5384693101)+0.2369268851*f(t_end+t_start*0.9061798459));
end