%Calculate 1D TEM analytical solutions for frequency domain fields
%Author: You Xiran    Time: October 2023
function[Hzt_jx_t]=Return_jx_Hz_t(I,res,t,u0,rxx,ry,x)
        rox1=rxx-x;
        roy1=ry;
        r=sqrt(rox1^2+roy1^2);
        sinf=roy1/r;
        u=r*sqrt(u0/2./t/res(1)); 
        fa=erf(u/sqrt(2.)); 
                Hzt_jx_t=fa-sqrt(2/pi)*u*(1+u^2/3)*exp(-u^2/2);
        Hzt_jx_t=Hzt_jx_t*3*I*sinf*res(1)/2/pi/r/r/r/r/u0;  
end