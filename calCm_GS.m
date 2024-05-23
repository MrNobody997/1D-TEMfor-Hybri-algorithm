%Calculate the coefficients of the G-S algorithm
%Author: You Xiran    Time: October 2023
function [Cmg]=calCm_GS(M)
Cmg=zeros(1,M);  
for m=1:M
    Csum=0.0;
    for k=floor((m+1.0)/2.0):min(m,M/2.0)
        A=k^(M/2.0)*factorial(2.0*k);
        B=factorial(M/2.0-k)*factorial(k)*factorial(k-1)*factorial(m-k)*factorial(2.0*k-m);
        Csum=Csum+A/B;
    end
    Cmg(m)=(-1)^(M/2.0+m)*Csum;
end
Cmg=Cmg';
end