%Calculate the coefficients of the Talbot algorithm
%Author: You Xiran    Time: October 2023
function [cTm]=calDeltam_Talbot(M)
    Deltam=zeros(1,M-1);
    cTm=zeros(1,M-1);
    Delta0=2*M/5;
    cT0=0.5*exp(Delta0);
    for m=1:M-1
        Deltam(m)=2*m*pi/5*(cot(m*pi/M)+1i);
        cTm(m)=(1+1i*m*pi/M*(1+cot(m*pi/M)^2)-1i*cot(m*pi/M))*exp(Deltam(m));
    end
    cTm=[cT0,cTm];
end