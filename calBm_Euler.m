function [cEm]=calBm_Euler(M)
    fm=zeros(1,2*M);
    cEm=zeros(1,2*M);
    f0=0.5;
    fm(end)=2^(-M);
    for m=1:M
       fm(m)=1; 
       if m<M
           fm(2*M-m)=fm(2*M-m+1)+2^(-M)*factorial(M)/(factorial(m)*factorial(M-m));
       end
    end
    for m=1:2*M
        cEm(m)=(-1)^m*fm(m);
    end
    cE0=(-1)^0*f0;
    cEm=[cE0,cEm];
end