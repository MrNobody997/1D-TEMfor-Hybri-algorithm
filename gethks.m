%Obtain the independent variables of time-frequency conversion and Hanker transformation
%Author: You Xiran    Time: October 2023
function [sg,hnk,sw,sc,se,st,sr]=gethks(gsflt,sinflt,cosflt,cEm,cTm,deltsin,deltcos,hankfit,a0,Gup,deltx,t,r)
sw=zeros(length(t),length(sinflt));  
sc=zeros(length(t),length(cosflt));  
sg=zeros(length(t),length(gsflt));   
se=zeros(length(t),length(cEm)); 
st=zeros(length(t),length(cTm)); 
sr=zeros(length(t),21);  
%Euler
Me=(length(cEm)-1)/2;
%Talbot
Deltam=zeros(1,length(cTm)-1);
Mt=length(cTm);
Delta0=2*Mt/5;
for m=1:Mt-1
    Deltam(m)=2*m*pi/5*(cot(m*pi/Mt)+1i);
end
Deltam=[Delta0,Deltam];
%
for ii=1:length(t)
    %sine independent variable
    for jj=1:length(sinflt)
        sw(ii,jj)=exp((jj-60.d0)*deltsin)/t(ii); % sine 160
        %             sw(ii,jj)=exp((jj-57.d0)*deltsin)/t(ii); % sine 123
        %             sw(ii,jj)=exp((jj-150.d0)*deltsin)/t(ii); % sine 250
    end
    %cosine independent variable
    for jj=1:length(cosflt)
        sc(ii,jj)=exp((jj-150.d0)*deltcos)/t(ii); % cosine 250
    end
    %G-S independent variable
    for jj=1:length(gsflt)
        sg(ii,jj)=log(2.0)*jj/t(ii); % ln2/t*m ,m=1,2,...,n
    end
    %Euler independent variable
    Bm=zeros(1,length(cEm)); %2M+1
    for jj=1:length(cEm) % m=[0,2M]
        Bm(jj)=Me*log(10)/3+1i*pi*(jj-1);
        se(ii,jj)=Bm(jj)/t(ii);
    end
    %Talbot independent variable
    for jj=1:length(cTm)
        st(ii,jj)=Deltam(jj)/t(ii);
    end
    %Guptasarma independent variable
    for jj=1:21
        sr(ii,jj)=10^(Gup(jj,1)-log10(t(ii)));
    end
end
hnk=10.^(a0+((1:length(hankfit))-1.d0)*deltx)/r; 
end