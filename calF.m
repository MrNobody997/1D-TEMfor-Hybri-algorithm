%Calculate the field value in frequency domain
%Author: You Xiran    Time: October 2023
function [fd]=calF(t,h0,r,dl,sintheta,I0,nturns,srx,sinflt,hankfit,hnk,U1,miu,sinv,sign1)
nt=length(t);
fd=zeros(1,nt);  
for it=1:nt
    vt=0;
    for ks=1:length(sinflt)
        vs=0;
        for kh=1:length(hankfit)
            mid=(hnk(kh)-U1(1,it,ks,kh))/(hnk(kh)+U1(1,it,ks,kh)); 
            mid=(1+mid)*exp(-h0*hnk(kh))*hnk(kh); 
            vs=vs+mid*hankfit(kh)/r; 
        end
    
        if (sinv==1)
            vt=vt+imag(vs)*sinflt(ks); %sine
        elseif (sinv==0)
            vt=vt+vs*sinflt(ks); %G-S
        elseif (sinv==2)
            vt=vt+real(vs)*sinflt(ks); %Euler
        elseif (sinv==3)
            vt=vt+real(vs*sinflt(ks)); %Talbot
        elseif (sinv==4)
            vt=vt+real(vs*sinflt(ks)); %cosine
        elseif (sinv==5)
            vt=vt+imag(vs)*sinflt(ks,2)*10^(sinflt(ks,1)-log10(t(it))); %Guptasarma
        end
    end
    if (sinv==1)
        fd(it)=0.25/pi*I0*sintheta*(-2.0)/pi*sqrt(pi/2.0)*vt/t(it); 
    elseif (sinv==0)
        fd(it)=0.25/pi*I0*sintheta*log(2.0)*vt/t(it);
    elseif (sinv==2)
        M=(length(sinflt)-1)/2;
        fd(it)=0.25/pi*I0*sintheta*10^(M/3)*vt/t(it);
    elseif (sinv==3)
        fd(it)=0.25/pi*I0*sintheta*0.4*vt/t(it);
    elseif (sinv==4)
        fd(it)=0.25/pi*I0*sintheta*(-2.0)/pi*sqrt(pi/2.0)*vt/t(it);
    elseif (sinv==5)
        fd(it)=0.25/pi*I0*sintheta*vt;
    end
    if (sign1==1)
        fd(it)=fd(it); % dHz/dt
    elseif (sign1==2)
        fd(it)=miu*fd(it); % dBz/dt
    else
        fd(it)=nturns*srx*miu*fd(it); % Vbz
    end
end
end