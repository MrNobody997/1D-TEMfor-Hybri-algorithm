%Calculate kernel function
function [U1,u]=calrte(t,hnk,sw,miu,p,h,am,tao,c)
nl=length(p);
thickness=h;
U1=zeros(nl,length(t),length(sw(1,:)),length(hnk));
u=zeros(nl,length(t),length(sw(1,:)),length(hnk));
resip=zeros(nl,length(t),length(sw(1,:)));
for it=1:length(t)
    for ks=1:length(sw(1,:))
        %cole-cole:
        for mm=1:nl
            resip(mm,it,ks)=p(mm)*(1-am(mm)*(1-1/(1+(sw(it,ks)*tao(mm))^c(mm))))/(1-am(mm));
        end
        for kh=1:length(hnk)
            u(nl,it,ks,kh)=sqrt(hnk(kh)^2 + sw(it,ks)*miu/resip(nl,it,ks));
            U1(nl,it,ks,kh)=u(nl,it,ks,kh);
            for ii=nl-1:-1:1
                u(ii,it,ks,kh)=sqrt(hnk(kh)^2 + sw(it,ks)*miu/resip(ii,it,ks));
                uii=u(ii,it,ks,kh);
                Uii1=U1(ii+1,it,ks,kh);
                tanhah=tanh(uii*thickness(ii));
                U1(ii,it,ks,kh)=uii*( (Uii1+uii*tanhah)/(uii+Uii1*tanhah) );
            end
        end
    end
end
end