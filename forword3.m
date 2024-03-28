%Calculate the response of a single dipole:
function [fd,resip]=forword3(gsflt,sinflt,cosflt,cEm,cTm,deltsin,deltcos,hankfit,a0,Gup,deltx,times,Mxyz,x,miu,p,h,am,tao,c,I0,nturns,srx,sinv,sign1)
rx=Mxyz(1);
ry=Mxyz(2);
h0=Mxyz(3);
r=sqrt((rx-x)^2+ry^2);
sintheta=ry/r;
[sg,hnk,sw,sc,se,st,sr]=gethks(gsflt,sinflt,cosflt,cEm,cTm,deltsin,deltcos,hankfit,a0,Gup,deltx,times,r);
if sinv==1
    [U1,~]=calrte(times,hnk,1i*sw,miu,p,h,am,tao,c);
    [fd]=calF(times,h0,r,1,sintheta,I0,nturns,srx,sinflt,hankfit,hnk,U1,miu,sinv,sign1);
elseif sinv==0
    [U1,~]=calrte(times,hnk,sg,miu,p,h,am,tao,c);
    [fd]=calF(times,h0,r,1,sintheta,I0,nturns,srx,gsflt,hankfit,hnk,U1,miu,sinv,sign1);
elseif sinv==2
    [U1,~]=calrte(times,hnk,se,miu,p,h,am,tao,c);
    [fd]=calF(times,h0,r,1,sintheta,I0,nturns,srx,cEm,hankfit,hnk,U1,miu,sinv,sign1);
elseif sinv==3
    [U1,~]=calrte(times,hnk,st,miu,p,h,am,tao,c);
    [fd]=calF(times,h0,r,1,sintheta,I0,nturns,srx,cTm,hankfit,hnk,U1,miu,sinv,sign1);
elseif sinv==4
    [U1,~]=calrte(times,hnk,1i*sc,miu,p,h,am,tao,c);
    [fd]=calF(times,h0,r,1,sintheta,I0,nturns,srx,cosflt,hankfit,hnk,U1,miu,sinv,sign1);
elseif sinv==5
    [U1,~]=calrte(times,hnk,1i*sr,miu,p,h,am,tao,c);
    [fd]=calF(times,h0,r,1,sintheta,I0,nturns,srx,Gup,hankfit,hnk,U1,miu,sinv,sign1);
end
end