%1D Ground sources TEM forward modeling
%Author: Yxr    Time: October 2023
clc
clear all
format short e
%Load filtering coefficients:
tic
%Setting forward modeling parameters:
I0=10; % Current
L=1000; % Length of the source wire
Mxyz=[0,500,0]; %Location of measuring points(x,y,z)
srx=1; % Effective area of receiving coil
nturns=1; % Turn Ratio 
Lx=[-L/2,L/2]; % The x-coordinates at both ends of the wire (along the x-axis direction)
Ta=-7; Tb=-0; % Time range(log10)
trace=100; % Time traces
times=logspace(Ta,Tb,trace);
%%
%Setting model parameters:
sinv=1; %Time-frequency conversion sign:0,G-S;1,sine;2,cosine;3,Euler;4,Talbot;5,Guptasarma.
sign1=2; % Output sign = 1.dHz|2.dBz|3.Vbz
p=[100 500 100]; % Resistivity
h=[100 10]; % Layer thickness
%Induced polarization parameters():
IP=0; % The sign of enabling IP effect, 1, Yes; 0, No.
am=[ 0.0 0.8 0.0 0]; % Polarizability
tao=[0.1 0.1 0.1 0]; % Time constant
c=[0.25 0.25 0.25 0]; % Frequency-dependant coefficient
%%
for M=12:1:34  %Different coefficients
[gsflt,sinflt,cosflt,cEm,cTm,deltsin,deltcos,hankfit,a0,Gup,deltx,miu]=loadsinhank(M);
%1D Forward modeling:
if IP==0
    am=zeros(1,length(p));
    tao=am; c=am;
end
f=@(x)forword3(gsflt,sinflt,cosflt,cEm,cTm,deltsin,deltcos,hankfit,a0,Gup,deltx,times,Mxyz,x,miu,p,h,am,tao,c,I0,nturns,srx,sinv,sign1);
[V]=IntGauss_Legendre(f,Lx(1),Lx(2)); %¶Ôdl»ý·Ö
V=abs(V);
toc

%Plotting:
figure(1)
loglog(times,abs(V),'-','LineWidth',1)
ylabel('dHz/dt');xlabel('Time');
hold on
legend
%%
%Half-space analytical formula:
for n=1:trace
    t=times(n);
    Hzt_jx=Return_sum_jx_field_t(Lx(1),Lx(2),I0,p,t,miu,Mxyz(1),Mxyz(2));
    dBz_jx(n)=abs(Hzt_jx)*miu;
end
% figure(1)
% loglog(times,dBz_jx,'k:','LineWidth',2.0);
% ylabel('dBz/dt');xlabel('Time');
% hold on
% legend
end