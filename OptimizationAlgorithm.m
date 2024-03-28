%1D Forward modeling of ground sources TEM based on adaptive hybrid algorithm
%Author: Yxr    Time: October 2023
clc
clear all
format short e
%Load filtering coefficients:
hankfit=load('hkf.txt'); % 140 points J1 Hankel filtering coefficient
a0=-7.91001919000d+00; % Initial value of Hankel
deltx=8.79671439570d-02; % Hanke sampling interval
miu=pi*4d-7; % μ0=π*4*10^-7
gsflt=calCm_GS(14); % load 14 points G-S coefficients
cTm=calDeltam_Talbot(20); % load 20 points Talbot coefficients
%%
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
sign1=2; % Output sign = 1.dHz|2.dBz|3.Vbz
p=[100 500 100]; % Resistivity
h=[100 10]; % Layer thickness
%Induced polarization parameters():
IP=0; % The sign of enabling IP effect, 1, Yes; 0, No.
am=[ 0.0 0.8 0.0 0]; % Polarizability
tao=[0.1 0.1 0.1 0]; % Time constant
c=[0.25 0.25 0.25 0]; % Frequency-dependant coefficient
%%
%1D Forward modeling:
tic 
if IP==0
    am=zeros(1,length(p));
    tao=am; c=am;
end
%First calculation:
f=@(x)forword3(gsflt,[],[],[],[],[],[],hankfit,a0,zeros(21,2),deltx,times,Mxyz,x,miu,p,h,am,tao,c,I0,nturns,srx,0,sign1);
[Vg]=IntGauss_Legendre(f,Lx(1),Lx(2));
V=abs(Vg);
%twice calculation:
iii=trace;
f=@(x)forword3([],[],[],[],cTm,[],[],hankfit,a0,zeros(21,2),deltx,times(iii),Mxyz,x,miu,p,h,am,tao,c,I0,nturns,srx,3,sign1);
[Vt]=IntGauss_Legendre(f,Lx(1),Lx(2));
Vt=abs(Vt);
while( abs(Vt-V(iii))>5e-6  )
    V(iii)=Vt;
    iii=iii-1;
    if iii<1
        break;
    end
    f=@(x)forword3([],[],[],[],cTm,[],[],hankfit,a0,zeros(21,2),deltx,times(iii),Mxyz,x,miu,p,h,am,tao,c,I0,nturns,srx,3,sign1);
    [Vt]=IntGauss_Legendre(f,Lx(1),Lx(2)); 
    Vt=abs(Vt);
end
toc
%Plotting:
figure(1)
loglog(times,abs(V),'-','LineWidth',1)
ylabel('dBz/dt');xlabel('Time');
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
%%
%Abs E
Err=abs(dBz_jx-V)./dBz_jx;
figure(2)
loglog(times,Err);
ylabel('相对误差百分比');xlabel('Time');
hold on
legend
%MSE
MSE_GS=mse(dBz_jx,V);