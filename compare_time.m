%Compare calculation time
%Author: You Xiran    Time: October 2023
clc
clear all
format short e
I0=10;     % Current
L=1000;     % Length of the source wire
Mxyz=[0 500 0];   %Location of measuring points(x,y,z)
srx=1;     % Effective area of receiving coil
nturns=1;    % Turn Ratio
Lx=[-L/2,L/2];  % The x-coordinates at both ends of the wire (along the x-axis direction)
Ta=-7; Tb=-0;    % Time range(log10)
trace=100;       % Time traces
times=logspace(Ta,Tb,trace);
%Set model parameters:
sinv=3;     %Time-frequency conversion sign:0,G-S;1,sine;2,cosine;3,Euler;4,Talbot;5,Guptasarma.
sign1=2;    % Output sign = 1.dHz;2.dBz;3.Vbz
p=[100 500 100];
h=[100 200 100];
am=[ 0.0 0.0 0.0 0];
tao=[0.1 0.1 0.1 0];
c=[0.25 0.25 0.25 0];
%%
% for m=1:3
% T=zeros(1,48);
% M=1;
T=zeros(1,34);
for tt=1:10
    tt
    for M=10:2:16
        %Load filtering coefficients:
        [gsflt,sinflt,cosflt,cEm,cTm,deltsin,deltcos,hankfit,a0,Gup,deltx,miu]=loadsinhank(M);
        tic
        %Forward modeling calculation:
        f=@(x)forword3(gsflt,sinflt,cosflt,cEm,cTm,deltsin,deltcos,hankfit,a0,Gup,deltx,times,Mxyz,x,miu,p,h,am,tao,c,I0,nturns,srx,sinv,sign1);
        [V]=IntGauss_Legendre(f,Lx(1),Lx(2)); 
        V=abs(V);
        
        %Plotting:
        % figure(1)
        % loglog(times,abs(V),'-','LineWidth',1)
        % ylabel('dBz/dt');xlabel('Time');
        % hold on
        % legend
        
        %Half-space analytical formula:
%         for n=1:trace
%             t=times(n);
%             Hzt_jx=Return_sum_jx_field_t(Lx(1),Lx(2),I0,p,t,miu,Mxyz(1),Mxyz(2));
%             dBz_jx(n)=abs(Hzt_jx)*miu;
%         end
        tsum=toc
        T(M)=T(M)+tsum;
        % figure(1)
        % loglog(times,dBz_jx,':','LineWidth',2.0);
        % ylabel('dBz/dt');xlabel('Time');
        % hold on
        % legend
    end
end

figure(1)
semilogy(T./10,'k-o','LineWidth',1.0);ylabel('Average time/s');xlabel('M');hold on
tmin=min(min(T(:,2:end)))
tmax=max(max(T(:,2:end)))
%%
%Batch Legend:
legendShow={ };j=1;
for i=4:2:16
    legendShow{1,j}=num2str(i);
    j=j+1;
end
legend(legendShow)