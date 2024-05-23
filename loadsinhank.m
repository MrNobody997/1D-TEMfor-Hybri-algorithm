%Load filtering coefficients
%Author: You Xiran    Time: October 2023
function [gsflt,sinflt,cosflt,cEm,cTm,deltsin,deltcos,hankfit,a0,Gup,deltx,miu]=loadsinhank(M)
  sinflt=load('sinf.txt'); % 160 points sine
%   sinflt=load('sin_xs_123.txt'); % 123 points sine
%   sinflt=load('sin_250.txt'); % 250 points sine
  cosflt=load('cos_250.txt'); %250 points cosine
  hankfit=load('hkf.txt'); % 140 points J1 Hankel filtering coefficient
  a0=-7.91001919000d+00; % Initial value of Hankel
  Gup=load('Guptasarma.txt'); % Guptasarma filtering coefficients
  miu=pi*4d-7; 
  deltx=8.79671439570d-02; % Hanke sampling interval
  deltsin=log(10.0d0)/10.0d0; % 160 points sine sampling interval(ln(10)/10)
%   deltsin=log(10.0d0)/20.0d0; % 250 points sine sampling interval(ln(10)/20)
  deltcos=log(10.0d0)/20.0d0; % 250 points sine sampling interval(ln(10)/20)
%Calculate the coefficients
  if mod(M,2)==0
      gsflt=calCm_GS(M); %G-S
  else
      gsflt=[];
  end
  cEm=calBm_Euler(M); %Euler
  cTm=calDeltam_Talbot(M); %Talbot
end