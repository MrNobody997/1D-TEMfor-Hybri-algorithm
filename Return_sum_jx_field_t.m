function[Hzt_jx_t]=Return_sum_jx_field_t(Aset,Bset,I,res,t,u0,rxx,ry)
              f_jx=@(x)Return_jx_Hz_t(I,res,t,u0,rxx,ry,x); %
              [Hzt1_jx_t]=IntGauss_Legendre(f_jx,Aset,Bset); %高斯-勒让德积分 电偶极子求积分
              Hzt_jx_t=Hzt1_jx_t;