%求解非线性电路
%《Nonlinear Axisymmetric Magnetostatic Analysis forElectromagnetic Device Using TLM-Based Finite-Element Method》
%By QzLancer
%2019/3/6

clear all
I0 = 1; R1 = 10; Y0 = 1;
Vi = 0;
Vr = 0;
V = zeros(50,1);
x = 0:0.01:1;
plot(x,10-10*x,x,10*x);
% ylim([0,1]);
hold on;
%二分法
e = 0.00001;
for i = 1:50
    %入射过程
    Va = (I0+2*Y0*Vi)/(1/R1+Y0);
    Ia = I0 - Va/R1;
    plot(Ia,Va,'o');
    Vr = Va - Vi;
    V(i) = Va;
    %反射过程，I20为左边界，I22为右边界，I21为中点
    I20 = 0; I22 = 2;
    while abs(I22-I20) >= e
        I21 = (I20+I22)/2;
        f = 10*I21*Y0+I21-2*Y0*Vr;
        if f > 0
            I22 = I21;
        else
            I20 = I21;
        end
    end
    Vc = 10*I21;
    Vi = Vc - Vr;
    plot(I21,Vc,'*');
    hold on;
end