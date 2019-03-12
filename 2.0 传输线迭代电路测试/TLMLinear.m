%求解简单电路
%《Nonlinear Axisymmetric Magnetostatic Analysis forElectromagnetic Device Using TLM-Based Finite-Element Method》
%by QzLancer
%2019/3/6
 
clear all;
I0 = 1; R1 = 10; R2 = 10; Y0 = 1;
Vi = 0; Vr = 0;
V = zeros(50, 1);
for i = 1:50
    %入射过程
    Va = (I0+2*Vi*Y0)/(Y0+1/R1);
    Vr = Va - Vi;
    V(i) = Va;
    %反射过程
    Vc = 2*Vr*Y0/(Y0+1/R2);
    Vi = Vc - Vr;
end