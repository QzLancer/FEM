%求解非线性电路
%《Nonlinear Axisymmetric Magnetostatic Analysis forElectromagnetic Device Using TLM-Based Finite-Element Method》
%参考nrtlm.m实现分步绘图过程
%二分法，牛顿法和割线法
%By QzLancer
%2019/3/6
close all;
clear all;

I0 = 1; R1 = 10; Y0 = 1;
Vi = 0;
Vr = 0;
V = zeros(20,1);
I2 = 0:0.0001:2;
Ua = 10-10*I2;
Uc = I2.^2;
f = @(I,Vr) I^2*Y0+I-2*Y0*Vr;
plot(I2,Ua,'k',I2,Uc,'k');
axis([0,2,0,8]);
% ylim([0,1]);
hold on;

e = 0.00001;
I21 = 0; Vc = 0;
for i = 1:length(V)
    %入射过程
    Va = (I0+2*Y0*Vi)/(1/R1+Y0);
    Ia = I0 - Va/R1;
    title(['Step:', num2str(i)]);
    plot(Ia,Va,'Color','k','Marker','.','MarkerSize',10);
    line([I21,Ia],[Vc,Va],'Color','k');
    Vr = Va - Vi;
    V(i) = Va;
    %反射过程，I20为左边界，I22为右边界，I21为中点
%     %二分法
%     I20 = 0; I22 = 2;
%     while abs(I22-I20) >= e
%         I21 = (I20+I22)/2;
%         f = (I21)^2*Y0+I21-2*Y0*Vr;
%         if f > 0
%             I22 = I21;
%         else
%             I20 = I21;
%         end
%     end

%     %牛顿迭代法
%     I20 = 0;
% %     I21 = I20-(I20^2*Y0+I20-2*Y0*Vr)/(2*I20*Y0+1);
%     I21 = I20-f(I20,Vr)/(2*I20*Y0+1);
%     while abs(I21-I20) >= e
%         I20 = I21;
%         I21 = I20-(I20^2*Y0+I21-2*Y0*Vr)/(2*I20*Y0+1);
%     end
    
    %割线法
    I20 = 0; I21 = 1;
    I22 = I21 - f(I21,Vr)/(f(I21,Vr)-f(I20,Vr))*(I21-I20);
    while abs(I21-I22) >= e
        I20 = I21;
        I21 = I22;
        I22 = I21 - f(I21,Vr)/(f(I21,Vr)-f(I20,Vr))*(I21-I20);
    end
        
    Vc = (I21)^2;
    Vi = Vc - Vr;
    plot(I21,Vc,'Color','k','Marker','.','MarkerSize',10);
    line([Ia,I21], [Va,Vc],'Color','k');
    hold on;
end
figure(2)
plot([1:length(V)],V,'s');
