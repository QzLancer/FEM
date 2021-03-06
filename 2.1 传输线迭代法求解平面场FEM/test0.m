%用传输线迭代法求解简单的有限元模型
%采用的模型是《Numerical techniques in electromagnetics Chapter6.2.4》
%By QzLancer
%2019/3/10
clear all;
%--------------------------------初始化几何参数
X = [0.8;1.4;2.1;1.2];
Y = [1.8;1.4;2.1;2.7];
TriElement = [1,2,4;2,3,4];
TriX = X(TriElement);
TriY = Y(TriElement);
Y0 = 1; %传输线导纳
e = 0.01;
%--------------------------------求解必要的几何参数
p(:,1) = TriX(:,2).*TriY(:,3) - TriY(:,2).*TriX(:,3);
p(:,2) = TriX(:,3).*TriY(:,1) - TriY(:,3).*TriX(:,1);
p(:,3) = TriX(:,1).*TriY(:,2) - TriY(:,1).*TriX(:,2);
q(:,1) = TriY(:,2) - TriY(:,3);
q(:,2) = TriY(:,3) - TriY(:,1);
q(:,3) = TriY(:,1) - TriY(:,2);
r(:,1) = TriX(:,3) - TriX(:,2);
r(:,2) = TriX(:,1) - TriX(:,3);
r(:,3) = TriX(:,2) - TriX(:,1);
Area = (q(:,1).*r(:,2) - q(:,2).*r(:,1))/2;
%---------------------------------单元分析和总体合成，使用传输线迭代法
V0 = zeros(4,1);
V1 = ones(4,1);
Ym = zeros(4,4);
% Vi = zeros(2,6);
% Vr = zeros(2,6);
Vi = zeros(3,3);
Vr = zeros(3,3);
Vc = zeros(3,3);
I = zeros(4,1);
Se = zeros(3,3);
%组装大型系数矩阵，单元内的所有元素都是Y0
for k = 1:2
    for i = 1:3
        for j = 1:3
            if i == j
                Ym(TriElement(k,i),TriElement(k,j)) = Ym(TriElement(k,i),TriElement(k,j))+3*Y0;
            else
                Ym(TriElement(k,i),TriElement(k,j)) = Ym(TriElement(k,i),TriElement(k,j))-Y0;
            end
        end
    end
end
%求解入射过程电压Va，直接法
Va1 = [0;0;10;0];
I = I - Ym*Va1;
Va1([2,4]) = Ym([2,4],[2,4])\I([2,4]); 
Va0 = [100;100;100;100];
%反射过程求解，装配，入射过程求解，循环
while norm(Va1-Va0) >= e %循环判定，范数？
    Va0 = Va1;
    I = zeros(4,1);
    for k = 1:2
        for i = 1:3
            for j = 1:3
                if i == j
                    Vr(i,j) = Va0(TriElement(k,i))-Vi(i,j);
                else
                    Vr(i,j) = Va0(TriElement(k,i))-Va0(TriElement(k,j))-Vi(i,j);
                end
                Se(i,j) = (q(k,i)*q(k,j)+r(k,i)*r(k,j))/(4*Area(k));
                if i == j
                    Vc(i,j) = 2*Vr(i,j)*Y0/(Y0+Se(i,1)+Se(i,2)+Se(i,3));
                else
                    Vc(i,j) = 2*Vr(i,j)*Y0/(Y0-Se(i,j));
                end
                Vi(i,j) = Vc(i,j) - Vr(i,j);
                I(TriElement(k,i)) = I(TriElement(k,i))+2*Y0*Vi(i,j);
            end
        end
    end
    Va1 = [0;0;10;0];
    I = I - Ym*Va1;
    Va1([2,4]) = Ym([2,4],[2,4])\I([2,4]); 
end
