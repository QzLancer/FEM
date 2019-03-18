%�ô����ߵ��������򵥵�����Ԫģ��
%���õ�ģ���ǡ�Numerical techniques in electromagnetics Chapter6.2.4��
%��Ч��������ʱ���淽��������
%V2 = 3.7075, V4 = 4.4386.
%By QzLancer
%2019/3/10
clear all;
%--------------------------------��ʼ�����β���
X = [0.8;1.4;2.1;1.2];
Y = [1.8;1.4;2.1;2.7];
TriElement = [1,2,4;2,3,4];
TriX = X(TriElement);
TriY = Y(TriElement);
Y0 = 1; %�����ߵ���
e = 0.0000001;
%--------------------------------����Ҫ�ļ��β���
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
%---------------------------------��Ԫ����������ϳɣ�ʹ�ô����ߵ�����
V0 = zeros(4,1); 
V1 = ones(4,1);
Ym = zeros(4,4);
Vi = zeros(3,3,2);
Vr = zeros(3,3,2);
Vc = zeros(3,3,2);
I = zeros(4,1);
Se = zeros(3,3,2);
S = zeros(4,4);
F = zeros(4,1);
%��װ����ϵ������
for k = 1:2
    for i = 1:3
        for j = 1:3
            Se(i,j,k) = (q(k,i)*q(k,j)+r(k,i)*r(k,j))/(4*Area(k));
            S(TriElement(k,i),TriElement(k,j)) = S(TriElement(k,i),TriElement(k,j))+Se(i,j,k);
            if i == j
                Ym(TriElement(k,i),TriElement(k,j)) = Ym(TriElement(k,i),TriElement(k,j))+2*Y0;
            else
                if Se(i,j,k) < 0
                    Ym(TriElement(k,i),TriElement(k,j)) = Ym(TriElement(k,i),TriElement(k,j))-Y0;
                else %���ڸ���������
                    Ym(TriElement(k,i),TriElement(k,i)) = Ym(TriElement(k,i),TriElement(k,i))-Y0;
                end
            end
        end
    end
end
%���
Solve = [0;0;10;0];
F = F-S*Solve;
Solve([2,4]) = S([2,4],[2,4])\F([2,4]);
%���������̵�ѹVa��ֱ�ӷ�
Va1 = [0;0;10;0];
I = I - Ym*Va1;
Va1([2,4]) = Ym([2,4],[2,4])\I([2,4]); 
Va0 = [100;100;100;100];
m = 1;
%���������⣬װ�䣬���������⣬ѭ��
while norm(Va1-Va0)/norm(Va1) >= e %ѭ���ж�������
    Va0 = Va1;
    I = zeros(4,1);
    for k = 1:2
        for i = 1:3
            for j = 1:3
                if Se(i,j,k) < 0
                    Vr(i,j,k) = Va0(TriElement(k,i))-Va0(TriElement(k,j))-Vi(i,j,k);
                    Vc(i,j,k) = 2*Vr(i,j,k)*Y0/(Y0-Se(i,j));
                    Vi(i,j,k) = Vc(i,j,k) - Vr(i,j,k);
                    I(TriElement(k,i)) = I(TriElement(k,i))+2*Y0*Vi(i,j,k);
                else %���ڸ�����
                    I(TriElement(k,i)) = I(TriElement(k,i))+Se(i,j,k)*(Va0(TriElement(k,i))-Va0(TriElement(k,j)));
                end
            end
        end
    end
    Va1 = [0;0;10;0];
    I = I - Ym*Va1;
    Va1([2,4]) = Ym([2,4],[2,4])\I([2,4]); 
    m = m+1;
end
