%通过直接法求解二维轴对称稳态热场
%这是一个线性问题，包含第一类边界条件,第二类边界条件还不会处理（逃
%材料的热导率为52W/(m.K),边界的温度为293.15K，求解域2的单位体积发热为10000000W/m^3
%参考颜威利《电气工程电磁场数值分析》 P57 轴对称磁场有限元法,没成
%（自己推导出来的泛函应该怎么离散？？？）
%2018/12/10
%By QzLancer
%-------------------------------读取分网文件
clear all;
[Coor,VtxElement,VtxEntity,EdgElement,EdgEntity,TriElement,TriEntity] = readcomsol('mesh1.mphtxt');
%-------------------------------根据插值要求求解各单元的几何参数，目前只考虑内部单元
%导出所有三角形单元RZ轴的坐标
R = Coor(:,1);
Z = Coor(:,2);
TriR = R(TriElement);
TriZ = Z(TriElement);
%计算出所有单元的p,q,r和Area
p(:,1) = TriR(:,2).*TriZ(:,3) - TriZ(:,2).*TriR(:,3);
p(:,2) = TriR(:,3).*TriZ(:,1) - TriZ(:,3).*TriR(:,1);
p(:,3) = TriR(:,1).*TriZ(:,2) - TriZ(:,1).*TriR(:,2);
q(:,1) = TriZ(:,2) - TriZ(:,3);
q(:,2) = TriZ(:,3) - TriZ(:,1);
q(:,3) = TriZ(:,1) - TriZ(:,2);
r(:,1) = TriR(:,3) - TriR(:,2);
r(:,2) = TriR(:,1) - TriR(:,3);
r(:,3) = TriR(:,2) - TriR(:,1);
Area = (q(:,1).*r(:,2) - q(:,2).*r(:,1))/2;
% 三角形所有单元的重心处的半径
% TriRadius = 1.5./(1./(TriR(:,1)+TriR(:,2))+1./(TriR(:,2)+TriR(:,3))+1./(TriR(:,3)+TriR(:,1)));
TriRadius = (TriR(:,1)+TriR(:,2)+TriR(:,3))./3;
%-------------------------------求出负载和热导率
Cond = 52;
%等效热导率
EquCond = Cond./TriRadius;
%找出热源所在单元
Domain2 = find(TriEntity==2);
%给所有单元附上热源
Source = zeros(length(TriElement),1);
Source(Domain2) = 10000000;
%-------------------------------单元分析和总体合成
S = zeros(length(Coor));
F = zeros(length(Coor),1);
for k = 1:length(TriElement)
    for i = 1:3
        for j = 1:3
            Se= EquCond(k) * (r(k,i)*r(k,j) + q(k,i)*q(k,j)) / (4*Area(k));
            S(TriElement(k,i),TriElement(k,j)) = S(TriElement(k,i),TriElement(k,j)) + Se;
        end
        Fe = Source(k)*Area(k)/3;
        F(TriElement(k,i)) = F(TriElement(k,i)) + Fe;
    end
end
%-------------------------------采用直接法求解
%检索出边界点
EquTemp = zeros(length(Coor),1);
Boundary = find(Z==0 | Z==0.14 | R==0.1);
FreeNodes = find(~(Z==0 | Z==0.14 | R==0.1));

EquTemp(Boundary) = 293.15.*R(Boundary);
F1 = S(FreeNodes,:)*EquTemp;
F2 = F(FreeNodes) - F1;
EquTemp(FreeNodes) = S(FreeNodes,FreeNodes)\F2;
Temp = EquTemp./R;
%-------------------------------后处理
Interp1 = scatteredInterpolant(R,Z,Temp);
% tx = 0:1e-3:0.08;
% ty = 0:1e-3:0.14;
tx = 0.02:1e-3:0.1;
ty = 0:1e-3:0.14;
[qx,qy] = meshgrid(tx,ty);
qz = Interp1(qx,qy);
% mesh(qx,qy,qz);
subplot(1,2,2);
contourf(qx,qy,qz,20);colorbar;
hold on;
EdgNode1 = EdgElement(:, 1);
EdgNode2 = EdgElement(:, 2);
EdgCoor1 = Coor(EdgNode1, :);
EdgCoor2 = Coor(EdgNode2, :);
% hold on;
% plot(EdgCoor1(:, 1), EdgCoor1(:, 2), 'x');
% hold on;
% plot(EdgCoor2(:, 1), EdgCoor2(:, 2), 'x');
Edglength = length(EdgElement); 
for i = 1:Edglength
    plot([EdgCoor1(i, 1), EdgCoor2(i, 1)], [EdgCoor1(i, 2), EdgCoor2(i, 2)],'Color','k','Linewidth',1);
    hold on;
end
%-------------------------读取COMSOL计算出来的结果并进行后处理
[fileID, Errmessage] = fopen('comsolsolution.txt', 'r');
if fileID == -1
    disp(Errmessage);
end
for i=1:9
    fgetl(fileID);
end
comsoldata = fscanf(fileID,'%lf %lf %lf\n',[3,length(Coor)]);
comsoldata = comsoldata';
fclose(fileID);
Interp2 = scatteredInterpolant(comsoldata(:,1),comsoldata(:,2),comsoldata(:,3));
qz = Interp2(qx,qy);
subplot(1,2,1);
contourf(qx,qy,qz,20);colorbar;
%----------------------------求解结果和仿真结果差距很大，一些验证工作
EquTemp1 = comsoldata(:,3).*R;
FF= S*EquTemp1;