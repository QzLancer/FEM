%实际接触器模型计算，C++求出奇异，用这个来debug
%2019/5/9
%By QzLancer
%-------------------------------读取分网文件
clear all;
[Coor,VtxElement,VtxEntity,EdgElement,EdgEntity,TriElement,TriEntity] = readcomsol('mesh_contactor.mphtxt');
tic;
%-------------------------------根据插值要求求解各单元的几何参数
%导出所有三角形单元RZ轴的坐标
R = Coor(:,1);
Z = Coor(:,2);
TriR = R(TriElement);
TriZ = Z(TriElement);
%计算出所有单元的p,q,r和Area出错 solvecontactor (line 14)
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
TriRadius = (TriR(:,1)+TriR(:,2)+TriR(:,3))./3;
%计算出所有线单元的坐标
EdgR = R(EdgElement);
EdgZ = Z(EdgElement);
d = sqrt((EdgR(:, 1) - EdgR(:, 2)).^2 + (EdgZ(:, 1) - EdgZ(:, 2)).^2);
%计算所有边界单元的平均半径r
EdgRadius = (EdgR(:, 1) + EdgR(:, 2)) / 2;
%--------------------------------------------------
%热传导边界定义
h = zeros(length(EdgEntity),1);
Text = zeros(length(EdgEntity),1);
HeatConvBound = find(EdgEntity==2 | EdgEntity==9 | EdgEntity == 49);
h(HeatConvBound) = 20;
Text(HeatConvBound) = 293.15;
a = EdgElement(HeatConvBound,:);
%负载和热导率定义
cond = zeros(length(TriEntity), 1);
Source = zeros(length(TriEntity), 1);
m1 = find(TriEntity==10 | TriEntity==11);
Source(m1) = 500000;
cond(m1) = 400;
m1 = find(TriEntity == 2 | TriEntity == 3 | TriEntity == 4 | TriEntity == 6 | TriEntity == 8 );
cond(m1) = 76.2;
m1 =  find(TriEntity == 1 | TriEntity == 7 | TriEntity == 9);
cond(m1) = 0.26;
m1 = find(TriEntity == 5 | TriEntity == 12 | TriEntity == 13);
cond(m1) = 0.26;
%------------------------------面单元分析和装配
S = zeros(length(Coor));
F = zeros(length(Coor),1);
Se = zeros(3,3,length(Coor));
for k = 1:length(TriElement)
    for i = 1:3
        for j = 1:3
            Se(i,j,k)= (pi*cond(k)*TriRadius(k)*(r(k,i)*r(k,j) + q(k,i)*q(k,j)))/(2*Area(k));
%             Se= Cond .* (r(k,i)*r(k,j) + q(k,i)*q(k,j)) / (4*Area(k));
%             Se= (pi*Cond*sqrt(RR(k))*(r(k,i)*r(k,j) + q(k,i)*q(k,j)))/(2*Area(k));
            S(TriElement(k,i),TriElement(k,j)) = S(TriElement(k,i),TriElement(k,j)) + Se(i,j,k);
        end
        Fe = pi*Source(k)*Area(k)*(R(TriElement(k,i))+3*TriRadius(k))/6;
%         Fe = pi*Source(k)*Area(k)*cofF1(k,i)/15;
        F(TriElement(k,i)) = F(TriElement(k,i)) + Fe;
    end
end
S3 = S;
%------------------------------线单元分析和装配
Fl = zeros(length(EdgElement),2);
for k = 1:length(EdgElement);
    for i = 1:2
       for j = 1:2
           if i == j
               Sl =  pi*h(k)*d(k)*(2*EdgRadius(k)+2*R(EdgElement(k,i)))/6;
           else
               Sl =  pi*h(k)*d(k)*(2*EdgRadius(k))/6;
           end
               S(EdgElement(k,i), EdgElement(k,j)) = S(EdgElement(k,i), EdgElement(k,j)) + Sl;
       end
       Fl(k,i) = pi*h(k)*Text(k)*d(k)*(2*EdgRadius(k)+R(EdgElement(k,i)))/3;
       F(EdgElement(k,i)) = F(EdgElement(k,i)) + Fl(k,i);
    end
end
% %-------------------------------采用直接法求解
    Temp = S\F;
toc;
%-------------------------------后处理
Interp1 = scatteredInterpolant(R,Z,Temp);
% tx = 0:1e-3:0.08;
% ty = 0:1e-3:0.14;
tx = 0:1e-3:0.026;
ty = -0.025:1e-3:0.024;
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
[fileID, Errmessage] = fopen('solve.txt', 'r');
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